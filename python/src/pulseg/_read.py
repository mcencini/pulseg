"""Reading a Pulseq sequence into a PulSeg representation through pulserver."""

from __future__ import annotations

import datetime
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pulserver.ir as ir
import pypulseqpp as pp

from ._model import (
    CONSTANT_DELAY_ID,
    GRADIENT_AXES,
    VARIABLE_DELAY_ID,
    BaseBlock,
    PulseqBlock,
    Representation,
    SegmentInstance,
    VirtualSegment,
)

_RF_USE = {0: "undefined", 1: "excitation", 2: "refocusing"}


def _same(a: np.ndarray, b: np.ndarray, tol: float = 1e-6) -> bool:
    return a.shape == b.shape and bool(np.allclose(a, b, rtol=0, atol=tol, equal_nan=True))


class _Position:
    """One block position of a virtual segment, as the instances play it."""

    def __init__(self):
        self.base: PulseqBlock | None = None
        self.shots: list[tuple[np.ndarray, np.ndarray]] = []

    def shot(self, tt: np.ndarray, shape: np.ndarray) -> int:
        for n, (t, s) in enumerate(self.shots):
            if _same(t, tt) and _same(s, shape):
                return n
        self.shots.append((tt, shape))
        return len(self.shots) - 1


def _rf_event(played, b, system) -> SimpleNamespace:
    first, stop = played["rf_span"][b]
    delay = played["rf_delay_us"][b] * 1e-6
    signal = np.asarray(played["rf_waveform_hz"][first:stop], dtype=complex)
    signal = signal / played["rf_amp_hz"][b]
    t = played["rf_time_us"][first:stop] * 1e-6 - delay
    return SimpleNamespace(
        type="rf",
        signal=signal,
        t=t,
        shape_dur=len(signal) * system.rf_raster_time,
        delay=delay,
        center=played["rf_center_us"][b] * 1e-6 - delay,
        dead_time=system.rf_dead_time,
        ringdown_time=system.rf_ringdown_time,
        use=_RF_USE.get(int(played["rf_use"][b]), "undefined"),
    )


def _adc_event(played, b, system) -> SimpleNamespace:
    return SimpleNamespace(
        type="adc",
        num_samples=int(played["adc_samples"][b]),
        dwell=played["adc_dwell_ns"][b] * 1e-9,
        delay=played["adc_delay_us"][b] * 1e-6,
        dead_time=system.adc_dead_time,
        phase_modulation=np.zeros(0),
    )


def _gradient_shapes(played, b) -> tuple[np.ndarray, np.ndarray]:
    """Corner times (s, from the block start) and normalized values, ``(3, n)`` each; NaN off-axis."""
    spans = played["gradient_span"][b]
    n = max(int(stop - first) for first, stop in spans)
    tt = np.full((3, n), np.nan)
    shape = np.full((3, n), np.nan)
    for a, (first, stop) in enumerate(spans):
        if stop > first:
            amplitude = played["gradient_hz_per_m"][b][a]
            values = played["gradient_waveform_hz_per_m"][first:stop]
            tt[a, : stop - first] = played["gradient_time_us"][first:stop] * 1e-6
            shape[a, : stop - first] = values / amplitude if amplitude != 0 else values
    return tt, shape


def _plays_gradient(played, b) -> bool:
    spans = played["gradient_span"][b]
    return bool(np.any(spans[:, 1] > spans[:, 0]))


def _gradient_events(tt, shape, shots) -> dict:
    events = {}
    for a, axis in enumerate(GRADIENT_AXES):
        valid = ~np.isnan(tt[a])
        if not valid.any():
            continue
        t = tt[a, valid]
        event = SimpleNamespace(type="grad", channel=axis[1], delay=float(t[0]), tt=t - t[0])
        event.shape_dur = float(event.tt[-1])
        event.waveform = shape[a, valid]
        event.first, event.last = float(event.waveform[0]), float(event.waveform[-1])
        if len(shots) > 1:
            event.shot_tt = [s_tt[a, ~np.isnan(s_tt[a])] - t[0] for s_tt, _ in shots]
            event.shot_waveforms = [s[a, ~np.isnan(s[a])] for _, s in shots]
            event.shot_first = np.array([w[0] for w in event.shot_waveforms])
            event.shot_last = np.array([w[-1] for w in event.shot_waveforms])
        events[axis] = event
    return events


_RASTERS = {
    "rf_raster_time": "RadiofrequencyRasterTime",
    "grad_raster_time": "GradientRasterTime",
    "adc_raster_time": "AdcRasterTime",
    "block_duration_raster": "BlockDurationRaster",
}


def sequence_system(seq) -> pp.Opts:
    """Return the system a sequence was built with.

    Parameters
    ----------
    seq : pypulseqpp.Sequence

    Returns
    -------
    pypulseqpp.Opts
        ``seq.system`` if set; otherwise the rasters the sequence's
        ``[DEFINITIONS]`` record, on pypulseqpp's defaults for the rest.
    """
    if seq.system is not None:
        return seq.system
    rasters = {}
    for field, name in _RASTERS.items():
        value = seq.get_definition(name)
        if value:
            rasters[field] = float(value[0])
    return pp.Opts(**rasters)


def read(seq, system=None, *, verbose: bool = False) -> Representation:
    """Read a Pulseq sequence or ``.seq`` file into a PulSeg representation.

    pulserver segments the sequence (``pulserver.ir.convert``) and plays its
    cache back (``pulserver.ir.play``); each of pulserver's segments is one
    virtual segment, and each of its block positions one base block. Segment
    boundaries need no ``TRID`` label; ``pulserver.ir.Grouping`` documents
    where they fall.

    A gradient that a position plays with another shape or timing in some
    instances, such as a spiral interleaf or a rotated readout, is a shot
    variant of that position's base block. Rotations reach the gradient
    shapes this way, so ``rotation_matrix`` is left unset.

    Parameters
    ----------
    seq : pypulseqpp.Sequence or str or pathlib.Path
        The sequence, or the path of a ``.seq`` file.
    system : pypulseqpp.Opts, optional
        Limits and rasters the sequence is segmented under; the sequence's
        own (:func:`sequence_system`) if omitted.
    verbose : bool, default=False
        Print a summary.

    Returns
    -------
    Representation
        Version 2.1-alpha, or 2.2-alpha where an instance sets ``NOROT`` or
        ``NOPOS``.

    Raises
    ------
    ValueError
        If pulserver cannot segment the sequence.
    """
    source_file = None
    if isinstance(seq, str | Path):
        source_file = str(seq)
        loaded = pp.Sequence()
        loaded.read(source_file)
    else:
        loaded = seq
    system = sequence_system(loaded) if system is None else system
    with tempfile.TemporaryDirectory() as folder:
        if source_file is not None:
            path = Path(source_file)
        else:
            path = Path(folder) / "sequence.seq"
            seq.write(str(path))
        cache = ir.convert(path, system, cache_ext=".pulseg-read")
        try:
            summary = ir.summary(path, system, cache_ext=".pulseg-read")
            played = ir.play(path, ".pulseg-read", waveforms=True)
        finally:
            cache.unlink(missing_ok=True)
    rep = _representation(summary["segments"], played, system)
    rep.source_file = source_file
    rep.system = system
    if verbose:
        print(
            f"read: {len(played['segment'])} blocks, {len(rep.base_blocks)} base blocks, "
            f"{len(rep.virtual_segments)} virtual segments, {len(rep.execution_stream)} segment instances"
        )
    return rep


def _representation(segments, played, system) -> Representation:
    n_blocks = len(played["segment"])
    positions = {s: [_Position() for _ in range(int(seg["num_blocks"]))] for s, seg in enumerate(segments)}
    instances = []
    b = 0
    while b < n_blocks:
        s = int(played["segment"][b])
        length = int(segments[s]["num_blocks"])
        blocks, shots = range(b, b + length), []
        for p, k in enumerate(blocks):
            position = positions[s][p]
            if position.base is None:
                position.base = PulseqBlock(block_duration=0.0)
                if played["rf_channels"][k]:
                    position.base.rf = _rf_event(played, k, system)
                if played["adc"][k]:
                    position.base.adc = _adc_event(played, k, system)
            if _plays_gradient(played, k):
                shots.append(position.shot(*_gradient_shapes(played, k)))
        instances.append((s, blocks, shots, bool(segments[s]["has_trigger"])))
        b += length

    base_blocks, ids = [], {}
    for s, segment_positions in positions.items():
        plays = [blocks for t, blocks, *_ in instances if t == s]
        durations = np.array([[played["duration_us"][k] for k in blocks] for blocks in plays])
        for p, position in enumerate(segment_positions):
            block = position.base
            if position.shots:
                for axis, event in _gradient_events(*position.shots[0], position.shots).items():
                    setattr(block, axis, event)
            if block is None or block.is_pure_delay:
                variable = durations.size and np.ptp(durations[:, p]) > 0
                ids[s, p] = VARIABLE_DELAY_ID if variable else CONSTANT_DELAY_ID
                continue
            block.block_duration = float(durations[0, p]) * 1e-6
            ids[s, p] = len(base_blocks) + 2
            base_blocks.append(BaseBlock(id=ids[s, p], block=block))

    virtual_segments = [
        VirtualSegment(id=s + 1, base_block_ids=[ids[s, p] for p in range(len(segment_positions))])
        for s, segment_positions in positions.items()
    ]
    uses_shots = any(len(p.shots) > 1 for ps in positions.values() for p in ps)
    stream = []
    for s, blocks, shots, physio in instances:
        rf = [k for k in blocks if played["rf_channels"][k]]
        grad = [k for k in blocks if _plays_gradient(played, k)]
        adc = [k for k in blocks if played["adc"][k]]
        stream.append(
            SegmentInstance(
                virtual_segment_id=s + 1,
                rf_amplitude=np.asarray(played["rf_amp_hz"][rf], dtype=float),
                rf_phase_offset=np.asarray(played["rf_phase_rad"][rf], dtype=float),
                rf_frequency_offset=np.asarray(played["rf_freq_hz"][rf], dtype=float),
                gradient_amplitude=np.asarray(played["gradient_hz_per_m"][grad], dtype=float).reshape(-1, 3),
                adc_phase_offset=np.asarray(played["adc_phase_rad"][adc], dtype=float),
                adc_frequency_offset=np.asarray(played["adc_freq_hz"][adc], dtype=float),
                block_duration=np.asarray(played["duration_us"][list(blocks)], dtype=float) * 1e-6,
                gradient_shot_index=np.asarray(shots, dtype=int) if uses_shots else None,
                physio_trigger=int(physio),
                norot=int(played["norot"][blocks.start]),
                nopos=int(played["nopos"][blocks.start]),
            )
        )
    uses_flags = any(i.norot or i.nopos for i in stream)
    return Representation(
        pulseg_version="2.2-alpha" if uses_flags else "2.1-alpha",
        base_blocks=base_blocks,
        virtual_segments=virtual_segments,
        execution_stream=stream,
        creation_date=datetime.date.today().isoformat(),
    )
