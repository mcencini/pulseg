"""Conversion of a Pulseq sequence into a PulSeg representation."""

from __future__ import annotations

import datetime
from pathlib import Path

import numpy as np

from ._blocks import gradient_timing, normalize_block
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

PRESCRIPTION_FLAGS = ("NOROT", "NOPOS", "PMC")
KEY_DECIMALS = 9


def _pypulseqpp():
    try:
        import pypulseqpp
    except ImportError as error:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "Converting Pulseq sequences needs pypulseqpp: pip install 'pulseg[pulseq]'"
        ) from error
    return pypulseqpp


def quaternion_to_matrix(q) -> np.ndarray:
    """Return the rotation matrix of a scalar-first quaternion ``(w, x, y, z)``."""
    w, x, y, z = np.asarray(q, dtype=float) / np.linalg.norm(q)
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def _round(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            value = np.round(value.real, KEY_DECIMALS) + 1j * np.round(value.imag, KEY_DECIMALS)
        else:
            value = np.round(value, KEY_DECIMALS) + 0.0
        return (value.dtype.str, value.shape, value.tobytes())
    if isinstance(value, float):
        return round(value, KEY_DECIMALS) + 0.0
    return value


def _event_key(event):
    if event is None:
        return None
    return tuple(sorted((k, _round(v)) for k, v in vars(event).items()))


def _timing_key(block: PulseqBlock) -> tuple:
    """Identity of a base block: RF and ADC by shape, gradients by timing."""
    return (
        _round(block.block_duration),
        _event_key(block.rf),
        *(
            gradient_timing(getattr(block, a)) if getattr(block, a) is not None else None
            for a in GRADIENT_AXES
        ),
        _event_key(block.adc),
        tuple(_event_key(t) for t in block.trig),
    )


def _shot_key(block: PulseqBlock) -> tuple:
    """The normalized arbitrary-gradient shapes of a block, one entry per axis."""
    return tuple(
        _round(np.asarray(g.waveform)) if g is not None and g.type == "grad" else None
        for g in (getattr(block, a) for a in GRADIENT_AXES)
    )


def _add_shot(base: PulseqBlock, block: PulseqBlock) -> None:
    for axis in GRADIENT_AXES:
        stored, new = getattr(base, axis), getattr(block, axis)
        if stored is None or stored.type != "grad":
            continue
        if not hasattr(stored, "shot_waveforms"):
            stored.shot_waveforms = np.asarray(stored.waveform)[None, :]
            stored.shot_first = np.array([stored.first])
            stored.shot_last = np.array([stored.last])
        stored.shot_waveforms = np.vstack([stored.shot_waveforms, np.asarray(new.waveform)[None, :]])
        stored.shot_first = np.append(stored.shot_first, new.first)
        stored.shot_last = np.append(stored.shot_last, new.last)


def _resolve_soft_delay(b, soft_delay_input_ms, index):
    if soft_delay_input_ms is None:
        raise ValueError(
            "This Pulseq sequence contains soft delay events. Call pulseg.import_(..., "
            "soft_delay_input_ms=value) to specify the soft delay input in milliseconds."
        )
    duration = soft_delay_input_ms * 1e-3 / b.soft_delay.factor + b.soft_delay.offset
    if duration < 0:
        raise ValueError(f"Soft delay on block {index} evaluates to a negative duration ({duration:g} s).")
    b.soft_delay = None
    b.block_duration = duration
    return b


def _prescription_flags(seq, n_blocks: int) -> dict:
    labels = seq.evaluate_labels(init={f: 0 for f in PRESCRIPTION_FLAGS}, evolution="blocks")
    return {f: np.broadcast_to(np.asarray(labels[f], dtype=int), (n_blocks,)) for f in PRESCRIPTION_FLAGS}


def _instance_starts(seq, decoded: list, flags: dict) -> list[int]:
    n_blocks = len(decoded)
    size, start = seq.repetition()
    starts = set(range(start - 1, n_blocks, max(size, 1))) | {0}
    for i, block in enumerate(decoded):
        if any(t.type == "trigger" for t in (getattr(block, "trig", None) or [])):
            starts.add(i)
    for i in range(1, n_blocks):
        if any(flags[f][i] != flags[f][i - 1] for f in PRESCRIPTION_FLAGS):
            starts.add(i)
    return sorted(starts)


def import_(seq, *, soft_delay_input_ms: float | None = None, verbose: bool = False) -> Representation:
    """Convert a Pulseq sequence or ``.seq`` file into a PulSeg representation.

    The Python counterpart of the MATLAB ``pulseg.import`` (``import`` is a
    Python keyword). Segment boundaries need no ``TRID`` label: a segment
    instance starts at every repetition ``pypulseqpp.Sequence.repetition``
    reports, at every block carrying a physiological trigger, and wherever
    the ``NOROT``, ``NOPOS`` or ``PMC`` label changes. Instances with the same
    sequence of base blocks play one virtual segment.

    A base block is identified by the shapes of its RF and ADC events and by
    the timing of its gradients. Arbitrary gradients that share a timing but
    differ in shape after normalization (independently designed spiral
    interleaves, for instance) are shot variants of one base block, selected
    per instance by ``gradient_shot_index``. A pure delay is base block 1
    where its duration differs between the instances of its virtual segment,
    and base block 0 otherwise.

    Parameters
    ----------
    seq : pypulseqpp.Sequence or str or pathlib.Path
        The sequence, or the path of a ``.seq`` file.
    soft_delay_input_ms : float, optional
        Value entered for every soft delay, in ms; required if the sequence
        has soft delays. Each resolves to ``input / factor + offset``.
    verbose : bool, default=False
        Print a summary.

    Returns
    -------
    Representation
        Version 2.1-alpha, or 2.2-alpha where the sequence uses RF shims or
        the ``NOROT``, ``NOPOS`` or ``PMC`` labels.

    Raises
    ------
    ValueError
        For soft delays without ``soft_delay_input_ms``, or one that resolves
        to a negative duration.
    NotImplementedError
        For RF or ADC offsets given in ppm.
    """
    pp = _pypulseqpp()
    source_file = None
    if isinstance(seq, str | Path):
        source_file = str(seq)
        loaded = pp.Sequence()
        loaded.read(source_file)
        seq = loaded

    n_blocks = seq.num_blocks
    decoded = []
    for i in range(n_blocks):
        b = seq.get_block(i + 1)
        if getattr(b, "soft_delay", None) is not None:
            b = _resolve_soft_delay(b, soft_delay_input_ms, i + 1)
        decoded.append(b)
    rotations = seq.block_rotations()
    uses_shim = any(getattr(b, "rf_shim", None) is not None for b in decoded)
    flags = _prescription_flags(seq, n_blocks)
    starts = _instance_starts(seq, decoded, flags) + [n_blocks]

    base_ids: dict[tuple, int] = {}
    shot_ids: dict[tuple, dict[tuple, int]] = {}
    base_blocks: list[BaseBlock] = []
    segment_ids: dict[tuple, int] = {}
    segments: list[VirtualSegment] = []
    raw = []

    for first, stop in zip(starts[:-1], starts[1:], strict=True):
        ids, rf, grad, shots, rot, adc, shim, durations = [], [], [], [], [], [], [], []
        physio = 0
        for i in range(first, stop):
            b = decoded[i]
            b0, scales = normalize_block(b)
            durations.append(b0.block_duration)
            physio |= any(t.type == "trigger" for t in b0.trig)
            if b0.is_pure_delay:
                ids.append(-1)
                continue
            key = _timing_key(b0)
            if key not in base_ids:
                base_ids[key] = len(base_blocks) + 2
                base_blocks.append(BaseBlock(id=base_ids[key], block=b0))
                shot_ids[key] = {_shot_key(b0): 0}
            elif _shot_key(b0) not in shot_ids[key]:
                shot_ids[key][_shot_key(b0)] = len(shot_ids[key])
                _add_shot(base_blocks[base_ids[key] - 2].block, b0)
            ids.append(base_ids[key])
            if scales["rf"] is not None:
                rf.append(scales["rf"])
                rf_shim = getattr(b, "rf_shim", None)
                shim.append(np.asarray(rf_shim.shim_vector, dtype=complex) if rf_shim is not None else None)
            if b0.has_gradient:
                grad.append(scales["grad"])
                shots.append(shot_ids[key][_shot_key(b0)])
                q = b.rotation.quaternion if b.rotation is not None else None
                rot.append(quaternion_to_matrix(q) if rotations[i] and q is not None else np.eye(3))
            if scales["adc"] is not None:
                adc.append(scales["adc"])
        if tuple(ids) not in segment_ids:
            segment_ids[tuple(ids)] = len(segments) + 1
            segments.append(VirtualSegment(id=segment_ids[tuple(ids)], base_block_ids=list(ids)))
        raw.append((segment_ids[tuple(ids)], rf, grad, shots, rot, adc, shim, durations, physio, first))

    variable = {s.id: np.zeros(len(s.base_block_ids), bool) for s in segments}
    reference: dict[int, np.ndarray] = {}
    for segment_id, *_, durations, _, _ in raw:
        durations = np.asarray(durations)
        variable[segment_id] |= np.abs(durations - reference.setdefault(segment_id, durations)) > 1e-12
    for segment in segments:
        segment.base_block_ids = [
            (VARIABLE_DELAY_ID if variable[segment.id][p] else CONSTANT_DELAY_ID) if b == -1 else b
            for p, b in enumerate(segment.base_block_ids)
        ]

    uses_rotation = bool(np.any(rotations))
    uses_shots = any(len(v) > 1 for v in shot_ids.values())
    uses_flags = any(np.any(flags[f]) for f in PRESCRIPTION_FLAGS)
    stream = []
    for segment_id, rf, grad, shots, rot, adc, shim, durations, physio, first in raw:
        rf = np.asarray(rf, dtype=float).reshape(-1, 3)
        adc = np.asarray(adc, dtype=float).reshape(-1, 2)
        inst = SegmentInstance(
            virtual_segment_id=segment_id,
            rf_amplitude=rf[:, 0],
            rf_phase_offset=rf[:, 1],
            rf_frequency_offset=rf[:, 2],
            gradient_amplitude=np.asarray(grad, dtype=float).reshape(-1, 3),
            adc_phase_offset=adc[:, 0],
            adc_frequency_offset=adc[:, 1],
            block_duration=np.asarray(durations, dtype=float),
            gradient_shot_index=np.asarray(shots, dtype=int) if uses_shots else None,
            rotation_matrix=np.asarray(rot, dtype=float).reshape(-1, 3, 3) if uses_rotation else None,
            physio_trigger=int(physio),
            norot=int(flags["NOROT"][first]),
            nopos=int(flags["NOPOS"][first]),
            pmc=int(flags["PMC"][first]),
        )
        if uses_shim:
            width = max((len(s) for s in shim if s is not None), default=1)
            weights = [s if s is not None else np.ones(width) for s in shim]
            inst.rf_shim = np.array(weights, dtype=complex).reshape(-1, width)
        stream.append(inst)

    rep = Representation(
        pulseg_version="2.2-alpha" if (uses_shim or uses_flags) else "2.1-alpha",
        base_blocks=base_blocks,
        virtual_segments=segments,
        execution_stream=stream,
        source_file=source_file,
        creation_date=datetime.date.today().isoformat(),
    )
    if verbose:
        print(
            f"import_: {n_blocks} blocks, {len(base_blocks)} base blocks, "
            f"{len(segments)} virtual segments, {len(stream)} segment instances"
        )
    return rep
