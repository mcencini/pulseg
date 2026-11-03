"""Conversion between Pulseq sequences (through pypulseqpp) and PulSeg."""

from __future__ import annotations

import datetime
from types import SimpleNamespace

import numpy as np

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
_RF_FIELDS = ("t", "shape_dur", "delay", "center", "dead_time", "ringdown_time", "use")
_ADC_FIELDS = ("num_samples", "dwell", "delay", "dead_time")
_TRAP_FIELDS = ("delay", "rise_time", "flat_time", "fall_time")
_GRAD_FIELDS = ("tt", "delay", "shape_dur")
_TRIG_FIELDS = ("type", "channel", "delay", "duration")


def _pypulseqpp():
    try:
        import pypulseqpp
    except ImportError as error:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "Converting Pulseq sequences needs pypulseqpp: pip install 'pulseg[pulseq]'"
        ) from error
    return pypulseqpp


def quaternion_to_matrix(q) -> np.ndarray:
    """Return the rotation matrix of a scalar-first unit quaternion.

    Parameters
    ----------
    q : array_like
        ``(w, x, y, z)``; normalised here.

    Returns
    -------
    numpy.ndarray
        ``(3, 3)`` rotation matrix.
    """
    w, x, y, z = np.asarray(q, dtype=float) / np.linalg.norm(q)
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def _copy(event, fields):
    return SimpleNamespace(
        **{
            f: np.array(getattr(event, f)) if isinstance(getattr(event, f), np.ndarray) else getattr(event, f)
            for f in fields
        }
    )


def normalize_block(block) -> tuple[PulseqBlock, dict]:
    """Split a decoded Pulseq block into its normalized shape and its scales.

    RF is divided by its peak magnitude. Each gradient axis is divided by its
    signed amplitude at the sample of largest magnitude, so a polarity flip is
    a negative scale on the same base block. RF and ADC offsets are moved out
    of the block.

    Parameters
    ----------
    block : SimpleNamespace
        A block as ``pypulseqpp.Sequence.get_block`` decodes it.

    Returns
    -------
    base : PulseqBlock
        The normalized block.
    scales : dict
        ``rf`` (amplitude, phase, frequency) or None, ``grad`` (three signed
        scales) and ``adc`` (phase, frequency) or None.

    Raises
    ------
    NotImplementedError
        For soft delays and for RF or ADC offsets given in ppm.
    """
    if getattr(block, "soft_delay", None) is not None:
        raise NotImplementedError(
            "soft delays are not converted; resolve them with Sequence.apply_soft_delay first"
        )
    base = PulseqBlock(block_duration=float(block.block_duration))
    scales = {"rf": None, "grad": np.zeros(3), "adc": None}

    if block.rf is not None:
        rf = block.rf
        if rf.freq_ppm or rf.phase_ppm:
            raise NotImplementedError("RF offsets in ppm are not converted")
        signal = np.asarray(rf.signal, dtype=complex)
        peak = float(np.max(np.abs(signal))) if signal.size else 0.0
        base.rf = _copy(rf, _RF_FIELDS)
        base.rf.type = "rf"
        base.rf.signal = signal / peak if peak > 0 else signal.copy()
        scales["rf"] = (peak, float(rf.phase_offset), float(rf.freq_offset))

    for a, axis in enumerate(GRADIENT_AXES):
        grad = getattr(block, axis)
        if grad is None:
            continue
        if grad.type == "trap":
            scale = float(grad.amplitude)
            event = _copy(grad, _TRAP_FIELDS)
            event.amplitude = 1.0
        else:
            waveform = np.asarray(grad.waveform, dtype=float)
            scale = float(waveform[np.argmax(np.abs(waveform))]) if waveform.size else 0.0
            event = _copy(grad, _GRAD_FIELDS)
            event.waveform = waveform / scale if scale != 0 else waveform.copy()
            event.first = float(grad.first) / scale if scale != 0 else 0.0
            event.last = float(grad.last) / scale if scale != 0 else 0.0
        event.type = grad.type
        event.channel = axis[1]
        setattr(base, axis, event)
        scales["grad"][a] = scale

    if block.adc is not None:
        adc = block.adc
        if adc.freq_ppm or adc.phase_ppm:
            raise NotImplementedError("ADC offsets in ppm are not converted")
        base.adc = _copy(adc, _ADC_FIELDS)
        base.adc.type = "adc"
        base.adc.num_samples = int(adc.num_samples)
        base.adc.phase_modulation = np.array(adc.phase_modulation, dtype=float)
        scales["adc"] = (float(adc.phase_offset), float(adc.freq_offset))

    base.trig = [_copy(t, _TRIG_FIELDS) for t in (getattr(block, "trig", None) or [])]
    return base, scales


def _round(value):
    if isinstance(value, np.ndarray):
        array = np.round(value, KEY_DECIMALS) + 0.0
        if np.iscomplexobj(array):
            array = np.round(array.real, KEY_DECIMALS) + 1j * np.round(array.imag, KEY_DECIMALS)
        return (array.dtype.str, array.shape, array.tobytes())
    if isinstance(value, float):
        return round(value, KEY_DECIMALS) + 0.0
    return value


def _event_key(event):
    if event is None:
        return None
    return tuple(sorted((k, _round(v)) for k, v in vars(event).items()))


def block_key(block: PulseqBlock) -> tuple:
    """Return a hashable key equal for blocks with the same normalized shape and timing."""
    return (
        _round(block.block_duration),
        _event_key(block.rf),
        *(_event_key(getattr(block, axis)) for axis in GRADIENT_AXES),
        _event_key(block.adc),
        tuple(_event_key(t) for t in block.trig),
    )


def _instance_starts(seq, n_blocks: int, decoded: list) -> list[int]:
    size, start = seq.repetition()
    starts = set(range(start - 1, n_blocks, max(size, 1)))
    starts.add(0)
    for i, block in enumerate(decoded):
        if any(t.type == "trigger" for t in (getattr(block, "trig", None) or [])):
            starts.add(i)
    flags = _prescription_flags(seq, n_blocks)
    for i in range(1, n_blocks):
        if any(flags[f][i] != flags[f][i - 1] for f in PRESCRIPTION_FLAGS):
            starts.add(i)
    return sorted(starts)


def _prescription_flags(seq, n_blocks: int) -> dict:
    labels = seq.evaluate_labels(init={f: 0 for f in PRESCRIPTION_FLAGS}, evolution="blocks")
    return {f: np.broadcast_to(np.asarray(labels[f], dtype=int), (n_blocks,)) for f in PRESCRIPTION_FLAGS}


def from_sequence(seq, *, source_file: str | None = None) -> Representation:
    """Convert a Pulseq sequence into a PulSeg representation.

    Segment boundaries follow this rule. A segment instance starts at every
    repetition that ``pypulseqpp.Sequence.repetition`` reports, at every
    block carrying a physiological trigger, and wherever the ``NOROT``,
    ``NOPOS`` or ``PMC`` label changes. Instances with the same sequence of
    base blocks play one virtual segment. A pure delay is base block 1 where
    its duration differs between the instances of its virtual segment, and
    base block 0 otherwise. No ``TRID`` label is needed.

    Parameters
    ----------
    seq : pypulseqpp.Sequence
        The sequence to convert.
    source_file : str, optional
        Recorded in the representation.

    Returns
    -------
    Representation
        Version 2.1-alpha, or 2.2-alpha where the sequence uses RF shims or
        the ``NOROT``, ``NOPOS`` or ``PMC`` labels.

    Raises
    ------
    NotImplementedError
        For soft delays and for offsets given in ppm.
    """
    _pypulseqpp()
    n_blocks = seq.num_blocks
    decoded = [seq.get_block(i + 1) for i in range(n_blocks)]
    rotations = seq.block_rotations()
    uses_shim = any(getattr(block, "rf_shim", None) is not None for block in decoded)
    flags = _prescription_flags(seq, n_blocks)
    starts = _instance_starts(seq, n_blocks, decoded) + [n_blocks]

    base_by_key: dict[tuple, int] = {}
    base_blocks: list[BaseBlock] = []
    segment_by_key: dict[tuple, int] = {}
    segments: list[VirtualSegment] = []
    raw_instances = []

    for first, stop in zip(starts[:-1], starts[1:], strict=True):
        ids, rf, grad, rot, adc, shim, durations = [], [], [], [], [], [], []
        physio = 0
        for i in range(first, stop):
            block = decoded[i]
            base, scales = normalize_block(block)
            durations.append(base.block_duration)
            physio |= any(t.type == "trigger" for t in base.trig)
            if base.is_pure_delay:
                ids.append(-1)
                continue
            key = block_key(base)
            if key not in base_by_key:
                base_by_key[key] = len(base_blocks) + 2
                base_blocks.append(BaseBlock(id=base_by_key[key], block=base))
            ids.append(base_by_key[key])
            if scales["rf"] is not None:
                rf.append(scales["rf"])
                rf_shim = getattr(block, "rf_shim", None)
                shim.append(np.asarray(rf_shim.shim_vector, dtype=complex) if rf_shim is not None else None)
            if base.has_gradient:
                grad.append(scales["grad"])
                q = block.rotation.quaternion if block.rotation is not None else None
                rot.append(quaternion_to_matrix(q) if rotations[i] and q is not None else np.eye(3))
            if scales["adc"] is not None:
                adc.append(scales["adc"])
        segment_key = tuple(ids)
        if segment_key not in segment_by_key:
            segment_by_key[segment_key] = len(segments) + 1
            segments.append(VirtualSegment(id=segment_by_key[segment_key], base_block_ids=list(ids)))
        raw_instances.append(
            (segment_by_key[segment_key], rf, grad, rot, adc, shim, durations, physio, first)
        )

    variable = {s.id: np.zeros(len(s.base_block_ids), bool) for s in segments}
    reference: dict[int, np.ndarray] = {}
    for segment_id, *_, durations, _, _ in raw_instances:
        durations = np.asarray(durations)
        ref = reference.setdefault(segment_id, durations)
        variable[segment_id] |= np.abs(durations - ref) > 1e-12
    for segment in segments:
        segment.base_block_ids = [
            (VARIABLE_DELAY_ID if variable[segment.id][p] else CONSTANT_DELAY_ID) if b == -1 else b
            for p, b in enumerate(segment.base_block_ids)
        ]

    uses_rotation = bool(np.any(rotations))
    uses_flags = any(np.any(flags[f]) for f in PRESCRIPTION_FLAGS)
    stream = []
    for segment_id, rf, grad, rot, adc, shim, durations, physio, first in raw_instances:
        rf = np.asarray(rf, dtype=float).reshape(-1, 3)
        adc = np.asarray(adc, dtype=float).reshape(-1, 2)
        instance = SegmentInstance(
            virtual_segment_id=segment_id,
            rf_amplitude=rf[:, 0],
            rf_phase_offset=rf[:, 1],
            rf_frequency_offset=rf[:, 2],
            gradient_amplitude=np.asarray(grad, dtype=float).reshape(-1, 3),
            adc_phase_offset=adc[:, 0],
            adc_frequency_offset=adc[:, 1],
            block_duration=np.asarray(durations, dtype=float),
            rotation_matrix=np.asarray(rot, dtype=float).reshape(-1, 3, 3) if uses_rotation else None,
            physio_trigger=int(physio),
            norot=int(flags["NOROT"][first]),
            nopos=int(flags["NOPOS"][first]),
            pmc=int(flags["PMC"][first]),
        )
        if uses_shim:
            width = max((len(s) for s in shim if s is not None), default=1)
            instance.rf_shim = np.array(
                [s if s is not None else np.ones(width) for s in shim], dtype=complex
            ).reshape(-1, width)
        stream.append(instance)

    return Representation(
        pulseg_version="2.2-alpha" if (uses_shim or uses_flags) else "2.1-alpha",
        base_blocks=base_blocks,
        virtual_segments=segments,
        execution_stream=stream,
        source_file=source_file,
        creation_date=datetime.date.today().isoformat(),
    )


def read_seq(path: str) -> Representation:
    """Read a Pulseq ``.seq`` file into a PulSeg representation.

    Parameters
    ----------
    path : str
        The ``.seq`` file.

    Returns
    -------
    Representation
        As :func:`from_sequence` converts it.
    """
    pp = _pypulseqpp()
    seq = pp.Sequence()
    seq.read(path)
    return from_sequence(seq, source_file=str(path))


def physical_blocks(rep: Representation):
    """Yield every block of the execution stream at its physical amplitudes.

    Applies each instance's RF amplitude, phase and frequency, gradient
    scales, rotation, ADC offsets and RF shim to its base blocks.

    Parameters
    ----------
    rep : Representation

    Yields
    ------
    block : PulseqBlock
        The block as played, gradients on the logical axes before rotation.
    rotation : numpy.ndarray or None
        ``(3, 3)`` rotation of a gradient block, else None.
    shim : numpy.ndarray or None
        Per-channel weights of an RF block, else None.
    instance : SegmentInstance
        The instance the block belongs to.
    """
    for inst in rep.execution_stream:
        segment = rep.virtual_segment(inst.virtual_segment_id)
        n_rf = n_grad = n_adc = 0
        for position, block_id in enumerate(segment.base_block_ids):
            duration = float(inst.block_duration[position])
            if block_id in (CONSTANT_DELAY_ID, VARIABLE_DELAY_ID):
                yield PulseqBlock(block_duration=duration), None, None, inst
                continue
            base = rep.base_block(block_id).block
            block = PulseqBlock(block_duration=duration, trig=[SimpleNamespace(**vars(t)) for t in base.trig])
            shim = rotation = None
            if base.rf is not None:
                rf = SimpleNamespace(**vars(base.rf))
                rf.signal = base.rf.signal * inst.rf_amplitude[n_rf]
                rf.phase_offset = float(inst.rf_phase_offset[n_rf])
                rf.freq_offset = float(inst.rf_frequency_offset[n_rf])
                rf.freq_ppm = rf.phase_ppm = 0.0
                if inst.rf_shim is not None:
                    shim = np.asarray(inst.rf_shim[n_rf])
                block.rf = rf
                n_rf += 1
            if base.has_gradient:
                scale = inst.gradient_amplitude[n_grad]
                for a, axis in enumerate(GRADIENT_AXES):
                    grad = getattr(base, axis)
                    if grad is None:
                        continue
                    event = SimpleNamespace(**vars(grad))
                    if grad.type == "trap":
                        event.amplitude = grad.amplitude * scale[a]
                    else:
                        event.waveform = grad.waveform * scale[a]
                        event.first = grad.first * scale[a]
                        event.last = grad.last * scale[a]
                    setattr(block, axis, event)
                if inst.rotation_matrix is not None:
                    rotation = np.asarray(inst.rotation_matrix[n_grad])
                n_grad += 1
            if base.adc is not None:
                adc = SimpleNamespace(**vars(base.adc))
                adc.phase_offset = float(inst.adc_phase_offset[n_adc])
                adc.freq_offset = float(inst.adc_frequency_offset[n_adc])
                adc.freq_ppm = adc.phase_ppm = 0.0
                block.adc = adc
                n_adc += 1
            yield block, rotation, shim, inst


def to_sequence(rep: Representation, system=None):
    """Write a PulSeg representation back to a Pulseq sequence.

    Every block is emitted at its physical amplitudes, with its rotation and
    RF shim as Pulseq extensions, and ``NOROT``, ``NOPOS`` and ``PMC`` set
    where an instance changes them.

    Parameters
    ----------
    rep : Representation
    system : pypulseqpp.Opts, optional
        Scanner limits; pypulseqpp's defaults if omitted.

    Returns
    -------
    pypulseqpp.Sequence
    """
    pp = _pypulseqpp()
    seq = pp.Sequence(system) if system is not None else pp.Sequence()
    state = {f: 0 for f in PRESCRIPTION_FLAGS}
    previous = None
    for block, rotation, shim, inst in physical_blocks(rep):
        events = [e for e in (block.rf, block.gx, block.gy, block.gz, block.adc) if e is not None]
        events += block.trig
        if inst is not previous:
            for flag in PRESCRIPTION_FLAGS:
                value = getattr(inst, flag.lower())
                if value != state[flag]:
                    events.append(pp.make_label(flag, "SET", bool(value)))
                    state[flag] = value
            previous = inst
        if rotation is not None and not np.allclose(rotation, np.eye(3)):
            events.append(pp.make_rotation(rotation))
        if shim is not None:
            events.append(pp.make_rf_shim(shim))
        if not events:
            events.append(pp.make_delay(block.block_duration))
        seq.add_block(*events)
        seq.block_durations[seq.num_blocks] = block.block_duration
    return seq
