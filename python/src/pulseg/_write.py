"""Flattening a PulSeg representation back into Pulseq blocks."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pypulseqpp as pp

from ._blocks import get_grad_energy, shot_shape
from ._model import CONSTANT_DELAY_ID, GRADIENT_AXES, VARIABLE_DELAY_ID, PulseqBlock, Representation
from ._validate import validate_ir

PRESCRIPTION_FLAGS = ("NOROT", "NOPOS", "PMC")


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
        n_rf = n_grad = n_adc = n_trig = 0
        if inst.adc_acquire is not None and np.any(np.asarray(inst.adc_acquire) == 0):
            raise NotImplementedError(
                "adc_acquire has a readout that is played but not stored, which a Pulseq file cannot express"
            )
        for position, block_id in enumerate(segment.base_block_ids):
            duration = float(inst.block_duration[position])
            if block_id in (CONSTANT_DELAY_ID, VARIABLE_DELAY_ID):
                yield PulseqBlock(block_duration=duration), None, None, inst
                continue
            base = rep.base_block(block_id).block
            block = PulseqBlock(block_duration=duration, trig=[SimpleNamespace(**vars(t)) for t in base.trig])
            if any(t.type == "output" for t in base.trig):
                if inst.trig_enable is not None and inst.trig_enable[n_trig] == 0:
                    block.trig = [t for t in block.trig if t.type != "output"]
                n_trig += 1
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
                shot = int(inst.gradient_shot_index[n_grad]) if inst.gradient_shot_index is not None else 0
                for a, axis in enumerate(GRADIENT_AXES):
                    grad = getattr(base, axis)
                    if grad is None:
                        continue
                    event = SimpleNamespace(**vars(grad))
                    for name in ("shot_tt", "shot_waveforms", "shot_first", "shot_last"):
                        event.__dict__.pop(name, None)
                    if grad.type == "trap":
                        event.amplitude = grad.amplitude * scale[a]
                    else:
                        tt, waveform, first, last = shot_shape(grad, shot)
                        event.tt = np.asarray(tt)
                        event.shape_dur = float(event.tt[-1])
                        event.waveform = waveform * scale[a]
                        event.first = first * scale[a]
                        event.last = last * scale[a]
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


def write(rep: Representation, seqfile=None, *, strict: bool = False, verbose: bool = False, system=None):
    """Flatten a PulSeg representation into a Pulseq sequence, and optionally write it.

    Every block is emitted at its physical amplitudes and selected shot, with
    its rotation and RF shim as Pulseq extensions, and ``NOROT``, ``NOPOS``
    and ``PMC`` set where an instance changes them.

    Parameters
    ----------
    rep : Representation
    seqfile : str or pathlib.Path, optional
        Where to write the ``.seq`` file; nothing is written if omitted.
    strict : bool, default=False
        Validate the representation first.
    verbose : bool, default=False
        Print a summary.
    system : pypulseqpp.Opts, optional
        Limits and rasters of the written sequence; ``rep.system`` if
        omitted, which :func:`read` sets to the source's system.

    Returns
    -------
    pypulseqpp.Sequence

    Raises
    ------
    ValidationError
        With ``strict``, for an invalid representation.
    NotImplementedError
        Where ``adc_acquire`` is 0: Pulseq has no readout that is played but
        not stored.
    """
    if strict:
        validate_ir(rep)
    system = rep.system if system is None else system
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
    if seqfile is not None:
        seq.write(str(seqfile))
    if verbose:
        print(f"write: {seq.num_blocks} blocks" + (f" written to {seqfile}" if seqfile is not None else ""))
    return seq


def stream2loop(rep: Representation, dt: float | None = None) -> np.ndarray:
    """Return the execution stream as one row per Pulseq block.

    Parameters
    ----------
    rep : Representation
    dt : float, optional
        Accepted for parity with the MATLAB toolbox; gradient energies are
        integrated over each gradient's own sample times.

    Returns
    -------
    numpy.ndarray
        ``(n_blocks, 24)``. Columns: virtual segment ID, base-block ID, RF
        amplitude, phase and frequency, then for x, y and z the gradient
        scale and energy ((Hz/m)^2 s), ADC phase and frequency, block
        duration (s), physiological trigger (first block of an instance),
        and the rotation matrix, row-major.
    """
    del dt
    validate_ir(rep)
    rows = []
    for inst in rep.execution_stream:
        segment = rep.virtual_segment(inst.virtual_segment_id)
        n_rf = n_grad = n_adc = 0
        for position, block_id in enumerate(segment.base_block_ids):
            row = np.zeros(24)
            row[0], row[1] = segment.id, block_id
            row[13] = inst.block_duration[position]
            row[14] = inst.physio_trigger if position == 0 else 0
            row[15:24] = np.eye(3).ravel()
            if block_id not in (CONSTANT_DELAY_ID, VARIABLE_DELAY_ID):
                base = rep.base_block(block_id).block
                if base.rf is not None:
                    row[2:5] = (
                        inst.rf_amplitude[n_rf],
                        inst.rf_phase_offset[n_rf],
                        inst.rf_frequency_offset[n_rf],
                    )
                    n_rf += 1
                if base.has_gradient:
                    shot = (
                        int(inst.gradient_shot_index[n_grad]) if inst.gradient_shot_index is not None else 0
                    )
                    for a, axis in enumerate(GRADIENT_AXES):
                        grad = getattr(base, axis)
                        if grad is not None:
                            scale = inst.gradient_amplitude[n_grad][a]
                            row[5 + 2 * a] = scale
                            row[6 + 2 * a] = get_grad_energy(grad, scale, shot)
                    if inst.rotation_matrix is not None:
                        row[15:24] = np.asarray(inst.rotation_matrix[n_grad]).ravel()
                    n_grad += 1
                if base.adc is not None:
                    row[11:13] = inst.adc_phase_offset[n_adc], inst.adc_frequency_offset[n_adc]
                    n_adc += 1
            rows.append(row)
    return np.array(rows).reshape(-1, 24)
