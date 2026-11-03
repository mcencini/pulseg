"""Normalizing, comparing and checking Pulseq blocks as PulSeg base blocks."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace

import numpy as np

from ._model import GRADIENT_AXES, PulseqBlock

NORMALIZATION_TOLERANCE = 1e-6
_RF_FIELDS = ("t", "shape_dur", "delay", "center", "dead_time", "ringdown_time", "use")
_ADC_FIELDS = ("num_samples", "dwell", "delay", "dead_time")
_TRAP_FIELDS = ("delay", "rise_time", "flat_time", "fall_time")
_GRAD_FIELDS = ("tt", "delay", "shape_dur")
_TRIG_FIELDS = ("type", "channel", "delay", "duration")


def _copy(event, fields):
    values = {}
    for name in fields:
        value = getattr(event, name)
        values[name] = np.array(value) if isinstance(value, np.ndarray) else value
    return SimpleNamespace(**values)


def normalize_block(b) -> tuple[PulseqBlock, dict]:
    """Normalize the RF and gradient amplitudes of a Pulseq block.

    RF is divided by its peak magnitude. Each gradient axis is divided by its
    signed amplitude at the sample of largest magnitude, so a polarity flip is
    a negative scale of the same base block; a trapezoid takes amplitude 1. RF
    and ADC phase and frequency offsets are moved out of the block.

    Parameters
    ----------
    b : SimpleNamespace
        A block as ``pypulseqpp.Sequence.get_block`` decodes it.

    Returns
    -------
    b0 : PulseqBlock
        The normalized block.
    scales : dict
        ``rf``: (amplitude, phase, frequency) or None; ``grad``: the three
        signed scales; ``adc``: (phase, frequency) or None.

    Raises
    ------
    NotImplementedError
        For RF or ADC offsets given in ppm.
    """
    b0 = PulseqBlock(block_duration=float(b.block_duration))
    scales = {"rf": None, "grad": np.zeros(3), "adc": None}

    if b.rf is not None:
        rf = b.rf
        if rf.freq_ppm or rf.phase_ppm:
            raise NotImplementedError("RF offsets in ppm are not converted")
        signal = np.asarray(rf.signal, dtype=complex)
        peak = float(np.max(np.abs(signal))) if signal.size else 0.0
        b0.rf = _copy(rf, _RF_FIELDS)
        b0.rf.type = "rf"
        b0.rf.signal = signal / peak if peak > 0 else signal.copy()
        scales["rf"] = (peak, float(rf.phase_offset), float(rf.freq_offset))

    for a, axis in enumerate(GRADIENT_AXES):
        grad = getattr(b, axis)
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
            divisor = scale if scale != 0 else 1.0
            event.waveform = waveform / divisor
            event.first = float(grad.first) / divisor
            event.last = float(grad.last) / divisor
        event.type = grad.type
        event.channel = axis[1]
        setattr(b0, axis, event)
        scales["grad"][a] = scale

    if b.adc is not None:
        adc = b.adc
        if adc.freq_ppm or adc.phase_ppm:
            raise NotImplementedError("ADC offsets in ppm are not converted")
        b0.adc = _copy(adc, _ADC_FIELDS)
        b0.adc.type = "adc"
        b0.adc.num_samples = int(adc.num_samples)
        b0.adc.phase_modulation = np.array(adc.phase_modulation, dtype=float)
        scales["adc"] = (float(adc.phase_offset), float(adc.freq_offset))

    b0.trig = [_copy(t, _TRIG_FIELDS) for t in (getattr(b, "trig", None) or [])]
    return b0, scales


def num_shots(block: PulseqBlock) -> int:
    """Number of gradient shot variants a base block carries (1 without variants)."""
    counts = [
        len(getattr(getattr(block, axis), "shot_waveforms", [None]))
        for axis in GRADIENT_AXES
        if getattr(block, axis) is not None
    ]
    return max(counts, default=1)


def shot_shape(event, shot: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Return the corner times, normalized waveform, first and last value of one shot.

    Parameters
    ----------
    event : SimpleNamespace
        An arbitrary gradient (``type == "grad"``) of a base block.
    shot : int
        Shot variant, 0 for the first.

    Returns
    -------
    tt, waveform : numpy.ndarray
    first, last : float
    """
    if hasattr(event, "shot_waveforms"):
        tt = event.shot_tt[shot] if hasattr(event, "shot_tt") else event.tt
        return tt, event.shot_waveforms[shot], float(event.shot_first[shot]), float(event.shot_last[shot])
    if shot != 0:
        raise IndexError(f"gradient has one shot variant, shot {shot} requested")
    return event.tt, event.waveform, float(event.first), float(event.last)


def _same(a, b, tol) -> bool:
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a, b = np.asarray(a), np.asarray(b)
        return a.shape == b.shape and (a.size == 0 or np.allclose(a, b, rtol=0, atol=tol))
    if isinstance(a, float) or isinstance(b, float):
        return isinstance(a, int | float) and isinstance(b, int | float) and abs(a - b) <= tol
    return a == b


def _same_event(e1, e2, tol) -> bool:
    if e1 is None or e2 is None:
        return e1 is e2
    v1, v2 = vars(e1), vars(e2)
    return v1.keys() == v2.keys() and all(_same(v1[k], v2[k], tol) for k in v1)


def compare_normalized_blocks(b1: PulseqBlock, b2: PulseqBlock, tol: float = 1e-12) -> bool:
    """Return True if two normalized blocks have the same shape and timing.

    RF and ADC offsets are not part of a normalized block, so they are not
    compared; gradient polarity is a signed scale, so it is not either.

    Parameters
    ----------
    b1, b2 : PulseqBlock
        Blocks as :func:`normalize_block` returns them.
    tol : float, default=1e-12
        Absolute tolerance on every number.

    Returns
    -------
    bool
    """
    if abs(b1.block_duration - b2.block_duration) > tol or len(b1.trig) != len(b2.trig):
        return False
    pairs = [(getattr(b1, f), getattr(b2, f)) for f in ("rf", *GRADIENT_AXES, "adc")]
    pairs += list(zip(b1.trig, b2.trig, strict=True))
    return all(_same_event(e1, e2, tol) for e1, e2 in pairs)


def get_grad_energy(g, scale: float, shot: int = 0) -> float:
    """Return the energy of one gradient channel at an instance's scale.

    Parameters
    ----------
    g : SimpleNamespace or None
        A normalized gradient event of a base block.
    scale : float
        The instance's signed scale for this axis (Hz/m).
    shot : int, default=0
        Shot variant of an arbitrary gradient.

    Returns
    -------
    float
        Integral of the squared gradient, in (Hz/m)^2 s; 0 without a gradient.
    """
    if g is None:
        return 0.0
    if g.type == "trap":
        a = g.amplitude
        energy = a * a * (g.rise_time / 3 + g.flat_time + g.fall_time / 3)
    else:
        tt, waveform, _, _ = shot_shape(g, shot)
        w, dt = np.asarray(waveform, dtype=float), np.diff(np.asarray(tt, dtype=float))
        energy = float(np.sum(dt * (w[:-1] ** 2 + w[:-1] * w[1:] + w[1:] ** 2) / 3))
    return float(scale) ** 2 * energy


@dataclass
class CheckReport:
    """What a base-block check found and changed.

    Attributes
    ----------
    ok : bool
        True when no error was found.
    errors, warnings, repairs : list of str
        One sentence each; ``repairs`` lists what hydration changed.
    """

    ok: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    repairs: list[str] = field(default_factory=list)

    def merge(self, other: CheckReport) -> None:
        """Add another report's findings to this one."""
        self.ok &= other.ok
        self.errors += other.errors
        self.warnings += other.warnings
        self.repairs += other.repairs


def _peak(values) -> float:
    values = np.asarray(values)
    return float(np.max(np.abs(values))) if values.size else 0.0


def check_pulseq_block(
    b: PulseqBlock,
    block_name: str = "block",
    *,
    hydrate: bool = False,
    require_normalized: bool = True,
    error_on_fail: bool = True,
    forbid_explicit_delay_blocks: bool = True,
    strip_metadata_extensions: bool = False,
) -> tuple[PulseqBlock, CheckReport]:
    """Check, and optionally hydrate, a base block's Pulseq block.

    Parameters
    ----------
    b : PulseqBlock
        The block to check; hydration edits it in place.
    block_name : str, default="block"
        Name used in messages.
    hydrate : bool, default=False
        Add missing writer fields (event ``type``, gradient ``channel``, ADC
        ``phase_modulation``) and correct a gradient ``channel`` that
        disagrees with its axis.
    require_normalized : bool, default=True
        Require RF peak magnitude 1 and gradient peak amplitude 1 (or 0).
    error_on_fail : bool, default=True
        Raise when an error is found.
    forbid_explicit_delay_blocks : bool, default=True
        Treat a block with no RF, gradient, ADC or trigger event as an error;
        pure delays use the reserved IDs 0 and 1.
    strip_metadata_extensions : bool, default=False
        With ``hydrate``, remove trigger events from the block.

    Returns
    -------
    b : PulseqBlock
    report : CheckReport

    Raises
    ------
    ValueError
        With ``error_on_fail``, listing the errors.
    """
    report = CheckReport()

    def error(message):
        report.errors.append(f"{block_name}: {message}")

    duration = getattr(b, "block_duration", None)
    if duration is None or not np.isfinite(duration) or duration < 0:
        error(f"block_duration is {duration}, expected finite and >= 0")
    if forbid_explicit_delay_blocks and b.is_pure_delay:
        error("is a pure delay; pure delays use the reserved IDs 0 and 1")

    if b.rf is not None:
        if hydrate and getattr(b.rf, "type", None) != "rf":
            b.rf.type = "rf"
            report.repairs.append(f"{block_name}.rf: set type")
        signal, t = getattr(b.rf, "signal", None), getattr(b.rf, "t", None)
        if signal is None or t is None or np.size(signal) != np.size(t):
            error("rf needs signal and t of equal length")
        elif require_normalized and _peak(signal) > 0 and abs(_peak(signal) - 1) > NORMALIZATION_TOLERANCE:
            error(f"rf peak magnitude is {_peak(signal):g}, expected 1")

    for axis in GRADIENT_AXES:
        grad = getattr(b, axis)
        if grad is None:
            continue
        if getattr(grad, "channel", axis[1]) != axis[1]:
            if hydrate:
                report.repairs.append(f"{block_name}.{axis}: channel {grad.channel!r} set to {axis[1]!r}")
            else:
                report.warnings.append(f"{block_name}.{axis}: channel {grad.channel!r}; the axis {axis} wins")
        if hydrate:
            grad.channel = axis[1]
        kind = getattr(grad, "type", None)
        if kind not in ("trap", "grad"):
            error(f"{axis} has type {kind!r}, expected 'trap' or 'grad'")
            continue
        if kind == "trap":
            if any(getattr(grad, f, None) is None for f in _TRAP_FIELDS):
                error(f"{axis} trapezoid needs delay, rise_time, flat_time and fall_time")
            peaks = [abs(getattr(grad, "amplitude", 0.0))]
        else:
            if np.size(getattr(grad, "waveform", [])) != np.size(getattr(grad, "tt", [])):
                error(f"{axis} needs waveform and tt of equal length")
            shots = getattr(grad, "shot_waveforms", [grad.waveform])
            peaks = [_peak(w) for w in shots]
        if require_normalized:
            for peak in peaks:
                if peak > 0 and abs(peak - 1) > NORMALIZATION_TOLERANCE:
                    error(f"{axis} peak amplitude is {peak:g}, expected 1")

    if b.adc is not None:
        if hydrate and getattr(b.adc, "type", None) != "adc":
            b.adc.type = "adc"
            report.repairs.append(f"{block_name}.adc: set type")
        if hydrate and not hasattr(b.adc, "phase_modulation"):
            b.adc.phase_modulation = np.zeros(0)
            report.repairs.append(f"{block_name}.adc: added an empty phase_modulation")
        if getattr(b.adc, "num_samples", 0) < 1 or getattr(b.adc, "dwell", 0) <= 0:
            error("adc needs num_samples >= 1 and dwell > 0")

    if hydrate and strip_metadata_extensions and b.trig:
        b.trig = []
        report.repairs.append(f"{block_name}: removed trigger events")

    report.ok = not report.errors
    if error_on_fail and report.errors:
        raise ValueError("\n".join(report.errors))
    return b, report


def check_base_blocks(rep, **options):
    """Check, and optionally hydrate, every explicit base block of a representation.

    Parameters
    ----------
    rep : Representation
    **options
        Passed to :func:`check_pulseq_block`; ``error_on_fail`` applies to
        the whole representation.

    Returns
    -------
    rep : Representation
    report : CheckReport

    Raises
    ------
    ValueError
        With ``error_on_fail`` (the default), listing every error.
    """
    error_on_fail = options.pop("error_on_fail", True)
    report = CheckReport()
    for base in rep.base_blocks:
        _, found = check_pulseq_block(
            base.block, f"base_blocks[id={base.id}].block", error_on_fail=False, **options
        )
        report.merge(found)
    if error_on_fail and report.errors:
        raise ValueError("\n".join(report.errors))
    return rep, report
