"""Structural and base-block validation of a PulSeg representation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._blocks import check_pulseq_block, num_shots
from ._model import (
    CONSTANT_DELAY_ID,
    SUPPORTED_VERSIONS,
    VARIABLE_DELAY_ID,
    PulseqBlock,
    Representation,
    SegmentInstance,
)

V22_FIELDS = ("rf_shim", "adc_acquire", "trig_enable")
V22_FLAGS = ("norot", "nopos", "pmc")


class ValidationError(ValueError):
    """A representation that breaks the specification.

    Attributes
    ----------
    problems : list of str
        Every problem found, one sentence each.
    """

    def __init__(self, problems: list[str]):
        self.problems = problems
        super().__init__(
            f"{len(problems)} problem(s) in the PulSeg representation:\n- " + "\n- ".join(problems)
        )


@dataclass(frozen=True)
class EventCounts:
    """Events of one virtual segment that the per-instance arrays index."""

    blocks: int
    rf: int
    grad: int
    adc: int
    trig: int
    shots: tuple[int, ...] = ()


def has_output_trigger(block: PulseqBlock) -> bool:
    """True if the block carries at least one output trigger."""
    return any(getattr(t, "type", None) == "output" for t in block.trig)


def count_events(rep: Representation, base_block_ids: list[int]) -> EventCounts:
    """Count the blocks and events of an ordered list of base blocks.

    Parameters
    ----------
    rep : Representation
        Holds the base blocks the IDs refer to.
    base_block_ids : list of int
        Base-block IDs, the reserved delays included.

    Returns
    -------
    EventCounts
    """
    rf = grad = adc = trig = 0
    shots = []
    for block_id in base_block_ids:
        if block_id in (CONSTANT_DELAY_ID, VARIABLE_DELAY_ID):
            continue
        block = rep.base_block(block_id).block
        rf += block.rf is not None
        grad += block.has_gradient
        if block.has_gradient:
            shots.append(num_shots(block))
        adc += block.adc is not None
        trig += has_output_trigger(block)
    return EventCounts(len(base_block_ids), rf, grad, adc, trig, tuple(shots))


def _check_array(problems, name, value, shape, finite=True):
    array = np.asarray(value)
    if array.shape != shape:
        problems.append(f"{name} has shape {array.shape}, expected {shape}.")
        return None
    if finite and array.size and not np.all(np.isfinite(array)):
        problems.append(f"{name} contains non-finite values.")
    return array


def _check_instance(problems, name, inst: SegmentInstance, counts: EventCounts, version: str):
    _check_array(problems, f"{name}.rf_amplitude", inst.rf_amplitude, (counts.rf,))
    _check_array(problems, f"{name}.rf_phase_offset", inst.rf_phase_offset, (counts.rf,))
    _check_array(problems, f"{name}.rf_frequency_offset", inst.rf_frequency_offset, (counts.rf,))
    _check_array(problems, f"{name}.gradient_amplitude", inst.gradient_amplitude, (counts.grad, 3))
    _check_array(problems, f"{name}.adc_phase_offset", inst.adc_phase_offset, (counts.adc,))
    _check_array(problems, f"{name}.adc_frequency_offset", inst.adc_frequency_offset, (counts.adc,))
    durations = _check_array(problems, f"{name}.block_duration", inst.block_duration, (counts.blocks,))
    if durations is not None and durations.size and np.any(durations < 0):
        problems.append(f"{name}.block_duration has negative entries.")
    if inst.gradient_shot_index is not None:
        shots = _check_array(
            problems, f"{name}.gradient_shot_index", inst.gradient_shot_index, (counts.grad,)
        )
        if (
            shots is not None
            and shots.size
            and (np.any(shots < 0) or np.any(shots >= np.array(counts.shots)))
        ):
            problems.append(
                f"{name}.gradient_shot_index selects a shot variant its base block does not carry."
            )
    if inst.rotation_matrix is not None:
        rotations = _check_array(
            problems, f"{name}.rotation_matrix", inst.rotation_matrix, (counts.grad, 3, 3)
        )
        if rotations is not None and rotations.size:
            gram = np.einsum("nij,nkj->nik", rotations, rotations)
            if not np.allclose(gram, np.eye(3), atol=1e-6):
                problems.append(f"{name}.rotation_matrix holds a matrix that is not a rotation.")
    if inst.physio_trigger not in (0, 1):
        problems.append(f"{name}.physio_trigger is {inst.physio_trigger}, expected 0 or 1.")
    if inst.label is not None and not isinstance(inst.label, str):
        problems.append(f"{name}.label must be a string.")

    used_v22 = [f for f in V22_FIELDS if getattr(inst, f) is not None]
    used_v22 += [f for f in V22_FLAGS if getattr(inst, f) != 0]
    if used_v22 and version != "2.2.0":
        problems.append(f"{name} sets {', '.join(used_v22)}, which require pulseg_version 2.2.0.")
    if inst.rf_shim is not None:
        shim = np.asarray(inst.rf_shim)
        if shim.ndim != 2 or shim.shape[0] != counts.rf:
            problems.append(f"{name}.rf_shim has shape {shim.shape}, expected ({counts.rf}, n_channels).")
    for flag_field, count in (("adc_acquire", counts.adc), ("trig_enable", counts.trig)):
        value = getattr(inst, flag_field)
        if value is not None:
            flags = _check_array(problems, f"{name}.{flag_field}", value, (count,))
            if flags is not None and not np.all(np.isin(flags, (0, 1))):
                problems.append(f"{name}.{flag_field} entries must be 0 or 1.")
    for flag in V22_FLAGS:
        if getattr(inst, flag) not in (0, 1):
            problems.append(f"{name}.{flag} is {getattr(inst, flag)}, expected 0 or 1.")


def find_problems(rep: Representation) -> list[str]:
    """List every way a representation breaks the specification.

    Checks the version, the uniqueness and ranges of IDs, the references
    between base blocks, virtual segments and instances, the length of every
    per-instance array against the events of its virtual segment, the
    normalization of every base block, and that a pure delay keeps one
    duration across the instances of a constant-delay position.

    Parameters
    ----------
    rep : Representation

    Returns
    -------
    list of str
        Empty for a valid representation.
    """
    problems: list[str] = []
    if rep.pulseg_version not in SUPPORTED_VERSIONS:
        problems.append(f"pulseg_version is {rep.pulseg_version!r}, expected one of {SUPPORTED_VERSIONS}.")
    for field_name in ("base_blocks", "virtual_segments", "execution_stream"):
        if not getattr(rep, field_name):
            problems.append(f"{field_name} is empty.")

    base_ids = [b.id for b in rep.base_blocks]
    if len(set(base_ids)) != len(base_ids):
        problems.append("base block IDs are not unique.")
    for base in rep.base_blocks:
        if base.id < 2:
            problems.append(f"base block ID {base.id} is below 2; 0 and 1 are reserved pure delays.")
        _, report = check_pulseq_block(base.block, f"base_blocks[id={base.id}].block", error_on_fail=False)
        problems += report.errors

    segment_ids = [s.id for s in rep.virtual_segments]
    if len(set(segment_ids)) != len(segment_ids):
        problems.append("virtual segment IDs are not unique.")
    known_blocks = set(base_ids) | {CONSTANT_DELAY_ID, VARIABLE_DELAY_ID}
    counts = {}
    for segment in rep.virtual_segments:
        name = f"virtual_segments[id={segment.id}]"
        if segment.id < 1:
            problems.append(f"{name} has a non-positive ID.")
        if not segment.base_block_ids:
            problems.append(f"{name} has no base blocks.")
        unknown = sorted(set(segment.base_block_ids) - known_blocks)
        if unknown:
            problems.append(f"{name} references unknown base blocks {unknown}.")
            continue
        counts[segment.id] = count_events(rep, segment.base_block_ids)

    constant_durations: dict[tuple[int, int], float] = {}
    for n, inst in enumerate(rep.execution_stream):
        name = f"execution_stream[{n}]"
        if inst.virtual_segment_id not in counts:
            problems.append(f"{name} references unknown virtual segment {inst.virtual_segment_id}.")
            continue
        _check_instance(problems, name, inst, counts[inst.virtual_segment_id], rep.pulseg_version)
        segment = rep.virtual_segment(inst.virtual_segment_id)
        durations = np.asarray(inst.block_duration)
        if durations.shape != (len(segment.base_block_ids),):
            continue
        for position, block_id in enumerate(segment.base_block_ids):
            if block_id != CONSTANT_DELAY_ID:
                continue
            key = (segment.id, position)
            first = constant_durations.setdefault(key, float(durations[position]))
            if abs(first - durations[position]) > 1e-12:
                problems.append(
                    f"{name} plays constant delay {position} of virtual segment {segment.id} for "
                    f"{durations[position]:g} s, where an earlier instance plays {first:g} s."
                )
    return problems


def validate_ir(rep: Representation) -> None:
    """Raise if a representation breaks the specification.

    Parameters
    ----------
    rep : Representation

    Raises
    ------
    ValidationError
        Listing every problem :func:`find_problems` finds.
    """
    problems = find_problems(rep)
    if problems:
        raise ValidationError(problems)
