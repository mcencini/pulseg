"""Structural validation catches each way a representation can break the specification."""

import numpy as np
import pytest

import pulseg


def test_a_well_formed_representation_has_no_problems(representation):
    assert pulseg.find_problems(representation) == []
    pulseg.validate_ir(representation)


def test_validate_reports_every_problem_at_once(representation):
    representation.pulseg_version = "9.9"
    representation.base_blocks[0].id = 1
    with pytest.raises(pulseg.ValidationError) as error:
        pulseg.validate_ir(representation)
    assert len(error.value.problems) >= 2


def _break_version(rep):
    rep.pulseg_version = "2.0-alpha"


def _break_reserved_id(rep):
    rep.base_blocks[0].id = 0


def _break_duplicate_id(rep):
    rep.base_blocks[1].id = 2


def _break_reference(rep):
    rep.virtual_segments[0].base_block_ids[1] = 99


def _break_unknown_segment(rep):
    rep.execution_stream[0].virtual_segment_id = 7


def _break_rf_length(rep):
    rep.execution_stream[0].rf_amplitude = np.array([1.0, 2.0])


def _break_gradient_shape(rep):
    rep.execution_stream[0].gradient_amplitude = np.array([1e5, 0.0, 0.0])


def _break_negative_duration(rep):
    rep.execution_stream[0].block_duration = np.array([1e-3, -2e-3, 1e-3])


def _break_rf_normalization(rep):
    rep.base_blocks[0].block.rf.signal = rep.base_blocks[0].block.rf.signal * 2


def _break_gradient_normalization(rep):
    rep.base_blocks[1].block.gx.amplitude = 3.0


def _break_explicit_pure_delay(rep):
    rf = rep.base_blocks[0].block
    rf.rf = None


def _break_v22_field_under_v21(rep):
    rep.execution_stream[0].norot = 1


def _break_constant_delay(rep):
    rep.virtual_segments[0].base_block_ids[2] = pulseg.CONSTANT_DELAY_ID


def _break_rotation(rep):
    rep.execution_stream[0].rotation_matrix = np.array([2 * np.eye(3)])


def _break_physio_flag(rep):
    rep.execution_stream[0].physio_trigger = 2


def _break_shot_index(rep):  # the readout carries one shot
    rep.execution_stream[0].gradient_shot_index = np.array([1])


BREAKS = [
    _break_version,
    _break_reserved_id,
    _break_duplicate_id,
    _break_reference,
    _break_unknown_segment,
    _break_rf_length,
    _break_gradient_shape,
    _break_negative_duration,
    _break_rf_normalization,
    _break_gradient_normalization,
    _break_explicit_pure_delay,
    _break_v22_field_under_v21,
    _break_constant_delay,
    _break_rotation,
    _break_physio_flag,
    _break_shot_index,
]


@pytest.mark.parametrize("breaker", BREAKS, ids=[b.__name__.removeprefix("_break_") for b in BREAKS])
def test_each_violation_is_reported(representation, breaker):
    breaker(representation)
    assert pulseg.find_problems(representation)


def test_v22_fields_are_accepted_under_v22(representation):
    representation.pulseg_version = "2.2-alpha"
    inst = representation.execution_stream[0]
    inst.norot = 1
    inst.adc_acquire = np.array([0])
    inst.rf_shim = np.array([[1.0, 1j]])
    assert pulseg.find_problems(representation) == []


def test_num_blocks_counts_every_played_block(representation):
    assert representation.num_blocks == 6
