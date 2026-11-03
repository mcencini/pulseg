"""Conversion from and back to Pulseq through pypulseqpp."""

import numpy as np
import pytest

import pulseg

pp = pytest.importorskip("pypulseqpp")


def gradient_echo(n_lines=4, rotate=False, trid=False, tr_fill=None):
    """A small 2D gradient echo: excitation, phase encode, readout, TR fill."""
    seq = pp.Sequence(pp.Opts())
    rf, gz, _ = pp.make_sinc_pulse(np.pi / 6, duration=2e-3, slice_thickness=5e-3, return_gz=True)
    gx = pp.make_trapezoid("x", flat_area=64 / 0.25, flat_time=3e-3)
    adc = pp.make_adc(num_samples=64, duration=3e-3, delay=gx.rise_time)
    gy_max = pp.make_trapezoid("y", area=n_lines * 2, duration=1e-3)
    for line in range(n_lines):
        rf.phase_offset = (117 * line * (line + 1) / 2 % 360) * np.pi / 180
        adc.phase_offset = rf.phase_offset
        excitation = [rf, gz]
        if rotate:
            excitation.append(pp.make_rotation(np.pi * line / n_lines))
        if trid:
            excitation.append(pp.make_label("TRID", "SET", 1))
        seq.add_block(*excitation)
        gy = pp.make_trapezoid(
            "y",
            amplitude=gy_max.amplitude * (line - n_lines / 2 + 0.5) / (n_lines / 2),
            rise_time=gy_max.rise_time,
            flat_time=gy_max.flat_time,
            fall_time=gy_max.fall_time,
        )
        seq.add_block(gy)
        seq.add_block(gx, adc)
        fill = tr_fill(line) if tr_fill else 5e-3
        seq.add_block(pp.make_delay(fill))
    return seq


def assert_same_waveforms(a, b):
    wa, wb = a.waveforms_and_times()[0], b.waveforms_and_times()[0]
    assert len(wa) == len(wb)
    for ca, cb in zip(wa, wb, strict=True):
        np.testing.assert_allclose(np.asarray(cb), np.asarray(ca), rtol=1e-9, atol=1e-6)


def test_a_converted_sequence_is_a_valid_representation():
    rep = pulseg.from_sequence(gradient_echo())
    pulseg.validate(rep)


def test_conversion_needs_no_trid_label():
    without = pulseg.from_sequence(gradient_echo(trid=False))
    with_trid = pulseg.from_sequence(gradient_echo(trid=True))
    assert len(without.execution_stream) == len(with_trid.execution_stream) == 4


def test_every_repetition_is_one_instance_of_one_virtual_segment():
    rep = pulseg.from_sequence(gradient_echo(n_lines=6))
    assert len(rep.execution_stream) == 6
    assert {inst.virtual_segment_id for inst in rep.execution_stream} == {1}


def test_the_execution_stream_plays_every_block():
    seq = gradient_echo()
    assert pulseg.from_sequence(seq).num_blocks == seq.num_blocks


def test_a_round_trip_plays_the_same_waveforms():
    seq = gradient_echo()
    assert_same_waveforms(seq, pulseg.to_sequence(pulseg.from_sequence(seq)))


def test_a_round_trip_keeps_rf_and_adc_phase_offsets():
    seq = gradient_echo()
    back = pulseg.to_sequence(pulseg.from_sequence(seq))
    for i in range(1, seq.num_blocks + 1):
        a, b = seq.get_block(i), back.get_block(i)
        if a.rf is not None:
            assert b.rf.phase_offset == pytest.approx(a.rf.phase_offset)
        if a.adc is not None:
            assert b.adc.phase_offset == pytest.approx(a.adc.phase_offset)


def test_a_round_trip_keeps_block_rotations():
    seq = gradient_echo(rotate=True)
    rep = pulseg.from_sequence(seq)
    assert rep.execution_stream[0].rotation_matrix is not None
    back = pulseg.to_sequence(rep)
    for i in range(1, seq.num_blocks + 1):
        qa, qb = seq.get_block(i).rotation, back.get_block(i).rotation
        if qa is None or np.allclose(qa.quaternion, [1, 0, 0, 0]):
            continue
        assert abs(np.dot(qa.quaternion, qb.quaternion)) == pytest.approx(1.0)


def test_a_tr_fill_that_varies_is_a_variable_delay():
    rep = pulseg.from_sequence(gradient_echo(tr_fill=lambda line: (5 + line) * 1e-3))
    assert all(
        s.base_block_ids[-1] == pulseg.VARIABLE_DELAY_ID
        for s in rep.virtual_segments
        if len(s.base_block_ids) == 4
    )


def test_a_fixed_tr_fill_is_a_constant_delay():
    rep = pulseg.from_sequence(gradient_echo())
    assert all(s.base_block_ids[-1] == pulseg.CONSTANT_DELAY_ID for s in rep.virtual_segments)


def test_a_physiological_trigger_starts_an_instance_that_waits_for_it():
    seq = gradient_echo(n_lines=2)
    seq.add_block(pp.make_trigger("physio1", duration=2e-3))
    seq.add_block(pp.make_delay(1e-3))
    rep = pulseg.from_sequence(seq)
    assert rep.execution_stream[-1].physio_trigger == 1
    assert all(inst.physio_trigger == 0 for inst in rep.execution_stream[:-1])


def test_a_norot_block_sets_the_instance_flag_and_survives_a_round_trip():
    seq = gradient_echo(n_lines=2)
    seq.add_block(pp.make_label("NOROT", "SET", True), pp.make_trapezoid("z", area=100, duration=1e-3))
    rep = pulseg.from_sequence(seq)
    assert rep.pulseg_version == "2.2-alpha"
    assert rep.execution_stream[-1].norot == 1
    pulseg.validate(rep)
    back = pulseg.to_sequence(rep)
    labels = back.evaluate_labels(init={"NOROT": 0}, evolution="blocks")["NOROT"]
    assert labels[-1] == 1 and labels[0] == 0


def test_an_rf_shim_is_carried_per_rf_event_and_survives_a_round_trip():
    seq = pp.Sequence(pp.Opts())
    pulse = pp.make_block_pulse(np.pi / 2, duration=1e-3)
    weights = np.array([1.0, 0.5j])
    seq.add_block(pulse, pp.make_rf_shim(weights))
    seq.add_block(pp.make_delay(1e-3))
    rep = pulseg.from_sequence(seq)
    assert rep.pulseg_version == "2.2-alpha"
    np.testing.assert_allclose(rep.execution_stream[0].rf_shim[0], weights)
    back = pulseg.to_sequence(rep)
    assert back.block_shims()[0] == 1


def test_offsets_in_ppm_are_refused():
    seq = pp.Sequence(pp.Opts())
    pulse = pp.make_block_pulse(np.pi / 2, duration=1e-3, freq_ppm=-3.4)
    seq.add_block(pulse)
    with pytest.raises(NotImplementedError):
        pulseg.from_sequence(seq)


def test_read_seq_records_the_source_file(tmp_path):
    path = tmp_path / "gre.seq"
    gradient_echo().write(str(path))
    rep = pulseg.read_seq(str(path))
    assert rep.source_file == str(path)
    pulseg.validate(rep)
