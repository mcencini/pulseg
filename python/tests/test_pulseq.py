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


def assert_same_waveforms(a, b, rtol=1e-9):
    wa, wb = a.waveforms_and_times()[0], b.waveforms_and_times()[0]
    assert len(wa) == len(wb)
    for ca, cb in zip(wa, wb, strict=True):
        np.testing.assert_allclose(np.asarray(cb), np.asarray(ca), rtol=rtol, atol=1e-6)


def test_a_converted_sequence_is_a_valid_representation():
    rep = pulseg.import_(gradient_echo())
    pulseg.validate_ir(rep)


def test_conversion_needs_no_trid_label():
    without = pulseg.import_(gradient_echo(trid=False))
    with_trid = pulseg.import_(gradient_echo(trid=True))
    assert len(without.execution_stream) == len(with_trid.execution_stream) == 4


def test_every_repetition_is_one_instance_of_one_virtual_segment():
    rep = pulseg.import_(gradient_echo(n_lines=6))
    assert len(rep.execution_stream) == 6
    assert {inst.virtual_segment_id for inst in rep.execution_stream} == {1}


def test_the_execution_stream_plays_every_block():
    seq = gradient_echo()
    assert pulseg.import_(seq).num_blocks == seq.num_blocks


def test_a_round_trip_plays_the_same_waveforms():
    seq = gradient_echo()
    assert_same_waveforms(seq, pulseg.write_seq(pulseg.import_(seq), system=pp.Opts()))


def test_a_round_trip_keeps_rf_and_adc_phase_offsets():
    seq = gradient_echo()
    back = pulseg.write_seq(pulseg.import_(seq), system=pp.Opts())
    for i in range(1, seq.num_blocks + 1):
        a, b = seq.get_block(i), back.get_block(i)
        if a.rf is not None:
            assert b.rf.phase_offset == pytest.approx(a.rf.phase_offset)
        if a.adc is not None:
            assert b.adc.phase_offset == pytest.approx(a.adc.phase_offset)


def test_a_round_trip_keeps_block_rotations():
    seq = gradient_echo(rotate=True)
    rep = pulseg.import_(seq)
    assert rep.execution_stream[0].rotation_matrix is not None
    back = pulseg.write_seq(rep, system=pp.Opts())
    for i in range(1, seq.num_blocks + 1):
        qa, qb = seq.get_block(i).rotation, back.get_block(i).rotation
        if qa is None or np.allclose(qa.quaternion, [1, 0, 0, 0]):
            continue
        assert abs(np.dot(qa.quaternion, qb.quaternion)) == pytest.approx(1.0)


def test_a_tr_fill_that_varies_is_a_variable_delay():
    rep = pulseg.import_(gradient_echo(tr_fill=lambda line: (5 + line) * 1e-3))
    assert all(
        s.base_block_ids[-1] == pulseg.VARIABLE_DELAY_ID
        for s in rep.virtual_segments
        if len(s.base_block_ids) == 4
    )


def test_a_fixed_tr_fill_is_a_constant_delay():
    rep = pulseg.import_(gradient_echo())
    assert all(s.base_block_ids[-1] == pulseg.CONSTANT_DELAY_ID for s in rep.virtual_segments)


def test_a_physiological_trigger_starts_an_instance_that_waits_for_it():
    seq = gradient_echo(n_lines=2)
    seq.add_block(pp.make_trigger("physio1", duration=2e-3))
    seq.add_block(pp.make_delay(1e-3))
    rep = pulseg.import_(seq)
    assert rep.execution_stream[-1].physio_trigger == 1
    assert all(inst.physio_trigger == 0 for inst in rep.execution_stream[:-1])


def test_a_norot_block_sets_the_instance_flag_and_survives_a_round_trip():
    seq = gradient_echo(n_lines=2)
    seq.add_block(pp.make_label("NOROT", "SET", True), pp.make_trapezoid("z", area=100, duration=1e-3))
    rep = pulseg.import_(seq)
    assert rep.pulseg_version == "2.2-alpha"
    assert rep.execution_stream[-1].norot == 1
    pulseg.validate_ir(rep)
    back = pulseg.write_seq(rep, system=pp.Opts())
    labels = back.evaluate_labels(init={"NOROT": 0}, evolution="blocks")["NOROT"]
    assert labels[-1] == 1 and labels[0] == 0


def test_an_rf_shim_is_carried_per_rf_event_and_survives_a_round_trip():
    seq = pp.Sequence(pp.Opts())
    pulse = pp.make_block_pulse(np.pi / 2, duration=1e-3)
    weights = np.array([1.0, 0.5j])
    seq.add_block(pulse, pp.make_rf_shim(weights))
    seq.add_block(pp.make_delay(1e-3))
    rep = pulseg.import_(seq)
    assert rep.pulseg_version == "2.2-alpha"
    np.testing.assert_allclose(rep.execution_stream[0].rf_shim[0], weights)
    back = pulseg.write_seq(rep, system=pp.Opts())
    assert back.block_shims()[0] == 1


def test_offsets_in_ppm_are_refused():
    seq = pp.Sequence(pp.Opts())
    pulse = pp.make_block_pulse(np.pi / 2, duration=1e-3, freq_ppm=-3.4)
    seq.add_block(pulse)
    with pytest.raises(NotImplementedError):
        pulseg.import_(seq)


def test_import_reads_a_seq_file_and_records_its_path(tmp_path):
    path = tmp_path / "gre.seq"
    gradient_echo().write(str(path))
    rep = pulseg.import_(str(path))
    assert rep.source_file == str(path)
    pulseg.validate_ir(rep)


def test_write_seq_writes_a_file_that_reads_back_as_the_same_scan(tmp_path):
    seq = gradient_echo()
    path = tmp_path / "out.seq"
    pulseg.write_seq(pulseg.import_(seq), path, strict=True, system=pp.Opts())
    back = pp.Sequence()
    back.read(str(path))
    assert_same_waveforms(seq, back, rtol=1e-5)  # a text .seq file stores six significant digits


def spiral_interleaves(n_shots=3, n=40):
    """Interleaves drawn as distinct arbitrary shapes on one timing."""
    seq = pp.Sequence(pp.Opts())
    t = np.linspace(0, 1, n)
    for shot in range(n_shots):
        phase = 2 * np.pi * shot / n_shots
        envelope = np.sin(np.pi * t) * (1 + 0.3 * shot)  # not a rotation of one another
        gx = pp.make_arbitrary_grad("x", 1e5 * envelope * np.cos(6 * t + phase), first=0, last=0)
        gy = pp.make_arbitrary_grad("y", 1e5 * envelope * np.sin(6 * t + phase), first=0, last=0)
        seq.add_block(gx, gy)
        seq.add_block(pp.make_delay(1e-3))
    return seq


def test_interleaves_on_one_timing_are_shot_variants_of_one_base_block():
    rep = pulseg.import_(spiral_interleaves())
    pulseg.validate_ir(rep)
    gradient_blocks = [b for b in rep.base_blocks if b.block.has_gradient]
    assert len(gradient_blocks) == 1
    assert len({i.virtual_segment_id for i in rep.execution_stream}) == 1
    shots = np.concatenate([i.gradient_shot_index for i in rep.execution_stream])
    assert sorted(shots.tolist()) == [0, 1, 2]


def test_shot_variants_play_back_as_the_shapes_they_were_drawn_with():
    seq = spiral_interleaves()
    assert_same_waveforms(seq, pulseg.write_seq(pulseg.import_(seq), system=pp.Opts()))


def test_a_scaled_copy_of_a_shape_is_not_a_new_shot():
    seq = pp.Sequence(pp.Opts())
    shape = np.sin(np.linspace(0, np.pi, 30))
    for scale in (1e5, -5e4):
        seq.add_block(pp.make_arbitrary_grad("x", scale * shape))
    rep = pulseg.import_(seq)
    assert rep.execution_stream[0].gradient_shot_index is None
    assert len(rep.base_blocks) == 1


def test_a_soft_delay_needs_an_input_value():
    seq = gradient_echo(n_lines=1)
    seq.add_block(pp.make_delay(5e-3), pp.make_soft_delay("TE", offset=1e-3, factor=2.0))
    with pytest.raises(ValueError, match="soft_delay_input_ms"):
        pulseg.import_(seq)


def test_a_soft_delay_resolves_to_input_over_factor_plus_offset():
    seq = gradient_echo(n_lines=1)
    seq.add_block(pp.make_delay(5e-3), pp.make_soft_delay("TE", offset=1e-3, factor=2.0))
    rep = pulseg.import_(seq, soft_delay_input_ms=8.0)
    assert rep.execution_stream[-1].block_duration[-1] == pytest.approx(8e-3 / 2 + 1e-3)


def test_stream2loop_has_one_row_per_block_with_its_scales():
    seq = gradient_echo()
    rep = pulseg.import_(seq)
    loop = pulseg.stream2loop(rep)
    assert loop.shape == (seq.num_blocks, 24)
    rf_rows = loop[loop[:, 2] != 0]
    assert len(rf_rows) == 4
    np.testing.assert_allclose(
        loop[:, 13], [seq.get_block(i).block_duration for i in range(1, seq.num_blocks + 1)]
    )
