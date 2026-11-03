"""Reading Pulseq sequences through pulserver and writing them back."""

from types import SimpleNamespace

import numpy as np
import pypulseqpp as pp
import pytest

import pulseg


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


def _sampled(seq):
    """Gradients along x, y and z, as (corner times, values) per axis."""
    waves = seq.waveforms_and_times()[0]
    return [(np.asarray(w[0]), np.asarray(w[1])) for w in waves[:3]]


def assert_same_waveforms(a, b, rtol=1e-5):
    """Same gradients on every axis, compared on the union of both corner grids."""
    for (ta, ga), (tb, gb) in zip(_sampled(a), _sampled(b), strict=True):
        t = np.union1d(ta, tb)
        if t.size == 0:
            continue
        on_a = np.interp(t, ta, ga, left=0, right=0) if ta.size else np.zeros_like(t)
        on_b = np.interp(t, tb, gb, left=0, right=0) if tb.size else np.zeros_like(t)
        scale = max(np.max(np.abs(on_a)), 1.0)
        np.testing.assert_allclose(on_b, on_a, rtol=0, atol=rtol * scale)  # pulserver plays single precision


def rf_instances(rep):
    return [i for i in rep.execution_stream if i.rf_amplitude.size]


def test_a_converted_sequence_is_a_valid_representation():
    rep = pulseg.read(gradient_echo())
    pulseg.validate_ir(rep)


def test_conversion_needs_no_trid_label():
    without = pulseg.read(gradient_echo(trid=False))
    with_trid = pulseg.read(gradient_echo(trid=True))
    assert len(rf_instances(without)) == len(rf_instances(with_trid)) == 4


def test_every_repetition_plays_one_instance_of_one_virtual_segment():
    rep = pulseg.read(gradient_echo(n_lines=6))
    assert len(rf_instances(rep)) == 6
    assert len({inst.virtual_segment_id for inst in rf_instances(rep)}) == 1


def test_the_execution_stream_plays_every_block():
    seq = gradient_echo()
    assert pulseg.read(seq).num_blocks == seq.num_blocks


def test_a_round_trip_plays_the_same_waveforms():
    seq = gradient_echo()
    assert_same_waveforms(seq, pulseg.write(pulseg.read(seq)))


def test_a_round_trip_keeps_rf_and_adc_phase_offsets():
    seq = gradient_echo()
    back = pulseg.write(pulseg.read(seq))
    for i in range(1, seq.num_blocks + 1):
        a, b = seq.get_block(i), back.get_block(i)
        if a.rf is not None:
            assert b.rf.phase_offset == pytest.approx(a.rf.phase_offset, abs=1e-5)
        if a.adc is not None:
            assert b.adc.phase_offset == pytest.approx(a.adc.phase_offset, abs=1e-5)


def test_a_rotated_block_plays_back_as_the_same_physical_gradients():
    seq = gradient_echo(rotate=True)
    assert_same_waveforms(seq, pulseg.write(pulseg.read(seq)))


def _delay_ids(rep):
    return [b for s in rep.virtual_segments for b in s.base_block_ids if b < 2]


def test_a_tr_fill_that_varies_is_a_variable_delay():
    rep = pulseg.read(gradient_echo(tr_fill=lambda line: (5 + line) * 1e-3))
    assert _delay_ids(rep) == [pulseg.VARIABLE_DELAY_ID]


def test_a_fixed_tr_fill_is_a_constant_delay():
    rep = pulseg.read(gradient_echo())
    assert _delay_ids(rep) == [pulseg.CONSTANT_DELAY_ID]


def test_a_physiological_trigger_starts_an_instance_that_waits_for_it():
    seq = gradient_echo(n_lines=2)
    seq.add_block(pp.make_trigger("physio1", duration=2e-3))
    seq.add_block(pp.make_delay(1e-3))
    rep = pulseg.read(seq)
    waiting = [i for i in rep.execution_stream if i.physio_trigger]
    assert len(waiting) == 1
    assert rep.execution_stream.index(waiting[0]) == len(rep.execution_stream) - 2


def test_a_norot_block_sets_the_instance_flag_and_survives_a_round_trip():
    seq = gradient_echo(n_lines=2)
    seq.add_block(pp.make_label("NOROT", "SET", True), pp.make_trapezoid("z", area=100, duration=1e-3))
    rep = pulseg.read(seq)
    assert rep.pulseg_version == "2.2.0"
    assert rep.execution_stream[-1].norot == 1
    pulseg.validate_ir(rep)
    back = pulseg.write(rep)
    labels = back.evaluate_labels(init={"NOROT": 0}, evolution="blocks")["NOROT"]
    assert labels[-1] == 1 and labels[0] == 0


def test_read_records_the_path_of_a_seq_file(tmp_path):
    path = tmp_path / "gre.seq"
    gradient_echo().write(str(path))
    rep = pulseg.read(str(path))
    assert rep.source_file == str(path)
    pulseg.validate_ir(rep)


def test_write_writes_a_file_that_reads_back_as_the_same_scan(tmp_path):
    seq = gradient_echo()
    path = tmp_path / "out.seq"
    pulseg.write(pulseg.read(seq), path, strict=True)
    back = pp.Sequence()
    back.read(str(path))
    assert_same_waveforms(seq, back)


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
    rep = pulseg.read(spiral_interleaves())
    pulseg.validate_ir(rep)
    gradient_blocks = [b for b in rep.base_blocks if b.block.has_gradient]
    assert len(gradient_blocks) == 1
    shots = [int(i.gradient_shot_index[0]) for i in rep.execution_stream if i.gradient_amplitude.size]
    assert shots == [0, 1, 2]


def test_shot_variants_play_back_as_the_shapes_they_were_drawn_with():
    seq = spiral_interleaves()
    assert_same_waveforms(seq, pulseg.write(pulseg.read(seq)))


def test_a_scaled_copy_of_a_shape_is_not_a_new_shot():
    seq = pp.Sequence(pp.Opts())
    shape = np.sin(np.linspace(0, np.pi, 30))
    for scale in (1e5, -5e4):
        seq.add_block(pp.make_arbitrary_grad("x", scale * shape, first=0, last=0))
    rep = pulseg.read(seq)
    assert rep.execution_stream[0].gradient_shot_index is None
    assert len(rep.base_blocks) == 1


def test_stream2loop_has_one_row_per_block_with_its_scales():
    seq = gradient_echo()
    loop = pulseg.stream2loop(pulseg.read(seq))
    assert loop.shape == (seq.num_blocks, 24)
    rf_rows = loop[loop[:, 2] != 0]
    assert len(rf_rows) == 4
    np.testing.assert_allclose(
        loop[:, 13], [seq.get_block(i).block_duration for i in range(1, seq.num_blocks + 1)]
    )


def test_read_and_write_follow_the_rasters_a_seq_file_records(tmp_path):
    system = pp.Opts(grad_raster_time=1e-5, block_duration_raster=1e-5, rf_raster_time=1e-6)
    seq = pp.Sequence(system)
    shape = np.sin(np.linspace(0, np.pi, 30))
    seq.add_block(pp.make_arbitrary_grad("x", 1e5 * shape, first=0, last=0, system=system))
    path = tmp_path / "fine.seq"
    seq.write(str(path))
    rep = pulseg.read(path)
    assert rep.system.grad_raster_time == pytest.approx(1e-5)
    source, back = seq.get_block(1).gx, pulseg.write(rep).get_block(1).gx
    on_back = np.interp(source.tt, back.tt, back.waveform)
    np.testing.assert_allclose(on_back, source.waveform, rtol=0, atol=1e-5 * 1e5)


def test_write_defaults_to_the_system_the_sequence_was_read_under():
    seq = spiral_interleaves()
    rep = pulseg.read(seq)
    assert rep.system is seq.system


def test_an_rf_offset_in_ppm_reads_as_hz_at_the_systems_field():
    system = pp.Opts(B0=7.0)
    seq = pp.Sequence(system)
    seq.add_block(pp.make_block_pulse(np.pi / 2, duration=1e-3, freq_ppm=-3.4, system=system))
    seq.add_block(pp.make_delay(1e-3))
    rep = pulseg.read(seq)
    expected = -3.4e-6 * system.gamma * 7.0
    assert rep.execution_stream[0].rf_frequency_offset[0] == pytest.approx(expected, rel=1e-5)


def test_a_soft_delay_reads_as_its_designed_duration():
    seq = gradient_echo(n_lines=1)
    seq.add_block(pp.make_delay(5e-3), pp.make_soft_delay("TE", offset=1e-3, factor=2.0))
    rep = pulseg.read(seq)
    assert rep.execution_stream[-1].block_duration[-1] == pytest.approx(5e-3)


def test_a_disabled_output_trigger_is_not_written(representation):
    trigger = SimpleNamespace(type="output", channel="osc0", delay=0.0, duration=1e-4)
    representation.base_block(2).block.trig = [trigger]
    representation.pulseg_version = "2.2.0"
    representation.execution_stream[1].trig_enable = np.array([0])
    pulseg.validate_ir(representation)
    back = pulseg.write(representation, system=pp.Opts())
    fired = [i for i in range(1, back.num_blocks + 1) if getattr(back.get_block(i), "trig", None)]
    assert fired == [1]


def test_a_readout_played_but_not_stored_is_refused_by_write(representation):
    representation.pulseg_version = "2.2.0"
    representation.execution_stream[0].adc_acquire = np.array([0])
    pulseg.validate_ir(representation)
    with pytest.raises(NotImplementedError, match="adc_acquire"):
        pulseg.write(representation, system=pp.Opts())


def test_an_rf_shim_is_read_per_rf_event_and_survives_a_round_trip():
    seq = pp.Sequence(pp.Opts())
    weights = np.array([1.0, 0.5j])
    seq.add_block(pp.make_block_pulse(np.pi / 2, duration=1e-3), pp.make_rf_shim(weights))
    seq.add_block(pp.make_delay(1e-3))
    rep = pulseg.read(seq)
    rf_instance = rf_instances(rep)[0]
    np.testing.assert_allclose(rf_instance.rf_shim[0], weights, atol=1e-5)
    back = pulseg.write(rep)
    np.testing.assert_allclose(back.get_block(1).rf_shim.shim_vector, weights, atol=1e-5)


def test_an_output_trigger_is_read_into_its_base_block_and_written_back():
    seq = gradient_echo(n_lines=2)
    seq.add_block(pp.make_delay(2e-3), pp.make_digital_output_pulse("ext1", delay=20e-6, duration=100e-6))
    rep = pulseg.read(seq)
    triggered = [b for b in rep.base_blocks if b.block.trig]
    assert len(triggered) == 1
    assert triggered[0].block.trig[0].channel == "ext1"
    back = pulseg.write(rep)
    fired = [i for i in range(1, back.num_blocks + 1) if getattr(back.get_block(i), "trig", None)]
    assert fired == [back.num_blocks]
