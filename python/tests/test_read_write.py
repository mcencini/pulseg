"""Reading Pulseq sequences through pulserver and writing them back."""

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
    assert_same_waveforms(seq, pulseg.write(pulseg.read(seq), system=pp.Opts()))


def test_a_round_trip_keeps_rf_and_adc_phase_offsets():
    seq = gradient_echo()
    back = pulseg.write(pulseg.read(seq), system=pp.Opts())
    for i in range(1, seq.num_blocks + 1):
        a, b = seq.get_block(i), back.get_block(i)
        if a.rf is not None:
            assert b.rf.phase_offset == pytest.approx(a.rf.phase_offset, abs=1e-5)
        if a.adc is not None:
            assert b.adc.phase_offset == pytest.approx(a.adc.phase_offset, abs=1e-5)


def test_a_rotated_block_plays_back_as_the_same_physical_gradients():
    seq = gradient_echo(rotate=True)
    assert_same_waveforms(seq, pulseg.write(pulseg.read(seq), system=pp.Opts()))


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
    assert rep.pulseg_version == "2.2-alpha"
    assert rep.execution_stream[-1].norot == 1
    pulseg.validate_ir(rep)
    back = pulseg.write(rep, system=pp.Opts())
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
    pulseg.write(pulseg.read(seq), path, strict=True, system=pp.Opts())
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
    assert_same_waveforms(seq, pulseg.write(pulseg.read(seq), system=pp.Opts()))


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
