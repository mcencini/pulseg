"""Base-block helpers: normalization, comparison, checks and gradient energy."""

from types import SimpleNamespace

import numpy as np
import pypulseqpp as pp
import pytest
from conftest import make_readout_block, make_rf_block

import pulseg


def test_a_well_formed_block_passes_its_check():
    _, report = pulseg.check_pulseq_block(make_rf_block())
    assert report.ok and not report.errors


def test_a_block_failing_its_check_raises_unless_asked_not_to():
    block = make_rf_block()
    block.rf.signal = block.rf.signal * 3
    with pytest.raises(ValueError, match="peak magnitude"):
        pulseg.check_pulseq_block(block)
    _, report = pulseg.check_pulseq_block(block, error_on_fail=False)
    assert not report.ok


def test_hydration_sets_a_gradient_channel_from_its_axis():
    block = make_readout_block()
    block.gx.channel = "z"
    block.adc.__dict__.pop("phase_modulation")
    block, report = pulseg.check_pulseq_block(block, hydrate=True)
    assert block.gx.channel == "x"
    assert hasattr(block.adc, "phase_modulation")
    assert len(report.repairs) == 2


def test_hydration_can_strip_trigger_events():
    block = make_readout_block()
    block.trig = [SimpleNamespace(type="output", channel="osc0", delay=0.0, duration=1e-4)]
    block, _ = pulseg.check_pulseq_block(block, hydrate=True, strip_metadata_extensions=True)
    assert block.trig == []


def test_check_base_blocks_collects_errors_from_every_block(representation):
    for base in representation.base_blocks:
        base.block.block_duration = -1.0
    _, report = pulseg.check_base_blocks(representation, error_on_fail=False)
    assert len(report.errors) == 2


def test_compare_normalized_blocks_ignores_nothing_but_offsets():
    a, b = make_rf_block(), make_rf_block()
    assert pulseg.compare_normalized_blocks(a, b)
    b.rf.signal = np.roll(b.rf.signal, 1)
    assert not pulseg.compare_normalized_blocks(a, b)
    assert not pulseg.compare_normalized_blocks(make_rf_block(1e-3), make_rf_block(2e-3))


def test_trapezoid_energy_matches_its_sampled_integral():
    trap = make_readout_block().gx
    t = np.linspace(0, trap.rise_time + trap.flat_time + trap.fall_time, 200001)
    g = np.interp(t, np.cumsum([0, trap.rise_time, trap.flat_time, trap.fall_time]), [0, 1, 1, 0])
    sampled = np.sum(g[:-1] ** 2 * np.diff(t))
    assert pulseg.get_grad_energy(trap, 2.0) == pytest.approx(4 * sampled, rel=1e-4)


def test_no_gradient_has_no_energy():
    assert pulseg.get_grad_energy(None, 5.0) == 0.0


def test_normalize_block_moves_polarity_into_a_signed_scale():
    seq = pp.Sequence(pp.Opts())
    seq.add_block(pp.make_arbitrary_grad("x", -2e4 * np.sin(np.linspace(0, np.pi, 20))))
    block = seq.get_block(1)
    b0, scales = pulseg.normalize_block(block)
    assert scales["grad"][0] == pytest.approx(-np.max(np.abs(block.gx.waveform)))
    assert np.max(np.abs(b0.gx.waveform)) == pytest.approx(1.0)
