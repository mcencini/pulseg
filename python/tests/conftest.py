"""Shared fixtures: small representations built by hand."""

from types import SimpleNamespace

import numpy as np
import pytest

import pulseg


def make_rf_block(duration=1e-3):
    t = np.arange(10) * 1e-4 + 5e-5
    rf = SimpleNamespace(
        type="rf",
        signal=np.hanning(10).astype(complex) / np.hanning(10).max(),
        t=t,
        shape_dur=1e-3,
        delay=0.0,
        center=5e-4,
        dead_time=0.0,
        ringdown_time=0.0,
        use="excitation",
    )
    return pulseg.PulseqBlock(block_duration=duration, rf=rf)


def make_readout_block(duration=2e-3):
    gx = SimpleNamespace(
        type="trap", channel="x", amplitude=1.0, delay=0.0, rise_time=1e-4, flat_time=1.8e-3, fall_time=1e-4
    )
    adc = SimpleNamespace(
        type="adc", num_samples=64, dwell=2e-5, delay=1e-4, dead_time=0.0, phase_modulation=np.zeros(0)
    )
    return pulseg.PulseqBlock(block_duration=duration, gx=gx, adc=adc)


@pytest.fixture
def representation():
    """Two instances of an excitation, a readout and a variable delay."""
    instances = [
        pulseg.SegmentInstance(
            virtual_segment_id=1,
            rf_amplitude=np.array([250.0]),
            rf_phase_offset=np.array([phase]),
            rf_frequency_offset=np.array([0.0]),
            gradient_amplitude=np.array([[1e5, 0.0, 0.0]]),
            adc_phase_offset=np.array([phase]),
            adc_frequency_offset=np.array([0.0]),
            block_duration=np.array([1e-3, 2e-3, delay]),
        )
        for phase, delay in ((0.0, 1e-3), (np.pi / 2, 3e-3))
    ]
    return pulseg.Representation(
        pulseg_version="2.1-alpha",
        base_blocks=[pulseg.BaseBlock(2, make_rf_block()), pulseg.BaseBlock(3, make_readout_block())],
        virtual_segments=[pulseg.VirtualSegment(1, [2, 3, pulseg.VARIABLE_DELAY_ID])],
        execution_stream=instances,
    )
