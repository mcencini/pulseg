"""PulSeg data model: base blocks, virtual segments, segment instances."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace

import numpy as np

SUPPORTED_VERSIONS = ("2.1-alpha", "2.2-alpha")
"""Specification versions this package reads and validates."""

CONSTANT_DELAY_ID = 0
"""Reserved base-block ID of a pure delay whose duration is fixed at its position."""

VARIABLE_DELAY_ID = 1
"""Reserved base-block ID of a pure delay whose duration varies between instances."""

GRADIENT_AXES = ("gx", "gy", "gz")


@dataclass
class PulseqBlock:
    """A Pulseq block at normalized amplitude.

    Events are :class:`types.SimpleNamespace` objects in the PyPulseq field
    vocabulary (``signal``, ``t``, ``amplitude``, ``rise_time``, ...), so they
    can be handed to ``pypulseqpp.Sequence.add_block`` after the instance
    parameters are applied.

    Attributes
    ----------
    block_duration : float
        Block duration in seconds.
    rf, gx, gy, gz, adc : SimpleNamespace or None
        At most one event of each kind.
    trig : list of SimpleNamespace
        Trigger events: output triggers (``type == "output"``) and
        physiological triggers (``type == "trigger"``).
    """

    block_duration: float
    rf: SimpleNamespace | None = None
    gx: SimpleNamespace | None = None
    gy: SimpleNamespace | None = None
    gz: SimpleNamespace | None = None
    adc: SimpleNamespace | None = None
    trig: list = field(default_factory=list)

    @property
    def has_gradient(self) -> bool:
        """True if the block plays a gradient on any axis."""
        return any(getattr(self, axis) is not None for axis in GRADIENT_AXES)

    @property
    def is_pure_delay(self) -> bool:
        """True if the block plays no RF, gradient, ADC or trigger event."""
        return self.rf is None and not self.has_gradient and self.adc is None and not self.trig


@dataclass
class BaseBlock:
    """A Pulseq block with RF and gradient amplitudes normalized to one.

    Attributes
    ----------
    id : int
        Unique identifier, at least 2.
    block : PulseqBlock
        The normalized block.
    name : str or None
        Optional descriptive label.
    """

    id: int
    block: PulseqBlock
    name: str | None = None


@dataclass
class VirtualSegment:
    """An ordered list of base blocks forming one reusable unit of the sequence.

    Attributes
    ----------
    id : int
        Unique positive identifier.
    base_block_ids : list of int
        Base-block IDs in play order; 0 and 1 are the reserved pure delays.
    name : str or None
        Optional descriptive label.
    """

    id: int
    base_block_ids: list[int]
    name: str | None = None


@dataclass
class SegmentInstance:
    """One occurrence of a virtual segment in the execution stream.

    Per-event arrays are indexed by the RF, gradient and ADC events of the
    virtual segment in play order. A gradient event is one block with a
    gradient on any axis, so it consumes one row of ``gradient_amplitude``.

    Attributes
    ----------
    virtual_segment_id : int
        The virtual segment this instance plays.
    rf_amplitude, rf_phase_offset, rf_frequency_offset : numpy.ndarray
        ``(n_rf,)``: scale factor (Hz), phase (rad) and frequency (Hz).
    gradient_amplitude : numpy.ndarray
        ``(n_grad, 3)``: signed scale factors on x, y and z (Hz/m).
    adc_phase_offset, adc_frequency_offset : numpy.ndarray
        ``(n_adc,)``: receiver phase (rad) and frequency (Hz).
    block_duration : numpy.ndarray
        ``(n_blocks,)``: duration of each block in seconds.
    gradient_shot_index : numpy.ndarray or None
        ``(n_grad,)`` selected shot variant; None means 0 everywhere.
    rotation_matrix : numpy.ndarray or None
        ``(n_grad, 3, 3)`` rotation applied after scaling; None means identity.
    rf_shim : numpy.ndarray or None
        ``(n_rf, n_channels)`` complex weight per transmit channel (2.2-alpha);
        None means unit weights.
    adc_acquire : numpy.ndarray or None
        ``(n_adc,)``, 0 where the readout is played but not stored (2.2-alpha);
        None means 1 everywhere.
    trig_enable : numpy.ndarray or None
        ``(n_trig_blocks,)``, 0 where the block's output triggers do not fire
        (2.2-alpha); None means 1 everywhere.
    norot, nopos, pmc : int
        Pulseq ``NOROT``, ``NOPOS`` and ``PMC`` flags of the instance (2.2-alpha).
    physio_trigger : int
        1 if the instance waits for a physiological trigger.
    label : str or None
        Optional label.
    """

    virtual_segment_id: int
    rf_amplitude: np.ndarray
    rf_phase_offset: np.ndarray
    rf_frequency_offset: np.ndarray
    gradient_amplitude: np.ndarray
    adc_phase_offset: np.ndarray
    adc_frequency_offset: np.ndarray
    block_duration: np.ndarray
    gradient_shot_index: np.ndarray | None = None
    rotation_matrix: np.ndarray | None = None
    rf_shim: np.ndarray | None = None
    adc_acquire: np.ndarray | None = None
    trig_enable: np.ndarray | None = None
    norot: int = 0
    nopos: int = 0
    pmc: int = 0
    physio_trigger: int = 0
    label: str | None = None


@dataclass
class Representation:
    """A complete PulSeg representation.

    Attributes
    ----------
    pulseg_version : str
        Specification version, one of :data:`SUPPORTED_VERSIONS`.
    base_blocks : list of BaseBlock
        Explicit base blocks (IDs at least 2).
    virtual_segments : list of VirtualSegment
    execution_stream : list of SegmentInstance
    source_file : str or None
    creation_date : str or None
        ISO 8601 date.
    system : pypulseqpp.Opts or None
        Limits and rasters the source sequence was built with; not part of
        the specification.
    """

    pulseg_version: str
    base_blocks: list[BaseBlock]
    virtual_segments: list[VirtualSegment]
    execution_stream: list[SegmentInstance]
    source_file: str | None = None
    creation_date: str | None = None
    system: object | None = None

    def base_block(self, block_id: int) -> BaseBlock:
        """Return the explicit base block with this ID.

        Raises
        ------
        KeyError
            If no explicit base block has this ID.
        """
        for base in self.base_blocks:
            if base.id == block_id:
                return base
        raise KeyError(f"no base block with id {block_id}")

    def virtual_segment(self, segment_id: int) -> VirtualSegment:
        """Return the virtual segment with this ID.

        Raises
        ------
        KeyError
            If no virtual segment has this ID.
        """
        for segment in self.virtual_segments:
            if segment.id == segment_id:
                return segment
        raise KeyError(f"no virtual segment with id {segment_id}")

    @property
    def num_blocks(self) -> int:
        """Number of Pulseq blocks the execution stream plays."""
        lengths = {s.id: len(s.base_block_ids) for s in self.virtual_segments}
        return sum(lengths[i.virtual_segment_id] for i in self.execution_stream)
