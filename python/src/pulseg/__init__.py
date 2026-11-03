"""PulSeg: a segment-based intermediate representation of MRI pulse sequences."""

from ._model import (
    CONSTANT_DELAY_ID,
    SUPPORTED_VERSIONS,
    VARIABLE_DELAY_ID,
    BaseBlock,
    PulseqBlock,
    Representation,
    SegmentInstance,
    VirtualSegment,
)
from ._pulseq import from_sequence, physical_blocks, read_seq, to_sequence
from ._validate import ValidationError, find_problems, validate

__version__ = "0.1.0.dev0"

__all__ = [
    "CONSTANT_DELAY_ID",
    "SUPPORTED_VERSIONS",
    "VARIABLE_DELAY_ID",
    "BaseBlock",
    "PulseqBlock",
    "Representation",
    "SegmentInstance",
    "ValidationError",
    "VirtualSegment",
    "find_problems",
    "from_sequence",
    "physical_blocks",
    "read_seq",
    "to_sequence",
    "validate",
]
