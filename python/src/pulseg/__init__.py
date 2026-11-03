"""PulSeg: a segment-based intermediate representation of MRI pulse sequences.

:func:`read` segments a Pulseq sequence with pulserver; :func:`write` plays a
representation back as a Pulseq sequence. The other functions mirror the
MATLAB toolbox in ``matlab/+pulseg``.
"""

from ._blocks import (
    CheckReport,
    check_base_blocks,
    check_pulseq_block,
    compare_normalized_blocks,
    get_grad_energy,
    normalize_block,
)
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
from ._read import read
from ._validate import ValidationError, find_problems, validate_ir
from ._write import stream2loop, write

__version__ = "0.1.0.dev0"

__all__ = [
    "CONSTANT_DELAY_ID",
    "SUPPORTED_VERSIONS",
    "VARIABLE_DELAY_ID",
    "BaseBlock",
    "CheckReport",
    "PulseqBlock",
    "Representation",
    "SegmentInstance",
    "ValidationError",
    "VirtualSegment",
    "check_base_blocks",
    "check_pulseq_block",
    "compare_normalized_blocks",
    "find_problems",
    "get_grad_energy",
    "normalize_block",
    "read",
    "stream2loop",
    "validate_ir",
    "write",
]
