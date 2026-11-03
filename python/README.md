# pulseg (Python)

Python tools for the [PulSeg](https://github.com/HarmonizedMRI/pulseg) intermediate
representation of MRI pulse sequences: base blocks, virtual segments, and the execution
stream of segment instances.

## Install

```bash
pip install pulseg
```

pulseg depends on [pulserver](https://pypi.org/project/pulserver/), which segments
Pulseq sequences, and on [pypulseqpp](https://pypi.org/project/pypulseqpp/).

From a clone, for development:

```bash
cd python
pip install -e ".[dev]"
pytest
```

## Usage

```python
import pulseg

rep = pulseg.read("gre.seq")  # Pulseq file or pypulseqpp.Sequence -> PulSeg
pulseg.validate_ir(rep)  # raises ValidationError listing every problem
seq = pulseg.write(rep, "gre_roundtrip.seq")
loop = pulseg.stream2loop(rep)  # one row per block, the MATLAB column layout
```

`read` hands the sequence to pulserver (`pulserver.ir.convert`) and plays the result back
(`pulserver.ir.play`), so pulseg and the pulserver interpreter see the same segmentation.
Each pulserver segment is one virtual segment, and each of its block positions one base
block. `write` plays a representation back as a Pulseq sequence.

Both use the sequence's own system unless given one: `seq.system` for an in-memory
sequence, or the rasters a `.seq` file's `[DEFINITIONS]` record. `read` keeps it as
`rep.system`, which `write` then uses.

The remaining functions mirror the MATLAB toolbox in `matlab/+pulseg`:
`validate_ir`, `check_pulseq_block`, `check_base_blocks`, `normalize_block`,
`compare_normalized_blocks`, `get_grad_energy` and `stream2loop`.

## Segment boundaries

The source file needs no `TRID` label. `pulserver.ir.Grouping` documents where
boundaries fall. A pure delay is base block 1 where its duration differs between the
instances of its virtual segment, and base block 0 otherwise.

## Gradient shot variants

A gradient that a block position plays with another normalized shape or timing in some
instances, such as a spiral interleaf or a rotated readout, is a shot variant of that
position's base block, and each segment instance selects one through
`gradient_shot_index`. A scaled copy of a waveform is the same shot with a different
gradient amplitude. Rotations reach the gradient shapes this way, so `read` leaves
`rotation_matrix` unset.

Amplitudes and shapes come from pulserver in single precision.

## Specification versions

`pulseg.SUPPORTED_VERSIONS` lists the versions the validator accepts: 2.1-alpha and the
2.2-alpha draft. `read` writes 2.1-alpha unless an instance sets `NOROT` or `NOPOS`,
which only 2.2-alpha defines.

## Not yet supported

- RF shims and output triggers in `read`.
- `translateFOVrf`.
- A stored file format: representations live in memory.

## Release

A tag `python-vX.Y.Z` on the upstream repository builds the sdist and wheel and
publishes them to PyPI through trusted publishing (`.github/workflows/python-publish.yml`).
Set `__version__` in `src/pulseg/__init__.py` to `X.Y.Z` before tagging.
