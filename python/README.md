# pulseg (Python)

Python tools for the [PulSeg](https://github.com/HarmonizedMRI/pulseg) intermediate
representation of MRI pulse sequences: base blocks, virtual segments, and the execution
stream of segment instances.

## Install

```bash
pip install pulseg              # data model and validation (numpy only)
pip install "pulseg[pulseq]"    # plus conversion from and to Pulseq, via pypulseqpp
```

From a clone, for development:

```bash
cd python
pip install -e ".[dev]"
pytest
```

## Usage

```python
import pulseg

rep = pulseg.read_seq("gre.seq")  # Pulseq file -> PulSeg representation
pulseg.validate(rep)  # raises ValidationError listing every problem
seq = pulseg.to_sequence(rep)  # back to a pypulseqpp Sequence
seq.write("gre_roundtrip.seq")
```

`pulseg.from_sequence(seq)` converts an in-memory `pypulseqpp.Sequence`.

## Segment boundaries

The source file needs no `TRID` label. A segment instance starts at every repetition
`pypulseqpp.Sequence.repetition()` reports, at every block carrying a physiological
trigger, and wherever the `NOROT`, `NOPOS` or `PMC` label changes. Instances with the
same sequence of base blocks play one virtual segment. A pure delay is base block 1 where
its duration differs between the instances of its virtual segment, and base block 0
otherwise.

## Specification versions

`pulseg.SUPPORTED_VERSIONS` lists the versions the validator accepts: 2.1-alpha and the
2.2-alpha draft. A converted sequence is written as 2.1-alpha unless it uses a field
only 2.2-alpha defines (`rf_shim`, `norot`, `nopos`, `pmc`, `adc_acquire`, `trig_enable`).

## Not yet supported

- Gradient shot variants: every gradient shape is its own base block, so
  `gradient_shot_index` is always 0.
- Soft delays and RF or ADC offsets given in ppm; resolve them in the sequence first.
- A stored file format: representations live in memory.

## Release

A tag `python-vX.Y.Z` on the upstream repository builds the sdist and wheel and
publishes them to PyPI through trusted publishing (`.github/workflows/python-publish.yml`).
Set `__version__` in `src/pulseg/__init__.py` to `X.Y.Z` before tagging.
