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

The functions mirror the MATLAB toolbox in `matlab/+pulseg`, with Python naming.
`import` is a Python keyword, so the MATLAB `pulseg.import` is `pulseg.import_`.

```python
import pypulseqpp as pp
import pulseg

rep = pulseg.import_("gre.seq")  # Pulseq file or pypulseqpp.Sequence -> PulSeg
pulseg.validate_ir(rep)  # raises ValidationError listing every problem
seq = pulseg.write_seq(rep, "gre_roundtrip.seq", system=pp.Opts())
loop = pulseg.stream2loop(rep)  # one row per block, the MATLAB column layout
```

| MATLAB | Python |
| --- | --- |
| `pulseg.import` | `pulseg.import_` |
| `pulseg.validate_ir` | `pulseg.validate_ir`, `pulseg.find_problems` |
| `pulseg.check_pulseq_block` | `pulseg.check_pulseq_block` |
| `pulseg.check_base_blocks` | `pulseg.check_base_blocks` |
| `pulseg.normalize_block` | `pulseg.normalize_block` |
| `pulseg.compare_normalized_blocks` | `pulseg.compare_normalized_blocks` |
| `pulseg.get_grad_energy` | `pulseg.get_grad_energy` |
| `pulseg.write_seq` | `pulseg.write_seq` |
| `pulseg.stream2loop` | `pulseg.stream2loop` |

`write_seq` rebuilds arbitrary gradients on the gradient raster of `system`; pass the
system the source sequence was built with.

## Segment boundaries

The source file needs no `TRID` label. A segment instance starts at every repetition
`pypulseqpp.Sequence.repetition()` reports, at every block carrying a physiological
trigger, and wherever the `NOROT`, `NOPOS` or `PMC` label changes. Instances with the
same sequence of base blocks play one virtual segment. A pure delay is base block 1 where
its duration differs between the instances of its virtual segment, and base block 0
otherwise.

## Gradient shot variants

A base block is identified by the shapes of its RF and ADC events and by the timing of
its gradients. Arbitrary gradients that share a timing but differ in shape after
normalization, such as independently designed spiral interleaves, become shot variants
of one base block, and each segment instance selects one through `gradient_shot_index`.
A scaled copy of a waveform is the same shot with a different gradient amplitude.

## Soft delays

Pass the value entered at the console, in ms, as `import_(..., soft_delay_input_ms=...)`.
Every soft delay resolves to `input / factor + offset`.

## Specification versions

`pulseg.SUPPORTED_VERSIONS` lists the versions the validator accepts: 2.1-alpha and the
2.2-alpha draft. A converted sequence is written as 2.1-alpha unless it uses a field
only 2.2-alpha defines (`rf_shim`, `norot`, `nopos`, `pmc`, `adc_acquire`, `trig_enable`).

## Not yet supported

- RF or ADC offsets given in ppm.
- `translateFOVrf`.
- A stored file format: representations live in memory.

## Release

A tag `python-vX.Y.Z` on the upstream repository builds the sdist and wheel and
publishes them to PyPI through trusted publishing (`.github/workflows/python-publish.yml`).
Set `__version__` in `src/pulseg/__init__.py` to `X.Y.Z` before tagging.
