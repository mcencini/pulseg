![pulseg logo](../assets/logo.svg)

# PulSeg Intermediate Representation Specification

**Version:** 2.2.0  
**Date:** 2026-11-04  
**Status:** Draft proposal  
**Authors:** Jon-Fredrik Nielsen, Matteo Cencini  
**Repository:** https://github.com/HarmonizedMRI/pulseg

---

## 1. Overview

PulSeg is a vendor-neutral intermediate representation (IR) for MRI pulse sequences. It is
designed to sit between a high-level sequence description format, such as
[Pulseq](https://pulseq.github.io/), and a vendor-specific hardware execution layer.

By making explicit the repeating structural units present in most MRI sequences, PulSeg
enables efficient and unambiguous mapping onto the internal execution models of scanner
hardware, simulation frameworks, and sequence analysis tools.

This document defines the PulSeg data structures, terminology, and conversion requirements.
It is intended for developers of sequence conversion tools, scanner interpreters, and
simulation frameworks.

This specification is stable at version `2.2.0`. Any change affecting data structure
definitions, field semantics, required fields, or normalization rules must increment the
version number and include migration notes in the changelog.

Implementations MAY include additional metadata fields not defined in this specification.
Compliant readers SHOULD ignore unknown fields unless they conflict with required PulSeg
semantics. Such fields are non-normative and MUST NOT be required for interpreting the core
IR.

---

## 2. Definitions

- **Pulseq block:**  
  An atomic time unit in the [Pulseq](https://pulseq.github.io/) sequence format,
  containing at most one waveform per gradient axis, one RF waveform, and one ADC window.
  Waveform amplitudes in a physical Pulseq block reflect the actual amplitudes used during
  that block.

- **PulSeg PulseqBlock:**  
  A Pulseq-compatible block representation used inside a PulSeg `BaseBlock`. A PulSeg
  `PulseqBlock` stores waveform and timing shape information. It is not required to preserve
  every implementation-specific field used internally by a particular Pulseq library.

- **Base block:**  
  A PulSeg `PulseqBlock` in which all RF and gradient waveform amplitudes have been
  normalized to a canonical value. A base block defines waveform/timing *shape*, not the
  physical amplitudes, phases, frequency offsets, or spatial rotations used in any particular
  execution. Base blocks are the atomic elements of a PulSeg representation.

- **Virtual segment:**  
  An ordered, finite sequence of base block IDs representing a reusable unit of the MRI
  sequence — for example, a readout module, a preparation module, or a spoiler module (it
  need not be periodic). A virtual segment defines sequence *structure*, not the specific amplitudes,
  phases, frequency offsets, or durations used in any particular execution.

- **Segment instance:**  
  A concrete realization of a virtual segment within the execution stream. A segment instance
  associates a virtual segment with physical RF amplitude scales, RF phase offsets, RF
  frequency offsets, gradient amplitude scales, ADC receiver offsets, block durations, and
  optional execution metadata.

- **Execution stream:**  
  The ordered sequence of segment instances that defines the complete execution of the MRI
  sequence.

- **Hydration:**  
  The process of adding or correcting implementation-specific Pulseq event fields that are
  redundant with PulSeg semantics and required by a particular Pulseq software library. For
  example, a MATLAB Pulseq writer may require an event `type`, gradient `channel`, or ADC
  `phaseModulation` field. Hydration MUST NOT alter waveform shape, physical timing,
  physical amplitude, phase, frequency offset, event ordering, or other PulSeg execution
  semantics.

---

## 3. Data Structures

All field names use `snake_case`. Fields marked **required** must be present in any compliant
PulSeg representation. Fields marked **optional** may be omitted.

Unless otherwise stated, durations are expressed in seconds, phase offsets in radians,
frequency offsets in Hz, and gradient/RF amplitude scaling factors in the native units implied
by the source sequence or target backend. A stored representation MAY use other units and
numeric formats if it declares them (Section 5.5).

### 3.0 PulSeg PulseqBlock

Within PulSeg, a `PulseqBlock` is a Pulseq-compatible block representation containing
waveform and timing shape information.

A PulSeg `PulseqBlock` MUST be a block-like structure containing:

| Field | Type | Required | Description |
|---|---:|---:|---|
| `blockDuration` | float | required | Finite non-negative duration of the block in seconds. |

A PulSeg `PulseqBlock` MAY contain at most one of each event field:

| Field | Description |
|---|---|
| `rf` | RF transmit event. |
| `gx` | Gradient event on the physical x axis. |
| `gy` | Gradient event on the physical y axis. |
| `gz` | Gradient event on the physical z axis. |
| `adc` | ADC receive event. |
| `trig` | Output trigger event: channel, delay and duration. |

Explicit PulSeg base blocks with IDs ≥ 2 MUST contain at least one RF, gradient, ADC or
output trigger event. Pure delays MUST be represented using the reserved implicit base block IDs `0` and
`1`, rather than explicit `BaseBlock` entries.

#### 3.0.1 Gradient axis convention

For PulSeg semantics, the gradient axis is determined by the containing block field name:

- `gx` denotes the physical x-gradient channel
- `gy` denotes the physical y-gradient channel
- `gz` denotes the physical z-gradient channel

If an implementation stores a nested gradient event field such as `channel`, that field MUST
be consistent with the containing block field. If inconsistent, the containing block field is
authoritative.

For example, a gradient event stored as `block.gz` is interpreted as a z-gradient even if an
implementation-specific nested field says `block.gz.channel = "x"`. A compliant
implementation MAY hydrate this event by setting `block.gz.channel = "z"`.

#### 3.0.2 Event-level shape fields

The exact internal representation of RF, gradient, and ADC events may vary between Pulseq
software libraries. However, a PulSeg `PulseqBlock` MUST contain enough information to
recover the waveform/timing shape of each event.

At minimum:

- RF events SHOULD contain RF signal samples and sample timing information.
- Arbitrary/sample-based gradient events SHOULD contain gradient waveform samples and sample
  timing information.
- Trapezoid/scalar gradient events SHOULD contain amplitude and timing fields sufficient to
  reconstruct the trapezoid shape.
- ADC events SHOULD contain sample count, dwell time, and delay.

Implementation-specific fields such as event `type`, gradient `channel`, ADC
`phaseModulation`, or Pulseq writer-specific defaults MAY be added during hydration. Such
fields are not part of PulSeg semantics unless explicitly specified here.

#### 3.0.3 Base-block metadata and Pulseq extensions

A `BaseBlock.block` defines reusable waveform and timing shape. Pulseq labels, rotation
extensions, comments, and vendor-specific extensions SHOULD NOT be stored as part of
normalized base-block shape unless explicitly required by a future version of this
specification.

An output trigger is the exception: its delay and duration are part of the block's timing,
which a backend prepares with the segment, so it is stored as the `trig` event. Whether it
fires in a given instance is set by `trig_enable` (Section 3.3).

Execution-dependent metadata should be represented through `SegmentInstance` fields such as
`rotation_matrix`, `norot`, `nopos`, `pmc`, `physio_trigger`, or `label`.

---

### 3.1 BaseBlock

A base block wraps a single normalized PulSeg `PulseqBlock`. A base block is identified by
the normalized *shapes* of its RF and ADC events and by the *timing* of its gradient events.
Because gradient identity is timing-based, a single gradient event may carry more than one
co-timed waveform *shape* — a *shot* dimension — much as a multi-channel (pTx) RF event
carries more than one transmit channel (see the normalization rules).

| Field | Type | Required | Description |
|---|---:|---:|---|
| `id` | int | required | Unique identifier for this base block. Must be a non-negative integer. IDs `0` and `1` are globally reserved. User-defined explicit base blocks must use IDs ≥ 2. |
| `block` | PulseqBlock | required | A PulSeg `PulseqBlock` with RF and gradient waveform amplitudes normalized according to the rules below. |
| `name` | string | optional | Human-readable descriptive label, e.g. `"rf_prep"`, `"readout"`, `"gx_spoil"`. |

A `BaseBlock` entry with `id >= 2` MUST NOT represent a pure delay block. Explicit base
blocks MUST contain at least one non-empty `rf`, `gx`, `gy`, `gz`, or `adc` event.

`base_blocks` contains only explicit base blocks with IDs ≥ 2. Reserved delay blocks with IDs
`0` and `1` are implicit and MUST NOT appear as entries in `base_blocks`.

For PulSeg `2.2.0`, `base_blocks` MUST be non-empty.

#### Reserved identifiers

To support real-time timing optimizations on scanner hardware, the following base block IDs
are globally reserved and implicitly defined across all PulSeg implementations:

- `id == 0` — **Implicit Constant Delay Block**  
  Represents a pure delay block whose duration remains fixed and invariant across all
  instances of the corresponding virtual-segment position in the execution stream. This
  allows hardware interpreters to pre-compute structural timing gaps.

- `id == 1` — **Implicit Variable Delay Block**  
  Represents a pure delay block whose duration may vary dynamically between different
  instances of the corresponding virtual-segment position, for example to extend TE or TR.
  The runtime duration is dictated explicitly by the corresponding entry in the segment
  instance's `block_duration` array.

Whether a pure delay is constant or variable follows from the execution stream: it is
variable if its `block_duration` differs between instances of its virtual segment. A stored
representation MAY therefore omit the reserved IDs and let the reader assign them on loading
(Section 5.5).

#### Normalization rules

- **Single-channel RF waveforms:**  
  Normalize by peak magnitude such that:

  ```text
  max(abs(rf.signal)) == 1.0
  ```

  for nonzero RF waveforms.

- **Multi-channel RF waveforms / pTx:**  
  Normalize all transmit channels by a single global scaling factor derived from the maximum
  peak magnitude across all channels combined:

  ```text
  max_c(max(abs(rf.signal_c))) == 1.0
  ```

  This preserves relative amplitude between physical transmit channels.

- **Gradient waveforms:**  
  Normalize each physical gradient channel independently by peak absolute amplitude such that:

  ```text
  max(abs(gx.waveform)) == 1.0
  max(abs(gy.waveform)) == 1.0
  max(abs(gz.waveform)) == 1.0
  ```

  for nonzero arbitrary/sample-based gradient waveforms.

  For trapezoid/scalar gradients, the scalar amplitude is normalized such that:

  ```text
  abs(grad.amplitude) == 1.0
  ```

  for nonzero trapezoid amplitudes.

  A gradient event MAY carry multiple co-timed *shot variants* — gradient shapes sharing one
  timing structure but differing in sample values, e.g. independently optimized spiral
  interleaves — and each shot variant is normalized independently (peak `1.0` per variant).
  Two gradient waveforms share a *timing structure* when they have identical
  `(delay, rise, flat, fall)` for trapezoids, `(delay, time_shape_id)` for extended
  trapezoids, or `(delay, num_samples)` on the gradient raster for uniform-raster arbitrary
  gradients. Shot variants are reserved for shapes NOT related by scaling or rotation;
  families related by scaling or rotation (phase-encode blips, radial spokes, rotated
  interleaves) MUST use a single shot variant scaled via `gradient_amplitude` /
  `rotation_matrix`. The single-shot case is the common one.

- **ADC windows:**  
  ADC sample count, dwell time, delay, and sampling geometry are not amplitude-normalized and
  are copied directly into the base block.

- **Output triggers:**  
  Channel, delay and duration are not normalized and are copied directly into the base block.

- **Unused channels:**  
  A channel with no waveform or event in the source block MUST have no corresponding event in
  the base block.

- **Execution-dependent offsets:**  
  RF phase/frequency offsets and ADC phase/frequency offsets are execution-dependent and are
  stored in `SegmentInstance`, not in the normalized base-block semantics. Implementations
  SHOULD set such offsets in normalized base blocks to zero or ignore them when
  reconstructing physical segment instances.

---

### 3.2 VirtualSegment

A virtual segment defines an ordered sequence of base block IDs constituting a reusable
sequence unit.

| Field | Type | Required | Description |
|---|---:|---:|---|
| `id` | int | required | Unique identifier for this virtual segment. Must be a positive integer. |
| `base_block_ids` | int[] | required | Ordered list of base block IDs comprising this segment. Must be non-empty. IDs may reference explicit base blocks or the reserved implicit delay IDs `0` and `1`. |
| `name` | string | optional | Human-readable descriptive label, e.g. `"tr"`, `"inversion_prep"`, `"readout"`. |

Constraints:

- `base_block_ids` MUST contain at least one entry.
- All IDs in `base_block_ids` MUST reference valid explicit base blocks or reserved IDs `0`
  and `1`.

---

### 3.3 SegmentInstance

A segment instance associates a virtual segment with concrete parameters for a single
execution in the execution stream.

| Field | Type | Required | Description |
|---|---:|---:|---|
| `virtual_segment_id` | int | required | ID of the virtual segment being instantiated. Must reference a valid virtual segment. |
| `rf_amplitude` | float[] | required | RF waveform amplitude scale factors, one per RF event in the virtual segment. Multiply by the normalized base-block RF waveform to recover physical RF amplitude. |
| `rf_phase_offset` | float[] | required | RF phase offsets in radians, one per RF event in the virtual segment. |
| `rf_frequency_offset` | float[] | required | RF transmit frequency offsets in Hz, one per RF event in the virtual segment. |
| `gradient_amplitude` | float[N_grad][3] | required | Signed gradient amplitude scale factors, one row per gradient event in the virtual segment. Each row is `[sx, sy, sz]` for physical axes `[Gx, Gy, Gz]`. Negative values indicate polarity inversion relative to the normalized base block. |
| `gradient_shot_index` | int[] | optional | Selected shot variant for each gradient event (see Section 3.1), one per gradient event. Defaults to 0 — the first/only shot — if omitted. |
| `adc_phase_offset` | float[] | required | ADC receiver phase offsets in radians, one per ADC event in the virtual segment. |
| `adc_frequency_offset` | float[] | required | ADC receiver frequency offsets in Hz, one per ADC event in the virtual segment. |
| `block_duration` | float[] | required | Pulseq block duration in seconds, one per block in the virtual segment. |
| `rf_shim` | float[N_rf][N_ch][2] | optional | Per-channel shim (magnitude, phase in radians), one set per RF event. Defaults to 1 and 0. |
| `adc_acquire` | int[] | optional | 0 if the readout is played but not stored (e.g. dummy scans), one per ADC event. Defaults to 1. |
| `trig_enable` | int[] | optional | 0 if the output trigger does not fire, one per `trig` event. Defaults to 1. |
| `rotation_matrix` | float[3][3][N_grad] | optional | One 3D spatial rotation matrix per gradient event in the virtual segment. Defaults to identity if omitted. |
| `norot` | int | optional | 1 if the prescription's rotation is not applied (Pulseq `NOROT`). Defaults to 0. |
| `nopos` | int | optional | 1 if the prescription's offset is not applied (Pulseq `NOPOS`). Defaults to 0. |
| `pmc` | int | optional | 1 if motion updates apply (Pulseq `PMC`). Defaults to 0. |
| `physio_trigger` | int | optional | Binary hardware flag, `1` or `0`, indicating whether execution must pause to await a physiological gating event before playing out this instance. Defaults to `0` if omitted. |
| `label` | string | optional | Optional execution-stream label for this instance, e.g. for slice, contrast, or repetition indexing. |

#### 3.3.1 Event counting convention

For purposes of indexing the per-event arrays in `SegmentInstance`:

- An RF event is one `PulseqBlock` containing a non-empty `rf` event.
- An ADC event is one `PulseqBlock` containing a non-empty `adc` event.
- A gradient event is one `PulseqBlock` containing at least one non-empty gradient event among
  `gx`, `gy`, and `gz`.

Thus, a `PulseqBlock` containing both `gx` and `gy` consumes one row of
`gradient_amplitude`, not two. The corresponding row contains scale factors for all three
physical axes `[sx, sy, sz]`. Missing axes are treated as absent waveforms, regardless of the
corresponding scale value.

If a virtual segment contains no RF, gradient, or ADC events, the associated arrays
(`rf_amplitude`, `gradient_amplitude`, `adc_phase_offset`, etc.) may be empty arrays but MUST
still be present as fields.

#### 3.3.2 Gradient scaling and rotation convention

For a gradient-containing base block, let the normalized base waveform vector be:

```text
G_base(t) = [Gx_base(t), Gy_base(t), Gz_base(t)]^T
```

where missing axes are treated as zero.

Let:

```text
s = [sx, sy, sz]^T
```

be the corresponding row of `gradient_amplitude`, and let `R` be the corresponding
`rotation_matrix`, or identity if `rotation_matrix` is omitted.

The physical gradient vector is:

```text
G_phys(t) = R * diag(s) * G_base(t)
```

Thus, signed per-axis amplitude scaling is applied before spatial rotation.

For a gradient event with shot variants, `G_base(t)` is the variant selected by
`gradient_shot_index`.

#### 3.3.3 Block duration convention

`block_duration` values are expressed in seconds and MUST be finite and non-negative. For
representations imported from Pulseq, these durations SHOULD correspond to the source Pulseq
block durations.

When exporting PulSeg back to a Pulseq `.seq` file, implementations MUST ensure that emitted
Pulseq block durations are compatible with the target Pulseq block-duration raster. This may
require raster-aligning generated delay blocks or splitting/padding emitted Pulseq blocks
without changing intended execution timing.

---

### 3.4 PulSeg Representation

A complete PulSeg representation consists of the following top-level fields:

| Field | Type | Required | Description |
|---|---:|---:|---|
| `pulseg_version` | string | required | Version of this specification. Must be `"2.2.0"` for representations compliant with this document. |
| `base_blocks` | BaseBlock[] | required | List of all explicit base blocks. Must be non-empty. IDs must be unique and ≥ 2. |
| `virtual_segments` | VirtualSegment[] | required | List of all virtual segments. Must be non-empty. IDs must be unique. |
| `execution_stream` | SegmentInstance[] | required | Ordered list of segment instances defining the complete scan execution. Must be non-empty. |
| `source_file` | string | optional | Path or filename of the source Pulseq `.seq` file from which this representation was generated. |
| `creation_date` | string | optional | ISO 8601 date string, e.g. `"2025-02-20"`, indicating when this representation was created. |

---

## 4. Conversion from Pulseq

### 4.1 Overview

Conversion from a Pulseq `.seq` file to PulSeg proceeds in three steps:

1. **Parse and normalize**  
   Read all Pulseq blocks and normalize RF and gradient waveform amplitudes to produce base
   blocks according to Section 3.1.

2. **Identify segments**  
   Group consecutive base blocks into virtual segments. Segment boundaries are chosen by the
   converter according to Section 4.2.

3. **Build the execution stream**  
   For each instance of each virtual segment in the Pulseq block stream, record the physical
   amplitudes, phases, frequency offsets, block durations, and supported execution metadata as
   a `SegmentInstance`.

---

### 4.2 Segment Boundaries

Segment boundaries are chosen by the converter. The source Pulseq file is not required to
carry any label for this purpose.

- A converter MUST document the rule by which it places boundaries.
- Boundaries fall between Pulseq blocks; a block is never split across segment instances.
- All instances of a virtual segment MUST have the same number of Pulseq blocks and the same
  normalized base-block structure. For gradient events, *the same structure* is defined by
  timing structure (Section 3.1): instances MAY select different shot variants via
  `gradient_shot_index`, provided the gradient timing structure is identical across instances.
- A converter SHOULD place boundaries only where every gradient is at or near zero, since a
  backend may set its gradients at the start of each instance.
- A converter MUST place a boundary wherever `norot` or `pmc` changes, since a backend sets
  one prescription transform per instance.

Different converters may segment the same file differently; each result is valid if it
satisfies Section 3. A future version may define an optional label by which a sequence
designer fixes the boundaries.

---

### 4.3 Pulseq-to-PulSeg Conversion Requirements

The following requirements apply to any compliant Pulseq-to-PulSeg conversion:

- All explicit base block waveforms MUST be normalized according to the rules in Section 3.1.
- Every explicit base block MUST have a unique ID ≥ 2.
- Reserved base block IDs `0` and `1` MUST NOT appear as entries in `base_blocks`.
- Every virtual segment MUST have a unique positive ID.
- Every virtual segment `base_block_ids` list MUST reference only valid explicit base block IDs
  or reserved IDs `0` and `1`.
- The execution stream MUST account for every Pulseq block in the source file.
- Conversion MUST be lossless with respect to the supported PulSeg execution model: waveform
  shape, RF/gradient/ADC event presence, amplitude scaling, RF/ADC phase and frequency
  offsets, block durations, virtual segment structure, and execution ordering.
- Pulseq metadata or extensions not represented by this specification, such as arbitrary
  labels, comments, definitions, or vendor-specific extensions, MAY be omitted, transformed,
  stripped, or stored as implementation-specific metadata.
- Gradient shapes that differ across instances of one virtual segment MUST be encoded per the
  shot-variant rules in Section 3.1.
- The `pulseg_version` field MUST be set to `"2.2.0"`.
- Amplitude scaling factors MUST satisfy:

  ```text
  physical amplitude = normalized amplitude × scale factor
  ```

---

### 4.4 Conversion from PulSeg to Pulseq

Conversion from PulSeg to a flattened Pulseq `.seq` stream proceeds by iterating through the
`execution_stream`.

For each `SegmentInstance`:

1. Resolve `virtual_segment_id` to a `VirtualSegment`.
2. For each `base_block_id` in `virtual_segment.base_block_ids`:
   - If the ID is `0` or `1`, emit a pure delay block with duration from the corresponding
     `block_duration` entry.
   - Otherwise, copy the referenced normalized `BaseBlock.block`.
   - Apply RF amplitude, phase offset, and frequency offset to each RF event.
   - Apply gradient amplitude scaling and optional rotation to each gradient event.
   - Apply ADC phase and frequency offsets to each ADC event.
   - Emit each output trigger whose `trig_enable` entry is `1`.
   - Emit a Pulseq block or sequence of Pulseq blocks whose total duration equals the
     corresponding `block_duration`.

A PulSeg-to-Pulseq converter MAY emit additional pure delay blocks or add duration-padding
events to satisfy Pulseq writer or raster constraints, provided the physical timing, event
ordering, and intended execution semantics are preserved.

---

## 5. Validation and Compatibility

PulSeg distinguishes between several validation and compatibility layers.

### 5.1 IR structural validation

IR structural validation verifies required fields, field types, ID uniqueness, references,
array lengths, and event counts according to this specification.

Examples include:

- `pulseg_version` is present and equals `"2.2.0"`.
- `base_blocks`, `virtual_segments`, and `execution_stream` are non-empty.
- Explicit base block IDs are unique and ≥ 2.
- Virtual segment IDs are unique and positive.
- Segment instances reference valid virtual segments.
- Per-instance arrays have lengths consistent with event counts in the referenced virtual
  segment.

### 5.2 Base-block Pulseq validation

Base-block Pulseq validation verifies that each explicit `BaseBlock.block` is a normalized
Pulseq-compatible block containing valid RF, gradient, and ADC event structures.

Such validation SHOULD check:

- presence and validity of `blockDuration`
- RF signal and timing fields
- gradient waveform or trapezoid timing fields
- ADC sample count, dwell time, and delay
- normalization of RF and gradient amplitudes
- consistency of gradient axes with `gx`, `gy`, and `gz` field names

### 5.3 Pulseq library compatibility and hydration

Different Pulseq software libraries may require additional implementation-specific fields when
constructing or writing `.seq` files. Examples include:

- event `type` identifiers
- gradient `channel` fields
- ADC `phaseModulation` fields
- library-specific default timing fields
- block-duration raster alignment

A PulSeg implementation MAY hydrate PulseqBlock events by adding harmless default fields
required by a target Pulseq library. Hydration MUST NOT change waveform shape, physical
timing, physical amplitude, phase, frequency offset, event ordering, or other PulSeg execution
semantics.

For example, setting `block.gy.channel = "y"` when the event is stored in `block.gy` is a
valid hydration step because the field name `gy` is already authoritative for PulSeg
semantics.

### 5.4 Scanner/backend validation

Scanner or backend validation verifies hardware limits, safety constraints, SAR, gradient
constraints, real-time behavior, and vendor-specific execution requirements.

Such validation is outside the scope of the PulSeg IR specification and MUST be performed by
the relevant scanner interpreter, backend compiler, or execution environment.

### 5.5 Stored representations

A stored PulSeg representation is any file from which a reader can reconstruct Section 3
without loss. It MUST state the PulSeg version it encodes (directly or through its own format
version) and declare any units or numeric formats that differ from Section 3. It MAY number
base blocks and virtual segments by its own scheme and omit the reserved IDs 0 and 1, which
the reader then assigns (Section 3.1); select shot variants by shape reference; encode the
execution stream compactly, provided it expands to the same stream; and carry additional
sections that a reader may skip. The pulserver IR cache (`.pseg`) at layout revision 0.0.7
is a stored representation of PulSeg 2.2.0.

### 5.6 Conformance

A **converter** produces from Pulseq a representation satisfying Section 4. A **consumer**
plays any representation satisfying Section 3, or refuses it with a stated reason (e.g. a
hardware limit); it MAY ignore an optional field only where that field is at its default
throughout.

---

## 6. Diagram

![Intermediate Representation](./spec-diagram.png)

*Figure 1. Schematic of the PulSeg intermediate representation, showing the relationship
between base blocks, virtual segments, and segment instances in the execution stream.*

---

## 7. Versioning and Changelog

**Current version:** `2.2.0`

Any change to this specification that affects data structure definitions, field names, field
types, required/optional status, or normalization rules must:

1. Increment the version number.
2. Add an entry to the changelog below.
3. Update the `pulseg_version` field description in Section 3.4.

Clarifications that do not change data structures or required semantics MAY be incorporated
without changing the `pulseg_version` string, but SHOULD be documented in the changelog or
repository release notes.

### Changelog

| Version | Date | Description |
|---|---:|---|
| 1.0 | 2025-02-20 | Initial release |
| 2.0-alpha | 2026-06-11 | Class definition structural upgrade; minor naming standard alignment to PyPulseq variable guide |
| 2.0-alpha clarification | 2026-06-14 | Clarified PulSeg PulseqBlock semantics, event counting, base-block validation/hydration, Pulseq export requirements, and validation layers |
| 2.1-alpha | 2026-06-18 | Clarified that a virtual segment is a reusable structural unit (removed the "TR period" example; it need not be periodic). Folded multishot gradients into the base-block model: a gradient event may carry co-timed shot variants identified by timing structure, selected per instance via the optional `gradient_shot_index` field. Additive and backward-compatible; single-shot representations are unaffected. |
| 2.2.0 | 2026-11-04 | `TRID` labels no longer required; segment boundaries are chosen by the converter. Added the optional `trig` event and instance fields `rf_shim`, `adc_acquire`, `trig_enable`, `norot`, `nopos`, `pmc`; stored representations and conformance (Sections 5.5, 5.6). A 2.1-alpha representation is valid 2.2.0 once its version string is updated. |

---

## 8. References

- Layton KJ et al. Pulseq: A rapid and hardware-independent pulse sequence prototyping
  framework. *Magn Reson Med.* 2017;77(4):1544–1552.
- [Pulseq specification](https://pulseq.github.io/)
- [PulSeg GitHub repository](https://github.com/HarmonizedMRI/pulseg)

---

## 9. Contact and Contributions

For questions, bug reports, or change requests, please open a GitHub issue:  
https://github.com/HarmonizedMRI/pulseg/issues

For correspondence regarding this specification:  
*jfnielsen@gmail.com*
