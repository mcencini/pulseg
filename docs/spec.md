![pulseg logo](../assets/logo.svg)

# PulSeg Intermediate Representation Specification

**Version:** 2.2.0  
**Date:** 2026-11-04  
**Status:** Draft proposal  
**Authors:** Jon-Fredrik Nielsen, Matteo Cencini  
**Repository:** https://github.com/HarmonizedMRI/pulseg

---

## 1. Overview

PulSeg is a vendor-neutral intermediate representation (IR) for MRI pulse sequences. 
It is designed to sit between a high-level sequence description format (such as Pulseq) 
and a vendor-specific hardware execution layer. By making explicit the repeating structural 
units present in most MRI sequences, PulSeg enables efficient and unambiguous mapping onto 
the internal execution models of real scanner hardware.

This document defines the PulSeg data structures, terminology, and conversion requirements. 
It is intended for developers of sequence conversion tools, scanner interpreters, and 
simulation frameworks.

This specification is stable at version 2.2.0. Any change affecting data structure definitions, 
field semantics, or required fields must increment the version number and include 
migration notes in the changelog (see Section 6).

Implementations MAY include additional fields not defined in this specification. Compliant
readers SHOULD ignore unknown fields unless they conflict with required PulSeg semantics.
Such fields MUST NOT be required for interpreting the core IR.

---

## 2. Definitions

- **Pulseq block:**  
  An atomic time unit in the [Pulseq](https://pulseq.github.io/) sequence format,
  containing at most one waveform per gradient axis, one RF waveform, and one ADC window. 
  Waveform amplitudes in a Pulseq block
  reflect the actual physical amplitudes used during that block.

- **Base block:**  
  A Pulseq block in which all waveform amplitudes have been normalized to a canonical value 
  (see Section 3.1). A base block defines the *shape* of the waveforms but not their 
  physical amplitudes. Base blocks are the atomic elements of a PulSeg representation.

- **Virtual segment:**  
  An ordered, finite sequence of base blocks representing a generic, reusable unit of the 
  MRI sequence — for example, a readout module, a spoiler module, or a contrast preparation
  module (it need not be periodic). A virtual segment 
  defines the *structure* of a sequence unit but not the specific amplitudes, phases, or 
  frequency offsets used in any particular execution.

- **Segment instance:**  
  A concrete realization of a virtual segment within the execution stream. A segment instance 
  associates a virtual segment with specific waveform amplitudes, RF phase offsets, and 
  frequency offsets for a single occurrence in the scan.

- **Execution stream:**  
  The ordered sequence of segment instances that defines the complete execution of the 
  MRI sequence on the scanner.

---

## 3. Data Structures

All field names use snake_case. Fields marked **(required)** must be present in any 
compliant PulSeg representation. Fields marked **(optional)** may be omitted.

Durations are in seconds, phase offsets in radians, frequency offsets in Hz, and RF and
gradient amplitude scale factors in the native units of the source sequence. A stored
representation MAY use other units and numeric formats if it declares them (Section 3.5).

### 3.1 BaseBlock

A base block wraps a single Pulseq block with normalized waveform amplitudes. A base block is
identified by the normalized *shapes* of its RF and ADC events and by the *timing* of its
gradient events. Because gradient identity is timing-based, a single gradient event may carry
more than one co-timed waveform *shape* — a *shot* dimension — much as a multi-channel (pTx) RF
event carries more than one transmit channel (see the normalization rules).

| Field | Type | Required | Description |
|---|---|---|---|
| `id` | int | required | Unique identifier for this base block. Must be a non-negative integer. IDs 0 and 1 are globally reserved (see Reserved Identifiers below). User-defined blocks must use IDs ≥ 2. |
| `block` | PulseqBlock | required | A Pulseq block with all waveform amplitudes normalized to 1.0 (or 0.0 for unused channels). See normalization rules below. |
| `name` | string | optional | Human-readable descriptive label (e.g., `"rf_prep"`, `"gx_spoil"`). |

A base block contains at most one each of `rf`, `gx`, `gy`, `gz`, `adc` and `trig`, where
`trig` is an output trigger (channel, delay, duration). Its delay and duration are part of the
block's timing, which a backend prepares with the segment; whether it fires is set per
instance by `trig_enable`. A base block with ID ≥ 2 MUST contain at least one event; pure
delays are represented only by the reserved IDs below and do not appear in `base_blocks`.

**Reserved Identifiers:**  
To support real-time timing optimizations on scanner hardware, the following base block IDs are globally reserved and implicitly defined across all PulSeg implementations:

- id == 0 (Implicit Constant Delay Block): Represents a pure delay block whose duration remains fixed and invariant across all instances in the execution stream. This allows hardware interpreters to pre-compute structural timing gaps.

- id == 1 (Implicit Variable Delay Block): Represents a pure delay block whose duration can vary dynamically between different instances in the execution stream (e.g., for extending TE or TR). The runtime duration is dictated explicitly by the corresponding entry in the instance's block_duration array.

Whether a pure delay is constant or variable follows from the execution stream: it is
variable if its `block_duration` differs between instances of its virtual segment. A stored
representation MAY therefore omit the reserved IDs and let the reader assign them on loading
(Section 3.5).

**Normalization rules:**
- Single-channel RF waveforms: Normalize by peak magnitude, such that `max(|rf.signal|) == 1.0`.
- Multi-channel RF waveforms (Parallel Transmit / pTx): Normalize all transmit channels by a single global scaling factor 
derived from the maximum peak magnitude across all channels combined, such that 
$\max_{c}(\max(\lvert\text{rf.signal}_c\rvert)) == 1.0$. 
This ensures that the relative amplitude between distinct physical transmit coils are strictly preserved.
- Gradient waveforms: normalize each channel independently by peak absolute amplitude, such that `max(|grad.waveform|) == 1.0` for each channel. A gradient event MAY carry multiple co-timed *shot variants* — gradient shapes sharing one timing structure but differing in sample values, e.g. independently optimized spiral interleaves — and each shot variant is normalized independently (peak `1.0` per variant). Two gradient waveforms share a *timing structure* when they have identical `(delay, rise, flat, fall)` for trapezoids, `(delay, time_shape_id)` for extended trapezoids, or `(delay, num_samples)` on the gradient raster for uniform-raster arbitrary gradients. Shot variants are reserved for shapes NOT related by scaling or rotation; families related by scaling or rotation (phase-encode blips, radial spokes, rotated interleaves) MUST use a single shot variant scaled via `gradient_amplitude` / `rotation_matrix`. The single-shot case is the common one.
- ADC windows: not normalized; copied directly from the Pulseq block
- Output triggers: not normalized; copied directly from the Pulseq block
- A channel with no waveform in the original block must have no waveform in the base block

### 3.2 VirtualSegment

A virtual segment defines an ordered sequence of base blocks constituting a reusable 
sequence unit.

| Field | Type | Required | Description |
|---|---|---|---|
| `id` | int | required | Unique identifier for this virtual segment. Must be a positive integer. |
| `base_block_ids` | int[] | required | Ordered list of base block IDs comprising this segment. Must be non-empty. All referenced IDs must exist in the base block list. |
| `name` | string | optional | Human-readable descriptive label (e.g., `"tr_delay"`, `"inversion_prep"`). |

**Constraints:**
- `base_block_ids` must contain at least one entry
- All IDs in `base_block_ids` must reference valid base blocks

### 3.3 SegmentInstance

A segment instance associates a virtual segment with the concrete parameters for a 
single execution in the execution stream.

| Field | Type | Required | Description |
|---|---|---|---|
| `virtual_segment_id` | int | required | ID of the virtual segment being instantiated. Must reference a valid virtual segment. |
| `rf_amplitude` | float[] | required | Scaling factors for RF waveform amplitudes, one per RF event in the virtual segment. Multiply by the normalized base block RF amplitude to recover the physical amplitude. |
| `rf_phase_offset` | float[] | required | RF phase offsets in radians, one per RF event in the virtual segment. |
| `rf_frequency_offset` | float[] | required | RF transmit frequency offset in Hz, one per RF event in the virtual segment. |
| `gradient_amplitude` | float[N_grad][3] | required | Signed scaling factors for gradient amplitudes, one row `[sx, sy, sz]` per gradient event in the virtual segment. Negative values indicate polarity inversion relative to the normalized base block. |
| `gradient_shot_index` | int[] | optional | Selected shot variant for each gradient event (see BaseBlock, Section 3.1), one per gradient event. Defaults to 0 — the first/only shot — if omitted. |
| `rf_shim` | float[N_rf][N_ch][2] | optional | Per-channel shim (magnitude, phase in radians), one set per RF event. Defaults to 1 and 0. |
| `adc_phase_offset` | float[] | required | ADC receiver phase offsets in radians, one per ADC event in the virtual segment. |
| `adc_frequency_offset` | float[] | required | ADC receiver frequency offset in Hz, one per ADC event in the virtual segment. |
| `adc_acquire` | int[] | optional | 0 if the readout is played but not stored (e.g. dummy scans), one per ADC event. Defaults to 1. |
| `trig_enable` | int[] | optional | 0 if the output trigger does not fire, one per `trig` event. Defaults to 1. |
| `block_duration` | float[] | required | Pulseq block duration in seconds, one per block in the virtual segment. |
| `rotation_matrix` | float[3][3][N_grad] | optional | 3D spatial rotation matrices applied to the gradient axes, one per gradient event in the virtual segment. Defaults to identity if omitted.|
| `norot` | int | optional | 1 if the prescription's rotation is not applied (Pulseq `NOROT`). Defaults to 0. |
| `nopos` | int | optional | 1 if the prescription's offset is not applied (Pulseq `NOPOS`). Defaults to 0. |
| `pmc` | int | optional | 1 if motion updates apply (Pulseq `PMC`). Defaults to 0. |
| `physio_trigger` | int | optional | Binary hardware flag (1 or 0) indicating whether execution must pause to await a physical gating event (e.g., ECG R-wave or respiratory trigger) before playing out this instance. Defaults to 0 if omitted. |
| `label` | string | optional | Optional execution stream label for this instance (e.g., for slice or contrast indexing). |

**Notes:**
- If a virtual segment contains no RF/gradient/adc events, the associated columns (e.g., `rf_amplitude`, `rf_phase_offset`, etc)
  may be empty arrays but must still be present as fields
- Physical amplitude = base block normalized amplitude × scaling factor
- For multishot base gradients, the per-instance physical gradient uses the shot variant
  selected by `gradient_shot_index`: physical gradient = rotation × signed-scale ×
  (selected-variant normalized waveform)
- Event counting: an RF event is one block with an `rf` event, an ADC event one block with
  an `adc` event, and a gradient event one block with at least one of `gx`, `gy`, `gz`. A
  block with both `gx` and `gy` consumes one row of `gradient_amplitude`.
- Scaling is applied before rotation: `G_phys(t) = R * diag(s) * G_base(t)`, with `s` the row
  of `gradient_amplitude`, `R` the rotation (identity if omitted) and missing axes of
  `G_base` zero.

### 3.4 PulSeg Representation (Top-Level Structure)

A complete PulSeg representation consists of the following top-level fields:

| Field | Type | Required | Description |
|---|---|---|---|
| `pulseg_version` | string | required | Version of this specification. Must be `"2.2.0"` for representations compliant with this document. |
| `base_blocks` | BaseBlock[] | required | List of all base blocks. Must be non-empty. IDs must be unique. |
| `virtual_segments` | VirtualSegment[] | required | List of all virtual segments. Must be non-empty. IDs must be unique. |
| `execution_stream` | SegmentInstance[] | required | Ordered list of segment instances defining the complete scan execution. Must be non-empty. |
| `source_file` | string | optional | Path or filename of the source Pulseq `.seq` file from which this representation was generated. |
| `creation_date` | string | optional | ISO 8601 date string (e.g., `"2025-02-20"`) indicating when this representation was created. |

### 3.5 Stored Representations

A stored PulSeg representation is any file from which a reader can reconstruct Section 3
without loss. It MUST state the PulSeg version it encodes (directly or through its own format
version) and declare any units or numeric formats that differ from Section 3. It MAY number
base blocks and virtual segments by its own scheme and omit the reserved IDs 0 and 1, which
the reader then assigns (Section 3.1); select shot variants by shape reference; encode the
execution stream compactly, provided it expands to the same stream; and carry additional
sections that a reader may skip. The pulserver IR cache (`.pseg`) at layout revision 0.0.7
is a stored representation of PulSeg 2.2.0.

### 3.6 Conformance

A **converter** produces from Pulseq a representation satisfying Section 4. A **consumer**
plays any representation satisfying Section 3, or refuses it with a stated reason (e.g. a
hardware limit); it MAY ignore an optional field only where that field is at its default
throughout.

---

## 4. Conversion from Pulseq

### 4.1 Overview

Conversion from a Pulseq `.seq` file to PulSeg proceeds in three steps:

1. **Parse and normalize** — Read all Pulseq blocks; normalize waveform amplitudes to 
   produce base blocks (see Section 3.1 normalization rules).
2. **Identify segments** — Group consecutive base blocks into virtual segments. Segment 
   boundaries are chosen by the converter (see Section 4.2).
3. **Build the execution stream** — For each instance of each virtual segment in the Pulseq 
   block stream, record the physical amplitude, phase, and frequency parameters as a 
   segment instance.


### 4.2 Segment Boundaries

Segment boundaries are chosen by the converter. The source Pulseq file is not required to
carry any label for this purpose.

- A converter MUST document the rule by which it places boundaries.
- Boundaries fall between Pulseq blocks; a block is never split across segment instances.
- All instances of a virtual segment MUST have the same number of Pulseq blocks and the same normalized base-block structure. For gradient events, *the same structure* is defined by timing structure (Section 3.1): instances MAY select different shot variants via `gradient_shot_index`, provided the gradient timing structure is identical across instances.
- A converter SHOULD place boundaries only where every gradient is at or near zero, since a
  backend may set its gradients at the start of each instance.
- A converter MUST place a boundary wherever `norot` or `pmc` changes, since a backend sets
  one prescription transform per instance.

Different converters may segment the same file differently; each result is valid if it
satisfies Section 3. A future version may define an optional label by which a sequence
designer fixes the boundaries.

### 4.3 Conversion Requirements

The following requirements apply to any compliant Pulseq-to-PulSeg conversion:

- All base block waveforms MUST be normalized according to the rules in Section 3.1
- Every virtual segment MUST have a unique ID
- Every base block MUST have a unique ID
- The execution stream MUST account for every block in the source Pulseq file (conversion is lossless)
- The `pulseg_version` field MUST be set to the version of this specification
- Amplitude scaling factors MUST be such that: physical amplitude = normalized amplitude × scale factor
- Gradient shapes that differ across instances of one virtual segment MUST be encoded per the shot-variant rules in Section 3.1 (shot variants only for shapes not related by scaling or rotation; scale/rotation families use a single shot variant).
- Pulseq labels and extensions not represented by this specification MAY be omitted or kept
  as additional fields.

---

## 5. Diagram

![Intermediate Representation](./spec-diagram.png)

*Figure 1. Schematic of the PulSeg intermediate representation, showing the relationship 
between base blocks, virtual segments, and segment instances in the execution stream.*

---

## 6. Versioning and Changelog

**Current version:** 2.2.0

Any change to this specification that affects data structure definitions, field names, 
field types, required/optional status, or normalization rules must:
1. Increment the version number (patch increment for clarifications, minor increment for 
   additive changes, major increment for breaking changes)
2. Add an entry to the changelog below
3. Update the `pulseg_version` field description in Section 3.4

### Changelog

| Version | Date | Description |
|---|---|---|
| 1.0 | 2025-02-20 | Initial release |
| 2.0-alpha | 2026-06-11 | Class definition structural upgrade; minor naming standard alignment to PyPulseq variable guide |
| 2.1-alpha | 2026-06-18 | Clarified that a virtual segment is a reusable structural unit (replaced the "TR period" example; it need not be periodic). Folded multishot gradients into the base-block model: a gradient event may carry co-timed shot variants identified by timing structure, selected per instance via the optional `gradient_shot_index` field. Additive and backward-compatible; single-shot representations are unaffected. |
| 2.2.0 | 2026-11-04 | `TRID` labels no longer required; segment boundaries are chosen by the converter. Added the optional `trig` event and instance fields `rf_shim`, `adc_acquire`, `trig_enable`, `norot`, `nopos`, `pmc`; stored representations and conformance (Sections 3.5, 3.6). Added the unknown-fields clause, units, event counting, scaling before rotation and array shapes. A 2.1-alpha representation is valid 2.2.0 once its version string is updated. |

---

## 7. References

- Layton KJ et al. Pulseq: A rapid and hardware-independent pulse sequence prototyping 
  framework. *Magn Reson Med.* 2017;77(4):1544–1552.
- [Pulseq specification](https://pulseq.github.io/)
- [PulSeg GitHub repository](https://github.com/HarmonizedMRI/pulseg)

---

## 8. Contact and Contributions

For questions, bug reports, or change requests, please open a GitHub issue:  
https://github.com/HarmonizedMRI/pulseg/issues

For correspondence regarding this specification:  
*jfnielsen@gmail.com*

