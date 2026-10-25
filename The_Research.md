# Project Solar — The Research

## 1. Research Objectives

### 1.1 Identification of Magnetic Switchbacks

Identify magnetic switchback candidates using magnetic-field measurements from Parker Solar Probe's FIELDS instrument.

Examine large rotations of the magnetic-field vector, including reversals of its radial component ($B_R$), and investigate whether the total magnetic-field magnitude ($|B|$) remains approximately constant.

Develop quantitative criteria to distinguish switchbacks from ordinary magnetic fluctuations and other structures.

### 1.2 Characterization of Magnetic Switchbacks

Quantitatively investigate the physical and statistical properties of detected switchbacks:

- **Duration:** How long individual events last.
- **Rotation angle:** The angular deflection of the magnetic field.
- **Magnetic-field magnitude:** Changes in $|B|$ during events.
- **Occurrence rate:** How frequently switchbacks are detected.
- **Event variability:** Differences between individual switchbacks.

Identify recurring patterns and statistical distributions across events.

### 1.3 Relationships with Other Physical Measurements

Investigate correlations between switchback properties and other measurements recorded by Parker Solar Probe.

The investigation may include:

**Spacecraft position and environment**
- Heliocentric distance
- Spacecraft position and trajectory
- Changes across different solar encounters

**Plasma properties (SWEAP)**
- Proton and electron density
- Solar-wind velocity
- Proton and electron temperatures
- Temperature anisotropy
- Ion composition and distribution functions
- Electron heat flux

**Electromagnetic properties (FIELDS)**
- Magnetic-field fluctuations and turbulence
- Electric-field measurements
- Plasma-wave activity and spectral power
- Magnetic and electric-field variability

**Energetic particles (IS☉IS)**
- Energetic ion and electron intensities
- Particle energy spectra
- Directional particle distributions
- Changes in energetic-particle populations near switchbacks

**Derived physical quantities**
- Plasma beta
- Alfvén speed
- Alfvénicity
- Magnetic and plasma fluctuation statistics

The objective is to identify statistically meaningful relationships between magnetic switchbacks and their surrounding physical environment.

Measurements will be compared only where sufficient temporal coverage, resolution, and data quality are available.

Correlations will be interpreted cautiously, distinguishing statistical association from physical causation.

## 2. Scientific Background

### 2.1 The Solar Wind and Its Magnetic Field

The solar wind is a continuous flow of plasma, primarily consisting of protons, electrons, and helium ions, escaping from the solar corona into interplanetary space.

Because the solar wind is electrically conducting, its motion is strongly coupled to the interplanetary magnetic field. As the plasma travels outward, it carries magnetic structures and disturbances away from the Sun.

Parker Solar Probe measures the magnetic field in the RTN coordinate system:

- **R (Radial):** Points outward from the Sun toward the spacecraft.
- **T (Tangential):** Lies perpendicular to R, broadly in the direction of solar rotation.
- **N (Normal):** Completes the right-handed coordinate system.

The magnetic-field magnitude is:

$$
|\mathbf{B}| = \sqrt{B_R^2 + B_T^2 + B_N^2}
$$

These measurements allow us to distinguish changes in magnetic-field direction from changes in magnetic-field strength.

### 2.2 Magnetic Switchbacks

Magnetic switchbacks are large, localized deflections of the magnetic-field vector observed in the solar wind.

During a switchback, the magnetic field can rotate so dramatically that its radial component temporarily reverses relative to the surrounding magnetic-field polarity.

One particularly interesting property is that **the total magnetic-field magnitude often remains approximately constant despite the large directional change**.

Switchbacks are frequently associated with enhancements in solar-wind velocity and other Alfvénic signatures.

The observations suggest that switchbacks represent strongly deflected magnetic structures propagating through or being carried by the solar wind.

However, a reversal of $B_R$ alone does not necessarily establish the presence of a switchback.

### 2.3 Alfvén Waves and Plasma Disturbances

Alfvén waves are disturbances in magnetized plasma involving coupled fluctuations in magnetic fields and plasma velocity.

Their characteristic propagation speed is the Alfvén speed:

$$
v_A = \frac{B}{\sqrt{\mu_0 \rho}}
$$

where:

- $B$ is the magnetic-field magnitude.
- $\rho$ is the plasma mass density.
- $\mu_0$ is the vacuum permeability.

Many switchbacks exhibit Alfvénic behavior: changes in magnetic-field direction are accompanied by correlated plasma-velocity fluctuations.

Studying these relationships may help determine whether particular switchbacks behave like propagating Alfvénic structures and how they interact with the surrounding solar wind.

### 2.4 Proposed Formation Mechanisms

The physical origin of magnetic switchbacks remains an active research topic.

Several mechanisms have been proposed:

#### Interchange Magnetic Reconnection

Magnetic reconnection between open and closed magnetic-field lines near the Sun may release disturbances into the solar wind, potentially providing the initial conditions for switchback formation.

#### Growing Alfvénic Fluctuations

Magnetic fluctuations originating near the Sun may increase in relative amplitude as the solar wind expands, eventually developing into large magnetic deflections or reversals.

#### Solar-Wind Turbulence

Nonlinear interactions between magnetic and velocity fluctuations may produce strongly distorted magnetic structures.

#### Velocity Shear and Plasma Instabilities

Differences in velocity between neighboring solar-wind streams may distort magnetic structures and contribute to the development of switchbacks.

These mechanisms are not necessarily mutually exclusive. Disturbances generated near the Sun may evolve through several processes as they propagate outward.

Current research therefore investigates both the solar origins of magnetic disturbances and their subsequent evolution within the solar wind.

### 2.5 Scientific References

1. [Bale et al. (2019) — Highly structured slow solar wind emerging from an equatorial coronal hole](https://doi.org/10.1038/s41586-019-1818-7)
2. [Kasper et al. (2019) — Alfvénic velocity spikes and rotational flows in the near-Sun solar wind](https://doi.org/10.1038/s41586-019-1813-z)

## 3. Observations and Figures

This section presents observational results obtained from Parker Solar Probe measurements, including magnetic switchback detection, event characterization, statistical distributions, and relationships with other physical quantities.

### 3.1 Preliminary Observations

This section presents our initial observations from Parker Solar Probe FIELDS magnetometer measurements.

The investigation begins with one-minute magnetic-field measurements in RTN coordinates, followed by higher-resolution examination of selected magnetic-rotation candidates.

These observations are preliminary. Candidate events have not yet undergone complete physical and statistical validation.

#### Figure 1 — Initial Magnetic Field Measurements

![FIELDS magnetic-field observations](figures/fields_20181106_1min.png)

**Observation date:** November 6, 2018.

The figure shows the radial, tangential, and normal magnetic-field components, together with the total magnetic-field magnitude.

The individual components exhibit substantial temporal variations, including directional reversals.

This overview demonstrates the need for automated screening to identify interesting magnetic structures across extended observation periods.

#### Figure 2 — Promising Magnetic Switchback Candidate

![November 7 magnetic-rotation candidate](figures/highres_strong_20181107_0332.png)

**Observation date:** November 7, 2018, approximately 03:32 UTC.

The measurements reveal a substantial reversal of the radial magnetic-field component, followed by a return toward its previous direction.

The total magnetic-field magnitude appears comparatively stable over much of the interval.

This behavior is consistent with a potential magnetic switchback, although the event requires further validation.

#### Figure 3 — Additional Magnetic Switchback Candidate

![November 20 magnetic-rotation candidate](figures/highres_strong_20181120_2158.png)

**Observation date:** November 20, 2018, approximately 21:58 UTC.

This interval exhibits another pronounced magnetic-field rotation, with a reversal of the radial component.

The event provides an additional example for comparing magnetic-field behavior across different observation periods.

Its physical classification remains preliminary.

#### Figure 4 — Magnetic Rotation with Significant Magnitude Variations

![November 8 magnetic-rotation candidate](figures/highres_strong_20181108_0944.png)

**Observation date:** November 8, 2018, approximately 09:44 UTC.

Unlike the previous examples, this interval exhibits substantial variations in the total magnetic-field magnitude alongside changes in field direction.

This example illustrates why a large magnetic rotation alone is insufficient to establish the presence of a switchback.

Further examination is necessary to determine whether the variations represent physical plasma structures, instrumental effects, or other magnetic disturbances.

#### Preliminary Findings

Our initial investigation demonstrates that:

1. Parker Solar Probe FIELDS measurements contain substantial magnetic-field rotations and radial reversals.
2. Some candidate events exhibit relatively stable magnetic-field magnitude despite large directional changes.
3. Other events exhibit significant magnetic-field magnitude variations.
4. Automated screening can identify intervals for further investigation, but its classifications require validation.

These observations motivate the development of more rigorous switchback detection and characterization methods in Section 3.2.

No population-wide conclusions are drawn from these selected examples.

### 3.2 Switchback Identification and Characterization

#### 3.2.1 Detector and reference definitions

Our benchmark retains the angular detection thresholds used by Laker et al. (2022):
45° trigger, 30° boundary, 10-second minimum duration, and merging separated
segments when the intervening valid-data gap is less than 20 seconds.

For this **local-background benchmark**, the magnetic deflection angle is
calculated against a 181-minute rolling median of the magnetic RTN components:

$$
\theta(t)=\cos^{-1}\left(\frac{\mathbf{B}(t)\cdot\mathbf{B}_{bg}(t)}{|\mathbf{B}(t)|\,|\mathbf{B}_{bg}(t)|}\right).
$$

The independent magnitude reference is the 181-minute rolling median of **scalar**
$|\mathbf{B}|$, not the magnitude of the component-wise median vector.
Magnitude deviations are reported as percentages of this scalar reference.

This benchmark is **not an exact reproduction** of Laker et al.'s Parker-spiral
method. The optional `parker-spc` mode uses good-quality SWEAP SPC Level 3
primary-proton radial velocity, spacecraft HCI position, an assumed 2.5 solar-radius
source surface and local magnetic-sector polarity. It is exploratory: the published
method primarily uses SPAN-i measurements and a filtered velocity reference.

#### 3.2.2 Event catalog and data-quality safeguards

The detector uses one-second median magnetic samples; missing samples are not
interpolated. Each event records a peak angle, estimated start and end, radial
reversal, scalar-field-magnitude deviation, and reference method.

An event that begins or ends at the analysis-window boundary (or adjacent
to a missing angular sample) is **censored**: its complete duration is unknown.
Its `duration_seconds` is left blank, while `observed_span_seconds` is retained.
Duplicated detections from overlapping seed windows are consolidated in the
corrected catalog. All detections remain **unvalidated magnetic rotations**.

#### 3.2.3 Corrected preliminary results

- Selected windows successfully processed: **20**.
- Original seed-associated detections after global deduplication: **57**.
- Duplicate observations consolidated: **4**.
- Boundary/gap-censored detections: **3**.
- Median duration among only uncensored detected intervals: **223 seconds**.
- Median absolute scalar $|B|$ deviation across detected intervals: **6.13%**.

These values describe a deliberately selected set of observation windows. They are
**not** population-wide occurrence rates, event-confirmation statistics, or an
estimate of false-positive rate.

An initial exploratory SWEAP-SPC Parker-spiral comparison covered three selected windows and returned **10 unvalidated intervals**; it was later expanded to 16 common windows. It is NOT the full Laker et al. method.

#### 3.2.4 Selected diagnostic figures

![Corrected local-background detection at 20181107_0332](figures/corrected_local_strong_20181107_0332_3.png)

![Corrected local-background detection at 20181108_0944](figures/corrected_local_strong_20181108_0944_4.png)

#### 3.2.5 Cross-reference disagreement analysis

A comparison of the **three observation windows processed successfully by both
reference methods** (November 1, 5 and 7, 2018) contained 10 local-reference
detections and 10 Parker-SPC-reference detections. When detections overlapping
by at least 50% of the shorter interval were assembled into connected groups,
there were **eight overlapping groups**: four one-to-one, two where the local
reference split an SPC event, and two where the SPC reference split a local event.
**No groups were exclusive to a single method** under that overlap rule.

The earlier greedy one-to-one comparison's two local-only and two Parker-only
labels were misleading: these intervals temporally overlapped detections from
the other method but could not be matched one-to-one. The actual discrepancies
primarily concern **event segmentation and boundary placement**.

This result shows timing overlap of magnetic rotations, **not confirmation of
physical switchbacks, algorithm accuracy, or a population-wide detection rate**.
See the [full disagreement audit](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/reference_disagreement_audit.md).
The revised `compare_detectors.py` also exports overlap groups and pairs.

#### 3.2.6 Reference-angle diagnostics and next validation stage

Overlays of identical FIELDS one-second measurements against both reference
directions were generated for the four groups with event-segmentation
disagreements. Median absolute angle differences were **16.93°** and **15.25°**
for the two cases where the local detector split a Parker-SPC event, and
**7.60°** and **4.23°** for the two cases where the Parker-SPC detector
split a local event. The sample fractions exceeding the 45° trigger
threshold differed by as much as 12.19 percentage points in these
selected diagnostic windows; these fractions are **not occurrence rates**.

Qualitative review of the plots finds that differences in background
direction cause changed threshold-crossing and boundary positions.
In particular, the November 1 example dips briefly below the local
30° boundary while the Parker-SPC angle remains elevated.
Other groups are sensitive to short gaps and the event-merging rule.
See the [diagnostic results](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/reference_angle_diagnostics.csv)
and the [visual audit](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/reference_angle_findings.md).

The [review worksheet](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/reference_review_sheet.csv) remains
unlabeled. The optional published 2019 example plot was not yet
available in the shared repository, and the independent
[Martinović et al. E1–E2 catalog](https://doi.org/10.5281/zenodo.4579821)
has not yet been cross-matched. **Neither method can be declared
scientifically more accurate at this stage.**

#### 3.2.7 Remaining validation

A full replication must reproduce the published plasma-based Parker-spiral
reference (including SPAN-i primary measurements and smoothing), apply consistent
quality flags, test threshold sensitivity, independently review candidates and
non-candidates, and quantify false positives and missed events before scaling to
multi-year statistical results.

An independent positive-catalog comparison tool has now been added:
[external_switchback_validation.py](external_switchback_validation.py), with
[reproducibility and interpretation protocol](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/external_catalog_validation_protocol.md).
It targets the Martinović et al. (2021) Encounter 1–2 visually selected catalog,
using common successfully processed seed windows and one-to-one temporal matching.
The external HDF5 schema is checksum-verified. Its 5,370 entries represent
five adjacent subregions per 1,074 listed events, not 5,370 independent
switchbacks. The updated parser identifies each event from its encounter/event
number, validates five chronological regions, and selects the central
switchback interval (the published LQR–LTR–SB–TTR–TQR sequence). It uses
nanosecond-resolution timestamps from `/basic_data/block3_values` without
deserializing the file's pickled region column.

The archived catalog cross-match was subsequently completed and is reported in Sections 3.2.8–3.2.9; none of its exploratory match fractions establish detector accuracy.
Published positives alone cannot establish a false-positive rate.

#### 3.2.8 Independent catalog matching and threshold sensitivity

The first cross-match completed against Martinović et al. (2021).
Seventeen published central switchback intervals fall fully within the
three shared detector windows (November 1, 5, and 7, 2018).
With the original rule (at least 50% overlap of the shorter interval),
the local reference matches **6/17** and Parker-SPC matches **7/17**.

However, the rule counts very short published switchbacks as matched
to candidate detections lasting many minutes. A sensitivity analysis
penalizing such boundary mismatch changes the comparison:

| Matching rule | Local | Parker-SPC |
|---|---:|---:|
| 50% of shorter interval | 6/17 | 7/17 |
| Intersection-over-union >= 0.25 | 5/17 | 5/17 |
| Intersection-over-union >= 0.50 | 4/17 | 3/17 |

Seven of the 17 published events last less than 30 seconds, including
three below the current detector's ten-second minimum duration. The
original 7-versus-6 advantage is therefore **not robust to the event
matching definition**, and no detector is selected at this stage.

These are match fractions for selected windows, not statistical estimates
of mission-wide accuracy. See the
[full independent-catalog audit](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/external_catalog_validation_findings.md)
and [reproducible sensitivity script](external_match_sensitivity.py).



#### 3.2.9 Expanded shared-window validation

We expanded the Parker-SPC run to all 20 original seed windows, successfully
processing **16** with both detectors. The remaining four lacked adequate
SPC reference coverage and were excluded from the head-to-head comparison.
The 16 shared seeds form **14** merged observation windows containing **27**
published Martinović et al. central switchback intervals.

The local detector reported **48** candidate intervals within shared coverage,
while Parker-SPC reported **44**. Independent-catalog matching yielded:

| Event matching rule | Local | Parker-SPC |
|---|---:|---:|
| 50% overlap of shorter event | 11/27 | 12/27 |
| IoU >= 0.25 | **9/27** | **9/27** |
| IoU >= 0.50 | 7/27 | 5/27 |

These results confirm that **there is no stable accuracy winner** under the
criteria tested. Good matches exist (for example, published event E1-457 is
closely bounded by both detectors), but long detector intervals also overlap
published switchbacks lasting only seconds, limiting the usefulness of
the looser rule.

This selected-window, positive-catalog comparison does **not** measure
mission-wide recall, false-positive rate, or specificity. We will not
choose or freeze either algorithm until we inspect boundary failures and
review a separate, independently sampled test set including non-events.
See the [expanded comparison and limitations](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/external_catalog_validation_findings.md).


**Reference:** [Laker et al. (2022)](https://academic.oup.com/mnras/article/517/1/1001/6692877).

#### 3.2.10 Independent holdout study (prepared, not yet run)

To address selection bias from the initial 20 screened windows, a new
[held-out experiment](independent_holdout.py) and
[preregistered protocol](https://github.com/Epic-coder69/project-solar/blob/archive/pre-cleanup-2026-10-25/results/holdout_protocol.md) have been added.

The experiment samples published Martinović et al. Encounter-1 switchback
intervals that do not overlap any original analysis exposure. It also
samples uniformly timed catalog-clear intervals, restricted to the published
encounter time range and to sufficiently observed dates. A catalog-clear
interval is **unlabeled**, never assumed to be a true negative.

Selection uses a fixed pseudorandom seed, excludes all original selected
windows (including their reference-estimation margins), forbids overlapping
holdout exposures, and records a frozen manifest hash. Both reference
methods will be evaluated on the identical holdout intervals using the
unchanged 30-degree boundary, 45-degree trigger, 10-second duration,
and less-than-20-second gap-merging settings. Insufficient FIELDS or
SWEAP coverage is audited separately.

Planned outputs are publication-positive IoU comparisons at 0.25 and
0.50, catalog-clear candidate flag counts, and **blinded raw FIELDS
plots with empty human-review classifications**. These are not estimates
of mission-wide sensitivity, specificity, or false-positive rate because
the positive catalog is not exhaustive and the sampled strata are selected.

**Status: code and protocol committed; the spacecraft-data holdout run and
manual classification are not yet complete.** A detector has not been frozen.

#### 3.2.11 Independent holdout results — completed analysis, manual review pending

The frozen independent holdout contains **10 cataloged Encounter-1 switchbacks**
and **10 catalog-clear, unlabeled intervals**, excluding the original
candidate-selection exposure windows. Both magnetic-reference methods
successfully processed **all 20** held-out intervals, with at least 98.72%
valid 1-second angle data in their respective evaluation regions.

Both detector methods produce *some temporally overlapping magnetic candidate*
for all **10/10 published positives**. However, **event-boundary agreement**
is substantially lower: at IoU >=0.25, both methods match **2/10** published
events (H003/E1-104 and H006/E1-129); at IoU >=0.50 the local method matches
**2/10**, versus **1/10** for Parker-SPC. Some detector intervals are
much broader than the cataloged switchback, and edge-censored intervals
have incomplete duration. Two of the published events last only about
3.50 seconds, below the 10-second minimum-event setting in both detectors.

Both methods also flag **7/10 identical catalog-clear intervals**. Since
the external event catalog is not exhaustive, these are *unlabeled candidate
detections*, **not demonstrated false positives**.

The raw holdout review graphs have been generated, but the
[manual blind review sheet](results/holdout_blind_review.csv)
contains no classifications. Hence **no reliable false-positive rate,
sensitivity/specificity estimate, or scientific winner** is established.
The sample is also too small to justify population-wide performance claims.
See the [holdout findings and complete caveats](results/holdout_findings.md).

**Next step:** blinded physical event review with FIELDS and, where possible,
SWEAP plasma velocity, followed by a fixed-threshold assessment of
complete, non-censored boundaries before selecting a detector for statistical
analysis.

### Repository provenance

The scripts, diagnostic catalogs and intermediate plots used during
method development have been archived on the
[pre-cleanup Git branch](https://github.com/Epic-coder69/project-solar/tree/archive/pre-cleanup-2026-10-25).
Links to intermediate evidence below refer to that preserved snapshot.
The maintained FIELDS code and a small number of final catalogs/figures
remain on the main branch. The replacement four-samples-per-cycle loader
requires verification against the original local analysis environment before
claiming full rerun equivalence.

### 3.3 Statistical Results

*To be added:*

- Distribution of switchback durations.
- Distribution of magnetic-field rotation angles.
- Switchback occurrence rates.
- Comparison of event properties across observation periods.
- Statistical patterns and variability.

### 3.4 Relationships with Other Physical Measurements

*To be added:*

- Relationships with heliocentric distance.
- Solar-wind velocity and Alfvénic behavior.
- Plasma density and temperature.
- Magnetic and electric-field fluctuations.
- Plasma-wave activity and turbulence.
- Energetic-particle intensity and energy spectra.
- Derived plasma quantities and statistical correlations.
- Assessment of statistical significance and measurement limitations.

## 4. Analysis and Discussion

## 5. Conclusions

## 6. References
