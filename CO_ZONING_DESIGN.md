> **Repository visitors:** start with `README.md`. This document is the internal design history and working notes — it records decisions, dead ends, and development context, and is not written as the public entry point.

# Paper 4: Where segregation follows workers to work — zoning, metropolitan position, and workplace segregation in the Denver region

*Target: Journal of the American Planning Association (Standard Article, 6,000 words). Formerly framed as a Denver case study for Paper 3; promoted to a standalone paper 2026-08-07.*

*Working design document, August 2026. Companion to the `6x_co_*.R` scripts in `~/Downloads/LODES/Colorado/`. Sibling to Paper 3 (`paper_pipeline/P3_DESIGN.md`), which is national and transit-focused; this one is Denver and zoning-focused.*

## Status: DESIGN + PIPELINE SCAFFOLD, GEOMETRY LOGIC VERIFIED (2026-08-06)

Scripts 60–65 written; models not yet run. The zoning×tract intersection (61's logic) was dry-run against the actual shapefile + TIGER 2024 CO tracts during scaffolding: **818 tracts intersect the layer, ~730 pass the 80% coverage screen** (Denver MSA ~660: Denver 178, Arapahoe 152, Jefferson 136, Adams 99, Douglas 72, Broomfield 24; Boulder ~69; Weld 3), fringe counties drop out at ~0% coverage as intended. Run order: 61 (zoning × tracts, ~10 min) → 62 (CO OD + wexp; fast if 51's co caches exist) → 63 → 64 → 65. Independent of the national 51/52 run — 62 downloads Colorado itself if 51's caches aren't there yet.

## The one-sentence contribution (protect this)

Paper 3 asks whether commuting infrastructure moderates the tract-level link between living in a segregated neighborhood and working in a segregated workhood; this case study asks whether the *land-use regime itself* — exclusionary residential zoning, ADU allowances, use mixing — predicts both sides of that link and moderates the connection between them, using a harmonized zoning map of the Denver region.

## Positioning

Paper 1 (MSA level) and Paper 3 (tract level) treat transport as the structural moderator. Zoning is the *upstream* structure: it fixes where housing of what density may exist and where jobs may locate, which is exactly the residential-sorting + job-geography machinery both papers gesture at. The zoning-atlas literature (Connecticut Zoning Atlas / National Zoning Atlas; Sara Bronin's work; Trounstine 2018 *Segregation by Design*; Rothwell & Massey 2009, 2010 on density zoning and segregation) links restrictive zoning to residential segregation — but not to *workplace* segregation or to the commuting bridge between the two spheres. That bridge is this study's novelty, and the LODES OD file is what makes it measurable.

## Research questions

- **RQ-Z1:** Do tracts under exclusionary zoning regimes (high share of land zoned low-density residential, low ADU allowance, low use mixing) have higher residential segregation? (Replication of the zoning–segregation literature at tract level, Denver.)
- **RQ-Z2 (headline):** Does zoning moderate the association between a tract's residential segregation and the OD-weighted workplace segregation its residents experience — is the res→workhood link tighter where zoning is exclusionary, and looser where residential land is mixed/ADU-permissive?
- **RQ-Z3 (converse):** What zoning hosts the workhoods that draw workers from segregated neighborhoods — is workplace segregation zoned into industrial land?
- **RQ-Z4 (trajectories):** Do 2023 zoning regimes predict 2011→2023 *changes* (or lock-in: persistence interactions) in residential segregation and workplace-seg exposure?
- **RQ-Z-mech:** Mechanism check via Paper 3's job-search-radius measures: does exclusionary zoning predict longer commutes, narrower destination portfolios (`eff_n_dest`), and worse jobs-housing balance?

## Data

| Piece | Source | Status |
|---|---|---|
| Zoning districts, Denver region | `~/Wellbeing/Zoning/ALL_DenverMSA_10.3.23.shp` — 100,837 polygons, 51 jurisdictions, harmonized `GenZone2` classes + ADU + GroupHome flags, Oct 2023, EPSG:32613 | **on disk** |
| Tract residential & workplace D (2011–2023) | Paper 2 panel `output/paper2/p2_tract_segregation_panel.rds` | **done, reuse** |
| Commuting flows → wexp | LODES 8.4 OD co main+aux JT01, tract-agg (reuses 51's co caches when present; else 62 downloads, ~1 min/yr) | script 62 |
| Tract covariates | Paper 2 `p2_tract_lodes_cov_panel.rds`, `p2_tract_income_panel.rds`, `tract_aland_2020.rds` | **done, reuse** |
| Tract geometry | TIGER 2024 CO tracts (2020-vintage GEOIDs), already unzipped locally | on disk |
| SLD transit measures (optional horse race) | national 53A output `p3_tract_sld.rds` | when 53 runs |

## Known data facts baked into the design

- **Coverage is the Denver region, not Colorado, and not complete.** Jurisdictions span Denver–Aurora–Lakewood (19740) plus Boulder-county (14500) and SW-Weld (24540) towns. Rural fringe counties of 19740 (Elbert, Park, Clear Creek, Gilpin) and unmapped jurisdictions are absent. Uncovered land = *unknown* zoning, never zero → tracts enter models only at ≥80% zoned coverage (`CO_MIN_ZONED_COVER`; sensitivity 50%). Main sample = Denver MSA; pooled-3-CBSA sensitivity with county FE.
- **Snapshot (Oct 2023).** Cross-sectional design at the 2023 anchor; the trajectory block (ZD) leans on zoning changing slowly, and says so. No causal claims anywhere: zoning both produces and is produced by segregation (Trounstine's point); language stays descriptive/moderation, exactly as papers 1–3.
- **GenZone2 classes** (verified against the DBF): Residential_Low (66% of polygons) /_Med /_MedHigh /_High, MixedUseRes_Low/_Med/_MedHigh/_High, MixedUse_Conditional, Commercial, Industrial, OpenSpace, Civic, MobileHome, Uncertain (570 polygons; excluded from entropy, kept as its own share). Script 61 **stops** on any unrecognized value.
- **ADU/GroupHome flags are messy free text** ("Yes", "YEs", "Yes'", "NO", "N", blank) → normalized; blanks are *not coded* (dropped from numerator and denominator), not "no".
- **Overlapping/stacked polygons are real** (verified): raw piece sums double-count up to 2.5× in the worst tract. 61 therefore *dissolves* (st_union) per tract×class before measuring, which kills within-class duplicates; cross-class conflicts (county and town zoning the same ground differently, e.g. Foxfield vs Arapahoe County) survive and are split proportionally (p99 excess ≈ 4%; conflicts >10% logged to diagnostics).
- **Coverage denominator is the tract POLYGON (land+water), not ALAND** (verified): jurisdictions zone their reservoirs (Standley Lake adds 5 km² of zoned water to one Westminster tract), so ALAND denominators push "coverage" to 2.5. Zoned water lands in the OpenSpace/Civic share.
- **Verified measure landmarks** (covered tracts): `pct_res_low` p10/50/90 ≈ 0/57/92 — wide variation, good for identification; `pct_adu_res` p10/50/90 ≈ 0/89/100 (Denver's recent ADU expansion shows; the interesting variation is in the suburbs).
- **Same-MSA OD restriction** (as 52): Denver↔Boulder cross-CBSA commuting is real, so the logged drop share will exceed the national 5–10% — reported, and consistent with Paper 3's construct.
- 2020–21 flagged pandemic years (WFH caveat), as in Papers 2–3.

## Tract-level zoning measures (script 61)

All area-weighted from the zoning×tract intersection in EPSG:5070: `pct_res_low` (share of zoned area Residential_Low — headline exclusionary proxy), `pct_reslow_of_res` (share of *residential* land that is Residential_Low — alternative denominator, zoning-atlas convention), `pct_adu_res` / `pct_grouphome_res` (share of residential land allowing ADUs / group homes), `pct_job_zone` (Commercial+Industrial+mixed-use — land where jobs can locate), `zoning_entropy` (normalized Shannon over 7 non-uncertain groups), plus per-group shares, coverage, and `jurisd_main` (largest jurisdiction — the clustering unit, because zoning is written by jurisdictions).

## Models (script 64, fixest, 55's conventions)

Cross-sections at 2023, county FE, **jurisdiction-clustered** SEs (51 clusters vs ~7 counties), z-scored variables, each zoning measure entered separately then a 3-var joint spec:

- **ZA:** res_seg ~ zoning + covs — RQ-Z1 (all three seg dimensions; plus own-tract WAC seg ~ industrial share).
- **ZB (headline):** wexp ~ res_seg × zoning + covs — RQ-Z2. Sensitivities: aspatial D, pooled 3 CBSAs, coverage ≥50%, min-commuters 50, **no-density** (worker density is plausibly a mediator of zoning — zoning caps density — so the density-conditional spec may over-control; `_nodens` variants bound the association from the other side, in both ZA and ZB), and a res_low-vs-transit horse race if SLD exists. **Exclusionary-measure comparison**: `pct_res_low` vs `pct_reslow_of_res` (r = 0.91) side-by-side across every variant + a joint horse race, written to `co_exclusionary_measure_comparison.csv` — first run showed the residential-land denominator is the sharper instrument (robust in the pooled sample where the all-zoned-area version is not).
- **ZC:** work-tract wres ~ work_seg + zoning-of-work-tract + industry-mix covs — RQ-Z3.
- **ZD:** Δres_seg and Δwexp (2011→2023) ~ 2023 zoning + baseline covs; persistence spec res23 ~ res11 × zoning — RQ-Z4.
- **ZM:** mean commute distance, `eff_n_dest`, same-tract share, jobs-housing ratio ~ zoning (± res_seg interaction) — RQ-Z-mech.

Descriptives: 2023 means by exclusionary tercile; zoning-measure correlation matrix (collinearity eyeball before reading the joint spec).

## Interpretation risks to think through before writing

- **Power**: Denver-MSA in-scope 2023 frame is a few hundred tracts. Report CIs, don't lean on stars; the jurisdiction figure (Fig 4) may carry the story better than tract regressions.
- **Sign ambiguity mirrors Paper 1's puzzle**: exclusionary zoning could *tighten* the res→workhood link (homogeneous suburbs exporting commuters into homogeneous job centers) or *loosen* it (residents of exclusive-but-integrated-workhood suburbs commute far and mix). Either is publishable as the zoning analogue of the transit surprise; the ZM block adjudicates mechanisms.
- **Res_Low ≠ single-family-exclusive exactly** — it's the harmonized low-density class. Say what it is; don't borrow "single-family zoning" language without checking the layer's coding notes.
- Denver is one region: frame as a case study demonstrating the linkage design that a National Zoning Atlas build-out could scale.

## Deliverable shells

Fig 1 three-panel map (exclusionary share | res D | wexp); Fig 2 binned res-vs-wexp scatter by exclusionary tercile; Fig 3 ZB interaction forest; Fig 4 jurisdiction scatter (% Res_Low vs res D, size = workers); Fig 5 2011–2023 trajectories by 2023 tercile. Tables: descriptives by tercile; `co_model_coefficients.csv`.


## RESULTS AND REFRAMING (2026-08-07) — read this before touching the draft

Scripts 61–67 have been RUN. Scripts 68–70 (accessibility, robustness, publication figures) were written and the analysis executed; a JAPA draft exists at `Colorado/Paper4_JAPA_draft.docx`.

### The finding that reframed the paper

Addressing the accessibility confound **changed the headline**. Job accessibility built from the 2023 WAC surface correlates **−0.928** with distance to the CBD — in monocentric Denver, "job access" and "centrality" are one construct. Letting metropolitan position moderate the res→workhood link alongside zoning collapses most of the zoning effect:

**Numbers below are from the R pipeline (authoritative), not the exploratory run.**

| moderator | naive | + metro position (**preferred**) | + local job surface (over-control) | % of baseline remaining |
|---|---|---|---|---|
| exclusionary (`reslow_of_res`) | +0.333 (p<.001) | **+0.103 (p=.035)** | +0.021 ns | 31% |
| ADU | +0.340 (p<.001) | +0.106 (p=.032) | +0.038 (p=.074) | 31% |
| use mix (entropy) | −0.224 (p<.001) | −0.063 (p=.060) **marginal** | −0.018 ns | 28% |
| job-permitting | −0.325 (p<.001) | −0.130 **ns** (p=.149) | −0.079 (p=.062) | 40% |
| exclusionary (`pct_res_low`, all zoned land) | +0.187 ns | +0.098 (p=.042) | +0.050 | 52% |

**Only the exclusionary measure remains individually significant in the position-adjusted specification.** Use mix is negative in every position specification and significant with CBD-distance alone (−0.069, p=.021) but marginal with both distance measures (p=.060) — describe its direction, do not claim an effect. `fixest`'s small-sample correction is why these p-values run slightly above the exploratory Python run; the manuscript uses the R values throughout.

Rung C conditions on the local job surface, which zoning itself partly produces → **an over-control sensitivity that may remove part of the pathway through which zoning operates; never the preferred estimate, and not a formal bound**. Rung B is the preferred specification for decomposing the association — with the explicit caveat (recorded in 68 and the manuscript) that position is not thereby a proven confounder: position, infrastructure, employment geography, and zoning co-evolved, and the ladder decomposes a cross-sectional association rather than identifying a causal zoning effect. Position-control choice barely matters for the exclusionary estimate (CBD-only +0.099 p=.086; employment-centre-only +0.147 p=.006; both +0.103 p=.035). It matters for use mix (CBD-only p=.021 vs both p=.060).

**Agreed framing (Priyanka, 2026-08-07, language tightened 2026-08-08): decomposition.** Naive zoning–segregation associations are roughly 3× the position-adjusted estimates; ~30% of the baseline moderation estimate remains after adjustment for metropolitan position (the position-adjusted zoning association — never call it a "zoning-specific component"). This is the paper's contribution and it is honest.

### Other resolved issues

- **Levels (ZA) is a clean null.** `reslow_of_res` → residential D: raw +0.175 (p=.038) → demographics +0.124 → income +0.108 → **density −0.049 ns** → position −0.053 ns. Density absorbs it entirely. Caveat in text: density is plausibly a *mediator* of zoning, so M0–M2 and M3–M4 bracket the truth. Report the whole path, not the endpoint.
- **ADU: demoted to a measurement caveat** (Priyanka's call). It is NOT a reform marker. Corr with CBD distance only +0.09; ADU-permissive tracts sit *farther* out (18.2 vs 15.4 km). Jurisdictional patchwork: Lakewood 99.6%, DouglasCo uninc. 98.3%, JeffCo 96.3%, Aurora 93.8% vs Centennial 3.7%, Westminster 9.0%, **Denver only 53.7%**. It marks suburban large-lot jurisdictions where an accessory unit is trivially accommodated. HB24-1152 (compliance 30 Jun 2025) makes a genuine before/after possible on a post-2025 vintage.
- **Effect magnitudes are small.** Preferred spec: ±1 SD residential segregation moves wexp by ~+2.6% / −3.5% of its mean. Lead with native-unit descriptives instead (below).
- **Coupling r = 0.602** in the Denver in-scope 2023 model frame. NOTE: 62's diagnostic reports ~0.67 because it covers all regional tracts in the 3 CBSAs. Do not mix the two.

### Native-unit descriptives (the practitioner-legible table = Table 1)

| | least | middle | most |
|---|---|---|---|
| % residential land low-density | 16.7 | 72.8 | 95.5 |
| residential D | 0.020 | 0.028 | 0.050 |
| workplace-seg exposure | 0.021 | 0.022 | 0.027 |
| jobs per resident worker | 3.04 | 1.14 | 1.80 |
| mean commute (km) | 13.5 | 15.1 | 19.2 |
| distance to CBD (km) | 11.8 | 16.2 | 22.7 |

### Group flows (66/67), kept as a short section

Low-earnings workers face more segregated workplaces everywhere (0.027 vs 0.022 for high earners) through far narrower destination sets (57 vs 81 effective destinations). But the zoning moderation itself is **uniform across groups** (~+0.30 to +0.36 for every earnings and industry group) — a place effect, not a class effect. Earnings thresholds are fixed nominal, so never plot SE-group trends over time.

#
## Corrections applied from the 2018-revision audit (2026-08-08)

Priyanka asked for the 2018-paper corrections to be applied here "in the same way". Per the alignment audit (`claude/paper2-2018-alignment-audit.md`): (1) the "White" group had been all NON-BLACK workers (CR01 never retained; `100 − pct_black` used as White) — fixed in 35 via the versioned cache, now `p2_tract_lodes_cov_panel_v3.rds` (adds `pct_white_{rac,wac}` from CR01 + raw group counts); (2) distance decay had been applied to unprojected degree coordinates — corrected convention is EPSG:5070, exp(−0.5·d_km), 10-km cutoff.

**Audit of Paper 4 (scripts 60–71):** the segregation indices and wexp are consumed from the corrected Paper 2 panel (script 32), so both fixes were already inherited; no Paper 4 model ever used a non-Black residual as "White" (models use `pct_black_rac` only, which was always CR02 and correct). Actions taken anyway, mirroring 36's defensive pattern: **60** now asserts the corrected conventions at load (CR01/CR02, CT02/CT01, β=0.5, 10 km, EPSG:5070) and defines `P4_LODES_COV_FILE` pointing at the v3 cache with a hard stop if only the pre-fix file exists; **63** reads v3 with a `stopifnot` on `pct_white_rac` (numerically identical panel — v3 only adds columns — so no model re-run is needed and the manuscript's numbers stand); **68** documents that the job-access gravity decay (0.10/km) is deliberately not the segregation decay. The manuscript's Methods now states the CR01/CR02 convention explicitly.

## SI live/work descriptives (script 71, added 2026-08-08)

Tables S1/S2 (in the manuscript SI + `output/models/p4_si_tabS1/S2*.csv`): where workers in each earnings bracket and industry supergroup lived and worked, Denver MSA 2023, from the same-MSA-restricted OD frame. (LODES has no occupation data — industry supergroups are the closest construct; noted in the SI text.) Headline: residence distributions are nearly identical across earnings groups (Denver county 22.0–24.8%), but workplace distributions diverge — high earners work in major employment centers at 2× the low-earnings rate (26.6% vs 13.3%); low-wage work is dispersed. The class gradient in workplace-segregation exposure is a workplace-geography fact, not a residential-sorting fact. Run 71 after 68.


## Major revision applied (2026-08-08, from Priyanka's full review)

**Status: v2 draft (`Paper4_JAPA_draft_v2.docx`) supersedes v1.** The review's verdict: novelty high, JAPA-relevant, but major-revision territory on identification language and spatial methods. All ten priority items applied:

1. **Decomposition, not confounding.** The ladder is now framed as decomposing a cross-sectional association into components tracked by regional geography vs jurisdictional regulation. All "confound," "70% of the effect disappears," "effect survives," "upper bound," "bracket the truth" language replaced ("reduces the estimated moderation by ~70%," "position-adjusted estimate ≈ one-third of baseline," "can be substantially overstated ... in Denver roughly three times larger"). Explicit statement that position/infrastructure/employment geography/zoning co-evolved and no causal subtraction is licensed.
2. **Measure sensitivity (script 72; decay part VERIFIED in cloud):** preferred exclusionary interaction is stable across decay — β=0.25: +0.096; β=0.5: +0.103; β=1.0: +0.125 — with 70–76% attenuation each time; aspatial +0.028 ns (discards the spatial structure). Maxdist 5/20 km requires local block-level recompute (72 Part B, minutes for CO; SI Table S3 has placeholder rows until run).
3. **Spatial inference (script 73; VERIFIED in cloud):** Moran's I on preferred-spec residuals is LARGE (kNN8 0.294, z=16.1; inv-dist≤10km 0.261, z=28.2). Conley Bartlett SEs: 10 km p=.042, 20 km p=.048, 50 km p=.006 vs jurisdiction-clustered p=.035. Inference stands but is marginal — manuscript says so and leans on the attenuation pattern (SI Table S4).
4. **Construct renamed** "workplace-location segregation exposure," defined against establishment-level segregation in the abstract and intro.
5. **2023 choice justified** (zoning observed once; retrospective application would assume stability) + **temporal descriptives added**: coupling r = .57–.67 across 2011–2023 with no trend (62's diagnostic); currently-exclusionary tercile already elevated in 2011 (Fig S4 = 65's trends figure, now embedded in SI).
6. **Metropolitan position reframed as a planning product** ("position is the sediment of past planning ... a finding about which scale of planning matters"), not a nuisance.
7. **Jobs mechanism disciplined:** eff_n_dest non-monotonicity (90.4/101.7/89.1) now flagged at Table 1 and used in Discussion — destination LOCATION, not number; commuting geography promoted to conceptual protagonist in title/intro.
8. **ADU section shortened** to Denver vs Douglas-uninc./Centennial contrast; full jurisdiction table pointed to SI.
9. Flagged sentences deleted/replaced ("Nothing about the naive specification is careless," "the result we did not expect...").
10. **Fig 3 is the centerpiece** ("the paper's central result" in caption); its title in 70 updated to "Estimated zoning moderation shrinks by 60–72% once metropolitan position is included."

Equation for D̃ + β/cutoff justification added to Methods. Race-flow estimand stated in the intro (destination-portfolio estimand, not group-specific destinations; earnings flows as the check). Title now leads with commuting geography. v2 ≈ 5,200 words incl. abstract/captions/tables.

**72/73 RUN AND VERIFIED (2026-08-08):** truncation sensitivity confirms stability — maxdist 5 km: naive +0.363 → preferred +0.107 (p=.040); 20 km: +0.331 → +0.104 (p=.033). Across ALL measure variants (β 0.25/0.5/1.0; cutoff 5/10/20 km) the preferred estimate spans +0.096 to +0.125 with ~58–76% attenuation. 73's Moran/Conley matched the cloud verification exactly. SI Tables S3/S4 in v2 are now fully pipeline-sourced. **Remaining before submission:** references; rerun 70 once for Fig 3's new title; the national-extension timing decision.


## Network-accessibility extension (script 75, designed 2026-08-08; reviewer-requested)

Reviewer's framing, adopted: the richer conceptual sequence is *residential segregation → metropolitan location + transportation network → accessible labor market → realized commuting destinations → workplace-location segregation*. One carefully designed analysis, not a transportation-variable expansion. **All of it is mechanism-consistent decomposition — never mediation** (infrastructure, zoning, sorting, and employment location co-evolved).

**Part A (runs now):** opportunity-set vs realized-destination decomposition. accD = impedance-weighted D̃ of each tract's *reachable* labor market (Euclidean-gravity, labeled as such); sorting gap = wexp − accD. Distinguishes the two planning stories the exemplar-pair figure raises: opportunity structure (the reachable labor market is itself segregated) vs sorting within opportunities (integrated work is reachable but actual jobs are elsewhere). The flat effective-destination counts (90.4/101.7/89.1) already rule out destination *quantity*.

**Part B (needs SLD from paper_pipeline 53A — the one remaining data download):** the accessibility ladder A (zoning only, +0.333) → B (geographic position, +0.103) → C (**network** accessibility instead: log D5AR auto-45-min jobs, log D5BR transit-45-min jobs, with res_seg interactions) → D (position + network together: does connectivity explain variation *among* neighborhoods at comparable positions?). Plus the headline test the exemplar figure invites: accessibility × res_seg on wexp, separately by mode — negative = accessibility loosens the home-to-work coupling. Falls back to the Euclidean gravity placeholder with an explicit warning if SLD is absent. Cross-link: this is Paper 3's home turf (transit moderation, national); if the Denver result is strong, it strengthens the case for the national Paper 3 rather than expanding Paper 4.

## Main-text display budget (2026-08-08)

Per Priyanka: **max 4 tables/figures in the main text.** Now: Figure 1 (zoning layer), Figure 2 (commute-flow "whose commutes end in segregated workplaces" — promoted from SI), Figure 3 (specification ladder), Table 1 (tercile descriptives). Moved to SI: three-spheres maps (S5), unadjusted coupling + native gradients (S6), exemplar pair (S7); S1–S4 unchanged; Tables S1–S4 unchanged. All in-text references renumbered and verified.

## Known gaps before submission

0. **Scripts 61–70 have all been run end-to-end and the manuscript's statistics were re-verified against their output (2026-08-07).** A bug in 68's survival pivot (rung C carries extra interaction terms → duplicate keys) was fixed; the same bug was caught in 70 before first run.
1. **No reference list.** The draft deliberately cites literature generically — no fabricated citations. Needs a real bibliography (Trounstine 2018; Rothwell & Massey 2009/2010; National Zoning Atlas / Urban Institute Connecticut work; Benner & Karner on low-wage jobs–housing fit; Reardon & O'Sullivan 2004; Hong et al. 2014; Tan & deSouza 2026).
2. Draft is ~4,300 words incl. abstract/captions — under JAPA's 6,000-word main-text ceiling, so there is room for a fuller background section.
3. Transit access (SLD, script 53A) still not run; the paper controls position via CBD/employment-centre distance and job accessibility instead. A reviewer may still ask for transit specifically.
4. Consider whether the national extension (Paper 3's OD run + more zoning atlases) should precede or follow submission — across-metro variation with metro fixed effects would identify zoning far better than one monocentric region.

## Add-on: group flows — who travels farther, into what kind of workhood (scripts 66–67, added 2026-08-06)

**LODES data facts constraining this design (state these in any writeup):** the OD file's only group flows are earnings (SE01 ≤$1,250/mo, SE02, SE03 >$3,333/mo), industry supergroup (SI01 goods, SI02 trade/transport/utilities, SI03 services), and age (SA01 ≤29, SA02 30–54, SA03 55+). **No occupation anywhere in LODES** (CTPP is the flow source for occupation — possible complement/validation, coarser geography, pooled years). **No education or race in OD**; education (CD01–04, workers 30+), full earnings (CE01–03), and detailed industry (CNS01–20) exist only as RAC/WAC margins.

Two assumption-free uses of the margins: (1) **destination-mix exposure** — flow-weighted normalized entropy of destination-tract WAC earnings/education/industry composition ("does this tract's workforce commute into class-diverse or class-monoculture workhoods?"); (2) home-tract RAC weights for education *descriptive* rows only (shared-destination assumption, labeled, same handling as Paper 3's race rows — education never enters tract-level models).

**Blocks (67):** Table G1 — distance, >24 km share, destination portfolio, wexp, destination class/edu/industry mix by group × exclusionary tercile (the "who travels farther" table). GA — within-tract class gaps (mean_dist SE01−SE03, wexp gap, portfolio gap; also age and goods−services gaps) ~ zoning: does exclusionary zoning widen the gap between the low-wage workers who *serve* a place and the high-wage workers who *live* it ("drive till you qualify" / service-worker exile)? Levels alongside gaps to see which end moves. GB — destination class-mixing ~ res_seg × zoning (the mixing analogue of ZB). GC — group-specific ZB: is the zoning amplification of the res→workhood link borne by low-earnings workers? Group cells floor at ≥10 commuters (suppression noise); partition diagnostics check SE/SI/SA families each sum to ~100% of S000.

**National pathway:** 51's national OD download already carries SE/SI/SA columns, so 66's aggregation lifts into 52 unchanged. The margin pulls generalize by looping states with tract-agg WAC/RAC downloads (~640 files, small) — or by widening `KEEP_BLOCK` in 31 and re-aggregating blocks. Denver first; if GA/GB show signal, the national version is a natural Paper 3 section or standalone equity paper.

## Run order & rough costs

61 zoning×tract (sf intersection, ~10 min) → 62 OD+wexp (fast; co downloads ~1 min/yr if 51 hasn't cached them) → 63 build (fast) → 64 models (fast) → 65 figures → 66 group flows (needs 62's OD caches; adds co WAC+RAC tract downloads, ~5 MB/yr, cached) → 67 group models/figures (needs 66 + 63) → **68 accessibility ladder** (the confound test) → **69 robustness** (levels ladder, ADU validity, effect translation) → **70 publication figures/maps** → **71 SI live/work descriptives** → **74 flow maps** → **75 network accessibility** (Part B needs 53A's SLD). Needs R packages beyond the pipeline's: `lehdr`, `fixest`, `patchwork`, `ggrepel` (+ `tigris` only if the TIGER folder moves). Nothing national is recomputed; nothing here blocks or is blocked by the national 51 run currently in progress.
