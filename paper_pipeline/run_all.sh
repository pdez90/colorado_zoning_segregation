#!/bin/bash
# =============================================================================
# run_all.sh  -- rebuild the upstream pipeline, run the case study, validate.
#
#   bash paper_pipeline/run_all.sh everything    # the whole analysis, in order
#
# Works from either layout: a clone of the repository (paper_pipeline/ inside
# the case-study folder) or paper_pipeline/ and Colorado/ side by side under
# the data root. Data and caches live under LODES_ROOT (default
# ~/Downloads/LODES); set CO_ZONING_SHP to the zoning shapefile.
#
# Resumable: every step caches its outputs, so if this stops (or you stop it)
# just run it again and it picks up where it left off.
#
# Run one stage at a time instead:
#   bash run_all.sh upstream     # 10, 20, 32, 35, 53   (~1h45 cold)
#   bash run_all.sh patch        # point 60 at this machine's paths
#   bash run_all.sh casestudy    # 61, 62, 63, 64       (~30 min)
#   bash run_all.sh validate     # compare against a reference clone of the repository
#                                # (REF_DIR, default \$LODES_ROOT/colorado_zoning_segregation)
#   bash run_all.sh fresh        # rebuild SLD, refit 63-87, influence checks, number dump
# =============================================================================
set -o pipefail

LODES="${LODES_ROOT:-$HOME/Downloads/LODES}"
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "$CO_DIR" ]; then CO="$CO_DIR"
elif [ -f "$PIPE/../60_co_setup.R" ]; then CO="$(cd "$PIPE/.." && pwd)"
else CO="$LODES/Colorado"; fi
export LODES_ROOT="$LODES" PIPE_DIR="$PIPE" CO_DIR="$CO"
mkdir -p "$LODES" "$CO/clean" "$CO/output/models" "$CO/output/figures" "$CO/diagnostics"
LOGS="$PIPE/logs"
mkdir -p "$LOGS"

# --- find Rscript ------------------------------------------------------------
RS="$(command -v Rscript)"
if [ -z "$RS" ]; then
  for c in /usr/local/bin/Rscript /opt/homebrew/bin/Rscript \
           /Library/Frameworks/R.framework/Resources/bin/Rscript; do
    [ -x "$c" ] && RS="$c" && break
  done
fi
if [ -z "$RS" ]; then
  echo "Rscript not found. Install R from https://cran.r-project.org, or open"
  echo "RStudio and run:  setwd('$PIPE'); source('RUN_ALL.R')"
  exit 1
fi
echo "Using $RS"

# keep the Mac awake for the long steps
CAF=""; command -v caffeinate >/dev/null && CAF="caffeinate -i"

stage="${1:-all}"
ts() { date "+%H:%M:%S"; }

run_upstream() {
  echo "[$(ts)] upstream rebuild -> $LOGS/upstream.log"
  $CAF "$RS" "$PIPE/RUN_ALL.R" 2>&1 | tee "$LOGS/upstream.log"
}

run_patch() {
  echo "[$(ts)] patching 60_co_setup.R paths"
  "$RS" -e "source('$PIPE/patch_60_paths.R')" 2>&1 | tee "$LOGS/patch.log"
}

run_casestudy() {
  echo "[$(ts)] case study 61-64 -> $LOGS/casestudy.log"
  ( cd "$CO" && $CAF "$RS" -e "
      for (s in c('61_co_zoning_tract.R','62_co_od_wexp.R',
                  '63_co_build_panel.R','64_co_models.R')) {
        message('\n===== ', s, ' =====')
        source(s)
      }" ) 2>&1 | tee "$LOGS/casestudy.log"
}

run_sld() {
  echo "[$(ts)] Smart Location Database -> $LOGS/sld.log"
  # 53 skips when its output exists; drop a stale/empty one so it rebuilds
  rm -f "$LODES/clean/p3_tract_sld.rds" "$LODES/clean/p3_tract_sld_variants.rds"
  $CAF "$RS" -e "setwd('$PIPE'); source('53_sld.R')" 2>&1 | tee "$LOGS/sld.log"
}

run_income() {
  echo "[$(ts)] rebuilding the income panel (ACS downloads are cached)"
  rm -f "$LODES/clean/p2_tract_income_panel.rds"
  "$RS" -e "setwd('$PIPE'); source('35_covariates.R')" 2>&1 | tee "$LOGS/income.log"
}

# Recompute the segregation panel and everything downstream of it. Needed
# after any change to compute_spatial_D. ~15 min: the index is the slow part.
run_reseg() {
  echo "[$(ts)] rebuilding the segregation panel -> $LOGS/reseg.log"
  rm -f "$LODES/output/paper2"/p2_seg_*.rds \
        "$LODES/output/paper2/p2_tract_segregation_panel.rds"
  # 62 caches per-year exposures keyed off the old panel; drop them too
  rm -f "$CO/clean"/co_wexp_*.rds "$CO/clean/co_wres_work_panel.rds" \
        "$CO/clean"/co_seg_maxdist_*.rds "$CO/clean/co_group_flows_panel.rds"
  $CAF "$RS" -e "setwd('$PIPE'); source('32_segregation_panel.R')" 2>&1 \
    | tee "$LOGS/reseg.log"
  ( cd "$CO" && $CAF "$RS" -e "
      for (s in c('62_co_od_wexp.R','63_co_build_panel.R','64_co_models.R')) {
        message('\n===== ', s, ' =====')
        source(s)
      }" ) 2>&1 | tee -a "$LOGS/reseg.log"
}

# The OD caches were built without the SA01-03 age columns, which 66 and 67
# need. Rebuild them, then redo the group-flow scripts that depend on them.
run_odcols() {
  echo "[$(ts)] rebuilding OD caches with the age columns -> $LOGS/odcols.log"
  rm -f "$LODES/raw"/od_tract_*_co.rds
  rm -f "$CO/clean"/co_od_*.rds "$CO/clean/co_group_flows_panel.rds" \
        "$CO/clean/co_wac_diversity_panel.rds" "$CO/clean/co_rac_weights_panel.rds"
  $CAF "$RS" -e "setwd('$PIPE'); source('20_lodes_blocks.R')" 2>&1 \
    | tee "$LOGS/odcols.log"
  ( cd "$CO" && $CAF "$RS" -e "
      for (s in c('66_co_group_flows.R','67_co_group_models.R')) {
        message('\n===== ', s, ' =====')
        source(s)
      }" ) 2>&1 | tee -a "$LOGS/odcols.log"
}

run_models() {
  echo "[$(ts)] rebuilding panel + models (63, 64) -> $LOGS/models.log"
  ( cd "$CO" && $CAF "$RS" -e "
      for (s in c('63_co_build_panel.R','64_co_models.R')) {
        message('\n===== ', s, ' =====')
        source(s)
      }" ) 2>&1 | tee "$LOGS/models.log"
}

# income fix + SLD + refit + revalidate, in the right order
run_repair() {
  run_income && run_sld; run_models && run_validate
}

# Scripts 65-87: figures, SI tables, robustness, accessibility, decentralization,
# excess commuting, the lot. Every script is attempted so that all failures are
# listed together, but any failure makes the stage exit non-zero, which stops
# the run before the number dump: stale outputs are never passed off as new.
run_rest() {
  echo "[$(ts)] scripts 65-87 -> $LOGS/rest.log"
  ( cd "$CO" && $CAF "$RS" -e "
      # NB: no backslash escapes in this pattern. Bash collapses \\\\ to a
      # single backslash inside double quotes, so R would receive '\.R\$' and
      # reject it as an unrecognized escape. [.] does the same job untouched.
      scripts <- sort(list.files('.', pattern = '^(6[5-9]|7[0-9]|8[0-7])_.*[.]R\$'))
      message(length(scripts), ' scripts to run')
      failed <- character(0); t0 <- Sys.time()
      for (s in scripts) {
        message('\n===== ', s, ' =====')
        ok <- tryCatch({ source(s); TRUE },
                       error = function(e) {
                         message('!! FAILED ', s, ': ', conditionMessage(e))
                         FALSE })
        if (!ok) failed <- c(failed, s)
      }
      message('\n', strrep('=', 70))
      message(length(scripts) - length(failed), ' of ', length(scripts),
              ' completed in ',
              round(difftime(Sys.time(), t0, units = 'mins'), 1), ' min')
      if (length(failed)) {
        message('failed:'); for (f in failed) message('  ', f)
        quit(status = 1)
      }" ) 2>&1 | tee "$LOGS/rest.log"
}

run_validate() {
  echo "[$(ts)] validating against the committed results"
  "$RS" -e "
    Sys.setenv(CO_DIR  = '$CO',
               REF_DIR = '${REF_DIR:-$LODES/colorado_zoning_segregation}')
    source('$PIPE/90_validate_rebuild.R')
    source('$PIPE/91_validate_all_outputs.R')" 2>&1 | tee "$LOGS/validate.log"
}

# Influence and functional-form checks on the network-accessibility and
# headline interactions (75b), the tract-size robustness check (46), and the
# one-file dump of every number the pipeline produces (45).
run_influence() {
  echo "[$(ts)] influence checks, land-area robustness, number dump -> $LOGS/influence.log"
  ( cd "$CO" && $CAF "$RS" -e "source('75b_co_access_influence.R'); source('88_co_manuscript_statistics.R')" ) 2>&1 | tee "$LOGS/influence.log" || return 1
  "$RS" -e "setwd('$PIPE'); source('46_land_area_robustness.R'); source('45_manuscript_number_audit.R')" 2>&1 | tee -a "$LOGS/influence.log"
}

# Everything downstream of the SLD, from a clean SLD build. Upstream caches
# (LODES, segregation panel, income) are reused; run `all` first on a new machine.
run_fresh() {
  run_sld && run_models && run_rest && run_influence
}

# Everything that depends on the segregation index or the SLD table, from
# scratch: SLD, segregation panel, exposures, panel, models, figures, checks.
# LODES downloads, geography and ACS income are reused. ~30-40 min.
run_full() {
  run_sld && run_reseg && run_rest && run_influence
}

# The whole analysis from the raw inputs, in dependency order. Downloads are
# and the upstream segregation panel are reused if present (use `full` to
# rebuild the panel); the SLD table and every derived cache of the case study
# are cleared first so that nothing stale survives. Validation against a reference copy is separate
# (`validate`) because it needs a second clone of the repository.
run_everything() {
  rm -f "$CO/clean"/co_tract_zoning.rds "$CO/clean"/co_tract_geom.rds \
        "$CO/clean"/co_zoning_jurisd.rds "$CO/clean"/co_accessibility_2023.rds \
        "$CO/clean"/co_group_flows_panel.rds "$CO/clean"/co_wac_diversity_panel.rds \
        "$CO/clean"/co_rac_weights_panel.rds "$CO/clean"/co_wexp_*.rds \
        "$CO/clean"/co_wres_work_panel.rds "$CO/clean"/co_seg_maxdist_*.rds
  rm -f "$LODES/clean/p3_tract_sld.rds" "$LODES/clean/p3_tract_sld_variants.rds"
  run_upstream && run_casestudy && run_rest && run_influence
}

case "$stage" in
  upstream)   run_upstream ;;
  patch)      run_patch ;;
  casestudy)  run_casestudy ;;
  sld)        run_sld ;;
  income)     run_income ;;
  models)     run_models ;;
  reseg)      run_reseg && run_validate ;;
  rest)       run_rest && run_validate ;;
  odcols)     run_odcols ;;
  repair)     run_repair ;;
  validate)   run_validate ;;
  influence)  run_influence ;;
  fresh)      run_fresh ;;
  full)       run_full ;;
  refit)      run_models && run_rest && run_influence ;;
  everything) run_everything ;;
  all)        run_upstream && run_patch && run_casestudy && run_validate ;;
  *) echo "usage: bash run_all.sh [everything|full|fresh|refit|influence|upstream|patch|casestudy|sld|income|models|reseg|rest|odcols|repair|validate|all]"
     exit 1 ;;
esac

rc=$?
if [ $rc -ne 0 ]; then echo "[$(ts)] STOPPED: a step failed (see $LOGS)"; exit $rc; fi
echo "[$(ts)] done. Logs in $LOGS"
