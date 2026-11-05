# Needle litter timing fix: implementation plan

Status: **validated in R, not yet implemented in C++**. Picking this up fresh?
Start at "What to implement" below — everything above it is context/evidence
already established, no need to re-derive.

## The problem this solves

`leaf_mass_correction`'s needle litter (`actual_growth.cpp`) currently drives
loss off `(1 - fS)`, gated by a sticky `fS_reached_one_leafmass` flag. This
has two confirmed problems:

1. **Duration**: `fS` (PRELES's smoothed-temperature photosynthetic-capacity
   state, `fS_model()` in `external/Rprebasso/src/gpp.c`) sits pinned at
   exactly `0.0` for the entire dormant winter (~124 days/year at Hyytiälä),
   not just a short autumn transition. Using `(1-fS)` as the senescence
   "degree" signal means litter runs near its maximum rate for ~1/3 of the
   year, not a few autumn weeks.
2. **False triggers**: the sticky flag arms once `fS` first reaches `1.0`
   (early-mid summer), then stays armed. Any later wobble in `fS` - even a
   normal, non-autumn mid-summer cool spell (`fS` dipping to ~0.83-0.87 for
   2-3 weeks, confirmed directly in the data for 2014/2016) - produces a
   second, spurious litter pulse. Confirmed via direct `fS` trace:
   `fS` reaches 1.0 on day 1 (an initialization artifact, correctly excluded
   by the `day>50` gate), dips genuinely low (0.03-0.5) through a normal
   spring ramp to day ~130, climbs to 1.0 by day ~155-160, wobbles down to
   ~0.83 for ~15-20 days (days 165-185, NOT autumn), recovers to ~1.0
   through day ~250, then genuinely declines (real autumn) from day ~250
   onward.

## What was tried and ruled out

- **Raising the `day>50` gate to `day>100`**: no effect. `fS` doesn't reach
  1.0 until ~day 155-160 regardless, so both thresholds are moot in the
  years tested.
- **A free-parameter `g_cold(T)` sigmoid + `f_L` bell**, fit from scratch
  (own `a_cold`/`b_cold`/`sLc`): kept degenerating toward "`g_cold`
  roughly constant all year, `sLc` as large as possible" under pooled-year
  SSE minimization - a real finding about the fitting method (a slowly
  rising curve can match sparse, mostly-near-zero-then-one-late-jump real
  data about as well as a genuine short pulse, in raw SSE terms), not a
  tuning failure. Confirmed via a dedicated diagnostic: this g_cold-based
  window, when window-started in **April**, picked up genuine spring cold
  snaps as false "cooling" signal (April is still cold at Hyytiälä) -
  wrongly front-loading the accumulator and putting the resulting pulse
  2+ months too early.
- **A cumulative-fraction normalization anchored to the real trap data's
  own first observation**: flagged by the user as self-referential
  ("you shouldn't be scaling the model nor the data") - fixed by
  normalizing real data against `0`, not the first observation, and (for
  magnitude checks) scaling the model by an *independent* target
  (`N_repo/n_age` from real SMEAR D,H), never by the trap data's own total.

## The validated fix: `fS^2`, not `fS`

`external/Rprebasso/src/gpp.c`'s `fPheno_model()` already contains, for
`DOY > 212`:

```c
fPheno = fS * fS;
if (fPheno < 0.5) fPheno = 0;
```

with the authors' own comment: *"Assume leaves drop much faster than fS."*
**Important caveat, confirmed by tracing all the way through**: this branch
is DEAD CODE for CASSIA's actual Scots pine parameterization. `make_p2()`
(`src/structures_read.cpp:504-506`) sets `t0 = pPREL[27]` (0-indexed), and
CASSIA's real `pPREL` vector has that value at `-999.9` - `fPheno_model`'s
`if (t0 > -998)` check fails, so it takes the `else { fPheno = 1; }`
("it is an evergreen") branch and the `fS^2` logic never actually runs in
any CASSIA simulation. **So this is "a formula sitting in the codebase that
turned out to work well when tested," not "an already-active, validated
conifer mechanism we can just hook up."** Don't claim the latter.

Despite that, using `degree = max(0, 1 - fS^2)` (instead of `1-fS`), gated
to `day >= 212` (0-indexed; the R exploration used 1-indexed `day_vec>212`),
**with no sticky flag needed at all** (fS^2's steeper decline already
discriminates the mid-summer wobble on its own - confirmed: the wobble's
`fS~0.83-0.87` gives `fS^2~0.69-0.76`, both safely above any reasonable
trigger threshold), **fixed both problems**:

- Double peaks gone in every one of 11 real years tested (2006-2016), train
  and held-out alike.
- Validation RMSE improved: 0.219 (old `fS`-based, `day>100` gate) ->
  0.175 (`fS^2`, day>=212 gate).
- Visually confirmed: exactly one clean sharp pulse per year, landing in
  the literature-documented day 250-290 (September-October) window in
  every single training AND validation year. See
  `/tmp/.../scratchpad/daily_and_cum_fPheno_{train,valid}.png` from this
  session if the scratchpad still exists (session-specific path, may be
  gone - the finding itself is what matters, re-derivable from the R
  scripts below if needed).

Fitted shape parameter: **`sLc = 2.163312`** (width of the `f_L` bell,
same functional form as `fN`: `f_L(s) = (sLc*sqrt(s) - s)/(sLc^2/4)` for
`0 < s < sLc^2`, else `0`). Fit via `optimize()` on pooled training years
(2006/08/10/12/14/16), validated on held-out years (2007/09/11/13/15).
Treat this the same way `alpha_R`/`tau_R` are treated: a first-pass,
site-specific fitted value, not independently calibrated - flag it the
same way in any new R-level argument/roxygen doc.

## Literature grounding (don't need to re-search, just cite)

- Scots pine litterfall "typically peaks during September" in Nordic
  countries, max values September-October - matches the fitted pulse
  timing directly.
  [Models of litterfall production for Scots pine in Finland](https://www.sciencedirect.com/science/article/abs/pii/S0378112704007790)
- A documented **secondary** June/July peak exists in real data too - so
  the *first* (now-eliminated) spurious peak wasn't necessarily 100%
  artifact; worth keeping in mind if future validation shows a real
  early-summer signal being lost.
  [Litterfall dynamics in Scots pine, Norway spruce and birch stands in Estonia](https://www.researchgate.net/publication/362079205_Litterfall_dynamics_in_Scots_pine_Pinus_sylvestris_Norway_spruce_Picea_abies_and_birch_Betula_stands_in_Estonia)
- Senescence, once triggered, is fast: "change from green to yellow
  occurred within less than 4 days in central Finland" - supports a sharp
  pulse over a gradual decline.
  [Long-term needle litterfall of a Scots pine stand: relation to temperature factors](https://link.springer.com/article/10.1007/BF00317216)

## PREBAS comparison (separate finding, informs magnitude not timing)

PREBAS (`external/Rprebasso`, the full forest model PRELES sits inside)
computes foliage litterfall as, annually, no within-year timing at all:

```fortran
wf_STKG = N * wf_treeKG                                    ! stand foliage mass
par_vf = par_vf0 / (1 + par_aETS*(ETS-ETS_ref)/ETS_ref)     ! turnover time, warmth-adjusted
S_fol = S_fol + wf_STKG / par_vf                            ! foliage litterfall
```

i.e. `mass / turnover_time`, same structural form as our own `N_repo/n_age`
target. Confirms the magnitude formula is the mainstream approach (two
independent sources now, not just our own derivation) - but PREBAS never
had to solve the within-year timing question CASSIA needs, since it's an
annual model. No existing timing mechanism to borrow from there either.

## Open, unresolved: magnitude overshoot (separate from timing, NOT fixed by this)

With `cum_frac` (the fitted shape) scaled by `N_repo/n_age` computed from
**real SMEAR D,H** (independent of the trap data - not a circular check),
the model consistently **overshoots** real trap totals by ~20-40% across
both complete training and validation years (e.g. 2006: 1.91 vs ~1.43 kg;
2009: 1.94 vs ~1.35 kg; 2015: 2.28 vs ~1.57 kg). 2010 is a likely
incomplete-data outlier (model more than doubles real total), not checked
further.

**This has nothing to do with `litter_boost_year` or the Makela
correction** - confirmed directly: the R exploration scaling this check
doesn't even use that mechanism. The overshoot comes purely from
`N_repo/n_age` itself (computed from Repola's biomass prediction divided
by `n_age=3`) being bigger than real observed totals. Two candidate causes,
not yet distinguished:
1. `n_age=3` may be too short a turnover time for real Scots pine needle
   retention (commonly cited 2-4 years depending on site/latitude -
   plausible either way, not obviously wrong).
2. Repola's own `N_repo` prediction may run high for this specific stand.

**Next step, not started**: test both independently (e.g. recompute the
magnitude check with `n_age=4`; separately sanity-check `N_repo` against
an independent real foliage-biomass estimate if one exists in
`smearII_data.xlsx`, e.g. `pine_foliage_biomass`, mentioned in
`paper_plan.Rmd`'s working notes as already having real inventory years).

## What to implement

Translate the validated `fS^2`/`f_L`/`sLc` mechanism into `actual_growth.cpp`
+ `growth.cpp`, replacing the current `degree=(1-fS)` / `n_age*124`
linear-rate formula. **This is not a one-line swap** - the R version's
normalization (`cum_frac = cum/max(cum)`, dividing by that YEAR'S OWN
eventual total) is retrospective and can't be computed day-by-day in a
forward simulation. The fix: mirror the **exact pattern already used for
`leaf_mass_growth`'s own `fN_normalized`** (a one-year-LAGGED
self-normalization: accumulate this year's running total, freeze it at
year-end, use it as *next* year's normalizer) rather than inventing a new
pattern.

### 1. `inst/include/function_structures.h` (`growth_state` struct)

Remove (no longer needed - `fS^2` self-discriminates without a sticky
flag):
```cpp
bool fS_reached_one_leafmass = false;
```

Add:
```cpp
// fS^2-based litter pulse (replaces the old (1-fS)/n_age*124 linear
// rate - see needle_litter_fS2_implementation_plan.md for why). sL_litter
// is the bell-curve accumulator (s_L in f_L(s_L,sLc)), reset at day==0.
// sum_degree_fL_running/last_year mirror sum_fN_running/sum_fN_last_year
// EXACTLY (same freeze-at-day-0, same "only overwrite if running>epsilon"
// bootstrap-year guard, same reasoning) - a one-year-lagged
// self-normalizer, since the R version's same-year retrospective
// normalization (divide by that year's own eventual max) can't be
// computed day-by-day in a forward simulation.
double sL_litter = 0.0;
double sum_degree_fL_running = 0.0;
double sum_degree_fL_last_year = 0.0;
```

### 2. `src/growth.cpp`

Add a bell-shape helper next to `makela_control_factor`/
`makela_loss_boost_factor` (same file, same style):
```cpp
// Same functional form as fN (growth.cpp's needle phenology curve) -
// zero at s=0, peaks at 1 when s=sLc^2/4, returns to zero past s=sLc^2.
// Self-terminating by construction, no external day-count cutoff needed.
double makela_bell_shape(double s, double sLc) {
  if (s <= 0 || s >= sLc * sLc) return 0.0;
  return (sLc * std::sqrt(s) - s) / (sLc * sLc / 4.0);
}
```

In the `if (boolsettings.leaf_mass_correction && day == 0)` block (where
`sum_fN_last_year` is frozen/reset), add the matching freeze/reset for
litter's own normalizer, SAME guard pattern:
```cpp
if (state.sum_degree_fL_running > 1e-6) {
  state.sum_degree_fL_last_year = state.sum_degree_fL_running;
}
state.sum_degree_fL_running = 0.0;
state.sL_litter = 0.0;  // fresh accumulator for the new year
```

### 3. `src/actual_growth.cpp`

Replace the whole current block (`needle_senescence_degree` through
`needle_litter_today`, roughly lines 176-221 as of this session) with:

```cpp
if (boolsettings.leaf_mass_correction) {
  // day is 0-indexed here (day==0 is Jan 1) - the R exploration used
  // 1-indexed day_vec>212, i.e. day>=212 in 0-indexed terms. No sticky
  // flag needed: fS^2's steeper decline already discriminates real
  // autumn from a mid-summer wobble on its own (confirmed: the wobble's
  // fS~0.83-0.87 gives fS^2~0.69-0.76, nowhere near triggering).
  double fS_sq = photosynthesis.fS * photosynthesis.fS;
  double degree = (day >= 212) ? std::max(0.0, 1.0 - fS_sq) : 0.0;

  tree_state.sL_litter += degree;
  double rate_unnormalized = degree * makela_bell_shape(tree_state.sL_litter, /* sLc */ boolsettings.sLc_litter);
  tree_state.sum_degree_fL_running += rate_unnormalized;

  // One-year-lagged self-normalization, same pattern as leaf_mass_growth's
  // fN_normalized. First simulated year: no prior-year total yet, hold at
  // 0 (same bootstrap treatment growth already gets via sum_fN_last_year).
  double needle_litter_today = 0.0;
  if (tree_state.sum_degree_fL_last_year > 1e-6) {
    needle_litter_today = (tree_state.N_repo_current / parameters.n_age)
                           * rate_unnormalized / tree_state.sum_degree_fL_last_year
                           * tree_state.litter_boost_year;
  }

  all_out.culm_growth.leaf_mass[days_gone+day] = std::max(0.0, all_out.culm_growth.leaf_mass[index_ref] - needle_litter_today + tree_state.leaf_mass_growth);

  all_out.culm_growth.needle_litter[days_gone+day] = all_out.culm_growth.needle_litter[index_ref] + needle_litter_today;
  all_out.needle_litter[days_gone+day] = needle_litter_today;

  tree_state.max_leaf_mass_running = std::max(tree_state.max_leaf_mass_running, all_out.culm_growth.leaf_mass[days_gone+day]);
}
```

Keep `litter_boost_year` multiplication - orthogonal, still wanted (links
litter to the same over/undershoot signal that damps growth).

### 4. New tunable parameter: `sLc_litter`

Mirror `tau_R`'s treatment exactly (`inst/include/function_structures.h`
`Settings` struct, `src/reading_settings.cpp`, `R/CASSIA_cpp.R` roxygen +
default args + settings-list entry): add `double sLc_litter = 2.163312;`
default, with a roxygen doc flagging it as fitted-but-not-independently-
calibrated, same caveat style as `alpha_R`.

### 5. Testing checklist before considering this done

- [ ] `pkgbuild::compile_dll(force = TRUE)` builds clean
- [ ] `devtools::test(".")` still `FAIL 0`
- [ ] Rerun the 15-year `leaf_mass` trajectory plot (pattern established
      earlier this session: `plot_leaf_mass_correction.R`-style script) -
      confirm still bounded, no NaN, no runaway
- [ ] Rerun the daily-rate hump check against real 2006-2016 trap data
      (reuse `fit_litter_fS_bell.R`/`plot_daily_and_cum_fPheno.R` patterns,
      but now calling the actual compiled C++ via `CASSIA_cpp()` instead
      of a parallel R reimplementation) - confirm single clean pulse per
      year, no double peaks, still lands ~day 250-290
- [ ] Check year-1 bootstrap behaviour explicitly (litter should hold
      flat, not decline, same as growth's own year-1 - verify no
      regression of the earlier "why does it go down in year 1" fix)
- [ ] Note the known ~20-40% magnitude overshoot is NOT expected to
      change from this alone (it's an `N_repo`/`n_age` question, separate
      from timing) - don't be surprised if it's still there after this
      lands; that's the next, separate piece of work
