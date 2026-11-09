# Make `fS` available to all photosynthesis models, not just PRELES

Status: **planned, not started**. Picking this up fresh? Start at
"What to implement" below.

## The problem this solves

`leaf_mass_correction`'s needle litter mechanism (`actual_growth.cpp`)
reads `photosynthesis.fS` - PRELES's smoothed-temperature
photosynthetic-capacity state (`fS_model()`, `external/Rprebasso/src/gpp.c`)
- as its senescence-timing signal (`degree = max(0, 1 - fS^2)`, see
`needle_litter_fS2_implementation_plan.md`). But `fS` is only ever
actually computed inside `preles_cpp()`'s own call chain.

Confirmed in `src/CASSIA.cpp`'s year loop:

```cpp
photosynthesis.fS = 0.0;   // reset once per YEAR (not per day - itself
                            // worth double-checking while touching this)
...
if (boolsettings.photosynthesis_as_input) {
  photosynthesis.GPP = climate.Photosynthesis_IN[day + days_gone];
  photosynthesis.ET = 0.0;
  photosynthesis.SoilWater = 0.0;
  // fS never touched here - stays 0.0 all year
} else if (boolsettings.preles) {
  photosynthesis = preles_cpp(day, reset_state, climate.PAR[...], ...);
  // preles_cpp() computes fS internally and overwrites the whole struct
} else if (boolsettings.phydro) {
  // fS never touched here either - stays 0.0 all year
}
```

**Practical consequence**: `leaf_mass_correction = TRUE` only produces a
real needle-litter pulse when `preles = TRUE` is also set. Combined with
`photosynthesis_as_input` (the project's long-standing default) or
`phydro`, the litter mechanism silently does nothing (`degree` is always
`max(0, 1-0^2) = 1` once `day > 212`, so the bell-curve accumulator `s`
still runs, but driven by a constant `degree=1` rather than a real
temperature-acclimation signal - almost certainly wrong, not just
"disabled"). This hasn't been checked/confirmed empirically yet - do that
first (see below) before assuming it's purely cosmetic.

## Why this is tractable

`fS_model()` (`external/Rprebasso/src/gpp.c:4`) depends on **only**:

```c
double fS_model(double *S, double T, p2 GPP_par) {
  *S = *S + (T - *S) / GPP_par.tau;
  double fS = (*S - GPP_par.S0 > 0) ? (*S - GPP_par.S0) : 0;
  return (fS / GPP_par.Smax > 1) ? 1 : fS / GPP_par.Smax;
}
```

- `T`: today's air temperature - already available every day regardless
  of photosynthesis model (`climate.TAir[days_gone+day]`, read once per
  day in `CASSIA.cpp` before the photosynthesis dispatch).
- `GPP_par` (`p2`, `tau`/`S0`/`Smax`): already built **unconditionally**
  near the top of `CASSIA.cpp` - `p2 parGPP = make_p2(pPREL);` (line ~92,
  before the year loop) - not gated behind `boolsettings.preles` at all.
- `S` (`*S`, the running temperature-sum state): currently a `static`
  local inside `preles_cpp()` (`src/preles_test.cpp`), separate from
  PRELES's photosynthesis state proper. Needs its own home if decoupled.

So nothing about `fS_model()` itself is PRELES-specific - it's only
wired to PRELES by where it happens to be *called from*, not by any
real data dependency.

## What to implement

1. **Give the temperature-sum state (`S`, for `fS_model`) a persistent
   home outside `preles_cpp()`'s static.** Add a field to `growth_state`
   (`inst/include/function_structures.h:39`), e.g. `double fS_temp_sum =
   0.0;` (name it to avoid confusion with `preles_cpp()`'s own internal
   `S_state`, which should stay as-is for the `preles` pathway itself -
   don't just reuse/alias it, or `soil_water_reset`-style annual-reset
   bugs become possible again).
2. **Call `fS_model()` unconditionally, once per day, in `CASSIA.cpp`'s
   year loop** - regardless of which photosynthesis branch runs that day
   - using `climate.TAir[days_gone+day]`, `parGPP`, and the new
   `tree_state.fS_temp_sum`. Decide: does this *replace* `preles_cpp()`'s
   own internal `fS` computation for the `preles` branch (i.e. one
   shared `fS` computed the same way for every pathway), or does
   `preles_cpp()` keep computing its own copy independently (two
   parallel computations, same formula, same `T`/`GPP_par` inputs, only
   differing by which `S` state they carry - should be numerically
   identical if both are seeded/reset the same way, but confirm this
   rather than assuming it)? The cleaner option is almost certainly the
   former (one shared computation) - avoids two static-vs-struct copies
   of the same state silently drifting apart - but check
   `preles_cpp()`'s own use of `fS` internally (it's also an input to
   `GPPfun()`, not just an output) before ripping it out of there.
3. **Reset timing**: decide what "day==0" should mean for
   `tree_state.fS_temp_sum` - probably nothing (not reset at all,
   continuous like `soil_water_reset=FALSE`'s state, since `fS` is a
   genuine multi-year-continuous physiological signal, not an
   annually-reset one - same reasoning already applied to PRELES's own
   `S_state` in `soil_water_reset`'s fix). Don't reuse the `soil_water_reset`
   toggle for this though - it's a different state variable with its own
   (likely simpler) correctness story; decide fresh rather than assuming
   the same toggle should gate both.
4. **Check the annual `photosynthesis.fS = 0.0` reset** (`CASSIA.cpp`,
   top of year loop, "TODO: think about this bit!" comment already
   flags this as unexamined) - once `fS` is computed continuously
   outside PRELES, does this annual reset need to go too? Confirm
   whether anything currently relies on this per-year reset before
   removing it.
5. **Verify empirically, the same way `needle_litter_fS2_implementation_plan.md`
   did**: run `leaf_mass_correction = TRUE` with `photosynthesis_as_input
   = TRUE` (not `preles`) before and after, check `fS`/`needle_litter`
   traces for a sane seasonal pulse (not the constant-`degree=1`
   behavior described above), and confirm the `preles = TRUE` pathway's
   own `needle_litter`/`fS` output is **unchanged** (byte-identical) by
   this refactor - same regression-safety bar as every other change this
   session.

## Where to validate against real data

`development_notebooks/leaf_mass_correction.Rmd` already has everything
needed to re-run this check once implemented: real `karike_df_all`
litter-trap data (Sect. 2.3/4/5), the fitted `sLc_litter`/bell-shape
comparison (Sect. 6), and the `photosynthesis_as_input` vs. `preles`
distinction is already called out in the doc's own intro ("this matters
here specifically because `photosynthesis.fS`... is only ever actually
computed inside `preles_cpp()`'s daily call"). Re-run the notebook's
model call with `photosynthesis_as_input = TRUE` once this is wired up,
as an additional validation pass alongside the existing `preles = TRUE`
one.
