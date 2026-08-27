#include "CASSIA.h"

/*
 * Potential Growth
 */

void log_potential_growth(int day,
                          int days_gone,
                          const growth_state& tree_state,
                          output_vector& growth_out) {

  growth_out.potential_height[days_gone + day] = tree_state.height;
  growth_out.potential_needles[days_gone + day] = tree_state.needles;
  growth_out.potential_roots[days_gone + day] = tree_state.roots;
  growth_out.potential_diameter[days_gone + day] = tree_state.diameter;
  growth_out.potential_wall[days_gone + day] = tree_state.wall;
  growth_out.potential_bud[days_gone + day] = tree_state.bud;
  growth_out.g[days_gone + day] = tree_state.g;
  growth_out.en_pot_growth[days_gone + day] = tree_state.use;
  growth_out.potential_ring_width[days_gone + day] = tree_state.pot_mm;
}


/*
 * Actual Growth
 */

void log_actual_growth(int day,
                       int days_gone,
                       const growth_state& tree_state,
                       output_vector& growth_out) {
  growth_out.height[days_gone + day]   = tree_state.height;
  growth_out.needles[days_gone + day]  = tree_state.needles;
  growth_out.roots[days_gone + day]    = tree_state.roots;
  growth_out.diameter[days_gone + day] = tree_state.diameter;
  growth_out.wall[days_gone + day]     = tree_state.wall;
  growth_out.bud[days_gone + day]      = tree_state.bud;
}

/*
 * Sugar
 */

void log_sugar(int day,
               int days_gone,
               const carbo_tracker& sugar,
               const carbo_tracker& starch,
               const carbo_tracker& storage_term,
               const growth_out& nitrogen_capacity,
               const double respiration_growth,
               const double respiration_maintenance,
               const double nitrogen_balance,
               const uptake_structre& uptake,
               const bool& tree_alive,
               output_vector& out) {

  // Sugar vector
  out.sugar_vector.needles[days_gone + day]     = sugar.needles;
  out.sugar_vector.phloem[days_gone + day]      = sugar.phloem;
  out.sugar_vector.roots[days_gone + day]       = sugar.roots;
  out.sugar_vector.xylem_sh[days_gone + day]    = sugar.xylem_sh;
  out.sugar_vector.xylem_st[days_gone + day]    = sugar.xylem_st;
  out.sugar_vector.surplus[days_gone + day]     = sugar.surplus;
  out.sugar_vector.to_mycorrhiza[days_gone + day] = sugar.mycorrhiza;

  // Total sugar
  out.sugar[days_gone + day] = sugar.needles + sugar.phloem + sugar.roots + sugar.xylem_sh + sugar.xylem_st;

  // Starch vector
  out.starch_vector.needles[days_gone + day]   = starch.needles;
  out.starch_vector.phloem[days_gone + day]    = starch.phloem;
  out.starch_vector.roots[days_gone + day]     = starch.roots;
  out.starch_vector.xylem_sh[days_gone + day]  = starch.xylem_sh;
  out.starch_vector.xylem_st[days_gone + day]  = starch.xylem_st;

  // Total starch
  out.starch[days_gone + day] = starch.needles + starch.phloem + starch.roots + starch.xylem_sh + starch.xylem_st;

  // Storage terms
  out.storage_term_vector.needles[days_gone + day]   = storage_term.needles;
  out.storage_term_vector.phloem[days_gone + day]    = storage_term.phloem;
  out.storage_term_vector.roots[days_gone + day]     = storage_term.roots;
  out.storage_term_vector.xylem_sh[days_gone + day]  = storage_term.xylem_sh;
  out.storage_term_vector.xylem_st[days_gone + day]  = storage_term.xylem_st;

  // Respiration
  out.respiration_output.growth[days_gone + day]      = respiration_growth;
  out.respiration_output.maintenance[days_gone + day] = respiration_maintenance;

  // Nitrogen capacity
  out.nitrogen_capacity_vector.needles[days_gone + day]  = nitrogen_capacity.needles;
  out.nitrogen_capacity_vector.bud[days_gone + day]      = nitrogen_capacity.bud;
  out.nitrogen_capacity_vector.roots[days_gone + day]    = nitrogen_capacity.roots;
  out.nitrogen_capacity_vector.diameter[days_gone + day] = nitrogen_capacity.wall;
  out.nitrogen_capacity_vector.height[days_gone + day]   = nitrogen_capacity.height;

  // Uptake structure
  out.uptake_vector.ectomycorrhizal_transfer[days_gone + day] = uptake.ectomycorrhizal_transfer;
  out.uptake_vector.ectomycorrhizal_uptake[days_gone + day]   = uptake.ectomycorrhizal_uptake;
  out.uptake_vector.root_uptake[days_gone + day]              = uptake.root_uptake;
  out.uptake_vector.total_uptake[days_gone + day]             = uptake.total_uptake;

  // Misc
  out.nitrogen_balance[days_gone + day] = nitrogen_balance;
  out.tree_alive[days_gone + day]       = tree_alive;
}


/*
 * Photosynthesis
 */

void log_photosynthesis(int day,
                        int days_gone,
                        const photosynthesis_out& input,
                        output_vector& out) {

  out.photosynthesis.GPP[days_gone + day]       = input.GPP;
  out.photosynthesis.ET[days_gone + day]        = input.ET;
  out.photosynthesis.SoilWater[days_gone + day] = input.SoilWater;
}


