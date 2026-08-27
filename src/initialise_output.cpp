#include "CASSIA.h"

void initialize_output_vector(output_vector& out, int simulation_time) {
  // Years and Days
  out.year.resize(simulation_time);
  out.day.resize(simulation_time);

  // Potential growth
  out.potential_diameter.resize(simulation_time);
  out.potential_needles.resize(simulation_time);
  out.potential_height.resize(simulation_time);
  out.potential_wall.resize(simulation_time);
  out.potential_roots.resize(simulation_time);
  out.potential_bud.resize(simulation_time);
  out.potential_ecto.resize(simulation_time);
  out.potential_ring_width.resize(simulation_time);

  // Actual growth
  out.g.resize(simulation_time);
  out.diameter.resize(simulation_time);
  out.needles.resize(simulation_time);
  out.height.resize(simulation_time);
  out.wall.resize(simulation_time);
  out.roots.resize(simulation_time);
  out.ecto.resize(simulation_time);
  out.bud.resize(simulation_time);
  out.use.resize(simulation_time);
  out.release.resize(simulation_time);
  out.ring_width.resize(simulation_time);
  out.en_pot_growth.resize(simulation_time);

  // Respiration
  out.RmN.resize(simulation_time);
  out.RmS.resize(simulation_time);
  out.RmR.resize(simulation_time);
  out.Rm_a.resize(simulation_time);

  // Repola
  out.needle_mass.resize(simulation_time);
  out.m_N.resize(simulation_time);
  out.m_N_tot.resize(simulation_time);
  out.m_R_tot.resize(simulation_time);

  out.fAPAR.resize(simulation_time);
  out.LAI.resize(simulation_time);

  // Tree status
  out.tree_alive.resize(simulation_time, true); // initialize to true
  out.sugar.resize(simulation_time);
  out.starch.resize(simulation_time);
  out.nitrogen_balance.resize(simulation_time);

  // Nested structures (if needed, they must also be initialized similarly)
  out.culm_growth.resize(simulation_time);
  out.respiration_output.resize(simulation_time);
  out.sugar_vector.resize(simulation_time);
  out.starch_vector.resize(simulation_time);
  out.storage_term_vector.resize(simulation_time);
  out.nitrogen_capacity_vector.resize(simulation_time);
  out.uptake_vector.resize(simulation_time);
  out.photosynthesis.resize(simulation_time);
}
