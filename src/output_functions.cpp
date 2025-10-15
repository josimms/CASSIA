#include "CASSIA.h"

/*
 * Output Functions
 */

// df1

Rcpp::DataFrame createGrowthDataFrame(const output_vector& out) {
  return Rcpp::DataFrame::create(
    Rcpp::_["year"] = out.year,
    Rcpp::_["day"] = out.day,
    Rcpp::_["g"] = out.g,
    Rcpp::_["bud_growth_potential"] = out.potential_bud,
    Rcpp::_["wall_growth_potential"] = out.potential_wall,
    Rcpp::_["needle_growth_potential"] = out.potential_needles,
    Rcpp::_["roots_growth_potential"] = out.potential_roots,
    Rcpp::_["height_growth_potential"] = out.potential_height,
    Rcpp::_["bud_growth"] = out.bud,
    Rcpp::_["diameter_growth"] = out.diameter,
    Rcpp::_["needle_growth"] = out.needles,
    Rcpp::_["root_growth"] = out.roots,
    Rcpp::_["height_growth"] = out.height,
    Rcpp::_["respiration_growth"] = out.respiration_output.growth,
    Rcpp::_["respiration_maintenance"] = out.respiration_output.maintenance,
    Rcpp::_["ring_width"] = out.ring_width,
    Rcpp::_["mm_potenital"] = out.potential_ring_width
  );
}

// df2

Rcpp::DataFrame createSugarDataFrame(const output_vector& out) {
  return Rcpp::DataFrame::create(
    Rcpp::_["sugar"] = out.sugar,
    Rcpp::_["starch"] = out.starch,
    Rcpp::_["storage"] = out.storage_term_vector.needles, // adjust if total needed
    Rcpp::_["starch_needles"] = out.starch_vector.needles,
    Rcpp::_["starch_phloem"] = out.starch_vector.phloem,
    Rcpp::_["starch_xylem_sh"] = out.starch_vector.xylem_sh,
    Rcpp::_["starch_xylem_st"] = out.starch_vector.xylem_st,
    Rcpp::_["starch_roots"] = out.starch_vector.roots,
    Rcpp::_["sugar_needles"] = out.sugar_vector.needles,
    Rcpp::_["sugar_phloem"] = out.sugar_vector.phloem,
    Rcpp::_["sugar_xylem_sh"] = out.sugar_vector.xylem_sh,
    Rcpp::_["sugar_xylem_st"] = out.sugar_vector.xylem_st,
    Rcpp::_["sugar_roots"] = out.sugar_vector.roots,
    Rcpp::_["sugar_to_mycorrhiza"] = out.sugar_vector.to_mycorrhiza,
    Rcpp::_["sugar_out_surplus"] = out.sugar_vector.surplus,
    Rcpp::_["storage_term_needles"] = out.storage_term_vector.needles,
    Rcpp::_["storage_term_phloem"] = out.storage_term_vector.phloem,
    Rcpp::_["storage_term_xylem_sh"] = out.storage_term_vector.xylem_sh,
    Rcpp::_["storage_term_xylem_st"] = out.storage_term_vector.xylem_st,
    Rcpp::_["storage_term_roots"] = out.storage_term_vector.roots,
    Rcpp::_["n_E_pot"] = out.nitrogen_capacity_vector.height,   // adjust if needed
    Rcpp::_["n_W_pot"] = out.nitrogen_capacity_vector.diameter,
    Rcpp::_["n_M_pot"] = out.nitrogen_capacity_vector.mycorrhiza,
    Rcpp::_["nitrogen_balance"] = out.nitrogen_balance,
    Rcpp::_["nitrogen_capacity_needles"] = out.nitrogen_capacity_vector.needles,
    Rcpp::_["nitrogen_capacity_wall"] = out.nitrogen_capacity_vector.diameter,
    Rcpp::_["nitrogen_capacity_height"] = out.nitrogen_capacity_vector.height,
    Rcpp::_["nitrogen_capacity_bud"] = out.nitrogen_capacity_vector.bud,
    Rcpp::_["nitrogen_capacity_roots"] = out.nitrogen_capacity_vector.roots
  );
}

// df3

Rcpp::DataFrame createPrelesDataFrame(const output_vector& out) {
  return Rcpp::DataFrame::create(
    Rcpp::_["GPP"] = out.photosynthesis.GPP,
    Rcpp::_["ET"] = out.photosynthesis.ET,
    Rcpp::_["SoilWater"] = out.photosynthesis.SoilWater,
    Rcpp::_["fAPAR"] = out.photosynthesis.fAPAR
  );
}

// df4

Rcpp::DataFrame createCulmGrowthDataFrame(const output_vector& out) {
  return Rcpp::DataFrame::create(
    Rcpp::_["culm_growth_height"] = out.culm_growth.height,
    Rcpp::_["culm_growth_needles"] = out.culm_growth.needles,
    Rcpp::_["culm_growth_diameter"] = out.culm_growth.diameter,
    Rcpp::_["culm_growth_diameter_potential"] = out.culm_growth.diameter_potential,
    Rcpp::_["culm_growth_roots"] = out.culm_growth.roots,
    Rcpp::_["culm_growth_xylem_sh"] = out.culm_growth.xylem_sh,
    Rcpp::_["culm_growth_xylem_st"] = out.culm_growth.xylem_st,
    Rcpp::_["culm_growth_phloem"] = out.culm_growth.phloem,
    Rcpp::_["tree_alive"] = out.tree_alive,
    Rcpp::_["LAI"] = out.LAI,
    Rcpp::_["nitrogen_balance"] = out.nitrogen_balance,
    Rcpp::_["nitrogen_capacity_needles"] = out.nitrogen_capacity_vector.needles,
    Rcpp::_["nitrogen_capacity_wall"] = out.nitrogen_capacity_vector.diameter,
    Rcpp::_["nitrogen_capacity_height"] = out.nitrogen_capacity_vector.height,
    Rcpp::_["nitrogen_capacity_bud"] = out.nitrogen_capacity_vector.bud,
    Rcpp::_["nitrogen_capacity_roots"] = out.nitrogen_capacity_vector.roots
  );
}


/*
 * Tests
 */

void printOutputVectorSizes(const output_vector& out) {


  std::cout << "=== Growth DataFrame ===" << std::endl;
  std::cout << "year: " << out.year.size() << std::endl;
  std::cout << "day: " << out.day.size() << std::endl;
  std::cout << "g: " << out.g.size() << std::endl;
  std::cout << "bud_growth_potential: " << out.potential_bud.size() << std::endl;
  std::cout << "wall_growth_potential: " << out.potential_wall.size() << std::endl;
  std::cout << "needle_growth_potential: " << out.potential_needles.size() << std::endl;
  std::cout << "roots_growth_potential: " << out.potential_roots.size() << std::endl;
  std::cout << "height_growth_potential: " << out.potential_height.size() << std::endl;
  std::cout << "bud_growth: " << out.bud.size() << std::endl;
  std::cout << "diameter_growth: " << out.diameter.size() << std::endl;
  std::cout << "needle_growth: " << out.needles.size() << std::endl;
  std::cout << "root_growth: " << out.roots.size() << std::endl;
  std::cout << "height_growth: " << out.height.size() << std::endl;
  std::cout << "respiration_growth: " << out.respiration_output.growth.size() << std::endl;
  std::cout << "respiration_maintenance: " << out.respiration_output.maintenance.size() << std::endl;
  std::cout << "ring_width: " << out.ring_width.size() << std::endl;
  std::cout << "mm_potenital: " << out.potential_ring_width.size() << std::endl;

  std::cout << "\n=== Sugar DataFrame ===" << std::endl;
  std::cout << "sugar: " << out.sugar.size() << std::endl;
  std::cout << "starch: " << out.starch.size() << std::endl;
  std::cout << "storage_term_vector.needles: " << out.storage_term_vector.needles.size() << std::endl;
  std::cout << "starch_vector.needles: " << out.starch_vector.needles.size() << std::endl;
  std::cout << "starch_vector.phloem: " << out.starch_vector.phloem.size() << std::endl;
  std::cout << "starch_vector.xylem_sh: " << out.starch_vector.xylem_sh.size() << std::endl;
  std::cout << "starch_vector.xylem_st: " << out.starch_vector.xylem_st.size() << std::endl;
  std::cout << "starch_vector.roots: " << out.starch_vector.roots.size() << std::endl;
  std::cout << "sugar_vector.needles: " << out.sugar_vector.needles.size() << std::endl;
  std::cout << "sugar_vector.phloem: " << out.sugar_vector.phloem.size() << std::endl;
  std::cout << "sugar_vector.xylem_sh: " << out.sugar_vector.xylem_sh.size() << std::endl;
  std::cout << "sugar_vector.xylem_st: " << out.sugar_vector.xylem_st.size() << std::endl;
  std::cout << "sugar_vector.roots: " << out.sugar_vector.roots.size() << std::endl;
  std::cout << "sugar_vector.to_mycorrhiza: " << out.sugar_vector.to_mycorrhiza.size() << std::endl;
  std::cout << "sugar_vector.surplus: " << out.sugar_vector.surplus.size() << std::endl;
  std::cout << "storage_term_vector.phloem: " << out.storage_term_vector.phloem.size() << std::endl;
  std::cout << "storage_term_vector.xylem_sh: " << out.storage_term_vector.xylem_sh.size() << std::endl;
  std::cout << "storage_term_vector.xylem_st: " << out.storage_term_vector.xylem_st.size() << std::endl;
  std::cout << "storage_term_vector.roots: " << out.storage_term_vector.roots.size() << std::endl;
  std::cout << "nitrogen_balance: " << out.nitrogen_balance.size() << std::endl;
  // TODO: n_E_pot isn't the right output
  std::cout << "n_E_pot (height): " << out.nitrogen_capacity_vector.height.size() << std::endl;
  std::cout << "n_W_pot (diameter): " << out.nitrogen_capacity_vector.diameter.size() << std::endl;
  std::cout << "n_M_pot (mycorrhiza): " << out.nitrogen_capacity_vector.mycorrhiza.size() << std::endl;
  std::cout << "nitrogen_capacity_needles: " << out.nitrogen_capacity_vector.needles.size() << std::endl;
  std::cout << "nitrogen_capacity_bud: " << out.nitrogen_capacity_vector.bud.size() << std::endl;
  std::cout << "nitrogen_capacity_roots: " << out.nitrogen_capacity_vector.roots.size() << std::endl;

  std::cout << "\n=== Preles DataFrame ===" << std::endl;
  std::cout << "GPP: " << out.photosynthesis.GPP.size() << std::endl;
  std::cout << "ET: " << out.photosynthesis.ET.size() << std::endl;
  std::cout << "SoilWater: " << out.photosynthesis.SoilWater.size() << std::endl;
  std::cout << "fAPAR: " << out.photosynthesis.fAPAR.size() << std::endl;

  std::cout << "\n=== Culm Growth DataFrame ===" << std::endl;
  std::cout << "culm_growth_height: " << out.culm_growth.height.size() << std::endl;
  std::cout << "culm_growth_needles: " << out.culm_growth.needles.size() << std::endl;
  std::cout << "culm_growth_diameter: " << out.culm_growth.diameter.size() << std::endl;
  std::cout << "culm_growth_diameter_potential: " << out.culm_growth.diameter_potential.size() << std::endl;
  std::cout << "culm_growth_roots: " << out.culm_growth.roots.size() << std::endl;
  std::cout << "tree_alive: " << out.tree_alive.size() << std::endl;
  std::cout << "LAI: " << out.LAI.size() << std::endl;
  std::cout << "nitrogen_balance: " << out.nitrogen_balance.size() << std::endl;
  std::cout << "nitrogen_capacity_needles: " << out.nitrogen_capacity_vector.needles.size() << std::endl;
  std::cout << "nitrogen_capacity_wall: " << out.nitrogen_capacity_vector.diameter.size() << std::endl;
  std::cout << "nitrogen_capacity_height: " << out.nitrogen_capacity_vector.height.size() << std::endl;
  std::cout << "nitrogen_capacity_bud: " << out.nitrogen_capacity_vector.bud.size() << std::endl;
  std::cout << "nitrogen_capacity_roots: " << out.nitrogen_capacity_vector.roots.size() << std::endl;

  std::cout << "\n===========================\n" << std::endl;
}

