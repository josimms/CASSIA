#include "CASSIA.h"

xylogenesis_out xylogenesis(int no_day,
                            int day,
                            CASSIA_parameters parameters,
                            CASSIA_common common,
                            bool xylogenesis_option,
                            bool environmental_effect_xylogenesis,
                            double TAir,

                            double n_rows,
                            double max_ew_cells,

                            double n_E_pot_old,
                            double n_W_pot_old,
                            double n_M_pot_old,

                            double g,
                            std::vector<double> en_growth_vector,

                            // Iterations
                            //double wall_tot,
                            //double tau_E,
                            //double tau_W,
                            //double n_E_old,
                            //double n_W_old,
                            //double n_M_old,
                            double tau_W_old,

                            // Parameters
                            double carbon_daily_rate_ew, // TODO: where does this number come from
                            double carbon_daily_rate_lw)
{

  // TODO: should these be a iteration or a parameter?
  double CE_ew, CE_lw, DE_ew, DE_lw;
  double cell_l, M_suc, GD;

  double sD; // TODO: what is this?
  double tau_E, tau_W;
  double n_E, n_W, n_M;

  double tot_cells, fD, LD, storage_reduction, g_sD_T;

  if (xylogenesis_option) {
    GD = g_sD_T * fD * LD	* storage_reduction; // Daily potential number of new cells per day (in one radial cell row), used to be called division
    // tot_cells = cumsum(GD);

    if (environmental_effect_xylogenesis) {
      double temperature_effect_on_wall_formation = (g + 1)/1.5;

      // Assume that we are considering the total potential xylogenesis
      //if (day > last_dameter_growth_day) {
      // TODO: need to check the difference between cells and days
      //}

      // TODO: this requires a for loop, which needs the total cells, which is in the predicted code


    } else {

      // "Uggla" divides cells to early and late wood
      tau_E = 0;
      tau_W = 0;
      if (sD < pow(parameters.sDc, 2) / parameters.Uggla) {
        tau_E = parameters.tau_Ee;
        tau_W = parameters.tau_We;
      } else if (sD >= pow(parameters.sDc, 2) / parameters.Uggla) {
        tau_E = parameters.tau_El;
        tau_W = parameters.tau_Wl;
      }

      // Number of cells in enlargement, wall formation and mature phases on each day
      double n_E_pot = 0;
      double n_W_pot = 0;
      double n_M_pot = 0;
      if (day > 0) {
        n_E_pot = n_E_pot_old + GD - n_E_pot_old / tau_E;
        n_W_pot = n_W_pot_old + n_E_pot_old / tau_E - n_W_pot_old / tau_W_old;
        n_M_pot = n_M_pot_old + n_W_pot_old / tau_W_old;
      }

      // carbon to enlargement of one earlywood/latewood cell per one day and cell (kg C cell-1 day-1)
      double CE_ew, CE_lw;
      CE_ew = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_ew / 2, 2) * cell_l * M_suc / (common.gas_const * (TAir + common.abs_zero) * 1000) * 12 * common.M_C / M_suc / tau_E;
      CE_lw = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_ew / 2, 2) * cell_l * M_suc / (common.gas_const * (TAir + common.abs_zero) * 1000) * 12 * common.M_C / M_suc / tau_E;

      // carbon.wall <- en.growth <- en.release <- wall.growth <- d.growth <- NULL

      // The carbon needed to enlarge the cells  (kg C /day) (in one radial cell row)
      double carbon_enlargement = 0;
      double ew_cells;
      if (sD < pow(parameters.sDc, 2) / parameters.Uggla) {
        carbon_enlargement = CE_ew * n_E;
        ew_cells = tot_cells; // TODO: tot_cells, should the ew_max value be updated?
      } else if (sD >= pow(parameters.sDc, 2) / parameters.Uggla) {
        carbon_enlargement = CE_lw * n_E;
        ew_cells = max_ew_cells;
      }

      // Carbon to enlargement per day (kg C /day) (in the whole tree)
      double en_growth = n_rows * carbon_enlargement;

      double en_release = 0;
      if (day < parameters.tau_Ee) {
        en_release = 0; // the carbon used in enlargement is released after some days.
      } else if (day >= parameters.tau_Ee) {
        en_release = en_growth_vector[day - round(parameters.tau_Ee)]; // TODO: not sure how to index this
      }

      // Carbon to wall formation
      // Carbon.daily.rate determined in parameters_common.R but NOTE!!!! not used at the moment, replaced by a parameter set to result in density app. 200 kg C m-3!

      /*
       double CW;
       if (sD < pow(parameters.sDc, 2) / parameters.Uggla) {
       CW = carbon_daily_rate_ew;
       } else if (sD >= pow(parameters.sDc, 2) / parameters.Uggla) {
       CW = carbon_daily_rate_lw;
       }
       */

      double CW = 0.000000000029;

      // The use of carbon to wall growth kg C per day
      double wall_growth = n_rows * CW * n_W;
      double lw_cells = tot_cells - ew_cells;

      double ew_cells_max, cell_density_ew, cell_density_lw; // TODO: how do these parameters work?
      double max_lw_cells, max_ew_cells;

      if (day == no_day) {
        // Calculation of the daily annual ring width, starting from 16.12.2013.
        // Note! This does not include the size of enlarging cells, only wall forming and mature
        double diameter_cells = n_M + n_W;
        double diameter_ew_cells, cell_wall_thickness, cell_d_final, cell_density;
        if (diameter_cells <= ew_cells_max) {
          diameter_ew_cells = diameter_cells;
          cell_wall_thickness = parameters.wall_thickness_ew;
          cell_d_final = parameters.cell_d_ew;
          cell_density = cell_density_ew;
        } else {
          diameter_ew_cells = max_ew_cells;
          cell_wall_thickness = parameters.wall_thickness_lw;
          cell_d_final = parameters.cell_d_lw;
          cell_density = cell_density_lw;
        }
        double diameter_lw_cells = diameter_cells - diameter_ew_cells;

        double ew_width = diameter_ew_cells * parameters.cell_d_ew * 1000;   // + ifelse(sD < sDc^2 / Uggla, n.E * parameters[c("cell.d.ew"), c(site)] / 2)
        double lw_width = diameter_lw_cells * parameters.cell_d_ew * 1000; // + ifelse(sD >= sDc^2 / Uggla, n.E * parameters[c("cell.d.ew"), c(site)] / 2)
        double ring_width = ew_width + lw_width;

        // double ring_density = cumsum(CW * n_W)[365] / (ew_cells[365] * pow(parameters.cell_d_ew, 2) * parameters.cell_l_ew + lw_cells[365] * pow(parameters.cell_d_ew, 2) * parameters.cell_l_lw;

        double ew_cells = max_ew_cells;
        double lw_cells = max_lw_cells;
      }
    }

    double en_growth, en_release, wall_growth;
    double en_pot_growth = en_growth;
    double en_pot_release = en_release;
    double wall_pot_growth = wall_growth;

  } else { // Else from xylogenesis

    tau_W = 0;
    tau_E = 0;
    if (sD < pow(parameters.sDc, 2) / parameters.Uggla) {
      tau_E = parameters.tau_Ee;
      tau_W = parameters.tau_We;
    } else if (sD >= pow(parameters.sDc, 2) / parameters.Uggla) {
      tau_E = parameters.tau_El;
      tau_W = parameters.tau_Wl;
    }

    double n_E_pot, n_W_pot, n_M_pot;
    if (day > 1) {
      n_E_pot = n_E_pot_old + GD - n_E_pot_old / tau_E;
      n_W_pot = n_W_pot_old + n_E_pot_old / tau_E - n_W_pot_old / tau_W_old;
      n_M_pot = n_M_pot_old + n_W_pot_old / tau_W_old;
    }

    double M_C;
    double CE_ew = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_ew / 2, 2) * cell_l * M_suc / common.gas_const * (TAir + common.abs_zero * 1000) * 12 * M_C / M_suc / tau_E;
    double CE_lw = common.osmotic_sugar_conc * M_PI * pow(parameters.cell_d_ew / 2, 2) * cell_l * M_suc / common.gas_const * (TAir + common.abs_zero * 1000) * 12 * M_C / M_suc / tau_E;

    double carbon_enlargement_pot = 0;
    if (sD < pow(parameters.sDc, 2) / parameters.Uggla) {
      carbon_enlargement_pot = CE_ew * n_E_pot;
    } else if (sD >= pow(parameters.sDc, 2) / parameters.Uggla) {
      carbon_enlargement_pot = CE_lw * n_E_pot;
    }

    double en_pot_growth = n_rows * carbon_enlargement_pot;
    double en_pot_release = 0;
    if (day < parameters.tau_Ee) {
      en_pot_release = 0; // the carbon used in enlargement is released after some days.
    } else if (day >= parameters.tau_Ee) {
      en_pot_release = en_growth_vector[day - round(parameters.tau_Ee)];
    }

    double CW = 0.000000000029;

    double DW;
    double wall_pot_growth = n_rows * DW * n_W_pot;

    double cells_pot = n_W_pot + n_M_pot;

    double ew_cells_pot = 0;
    double lw_cells_pot = 0;
    double max_cells_pot, max_ew_width_pot;
    if (sD < pow(parameters.sDc, 2) / parameters.Uggla) {
      ew_cells_pot = cells_pot;
      lw_cells_pot = cells_pot - ew_cells_pot - max_cells_pot;
    }
    double ew_width_pot = ew_cells_pot * parameters.cell_d_ew * 1000;
    double lw_width_pot = lw_cells_pot * parameters.cell_d_lw * 1000;
    double pot_mm = 0;
    if (sD < pow(parameters.sDc, 2) / parameters.Uggla) {
      pot_mm = ew_width_pot;
    } else {
      pot_mm = lw_width_pot + max_ew_width_pot;
    }


  } // End of the xylogenesis

  xylogenesis_out out;
  out.n_E = n_E;
  out.n_E_pot = 0.0;
  out.n_W_pot = 0.0;
  out.n_M_pot = 0.0;

  return out;
}

/*
 * Xylogenesis Wrapper
 */

// [[Rcpp::export]]
Rcpp::List xylogenesis_wrapper(int no_day,
                               int day,
                               Rcpp::DataFrame pCASSIA_parameters,
                               Rcpp::DataFrame pCASSIA_common,
                               Rcpp::DataFrame pCASSIA_sperling,
                               std::vector<double> extras_sperling,
                               bool xylogenesis_option,
                               bool environmental_effect_xylogenesis,
                               double TAir,
                               double n_rows,
                               double max_ew_cells,
                               double n_E_pot_old,
                               double n_W_pot_old,
                               double n_M_pot_old,
                               double g,
                               std::vector<double> en_growth_vector,
                               double tau_W_old,
                               double carbon_daily_rate_ew,
                               double carbon_daily_rate_lw) {

  // TODO: where does this number come from

  CASSIA_common common = make_common(pCASSIA_common);
  CASSIA_parameters parameters = make_CASSIA_parameters(pCASSIA_parameters, pCASSIA_sperling);

  xylogenesis_out out = xylogenesis(no_day,
                                    day,
                                    parameters,
                                    common,
                                    xylogenesis_option,
                                    environmental_effect_xylogenesis,
                                    TAir,

                                    n_rows,
                                    max_ew_cells,

                                    n_E_pot_old,
                                    n_W_pot_old,
                                    n_M_pot_old,

                                    g,
                                    en_growth_vector,

                                    // Iterations
                                    //double wall_tot,
                                    //double tau_E,
                                    //double tau_W,
                                    //double n_E_old,
                                    //double n_W_old,
                                    //double n_M_old,
                                    tau_W_old,

                                    // Parameters
                                    carbon_daily_rate_ew, // TODO: where does this number come from
                                    carbon_daily_rate_lw);

  return Rcpp::List::create(Rcpp::_["n_E"] = out.n_E,
                            Rcpp::_["n_E_pot"] = out.n_E_pot,
                            Rcpp::_["n_W_pot"] = out.n_W_pot,
                            Rcpp::_["n_M_pot"] = out.n_M_pot);
}
