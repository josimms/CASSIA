#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

struct carbo_tracker{
  std::vector<double> needles; // Inital value from inputs
  std::vector<double> phloem; // Inital value from inputs
  std::vector<double> xylem_sh; // Inital value from inputs
  std::vector<double> xylem_st; // Inital value from inputs
  std::vector<double> roots; // Inital value from inputs
  std::vector<double> mycorrhiza; // Initial value 0
  double initial_amount; // Inital value from inputs
  double B;
};

struct carbo_balance{
  carbo_tracker sugar;
  carbo_tracker statch;
};

struct conc_gradient{
  std::vector<double> needles_to_phooem;
  std::vector<double> phloem_to_xylem_sh;
  std::vector<double> phloem_to_xylem_st;
  std::vector<double> phloem_to_roots;
  std::vector<double> roots_to_myco;
  double ratio_needles_to_phloem;
  double ratio_phloem_to_xylem_sh;
  double ratio_phloem_to_xylem_st;
  double ratio_phloem_to_roots;
};

struct interaction {
  double needles_to_phloem;
  double phloem_to_xylem_sh;
  double phloem_to_xylem_st;
  double phloem_to_roots;
};

struct parameter_sperling {
  double needles;
  double phloem;
  double xylem_sh;
  double xylem_st;
  double roots;
};

struct common_str {
  double Rg_N;
  double Rg_S;
  double Rg_R;
};

struct organ_values {
  std::vector<double> needles;
  std::vector<double> height;
  std::vector<double> wall;
  std::vector<double> roots;
  std::vector<double> use;
  std::vector<double> release;
};

parameter_sperling storage_carbohydrate(parameter_sperling alfa, parameter_sperling critical_W, parameter_sperling Wala) {
  parameter_sperling out;

  out.needles = 1 / (1 - 1/std::exp(alfa.needles) * (critical_W.needles - Wala.needles)); // TODO: what is exp in c // TODO: check what the exp should cover!
  out.phloem = 1 / (1 - 1/std::exp(alfa.phloem) * (critical_W.phloem - Wala.phloem));
  out.xylem_sh = 1 / (1 - 1/std::exp(alfa.xylem_sh) * (critical_W.xylem_sh - Wala.xylem_sh));
  out.xylem_st = 1 / (1 - 1/std::exp(alfa.xylem_st) * (critical_W.xylem_st - Wala.xylem_st));
  out.roots = 1 / (1 - 1/std::exp(alfa.roots) * (critical_W.roots - Wala.roots));

  return out;
}

carbo_tracker As_initiliser(carbo_tracker As, carbo_tracker Ad, double equilibrium_temperature, double Bd, double Bs, double mycorhizza_threshold)
{
  As.needles.push_back(Ad.needles[1] * std::exp(equilibrium_temperature * (Bd - Bs)));
  As.phloem.push_back(Ad.phloem[1] * std::exp(equilibrium_temperature * (Bd - Bs)));
  As.xylem_sh.push_back(Ad.xylem_sh[1] * std::exp(equilibrium_temperature * (Bd - Bs)));
  As.xylem_st.push_back(Ad.xylem_st[1] * std::exp(equilibrium_temperature * (Bd - Bs)));
  As.roots.push_back(Ad.roots[1] * std::exp(equilibrium_temperature * (Bd - Bs)));

  return As;
}

double storage_update(double ak, double alfa, double sugar, double starch, double Wala) {
  double comparison = std::min(1.0, ak * (1 - 1/std::exp(alfa * (sugar + starch - Wala))));
  double out = std::max(0.0, comparison);
  return out;
}

double emergancy(double sugar, double starch, double tau_emergancy, double lower_bound) {
  double out;
  if (sugar < lower_bound) {
    double comaprison = std::max((lower_bound - sugar) / tau_emergancy, 0.0);
    out = std::min(starch, comaprison);
  } else {
    out = 0;
  }
  return out;
}

/*
 * Some of the parameters are site dependent, could define the site side in the main function code so I just import the right variables into this code
 *
 * Needle mass per year, should deinfe the year before this function!
 *
 */


// [[Rcpp::export]]
carbo_balance sugar_model(int ndays,
                          // Weather data
                          std::vector<double> Temp,
                          std::vector<double> PF,

                          // Vectors from CASSIA model
                          double HN0,
                          double D0,
                          double D00,
                          double LR0,
                          std::vector<double> RmN,
                          std::vector<double> RmS,
                          std::vector<double> RmR,

                          organ_values pot_growth,

                          // Initial conditions
                          bool storage_grows,
                          double mycorhiza_threshold,
                          parameter_sperling Ad0,

                          // Parameters
                          double Q10d, double Q10s,
                          double carbon_sugar,
                          double needle_mass,
                          common_str common,
                          interaction resistance) {
  /*
   * MAKE THE PARAMETERS AND VARIABLES
   */

  // Make the structures

  carbo_tracker sugar;
  carbo_tracker starch;

  carbo_tracker storage_term;
  carbo_tracker Ad;
  carbo_tracker As;
  carbo_tracker Kd;
  carbo_tracker Ks;
  conc_gradient concentration_gradient;

  // TODO Work out if I need these considering the normal model or whether they should all be in the storage grows section!

  // Storage

  parameter_sperling critical_amount_carbo; // W_0 input
  parameter_sperling lower_bound_W; // Wala, input
  parameter_sperling critical_W; // W.crit.needles, input
  parameter_sperling alfa; // alfa, TODO: input?
  parameter_sperling ak = storage_carbohydrate(alfa, critical_W, lower_bound_W);
  parameter_sperling tau_emergancy; // Bayesian initiation!
  parameter_sperling emergancy_transfer;
  parameter_sperling lambda; // should be an input
  parameter_sperling delta; // should be an input

  /*
   * INITIALISE THE STRUCTURES
   */

  if (storage_grows) { // TODO: define or work out how to input these!
    // TODO: deine if the storage doesn't grow!

    lower_bound_W.needles = HN0 / D00 * lower_bound_W.needles;
    lower_bound_W.phloem =  D0 / D00 * lower_bound_W.phloem;
    lower_bound_W.roots =  LR0 / D00 * lower_bound_W.roots;
    lower_bound_W.xylem_sh =  LR0 / D00 * lower_bound_W.xylem_sh;
    lower_bound_W.xylem_st =  LR0 / D00 * lower_bound_W.xylem_st;

    critical_W.needles = HN0 / D00 * critical_W.needles;
    critical_W.phloem =  D0 / D00 * critical_W.phloem;
    critical_W.roots =  LR0 / D00 * critical_W.roots;
    critical_W.xylem_sh =  LR0 / D00 * critical_W.xylem_sh;
    critical_W.xylem_st =  LR0 / D00 * critical_W.xylem_st;
  }

  // Compute initial Te by the mean temperature for the first week of # October plus 3C (for the exponential nature of the curves)
  double equilibrium_temperature = std::accumulate(Temp.begin() + 273, Temp.begin() + 280, 0.0) / Temp.size() + 3; // 	# Compute initial Te by the mean temperature for the first week of # October plus 3C (for the exponential nature of the curves)

  // Initialise the sperling things
  Ad.needles.push_back(Ad0.needles);
  Ad.phloem.push_back(Ad0.phloem);
  Ad.xylem_sh.push_back(Ad0.xylem_sh);
  Ad.xylem_st.push_back(Ad0.xylem_st);
  Ad.roots.push_back(Ad0.roots);

  sugar.B = log(Q10s)/10;
  starch.B = log(Q10d)/10;

  As = As_initiliser(As, Ad, equilibrium_temperature, starch.B, sugar.B, mycorhiza_threshold); // TODO: does this work as well as it as a iterative call!

  Kd.needles.push_back(Ad.needles[1]*exp(starch.B*Temp[1]));
  Ks.needles.push_back(As.needles[1]*exp(sugar.B*Temp[1])); // TODO: for all other organs as well!

  sugar.needles.push_back(0); // Input, TODO: input this!
  sugar.phloem.push_back(0);
  sugar.xylem_sh.push_back(0);
  sugar.xylem_st.push_back(0);
  sugar.roots.push_back(0);
  sugar.mycorrhiza.push_back(0);

  starch.needles.push_back(0); // Input, TODO: input this!
  starch.phloem.push_back(0);
  starch.xylem_sh.push_back(0);
  starch.xylem_st.push_back(0);
  starch.roots.push_back(0);

  storage_term.needles.push_back(1); // TODO: should this be initialised with the equation?
  // TODO: do this for the rest of the organs

  concentration_gradient.needles_to_phooem.push_back(sugar.needles[1]+starch.needles[1] - concentration_gradient.ratio_needles_to_phloem*(sugar.phloem[1]+starch.phloem[1]));
  concentration_gradient.phloem_to_roots.push_back(sugar.phloem[1]+starch.phloem[1] - concentration_gradient.ratio_phloem_to_roots*(sugar.roots[1]+starch.roots[1]));
  concentration_gradient.phloem_to_xylem_sh.push_back(sugar.phloem[1]+starch.phloem[1] - concentration_gradient.ratio_phloem_to_xylem_sh*(sugar.xylem_sh[1]+starch.xylem_sh[1]));
  concentration_gradient.phloem_to_xylem_st.push_back(sugar.phloem[1]+starch.phloem[1] - concentration_gradient.ratio_phloem_to_xylem_st*(sugar.xylem_st[1]+starch.xylem_st[1]));
  concentration_gradient.roots_to_myco.push_back(sugar.phloem[1]+starch.phloem[1] - mycorhizza_threshold);
    //concentration_gradient.ratio_roots_to_myco*(sugar.myco[1]+starch.myco[1])));

  /*
   * DAILY LOOP!
   */

  for (int day = 1; day <= ndays; ++day) // Remember! Starts from day 2, indexing starts from 0!
  { // Wala and alfa are baised of older definitions, should I go back and correct this?

    /*
     * PARAMETER UPDATES!
     */

    storage_term.needles.push_back(storage_update(ak.needles, alfa.needles, sugar.needles[day-1], starch.needles[day-1], lower_bound_W.needles));
    storage_term.phloem.push_back(storage_update(ak.phloem, alfa.phloem, sugar.phloem[day-1], starch.phloem[day-1], lower_bound_W.phloem));
    storage_term.roots.push_back(storage_update(ak.roots, alfa.roots, sugar.roots[day-1], starch.roots[day-1], lower_bound_W.roots));
    storage_term.xylem_sh.push_back(storage_update(ak.xylem_sh, alfa.xylem_sh, sugar.xylem_sh[day-1], starch.xylem_sh[day-1], lower_bound_W.xylem_sh));
    storage_term.xylem_st.push_back(storage_update(ak.xylem_st, alfa.xylem_st, sugar.xylem_st[day-1], starch.xylem_st[day-1], lower_bound_W.xylem_st));

    Kd.needles.push_back(Ad.needles[day-1]*std::exp(starch.B*Temp[day])); // TODO: is this initalised correctly when thinking of the code?
    Ks.needles.push_back(As.needles[day-1]*std::exp(sugar.B*Temp[day])); // TODO: for all other organs as well!

    /*
     * if there is a surplus goes to next organ down - concentration driven model
     * The differences are normalised by a multiplier which represents the average difference in magnitude between the two stores
     * otherwise all of the sugar would just immediately go to the roots
     *
     * NOTE: the forces and therefore sugar transfered are worked out for all organs based on the amount of sugar there in the beginning
     * this could lead to a slight error, but should be corrected by the starch latter just have to imagine that all of the sugar
     * goes to the allocated organs simultaneously
     */

    concentration_gradient.needles_to_phloem.push_back(sugar.needles[day-1]+starch.needles[day-1] - concentration_gradient.ratio_needles_to_phloem*(sugar.phloem[day-1]+starch.phloem[day-1]));
    concentration_gradient.phloem_to_roots.push_back(sugar.phloem[day-1]+starch.phloem[day-1] - concentration_gradient.ratio_phloem_to_roots*(sugar.roots[day-1]+starch.roots[day-1]));
    concentration_gradient.phloem_to_xylem_sh.push_back(sugar.phloem[day-1]+starch.phloem[day-1] - concentration_gradient.ratio_phloem_to_xylem_sh*(sugar.xylem_sh[day-1]+starch.xylem_sh[day-1]));
    concentration_gradient.phloem_to_xylem_st.push_back(sugar.phloem[day-1]+starch.phloem[day-1] - concentration_gradient.ratio_phloem_to_xylem_st*(sugar.xylem_st[day-1]+starch.xylem_st[day-1]));
    concentration_gradient.roots_to_myco.push_back(sugar.phloem[day-1]+starch.phloem[day-1] - mycorhizza_threshold);

    /*
     * SUGAR TRANSFER WITH ALL PROCESSES BUT EMERGANCY
     */

    // # Rm.a maintenance respiration separated into organs
    sugar.needles.push_back(sugar.needles[day-1] + PF[day] -          // Last day sugar + daily photosynthesis
      RmN[day] * storage_term.needles[day] -                          // - maintenance respiration (altered by the carbon storage)
      (1 + common.Rg_N) * storage_term_needles[day] * (pot_growth.needle[day] + pot_growth.bud[day]) -          // - growth and growth respiration altered by the storage
      pot_growth.use[day] + pot_growth.release[day] -                                                           // - growth sugar use and + release and to the rest of the organs
      resistance.needles_to_phloem * concentration_gradient.needles_to_phooem[day] +                            // - transfer between organs
      (Kd.needles[day] - Ks.needles[day]) * carbon_sugar * 0.001 * needle_mass);                                 // + sperling processes with links to the needle growth process


    // coefficients are from mass ratio in starch and sugar 2015 xls

    sugar.phloem.push_back(sugar.phloem[day-1] -
      0.082179938 * RmS[day] * storage_term.phloem[day] -
      0.082179938 * (1 + common.Rg_S) * storage_term.phloem[day] * (pot_growth.wall[day] + pot_growth.height[day]) +  // growth
      resistance.needles_to_phloem * concentration_gradient.needles_to_phooem[day] -                             // transfer between organs
      resistance.phloem_to_roots * concentration_gradient.phloem_to_roots[day] -                                 // transfer between organs
      resistance.phloem_to_xylem_sh * concentration_gradient.phloem_to_xylem_sh[day] -
      resistance.phloem_to_xylem_st * concentration_gradient.phloem_to_xylem_st[day] +
      (Kd.phloem[day] - Ks.phloem[day]) * carbon_sugar * 0.001 * 7.4);

    sugar.roots.push_back(sugar.roots[day-1] +
      resistance.phloem_to_roots * concentration_gradient.needles_to_phooem[day] -        // transfer between organs
      concentration_gradient.roots_to_myco[day] +                                         // transfer between organs, no multiplier as this is for mycorhiza and the model just takes the extra sugar
      (Kd.roots[day] - Ks.roots[day]) * carbon_sugar * 0.001 * 2.8 -
      (1 + common.Rg_R) * storage_term.roots[day] * pot_growth.roots[day] -               // growth
      RmR[day] * storage_term.roots[day]);                                                // maintenance respiration);


    sugar.xylem_sh.push_back(sugar.xylem_sh[day-1] -
      0.096020683 * RmS[day] * storage_term.xylem_sh[day] -                                   // maintenance respiration
      0.096020683 * (1 + common.Rg_S) * storage_term.xylem_sh[day] * (pot_growth.wall[day] + pot_growth.height[day]) +    // growth
      resistance.phloem_to_xylem_sh * concentration_gradient.phloem_to_xylem_sh[day] +
      (Kd.xylem_sh[i] - Ks.xylem_sh[i]) * carbon_sugar * 0.001 * 2.8);

    sugar.xylem_st.push_back(sugar.xylem_st[day-1] -
      0.821799379 * RmS[day] * storage_term.xylem_st[day] -                                 // maintenance respiration
      0.821799379 * (1 + common.Rg_S) * storage_term.xylem_st[day] * (pot_growth.wall[day] + pot_growth.height[day]) +  // growth
      resistance.phloem_to_xylem_st * concentration_gradient.phloem_to_xylem_st[day] +
      (Kd.xylem_st[day] - Ks.xylem_st[day]) * carbon_sugar * 0.001 * 2.8);

    sugar.mycorrhiza.push_back(concentration_gradient.roots_to_myco[day]);

    /*
     * STARCH UPDATED SPERLING
     */

    // SPERLING MODEL

    starch.needles.push_back(starch.needles[day-1] + (- Kd.needles[day] + Ks.needles[day]) * carbon_sugar * 0.001 * needle_mass); // Subtract starch degradation and add synthase to ST
    starch.phloem.push_back(starch.phloem[day-1] + (- Kd.phloem[day] + Ks.phloem[day]) * carbon_sugar * 0.001 * 7.4); // Subtract starch degradation and add synthase to ST
    starch.roots.push_back(starch.roots[i-1] + (- Kd.roots[i] + Ks.roots[i]) * carbon_sugar * 0.001 * 2.8); // Subtract starch degradation and add synthase to ST
    // TOOD: are the densities right here?
    starch.xylem_sh.push_back(starch.xylem_sh[day-1] + (- Kd.xylem_sh[day] + Ks.xylem_sh[day]) * carbon_sugar * 0.001 * 2.8); // Subtract starch degradation and add synthase to starch
    starch.xylem_st.push_back(starch.xylem_st[day-1] + (- Kd.xylem_st[day] + Ks.xylem_st[day]) * carbon_sugar * 0.001 * 2.8); // Subtract starch degradation and add synthase to starch

    /*
     * STARCH AND SUGAR UPDATED EMERGANCY MODEL
     *
     * If sugar is below a certain value starch is released so the sugar doesn't go negative before starch
     *  This is a proxy for a starch metabolism system, which seems to be present under stress in literature
     *  but I can't find a mechanism for scots pine
     *  values are below the lowest recorded value
     */

    // Work out energy transfer

    emergancy_transfer.needles =  emergancy(sugar.needles[day], starch.needles[day], tau_emergancy.needles[day], lower_bound_W.needles); //TODO: check if this lower bound is okay here or whether I should use the one in the R code
    emergancy_transfer.phloem = emergancy(sugar.phloem[day], starch.phloem[day], tau_emergancy.phloem[day], lower_bound_W.phloem);
    emergancy_transfer.xylem_sh = emergancy(sugar.xylem_sh[day], starch.xylem_sh[day], tau_emergancy.xylem_sh[day], lower_bound_W.xylem_sh);
    emergancy_transfer.xylem_st = emergancy(sugar.xylem_st[day], starch.xylem_st[day], tau_emergancy.xylem_st[day], lower_bound_W.xylem_st);
    emergancy_transfer.roots = emergancy(sugar.roots[day], starch.roots[day], tau_emergancy.roots[day], lower_bound_W.roots);

    // sugar update

    sugar.needles[day] = sugar.needles[day] + emergancy_transfer.needles;
    sugar.phloem[day] = sugar.phloem[day] + emergancy_transfer.phloem;
    sugar.xylem_sh[day] = sugar.xylem_sh[day] + emergancy_transfer.xylem_sh;
    sugar.xylem_st[day] = sugar.xylem_st[day] + emergancy_transfer.xylem_st;
    sugar.roots[day] = sugar.roots[day] + emergancy_transfer.roots;

    // starch update

    starch.needles[day] = starch.needles[day] - emergancy_transfer.needles;
    starch.phloem[day] = starch.phloem[day] - emergancy_transfer.phloem;
    starch.xylem_sh[day] = starch.xylem_sh[day] - emergancy_transfer.xylem_sh;
    starch.xylem_st[day] = starch.xylem_st[day] - emergancy_transfer.xylem_st;
    starch.roots[day] = starch.roots[day] - emergancy_transfer.roots;

    /*
     * SPERLING PARAMETER UPDATE FOR NEXT ITERATION
     */


    As.needles.push_back((1-lamda.needles)*As.needles[day-1]);
    As.phloem.push_back((1-lamda.phloem)*As.needles[day-1]);
    As.roots.push_back((1-lamda.roots)*As.needles[day-1]);
    As.xylem_sh.push_back((1-lamda.xylem_sh)*As.needles[day-1]);
    As.xylem_st.push_back((1-lamda.xylem_st)*As.needles[day-1]);

    Ad.needles.push_back((1-lamda.needles)*Ad.needles[day-1]);
    Ad.phloem.push_back((1-lamda.phloem)*Ad.needles[day-1]);
    Ad.roots.push_back((1-lamda.roots)*Ad.needles[day-1]);
    Ad.xylem_sh.push_back((1-lamda.xylem_sh)*Ad.needles[day-1]);
    Ad.xylem_st.push_back((1-lamda.xylem_st)*Ad.needles[day-1]);

    /*
     * Induce starch synthase if SC is high or degradation if it is low
     * These numbers are from september 2018
     * xylem not changed as no data to support it
     *
     * TODO: should I make these numbers automatic somehow?
     */

    if  (sugar.needles[day]>0.12)
    {
      As.needles[day]=As.needles[day]+delta.needles;
    }
    else (sugar.needles[day]<0.12 && starch.needles[day]> 0)
    {
      Ad.needles[day]=Ad.needles[day]+delta.needles;
    }
    if  (sugar.phloem[day]>0.28)
    {
      As.phloem[day]=As.phloem[day]+delta.phloem;
    }
    else (sugar.phloem[day]<0.28 && starch.phloem[day]> 0)
    {
      Ad.phloem[day]=Ad.phloem[day]+delta.phloem;
    }
    if  (sugar.roots[day]>0.09)
    {
      As.roots[day]=As.roots[day]+delta.roots;
    }
    else (sugar.roots[day]<0.09 && starch.roots[day]> 0)
    {
      Ad.roots[day]=Ad.roots[day]+delta.roots;
    }
    if  (sugar.xylem_sh[day]>0.049)
    {
      As.xylem_sh[day]=As.xylem_sh[day]+delta.xylem_sh;
    }
    else (sugar.xylem_sh[day]<0.049 && starch.xylem_sh[day]> 0)
    {
      Ad.xylem_sh[day]=Ad.xylem_sh[day]+delta.xylem_sh;
    }
    if  (sugar.xylem_st[day]>0.32)
    {
      As.xylem_st[day]=As.xylem_st[day]+delta.xylem_st;
    }
    else (sugar.xylem_st[day]<0.32 && starch.xylem_sh[day]> 0)
    {
      Ad.xylem_st[day]=Ad.xylem_st[day]+delta.xylem_st;
    }
  }
  } // End of the daily timeloop

  /*
   * OUTPUT!
   */

  // Output structure
  carbo_balance carbohydrate_balance;
  carbohydrate_balance.sugar = sugar;
  carbohydrate_balance.statch = starch;

  return(carbohydrate_balance) //

}



/*** R
library(Rcpp)
sugar_model()

## TODO: again make sure that this makes sense in the context and also make sure that this is if hte sperling model is true!

if (phloem.trigger == T) {
  if (count%%2 != 0) {
    if (length(which(sugar.phloem < sperling[c("SCb"),c(site)])) == 0) {
      warning(paste("Never cold enough for the model to trigger bud burst, sugar never lower than", sperling[c("SCb"),c(site)], "kg C"))
    } else {
      sB0 <- which(sugar.phloem < sperling[c("SCb"),c(site)])[1]
    }
  }
} else {
  if (count%%2 != 0) {
    if (length(which(sugar.needles+sugar.phloem+sugar.roots+sugar.xylem.sh+sugar.xylem.st < sperling[c("SCb"),c(site)])) == 0) {
      warning(paste("Never cold enough for the model to trigger bud burst, sugar never lower than", sperling[c("SCb"),c(site)], "kg C"))
    } else {
      sB0 <- which(sugar.needles+sugar.phloem+sugar.roots+sugar.xylem.sh+sugar.xylem.st < sperling[c("SCb"),c(site)])[1]
    }
  }
}

## TODO: make sure that the variable names are changed here, or at least make sense with the output from the sugar model
if (storage.reset == FALSE) {
  if (n.year != 1) {
    sperling[c("starch.needles0"), c(site)] = starch.needles[n.days]
    sperling[c("starch.phloem0"), c(site)] = starch.phloem[n.days]
    sperling[c("starch.xylem.sh0"), c(site)] = starch.xylem.sh[n.days]
    sperling[c("starch.xylem.st0"), c(site)] = starch.xylem.st[n.days]
    sperling[c("starch.roots0"), c(site)] = starch.roots[n.days]

    sperling[c("sugar.needles0"), c(site)] = sugar.needles[n.days]
    sperling[c("sugar.phloem0"), c(site)] = sugar.phloem[n.days]
    sperling[c("sugar.xylem.sh0"), c(site)] = sugar.xylem.sh[n.days]
    sperling[c("sugar.xylem.st0"), c(site)] = sugar.xylem.st[n.days]
    sperling[c("sugar.roots0"), c(site)] = sugar.roots[n.days]

    sperling[c("Ad0.needles"),c(site)] = Ad.needles[n.days]
    sperling[c("Ad0.phloem"),c(site)] = Ad.phloem[n.days]
    sperling[c("Ad0.roots"),c(site)] = Ad.roots[n.days]
    sperling[c("Ad0.xylem.sh"),c(site)] = Ad.xylem.sh[n.days]
    sperling[c("Ad0.xylem.st"),c(site)] = Ad.xylem.st[n.days]
  }
}
*/
