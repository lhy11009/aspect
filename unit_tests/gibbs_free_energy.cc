/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include "common.h"

#include <aspect/material_model/reaction_model/gibbs_free_energy.h>
#include <aspect/material_model/reaction_model/gibbs_free_energy/constant.h>
#include <aspect/material_model/reaction_model/gibbs_free_energy/linearized.h>

#include <deal.II/base/parameter_handler.h>

TEST_CASE("Constant Gibbs Free Energy")
{
  using namespace aspect::MaterialModel::ReactionModel;

  dealii::ParameterHandler prm;
  GibbsFreeEnergy<2> gibbs_free_energy;
  GibbsFreeEnergy<2>::declare_parameters(prm);

  prm.enter_subsection("Gibbs free energy");
  {
    prm.set("Models", "Constant Gibbs free energy|Constant Gibbs free energy");
    prm.enter_subsection("Constant Gibbs free energy");
    prm.set("Gibbs free energy differences", "-1000|2000");
    prm.leave_subsection();
  }
  prm.leave_subsection();

  gibbs_free_energy.parse_parameters(prm, 2);
  const GibbsFreeEnergyInputs reference_state {};

  SECTION("A to B is favored")
  {
    CHECK(gibbs_free_energy.delta_gibbs_free_energy(1500, 1e9, reference_state, 0) == Approx(-1000));
  }

  SECTION("B to A is favored")
  {
    CHECK(gibbs_free_energy.delta_gibbs_free_energy(1500, 1e9, reference_state, 1) == Approx(2000));
  }
}

TEST_CASE("Linearized Gibbs Free Energy")
{
  using namespace aspect::MaterialModel::ReactionModel;

  dealii::ParameterHandler prm;
  GibbsFreeEnergy<2> gibbs_free_energy;
  GibbsFreeEnergy<2>::declare_parameters(prm);

  prm.enter_subsection("Gibbs free energy");
  prm.set("Models", "Linearized Gibbs free energy");
  prm.leave_subsection();

  gibbs_free_energy.parse_parameters(prm, 1);

  GibbsFreeEnergyInputs reference_state;
  reference_state.reference_temperature = 1500;
  reference_state.reference_pressure = 1e9;
  reference_state.reference_delta_gibbs_free_energy = -500;
  reference_state.delta_entropy = 10;
  reference_state.delta_volume = 2e-6;

  CHECK(gibbs_free_energy.delta_gibbs_free_energy(1500, 1e9, reference_state, 0) == Approx(-500));
  CHECK(gibbs_free_energy.delta_gibbs_free_energy(1600, 1.2e9, reference_state, 0) == Approx(-1100));
}

TEST_CASE("Gibbs Free Energy Model Count")
{
  using namespace aspect::MaterialModel::ReactionModel;

  dealii::ParameterHandler prm;
  GibbsFreeEnergy<2> gibbs_free_energy;
  GibbsFreeEnergy<2>::declare_parameters(prm);

  prm.enter_subsection("Gibbs free energy");
  prm.set("Models", "Constant Gibbs free energy");
  prm.leave_subsection();

  CHECK_THROWS_AS(gibbs_free_energy.parse_parameters(prm, 2), dealii::ExceptionBase);
}

TEST_CASE("Clapeyron Gibbs Free Energy")
{
  using namespace aspect::MaterialModel::ReactionModel;

  dealii::ParameterHandler prm;
  GibbsFreeEnergy<2> gibbs_free_energy;
  GibbsFreeEnergy<2>::declare_parameters(prm);

  prm.enter_subsection("Gibbs free energy");
  prm.set("Models", "Clapeyron Gibbs free energy|Constant Gibbs free energy|Clapeyron Gibbs free energy|Clapeyron Gibbs free energy");
  prm.enter_subsection("Clapeyron Gibbs free energy");
  prm.set("Clapeyron slopes", "4e6|-2e6|0");
  prm.set("Reference pressures", "1e9|1e9|1e9");
  prm.set("Reference temperatures", "1500|1500|1500");
  prm.set("Molar volume changes", "-2e-6|-2e-6|-2e-6");
  prm.leave_subsection();
  prm.leave_subsection();

  SECTION("Equilibrium and Gibbs energy on either side")
  {
    gibbs_free_energy.parse_parameters(prm, 4);

    CHECK(gibbs_free_energy.delta_gibbs_free_energy(1600, 1.4e9, {}, 0) == Approx(0));
    CHECK(gibbs_free_energy.delta_gibbs_free_energy(1600, 1.5e9, {}, 0) == Approx(-200));
    CHECK(gibbs_free_energy.delta_gibbs_free_energy(1600, 1.3e9, {}, 0) == Approx(200));
    CHECK(gibbs_free_energy.delta_gibbs_free_energy(1600, 0.8e9, {}, 2) == Approx(0));
    CHECK(gibbs_free_energy.delta_gibbs_free_energy(1600, 1e9, {}, 3) == Approx(0));
  }

  SECTION("Incorrect number of slopes")
  {
    prm.enter_subsection("Gibbs free energy");
    prm.enter_subsection("Clapeyron Gibbs free energy");
    prm.set("Clapeyron slopes", "4e6");
    prm.leave_subsection();
    prm.leave_subsection();
    CHECK_THROWS_AS(gibbs_free_energy.parse_parameters(prm, 4), dealii::ExceptionBase);
  }

  SECTION("Zero molar volume change")
  {
    prm.enter_subsection("Gibbs free energy");
    prm.enter_subsection("Clapeyron Gibbs free energy");
    prm.set("Molar volume changes", "-2e-6|0|-2e-6");
    prm.leave_subsection();
    prm.leave_subsection();
    CHECK_THROWS_AS(gibbs_free_energy.parse_parameters(prm, 4), dealii::ExceptionBase);
  }
}
