import warnings
from colorama import Fore

import numpy as np
from numpy.testing import assert_allclose

from omnisoot import SootGas, SootWrapper, constants

from constants import (
    PRECURSOR_NAMES,
    ALL_SURFACE_REACTIONS_MODEL_TYPES,
    HOT_GAS_COMPOSITION,
    HOT_GAS_COMPOSITION_WITH_OXIDIZERS,
    CONSERVED_ATOL,
    CONSERVED_RTOL,
    PREDEFINED_SOOT_VAR_LIST,
    N_SECS,
    STANDALONE_PARTICLE_DYNAMIC_MODEL_TYPES,
)


def test_surface_reactions_loading_updating(gas):
    print(Fore.CYAN + "Checking loading and updating of surface reactions models...")
    soot_gas = SootGas(gas);
    soot = SootWrapper(soot_gas);
    
    soot.particle_dynamics_model_type = "Monodisperse";
    # Getting the handle to
    particle_dynamics_model = soot.particle_dynamics_model;
    # Initializing particle dynamics model
    particle_dynamics_model.set_vars(np.array(PREDEFINED_SOOT_VAR_LIST["Monodisperse"]));
    # Setting the PAH growth model
    surface_reactions_model = soot.surface_reactions_model
    for surface_reactions_model_type in ALL_SURFACE_REACTIONS_MODEL_TYPES:
        soot.surface_reactions_model_type = surface_reactions_model_type;
        particle_dynamics_model.update();


def test_surface_reactions_mass_balance(gas):
    print(Fore.CYAN + "Checking elemental mass balance in surface reactions models...");
    print(Fore.CYAN + "Surface growth...")
    soot_gas = SootGas(gas);
    soot_gas.TPX = 1800, 101350, HOT_GAS_COMPOSITION;
    soot = SootWrapper(soot_gas);
    
    # Enabling PAH growth model flags and disabling the rest
    soot.inception_enabled = False;
    soot.PAH_adsorption_enabled = False;
    soot.surface_growth_enabled = True;
    soot.oxidation_enabled = False;
    soot.coagulation_enabled = False;
    soot.coalescence_enabled = False;
    soot.complementary_growth_enabled = False;

    for particle_dynamics_model_type in STANDALONE_PARTICLE_DYNAMIC_MODEL_TYPES:
        print(Fore.CYAN + f"{particle_dynamics_model_type}...")
        soot.particle_dynamics_model_type = particle_dynamics_model_type;
        # Getting the handle to particle dynamics model
        particle_dynamics_model = soot.particle_dynamics_model;
        if "Sectional" in particle_dynamics_model_type:
            particle_dynamics_model.number_of_sections = N_SECS;
            particle_dynamics_model.spacing_factor = 1.5
        # Initializing particle dynamics model
        particle_dynamics_model.set_vars(np.array(PREDEFINED_SOOT_VAR_LIST[particle_dynamics_model_type]));
        surface_reactions_model = soot.surface_reactions_model
        for surface_reactions_model_type in ALL_SURFACE_REACTIONS_MODEL_TYPES:
            soot.surface_reactions_model_type = surface_reactions_model_type;
            particle_dynamics_model.update();

            # Carbon & hydrogen mass addition/removal to/from soot
            dC_tot_dt_mass_particle = particle_dynamics_model.C_tot_growth() # [mol_c/kg-s]
            dH_tot_dt_mass_particle = particle_dynamics_model.H_tot_growth() # [mol_h/kg-s]
            # Scrub rate
            dC_tot_dt_mass_gas = particle_dynamics_model.w_scrub_elemental("C")
            dH_tot_dt_mass_gas = particle_dynamics_model.w_scrub_elemental("H")
            # Carbon mass balance
            dC_tot_dt_residual = np.abs(dC_tot_dt_mass_particle + dC_tot_dt_mass_gas);
            if dC_tot_dt_residual > CONSERVED_ATOL:
                print(
                    Fore.YELLOW + f"Carbon mass residual of {surface_reactions_model_type} is {dC_tot_dt_residual:0.3e} > atol ({CONSERVED_ATOL:0.3})"
                )
            # Hydrogen mass balance
            dH_tot_dt_residual = np.abs(dH_tot_dt_mass_particle + dH_tot_dt_mass_gas);
            if dH_tot_dt_residual > CONSERVED_ATOL:
                print(
                    Fore.YELLOW + f"Hydrogen mass residual of {surface_reactions_model_type} is {dH_tot_dt_residual:0.3e} > atol ({CONSERVED_ATOL:0.3})"
                );
    #----------------------------------------------------------------------------------
    print(Fore.CYAN + "Oxidation...")
    soot_gas.TPX = 1800, 101350, HOT_GAS_COMPOSITION_WITH_OXIDIZERS;
    soot.surface_growth_enabled = False;
    soot.oxidation_enabled = True;    

    for particle_dynamics_model_type in STANDALONE_PARTICLE_DYNAMIC_MODEL_TYPES:
        print(Fore.CYAN + f"{particle_dynamics_model_type}...")
        soot.particle_dynamics_model_type = particle_dynamics_model_type;
        # Getting the handle to particle dynamics model
        particle_dynamics_model = soot.particle_dynamics_model;
        if "Sectional" in particle_dynamics_model_type:
            particle_dynamics_model.number_of_sections = N_SECS;
            particle_dynamics_model.spacing_factor = 1.5
        particle_dynamics_model.set_vars(np.array(PREDEFINED_SOOT_VAR_LIST[particle_dynamics_model_type]));
        for surface_reactions_model_type in ALL_SURFACE_REACTIONS_MODEL_TYPES:
            soot.surface_reactions_model_type = surface_reactions_model_type;
            particle_dynamics_model.update();
            # Carbon & hydrogen mass addition/removal to/from soot
            dC_tot_dt_mass_particle = particle_dynamics_model.C_tot_oxidation() # [mol_c/kg-s]
            dH_tot_dt_mass_particle = particle_dynamics_model.H_tot_oxidation() # [mol_h/kg-s]
            # Scrub rate
            dC_tot_dt_mass_gas = particle_dynamics_model.w_scrub_elemental("C")
            dH_tot_dt_mass_gas = particle_dynamics_model.w_scrub_elemental("H")
            dO_tot_dt_mass_gas = particle_dynamics_model.w_scrub_elemental("O")


            # Carbon mass balance
            dC_tot_dt_residual = np.abs(-dC_tot_dt_mass_particle + dC_tot_dt_mass_gas);
            if dC_tot_dt_residual > CONSERVED_ATOL:
                print(
                    Fore.YELLOW + f"Carbon mass residual of {surface_reactions_model_type} is {dC_tot_dt_residual:0.3e} > atol ({CONSERVED_ATOL:0.3})"
                )
            # Hydrogen mass balance
            dH_tot_dt_residual = np.abs(-dH_tot_dt_mass_particle + dH_tot_dt_mass_gas);
            if dH_tot_dt_residual > CONSERVED_ATOL:
                print(
                    Fore.YELLOW + f"Hydrogen mass residual of {surface_reactions_model_type} is {dH_tot_dt_residual:0.3e} > atol ({CONSERVED_ATOL:0.3})"
                );

            # Oxygen mass balance
            dO_tot_dt_residual = np.abs(dO_tot_dt_mass_gas);
            if dO_tot_dt_residual > CONSERVED_ATOL:
                print(
                    Fore.YELLOW + f"Oxygen mass residual of {surface_reactions_model_type} is {dO_tot_dt_residual:0.3e} > atol ({CONSERVED_ATOL:0.3})"
                );