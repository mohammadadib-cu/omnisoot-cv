import warnings
from colorama import Fore

import numpy as np
from numpy.testing import assert_allclose

from omnisoot import SootGas, SootWrapper, constants

from constants import (
    PRECURSOR_NAMES,
    ALL_PAH_GROWTH_MODEL_TYPES,
    STANDALONE_PAH_GROWTH_MODEL_TYPES,
    STANDALONE_PARTICLE_DYNAMIC_MODEL_TYPES,
    HOT_GAS_COMPOSITION,
    CONSERVED_ATOL,
    CONSERVED_RTOL,
    PREDEFINED_SOOT_VAR_LIST,
    N_SECS,
)


def test_pah_growth_loading_updating(gas):
    print(Fore.CYAN + "Checking loading and updating of PAH growth models...")
    soot_gas = SootGas(gas);
    soot = SootWrapper(soot_gas);

    particle_dynamics_model_type = "Monodisperse";
    soot.particle_dynamics_model_type = "Monodisperse";
    # Getting the handle to
    particle_dynamics_model = soot.particle_dynamics_model;
    # Initializing particle dynamics model
    particle_dynamics_model.set_vars(np.array(PREDEFINED_SOOT_VAR_LIST[particle_dynamics_model_type]));
    # Setting the PAH growth model
    PAH_growth_model = soot.PAH_growth_model
    for PAH_growth_model_type in ALL_PAH_GROWTH_MODEL_TYPES:
        soot.PAH_growth_model_type = PAH_growth_model_type;
        soot.set_precursor_names(PRECURSOR_NAMES);
        particle_dynamics_model.update();



def test_pah_growth_mass_balance(gas):
    print(Fore.CYAN + "Checking elemental mass balance in PAH growth models...")
    soot_gas = SootGas(gas);
    soot_gas.TPX = 1800, 101350, HOT_GAS_COMPOSITION;
    soot = SootWrapper(soot_gas);
    
    # Enabling PAH growth model flags and disabling the rest
    soot.inception_enabled = True;
    soot.PAH_adsorption_enabled = True;
    soot.surface_growth_enabled = False;
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
        # Setting the PAH growth model
        for PAH_growth_model_type in ALL_PAH_GROWTH_MODEL_TYPES:
            soot.PAH_growth_model_type = PAH_growth_model_type;
            PAH_growth_model = soot.PAH_growth_model
            PAH_growth_model.set_precursor_names(PRECURSOR_NAMES);
            particle_dynamics_model.update();
            # Carbon & hydrogen mass addition/removal to/from soot
            dC_tot_dt_mass_particle = particle_dynamics_model.C_tot_inception() + particle_dynamics_model.C_tot_PAH(); # [mol_c/kg-s]
            dH_tot_dt_mass_particle = particle_dynamics_model.H_tot_inception() + particle_dynamics_model.H_tot_PAH(); # [mol_h/kg-s]
            # Scrub rate
            w_scrub = particle_dynamics_model.w_scrub_debug(); # [mol_species/m3-s]
            dC_tot_dt_mass_gas = particle_dynamics_model.w_scrub_elemental("C")
            dH_tot_dt_mass_gas = particle_dynamics_model.w_scrub_elemental("H")

            # Carbon mass balance
            dC_tot_dt_residual = np.abs(dC_tot_dt_mass_particle + dC_tot_dt_mass_gas);
            if dC_tot_dt_residual > CONSERVED_ATOL:
                print(
                    Fore.YELLOW + f"Carbon mass residual of {PAH_growth_model_type} is {dC_tot_dt_residual:0.3e} > atol ({CONSERVED_ATOL:0.3})"
                )

            # Hydrogen mass balance
            dH_tot_dt_residual = np.abs(dH_tot_dt_mass_particle + dH_tot_dt_mass_gas);
            if dH_tot_dt_residual > CONSERVED_ATOL:
                print(
                    Fore.YELLOW + f"Hydrogen mass residual of {PAH_growth_model_type} is {dH_tot_dt_residual:0.3e} > atol ({CONSERVED_ATOL:0.3})"
                )


def test_pah_growth_hybrid_model(gas): 
    print(Fore.CYAN + "Checking hybrid model...")
    soot_gas = SootGas(gas);
    soot_gas.TPX = 1800, 101350, HOT_GAS_COMPOSITION;
    soot = SootWrapper(soot_gas);

    # Enabling PAH growth model flags and disabling the rest
    soot.inception_enabled = True;
    soot.PAH_adsorption_enabled = True;
    soot.surface_growth_enabled = False;
    soot.oxidation_enabled = False;
    soot.coagulation_enabled = False;
    soot.coalescence_enabled = False;
    soot.complementary_growth_enabled = False;

    print(Fore.CYAN + "Checking if each submodel is working correctly...")
    for particle_dynamics_model_type in STANDALONE_PARTICLE_DYNAMIC_MODEL_TYPES:
        print(Fore.CYAN + f"{particle_dynamics_model_type}...")
        soot.particle_dynamics_model_type = particle_dynamics_model_type;
        # Getting the handle to particle dynamics model
        particle_dynamics_model = soot.particle_dynamics_model;
        # Initializing particle dynamics model
        particle_dynamics_model.set_vars(np.array(PREDEFINED_SOOT_VAR_LIST[particle_dynamics_model_type]));
        # Setting the PAH growth model
        standalone_dC_tot_dt_mass = {}
        for PAH_growth_model_type in STANDALONE_PAH_GROWTH_MODEL_TYPES:
            soot.PAH_growth_model_type = PAH_growth_model_type;
            PAH_growth_model = soot.PAH_growth_model
            PAH_growth_model.set_precursor_names(PRECURSOR_NAMES);
            particle_dynamics_model.update();       
            standalone_dC_tot_dt_mass[PAH_growth_model_type] = {
                "inception" : particle_dynamics_model.C_tot_inception(),
                "adsorption": particle_dynamics_model.C_tot_PAH(),
            };
        for PAH_growth_model_type in STANDALONE_PAH_GROWTH_MODEL_TYPES:
            soot.PAH_growth_model_type = "HybridDimerization";
            PAH_growth_model = soot.PAH_growth_model;
            PAH_growth_model.submodel_type_names = [PAH_growth_model_type];
            PAH_growth_model.submodels[0].set_precursor_names(PRECURSOR_NAMES);
            particle_dynamics_model.update();
            if (
                (particle_dynamics_model.C_tot_inception() == standalone_dC_tot_dt_mass[PAH_growth_model_type]["inception"]) and
                (particle_dynamics_model.C_tot_PAH() == standalone_dC_tot_dt_mass[PAH_growth_model_type]["adsorption"])
            ):
                print(Fore.GREEN + f"{PAH_growth_model_type} is successfully tested with the hybrid model.")
            else:
                raise Exception(f"{PAH_growth_model_type} is failed in testing the hybrid model");
    print(Fore.CYAN + "Checking access to partial source terms...");
    soot.PAH_growth_model_type = "HybridDimerization";
    PAH_growth_model = soot.PAH_growth_model;
    # Getting the handle to particle dynamics model
    particle_dynamics_model = soot.particle_dynamics_model;
    PAH_growth_model.submodel_type_names = STANDALONE_PAH_GROWTH_MODEL_TYPES;
    for submodel_index in range(len(PAH_growth_model.submodels)):
        PAH_growth_model.submodels[submodel_index].set_precursor_names(PRECURSOR_NAMES);
    particle_dynamics_model.update();
    for submodel_index, submodel_type in enumerate(PAH_growth_model.submodel_type_names):
        print(Fore.CYAN + f"Checking {submodel_type}...");
        print(Fore.CYAN + f"Checking inception all...");
        dump_array = PAH_growth_model.N_agg_inception_all(submodel_index);
        dump_array = PAH_growth_model.C_tot_inception_all(submodel_index);
        dump_array = PAH_growth_model.H_tot_inception_all(submodel_index);
        print(Fore.CYAN + f"Checking adsorption all...");
        dump_array = PAH_growth_model.C_tot_PAH_all(submodel_index);
        dump_array = PAH_growth_model.H_tot_PAH_all(submodel_index);

        print(Fore.CYAN + f"Checking inception each for first dimer...");
        dump_array = PAH_growth_model.N_agg_inception_each(submodel_index, 0);
        dump_array = PAH_growth_model.C_tot_inception_each(submodel_index, 0);
        dump_array = PAH_growth_model.H_tot_inception_each(submodel_index, 0);
        print(Fore.CYAN + f"Checking adsorption each for first PAH...");
        dump_array = PAH_growth_model.C_tot_PAH_each(submodel_index, 0);
        dump_array = PAH_growth_model.H_tot_PAH_each(submodel_index, 0)