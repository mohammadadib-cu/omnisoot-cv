import cantera as ct
from omnisoot import PressureReactor, SootGas, SootThermo
import numpy as np
import pandas as pd
import os
import time

# Constants
Av = 6.0221408e+23; #1/mol;
# Axillary Functions
def make_extra_dict(reactor):
    soot = reactor.soot;
    extra_dict_base = dict(
        # Gas
        restime = reactor.restime,
        reactor_volume = reactor.reactor_volume,
        gas_volume = reactor.gas_volume,
        gas_mass = reactor.gas_mass,
        # Main Soot Variables
        N_pri = soot.N_pri,
        N_agg = soot.N_agg,
        C_tot = soot.C_tot,
        H_tot = soot.H_tot,
        # Derived Soot Variables
        A_tot = soot.A_tot,
        d_p = soot.d_p,
        d_m = soot.d_m,
        d_v = soot.d_v,
        d_g = soot.d_g,
        n_p = soot.n_p,
        total_soot_mass = soot.total_mass,
        volume_fraction = soot.volume_fraction,
        SSA = soot.SSA,
        inception_mass = soot.inception_mass,
        C_tot_inception = soot.C_tot_inception,
        C_tot_PAH = soot.C_tot_PAH,
        C_tot_growth = soot.C_tot_growth,
        total_carbon_mass = reactor.total_carbon_mass,
        total_hydrogen_mass = reactor.total_hydrogen_mass,
        soot_carbon_mass = reactor.soot_carbon_mass
    );
    return extra_dict_base;

def simulate_shocktube(
    mech_name, 
    gas,
    PAH_growth_model_type, 
    particle_dynamics_model_type,
    T5_K, 
    P5_atm, 
    initial_composition, 
    max_res_time, 
    precursors,
    inception_prefactor = 1,
    adsorption_prefactor = 1,
    results_dir = "results",
):
    ### ----------------------------------------------------------------------------------------------------------------------------
    ### This part of the code configures the reactor and soot model
    ### Setting up soot gas
    soot_gas = SootGas(gas);
    soot_gas.TPX = T5_K, P5_atm * ct.one_atm, initial_composition;
    ### Building and configruing the reactor
    reactor = PressureReactor(soot_gas);
    reactor.pressure_type = "constant";

    ### Getting the handle to soot wrapper in the reactor
    soot = reactor.soot;
    ### Setting the particle dynamics and PAH growth model
    soot.particle_dynamics_model_type = particle_dynamics_model_type
    soot.PAH_growth_model_type = PAH_growth_model_type;
    ### Setting soot precursors
    soot.set_precursor_names(precursors);
    ### Getting the handle to the "Particle Dynamics Model"
    pdynamics = soot.particle_dynamics_model;
    ### Getting the handle to "PAH Growth Model"
    pahgrowth = soot.PAH_growth_model;
    ### Applying prefactors to adjust inception and PAH adsorption rates
    pahgrowth.inception_prefactor = inception_prefactor;
    pahgrowth.adsorption_prefactor = adsorption_prefactor;
    ### Solver configurations
    reactor.max_step = 1e-3;
    reactor.solver_type = "BDF";
    reactor.atol = 1e-15;
    reactor.rtol = 1e-9;
    ### Starting the reactor
    reactor.start();
    ### ----------------------------------------------------------------------------------------------------------------------------
    ### This part of the code creates the Solution Array object of cantera to keep record of the solution during simulation
    extra_dict = make_extra_dict(reactor);
    states = ct.SolutionArray(gas, extra = list(extra_dict.keys()));
    ### Running the simulation
    step = 0;
    while reactor.restime < max_res_time:
        reactor.step();
        step += 1;
        if step % 10 == 0:
            extra_dict = make_extra_dict(reactor);
            states.append(gas.state, **extra_dict);
    else:
        extra_dict = make_extra_dict(reactor);
        states.append(gas.state, **extra_dict);
   
    ### Building the output data
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];
    species_columns = [];
    species_data = [];

    # Adding soot data
    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/kg]", "N_pri[#/kg]"]
    soot_data += [states.A_tot, states.N_agg * Av, states.N_pri * Av];

    soot_columns += ["N_agg[#/m3]", "N_pri[#/m3]"]
    soot_data += [states.N_agg * Av * states.density, states.N_pri * Av * states.density];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "d_v[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.d_v*1e9, states.n_p];

    soot_columns += ["total_soot_mass[kg]", "volume_fraction", "SSA[m2/g]"]
    soot_data += [states.total_soot_mass * states.gas_mass, states.volume_fraction, states.SSA];

    soot_columns += ["total_carbon_mass[kg]", "carbon_yield"];
    soot_data += [states.total_carbon_mass * states.gas_mass, states.soot_carbon_mass/states.total_carbon_mass[0]]

    soot_columns += ["total_hydrogen_mass[kg]"];
    soot_data += [states.total_hydrogen_mass * states.gas_mass]

    ## Total Carbon & Hydrogen Mass
    total_carbon_mass = states.total_carbon_mass * states.gas_mass;
    total_hydrogen_mass = states.total_hydrogen_mass * states.gas_mass;

    soot_columns += ["total_carbon_mass[kg]", "total_hydrogen_mass[kg]"];
    soot_data += [total_carbon_mass, total_hydrogen_mass]
    
    soot_columns += ["C_tot_inception[mol/kg-s]", "C_tot_PAH[mol/kg-s]", "C_tot_growth[mol/kg-s]"];
    soot_data += [states.C_tot_inception, states.C_tot_PAH, states.C_tot_growth];
    

    # Adding Enthalpy to the vectors
    H_gas = states.gas_mass * states.enthalpy_mass;
    H_soot = (
        states.total_soot_mass * states.gas_mass *
        np.array(list(SootThermo.h_mass_soot(T, P) for T, P in zip(states.T, states.P)))
    );

    soot_columns += ["H_soot[J]"];
    soot_data += [H_soot];

    flow_columns += ["t[s]", "P[Pa]" , "T[K]", "density[kg/m3]", "mean_molecular_weight[kg/mol]", "gas_mass[kg]", "gas_volume[m3]",
                     "reactor_volume[m3]"];
    flow_data += [states.restime, states.P, states.T, states.density, states.mean_molecular_weight/1000,
             states.gas_mass, states.gas_volume, states.reactor_volume];
    
    flow_columns += ["H_gas[J]"];
    flow_data += [H_gas];

    species_columns = [f"{sp}" for sp in gas.species_names];
    species_data += [states.X[:,i] for i in range(len(gas.species_names))];
    
    ### Combining flow, soot & species columns 
    columns = flow_columns + soot_columns + species_columns;
    data = np.vstack((
        np.array(flow_data + soot_data + species_data),
    )).transpose()

    ### creating data frame
    df = pd.DataFrame(data = data, columns = columns);
    ### making the directory
    output_dir = f"{results_dir}/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    ### Writing the file on disk
    filename = f"{output_dir}/sim_results.csv"
    df.to_csv(filename);
    print(f"{filename} was generated!")