# Declarations
import os
import time
import cantera as ct
import numpy as np
import pandas as pd
from omnisoot import PlugFlowReactor, SootGas
# Constants
Av = 6.0221408e+23; #1/mol
MECH_PATH = {
    "Caltech":"./data/Caltech.yaml",
}

# -------------------------------------------------------------------------
# This function simulates the flow using PFR reactor of omnisoot and
# writes the simulation results on the disk
# -------------------------------------------------------------------------

def simulate_flowreactor(
    mech_name,
    PAH_growth_model_type,
    particle_dynamics_model_type,
    precursors,
    inlet_composition,
    mdot,
    P,
    temperature_profile,
    reactor_length, 
):
    # Creating the solution object and passing it to SootGas which is wrapper connecting chemistry to soot model
    gas = ct.Solution(MECH_PATH[mech_name]);
    soot_gas = SootGas(gas);
    # Configuring the reactor
    pfr = PlugFlowReactor(soot_gas);
    pfr.inlet.TPX = temperature_profile.T[0], P, inlet_composition;
    pfr.inlet.mdot = mdot;
    fixed_profile = np.vstack((temperature_profile.z, temperature_profile.T))
    pfr.set_fix_temperature_profile(fixed_profile);
    pfr.max_step = 1e-4;
    # Configuring soot model
    soot = soot = pfr.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);
    # Getting the handle to particle dynamics model
    particle_dynamics = soot.particle_dynamics_model
    # Starting the reactor
    pfr.start()
    # Creating the states object to snapshot the solution vector during the simulation
    states = ct.SolutionArray(gas, 1, extra={'restime': [0.0] , 'z':[0.0], 'velocity': [pfr.u],
                                     'N_agg': [soot.N_agg], 'N_pri': [soot.N_pri], 
                                     'C_tot': [soot.C_tot], 'H_tot': [soot.H_tot],
                                     'A_tot': [soot.A_tot],
                                     'd_p': [soot.d_p], 'd_m': [soot.d_m],
                                     'd_v': [soot.d_v], 'd_g': [soot.d_g],
                                     'n_p': [soot.n_p],
                                     'total_mass' : [soot.total_mass],
                                     'volume_fraction':[soot.volume_fraction],
                                     'SSA': [soot.SSA],
                                     'inception_mass': [soot.inception_mass],
                                     'coagulation_mass': [soot.coagulation_mass],
                                     'PAH_adsorption_mass' : [soot.PAH_adsorption_mass],
                                     'surface_growth_mass' : [soot.surface_growth_mass],
                                     'oxidation_mass': [soot.oxidation_mass],
                                     'total_carbon_flux' : [pfr.total_carbon_flux],
                                     'soot_carbon_flux' : [pfr.soot_carbon_flux]
                                        });

    # Marching forward along the reactor until the integrated length reaches the reactor length
    start_time = time.time()
    step = 0
    while pfr.z < reactor_length:
        pfr.step();
        step +=1;
        if step % 20 == 0:
            states.append(gas.state ,z = pfr.z, restime=pfr.restime, velocity = pfr.u,
                N_agg = soot.N_agg, N_pri =soot.N_pri,
                C_tot = soot.C_tot, H_tot = soot.H_tot, A_tot = soot.A_tot,
                d_p = soot.d_p, d_m = soot.d_m, d_v = soot.d_v, d_g = soot.d_g,
                n_p = soot.n_p,
                total_mass = soot.total_mass, volume_fraction = soot.volume_fraction, SSA = soot.SSA,
                inception_mass = soot.inception_mass,
                coagulation_mass = soot.coagulation_mass,
                PAH_adsorption_mass = soot.PAH_adsorption_mass,
                surface_growth_mass = soot.surface_growth_mass,
                oxidation_mass = soot.oxidation_mass,
                total_carbon_flux = pfr.total_carbon_flux,
                soot_carbon_flux = pfr.soot_carbon_flux
            );

    end_time = time.time();
    print(f"The simulation time took {end_time-start_time} s");
    
    # Creating the dataframe from solution vector
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];
    species_columns = [];
    species_data = [];

    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]", "N_agg[#/g]", "N_pri[#/g]"]
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6,
                 states.N_agg*Av/1e3, states.N_pri*Av/1e3];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "d_v[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.d_v*1e9, states.n_p];

    soot_columns += ["soot_mass[ug/g]", "volume_fraction[-]", "SSA[m2/g]",
                     "total_carbon_flux[kg/m2-s]", "carbon_yield[-]"];
    soot_data += [states.total_mass*1e6, states.volume_fraction, states.SSA,
                  states.total_carbon_flux, states.soot_carbon_flux/states.total_carbon_flux[0]];

    soot_columns += ['inception_mass[mol/kg-s]', "coagulation_mass[mol/kg-s]",
                     "PAH_adsorption_mass[mol/kg-s]", "surface_growth_mass[mol/kg-s]"];
    soot_data += [states.inception_mass, states.coagulation_mass,
                  states.PAH_adsorption_mass, states.surface_growth_mass];

    flow_columns += ["t[s]", "z[m]", "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "velocity[m/s]"];
    flow_data += [states.restime, states.z ,states.T, states.density, states.mean_molecular_weight/1000, states.velocity];

    species_columns = [f"{sp}" for sp in gas.species_names];
    species_data += [states.X[:,i] for i in range(len(gas.species_names))];
    
    columns = flow_columns + soot_columns + species_columns;
    data = (np.array(flow_data + soot_data + species_data)).transpose();
    df = pd.DataFrame(data = data, columns = columns);
    
    # Creating the results directory and storing the dataframe on disk
    output_dir = f"results/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);

    file_name = f"{output_dir}/sim_results.csv";
    df.to_csv(file_name, index=False);
    print(f"{file_name} was written!")