# Declarations
import os
import time
import cantera as ct
import numpy as np
import pandas as pd
import math
from omnisoot import PlugFlowReactor, SootGas
# Constants
Av = 6.0221408e+23; #1/mol
MECH_PATH = {
    "Caltech":"data/Caltech_basic.yaml",
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
    reactor_diameter,
    solver_type
):
    # Creating the solution object and passing it to SootGas which is wrapper connecting chemistry to soot model
    gas = ct.Solution(MECH_PATH[mech_name]);
    soot_gas = SootGas(gas);
    # Configuring the reactor
    pfr = PlugFlowReactor(soot_gas);
    pfr.inlet.TPX = temperature_profile[1, 0], P, inlet_composition;
    pfr.inlet.mdot = mdot;
    pfr.inlet_area = math.pi / 4.0 * reactor_diameter * reactor_diameter;
    pfr.temperature_solver_type = "profile_length"
    pfr.set_fix_temperature_profile(temperature_profile);
    pfr.max_step = 1e-3;
    pfr.solver_type = solver_type
    # Configuring soot model
    soot = pfr.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);    
    # Getting the handle to particle dynamics model
    particle_dynamics = soot.particle_dynamics_model
    particle_dynamics.eta_coag_type = "repulsion_d_based";
    # Getting the handle to surface reactions model
    sreacs = soot.surface_reactions_model;
    # Starting the reactor
    pfr.start();
    # Creating the states object to snapshot the solution vector during the simulation
    extras = {'restime': [0.0] , 'z':[0.0], 'velocity': [pfr.u], 'area' : [pfr.area],
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
             'total_carbon_mass_flow' : [pfr.total_carbon_mass_flow],
             'total_hydrogen_mass_flow' : [pfr.total_hydrogen_mass_flow],
             'soot_carbon_mass_flow' : [pfr.soot_carbon_mass_flow],
             'C_tot_inception': [soot.C_tot_inception],
             'C_tot_PAH': [soot.C_tot_PAH],
             'C_tot_growth': [soot.C_tot_growth],
        }
    
    states = ct.SolutionArray(gas, 1, extra=extras);

    def append():
        states.append(gas.state ,z = pfr.z, restime=pfr.restime, velocity = pfr.u, area = pfr.area,
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
            total_carbon_mass_flow = pfr.total_carbon_mass_flow,
            total_hydrogen_mass_flow = pfr.total_hydrogen_mass_flow,
            soot_carbon_mass_flow = pfr.soot_carbon_mass_flow,
            C_tot_inception = soot.C_tot_inception,
            C_tot_PAH = soot.C_tot_PAH,
            C_tot_growth = soot.C_tot_growth,
        );
    # Marching forward along the reactor until the integrated length reaches the reactor length
    step = 0
    while pfr.z < reactor_length:
        pfr.step();
        step +=1;
        if step % 5 == 0:
            append();
    else:
        append();    
    # Creating the dataframe from solution vector
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];
    species_columns = [];
    species_data = [];

    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]"];
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "d_v[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.d_v*1e9, states.n_p];

    soot_columns += ["total_carbon_mass_flow[kg/s]", "total_hydrogen_mass_flow", "soot_carbon_mass_flow[kg/s]"];
    soot_data += [states.total_carbon_mass_flow, states.total_hydrogen_mass_flow, states.soot_carbon_mass_flow];

    soot_columns += ["soot_mass[kg/kg]", "volume_fraction[-]", "SSA[m2/g]", 
                     "carbon_yield[-]"];
    soot_data += [states.total_mass, states.volume_fraction, states.SSA,
                  states.soot_carbon_mass_flow/states.total_carbon_mass_flow[0]];
    
    soot_columns += ['inception_mass[mol/kg-s]', "coagulation_mass[mol/kg-s]",
                     "PAH_adsorption_mass[mol/kg-s]", "surface_growth_mass[mol/kg-s]"];
    soot_data += [states.inception_mass, states.coagulation_mass,
                  states.PAH_adsorption_mass, states.surface_growth_mass];
    
    soot_columns += ['C_tot_inception[mol/kg-s]', "C_tot_PAH[mol/kg-s]", "C_tot_growth[mol/kg-s]"];
    soot_data += [states.C_tot_inception, states.C_tot_PAH, states.C_tot_growth];

    flow_columns += ["t[s]", "z[m]", "T[K]", "density[kg/m3]", "area[m2]", "mean_MW[kg/mol]", "velocity[m/s]"];
    flow_data += [states.restime, states.z ,states.T, states.density, states.area, states.mean_molecular_weight/1000, states.velocity];

    species_columns = [f"{sp}" for sp in gas.species_names];
    species_data += [states.X[:,i] for i in range(len(gas.species_names))];
    
    columns = flow_columns + soot_columns + species_columns;
    data = (np.array(flow_data + soot_data + species_data)).transpose();
    df = pd.DataFrame(data = data, columns = columns);
    
    # Creating the results directory and storing the dataframe on disk
    output_dir = f"results/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);

    file_name = f"{output_dir}/sim_results.csv";
    df.to_csv(file_name, index=False);