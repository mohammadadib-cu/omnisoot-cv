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
# -------------------------------------------------------------------------
# This function simulates the flow using PFR reactor of omnisoot and
# writes the simulation results on the disk
# -------------------------------------------------------------------------

def simulate_flowreactor(
    fl_category,
    gas,
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
    inception_prefactor = 1,
    adsorption_prefactor = 1,
    max_step = 1e-3,
    solver_type = "LSODA",
):

    soot_gas = SootGas(gas);
    # Configuring the reactor
    pfr = PlugFlowReactor(soot_gas);
    pfr.inlet.TPX = temperature_profile.T[0], P, inlet_composition;
    pfr.inlet.mdot = mdot;
    pfr.area = math.pi / 4.0 * reactor_diameter * reactor_diameter;
    fixed_profile = np.vstack((temperature_profile.z, temperature_profile.T));
    pfr.temperature_solver_type = "profile_length";
    pfr.set_fix_temperature_profile(fixed_profile);
    # Configuring soot model
    soot = pfr.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);
    soot.PAH_growth_model.inception_prefactor = inception_prefactor;
    soot.PAH_growth_model.adsorption_prefactor = adsorption_prefactor;
    
    # Getting the handle to particle dynamics model
    particle_dynamics = soot.particle_dynamics_model;
    particle_dynamics.eta_coag_type = "repulsion_d_based";
    # Starting the reactor
    pfr.max_step = max_step;
    pfr.solver_type = solver_type;
    pfr.start();
    # Creating the solution array for saving the solution
    sectional_data = None;
    
    def append_sectional():
        pdynamics = soot.particle_dynamics_model;
        n_secs = pdynamics.number_of_sections;
        params = [pdynamics.N_agg_sec_all(), pdynamics.d_m_sec_all()];
        for pi, param in enumerate(params):
            for si in range(n_secs): 
                sectional_data[pi*n_secs+si].append(param[si]);
                
    if particle_dynamics_model_type == "Sectional":
        n_secs = particle_dynamics.number_of_sections;
        sectional_header = (
            [f"N_agg_{i}" for i in range(n_secs)] + [f"d_m_{i}" for i in range(n_secs)]
        );
        sectional_data = [[] for _ in sectional_header];
        append_sectional();
        
    # Creating the states object to snapshot the solution vector during the simulation
    extras = {'restime': [0.0] , 'z':[0.0], 'velocity': [pfr.u], 'area' : [pfr.area],
             'N_agg': [soot.N_agg], 'N_pri': [soot.N_pri], 
             'C_tot': [soot.C_tot], 'H_tot': [soot.H_tot],
             'A_tot': [soot.A_tot],
             'd_p': [soot.d_p], 'd_m': [soot.d_m],
             'd_v': [soot.d_v], 'd_g': [soot.d_g],
             'n_p': [soot.n_p],
             'inception_mass' : [soot.inception_mass],
             'volume_fraction':[soot.volume_fraction],
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
            inception_mass = soot.inception_mass,
            volume_fraction = soot.volume_fraction,
            total_carbon_mass_flow = pfr.total_carbon_mass_flow,
            total_hydrogen_mass_flow = pfr.total_hydrogen_mass_flow,
            soot_carbon_mass_flow = pfr.soot_carbon_mass_flow,
            C_tot_inception = soot.C_tot_inception,
            C_tot_PAH = soot.C_tot_PAH,
            C_tot_growth = soot.C_tot_growth,
        );
        if particle_dynamics_model_type == "Sectional": 
            append_sectional();

    # Marching forward along the reactor until the integrated length reaches the reactor length
    start_time = time.time()
    step = 0
    while pfr.z < reactor_length:
        pfr.step();
        step +=1;
        if step % 20 == 0:
            append();
        if step % 1000 == 0:
            print(f"Step {step}, z = {pfr.z:.3f} m, restime = {pfr.restime:.3f} s");
    else:
        append();

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

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]"];
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "d_v[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.d_v*1e9, states.n_p];

    soot_columns += ["total_carbon_mass_flow[kg/s]", "total_hydrogen_mass_flow", "soot_carbon_mass_flow[kg/s]"];
    soot_data += [states.total_carbon_mass_flow, states.total_hydrogen_mass_flow, states.soot_carbon_mass_flow];

    soot_columns += ["volume_fraction[-]", "carbon_yield[-]"];
    soot_data += [states.volume_fraction, states.soot_carbon_mass_flow/states.total_carbon_mass_flow[0]];
    
    soot_columns += ['inception_mass[mol/kg-s]'];
    soot_data += [states.inception_mass];
    
    soot_columns += ['C_tot_inception[mol/kg-s]', "C_tot_PAH[mol/kg-s]", "C_tot_growth[mol/kg-s]"];
    soot_data += [states.C_tot_inception, states.C_tot_PAH, states.C_tot_growth];

    flow_columns += ["t[s]", "z[m]", "T[K]", "density[kg/m3]", "area[m2]", "mean_MW[kg/mol]", "velocity[m/s]"];
    flow_data += [states.restime, states.z ,states.T, states.density, states.area, states.mean_molecular_weight/1000, states.velocity];

    
    species_columns = [f"{sp}" for sp in precursors];
    species_data += [states.X[:,i] for i in range(len(precursors))];

    if sectional_data:
        soot_columns += sectional_header;
        soot_data += sectional_data;
    
    columns = flow_columns + soot_columns + species_columns;
    data = (np.array(flow_data + soot_data + species_data)).transpose();
    df = pd.DataFrame(data = data, columns = columns);
    
    # Creating the results directory and storing the dataframe on disk
    output_dir = f"results/{fl_category}/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);

    file_name = f"{output_dir}/sim_results.csv";
    df.to_csv(file_name, index=False);
    print(f"{file_name} was written!")