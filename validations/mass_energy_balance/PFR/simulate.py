# Declarations
import os
import time
import cantera as ct
import numpy as np
import pandas as pd
import math
from omnisoot import PlugFlowReactor, SootGas, SootThermo
# Constants
Av = 6.0221408e+23; #1/mol
# -------------------------------------------------------------------------
# This function simulates the flow using PFR reactor of omnisoot and
# writes the simulation results on the disk
# -------------------------------------------------------------------------

def simulat_reactor(
    gas,
    mech_name,
    PAH_growth_model_type,
    particle_dynamics_model_type,
    precursors,
    inlet_composition,
    T_inlet,
    P,
    mdot,
    reactor_length, 
    reactor_diameter,
    solver_type,
):
    # Creating the solution object and passing it to SootGas which is wrapper connecting chemistry to soot model
    #gas = ct.Solution(MECH_PATH[mech_name]);
    soot_gas = SootGas(gas);
    # Configuring the reactor
    pfr = PlugFlowReactor(soot_gas);
    pfr.inlet.TPX = T_inlet, P, inlet_composition;
    pfr.inlet.mdot = mdot;
    pfr.inlet_area = math.pi / 4.0 * reactor_diameter * reactor_diameter;
    pfr.temperature_solver_type = "energy_equation";
    pfr.atol = 1e-15;
    pfr.rtol = 1e-10;
    pfr.first_step = 1e-23;
    pfr.solver_type = solver_type;
    # Configuring soot model
    soot = pfr.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);    
    # Getting the handle to particle dynamics model
    particle_dynamics = soot.particle_dynamics_model;
    # Starting the reactor
    pfr.start();
    # Creating the solution array for saving the solution
    sectional_data = None;
    
    def append_sectional():
        pdynamics = soot.particle_dynamics_model;
        n_secs = pdynamics.number_of_sections;
        params = [pdynamics.N_agg_sec_all(), pdynamics.N_pri_sec_all(), 
                      pdynamics.n_p_sec_all(), pdynamics.d_p_sec_all(),
                      pdynamics.d_m_sec_all(), pdynamics.d_g_sec_all()];
        for pi, param in enumerate(params):
            for si in range(n_secs): 
                sectional_data[pi*n_secs+si].append(param[si]);
                
    if particle_dynamics_model_type == "Sectional":
        n_secs = particle_dynamics.number_of_sections;
        sectional_header = (
            [f"N_agg_{i}" for i in range(n_secs)] + [f"N_pri_{i}" for i in range(n_secs)] 
            + [f"n_p_{i}" for i in range(n_secs)] + [f"d_p_{i}" for i in range(n_secs)] 
            + [f"d_m_{i}" for i in range(n_secs)] + [f"d_g_{i}" for i in range(n_secs)]
        );
        sectional_data = [[] for _ in sectional_header];
        append_sectional();
        
    # Creating the states object to snapshot the solution vector during the simulation
    extras = {'restime': [0.0] , 'z':[0.0], 'velocity': [pfr.u], 'area' : [pfr.area],
              'mdot' : [pfr.mdot],
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
             'soot_mass_flow' : [pfr.soot_mass_flow],
             'C_tot_inception': [soot.C_tot_inception],
             'C_tot_PAH': [soot.C_tot_PAH],
             'C_tot_growth': [soot.C_tot_growth],
        }
    
    states = ct.SolutionArray(gas, 1, extra=extras);

    def append():
        states.append(gas.state ,z = pfr.z, restime=pfr.restime, velocity = pfr.u, area = pfr.area, mdot = pfr.mdot,
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
            soot_mass_flow = pfr.soot_mass_flow,
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
        if step % 25 == 0:
            append();
        if step % 400 == 0:
            log = f"z = {pfr.z:0.5f}";
    else:
        append();

    end_time = time.time();
    print(f"The simulation time took {end_time-start_time} s");
    
    # Creating the dataframe from solution vector
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];

    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]"];
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "d_v[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.d_v*1e9, states.n_p];

    soot_columns += ["total_carbon_mass_flow[kg/s]", "total_hydrogen_mass_flow[kg/s]", "soot_carbon_mass_flow[kg/s]"];
    soot_data += [states.total_carbon_mass_flow, states.total_hydrogen_mass_flow, states.soot_carbon_mass_flow];

    soot_columns += ['soot_h[J/s]'];
    TP = list(zip(states.T, states.P))
    soot_h = states.soot_mass_flow * np.array([SootThermo.h_mass_soot(T,P) for (T,P) in TP]);
    soot_data += [soot_h];


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

   # gas_h = states.velocity * states.density * states.area * (1.0 - states.volume_fraction) * states.enthalpy_mass;
    gas_h = states.mdot * states.enthalpy_mass;

    flow_columns += ["t[s]", "z[m]", "T[K]", "density[kg/m3]", "area[m2]", "mean_MW[kg/mol]", "gas_mdot[kg/s]" , "velocity[m/s]", "gas_h[J/s]"];
    flow_data += [states.restime, states.z ,states.T, states.density, states.area,
                   states.mean_molecular_weight/1000, states.mdot, states.velocity, gas_h];
    
    flow_columns += ["total_h"];
    flow_data += [gas_h+soot_h]

    if sectional_data:
        soot_columns += sectional_header;
        soot_data += sectional_data;
    
    columns = flow_columns + soot_columns;
    data = (np.array(flow_data + soot_data)).transpose();
    df = pd.DataFrame(data = data, columns = columns);
    
    # Creating the results directory and storing the dataframe on disk
    output_dir = f"results/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);

    file_name = f"{output_dir}/sim_results.csv";
    df.to_csv(file_name, index=False);
    print(f"{file_name} was written!")