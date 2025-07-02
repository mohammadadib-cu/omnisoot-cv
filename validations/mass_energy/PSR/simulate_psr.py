import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import pandas as pd
import os
import time
from omnisoot import PerfectlyStirredReactor, SootGas, SootThermo
# Constants
rho_soot = 1800; #kg/m3
Av = 6.0221408e+23;
MW_carbon = 12.011e-3 # [kg/mol]
MW_hydrogen = 1.008e-3 # [kg/mol]

simulation_time = 1;

def simulate_psr(gas,
                 reactor_volume, 
                 pressure, 
                 residence_time, 
                 reactor_temperature,
                 inlet_temperature,
                 composition,
                 precursors,
                 particle_dynamics_model_type,
                 PAH_growth_model_type,
                 solver_type = "BDF",
                ):
    
    print(f"PAH_growth_model_type = {PAH_growth_model_type}")
    print(f"particle_dynamics_model_type = {particle_dynamics_model_type}")
    
    gas.TPX = reactor_temperature, pressure, composition
    reactor_density = gas.density;
    
    gas.TPX = inlet_temperature, pressure, composition
    
    soot_gas = SootGas(gas)
    psr = PerfectlyStirredReactor(soot_gas);
    psr.mdot_in = reactor_density * reactor_volume / residence_time;
    psr.reactor_volume = reactor_volume;
    psr.P_outlet = ct.one_atm;
    psr._default_high_initial_temperature = 1800;
    psr.max_step = 1e-3;
    psr.first_step = None;
    psr.rtol = 1e-10;
    psr.atol = 1e-15;
    psr.solver_type = solver_type;
    
    soot = psr.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    pdynamics = soot.particle_dynamics_model
    if particle_dynamics_model_type == "Sectional":
        pdynamics.spacing_factor = 1.5;
        pdynamics.number_of_sections = 40;
    
    soot.PAH_growth_model_type = PAH_growth_model_type
    soot.set_precursor_names(precursors);
    
    psr.temperature_solver_type = "energy_equation";
    psr.start()
    
    states = ct.SolutionArray(gas, 1, extra={'restime': [0.0] ,
                                         'N_agg': [soot.N_agg], 'N_pri': [soot.N_pri], 
                                         'C_tot': [soot.C_tot], 'H_tot': [soot.H_tot],
                                         'A_tot': [soot.A_tot],
                                         'd_p': [soot.d_p], 'd_m': [soot.d_m],
                                         'd_v': [soot.d_v], 'd_g': [soot.d_g],
                                         'n_p': [soot.n_p],
                                         'total_soot_mass' : [soot.total_mass],
                                         'volume_fraction':[soot.volume_fraction],
                                         'SSA': [soot.SSA],
                                         'inception_mass': [soot.inception_mass],
                                         'coagulation_mass': [soot.coagulation_mass],
                                         'PAH_adsorption_mass' : [soot.PAH_adsorption_mass],
                                         'surface_growth_mass' : [soot.surface_growth_mass],
                                         'oxidation_mass': [soot.oxidation_mass],
                                         'total_carbon_mass_flow' : [psr.total_carbon_mass_flow_in],
                                         'total_hydrogen_mass_flow' : [psr.total_hydrogen_mass_flow_in],
                                         'soot_carbon_mass_flow' : [psr.soot_carbon_mass_flow_in],
                                         'soot_mass_flow' : [psr.soot_mass_flow_in],
                                         'mdot' : [psr.mdot_in],
                                         'removed_gas_mass' : [psr.removed_gas_mass],
                                            });
    
    def append():
        states.append(gas.state ,restime=psr.restime, 
                  N_agg = soot.N_agg, N_pri =soot.N_pri,
                  C_tot = soot.C_tot, H_tot = soot.H_tot, A_tot = soot.A_tot,
                  d_p = soot.d_p, d_m = soot.d_m, d_v = soot.d_v, d_g = soot.d_g,
                  n_p = soot.n_p,
                  total_soot_mass = soot.total_mass, volume_fraction = soot.volume_fraction, SSA = soot.SSA,
                  inception_mass = soot.inception_mass,
                  coagulation_mass = soot.coagulation_mass,
                  PAH_adsorption_mass = soot.PAH_adsorption_mass,
                  surface_growth_mass = soot.surface_growth_mass,
                  oxidation_mass = soot.oxidation_mass,
                  total_carbon_mass_flow = psr.total_carbon_mass_flow_out,
                  total_hydrogen_mass_flow = psr.total_hydrogen_mass_flow_out,
                  soot_carbon_mass_flow = psr.soot_carbon_mass_flow_out,
                  soot_mass_flow = psr.soot_mass_flow_out,
                  mdot = psr.mdot_out,
                  removed_gas_mass = psr.removed_gas_mass
                 );
        
    start_time = time.time()
    step = 0;
    while psr.restime < simulation_time:
        step += 1;
        psr.step();
        if step%20==0:
            append();
        if step % 1000 == 0:
            log = f"t={psr.restime:0.4e} s";
            print(log);
    append();
    end_time = time.time();
    print(f"The simulation took {(end_time-start_time):0.4f} seconds");


    ### Building the data frame of the results
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];

    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]", "N_agg[#/g]", "N_pri[#/g]"]
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6,
                 states.N_agg*Av/1e3, states.N_pri*Av/1e3];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_v[nm]", "d_g[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_v*1e9, states.d_g*1e9, states.n_p];

    soot_columns += ["soot_mass[kg/kg]", "volume_fraction", "SSA",
                     "total_carbon_mass_flow[kg/s]", "carbon_yield",
                    "total_hydrogen_mass_flow[kg/s]", "soot_mass_flow[kg/s]"];
    soot_data += [states.total_soot_mass, states.volume_fraction, states.SSA,
                  states.total_carbon_mass_flow, 
                  (states.soot_carbon_mass_flow)/(states.total_carbon_mass_flow[0]),
                 states.total_hydrogen_mass_flow, states.soot_mass_flow];

    soot_columns += ['inception_mass[mol/kg-s]', "coagulation_mass[mol/kg-s]",
                     "PAH_adsorption_mass[mol/kg-s]", "surface_growth_mass[mol/kg-s]"];
    soot_data += [states.inception_mass, states.coagulation_mass, states.PAH_adsorption_mass, states.surface_growth_mass];

    soot_columns += ['soot_h[J/m3]'];
    TP = list(zip(states.T, states.P))
    soot_h = np.array([SootThermo.h_mass_soot(T,P) for (T,P) in TP]);
    soot_data += [states.soot_mass_flow*soot_h];

    flow_columns += ["t[s]", "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "mdot[kg/s"];
    flow_data += [states.restime ,states.T, states.density, states.mean_molecular_weight/1000, states.mdot];
    
    
    C_flux_in = psr.mdot_in * (states.elemental_mass_fraction('C')[0]+states.C_tot[0]*MW_carbon);
    C_flux_out = (states.elemental_mass_fraction('C')+states.C_tot*MW_carbon)*states.mdot;
    C_flux_error = (C_flux_out-C_flux_in)/C_flux_in;

    
    H_flux_in = psr.mdot_in * states.elemental_mass_fraction('H')[0];
    H_flux_out = (states.elemental_mass_fraction('H')+states.H_tot*MW_hydrogen)*states.mdot;
    H_flux_error = abs(H_flux_out-H_flux_in)/H_flux_in;
    
    H1 = psr.mdot_in * psr.h_in ;
    H2 = states.mdot * states.enthalpy_mass;
    Hp = (states.total_soot_mass*states.mdot) * np.array(list(SootThermo.h_mass_soot(T, P) for T, P in zip(states.T, states.P)));
    H_error = (H2-H1+Hp)/H1;
    
    error_columns = ["C_mass_error", "H_mass_error", "enthalpy_error"];
    error_data = [C_flux_error, H_flux_error, H_error]

    columns = flow_columns + soot_columns + error_columns;
    data = (np.array(flow_data + soot_data + error_data)).transpose();
    
    output_dir = f"results/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    prop_df = pd.DataFrame(data=data, columns = columns)
    prop_df.to_csv(f"{output_dir}/results.csv");
