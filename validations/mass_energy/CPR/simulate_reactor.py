import os
import time
import cantera as ct
import numpy as np
import pandas as pd
from omnisoot import PressureReactor, SootGas, SootThermo

Av = 6.0221408e+23; #1/mol

mechanism_dict = {
    "Caltech":"../mechanisms/Caltech.yaml",
}

def simulate(
  mech_name,
  particle_dynamics_model_type,
  PAH_growth_model_type,
  precursors,
  P5,
  T5,
  initial_composition,
  max_res_time,
  solver_type,
):
    gas = ct.Solution(mechanism_dict[mech_name]);
    soot_gas = SootGas(gas);    
    soot_gas.TPX = T5, P5, initial_composition;
    reactor = PressureReactor(soot_gas);
    reactor.max_step = 1e-4;
    reactor.first_step = 1e-20;
    reactor.rtol = 1e-10;
    reactor.atol = 1e-15;
    reactor.solver_type = solver_type;
    reactor.pressure_type = "constant";
    soot = reactor.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);
    ### Starting the reactor
    reactor.start();
    # Creating the solution array for saving the solution
    states = ct.SolutionArray(gas, 1, extra={'restime': [0.0], "gas_volume" : [reactor.reactor_volume],
                                            "gas_mass" : [reactor.gas_mass], "W_enthalpy" : [reactor.W_enthalpy],
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
                                        'total_carbon_mass' : [reactor.total_carbon_mass],
                                        'soot_carbon_mass' : [reactor.soot_carbon_mass],
                                        'total_hydrogen_mass' : [reactor.total_hydrogen_mass],
                                            });
    def append():
        states.append(gas.state ,restime=reactor.restime, gas_volume = reactor.reactor_volume, gas_mass = reactor.gas_mass,
                W_enthalpy = reactor.W_enthalpy,
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
                total_carbon_mass = reactor.total_carbon_mass,
                soot_carbon_mass = reactor.soot_carbon_mass,
                total_hydrogen_mass = reactor.total_hydrogen_mass
                );

    start_time = time.time()
    step = 0;
    while reactor.restime < max_res_time:
        step += 1;
        reactor.step();
        if step % 25 == 0:
            append();
        if step % 500 == 0:
            log = f"time = {reactor.restime:0.5e} s,   ";
            print(log);

    append();

    end_time = time.time();
    print(f"The simulation time took {end_time-start_time} s");

    # Building the output
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];

    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]", "N_agg[#/g]", "N_pri[#/g]"]
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6,
                 states.N_agg*Av/1e3, states.N_pri*Av/1e3];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "d_v[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.d_v*1e9, states.n_p];

    soot_columns += ["soot_mass[ug/g]", "volume_fraction", "SSA[m2/g]"]

    soot_data += [states.total_mass*1e6, states.volume_fraction, states.SSA];

    soot_columns += ["total_carbon_mass_per_kggas[kg/kg]", "carbon_yield"];
    soot_data += [states.total_carbon_mass, states.soot_carbon_mass/states.total_carbon_mass[0]]

    soot_columns += ["total_hydrogen_mass_kggas[kg/kg]"];
    soot_data += [states.total_hydrogen_mass]

    ## Total Carbon & Hydrogen Mass
    total_carbon_mass = states.total_carbon_mass * states.gas_mass;
    total_hydrogen_mass = states.total_hydrogen_mass * states.gas_mass;

    soot_columns += ["total_carbon_mass[kg]", "total_hydrogen_mass[kg]"];
    soot_data += [total_carbon_mass, total_hydrogen_mass]

    soot_columns += ['inception_mass[mol/kg-s]', "coagulation_mass[mol/kg-s]",
                     "PAH_adsorption_mass[mol/kg-s]", "surface_growth_mass[mol/kg-s]"];
    soot_data += [states.inception_mass, states.coagulation_mass,
                  states.PAH_adsorption_mass, states.surface_growth_mass];

    flow_columns += ["t[s]", "P[kPa]" , "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "gas_mass[kg]"];
    flow_data += [states.restime, states.P, states.T, states.density, states.mean_molecular_weight/1000, states.gas_mass];


    # Adding Enthalpy to the vectors
    H_work = states.W_enthalpy
    H_gas = states.gas_mass * states.enthalpy_mass;
    H_soot = (
        states.total_mass * states.gas_mass *
        np.array(list(SootThermo.h_mass_soot(T, P) for T, P in zip(states.T, states.P)))
    );

    soot_columns += ["H_soot[J]"];
    soot_data += [H_soot];

    flow_columns += ["H_gas[J]", "W_work[J]"];
    flow_data += [H_gas, H_work];

    ### Combining flow, soot & species columns 
    columns = flow_columns + soot_columns;
    data = np.vstack((
        np.array(flow_data + soot_data),
    )).transpose();

    ### creating data frame
    df = pd.DataFrame(data = data, columns = columns);

    output_dir = f"results/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    ### Writing the file on disk
        
    filename = f"{output_dir}/results.csv"
    df.to_csv(filename);
    print(f"{filename} was stored on disk!")

