import cantera as ct
from omnisoot import PressureReactor, SootGas, SootThermo
import numpy as np
import pandas as pd
import os
import time

Av = 6.0221408e+23; #1/mol;

def create_extra_dict(reactor):
    soot = reactor.soot
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
    extra_dict = extra_dict_base.copy();
    return extra_dict;

def make_values_list(extra_dict):
    new_extra_dict = extra_dict.copy();
    for key in new_extra_dict.keys():
        new_extra_dict[key] = [new_extra_dict[key]];
        
    return new_extra_dict;

def simulate(
    gas,
    mech_name, 
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
    solver_type = "LSODA",
):
    soot_gas = SootGas(gas);
    soot_gas.TPX = T5_K, P5_atm * ct.one_atm, initial_composition;
    ### Building and configruing the reactor
    reactor = PressureReactor(soot_gas);
    reactor.pressure_type = "constant";
    reactor.first_step = 1e-20;
    reactor.max_step = 1e-4;
    reactor.solver_type = solver_type;
    soot = reactor.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);
    pdynamics = soot.particle_dynamics_model;
    pdynamics.eta_coag_type = "repulsion_d_based";
    if particle_dynamics_model_type == "Sectional":
        pdynamics.number_of_sections = 40;
    # PAH Growth Model
    soot.PAH_growth_model.inception_prefactor = inception_prefactor;
    soot.PAH_growth_model.adsorption_prefactor = adsorption_prefactor;
    ### Starting the reactor
    reactor.start();
    def append():
        extra_dict = create_extra_dict(reactor);
        states.append(gas.state, **extra_dict);
    
    extra_dict = create_extra_dict(reactor);
    states = ct.SolutionArray(gas, 1, extra = make_values_list(extra_dict));
    

    ### Running the simulation
    step = 0;
    while reactor.restime < max_res_time:
        reactor.step();
        step += 1;
        if step%10==0:
            append();
    else:
        append();
   
    ### Building the output data
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

    soot_columns += ["total_soot_mass[kg]", "volume_fraction", "SSA[m2/g]"]

    soot_data += [states.total_soot_mass * states.gas_mass, states.volume_fraction, states.SSA];

    soot_columns += ["total_carbon_mass[kg]", "carbon_yield"];
    soot_data += [states.total_carbon_mass * states.gas_mass,
                   states.soot_carbon_mass/states.total_carbon_mass[0]]

    soot_columns += ["total_hydrogen_mass[kg]"];
    soot_data += [states.total_hydrogen_mass * states.gas_mass]

    ## Total Carbon & Hydrogen Mass
    total_carbon_mass = states.total_carbon_mass * states.gas_mass;
    total_hydrogen_mass = states.total_hydrogen_mass * states.gas_mass;

    soot_columns += ["total_carbon_mass[kg]", "total_hydrogen_mass[kg]"];
    soot_data += [total_carbon_mass, total_hydrogen_mass]
    
    soot_columns += ['C_tot_inception[mol/kg-s]', "C_tot_PAH[mol/kg-s]",
                     "C_tot_growth[mol/kg-s]"];
    soot_data += [states.C_tot_inception, states.C_tot_PAH,
                  states.C_tot_growth];
    

    flow_columns += ["t[s]", "P[Pa]" , "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "gas_mass[kg]", "gas_volume[m3]",
                     "reactor_volume[m3]"];
    flow_data += [states.restime, states.P, states.T, states.density, states.mean_molecular_weight/1000,
             states.gas_mass, states.gas_volume, states.reactor_volume];    
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
    filename = f"{output_dir}/T{T5_K}K.csv"
    df.to_csv(filename);
    print(f"{filename} was written to the disk!")