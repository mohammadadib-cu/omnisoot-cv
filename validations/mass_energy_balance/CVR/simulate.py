import cantera as ct
from omnisoot import ConstantVolumeReactor, SootGas, SootThermo
import numpy as np
import pandas as pd
import os
import time

Av = 6.0221408e+23; #1/mol;
mech_dict = {
    "Caltech":"data/Caltech.yaml",
}


def create_extra_dict(reactor):
    soot = reactor.soot
    extra_dict_base = dict(
        # Gas
        restime = reactor.restime,
        reactor_volume = reactor.reactor_volume,
        gas_volume = reactor.gas_volume,
        gas_mass = reactor.gas_mass,
        total_soot_mass = soot.total_mass,
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
    mech_name, 
    PAH_growth_model_type, 
    particle_dynamics_model_type,
    T5_K, 
    P5_atm, 
    initial_composition, 
    max_res_time, 
    precursors,
    solver_type,
):
    gas = ct.Solution(mech_dict[mech_name]);
    soot_gas = SootGas(gas);
    soot_gas.TPX = T5_K, P5_atm * ct.one_atm, initial_composition;
    ### Building and configruing the reactor
    reactor = ConstantVolumeReactor(soot_gas);
    reactor.max_step = 1e-3;
    reactor.atol = 1e-15;
    reactor.rtol = 1e-10;
    reactor.solver_type = solver_type;
    
    soot = reactor.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    pdynamics = soot.particle_dynamics_model;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);
    ### Starting the reactor
    reactor.start();
    # Creating the solution array for saving the solution
    if particle_dynamics_model_type == "Sectional":
        n_secs = pdynamics.number_of_sections;
        sectional_header = [f"N_agg_{i}" for i in range(n_secs)] + [f"N_pri_{i}" for i in range(n_secs)] + [f"n_p_{i}" for i in range(n_secs)] + [f"d_p_{i}" for i in range(n_secs)] + [f"d_m_{i}" for i in range(n_secs)] + [f"d_g_{i}" for i in range(n_secs)];
        sectional_data = [pdynamics.N_agg_sec_all() + pdynamics.N_pri_sec_all() + pdynamics.n_p_sec_all() + 
                    pdynamics.d_p_sec_all() + pdynamics.d_m_sec_all() + pdynamics.d_g_sec_all()];
        
    
    def append_sectional():
        params = [pdynamics.N_agg_sec_all(), pdynamics.N_pri_sec_all(), 
                      pdynamics.n_p_sec_all(), pdynamics.d_p_sec_all(),
                      pdynamics.d_m_sec_all(), pdynamics.d_g_sec_all()];
        for pi, param in enumerate(params):
            for si in range(n_secs): 
                sectional_data[pi*n_secs+si].append(param[si]);

    def append():
        extra_dict = create_extra_dict(reactor);
        states.append(gas.state, **extra_dict);
        if particle_dynamics_model_type == "Sectional": 
            append_sectional();
    
    extra_dict = create_extra_dict(reactor);
    states = ct.SolutionArray(gas, 1, extra = make_values_list(extra_dict));

    sectional_data = None;
    if particle_dynamics_model_type == "Sectional":
        n_secs = pdynamics.number_of_sections;
        sectional_header = (
            [f"N_agg_{i}" for i in range(n_secs)] + [f"N_pri_{i}" for i in range(n_secs)] 
            + [f"n_p_{i}" for i in range(n_secs)] + [f"d_p_{i}" for i in range(n_secs)] 
            + [f"d_m_{i}" for i in range(n_secs)] + [f"d_g_{i}" for i in range(n_secs)]
        );
        sectional_data = [[] for _ in sectional_header];
        append_sectional();
    

    ### Running the simulation
    start_time = time.time();
    step = 0;
    while reactor.restime < max_res_time:
        reactor.step();
        step += 1;
        if step%20==0:
            append();
        if step % 400 == 0:
            log = f"time = {reactor.restime:0.4e}";
            print(log);
    else:
        append();

    end_time = time.time();
    print(f"The simulation is finished in {(end_time-start_time):0.2f} s")
    
    ### Building the output data
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];
    species_columns = [];
    species_data = [];

    soot_columns += ["total_soot_mass[kg]"]

    soot_data += [states.total_soot_mass * states.gas_mass];

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
    

    # Adding Enthalpy to the vectors
    U_gas = states.gas_mass * states.int_energy_mass;
    U_soot = (
        states.total_soot_mass * states.gas_mass *
        np.array(list(SootThermo.u_mass_soot(T) for T in states.T))
    );

    soot_columns += ["U_soot[J]"];
    soot_data += [U_soot];

    if sectional_data:
        soot_columns += sectional_header;
        soot_data += sectional_data;


    flow_columns += ["t[s]", "P[Pa]" , "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "gas_mass[kg]", "gas_volume[m3]",
                     "reactor_volume[m3]"];
    flow_data += [states.restime, states.P, states.T, states.density, states.mean_molecular_weight/1000,
             states.gas_mass, states.gas_volume, states.reactor_volume];

    flow_columns += ["U_gas[J]"];
    flow_data += [U_gas];

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
    output_dir = f"results/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    ### Writing the file on disk
    filename = f"{output_dir}/sim_results.csv"
    df.to_csv(filename);
    print(f"{filename} was written to the disk!")