import os
import time

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
from omnisoot import ConstantVolumeReactor, SootGas
import pandas as pd
### Constants
## Constant
rho_soot = 1800; #kg/m3
Av = 6.0221408e+23;

## Reaction Mechanisms
MECHANISMS = {
    "KAUST": "data/KAUST.yaml",
};

def create_extra_dict(reactor):
    soot = reactor.soot
    extra_dict_base = dict(
        # Flow
        restime = reactor.restime,
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

### Sectional Model Functions
def create_sectional_headers_data(reactor, number_of_sections):
    if reactor.soot.particle_dynamics_model_type == "Sectional":
        sectional_header = (
            [f"N_agg_{i}[mol/kg]" for i in range(number_of_sections)] + [f"N_pri_{i}[mol/kg]" for i in range(number_of_sections)] 
            + [f"n_p_{i}" for i in range(number_of_sections)] + [f"d_p_{i}[m]" for i in range(number_of_sections)] 
            + [f"d_m_{i}[m]" for i in range(number_of_sections)] + [f"d_g_{i}[m]" for i in range(number_of_sections)]
        );
        sectional_data = [[] for _ in sectional_header];
        return sectional_header, sectional_data;
    else:
        return [], [];

def append_sectional(reactor, **sectional_props):
    sectional_data = sectional_props.get("sectional_data");
    number_of_sections = sectional_props.get("number_of_sections");
    if reactor.soot.particle_dynamics_model_type == "Sectional":
        pdynamics = reactor.soot.particle_dynamics_model;
        params = [pdynamics.N_agg_sec_all(), pdynamics.N_pri_sec_all(), 
                      pdynamics.n_p_sec_all(), pdynamics.d_p_sec_all(),
                      pdynamics.d_m_sec_all(), pdynamics.d_g_sec_all()];
        for pi, param in enumerate(params):
            for si in range(number_of_sections): 
                sectional_data[pi*number_of_sections+si].append(param[si]);
    else:
        pass;

### Function to append data to SolutionArray
def append(states, gas, reactor, **sectional_props):
    extra_dict = create_extra_dict(reactor);
    states.append(gas.state, **extra_dict);
    append_sectional(reactor, **sectional_props);

### Building the data frame of the results
def create_dataframe(states, **sectional_args):
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];
    species_columns = [];
    species_data = [];

    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]"]
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_v[nm]", "d_g[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_v*1e9, states.d_g*1e9, states.n_p];

    soot_columns += ["soot_mass[kg/kg]",
                     "volume_fraction[-]",
                     "SSA[m2/g]",
                     "total_carbon_mass[kg/kg]",
                     "carbon_yield",
                     "total_carbon_mass[kg/m3]",
                     "total_hydrogen_mass[kg/kg]",
                     "total_hydrogen_mass[kg/m3]"];
    soot_data += [states.total_soot_mass,
                  states.volume_fraction, 
                  states.SSA,
                  states.total_carbon_mass, 
                  (states.soot_carbon_mass*states.density)/(states.total_carbon_mass[0]*states.density[0]),
                  states.total_carbon_mass*states.density,
                  states.total_hydrogen_mass, 
                  states.total_hydrogen_mass*states.density];

    soot_columns += ['inception_mass[mol/kg-s]'];
    soot_data += [states.inception_mass];
    
    soot_columns += ['C_tot_inception[mol/kg-s]', "C_tot_PAH[mol/kg-s]", "C_tot_growth[mol/kg-s]"];
    soot_data += [states.C_tot_inception, states.C_tot_PAH, states.C_tot_growth];
        
    soot_columns += sectional_args["sectional_headers"];
    soot_data += sectional_args["sectional_data"];
    
    flow_columns += ["t[s]", "T[K]", "density[kg/m3]", "mean_MW[kg/mol]"];
    flow_data += [states.restime, states.T, states.density, states.mean_molecular_weight/1000];

    species_columns = [f"{sp}" for sp in states.species_names];
    species_data += [states.X[:,i] for i in range(len(states.species_names))];

        
    columns = flow_columns + soot_columns + species_columns;
    data = (np.array(flow_data + soot_data + species_data)).transpose();
    return pd.DataFrame(data=data, columns = columns);


def simulate(
        mech_name,
        PAH_growth_model_type,
        particle_dynamics_model_type,
        precursors,
        max_res_time,
        initial_composition,
        initial_temperature,
        initial_pressure,
        **extra_props
    ):
    # Extracting extra props
    number_of_sections = extra_props.get("number_of_sections", 60);
    surface_reactivity_model_type = extra_props.get("surface_reactivity_model_type", "empirical");
    surface_reactivity_constant = extra_props.get("surface_reactivity_constant", 1.0);
    carbonization_enabled = extra_props.get("carbonization_enabled", False);
    eta_coag_type = extra_props.get("eta_coag_type", "unity");
    ## Chemistry
    gas = ct.Solution(MECHANISMS[mech_name])
    soot_gas = SootGas(gas);
    ### Reactor
    constUV = ConstantVolumeReactor(soot_gas);
    ### Initial conditions
    soot_gas.TPX = initial_temperature, initial_pressure, initial_composition;
    ### PAH growth model
    soot = constUV.soot;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    ### Particle dynamics model
    soot.particle_dynamics_model_type = particle_dynamics_model_type
    if particle_dynamics_model_type == "Sectional":
        soot.particle_dynamics_model.number_of_sections = number_of_sections;
    soot.particle_dynamics_model.eta_coag_type = eta_coag_type;
    ### Setting precursors
    soot.set_precursor_names(precursors);
    ### Surface reactiviy model
    soot.surface_reactions_model.surface_reactivity_model_type = surface_reactivity_model_type;
    if surface_reactivity_model_type == "constant":
        soot.surface_reactions_model.surface_reactivity_constant = surface_reactivity_constant;
    ### Carbonization
    soot.carbonization_enabled = carbonization_enabled;
    if carbonization_enabled:
        # Carbonization model's parameters were taken from Table 2 in Kholghy et al. (https://doi.org/10.1016/j.carbon.2016.01.022)
        soot.carbonization_model.A_carb = 1.85e9; #1/s
        soot.carbonization_model.Ea_carb = 242672; #J/mol
        soot.carbonization_model.HtoC_min = 0.005;
    ### Reactor integrator settings
    constUV.max_step = 5e-4;
    constUV.max_time = max_res_time;
    constUV.temperature_solver_type = "energy_equation";
    constUV.start();
    ### Soltuion array
    extra_dict = create_extra_dict(constUV);
    states = ct.SolutionArray(gas, 1, extra = make_values_list(extra_dict));
    ### Sectional data
    sectional_headers, sectional_data = create_sectional_headers_data(constUV, number_of_sections);
    sectional_props = {
        "number_of_sections" : number_of_sections,
        "sectional_headers" : sectional_headers, 
        "sectional_data" : sectional_data
    };
    append_sectional(constUV, **sectional_props);

    ### Running the reactor
    start_time = time.time()
    step = 0;
    while constUV.restime < max_res_time:
        step += 1;
        constUV.step();
        if step % 20 == 0:
            append(states, gas, constUV, **sectional_props);
    append(states, gas, constUV, **sectional_props);
    end_time = time.time();
    print(f"The simulation time is {end_time-start_time:0.2f} s");

    ## Saving the data on disk
    output_dir = f"results/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    
    prop_df = create_dataframe(states, **sectional_props);
    filename = f"{output_dir}/sim_results_surface_reactivity_{surface_reactivity_model_type}.csv";
    prop_df.to_csv(f"{filename}");
    print(f"{filename} was written!");