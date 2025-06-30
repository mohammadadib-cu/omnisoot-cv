## Declarations
import math
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import pandas as pd
from scipy import interpolate, integrate, optimize
from omnisoot import PerfectlyStirredReactor, PlugFlowReactor, SootGas, SootThermo

## Reaction Mechanisms
MECHANISMS = {
    "Caltech_basic":"Caltech_basic.yaml",
    "KAUST": "KAUST.yaml",
    "ABF1bar" : "abf_1bar.yaml",
    "CRECK" : "CRECK_2015.yaml",
};

## Constant
rho_soot = 1800; #kg/m3
Av = 6.0221408e+23;
MW_carbon = 12.011e-3 # [kg/mol]
MW_hydrogen = 1.008e-3 # [kg/mol]
PSR_SIMULATION_TIME = 0.5;
NUMBER_SECTIONS = 60;


def set_soot_params(reactor, particle_dynamics_model_type , PAH_growth_model_type, precursors, dimcoal_effs,
                    inception_prefactor, adsorption_prefactor):
    soot = reactor.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    if particle_dynamics_model_type == "Sectional":
        soot.particle_dynamics_model.number_of_sections = NUMBER_SECTIONS;
    soot.particle_dynamics_model.eta_coag_type = "repulsion_d_based";
    soot.PAH_growth_model_type = PAH_growth_model_type
    soot.set_precursor_names(precursors);
    soot.PAH_growth_model.inception_prefactor = inception_prefactor;
    soot.PAH_growth_model.adsorption_prefactor = adsorption_prefactor;
    soot.surface_reactions_model.alpha_model = "composition";
    if PAH_growth_model_type == "DimerCoalescence":
        soot.PAH_growth_model.stick_eff = dimcoal_effs;
    
    
def create_extra_dict(reactor, for_append = False):
    soot = reactor.soot
    extra_dict_base = dict(
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
        total_mass = soot.total_mass,
        volume_fraction = soot.volume_fraction,
        SSA = soot.SSA,
        inception_mass = soot.inception_mass,
        coagulation_mass = soot.coagulation_mass,
        PAH_adsorption_mass = soot.PAH_adsorption_mass,
        surface_growth_mass = soot.surface_growth_mass,
        oxidation_mass = soot.oxidation_mass,
        C_tot_inception = soot.C_tot_inception,
        C_tot_PAH = soot.C_tot_PAH,
        C_tot_growth = soot.C_tot_growth,
    );
    extra_dict_reactor = {};
    if reactor.reactor_type == "PerfectlyStirred":
        extra_dict_reactor = dict(
            restime = reactor.restime,
            mdot_out = reactor.mdot_in,
            total_carbon_mass_flow = reactor.total_carbon_mass_flow_in,
            total_hydrogen_mass_flow = reactor.total_hydrogen_mass_flow_in,
            soot_carbon_mass_flow = reactor.soot_carbon_mass_flow_in,
            soot_mass_flow = reactor.soot_mass_flow_in,
        );
        if for_append:
            extra_dict_reactor["total_carbon_mass_flow"] = reactor.total_carbon_mass_flow_out
            extra_dict_reactor["total_hydrogen_mass_flow"] = reactor.total_hydrogen_mass_flow_out
            extra_dict_reactor["soot_carbon_mass_flow"] = reactor.soot_carbon_mass_flow_out
            extra_dict_reactor["soot_mass_flow"] = reactor.soot_mass_flow_out
            extra_dict_reactor["mdot_out"] = reactor.mdot_out
    elif reactor.reactor_type == "PlugFlow":
        extra_dict_reactor = dict(
            z = reactor.z,
            restime = reactor.restime,
            velocity = reactor.u,
            total_carbon_mass_flow = reactor.total_carbon_mass_flow,
            total_hydrogen_mass_flow = reactor.total_hydrogen_mass_flow,
            soot_carbon_mass_flow = reactor.soot_carbon_mass_flow,
            soot_mass_flow = reactor.soot_mass_flow,
        );

    extra_dict = extra_dict_base.copy();
    extra_dict.update(extra_dict_reactor);
    return extra_dict;

def append_sectional(reactor, sectional_data):
    if reactor.soot.particle_dynamics_model_type == "Sectional":
        pdynamics = reactor.soot.particle_dynamics_model;
        params = [pdynamics.N_agg_sec_all(), pdynamics.N_pri_sec_all(), 
                      pdynamics.n_p_sec_all(), pdynamics.d_p_sec_all(),
                      pdynamics.d_m_sec_all(), pdynamics.d_g_sec_all()];
        for pi, param in enumerate(params):
            for si in range(NUMBER_SECTIONS): 
                sectional_data[pi*NUMBER_SECTIONS+si].append(param[si]);

def make_values_list(extra_dict):
    new_extra_dict = extra_dict.copy();
    for key in new_extra_dict.keys():
        new_extra_dict[key] = [new_extra_dict[key]];
        
    return new_extra_dict;

def create_sectional_headers_data():
    sectional_header = (
        [f"N_agg_{i}" for i in range(NUMBER_SECTIONS)] + [f"N_pri_{i}" for i in range(NUMBER_SECTIONS)] 
        + [f"n_p_{i}" for i in range(NUMBER_SECTIONS)] + [f"d_p_{i}" for i in range(NUMBER_SECTIONS)] 
        + [f"d_m_{i}" for i in range(NUMBER_SECTIONS)] + [f"d_g_{i}" for i in range(NUMBER_SECTIONS)]
    );
    sectional_data = [[] for _ in sectional_header];
    return sectional_header, sectional_data;

def append(states, gas, reactor, sectional_data):
    extra_dict = create_extra_dict(reactor, for_append = True);
    states.append(gas.state, **extra_dict);
    if reactor.soot.particle_dynamics_model_type == "Sectional":
        append_sectional(reactor, sectional_data);


### Building the data frame of the results
def create_dataframe(reactor, states, **kwargs):
    soot_columns = [];
    soot_data = [];
    flow_columns = [];
    flow_data = [];
    species_columns = [];
    species_data = [];
    error_columns = [];
    error_data = [];

    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]"]
    soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_v[nm]", "d_g[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_v*1e9, states.d_g*1e9, states.n_p];

    soot_columns += ["soot_mass[kg/kg]", "volume_fraction[-]", "SSA[m2/g]",
                     "total_carbon_mass_flow[kg/s]", "carbon_yield[-]",
                    "total_hydrogen_mass_flow[kg/s]", "soot_mass_flow[kg/s]"];
    soot_data += [states.total_mass, states.volume_fraction, states.SSA,
                  states.total_carbon_mass_flow, 
                  (states.soot_carbon_mass_flow)/(states.total_carbon_mass_flow[0]),
                 states.total_hydrogen_mass_flow, states.soot_mass_flow];

    soot_columns += ['inception_mass[mol/kg-s]', "coagulation_mass[mol/kg-s]",
                     "PAH_adsorption_mass[mol/kg-s]", "surface_growth_mass[mol/kg-s]"];
    soot_data += [states.inception_mass, states.coagulation_mass, states.PAH_adsorption_mass, states.surface_growth_mass];
    
    soot_columns += ['C_tot_inception[mol/kg-s]', "C_tot_PAH[mol/kg-s]", "C_tot_growth[mol/kg-s]"];
    soot_data += [states.C_tot_inception, states.C_tot_PAH, states.C_tot_growth];

    soot_columns += ['soot_h[J/m3]'];
    TP = list(zip(states.T, states.P))
    soot_h = np.array([SootThermo.h_mass_soot(T,P) for (T,P) in TP]);
    soot_data += [states.soot_mass_flow*soot_h];


    if reactor.soot.particle_dynamics_model_type == "Sectional":
        soot_columns += kwargs["sectional_headers"];
        soot_data += kwargs["sectional_data"];

    if reactor.reactor_type == "PlugFlow":
        flow_columns += ["z[m]", "velocity[m/s]"];
        flow_data += [states.z, states.velocity];
        
    flow_columns += ["t[s]", "T[K]", "density[kg/m3]", "mean_MW[kg/mol]"];
    flow_data += [states.restime ,states.T, states.density, states.mean_molecular_weight/1000];

    species_columns = [f"{sp}" for sp in states.species_names];
    species_data += [states.X[:,i] for i in range(len(states.species_names))];
        
    columns = flow_columns + soot_columns + species_columns;
    data = (np.array(flow_data + soot_data + species_data)).transpose();
    return pd.DataFrame(data=data, columns = columns);

def simulate_PSR_PFR(
        mech_name,
        PAH_growth_model_type,
        particle_dynamics_model_type,
        precursors,
        pfr_length,
        pfr_diameter,
        psr_volume,
        psr_residence_time,
        eq_ratio,
        fuel,
        oxidizer,
        psr_inlet_temperature,
        pressure,
        psr_inside_temperature = 1723,
        use_inside_for_init = True,
        use_constant_temperature = False,
        mdot_corrector = 1,
        dimcoal_effs = [0.002, 0.015, 0.025, 0.004],
        inception_prefactor = 1.0,
        adsorption_prefactor = 1.0,
        solver_type = "LSODA"
    ):

    pfr_area = math.pi / 4 * pfr_diameter * pfr_diameter;

    
    ## Soot Params
    soot_dict = dict(
        particle_dynamics_model_type = particle_dynamics_model_type,
        PAH_growth_model_type = PAH_growth_model_type, 
        precursors = precursors,
        dimcoal_effs = dimcoal_effs,
        inception_prefactor = inception_prefactor,
        adsorption_prefactor = adsorption_prefactor,
    )



    ## Chemistry
    gas = ct.Solution(MECHANISMS[mech_name])
#     gas.TP = psr_inlet_temperature, pressure
    gas.TP = psr_inside_temperature, pressure
    gas.set_equivalence_ratio(eq_ratio, fuel=fuel, oxidizer=oxidizer);
   
    reactor_inside_density = gas.density;
    
    if use_constant_temperature:
        pass
    else:
        gas.TP = psr_inlet_temperature, pressure
        gas.set_equivalence_ratio(eq_ratio, fuel=fuel, oxidizer=oxidizer);
        
    soot_gas = SootGas(gas);
    ## Partially Stirred Reactor
    psr = PerfectlyStirredReactor(soot_gas);
    if use_constant_temperature:
        psr.temperature_solver_type = "isothermal";
    else:
        psr.temperature_solver_type = "energy_equation";
#     psr.mdot_in = gas.density * psr_volume / psr_residence_time;
    psr.mdot_in = mdot_corrector * reactor_inside_density * psr_volume / psr_residence_time;
    psr._default_high_initial_temperature = psr_inside_temperature;
    psr.reactor_volume = psr_volume;
    psr.P_outlet = ct.one_atm
    psr.solver_type = solver_type
    psr.max_step = 1e-4;
    set_soot_params(psr, **soot_dict);
    psr.start()
    ## Extra Variable
    extra_dict = create_extra_dict(psr);
    states = ct.SolutionArray(gas, 1, extra=make_values_list(extra_dict));
    sectional_headers, sectional_data = create_sectional_headers_data();
    append_sectional(psr, sectional_data);


    ## Time integration for reaching the steady-state
    start_time = time.time()
    step = 0;
    while psr.restime < PSR_SIMULATION_TIME:
        step += 1;
        psr.step();
        if step % 50 == 0:
            append(states, gas, psr, sectional_data);
        if step % 500 == 0:
            print(f"psr step = {step}, psr time = {psr.restime:0.5e} s");
    else:
        append(states, gas, psr, sectional_data);
        
    end_time = time.time();
    print(f"The psr simulation time is {(end_time-start_time):0.1f} s");

    ## Saving the data on disk
    output_dir = f"results/{eq_ratio}/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    prop_df = create_dataframe(psr, states, sectional_headers = sectional_headers, sectional_data = sectional_data);
    filename = f"{output_dir}/sim_results_psr.csv";
    prop_df.to_csv(f"{filename}");
    print(f"{filename} was written!")

    ## Plug Flow Reactor
    pfr = PlugFlowReactor(soot_gas);
    set_soot_params(pfr, **soot_dict);
    pfr.inlet.TPX = soot_gas.T, soot_gas.P, soot_gas.X;
    pfr.inlet.mdot = psr.mdot_out;
    pfr.inlet.soot_inlet_type = "custom";
    pfr.inlet.soot_array = psr.soot_array;
    pfr.temperature_solver_type = "energy_equation";
    pfr.area = pfr_area;
    pfr.solver_type = solver_type
    pfr.max_step = 5e-4;
    pfr.start();

    ## Extra Variable
    extra_dict = create_extra_dict(pfr);
    states = ct.SolutionArray(gas, 1, extra=make_values_list(extra_dict));
    sectional_headers, sectional_data = create_sectional_headers_data();
    append_sectional(psr, sectional_data);

    ## Running the Simulation
    start_time = time.time()
    step = 0
    while pfr.z < pfr_length:
        pfr.step();
        step +=1;
        if step % 5 == 0:
            append(states, gas, pfr, sectional_data);
    else:
        append(states, gas, pfr, sectional_data);
    end_time = time.time();
    print(f"The pfr simulation time is {(end_time-start_time):0.1f} s");

    prop_df = create_dataframe(pfr, states, sectional_headers = sectional_headers, sectional_data = sectional_data)
    filename = f"{output_dir}/sim_results_pfr.csv";
    prop_df.to_csv(f"{filename}");
    print(f"{filename} was written!");
    