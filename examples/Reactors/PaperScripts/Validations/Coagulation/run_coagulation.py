## Declarations
from collections import defaultdict
from functools import reduce
import operator
import cantera as ct
from omnisoot import ConstantVolumeReactor, SootGas, SootThermo, constants
import numpy as np
import pandas as pd
import os
import json
import time


def run_coagulation(particle_dynamics_model_type):
    ## Cantera Solution Object
    gas = ct.Solution("data/Caltech.yaml");
    soot_gas = SootGas(gas);
    ## ---------------------------------------------------------------------------------
    ## Intial conditions of the reactor
    T = 1830; # K
    P = 101325; # Pa
    X = {'CH4':0.425, 'O2':0.435 , 'N2':0.14};
    ## ---------------------------------------------------------------------------------
    gas_density = 0.123;
    ## Reactor
    reactor = ConstantVolumeReactor(soot_gas);
    ### Initial conditions
    soot_gas.TPX = T, P, X;
    ## Simulation time
    max_res_time = 0.05;
    ## Getting the handle to soot
    soot = reactor.soot;
    ## Setting particle dynamics model
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    particle_dynamics_model = soot.particle_dynamics_model
    ## Modifying the soot options to deactivate inception, PAH adsorption, surface growth and oxidation
    soot.inception_enabled = False;
    soot.PAH_adsorption_enabled = False;
    soot.surface_growth_enabled = False;
    soot.oxidation_enabled = False;
    soot.coagulation_enabled = True;
    ## Modifying the initial soot of the reactor
    reactor.initial_soot_type = "custom";
    ## Calulating the soot variables based on the initial number of concentration of spherical particles 2 nm in diameters;
    N_0_number = 2626056626709460000; # 1/m3 from t=0 of DEM simulations
    N_0 = N_0_number / constants.Av / gas_density; # mol/kg
    # Calculating and constructing the soot array
    if particle_dynamics_model_type == "Monodisperse":
        N_agg = N_0;
        N_pri = N_0;
        C_tot = N_0 * constants.incepient_soot_C;
        H_tot = N_0 * constants.incepient_soot_H;
        soot_array = np.hstack((N_agg, N_pri, C_tot, H_tot));
    elif particle_dynamics_model_type == "Sectional":
        n_secs = particle_dynamics_model.n_secs;
        N_agg_sec = np.zeros((n_secs,));
        N_pri_sec = np.zeros((n_secs,));
        H_tot_sec = np.zeros((n_secs,));
        ## First section
        N_agg_sec[0] = N_0;
        N_pri_sec[0] = N_0;
        H_tot_sec[0] = N_0 * constants.incepient_soot_H / constants.Av;
        ## Other sections
        for sec in range(1, n_secs):
            N_agg_sec[sec] = 1.0 / constants.Av;    
            N_pri_sec[sec] = 1.0 / constants.Av;
            H_tot_sec[sec] = constants.incepient_soot_H / constants.Av;
            
        soot_array = np.hstack((N_agg_sec, N_pri_sec, H_tot_sec));
    
        ## Non-dimensionlized volume & number
        d_0 = 2e-9;
        v_p = (constants.Pi/6.0*d_0**3.0);
        v_agg_sec = np.zeros((n_secs,));
        v_agg_sec_diff = np.zeros((n_secs,));
        v_agg_sec[0] = v_p;
        for i in range(1, n_secs):
            v_agg_sec[i] = v_agg_sec[i-1] * particle_dynamics_model.spacing_factor;
            
        v_agg_sec_diff[:-1] = np.diff(v_agg_sec);
        v_agg_sec_diff[-1]= v_agg_sec_diff[-2]
            
    reactor.user_defined_initial_soot = soot_array;
    # The collision occurs in the free molecular regime, so the collision frequency interpolation type is set to "free_molecule"
    particle_dynamics_model.beta_interp_method_type = "free_molecule"
    ## Reactor integrator settings
    reactor.max_step = 2e-4;
    reactor.max_time = max_res_time;
    ## The energy equations is not solved, so the temperature is kept constants
    reactor.temperature_solver_type = "isothermal";
    ## Starting the reactor
    reactor.start();

    data_headers = ["t[s]", "lambda[m]", "N_agg[mol/kg]", "N_agg[#/m3]", "N_pri[mol/kg]", "N_pri[#/m3]", "n_p", "d_p[nm]", "d_m[nm]", "d_g[nm]"];
    data_list = defaultdict(list);
    sectional_variables = ["N_agg", "N_pri", "d_p", "d_m", "d_g"];
    sectional_data_list = defaultdict(list);
    sectional_headers = reduce(operator.add, [[f"{variable}_{i}" for i in range(particle_dynamics_model.n_secs)] for variable in sectional_variables])
    eta_psi_headers = reduce(operator.add, [[f"{variable}_{i}" for i in range(particle_dynamics_model.n_secs)] for variable in ["eta", "psi"]])
    extra_sectional_headers = ["d_mg[nm]", "sigma_mg"];
    
    def calc_sigma_mg(N_agg_sec, d_m_sec):
        # d_m_sec = sectional_data_list["d_m"][-1];
        # N_agg_sec = sectional_data_list["N_agg"][-1];
        d_mg = np.exp(
            np.sum(np.log(d_m_sec)*N_agg_sec)/np.sum(N_agg_sec)
        );
        log_sigma_cr = (
            (
                np.sum(
                    N_agg_sec * 
                    (
                        (np.log(d_m_sec/d_mg))**2.0
                    )
                )/(np.sum(N_agg_sec))
            ) ** 0.5
        );
        return d_mg, np.exp(log_sigma_cr);
    
    def calc_eta_psi(N_agg):
        v_tot_sec = N_agg * v_agg_sec;
        v_tot_avg = np.sum(v_tot_sec) / np.sum(N_agg)
        eta = v_agg_sec / v_tot_avg;
        psi = N_agg / v_agg_sec_diff * v_tot_avg / np.sum(N_agg);
    
        return eta.tolist(), psi.tolist()
    
    def append_common():
        data_list["t[s]"].append(reactor.restime);
        data_list["lambda[m]"].append(soot_gas.lambda_gas);
        data_list["N_agg[mol/kg]"].append(soot.N_agg);
        data_list["N_pri[mol/kg]"].append(soot.N_pri);
        data_list["N_agg[#/m3]"].append(soot.N_agg * constants.Av * gas_density);
        data_list["N_pri[#/m3]"].append(soot.N_pri * constants.Av * gas_density);
        data_list["n_p"].append(soot.n_p);
        data_list["d_p[nm]"].append(soot.d_p * 1e9);
        data_list["d_m[nm]"].append(soot.d_m * 1e9);
        data_list["d_g[nm]"].append(soot.d_g * 1e9);
    
    def append_sectional():
        for variable in sectional_variables:
            variable_secs = [getattr(particle_dynamics_model, f"{variable}_sec")(i) for i in range(particle_dynamics_model.n_secs)];
            sectional_data_list[variable].append(variable_secs);
        
        d_mg, sigma_mg = calc_sigma_mg(np.array(sectional_data_list["N_agg"][-1]), np.array(sectional_data_list["d_m"][-1]));
        sectional_data_list["sigma_mg"].append(sigma_mg);
        sectional_data_list["d_mg[nm]"].append(d_mg*1e9);
    
        eta, psi = calc_eta_psi(np.array(sectional_data_list["N_agg"][-1]));
        sectional_data_list["eta"].append(eta);
        sectional_data_list["psi"].append(psi);    
        
    
    append_call_list = [append_common];
    if particle_dynamics_model_type == "Sectional":
        append_call_list += [append_sectional]
    
    ## Adding initial values to the output data
    append_output = [append() for append in append_call_list];
    
    step = 0;
    ## Running the reactor
    start_time = time.time()
    while reactor.restime < max_res_time:
        step += 1;
        reactor.step();
        if step % 5 == 0:
            append_output = [append() for append in append_call_list];
    end_time = time.time();
    print(f"The simulation time is {(end_time-start_time):0.2f} s");

    ## Building the dataframe & writting the csv file
    all_data = np.array([data_list[header] for header in data_headers]).transpose();
    columns = [];
    columns += data_headers;
    if particle_dynamics_model_type == "Sectional":
        sectional_data = np.hstack(tuple(sectional_data_list[variable] for variable in sectional_variables));
        extra_sectional_data = np.vstack(tuple(sectional_data_list[headers] for headers in extra_sectional_headers));
        extra_sectional_data = extra_sectional_data.transpose();
        eta_psi_data = np.hstack(tuple(np.array(sectional_data_list[variable]) for variable in ["eta", "psi"]))
        
        all_data = np.hstack((all_data, sectional_data, extra_sectional_data, eta_psi_data));
        columns += sectional_headers + extra_sectional_headers + eta_psi_headers;
    df = pd.DataFrame(data = all_data, columns = columns);
    results_dir = "results";
    if not os.path.exists(results_dir):
        os.makedirs(results_dir);
    df.to_csv(f"{results_dir}/{particle_dynamics_model_type}.csv", index = False)