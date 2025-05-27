# Declarations
import os
import time
import cantera as ct
import numpy as np
import pandas as pd
from omnisoot import PerfectlyStirredReactor, SootGas, SootThermo
# Constants
Av = 6.0221408e+23; #1/mol
MW_carbon = 12.011e-3 # [kg/mol]
MW_hydrogen = 1.008e-3 # [kg/mol]
MECH_PATH = {
    "Caltech":"./data/Caltech.yaml",
}

MAX_SIMULATION_TIME = 3;
APPEND_PER_STEPS = 20;
# -------------------------------------------------------------------------
# This function simulates the reacting flow using PSR reactor of omnisoot and
# writes the simulation results on the disk
# -------------------------------------------------------------------------

def simulate_psr(
    mech_name,
    PAH_growth_model_type,
    particle_dynamics_model_type,
    energy_equation_mode,
    precursors,
    inlet_composition,
    reactor_pressure,
    reactor_temperature,
    reactor_volume,
    residence_time,
):
    # Creating the solution object and passing it to SootGas which is wrapper connecting chemistry to soot model
    gas = ct.Solution(MECH_PATH[mech_name]);
    gas.TPX = reactor_temperature, reactor_pressure, inlet_composition
    soot_gas = SootGas(gas);
    # Configuring the reactor
    psr = PerfectlyStirredReactor(soot_gas);
    if energy_equation_mode == "fixed_temperature":
        psr.temperature_solver_type = "isothermal";
    elif energy_equation_mode == "adiabatic":
        psr.temperature_solver_type = "energy_equation";

    psr.mdot_in = gas.density * reactor_volume / residence_time;
    psr.reactor_volume = reactor_volume;
    psr.P_outlet = ct.one_atm;
    psr.max_step = 1e-3;

    # Configuring soot model
    soot = psr.soot;
    soot.particle_dynamics_model_type = particle_dynamics_model_type;
    soot.PAH_growth_model_type = PAH_growth_model_type;
    soot.set_precursor_names(precursors);
    # Getting the handle to particle dynamics model
    particle_dynamics = soot.particle_dynamics_model
    # Starting the reactor
    psr.start();
    # Creating the states object to snapshot the solution vector during the simulation
    states = ct.SolutionArray(gas, 1, extra={
                                         'restime' : [psr.restime],
                                         'mdot' : [psr.mdot_in],
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
                                         'C_tot_inception': [soot.C_tot_inception],
                                         'C_tot_PAH_adsorption': [soot.C_tot_PAH],
                                         'C_tot_HACA_growth': [soot.C_tot_growth],
                                         'C_tot_HACA_oxidation': [soot.C_tot_oxidation],
                                         'coagulation_mass': [soot.coagulation_mass],
                                         'total_carbon_flux' : [psr.total_carbon_flux_in],
                                         'total_hydrogen_flux' : [psr.total_hydrogen_flux_in],
                                         'soot_carbon_flux' : [psr.soot_carbon_flux_in],
                                         'soot_mass_flux' : [psr.soot_mass_flux_in],
                                        });
    
    # Funtion to Append Data to States
    def append():
        states.append(
            gas.state,
            restime=psr.restime,
            mdot = psr.mdot_out,
            N_agg = soot.N_agg,
            N_pri = soot.N_pri,
            C_tot = soot.C_tot,
            H_tot = soot.H_tot,
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
            C_tot_inception=  soot.C_tot_inception,
            C_tot_PAH_adsorption = soot.C_tot_PAH,
            C_tot_HACA_growth = soot.C_tot_growth,
            C_tot_HACA_oxidation = soot.C_tot_oxidation,
            coagulation_mass = soot.coagulation_mass,
            total_carbon_flux = psr.total_carbon_flux_out,
            total_hydrogen_flux = psr.total_hydrogen_flux_out,
            soot_carbon_flux = psr.soot_carbon_flux_out,
            soot_mass_flux = psr.soot_mass_flux_out,
        );

    # Marching forward along the reactor until the integrated length reaches the reactor length
    start_time = time.time()
    step = 0
    while psr.restime < MAX_SIMULATION_TIME:
        psr.step();
        step +=1;
        if step % APPEND_PER_STEPS == 0:
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

    # Gas Flow Variables
    flow_columns += ["t[s]", "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "mdot[kg/s]"];
    flow_data += [states.restime, states.T, states.density, states.mean_molecular_weight/1000, states.mdot];

    # Soot Variables
    soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
    soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

    soot_columns += ["A_tot[m2/kg]", "N_agg[#/m3]", "N_pri[#/m3]"]
    soot_data += [states.A_tot, states.N_agg * Av * states.density, states.N_pri * Av * states.density];

    soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "d_v[nm]", "n_p"]
    soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.d_v*1e9, states.n_p];

    soot_columns += ["soot_mass[kg/kg]", "soot_volume_fraction[-]", "SSA[m2/g]",
                     "total_carbon_flux[kg/m2-s]", "carbon_yield[-]"];
    soot_data += [states.total_soot_mass, states.volume_fraction, states.SSA,
                  states.total_carbon_flux, states.soot_carbon_flux/states.total_carbon_flux[0]];

    soot_columns += ['inception_flux[#/m3-s]'];
    soot_data += [states.inception_mass * Av * states.density];

    soot_columns += [
        'C_tot_inception[mol/m3-s]', 'C_tot_PAH_adsorption[mol/m3-s]',
        'C_tot_HACA_growth[mol/m3-s]', 'C_tot_HACA_oxidation[mol/m3-s]',
    ];
    soot_data += [
        states.C_tot_inception * states.density, states.C_tot_PAH_adsorption * states.density,
        states.C_tot_HACA_growth * states.density, states.C_tot_HACA_oxidation * states.density,
    ];    

    
    # Species Variables
    species_columns = [f"{sp}" for sp in gas.species_names];
    species_data += [states.X[:,i] for i in range(len(gas.species_names))];
    
    # Mass Balance
    C_flux_in = psr.mdot_in * (states.elemental_mass_fraction('C')[0]+states.C_tot[0]*MW_carbon);
    C_flux_out = (states.elemental_mass_fraction('C')+states.C_tot*MW_carbon)*states.mdot;
    C_flux_error = (C_flux_out-C_flux_in)/C_flux_in;

    H_flux_in = psr.mdot_in * states.elemental_mass_fraction('H')[0];
    H_flux_out = (states.elemental_mass_fraction('H')+states.H_tot*MW_hydrogen)*states.mdot;
    H_flux_error = abs(H_flux_out-H_flux_in)/H_flux_in;
    
    # Energy Balance
    H1 = psr.mdot_in * psr.h_in;
    H2 = states.mdot * states.enthalpy_mass;
    Hp = (states.total_soot_mass*states.mdot) * np.array(list(SootThermo.h_mass_soot(T, P) for T, P in zip(states.T, states.P)));
    H_error = (H2-H1+Hp)/H1;
    
    error_columns = ["C_mass_error", "H_mass_error", "enthalpy_error"];
    error_data = [C_flux_error, H_flux_error, H_error]

    species_columns = [f"{sp}" for sp in gas.species_names];
    species_data = [states.X[:,i] for i in range(len(gas.species_names))];
    
    columns = flow_columns + soot_columns + species_columns + error_columns;
    data = (np.array(flow_data + soot_data + species_data + error_data)).transpose();
    df = pd.DataFrame(data = data, columns = columns);
    
    # Creating the results directory and storing the dataframe on disk
    output_dir = f"results/{particle_dynamics_model_type}/{PAH_growth_model_type}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);

    file_name = f"{output_dir}/sim_results.csv";
    df.to_csv(file_name, index=False);
    print(f"{file_name} was written!")