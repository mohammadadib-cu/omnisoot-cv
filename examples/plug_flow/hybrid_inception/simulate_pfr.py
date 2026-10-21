# --------------------------------------------------------------------------
# This example demonstrates the HybridDimerization model in Omnisoot (>0.3.2).
#
# HybridDimerization allows users to combine ReactiveDimerization,
# IrreversibleDimerization, DimerCoalescence, and EBridgeModified as
# submodels, each with independent precursors and inception and adsorption
# prefactors.
#
# The example also highlights several new features:
#   - Access to Omnisoot's data directory for shared files.
#   - Convenient and efficient access to simulation data through
#     reactor.solution_data_dict, which can be used to initialize and append
#     data to Cantera's SolutionArray.
#   - Configuration of soot submodels using a YAML file, reducing script
#     complexity and allowing settings to be reused across simulations.
#   - Initializing reactors using Cantera's Solution object instead of SootGas. 
# --------------------------------------------------------------------------
# Declarations
import time
import os
import cantera as ct
import omnisoot as oms
import numpy as np
import pandas as pd
# Simulation parameters
T5_K = 2300;
P5_atm = 5;
inlet_composition = {"CH4":0.2, "N2": 0.8};
reactor_inlet_mdot = 0.005; # kg/s
reactor_length = 1; # m
reactor_area = 3.14 * 15e-2**2; # m3
# Creating the solution object and passing it to SootGas which is wrapper connecting chemistry to soot model
gas = ct.Solution(oms.data_dir/"Caltech.yaml")
# Configuring the reactor
reactor = oms.PlugFlowReactor(gas);
reactor.inlet.TPX = T5_K, P5_atm*ct.one_atm, inlet_composition;
reactor.inlet.mdot = reactor_inlet_mdot;
reactor.inlet_area = reactor_area;
reactor.temperature_solver_type = "energy_equation"
# Configuring soot model
reactor.soot.configure_from_yaml("soot_config.yaml");
# Starting the reactor
reactor.start();
# Creating the states object to snapshot the solution vector during the simulation    
states = ct.SolutionArray(gas, 1, extra=reactor.solution_data_dict);
# Marching forward along the reactor until the integrated length reaches the reactor length
start_time = time.time()
step = 0
while reactor.length < reactor_length:
    reactor.step();
    step +=1;
    if step % 50 == 0:
        states.append(gas.state, **reactor.solution_data_dict)
else:
    states.append(gas.state, **reactor.solution_data_dict)
end_time = time.time();
print(f"The simulation time took {end_time-start_time} s");
# Building the dataframe and writing the results on the disk.
omnisoot_columns = reactor.solution_headers;
omnisoot_data = [getattr(states, key) for key in reactor.solution_keys]

species_columns = [f"{sp}" for sp in gas.species_names];
species_data = [states.X[:,i] for i in range(len(gas.species_names))];

data_array = (np.array(omnisoot_data + species_data)).transpose();
df = pd.DataFrame(data = data_array, columns = omnisoot_columns + species_columns);
df.to_csv("results.csv")