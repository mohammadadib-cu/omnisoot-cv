# This is an example of using the PlugFlowReactor model of omnisoot to simulate soot formation during ethylene pyrolysis in a flow reactor 
# The temperature profile and composition is based on Araki et al. (https://doi.org/10.1016/j.fuproc.2020.106673)

import time
import cantera as ct
import numpy as np
import pandas as pd
from simulate_flowreactor import simulate_flowreactor

mech_names = "Caltech";
PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence", "EBridgeModified", "IrreversibleDimerization"];
particle_dynamics_model_type = "Monodisperse";
solver_types = ["BDF", "LSODA", "Radau"];

df = pd.read_csv(f"data/T_1700K.csv");
temperature_profile = np.vstack((df["z[m]"].to_numpy(), df["T[K]"].to_numpy()));

arg_dict = dict(
    mech_name = "Caltech",
    PAH_growth_model_type = PAH_growth_model_types[0],
    particle_dynamics_model_type = particle_dynamics_model_type,
    precursors = ["A2", "A3R5", "A4R5", "A2R5"],
    inlet_composition = {"C2H4":0.01, "N2": 0.99},
    mdot = 3.793335e-05,
    P = 101350,
    #T_in = 2200,
    temperature_profile = temperature_profile,
    reactor_length = 0.64,
    reactor_diameter = 11e-3,
    solver_type = solver_types[0]
)

for solver_type in solver_types:
    t0 = time.time();
    arg_dict["solver_type"] = solver_type;
    for PAH_growth_model_type in PAH_growth_model_types:
        arg_dict["PAH_growth_model_type"] = PAH_growth_model_type;
        simulate_flowreactor(**arg_dict);
    t1 = time.time();
    print(f"The simulation finished in {(t1-t0):0.2f} seconds using {solver_type} solver.\n")