import cantera as ct
import numpy as np
from simulate_reactor import simulate
# from plot import plot

mech_name = "Caltech";
particle_dynamics_model_types = ["Monodisperse", "Sectional"];
PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];

arg_dict = dict(
        mech_name = mech_name,
        PAH_growth_model_type = "",
        particle_dynamics_model_type = "",
        P5 = 4.64 * ct.one_atm,
        T5 = 2355,
        initial_composition = {"CH4":0.1, "AR":0.9},
        max_res_time = 10e-3, 
        precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"],
        solver_type = "LSODA",
    )


for particle_dynamics_model_type in particle_dynamics_model_types:
    arg_dict["particle_dynamics_model_type"] = particle_dynamics_model_type;
    if particle_dynamics_model_type == "Monodisperse":
        arg_dict["solver_type"] = "BDF";
    else:
        arg_dict["solver_type"] = "LSODA";
    for PAH_growth_model_type in PAH_growth_model_types:
        arg_dict["PAH_growth_model_type"] = PAH_growth_model_type;
        simulate(**arg_dict)

# plot();