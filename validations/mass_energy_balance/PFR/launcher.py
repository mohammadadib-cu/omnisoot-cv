import cantera as ct
import numpy as np
from simulate import simulat_reactor

MECH_PATH = {
    "Caltech":"../mechanisms/Caltech.yaml",
}

mech_name = "Caltech";
particle_dynamics_model_types = ["Monodisperse", "Sectional"][1:];
PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"][2:];

arg_dict = dict(
        gas = ct.Solution(MECH_PATH[mech_name]),
        mech_name = mech_name,
        PAH_growth_model_type = "",
        particle_dynamics_model_type = "",
        precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"],
        inlet_composition = {"CH4":0.1, "AR":0.9},
        T_inlet = 2200,
        P = 101350,
        mdot = 0.0001,
        reactor_length = 0.8, 
        reactor_diameter = 11e-3,
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
        simulat_reactor(**arg_dict)

# plot();