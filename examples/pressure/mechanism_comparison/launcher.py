import cantera as ct
import numpy as np
from simulate import simulate
from simulate_nosoot import simulate_nosoot
from plot import plot_nosoot, plot_maxsoot

# Path to the mechanism files
MECH_DICT = {
    "Caltech":"./data/Caltech.yaml",
    "KAUST": "./data/KAUST.yaml",
    "ABF1bar" : "./data/abf_1bar.yaml",
    "CRECK" : "./data/CRECK_no_BINS.yaml",
    "ITV" : "./data/ITV_PAH.yaml",
    "NUIG" : "./data/NUIG.yaml",
    "FFCM2" : "./data/FFCM2.yaml",
}
# To suppress the long errors of loading mechanism in Cantera
import warnings
warnings.filterwarnings('ignore');

# No-soot simulation for different mechanisms
# Mechanisms are Caltech, KAUST, ABF1bar, ITV, CRECK, FFCM2
mech_names = ["Caltech", "KAUST", "ABF1bar", "ITV", "CRECK", "FFCM2"];
arg_dict = dict(
        gas = None,
        mech_name = "",
        T5_K = 2203,
        P5_atm = 4.5,
        initial_composition = {"CH4": 0.05, "AR": 0.95},
        max_res_time = 3.1e-3, 
        results_dir = "results",
    )

for mech_name in mech_names:
    arg_dict["gas"] = ct.Solution(MECH_DICT[mech_name]);
    arg_dict["mech_name"] = mech_name;
    simulate_nosoot(**arg_dict)

# Soot simulation for different mechanisms
# Mechanisms are Caltech, KAUST, ABF1bar, ITV, CRECK
precursors_dict = {
    "Caltech": ["A2", "A3", "A4", "A2R5"],
    "KAUST": ["A2", "A3", "A4", "A2R5"],
    "ABF1bar": ["A2", "A3", "A4", "A2R5"],
    "ITV": ["A2", "A3XC14H10", "A4XC16H10", "A2R5",],
    "CRECK": ["C10H8", "C14H10", "C16H10", "C12H8",],
}

mech_names = ["Caltech", "KAUST", "ABF1bar", "ITV", "CRECK"];
arg_dict = dict(
        gas = None,
        mech_name = "",
        PAH_growth_model_type = "IrreversibleDimerization",
        particle_dynamics_model_type = "Monodisperse",
        T5_K = 2203,
        P5_atm = 4.9,
        initial_composition = {"CH4":0.05, "AR": 0.95},
        max_res_time = 3.1e-3, 
        precursors = None, 
        results_dir = "results_Max",
        PAH_growth_model_config = None,
        solver_type = "BDF",
)

for mech_name in mech_names:
    gas = ct.Solution(MECH_DICT[mech_name]);
    arg_dict["mech_name"] = mech_name;
    arg_dict["gas"] = gas;
    precursors = precursors_dict[mech_name];
    arg_dict["precursors"] = precursors;
    arg_dict["PAH_growth_model_config"] = dict(
            stick_eff_dimerization =  [1.0] * len(precursors),
            stick_eff_adsorption = [1.0] * len(precursors),
    )
    simulate(**arg_dict)

plot_nosoot();
plot_maxsoot();