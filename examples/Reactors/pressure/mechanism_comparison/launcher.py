import numpy as np
from simulate_nosoot import simulate_nosoot
from plot import plot

# To suppress the long errors of loading mechanism in Cantera
import warnings
warnings.filterwarnings('ignore');

mech_names = ["Caltech", "KAUST", "ABF1bar", "ITV", "CRECK", "NUIG", "FFCM2"];
arg_dict = dict(
        mech_name = "",
        T5_K = 2200,
        P5_atm = 4.5,
        initial_composition = {"CH4": 0.05, "AR": 0.95},
        max_res_time = 5.1e-3, 
        results_dir = "results",
    )


for mech_name in mech_names:
    arg_dict["mech_name"] = mech_name;
    simulate_nosoot(**arg_dict)


plot();