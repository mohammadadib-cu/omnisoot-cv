from collections import namedtuple
import warnings

import cantera as ct
import numpy as np
import pandas as pd

from simulate_flowreactor import simulate_flowreactor

warnings.filterwarnings('ignore')

MECH_PATH = {
    "KAUST":"data/KAUST.yaml",
};

PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];
particle_dynamics_model_type = ["Monodisperse", "Sectional"][1];
TemperatureProfile = namedtuple("TemperatureProfile",["z","T"]);

df = pd.read_csv(f"data/temperature_profile.csv");
TemperatureProfile = namedtuple("TemperatureProfile",["z","T"]);

gas = ct.Solution(MECH_PATH["KAUST"]);
gas.TPX = 300, 101325, {'C2H4':0.006, 'N2':0.994};
mdot_1lmin = (1e-3 / 60) * gas.density

fl_categories = ["8.0lmin", "11.0lmin", "12.0lmin"];


prefactors_dict = {
    "ReactiveDimerization" : {
        "adsorption_prefactor" : 0.001,
        "inception_prefactor" : 0.01,
    },
    "DimerCoalescence" : {
        "adsorption_prefactor" : 0.01,
        "inception_prefactor" : 0.004,
    },
    "EBridgeModified" : {
        "adsorption_prefactor" : 0.1,
        "inception_prefactor" : 1e-2,
    },
    "IrreversibleDimerization" : {
        "adsorption_prefactor" : 1,
        "inception_prefactor" : 0.08,
    }, 
}

arg_dict = dict(
    fl_category = "",
    gas = gas,
    mech_name = "KAUST",
    PAH_growth_model_type = PAH_growth_model_types[0],
    particle_dynamics_model_type = particle_dynamics_model_type,
    precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"],
    inlet_composition = {'C2H4':0.006, 'N2':0.994},
    mdot = mdot_1lmin,
    P = 101350,
    temperature_profile = None,
    reactor_length = 1.4,
    reactor_diameter = 16e-3,
    adsorption_prefactor = 1,
    inception_prefactor = 1,
    solver_type = "BDF",
);

fl_correction_dict = {
    "8.0lmin" : 0.82, 
    "11.0lmin": 1.7, 
    "12.0lmin" : 2.4,
};

for fl_category in fl_categories:
    arg_dict["fl_category"] = fl_category;
    fl_num = float(fl_category[:-4]);
    arg_dict["mdot"] = mdot_1lmin * fl_num * fl_correction_dict[fl_category]
    arg_dict["temperature_profile"] = TemperatureProfile(z = df["z[m]"].to_numpy(), T = df[fl_category].to_numpy());
    for PAH_growth_model_type in PAH_growth_model_types:
        arg_dict["adsorption_prefactor"] = prefactors_dict[PAH_growth_model_type]["adsorption_prefactor"]
        arg_dict["inception_prefactor"] = prefactors_dict[PAH_growth_model_type]["inception_prefactor"]
        arg_dict["PAH_growth_model_type"] = PAH_growth_model_type;
        simulate_flowreactor(**arg_dict)