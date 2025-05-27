from collections import namedtuple
import pandas as pd
from simulate_psr import simulate_psr

PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];
particle_dynamics_model_types = ["Monodisperse", "Sectional"];
energy_equation_modes = ["fixed_temperature", "adiabatic"]


arg_dict = dict(
    mech_name = "Caltech",
    PAH_growth_model_type = PAH_growth_model_types[0],
    particle_dynamics_model_type = particle_dynamics_model_types[0],
    energy_equation_mode = energy_equation_modes[0],
    precursors = ["A2", "A3", "A4", "A2R5"],
    inlet_composition = {"C2H4":0.174, "O2": 0.174, "N2": 0.652},
    reactor_pressure = 101350,  # Pa
    reactor_temperature = 1800, # K
    reactor_volume = 250e-6,    # m3
    residence_time = 8.5e-3,    # s
);

# Simulate
for PAH_growth_model_type in PAH_growth_model_types:
    arg_dict["PAH_growth_model_type"] = PAH_growth_model_type;
    simulate_psr(**arg_dict);
    
