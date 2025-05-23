from collections import namedtuple
import pandas as pd
from simulate_flowreactor import simulate_flowreactor
from plot import plot

PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];
particle_dynamics_model_types = ["Monodisperse", "Sectional"];
df = pd.read_csv("./data/temperature_profile.csv");
TemperatureProfile = namedtuple("TemperatureProfile",["z","T"]);
temperature_profile = TemperatureProfile(z = df["z[cm]"].to_numpy(), T = df["T[K]"].to_numpy());

arg_dict = dict(
    mech_name = "Caltech",
    PAH_growth_model_type = PAH_growth_model_types[0],
    particle_dynamics_model_type = particle_dynamics_model_types[0],
    precursors = ["A2", "A3", "A4", "A2R5"],
    inlet_composition = {"C2H4":0.01, "N2": 0.99},
    mdot = 0.39915,
    P = 101350,
    temperature_profile = temperature_profile,
    reactor_length = 0.64
);

# Simulate
for particle_dynamics_model_type in particle_dynamics_model_types:
    arg_dict["particle_dynamics_model_type"] = particle_dynamics_model_type;
    simulate_flowreactor(**arg_dict);
    
# Plot
plot(particle_dynamics_model_types,  PAH_growth_model_types[0]);