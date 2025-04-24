import time
import cantera as ct
from plot import plot
from simulate_shocktube import simulate_shocktube

MECH_PATH = {
    "Caltech":"./data/Caltech.yaml",
}

mech_name = "Caltech";
results_dir = "results";
particle_dynamics_model_type = "Monodisperse";
PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];
initial_composition = {"CH4":0.05, "AR":0.95};
T5_K = 2200; # K
P5_atm = 4.5; # atm
max_res_time = 100e-3; # s
precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"];

gas = ct.Solution(MECH_PATH[mech_name])

t0 = time.time();
for PAH_growth_model_type in PAH_growth_model_types:
    simulate_shocktube(
        mech_name = mech_name, 
        gas = gas,
        PAH_growth_model_type = PAH_growth_model_type, 
        particle_dynamics_model_type = particle_dynamics_model_type,
        T5_K = T5_K, 
        P5_atm = P5_atm, 
        initial_composition = initial_composition, 
        max_res_time = max_res_time, 
        precursors = precursors,
        inception_prefactor = 0.1,
        adsorption_prefactor = 1,
        results_dir = results_dir,
    );
t1 = time.time();
print(f"The simulation took {(t1-t0):0.2f} s")


plot(
    mech_name = mech_name, 
    PAH_growth_model_types = PAH_growth_model_types, 
    particle_dynamics_model_type = particle_dynamics_model_type,
    results_dir = results_dir,
)