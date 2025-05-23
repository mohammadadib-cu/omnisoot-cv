from run_coagulation import run_coagulation
from plot import plot

particle_dynamics_models = ["Monodisperse", "Sectional"]

for particle_dynamics_model in particle_dynamics_models:
    run_coagulation(particle_dynamics_model);

plot();