###################################################################################################################
#  This is an example use-case of omnisoot to simulate soot formation during methane pyrolysis in a constant volume
#  reactor that allows using three different models to calculates surface reactivity (alpha): 
#  - 'constant' value 
#  - the 'empirical' equation by Appel et al. (https://doi.org/10.1016/S0010-2180(99)00135-2) 
#  - H/C ratio that relates soot surface reactivity to soot 'composition'
###################################################################################################################
import itertools
from simulate import simulate
# To suppress warning due to loading KAUST mechanism
import warnings
warnings.filterwarnings('ignore')

mech_names = ["KAUST"];
particle_dynamics_model_types = ["Monodisperse", "Sectional"][:1];
PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];
# Number of section in the sectional model
number_of_sections = 60;
# The model used to describe surface reactivity (alpha) in HACA 
surface_reactivity_model_type = ["empirical", "constant", "composition"][0];
# Carbonization removes H from soot particles to describe maturity
# Omnisoot was designed to accomodate different carbonization models,
# but currently only a simple Arrhenius equation is implemented
carbonization_enabled = True;
# Collision efficiency of particles are 1 by default. It can also be set to take the repulsion effect
# for small particles. The implementation is based on Eq.(16) of Hou et al.(https://doi.org/10.1016/j.jaerosci.2019.105478)
eta_coag_type = ["unity", "repulsion_d_based"][0];

extra_props = dict(
    number_of_sections = number_of_sections,
    surface_reactivity_model_type = surface_reactivity_model_type,
    surface_reactivity_constant = 0.01,
    carbonization_enabled = carbonization_enabled,
    eta_coag_type = eta_coag_type
)

main_props = dict(
    mech_name = "KAUST",
    PAH_growth_model_type = PAH_growth_model_types[0],
    particle_dynamics_model_type = PAH_growth_model_types[0],
    precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"],
    max_res_time = 10e-3,
    initial_composition = {"CH4": 0.1, "AR": 0.9},
    initial_temperature = 1800, # [K]
    initial_pressure = 101350, # [Pa]
)

# Running the simulations
for iterdata in itertools.product(particle_dynamics_model_types, PAH_growth_model_types):
    main_props["particle_dynamics_model_type"], main_props["PAH_growth_model_type"] = iterdata;
    simulate(**main_props, **extra_props);