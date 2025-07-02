from simulate_psr import simulate_psr
import cantera as ct

# Mechanisms
MECHANISMS = {
    "Caltech":"../mechanisms/Caltech.yaml",
};
# Mechanism list
mech_name = "Caltech";
# Model lists
PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];
particle_dynamics_model_types = ["Monodisperse", "Sectional"];

# Cantera Solution object
gas = ct.Solution(MECHANISMS[mech_name]);
eq_ratio = 2
fuel = {"C2H4":1.0}
oxidizer = {"N2":0.79, "O2":0.21};
gas = ct.Solution("abf_1bar.yaml");
gas.set_equivalence_ratio(eq_ratio, fuel=fuel, oxidizer=oxidizer);
composition = gas.X;

argdict = dict(
        gas = gas,
        mech_name = mech_name, 
        reactor_volume = 250e-6, 
        pressure = 101325, 
        residence_time = 10e-3, 
        reactor_temperature = 1800,
        inlet_temperature = 300,
        composition = composition,
        PAH_growth_model_type = PAH_growth_model_types[0],
        particle_dynamics_model_type = particle_dynamics_model_types[0],
        precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"],
        dimcoal_effs = [0.002, 0.015, 0.025, 0.004],
        solver_type = "LSODA",
);

for particle_dynamics_model_type in particle_dynamics_model_types:
    argdict["particle_dynamics_model_type"] = particle_dynamics_model_type;
    for PAH_growth_model_type in PAH_growth_model_types:
        argdict["PAH_growth_model_type"] = PAH_growth_model_type;
        simulate_psr(**argdict)