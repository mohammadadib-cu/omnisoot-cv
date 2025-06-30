# Declarations
import warnings
from simulate_PSR_PFR import simulate_PSR_PFR

# Suppress warnings
warnings.filterwarnings("ignore")

def main():
    # Main function to run the simulation
    PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence", "EBridgeModified", "IrreversibleDimerization"];
    eq_ratios = [1.9, 2.0, 2.1];

    prefactors_dict = {
        "ReactiveDimerization" : {
            "adsorption_prefactor" : 1e-2,
            "inception_prefactor" : 5e-4,
        },
        "DimerCoalescence" : {
            "adsorption_prefactor" : 1e-2,
            "inception_prefactor" : 1e-3,
        },
        "EBridgeModified" : {
            "adsorption_prefactor" : 1e-2,
            "inception_prefactor" : 2e-4,
        },
        "IrreversibleDimerization" : {
            "adsorption_prefactor" : 1,
            "inception_prefactor" : 1e-2,
        },
    }


    arg_dict = dict(
        mech_name = "KAUST",
        PAH_growth_model_type = PAH_growth_model_types[0],
        particle_dynamics_model_type = ["Monodisperse", "Sectional"][1],
        precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"],
        pfr_length = 0.7, #m
        pfr_diameter = 5.1e-2, #m
        psr_volume = 250*1e-6, #m3
        psr_residence_time = 11e-3, #s
        eq_ratio = 1.9,
        fuel = {"C2H4":1.0},
        oxidizer = {"N2":0.79, "O2":0.21},
        psr_inlet_temperature = 300,
        pressure = 101325,
        mdot_corrector = 1,
        use_constant_temperature = False,
        dimcoal_effs = [0.002, 0.015, 0.025, 0.004],
        inception_prefactor = 1,
        adsorption_prefactor = 1,
        solver_type = "LSODA",
    )


    for eq_ratio in eq_ratios:
        arg_dict["eq_ratio"] = eq_ratio;
        for PAH_growth_model_type in PAH_growth_model_types:
            arg_dict["PAH_growth_model_type"] = PAH_growth_model_type;
            arg_dict["inception_prefactor"] = prefactors_dict[PAH_growth_model_type]["inception_prefactor"];
            arg_dict["adsorption_prefactor"] = prefactors_dict[PAH_growth_model_type]["adsorption_prefactor"];
            simulate_PSR_PFR(**arg_dict)

if __name__ == "__main__":
    main();