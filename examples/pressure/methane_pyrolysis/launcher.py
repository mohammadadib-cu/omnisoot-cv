import cantera as ct

from C_tot_processor import main as C_tot_processor_main
from simulate import simulate

mech_dict = {
    "Caltech":"./data/Caltech.yaml",
}

P5_const = 0.002385;
mech_name = "Caltech";
particle_dynamics_model_type = ["Monodisperse", "Sectional"][1];
PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];

arg_dict = dict(
        gas = ct.Solution(mech_dict[mech_name]),
        mech_name = mech_name,
        PAH_growth_model_type = "",
        particle_dynamics_model_type = particle_dynamics_model_type,
        T5_K = 0,
        P5_atm = 0,
        initial_composition = {"CH4":0.05, "AR":0.95},
        max_res_time = 1.6e-3, 
        precursors = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"],
        inception_prefactor = 1,
        adsorption_prefactor = 1,
        results_dir = "results",
        solver_type = "BDF",
    )

prefactors_dict = {
    "ReactiveDimerization" : {
        "adsorption_prefactor" : 0.17,
        "inception_prefactor" : 0.17,
    },
    "DimerCoalescence" : {
        "adsorption_prefactor" : 0.21,
        "inception_prefactor" : 0.21,
    },
    "EBridgeModified" : {
        "adsorption_prefactor" : 0.36,
        "inception_prefactor" : 0.36,
    },
    "IrreversibleDimerization" : {
        "adsorption_prefactor" : 2.3,
        "inception_prefactor" : 2.3,
    },
};
T5_Ks = [1800, 1900, 2000, 2050, 2100, 2125, 2150, 2175, 2200, 2225, 2250, 2275, 2300, 2325, 2350, 2375, 2400, 2450, 2500, 2550, 2600, 2650, 2700, 2800, 2900, 3000];

for PAH_growth_model_type in PAH_growth_model_types:
    arg_dict["PAH_growth_model_type"] = PAH_growth_model_type;
    arg_dict["adsorption_prefactor"] = prefactors_dict[PAH_growth_model_type]["adsorption_prefactor"];
    arg_dict["inception_prefactor"] = prefactors_dict[PAH_growth_model_type]["inception_prefactor"];
    for T5 in T5_Ks:
        P5 = T5 * P5_const * 1e5;
        arg_dict["T5_K"] = T5;
        arg_dict["P5_atm"] = P5 / 101325;
        simulate(**arg_dict);

C_tot_processor_main();

