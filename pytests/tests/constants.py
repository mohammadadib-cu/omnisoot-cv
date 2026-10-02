from pathlib import Path

MECHANISM_PATH = Path(__file__).parent.parent / "data" / "Caltech.yaml"

PRECURSOR_NAMES = ["A2", "A3", "A4", "A2R5", "A3R5", "A4R5"];

INITIAL_COMPOSITION = {"CH4": 0.1, "AR": 0.9};

NUMBER_OF_TEST_STEPS = 10;

TEMPERATURE_PROFILE_TIME = [
    [0, 1e-3, 3e-3, 6e-3, 9e-3, 15e-3, 20e-3],
    [500, 800, 1500, 1800, 1800, 1200, 1000]
];

INITIAL_D_P = 2.0e-9;

SMALL = 1.0e-80;
CONSERVED_ATOL = 1e-18;
CONSERVED_RTOL = 1e-10;

HOT_GAS_COMPOSITION = {
    'N2': 0.9,
    'C2H2' : 0.05,
    'H' : 5e-5,
    'H2' : 0.04968,
    'A1' : 1e-4,
    'A2' : 3e-5,
    'A3' : 3e-5,
    'A4' : 3e-5,
    'A2R5' : 3e-5,
    'A3R5' : 1e-5,
    'A4R5' : 3e-5,
}

HOT_GAS_COMPOSITION_WITH_OXIDIZERS = {
    'N2': 0.599,
    'O2': 0.1,
    'H2O': 0.1,
    'CO2': 0.1,
    'OH' : 1e-3,
    'C2H2' : 0.05,
    'H' : 5e-5,
    'H2' : 0.04968,
    'A1' : 1e-4,
    'A2' : 3e-5,
    'A3' : 3e-5,
    'A4' : 3e-5,
    'A2R5' : 3e-5,
    'A3R5' : 1e-5,
    'A4R5' : 3e-5,
}

PREDEFINED_MONODISPERSE_SOOT_VAR_LIST = [
    5.76E-06*100, # N_agg
    2.42E-05*1000, # N_pri
    0.882520699*10, # C_tot
    0.509249656*10, # H_tot 
]

N_SECS = 45

PREDEFINED_SECTIONAL_SOOT_VAR_LIST = (
    [5.76E-06*10]*N_SECS+ # N_agg
    [2.42E-05*100]*N_SECS+ # N_pri
    [0.509249656]*N_SECS # H_tot 
)

PREDEFINED_SOOT_VAR_LIST = {
    "Monodisperse" : PREDEFINED_MONODISPERSE_SOOT_VAR_LIST,
    "Sectional" : PREDEFINED_SECTIONAL_SOOT_VAR_LIST,
}

STANDALONE_PARTICLE_DYNAMIC_MODEL_TYPES = [
    "Monodisperse",
    "Sectional",
]


ALL_PAH_GROWTH_MODEL_TYPES = [
    "ReactiveDimerization",
    "IrreversibleDimerization",
    "EBridgeModified",
    "DimerCoalescence",
    "HybridDimerization",
]

ALL_SURFACE_REACTIONS_MODEL_TYPES = [
    "FrenklachHACA",
    "FrenklachHACAModified",
    "BlanquartHACA",
]

STANDALONE_PAH_GROWTH_MODEL_TYPES = [
    "ReactiveDimerization",
    "IrreversibleDimerization",
    "EBridgeModified",
    "DimerCoalescence",
]

STANDALONE_PARTICLE_DYNAMIC_MODEL_TYPES = [
    "Monodisperse",
    "Sectional",
]

TEST_EXECUSION_ORDER = [
    "test_soot_gas",
    "test_particle_dynamics",
    "test_reactors",
]