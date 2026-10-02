# ----------------------------------------------------------------------------------------------------------------------------------
# This example illustrates the use of the HybridDimerization model in Omnisoot (>0.3.0).
# HybridDimerization allows users to combine ReactiveDimerization, IrriversibleDimerization, DimerCoalescence, and EBridgeModified
# as submodels with independent precursors, and inception and adsorption prefactors.
# ----------------------------------------------------------------------------------------------------------------------------------

# Declarations
import cantera as ct
from omnisoot import PressureReactor, SootGas, SootThermo
import numpy as np
import pandas as pd
import os

# Constants
Av = 6.0221408e+23; #1/mol;

# Helper functions
def create_extra_dict(reactor):
    soot = reactor.soot
    extra_dict_base = dict(
        # Gas
        restime = reactor.restime,
        reactor_volume = reactor.reactor_volume,
        gas_volume = reactor.gas_volume,
        gas_mass = reactor.gas_mass,
        # Main Soot Variables
        N_pri = soot.N_pri,
        N_agg = soot.N_agg,
        C_tot = soot.C_tot,
        H_tot = soot.H_tot,
        # Derived Soot Variables
        A_tot = soot.A_tot,
        d_p = soot.d_p,
        d_m = soot.d_m,
        d_g = soot.d_g,
        n_p = soot.n_p,
        total_soot_mass = soot.total_mass,
        volume_fraction = soot.volume_fraction,
        SSA = soot.SSA,
        total_carbon_mass = reactor.total_carbon_mass,
        total_hydrogen_mass = reactor.total_hydrogen_mass,
        soot_carbon_mass = reactor.soot_carbon_mass
    );
    extra_dict = extra_dict_base.copy();
    return extra_dict;

def make_values_list(extra_dict):
    new_extra_dict = extra_dict.copy();
    for key in new_extra_dict.keys():
        new_extra_dict[key] = [new_extra_dict[key]];
        
    return new_extra_dict;

# Loading the reaction mechanism
gas = ct.Solution("data/Caltech.yaml")
# Simulation parameters
T5_K = 2203;
P5_atm = 4.5;
initial_composition = {"CH4":0.5, "AR": 0.95};
max_res_time = 5e-3;
results_dir = "results";

# Setting up the reactor
soot_gas = SootGas(gas);
soot_gas.TPX = T5_K, P5_atm * ct.one_atm, initial_composition;
### Building and configruing the reactor
reactor = PressureReactor(soot_gas);
reactor.pressure_type = "constant";
# Particle dynamics model is set to Monodisperse by default with reactor is initiated.
soot = reactor.soot;
# Setting the PAH growth model to 'HybridDimerization'
soot.PAH_growth_model_type = 'HybridDimerization';
# Exploring the hybrid model
PAH_growth_model = soot.PAH_growth_model;
print(f"Default submodels: [" + ", ".join([f"'{stn}'" for stn in PAH_growth_model.submodel_type_names]) + "]")
# - Output: 
# - Default submodels: ['ReactiveDimerization', 'DimerCoalescence', 'EBridgeModified', 'IrreversibleDimerization']
# Setting submodels to ReactiveDimerization and EBridgeModified
PAH_growth_model.submodel_type_names = ['ReactiveDimerization', 'EBridgeModified']
print(f"Current submodels: [" + ", ".join([f"'{stn}'" for stn in PAH_growth_model.submodel_type_names]) + "]");
# - Output: 
# - Current submodels: ['ReactiveDimerization', 'EBridgeModified']
# Submodels can be accessed through index or type name
reacdim_0 = PAH_growth_model.submodels[0];
reacdim_1 = PAH_growth_model.get_submodel_by_type_name('ReactiveDimerization');
print(f"reacdim_0 & reacdim_1 refer to the same object: {reacdim_0 == reacdim_1}")
# Setting precursors for each submodel
PAH_growth_model.get_submodel_by_type_name('ReactiveDimerization').precursor_names = ['A2', 'A3', 'A4'];
PAH_growth_model.get_submodel_by_type_name('EBridgeModified').precursor_names = ['A2R5', 'A3R5', 'A4R5'];

## Starting the reactor
reactor.start();
## SolutionArray to store simulation history
extra_dict = create_extra_dict(reactor);
states = ct.SolutionArray(gas, 1, extra = make_values_list(extra_dict));

### Running the simulation
step = 0;
while reactor.restime < max_res_time:
    reactor.step();
    step += 1;
    if step%10==0:
        extra_dict = create_extra_dict(reactor);
        states.append(gas.state, **extra_dict);


### Building the output data
soot_columns = [];
soot_data = [];
flow_columns = [];
flow_data = [];
species_columns = [];
species_data = [];
soot_columns += ["N_agg[mol/kg]", "N_pri[mol/kg]", "C_tot[mol/kg]", "H_tot[mol/kg]"];
soot_data += [states.N_agg, states.N_pri, states.C_tot, states.H_tot];

soot_columns += ["A_tot[m2/kg]", "N_agg[#/cm3]", "N_pri[#/cm3]", "N_agg[#/g]", "N_pri[#/g]"]
soot_data += [states.A_tot, states.N_agg*Av*states.density/1e6, states.N_pri*Av*states.density/1e6,
             states.N_agg*Av/1e3, states.N_pri*Av/1e3];

soot_columns += ["d_p[nm]", "d_m[nm]", "d_g[nm]", "n_p"]
soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_g*1e9, states.n_p];

soot_columns += ["total_soot_mass[kg]", "volume_fraction", "SSA[m2/g]"]

soot_data += [states.total_soot_mass * states.gas_mass, states.volume_fraction, states.SSA];

soot_columns += ["total_carbon_mass[kg]", "carbon_yield"];
soot_data += [states.total_carbon_mass * states.gas_mass,
               states.soot_carbon_mass/states.total_carbon_mass[0]]

soot_columns += ["total_hydrogen_mass[kg]"];
soot_data += [states.total_hydrogen_mass * states.gas_mass]

## Total Carbon & Hydrogen Mass
total_carbon_mass = states.total_carbon_mass * states.gas_mass;
total_hydrogen_mass = states.total_hydrogen_mass * states.gas_mass;

soot_columns += ["total_carbon_mass[kg]", "total_hydrogen_mass[kg]"];
soot_data += [total_carbon_mass, total_hydrogen_mass];

flow_columns += ["t[s]", "P[Pa]" , "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "gas_mass[kg]", "gas_volume[m3]",
                 "reactor_volume[m3]"];
flow_data += [states.restime, states.P, states.T, states.density, states.mean_molecular_weight/1000,
         states.gas_mass, states.gas_volume, states.reactor_volume];    
### Combining flow, soot & species columns 
columns = flow_columns + soot_columns + species_columns;
data = np.vstack((
    np.array(flow_data + soot_data + species_data),
)).transpose()

### creating data frame
df = pd.DataFrame(data = data, columns = columns);
### making the directory
if not os.path.exists(results_dir):
    os.makedirs(results_dir);
### Writing the file on disk
filename = f"{results_dir}/T{T5_K}K.csv"
df.to_csv(filename);
print(f"{filename} was written to the disk!")