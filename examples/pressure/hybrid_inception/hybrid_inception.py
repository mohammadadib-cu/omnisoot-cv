# ----------------------------------------------------------------------------------------------------------------------------------
# This example illustrates the use of the HybridDimerization model in Omnisoot (>0.3.0).
# HybridDimerization allows users to combine ReactiveDimerization, IrriversibleDimerization, DimerCoalescence, and EBridgeModified
# as submodels with independent precursors, and inception and adsorption prefactors.
# ----------------------------------------------------------------------------------------------------------------------------------

# Declarations
import cantera as ct
import omnisoot as oms
import numpy as np
import pandas as pd
import os

# Loading the reaction mechanism
gas = ct.Solution(oms.data_dir/"Caltech.yaml")
# Simulation parameters
T5_K = 2203;
P5_atm = 4.5;
initial_composition = {"CH4":0.5, "AR": 0.95};
max_res_time = 5e-3;
results_dir = "results";


# Setting up the reactor
soot_gas = oms.SootGas(gas);
soot_gas.TPX = T5_K, P5_atm * ct.one_atm, initial_composition;
### Building and configruing the reactor
reactor = oms.PressureReactor(soot_gas);
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
states = ct.SolutionArray(gas, 1, extra = reactor.solution_data_dict);
### Running the simulation
step = 0;
while reactor.restime < max_res_time:
    reactor.step();
    step += 1;
    if step%10==0:
        states.append(gas.state, **reactor.solution_data_dict);

### Combining flow, soot & species columns 
oms_data = [getattr(states, key) for key in reactor.solution_keys];
species_columns = [f"{sp}" for sp in gas.species_names];
species_data = [states.X[:,i] for i in range(len(gas.species_names))];
data = np.vstack((
    np.array(oms_data + species_data),
)).transpose();
columns = reactor.solution_headers + species_columns

### creating data frame
df = pd.DataFrame(data = data, columns = columns);
### making the directory
if not os.path.exists(results_dir):
    os.makedirs(results_dir);
### Writing the file on disk
filename = f"{results_dir}/sim_results.csv"
df.to_csv(filename);
print(f"{filename} was written to the disk!")