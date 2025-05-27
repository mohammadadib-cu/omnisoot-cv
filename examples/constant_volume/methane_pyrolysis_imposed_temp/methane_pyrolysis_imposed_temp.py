import cantera as ct
from omnisoot import ConstantVolumeReactor, SootGas, SootThermo
import numpy as np
import os
import pandas as pd
import time
import matplotlib.pyplot as plt
### Constants
Av = 6.0221408e+23;
### gas
gas = ct.Solution("Caltech_basic.yaml");
soot_gas = SootGas(gas);
### ---------------------------------------------------------------------------------
### Temperature Time History
time_data = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1];
Temp_data = [473, 1405, 1613, 1659, 1670, 1672, 1672, 1671, 1670, 1669, 1666, 1662, 1655, 1643, 1623, 1590, 1536, 1448, 1301, 1060, 663];
fixed_profile = np.array([time_data, Temp_data])

### Intial conditions of the reactor
T = Temp_data[0]; # K
P = 10e5; # Pa
X = {'CH4':0.8, 'CO2': 0.1, 'H2': 0.1};
### ---------------------------------------------------------------------------------
### Reactor
constUV = ConstantVolumeReactor(soot_gas);
### Initial conditions
soot_gas.TPX = T, P, X;
### Simulation time
max_res_time = 0.1;
### PAH growht models
soot = constUV.soot;
PAH_growth_model_type = ["ReactiveDimerization", "DimerCoalescence",
                        "CrossLinkingModified"][0]
soot.PAH_growth_model_type = PAH_growth_model_type;
### Setting precursors
soot.set_precursor_names(['A2', 'A2R5', 'A3R5', 'A4R5']);
### Reactor integrator settings
constUV.max_step = 1e-3;
constUV.max_time = max_res_time;
constUV.temperature_solver_type = "profile_time";
constUV.set_fix_temperature_profile(fixed_profile)
constUV.start();
### Soltuion array
states = ct.SolutionArray(gas, 1, extra={'restime': [0.0] ,
                                         'N_agg': [soot.N_agg], 'N_pri': [soot.N_pri], 
                                         'C_tot': [soot.C_tot], 'H_tot': [soot.H_tot],
                                         'A_tot': [soot.A_tot],
                                         'd_p': [soot.d_p], 'd_m': [soot.d_m],
                                         'd_v': [soot.d_v], 'd_g': [soot.d_g],
                                         'n_p': [soot.n_p],
                                         'total_mass' : [soot.total_mass],
                                         'volume_fraction':[soot.volume_fraction],
                                         'SSA': [soot.SSA],
                                         'inception_mass': [soot.inception_mass],
                                         'coagulation_mass': [soot.coagulation_mass],
                                         'PAH_adsorption_mass' : [soot.PAH_adsorption_mass],
                                         'surface_growth_mass' : [soot.surface_growth_mass],
                                         'oxidation_mass': [soot.oxidation_mass],
                                         'total_carbon_mass' : [constUV.total_carbon_mass],
                                         'total_hydrogen_mass' : [constUV.total_hydrogen_mass],
                                         'soot_carbon_mass' : [constUV.soot_carbon_mass],
                                            });
def append():
    states.append(gas.state ,restime=constUV.restime, 
              N_agg = soot.N_agg, N_pri =soot.N_pri,
              C_tot = soot.C_tot, H_tot = soot.H_tot, A_tot = soot.A_tot,
              d_p = soot.d_p, d_m = soot.d_m, d_v = soot.d_v, d_g = soot.d_g,
              n_p = soot.n_p,
              total_mass = soot.total_mass, volume_fraction = soot.volume_fraction, SSA = soot.SSA,
              inception_mass = soot.inception_mass,
              coagulation_mass = soot.coagulation_mass,
              PAH_adsorption_mass = soot.PAH_adsorption_mass,
              surface_growth_mass = soot.surface_growth_mass,
              oxidation_mass = soot.oxidation_mass,
              total_carbon_mass = constUV.total_carbon_mass,
              total_hydrogen_mass = constUV.total_hydrogen_mass,
              soot_carbon_mass = constUV.soot_carbon_mass                
             );
    
### Running the reactor
start_time = time.time()
step = 0;
while constUV.restime < max_res_time:
    step += 1;
    constUV.step();
    if step%2==0:
        append();
else:
    append();
end_time = time.time();
print(f"the simulation took {(end_time-start_time):0.4f} s");

### Building the data frame of the results
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

soot_columns += ["d_p[nm]", "d_m[nm]", "d_v[nm]", "d_g[nm]", "n_p"]
soot_data += [states.d_p*1e9, states.d_m*1e9, states.d_v*1e9, states.d_g*1e9, states.n_p];

soot_columns += ["soot_mass[kg/kg]", "volume_fraction", "SSA",
                 "total_carbon_mass[kg/kg]", "carbon_yield", "total_carbon_mass[kg/m3]",
                "total_hydrogen_mass[kg/kg]", "total_hydrogen_mass[kg/m3]"];
soot_data += [states.total_mass, states.volume_fraction, states.SSA,
              states.total_carbon_mass, 
              (states.soot_carbon_mass*states.density)/(states.total_carbon_mass[0]*states.density[0]),
             states.total_carbon_mass*states.density,
             states.total_hydrogen_mass, 
             states.total_hydrogen_mass*states.density];

soot_columns += ['inception_mass[mol/kg-s]', "coagulation_mass[mol/kg-s]",
                 "PAH_adsorption_mass[mol/kg-s]", "surface_growth_mass[mol/kg-s]"];
soot_data += [states.inception_mass, states.coagulation_mass, states.PAH_adsorption_mass, states.surface_growth_mass];

soot_columns += ['u_soot[J/m3]'];
u_mass_soot = np.array([SootThermo.u_mass_soot(T) for T in states.T]);
U_soot = states.total_mass*states.density*u_mass_soot
soot_data += [U_soot];

flow_columns += ["t[s]", "P[kPa]" , "T[K]", "density[kg/m3]", "mean_MW[kg/mol]", "u_gas[J/m3]"];
U_gas = states.density*states.int_energy_mass;
flow_data += [states.restime, states.P ,states.T, states.density, states.mean_molecular_weight/1000,
             U_gas];

species_columns = [f"{sp}" for sp in gas.species_names];
species_data += [states.X[:,i] for i in range(len(gas.species_names))];

columns = flow_columns + soot_columns + species_columns;
data = (np.array(flow_data + soot_data + species_data)).transpose();

### Writing files to disk
output_dir = f"results/{PAH_growth_model_type}";
if not os.path.exists(output_dir):
    os.makedirs(output_dir);
prop_df = pd.DataFrame(data=data, columns = columns)
prop_df.to_csv(f"{output_dir}/results.csv");


### Elemental Balance
plt.rcParams.update({'font.family':'Times New Roman', 'font.size': 22});
fig, ax = plt.subplots(figsize=(9,8));


total_C_vol = states.total_carbon_mass*states.density
C_error = np.abs(total_C_vol - total_C_vol[0])/total_C_vol;
ax.plot(states.restime, C_error,
    color="#FF0000", linewidth = 2, label = "C", linestyle = "solid");

total_H_vol = states.total_hydrogen_mass*states.density
H_error = np.abs(total_H_vol - total_H_vol[0])/total_H_vol;
ax.plot(states.restime, H_error,
    color="#00FF00", linewidth = 2, label = "H", linestyle = "dashed");

# Y axis
ax.set_yscale('log');
ax.set_ylabel("Relative Error", {'fontsize' : 25})

# X axis
ax.set_xlabel("t [s]", {'fontsize' : 25})
ax.set_xlim([0, max_res_time]);

legend = ax.legend(loc = 'best', fontsize = 22, framealpha=0);

plot_dir = f"elemetal_balance";
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir);
plt.tight_layout()
plt.savefig(f"{plot_dir}/elemental_plot.jpg")