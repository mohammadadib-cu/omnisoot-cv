import cantera as ct
from omnisoot import PressureReactor, SootGas
import numpy as np
import pandas as pd
import os

Av = 6.0221408e+23; #1/mol;
MECH_DICT = {
    "Caltech":"./data/Caltech.yaml",
    "KAUST": "./data/KAUST.yaml",
    "ABF1bar" : "./data/abf_1bar.yaml",
    "CRECK" : "./data/CRECK_no_BINS.yaml",
    "ITV" : "./data/ITV_PAH.yaml",
    "NUIG" : "./data/NUIG.yaml",
    "FFCM2" : "./data/FFCM2.yaml",
}

def simulate_nosoot(
    mech_name, 
    T5_K, 
    P5_atm, 
    initial_composition, 
    max_res_time, 
    results_dir = "results"
):
    gas = ct.Solution(MECH_DICT[mech_name]);
    soot_gas = SootGas(gas);
    soot_gas.TPX = T5_K, P5_atm * ct.one_atm, initial_composition;
    ### Building and configruing the reactor
    reactor = PressureReactor(soot_gas);
    reactor.pressure_type = "constant";
    reactor.max_step = 1e-4;
    reactor.first_step = 1e-23;
    reactor.solver_type = "BDF";
    soot = reactor.soot;
    soot.soot_enabled = False;
    ### Starting the reactor
    reactor.start();

    def append():
        states.append(gas.state ,restime=reactor.restime
                 );
    
    states = ct.SolutionArray(gas, 1, extra={'restime': [0.0]
                                            });


    ### Running the simulation
    step = 0;
    while reactor.restime < max_res_time:
        reactor.step();
        step += 1;
        if step % 2 == 0:
            append();
    append();
   
    ### Building the output data
    flow_columns = [];
    flow_data = [];
    species_columns = [];
    species_data = [];

    flow_columns += ["t[s]", "P[Pa]" , "T[K]", "density[kg/m3]", "mean_MW[kg/mol]"];
    flow_data += [states.restime, states.P, states.T, states.density, states.mean_molecular_weight/1000];

    species_columns = [f"{sp}" for sp in gas.species_names];
    species_data += [states.X[:,i] for i in range(len(gas.species_names))];
    
    ### Combining flow, soot & species columns 
    columns = flow_columns + species_columns;
    data = np.vstack((
        np.array(flow_data + species_data),
    )).transpose()

    ### creating data frame
    df = pd.DataFrame(data = data, columns = columns);
    ### making the directory
    output_dir = f"{results_dir}/{mech_name}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    ### Writing the file on disk
    filename = f"{output_dir}/sim_results.csv"
    df.to_csv(filename);
    print(f"{filename} was written to the disk!")