from pathlib import Path
import os

import numpy as np
import pandas as pd
from scipy.integrate import simpson


def main():
    files_paths = Path(
         "results/"
    );
    for file_path in list(files_paths.glob("**/*.csv")):
        df = pd.read_csv(file_path);
        t_array = df["t[s]"].to_numpy();
        z_array = df["z[m]"].to_numpy();
        m_array = (1.0 / df["velocity[m/s]"]).to_numpy();
        C_tot_inception = df["C_tot_inception[mol/kg-s]"].to_numpy() * m_array;
        C_tot_PAH = df["C_tot_PAH[mol/kg-s]"].to_numpy() * m_array;
        C_tot_growth = df["C_tot_growth[mol/kg-s]"].to_numpy() * m_array;
        C_tot_tot = df["C_tot[mol/kg]"].to_numpy() * m_array;
        C_tot_rate = (C_tot_inception + C_tot_PAH + C_tot_growth);
        
        C_tot_int = np.zeros_like(C_tot_tot);
        C_tot_inception_int = np.zeros_like(C_tot_tot);
        C_tot_PAH_int = np.zeros_like(C_tot_tot);
        C_tot_growth_int = np.zeros_like(C_tot_tot);
        
        for i in range(1, len(t_array)):
            C_tot_int[i] = simpson(C_tot_rate[:(i+1)], x = z_array[:(i+1)]);
            C_tot_inception_int[i] = simpson(C_tot_inception[:(i+1)], x = z_array[:(i+1)]);
            C_tot_PAH_int[i] = simpson(C_tot_PAH[:(i+1)], x = z_array[:(i+1)]);
            C_tot_growth_int[i] = simpson(C_tot_growth[:(i+1)], x = z_array[:(i+1)]);
        
        df["C_tot_int[mol/kg]"] = C_tot_int;
        df["C_tot_inception_int[mol/kg]"] = C_tot_inception_int;
        df["C_tot_PAH_int[mol/kg]"] = C_tot_PAH_int;
        df["C_tot_growth_int[mol/kg]"] = C_tot_growth_int;
        
        df.to_csv(file_path, index = False)

if __name__ == "__main__":
    main();