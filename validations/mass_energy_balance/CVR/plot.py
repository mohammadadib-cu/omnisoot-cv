from collections import defaultdict, OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from pathlib import Path
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import shutil
import string

def save_plot(name, postprocess_dir, save_pdf = True, white_back = "None", pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    plt.tight_layout(pad = pad);
    facecolor = white_back;
    plt.savefig(f"{postprocess_dir}/{name}.svg", facecolor=facecolor);
    if save_pdf:
        filename = f"{postprocess_dir}/{name}.pdf";
        plt.savefig(filename);

def plot():
    if shutil.which('latex'):
        use_tex = True;
    else:
        use_tex = False;

    results_dir = "results";
    postprocess_dir = "post_process";

    linestyles = ["solid", "dashed", "dashdot", "dotted"];
    colors = ["#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#fc9d03", "#FF00FF", "#000000",
            "#ff0080", "#ff8000", "#4e5a00", "#a4805c", "#cb909a"]; 
    linewidth = 1.5;
    markersize = 10;
    ticklength = 5;
    tickwidth = 0.75;
    figsize = (6.4, 5);
    fontsize = 24;
    legend_fontsize = 20;
    axis_fontsize = 28;
    letter_fontsize = 36;

    PAH_growth_models = ["ReactiveDimerization", "DimerCoalescence", "EBridgeModified", "IrreversibleDimerization"];
    particle_dynamics_models = ["Monodisperse", "Sectional"];

    labels = {
        "ReactiveDimerization" : "Reac.Dim.",
        "DimerCoalescence" : "Dim.Coal.",
        "EBridgeModified" : "EBri.Mod.",
        "IrreversibleDimerization" : "Irrev.Dim."
    };

    results_df = defaultdict(
    lambda : defaultdict(
        dict
    )
    );
    for PAH_growth_model in PAH_growth_models:
        for particle_dynamics_model in particle_dynamics_models:
            results_df[PAH_growth_model][particle_dynamics_model] = pd.read_csv(
                f"{results_dir}/{particle_dynamics_model}/{PAH_growth_model}/sim_results.csv"
            );

    grid = (2,2);
    letters = string.ascii_lowercase;

    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][0];

    color_offset = 1
    step = 1;
    n_models = 3;
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)

    colors = [cmap(color_offset + i*step) for i in range(n_models)];

    plt.rcParams.update({'text.usetex' : use_tex, 'font.family':'Latin Modern Roman', 'font.size': fontsize});
    fig, axes = plt.subplots(nrows=grid[0], ncols=grid[1], figsize=(figsize[0]*grid[0], figsize[1]*grid[1]-1));
    axes = axes.reshape(grid[0]*grid[1]);
    for PAH_growth_model_index, PAH_growth_model in enumerate(PAH_growth_models):
        ax = axes[PAH_growth_model_index];
        for particle_dynamics_model_index, case in enumerate(particle_dynamics_models):
            df = results_df[PAH_growth_model][particle_dynamics_model]
            # Carbon Flux
            t_data = df["t[s]"].to_numpy()*1000;        
            C_flux = df["total_carbon_mass[kg]"].to_numpy();
            z_data = np.abs(C_flux - C_flux[0])/C_flux[0];
            label = " " if particle_dynamics_model_index == 0 else "Carbon";
            
            ax.plot(t_data, z_data,
                color=colors[0], linewidth = linewidth, label = label, linestyle = linestyles[particle_dynamics_model_index]);
            
            # Hydrogen Flux
            label = " " if particle_dynamics_model_index == 0 else "Hydrogen";
            H_flux = df["total_hydrogen_mass[kg]"];
            z_data = np.abs(H_flux - H_flux[0])/H_flux[0];
            ax.plot(t_data, z_data,
                color=colors[1], linewidth = linewidth, label = label, linestyle = linestyles[particle_dynamics_model_index]);
            
            # Energy Flux
            label = " " if particle_dynamics_model_index == 0 else "Energy";
            e_res = (
                (df["U_gas[J]"] - df["U_gas[J]"][0])
                + (df["U_soot[J]"] - df["U_soot[J]"][0])
            );
            z_data = np.abs(e_res)/ df["U_gas[J]"][0]
            
            ax.plot(t_data, z_data,
                color=colors[2], linewidth = linewidth, label = label, linestyle = linestyles[particle_dynamics_model_index]);
            
        # Y axis
        ax.set_yscale('log');
        if PAH_growth_model_index % 2 == 0:
            ax.set_ylabel("Relative error", {'fontsize' : 25})
        ax.set_ylim([1e-15, 3e-11]);

        # X axis
        if PAH_growth_model_index > 1:
            ax.set_xlabel("t [ms]", {'fontsize' : 25})
        ax.set_xlim([0, 10]);
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_minor_locator(MultipleLocator(0.4))
        
        if PAH_growth_model_index == 0:
            x_offset = -0.35;
            y_offset = 0.2;
            legend = ax.legend(loc = (0.4 + x_offset, 0.08 + y_offset), ncol = 2, fontsize = legend_fontsize,
                            labelspacing = 0.3, framealpha=0, columnspacing = 0.3);
            ax.text(0.41 + x_offset, 0.38 + y_offset, 'MPBM', transform=ax.transAxes, fontsize = 17);
            ax.text(0.6 + x_offset, 0.38 + y_offset, 'SPBM', transform=ax.transAxes, fontsize = 17);
            

        
        ax.text(0.05, 0.88, labels[PAH_growth_model], transform=ax.transAxes, fontsize = 22);
        
    save_plot("relerr_constuv", postprocess_dir = postprocess_dir)


if __name__ == "__main__":
    plot();