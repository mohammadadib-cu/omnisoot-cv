from collections import defaultdict
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator);

# Function to save plots to disk
def _save_plot(name, save_pdf = True):
    output_dir = f"{POSTPROCESS_DIR}";
    if not os.path.exists(output_dir):
        os.makedirs(output_dir);
    plt.tight_layout(pad=0.3);
    plt.savefig(f"{output_dir}/{name}.svg");
    if save_pdf:
        plt.savefig(f"{output_dir}/{name}.pdf")
        

# Plot style constant
RESULTS_DIR = "results";
DATA_DIR = "data";
POSTPROCESS_DIR = "post_process";

LINESTYLES = ["solid", "dashed", "dashdot", "dotted"];
COLORS = ["#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#fc9d03", "#FF00FF", "#000000",
         "#ff0080", "#ff8000", "#4e5a00", "#a4805c", "#cb909a"];

LINEWIDTH = 1.5;
MARKERSIZE = 10;
TICKLENGTH = 5;
TICKWIDTH = 0.75;
FIGSIZE = (6.5, 5.7);
FONTSIZE = 22;
LEGEND_FONTSIZE = 22;
AXIS_FONTSIZE = 24;
LETTER_FONTSIZE = 36;

LABELS = {
    "Monodisperse" : "Monodisperse",
    "Sectional" : "Sectional",
    "Data" : "Araki et al. (2019)"
}

def plot(
    particle_dynamics_model_types,
    PAH_growth_model_type
):

    results_df = defaultdict(
        dict
    );

    # Reading the results of simulations for a given PAH growth model with
    # Monodisperse and Sectional population balance models
    for particle_dynamics_model_type in particle_dynamics_model_types:
        filepath = f"{RESULTS_DIR}/{particle_dynamics_model_type}/{PAH_growth_model_type}/sim_results.csv";
        results_df[particle_dynamics_model_type] = pd.read_csv(filepath);
        
    # Experimental Data
    exp_df = pd.read_csv(f"{DATA_DIR}/experimental_data.csv");
    
    param_dict = {
        "soot_mass"  : "soot_mass[ug/g]",
        "N_agg"      : "N_agg[#/g]",
        "d_p"        : "d_p[nm]",
        "d_m"        : "d_m[nm]"
    };
    for param in param_dict.keys():
        plt.rcParams.update({'font.family':'Times New Roman', 'font.size': FONTSIZE});
        fig, ax = plt.subplots(figsize=FIGSIZE);

        # Numerical
        for particle_dynamics_model_index, particle_dynamics_model in enumerate(particle_dynamics_model_types):
            z_data = results_df[particle_dynamics_model]["z[m]"]*100;
            f_data = results_df[particle_dynamics_model][param_dict[param]];
            ax.plot(
                z_data, f_data,
                color= COLORS[particle_dynamics_model_index], 
                linewidth = LINEWIDTH, label = LABELS[particle_dynamics_model],
                linestyle = "solid"
            );
            
        # Experimental
        z_data = exp_df["z[m]"]*100;
        f_data = exp_df[param_dict[param]];
        ax.plot(
            z_data, f_data,
            color= "#000000",
            markersize = MARKERSIZE,
            marker = "s",
            markerfacecolor = "#FFFFFF",
            label = LABELS["Data"],
            linestyle = "None"
        );
        
        # Legend
        legend = ax.legend(loc = 'best', fontsize = LEGEND_FONTSIZE, framealpha=0);
        
        # X axis
        ax.set_xlim([0, 65]);
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(2))
        ax.set_xlabel("z [cm]", {'fontsize' : AXIS_FONTSIZE})
        
        # Y Axis
        if param == "soot_mass":
            ax.set_ylabel("Soot mass $\mathregular{[\mu g/g]}$", {'fontsize' : AXIS_FONTSIZE})
            ax.set_ylim([0, 210]);
        elif param == "N_agg":
            ax.set_yscale('log');
            ax.set_ylabel("$\mathregular{N_{{agg}} [1/g]}$", {'fontsize' : AXIS_FONTSIZE})
            ax.set_ylim([1e11, 3e14]);
        elif param == "d_p":
            ax.set_ylabel("Primary Particle Diameter $\mathregular{[nm]}$", {'fontsize' : AXIS_FONTSIZE})
            ax.set_ylim([1, 15]);
        elif param == "d_m":
            ax.set_ylabel("Mobility Diameter $\mathregular{[nm]}$", {'fontsize' : AXIS_FONTSIZE})
            ax.set_ylim([2, 80]);
        _save_plot(param)