import os
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import pandas as pd


LINEWIDTH = 2;
MARKERSIZE = 10;
TICKLENGTH = 5;
TICKWIDTH = 0.75;
PLOT_FONTSIZE = 22;
LEGEND_FONTSIZE = 18;
AXIS_FONTSIZE = 26;
LETTER_FONTSIZE = 25;
TITLE_FONTSIZE = 24;
TIME_FACTOR = 1000;
FIGURE_PATH = "figures"

def save_plot(filename, pad = 1.5):
    if not os.path.exists(FIGURE_PATH):
        os.makedirs(FIGURE_PATH);
    plt.tight_layout(pad = pad);
    plt.savefig(f"{FIGURE_PATH}/{filename}.jpg");
    plt.savefig(f"{FIGURE_PATH}/{filename}.pdf");


def plot(
    mech_name, 
    PAH_growth_model_types, 
    particle_dynamics_model_type,
    results_dir
):
    color_offset = 1
    step = 1;
    n_models = 4;
    cmap = plt.get_cmap('hsv', len(PAH_growth_model_types)*step + color_offset)
    colors = [cmap(color_offset + i*step) for i in range(n_models)];

    ### Reading files from disk
    results_df = defaultdict(dict);
    for PAH_growth_model_type in PAH_growth_model_types:
        results_df[PAH_growth_model_type] = pd.read_csv(f"{results_dir}/{mech_name}/{particle_dynamics_model_type}/{PAH_growth_model_type}/sim_results.csv");

    ### ----------------------------------------------------------------------------------------------------------------------------------------
    ### Energy and mass residual
    plt.rcParams.update({'font.family':'Times New Roman', 'font.size': PLOT_FONTSIZE});
    fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize = (12, 10));   
    axes = axes.reshape(-1,);
    for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
        ax = axes[PAH_growth_model_type_index];
        df = results_df[PAH_growth_model_type];
        
        t_data = df["t[s]"];
        # Residual of elemental carbon mass
        val = df["total_carbon_mass[kg]"].to_numpy();
        f_data = np.abs(val - val[0])/val[0]
        ax.plot(t_data * TIME_FACTOR, f_data, linewidth = LINEWIDTH, color = colors[0], label = "Carbon");

        # Residual of elemental hydrogen mass
        val = df["total_hydrogen_mass[kg]"].to_numpy();
        f_data = np.abs(val - val[0])/val[0]
        ax.plot(t_data * TIME_FACTOR, f_data, linewidth = LINEWIDTH, color = colors[1], label = "Hydrogen");

        # Residual of elemental energy
        val = (df["H_gas[J]"] + df["H_soot[J]"]).to_numpy();
        f_data = np.abs(val - val[0])/val[0]
        ax.plot(t_data * TIME_FACTOR, f_data, linewidth = LINEWIDTH, color = colors[2], label = "Energy");

        if PAH_growth_model_type_index % 2 == 0:
            ax.set_ylabel("Residual");
        ax.set_ylim([0, 3e-11]);
        
        
        if PAH_growth_model_type_index > 1:
            ax.set_xlabel("t [ms]");
        ax.xaxis.set_major_locator(MultipleLocator(20))
        ax.xaxis.set_minor_locator(MultipleLocator(4))
        ax.set_xlim([-0.1, 100.1]);

        legend = ax.legend(loc = "best", fontsize = LEGEND_FONTSIZE, framealpha=0, handletextpad = 0.3, handlelength = 1.5);

    save_plot("mass_energy_residual");

    ### ----------------------------------------------------------------------------------------------------------------------------------------
    ### Number of agglomerates
    plt.rcParams.update({'font.family':'Times New Roman', 'font.size': PLOT_FONTSIZE});
    fig, ax = plt.subplots(ncols = 1, nrows = 1, figsize = (6, 5)); 
    for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
        df = results_df[PAH_growth_model_type];
        t_data = df["t[s]"];
        f_data = df["N_agg[#/m3]"].to_numpy();
        ax.plot(t_data * TIME_FACTOR, f_data, linewidth = LINEWIDTH, color = colors[PAH_growth_model_type_index], label = PAH_growth_model_type);

    legend = ax.legend(loc = "best", fontsize = LEGEND_FONTSIZE, framealpha=0, handletextpad = 0.3, handlelength = 1.5);

    ax.set_ylabel("$N_{agg}$ [1/$m^3$]");
    ax.set_yscale("log")
    ax.set_ylim([1e16, 1e19]);

    ax.set_xlabel("t [ms]");
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    ax.set_xlim([-1, 50.1]);

    save_plot("N_agg");

    ### ----------------------------------------------------------------------------------------------------------------------------------------
    ### Primary Particle Diameter
    plt.rcParams.update({'font.family':'Times New Roman', 'font.size': PLOT_FONTSIZE});
    fig, ax = plt.subplots(ncols = 1, nrows = 1, figsize = (6, 5)); 
    for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
        df = results_df[PAH_growth_model_type];
        t_data = df["t[s]"];
        f_data = df["d_p[nm]"].to_numpy();
        ax.plot(t_data * TIME_FACTOR, f_data, linewidth = LINEWIDTH, color = colors[PAH_growth_model_type_index], label = PAH_growth_model_type);

    legend = ax.legend(loc = "best", fontsize = LEGEND_FONTSIZE, framealpha=0, handletextpad = 0.3, handlelength = 1.5);

    ax.set_ylabel("$d_p$ [nm]");
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.set_ylim([2, 20]);

    ax.set_xlabel("t [ms]");
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    ax.set_xlim([-1, 50.5]);

    save_plot("d_p");

    ### ----------------------------------------------------------------------------------------------------------------------------------------
    ### Specific surface area
    plt.rcParams.update({'font.family':'Times New Roman', 'font.size': PLOT_FONTSIZE});
    fig, ax = plt.subplots(ncols = 1, nrows = 1, figsize = (6, 5)); 
    for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
        df = results_df[PAH_growth_model_type];
        t_data = df["t[s]"];
        f_data = df["SSA[m2/g]"].to_numpy();
        ax.plot(t_data * TIME_FACTOR, f_data, linewidth = LINEWIDTH, color = colors[PAH_growth_model_type_index], label = PAH_growth_model_type);

    legend = ax.legend(loc = "best", fontsize = LEGEND_FONTSIZE, framealpha=0, handletextpad = 0.3, handlelength = 1.5);

    ax.set_ylabel("SSA [$m^2$/g]");
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(20))
    ax.set_ylim([100, 800]);

    ax.set_xlabel("t [ms]");
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    ax.set_xlim([-1, 50.5]);

    save_plot("SSA");