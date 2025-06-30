from collections import defaultdict
from pathlib import Path
import os
import string
import shutil

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, LogLocator)
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline


# Constants
Av = 6.0221408e+23;
MW_carbon = 0.012;
soot_density = 1800;

# Helper Functions
def save_plot(name, postprocess_dir, pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    # plt.tight_layout(pad = pad);
    filename = f"{postprocess_dir}/{name}.pdf";
    plt.savefig(filename);


def get_PSD(N_agg, d_m):
    dlog_dm = np.diff(np.log(d_m));
    N_agg_avg = (N_agg[1:] + N_agg[:-1])/2.0;
    dN_agg_dlogd_m = N_agg_avg/dlog_dm;
    d_m_avg = (d_m[1:] + d_m[:-1])/2.0;
    return d_m_avg, dN_agg_dlogd_m;

def get_d_g_sigma_d(N_agg, d_m):
    log_d_m = np.log(d_m);
    log_d_g = np.sum(N_agg * log_d_m) / np.sum(N_agg);
    d_g = np.exp(log_d_g)
    log_sigma_g = (
        np.sum(
            N_agg*
            (log_d_g - log_d_m) ** 2.0
        ) / np.sum(N_agg)
    )**0.5;
    sigma_g = np.exp(log_sigma_g)
    return d_g, sigma_g

def use_tex():
    """Check if LaTeX is installed on the system."""
    if shutil.which('latex'):
        return True;
    else:
        return False;

# Main Plotting Function
def main():    
    mech_names = ["FFCM2", "Caltech", "KAUST", "ABF1bar", "ITV", "CRECK"];
    # Plot configs
    result_dir = f"results/";
    postprocess_dir = f"post_process";
    
    linestyles = ["solid", "dashed", "dashdot", "dotted"];

    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][0];
    color_offset = 1
    step = 1;
    n_models = 4;
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)
    colors = [cmap(color_offset + i*step) for i in range(n_models)];
    
    title_fontsize = 26;
    legend_fontsize = 20;
    linewidth = 2;
    figsize = (6.5, 6.2);
    fontsize = 24;
    axis_fontsize = 26;
    plt.rcParams.update({"text.usetex": use_tex(), 'font.size': fontsize, 'figure.constrained_layout.use': True});
    
    labels = {
        "DimerCoalescence" : "Dim.Coal.",
        "ReactiveDimerization":"Reac.Dim.",
        "EBridgeModified" : "EBri.Mod.",
        "IrreversibleDimerization" : "Irrev.Dim.",
    }

    eq_ratios = [1.9, 2.0, 2.1];
    PAH_growth_model_types = [ "ReactiveDimerization", "DimerCoalescence", "EBridgeModified", "IrreversibleDimerization"];

    # Reading the simulation results
    results_df = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                dict
            )
        )
    );

    for eq_ratio in eq_ratios:
        for PAH_growth_model_type in PAH_growth_model_types:
            filename = f"{result_dir}/{eq_ratio}/KAUST/Sectional/{PAH_growth_model_type}/sim_results_pfr.csv";
            results_df[eq_ratio][PAH_growth_model_type] = pd.read_csv(filename);

    # Reading the experimental data
    exp_df = {};
    exp_df[1.9] = pd.read_csv("data/phi_1.9.csv");
    exp_df[2.0] = pd.read_csv("data/phi_2.0.csv");
    exp_df[2.1] = pd.read_csv("data/phi_2.1.csv");

    # ---------------------------------------------------------------------------------------------------------------
    # Ploting Particle Size Distribution (PSD)
    # ---------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 3, nrows = 1, figsize=(figsize[0] * 3, figsize[1]));

    for eq_index, eq_ratio in enumerate(eq_ratios):
        ax = axes[eq_index];
        title = f"({string.ascii_lowercase[eq_index]}) "+"$\\mathrm{\\phi}=$"+f"{eq_ratio}";
        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[eq_ratio][PAH_growth_model_type];

            density = df["density[kg/m3]"].iloc[-1];

            params = [f"d_m_{i}" for i in range(60)];
            sub_df = df[params];
            d_m = sub_df.iloc[-1, :].to_numpy()*1e9;
            params = [f"N_agg_{i}" for i in range(60)];
            sub_df = df[params];
            N_agg = sub_df.iloc[-1, :].to_numpy() * Av * density / 1e6;
            x_data, y_data = get_PSD(N_agg, d_m);
            

            d_g, sigma_g = get_d_g_sigma_d(N_agg, d_m);

            pos_index = np.where(np.diff(x_data)<0)[0][0];
            x_data = x_data[:(pos_index+1)]
            x_data_log = np.log(x_data);

            
            y_data = y_data * (y_data>0) + 10 * (y_data<0);
            y_data = y_data[:(pos_index+1)]
            y_data_log = np.log(y_data);
            
            x_data_log_ext = np.linspace(x_data_log[0], x_data_log[-1], 201);
            cs = CubicSpline(x_data_log, y_data_log);
            
            y_data_log_ext = cs(x_data_log_ext);
            
            x_data_ext = np.exp(x_data_log_ext);
            y_data_ext = np.exp(y_data_log_ext);
            
            ax.plot(x_data_ext, y_data_ext,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = labels[PAH_growth_model_type], linestyle = "solid");

        # Experimental
        x_data_exp = exp_df[eq_ratio]["dm"];
        y_data_exp = exp_df[eq_ratio]["dN/dm"];
        ax.plot(x_data_exp, y_data_exp, marker = "o", markersize = 10, markerfacecolor = "#FFFFFF",
            color="#000000", linewidth = linewidth, label = "Data", linestyle = "None"); 
        # Y axis
        if eq_index == 0:
            ax.set_ylabel("${dN_{agg}/d\\mathrm{log}(d_m)}$", {'fontsize' : axis_fontsize})

        ax.set_yscale('log');
        ax.set_ylim([1e5, 2e10]);

        # X axis
        ax.set_xlabel("$d_m$ nm", {'fontsize' : axis_fontsize})

        ax.set_xscale('log');
        ax.set_xlim([1, 700]);
        
        legend_loc = [
            (0.5, 0.5),
            (0.55, 0.55),
            (0.23, 0.0)
        ];
        
        if eq_index == 0:  
            ax.legend(loc = legend_loc[eq_index], ncol = 1, fontsize = legend_fontsize,
                            framealpha=0, labelspacing=0.3, columnspacing = 0,
                            handletextpad=0.52);
        title_pos = [
            {"x": 0.05, "y": 0.90},
            {"x": 0.65, "y": 0.9},
            {"x": 0.65, "y": 0.90},
        ]
        ax.text(title_pos[eq_index]["x"], title_pos[eq_index]["y"], title, transform=ax.transAxes, fontsize = title_fontsize);

    save_plot("PSD_eq_ratio_log", postprocess_dir = postprocess_dir);


    # ---------------------------------------------------------------------------------------------------------------
    # Primary Particle Diameter (along PFR)
    # ---------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 3, nrows = 1, figsize=(figsize[0] * 3, figsize[1]));

    param = "d_p[nm]";

    for eq_index, eq_ratio in enumerate(eq_ratios):
        title = f"({string.ascii_lowercase[eq_index]}) "+"$\\mathrm{\\phi}=$"+f"{eq_ratio}";
        ax = axes[eq_index];
        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[eq_ratio][PAH_growth_model_type];

            z_data = df['z[m]'].to_numpy();
            f_data = df[param];


            ax.plot(z_data, f_data,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = labels[PAH_growth_model_type], linestyle = "solid"); 
        # Y axis
        if eq_index == 0:
            ax.set_ylabel("Primary particle diameter, $d_p$ [nm]", {'fontsize' : axis_fontsize})

        ax.set_ylim([2.2, 13]);
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.4));

        # X axis
        ax.set_xlabel("z [m]", {'fontsize' : axis_fontsize})

        ax.set_xlim([0, 0.7]);
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02));

        legend = ax.legend(loc = "lower right", ncol = 1, fontsize = legend_fontsize,
                        framealpha=0, labelspacing=0.3, columnspacing = 0,
                        handletextpad=0.5);
        
        ax.text(0.04, 0.9, title, transform=ax.transAxes, fontsize = title_fontsize);

    save_plot("d_p_eq_ratio_all_single_mech", postprocess_dir = postprocess_dir)

    # ---------------------------------------------------------------------------------------------------------------
    # Soot volume fraction (along the PFR)
    # ---------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 3, nrows = 1, figsize=(figsize[0] * 3, figsize[1]));

    param = "volume_fraction[-]"

    for eq_index, eq_ratio in enumerate(eq_ratios):
        title = f"({string.ascii_lowercase[eq_index]}) "+"$\\mathrm{\\phi}=$"+f"{eq_ratio}"
        ax = axes[eq_index];
        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[eq_ratio][PAH_growth_model_type];
            z_data = df['z[m]'].to_numpy();
            f_data = df[param];
            ax.plot(z_data, f_data,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = labels[PAH_growth_model_type], linestyle = "solid"); 
        # Y axis
        if eq_index == 0:
            ax.set_ylabel("Soot volume fraction, $f_v$", {'fontsize' : axis_fontsize})


        # X axis
        ax.set_xlabel("z [m]", {'fontsize' : axis_fontsize})

        ax.set_xlim([0, 0.7]);
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02));
        
        if eq_index == 0 or eq_index == 1:
            loc = (0.01, 0.5);
        elif eq_index == 2:
            loc = "center left";
        
        if eq_index == 2:
            ax.legend(loc = loc, ncol = 1, fontsize = legend_fontsize,
                        framealpha=0, labelspacing=0.3, columnspacing = 0,
                        handletextpad=0.5);
        
        ax.text(0.35, 0.9, title, transform=ax.transAxes, fontsize = title_fontsize);

    save_plot("f_v_eq_ratio_all_single_mech", postprocess_dir = postprocess_dir);

    # ---------------------------------------------------------------------------------------------------------------
    # A2R5 & Inception Flux
    # ---------------------------------------------------------------------------------------------------------------
    
    fig, axes = plt.subplots(ncols = 3, nrows = 2, figsize=(figsize[0] * 3, figsize[1]+2));
    param = "A2R5";

    for eq_index, eq_ratio in enumerate(eq_ratios):
        title = f"({string.ascii_lowercase[eq_index]}) "+"$\\mathrm{\\phi}=$"+f"{eq_ratio}";
        ax = axes[0, eq_index];

        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[eq_ratio][PAH_growth_model_type];
            z_data = df["z[m]"].to_numpy();
            f_data = df[param];

            ax.plot(z_data, f_data,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = labels[PAH_growth_model_type], linestyle = "solid"); 
        # Y axis
        ax.set_ylim([1e-7, 1e-4])
        ax.set_yscale("log");
        ax.yaxis.set_major_locator(LogLocator(numticks = 6))
        ax.yaxis.set_minor_locator(LogLocator(subs = "all", numticks = 6))
        
        if eq_index == 0:
            ax.set_ylabel("Mole fraction of A2R5", {'fontsize' : axis_fontsize})

        # X axis
        ax.set_xlim([0, 0.7]);
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02));
        
        if eq_index == 0:
            legend = ax.legend(loc = (0.5, 0.05), ncol = 1, fontsize = legend_fontsize,
                        framealpha=0, labelspacing=0.3, columnspacing = 0,
                        handletextpad=0.5);
        
        ax.text(0.45, 0.85, title, transform=ax.transAxes, fontsize = title_fontsize);


    param = "inception_mass[mol/kg-s]";
    for eq_index, eq_ratio in enumerate(eq_ratios):
        title = f"({string.ascii_lowercase[eq_index+3]}) "+"$\\mathrm{\\phi}=$"+f"{eq_ratio}";
        ax = axes[1, eq_index];

        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[eq_ratio][PAH_growth_model_type];
            z_data = df["z[m]"].to_numpy();
            f_data = df[param] * Av * df["density[kg/m3]"] / 1e6;

            ax.plot(z_data, f_data,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = labels[PAH_growth_model_type], linestyle = "solid"); 
    
    # Y axis
    if eq_index == 0:
        ax.set_ylabel("Soot inception flux\n$I_{inc}$ [1/$\\mathrm{cm^3}$-s]", {'fontsize' : axis_fontsize})

    ax.set_ylim([1e8, 1e13]);
    ax.set_yscale("log");
    ax.yaxis.set_major_locator(LogLocator(numticks = 6))
    ax.yaxis.set_minor_locator(LogLocator(subs = "all", numticks = 6))
    
    # X axis
    ax.set_xlim([0, 0.7]);
    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.02));
    ax.set_xlabel("z [m]", {'fontsize' : axis_fontsize})
    
    if eq_index == 0:
        ax.legend(loc = (0.5, 0.05), ncol = 1, fontsize = legend_fontsize,
                       framealpha=0, labelspacing=0.3, columnspacing = 0,
                      handletextpad=0.5);

if __name__ == "__main__":
    main();

