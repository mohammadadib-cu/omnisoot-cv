from collections import defaultdict
from pathlib import Path
import os
import shutil
import string

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, LogLocator)
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

def save_plot(name, postprocess_dir, pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    # plt.tight_layout(pad = pad);
    filename = f"{postprocess_dir}/{name}.pdf";
    plt.savefig(filename);

def use_tex():
    """Check if LaTeX is installed on the system."""
    if shutil.which('latex'):
        return True;
    else:
        return False;


def custom_skewed_function(x, mu=0, sigma=1, alpha=0):
    """
    Custom function with adjustable asymmetry.
    """
    left = np.exp(-((x - mu) / (sigma * (1 + alpha)))**2)
    right = np.exp(-((x - mu) / (sigma * (1 - alpha)))**2)
    return np.where(x <= mu, left, right)

def plot():

    # Plot configurations
    postprocess_dir = "post_process";
    linestyles = ["solid", "dashed", "dashdot", "dotted"];
    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][1];

    color_offset = 1
    step = 1;
    n_models = 4;
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)

    colors = [cmap(color_offset + i*step) for i in range(n_models)];

    linewidth = 2;
    figsize = (6.1, 5.3);
    fontsize = 24;
    axis_fontsize = 26;
    letter_fontsize = 26;
    legend_fontsize = 20;

    plt.rcParams.update({"text.usetex": use_tex(), 'font.size': fontsize, 'figure.constrained_layout.use': True});

    labels = {
        "DimerCoalescence" : "Dim.Coal.",
        "ReactiveDimerization":"Reac.Dim.",
        "EBridgeModified" : "EBri.Mod.",
        "IrreversibleDimerization" : "Irrev.Dim.",
        "data" : "Data",
    };


    column_dict = {
        "mean_MW" : "mean_MW[kg/mol]",
        # Soot
        "N_agg" : "N_agg[#/cm3]",
        "n_p" : "n_p",
        "d_p" : "d_p[nm]",
        "d_m" : "d_m[nm]",
        "carbon_yield" : "carbon_yield",
        "N_pri" : "N_pri[#/cm3]",
        "C_tot_int": "C_tot_int[mol]",
        "C_tot_inception_int": "C_tot_inception_int[mol]",
        "C_tot_PAH_int": "C_tot_PAH_int[mol]",
        "C_tot_growth_int": "C_tot_growth_int[mol]",
    };


    headers = ["carbon_yield", "d_p", "mean_MW",
             "N_agg", "N_pri", "n_p",
             "C_tot_int", "C_tot_inception_int", "C_tot_PAH_int", "C_tot_growth_int"
            ];
    
    T_interp = np.linspace(1800, 2900, 201);
    
    PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];

    results_df = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(
                    dict
                )
            )
        )
    );

    # Reading numerical results from CSV files
    for PAH_growth_model_type in PAH_growth_model_types:
        files_paths = Path(
            f"results/Caltech/Sectional/{PAH_growth_model_type}/"
        );
        for file_path in files_paths.glob("**/*.csv"):
            T5 = file_path.name[1:-5];
            df = pd.read_csv(file_path);
            t_array = df["t[s]"];
            for header in headers:
                param = column_dict[header];
                f_array = df[param];
                results_df[PAH_growth_model_type][header][float(T5)] = np.interp(1.5e-3, t_array, f_array);

    # Reading experimental results from CSV files
    Em_s = [0.174, 0.37]
    df_exp =  defaultdict(
        lambda: defaultdict(
            dict
        )
    );

    df = pd.read_csv(f"data/yield_Em_5CH4.csv")
    df_exp["carbon_yield"]["carbon_yield_max"] = df["carbon_yield_Em"].to_numpy()/Em_s[0];
    df_exp["carbon_yield"]["carbon_yield_min"] = df["carbon_yield_Em"].to_numpy()/Em_s[1];
    df_exp["carbon_yield"]["carbon_yield_mean"] = (
        df_exp["carbon_yield"]["carbon_yield_max"] 
        + df_exp["carbon_yield"]["carbon_yield_min"]
    ) / 2.0;
    df_exp["carbon_yield"]["T"] = df["T"];


    # ----------------------------------------------------------------------------------------------------------------------------------
    # Plotting Yield & Primary particle diameter
    # ----------------------------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 2, figsize=(figsize[0]*2+0.5, figsize[1]));

    headers = ["carbon_yield", "d_p"];
    scale_factors = [100, 1];

    limiters = [
        (0, 1),
        (2, 1000)
    ]

    for header_index, header in enumerate(headers):
        title = f"({string.ascii_lowercase[header_index]})";
        ax = axes[header_index];
        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            val_dict = results_df[PAH_growth_model_type][header];
            x_data = np.array(list(val_dict.keys()));
            f_data = np.array(list(val_dict.values()));
            cs = CubicSpline(x_data, f_data);
            f_data = np.clip(cs(T_interp), *limiters[header_index]);
            x_data = T_interp
            ax.plot(x_data, f_data * scale_factors[header_index],
                linewidth = linewidth, linestyle = linestyles[0],
                    label = labels[PAH_growth_model_type], color = colors[PAH_growth_model_type_index]);


        if header == "carbon_yield":
            ### Experimental
            exp_color = "#002038"
            # Trendline
            mu = 2210;
            sigma = 190;
            alpha = -0.1;
            factor = 12   

            f_data_exp_ext = custom_skewed_function(T_interp, mu = mu, sigma = sigma, alpha = alpha)*factor;
            ax.plot(T_interp, f_data_exp_ext, linestyle = "dotted",  color=exp_color, linewidth = 1.5,);

            x_data_exp = df_exp[header]["T"];
            f_data_exp_mean = df_exp[header][f"{header}_mean"];
            f_data_exp_max = df_exp[header][f"{header}_max"];
            f_data_exp_min = df_exp[header][f"{header}_min"];

            ax.errorbar(x_data_exp, f_data_exp_mean, 
                    (f_data_exp_mean - f_data_exp_min,  f_data_exp_max - f_data_exp_mean), linestyle = "none",
                    marker = None , elinewidth=0.5,
                    color=exp_color, capsize=5, alpha = 0.8)

            ax.plot(x_data_exp, f_data_exp_mean, linestyle = "none",
                    marker="o", markersize = 9, markerfacecolor="#ffff",  color=exp_color,
                    label = labels["data"], alpha = 0.9);
        


        # Y axis
        if header_index == 0:
            ax.set_ylabel("Carbon yield, CY [\\%]", {'fontsize' : axis_fontsize});
            ax.set_ylim([-0.2, 18.1]);
            ax.yaxis.set_major_locator(MultipleLocator(4))
            ax.yaxis.set_minor_locator(MultipleLocator(0.8));
        else:
            ax.set_ylabel("Primary particle diameter\n$d_p$ [nm]", {'fontsize' : axis_fontsize});
            ax.set_ylim([1.5, 14.5]);
            ax.yaxis.set_major_locator(MultipleLocator(2))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            
        # X axis
        ax.set_xlabel("$\\mathrm{T_5}$ [K]", {'fontsize' : axis_fontsize})
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(40))
        ax.set_xlim([1800, 2900]);     
        
        text_pos = [
            (0.05, 0.88),
            (0.85, 0.88),
        ]

        ax.text(*text_pos[header_index], title, transform=ax.transAxes, fontsize = letter_fontsize)

        legend_pos = [
            (0.58, 0.45),
            (0.05, 0.55),
        ]
        
        ax.legend(loc = legend_pos[header_index], ncol = 1, fontsize = legend_fontsize, columnspacing = 0.1,
                    framealpha=0, handletextpad = 0.3, handlelength = 1.5);

    save_plot("carbon_yield_d_p_5CH4", postprocess_dir = postprocess_dir, pad = 0.4);


    # ----------------------------------------------------------------------------------------------------------------------------------
    # N_agg, N_pri & n_p
    # ----------------------------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 2, figsize=(figsize[0]*2+0.5, figsize[1]));

    headers = [["N_agg", "N_pri"], ["n_p"]];
    scale_factors = [1, 1];

    limiters = [
        (0, 1),
        (2, 1000)
    ]

    for header_index, header in enumerate(headers):
        title = f"({string.ascii_lowercase[header_index]})";
        ax = axes[header_index];
        for sub_header_index, sub_header in enumerate(header):
            for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
                val_dict = results_df[PAH_growth_model_type][sub_header];
                x_data = np.array(list(val_dict.keys()));
                f_data = np.array(list(val_dict.values()));
                if header_index == 0:
                    cs = CubicSpline(x_data, np.log10(f_data));
                    f_data = 10**cs(T_interp);
                else:
                    cs = CubicSpline(x_data, f_data);
                    f_data = np.clip(cs(T_interp), *limiters[header_index]);
                x_data = T_interp
                if sub_header_index == 0 and len(header) > 1:
                    label = " ";
                else:
                    label = labels[PAH_growth_model_type];
                ax.plot(x_data, f_data * scale_factors[header_index],
                    linewidth = linewidth, linestyle = linestyles[sub_header_index],
                        label = label, color = colors[PAH_growth_model_type_index]);
            
        # Y axis
        if header_index == 0:
            ax.set_ylabel("Number of particles [1/$\\mathrm{cm^3}$]", {'fontsize' : axis_fontsize});
            ax.set_ylim([1e9, 1.2e14]);
            ax.set_yscale("log");
            ax.yaxis.set_major_locator(LogLocator(base = 10, subs = (1.0,), numticks = 10))
        else:
            ax.set_ylabel("Number of primary particles\n per agglomerate, $n_p$", {'fontsize' : axis_fontsize-1});
            ax.set_ylim([0, 22]);
            ax.yaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
            
        # X axis
        ax.set_xlabel("$\\mathrm{T_5}$ [K]", {'fontsize' : axis_fontsize-1})
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(40))
        ax.set_xlim([1800, 2900]);     
        
        text_pos = [
            (0.85, 0.85),
            (0.85, 0.85),
        ]

        ax.text(*text_pos[header_index], title, transform=ax.transAxes, fontsize = letter_fontsize)

        legend_pos = [
            (0.03, 0.0),
            (0.55, 0.15),
        ]

        legend_col = [
            2,
            1,
        ];

        if header_index == 0:
            ax.text(legend_pos[header_index][0]+0.025, legend_pos[header_index][0]+0.39, "$N_{agg}$", transform=ax.transAxes, fontsize = 18);
            ax.text(legend_pos[header_index][0]+0.18, legend_pos[header_index][0]+0.39, "$N_{pri}$", transform=ax.transAxes, fontsize = 18);
        
            ax.legend(loc = legend_pos[header_index], ncol = legend_col[header_index], fontsize = legend_fontsize, columnspacing = 0.7,
                    framealpha=0, handletextpad = 0.5, handlelength = 1.5);
    
    save_plot("N_agg_N_pri_n_p_5CH4", postprocess_dir = postprocess_dir, pad = 0.4);

    # ----------------------------------------------------------------------------------------------------------------------------------
    # Carbon Contribution from different processes
    # ----------------------------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize=(figsize[0]*2 - 0.5, figsize[1]*2 ));
    axes = axes.reshape(-1,)

    colors_set = ["#FE9D52", "#DDDDDD", "#2A6A99"];

    f_data_dict = {},
    for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
        letter_index = PAH_growth_model_type_index;
        title = f"({string.ascii_lowercase[letter_index]}) " + labels[PAH_growth_model_type];
        ax = axes[PAH_growth_model_type_index];
        ax.text(0.4, 0.87, title, transform=ax.transAxes, fontsize = letter_fontsize-1 , color = "#000000")

        header = "C_tot_int";
        val_dict = results_df[PAH_growth_model_type][header];
        f_data_tot = np.array(list(val_dict.values()));
        x_data = np.array(list(val_dict.keys()));

        header = "C_tot_inception_int";
        val_dict = results_df[PAH_growth_model_type][header];
        f_data_inc = np.array(list(val_dict.values()));

        header = "C_tot_PAH_int";
        val_dict = results_df[PAH_growth_model_type][header];
        f_data_PAH = np.array(list(val_dict.values()));

        header = "C_tot_growth_int";
        val_dict = results_df[PAH_growth_model_type][header];
        f_data_HACA = np.array(list(val_dict.values()));

        f_data_dict = {
            "Inception" : f_data_inc/f_data_tot,
            "PAH adsorption" : f_data_PAH/f_data_tot,
            "HACA" : f_data_HACA/f_data_tot,
        }

        ax.stackplot(x_data, f_data_dict.values(), labels = f_data_dict.keys(), colors = colors_set,
                    alpha = 0.7);


        # Y axis
        if PAH_growth_model_type_index % 2 == 0:
            ax.set_ylabel("$\\frac{\\mathrm{Soot\\;carbon\\;mass\\;from\\;each\\;source}}{\\mathrm{Soot\\;carbon\\;mass}}$",
                        {'fontsize' : axis_fontsize+2});
        ax.set_ylim([0, 1]);
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.04))

        # X axis
        if PAH_growth_model_type_index > 1:
            ax.set_xlabel("$\\mathrm{T_5}$ [K]", {'fontsize' : axis_fontsize+2})
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(40))
        ax.set_xlim([1800, 2900]);
        ax.tick_params(axis='x', pad = 10)


        legend_pos = [
            (0.027, 0.4),
            (0.027, 0.4)
        ];

        if PAH_growth_model_type_index == 3:
            ax.legend(loc = (0.36, 0.5), ncol = 1, fontsize = legend_fontsize-2, columnspacing = 0.1,
                            framealpha = 1, handletextpad = 0.3, handlelength = 0.8);

    save_plot("C_tot_distmap_5CH4", postprocess_dir = postprocess_dir, pad = 0.4);

if __name__ == "__main__":
    plot();