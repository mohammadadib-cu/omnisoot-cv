# Constant Pressure Reactor

This is an example for simulating a 5% CH<sub>4</sub>-Ar mixture using the pressure reactor model, the monodisperse population balance model, and four different inception models (aka `PAH_growth_model`s).

## File Structure

- **launcher.py**: This is the entry point of the example. Running the script loops through the provided list of PAH growth models calls `simulate` function to simulate reactor and writes the results in the `results` directory
- **simulate.py**: This script constains the `simulate` function that receives simulation parameters such as mechanism name, PAH growth model, particle dynamcis model as arguments, and then creates and configures the reactor model and steps forward in time until a given residence time.
- **plot.py**: This script contains the `plot` function, and plots (i) the residual of carbon and hydrogen mass and energy (ii) the number of agglomerates (iii) primary particle diameter, d<sub>p</sub> and specific surface area (SSA).
- **data**: The directory contains Caltech mechanism file.

## Prerequisites

- Python 3.7 or higher
- Required Python libraries:
    - `omnisoot`
    -  `cantera`
    - `matplotlib`
    - `scipy`
    - `numpy`

Install dependencies using:
```bash
pip install cantera omnisoot matplotlib scipy numpy
```