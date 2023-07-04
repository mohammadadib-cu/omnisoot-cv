# omnisoot

omnisoot (previously eptlsoot) is an object-oriented computational tool for the simulation of the gas phase synthesis of nanoparticles including soot (Carbon Black), graphene, nickel, and carbon nanotubes in flames and reactors. omnisoot can be used for scientific research or process design and optimization of commercial reactors. It accommodates a variety of particle dynamics, inception, growth, and oxidation models that can be combined and used for target simulations in Premixed and Counterflow diffusion flames as well as Plug Flow, Constant Volume, Partially Stirred, and Shock Tube reactors. Reaction mechanisms are processed by [Cantera](https://cantera.org/) that is used to compute thermodynamic and kinetics properties of the gas mixture, and [SciPy](https://scipy.org/) ode package is used in reactor. The code is designed to handle high particle volume fractions. omnisoot can be used from Python and a web interface which provides post-processing and visualization tools. A series of examples are also available in the form of Jupyter notebooks and Python scripts.
omnisoot has been registered for copyright on 31 May 2023 by [Canadian Intellectual Property Office (CIPO)](https://www.canada.ca/en/services/business/ip.html) under regirstration number 1203139.

### updates
v0.1.7:
 - Performance improvements for reactors
 - Temperature time history capability is enabled for Constatn Volume Reactor

v0.1.6:
 - The formation and sensible energy of soot is tracked in the energy equation to close the energy gap in adiabatic processes. 
 - A new class, SootThermo, is added to the package to calculate the heat capacity and internal energy of soot is based on the properties of  [graphite](https://github.com/Cantera/cantera/blob/main/data/graphite.yaml). 
 - Elemental mass balance is improved.