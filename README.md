
# Omnisoot

This repository contains the validation and uses cases of omnisoot in Python, which is an object-oriented computational tool for the simulation of the gas phase synthesis of carbonaceous nano-particles in flames and reactors. Omnisoot can be used for scientific research or process design and optimization of commercial reactors. It accommodates a variety of particle dynamics, inception, growth, and oxidation models that can be combined and used for target simulations in Plug Flow, Constant Volume, Partially Stirred, and Pressure reactors. Reaction mechanisms are processed by [Cantera](https://cantera.org/) that is used to compute thermodynamic and kinetics properties of the gas mixture, and [SciPy](https://scipy.org/) ode package is used in reactor. The code is designed to handle high particle volume fractions. omnisoot can be used from Python and a web interface which provides post-processing and visualization tools.

  

omnisoot has been registered for copyright on 31 May 2023 by [Canadian Intellectual Property Office (CIPO)](https://www.canada.ca/en/services/business/ip.html) under regirstration number 1203139.

  

# Installation

Omnisoot is a standard Python package published on [PyPi](https://pypi.org/project/omnisoot/) that can be installed by `pip install omnisoot`

# Theoretical Foundation
The mathematical basis of omnisoot including the transport equations of the reactors, particle dynamics models, inception and surface growth models are explained in the corresponding [paper](https://github.com/mohammadadib-cu/omnisootpaper).