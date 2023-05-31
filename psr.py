# Decelarations
import cantera as ct
import numpy as np
import scipy.integrate
from .lib._omnisoot import CSootModel


class Inlet:
    def __init__(self, name, phase):
        self.phase = phase;
        self.name = name;
        self.X = {};
        self.T = None; #[K]
        self.mdot = None;

    def check_T(self):
        if self.T is None:
            raise Exception("Temperature of inlet is not set");
        elif self.T < 0.0:
            raise Exception("Temperature of inlet is too low");
        elif self.T > 3000.0:
            raise Exception("Temperature of inlet is too high");   
        return True;