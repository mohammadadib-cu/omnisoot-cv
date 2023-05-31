import cantera as ct

class SolutionArray(ct.SolutionArray):
    def __init__(self, soot_gas, shape=(0,), states=None, extra=None, meta=None):
        phase = soot_gas.cantera_gas;
        super().__init__(phase, shape=(0,), states=None, extra=None, meta=None)