from ..lib._omnisoot import CFrenklachHACA

class FrenklachHACA(CFrenklachHACA):
    serialized_name = "FrenklachHACA"
    def __init__(self, soot_gas):
        super().__init__(soot_gas);

SURFACE_REACTIONS_MODELS = [FrenklachHACA]