from ..lib._omnisoot import CReactDim, CDimerCoal, CCrossLink, CCrossLinkMod, CCrossLinkMerge, CIrrevDim, CEBridge

class PAHGrowthAbstract:
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);


    def set_precursor_names(self, precursor_list):
        PAH_indices = [];
        PAH_n_C = [];
        PAH_n_H = [];
        cantera_gas = self.soot_gas.cantera_gas;
        for precursor in precursor_list:
            if precursor in cantera_gas.species_names:
                PAH_indices.append(cantera_gas.species_names.index(precursor));
                PAH_n_C.append(cantera_gas.n_atoms(precursor, 'C'));
                PAH_n_H.append(cantera_gas.n_atoms(precursor, 'H'));
            else:
                raise ValueError(f"{precursor} does not exist in cantera gas object!");

        self.set_precursors(PAH_indices, PAH_n_C, PAH_n_H);

class ReactDim(PAHGrowthAbstract, CReactDim):
    serialized_name = "ReactiveDimerization"
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);


class DimerCoal(PAHGrowthAbstract, CDimerCoal):
    serialized_name = "DimerCoalescence"
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);


class CrossLink(PAHGrowthAbstract, CCrossLink):
    serialized_name = "CrossLinking"
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);

class CrossLinkMod(PAHGrowthAbstract, CCrossLinkMod):
    serialized_name = "CrossLinkingModified"
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);


class CrossLinkMerge(PAHGrowthAbstract, CCrossLinkMerge):
    serialized_name = "CrossLinkMerge"
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);

class IrrevDim(PAHGrowthAbstract, CIrrevDim):
    serialized_name = "IrreversibleDimerization"
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);

class EBridge(PAHGrowthAbstract, CEBridge):
    serialized_name = "EBridgeFormation"
    def __init__(self, soot_wrapper):
        super().__init__(soot_wrapper);

PAH_GROWTH_MODELS = [ReactDim, DimerCoal, CrossLink, CrossLinkMod, CrossLinkMerge, IrrevDim, EBridge]