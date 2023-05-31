from ..lib._omnisoot import CSootWrapper
from .sootmodels import SOOT_MODELS
from .pahgrowth import PAH_GROWTH_MODELS
from .surfacereactions import SURFACE_REACTIONS_MODELS

class SootWrapper(CSootWrapper):
    def __init__(self, soot_gas):
        super().__init__(soot_gas);
        self._soot_model_dict = {};
        self._PAH_growth_model_dict = {};
        self._surface_reactions_model_dict = {};
        self.set_default_soot();
    
    def set_default_soot(self):
        # Default soot model
        soot_model = self._get_or_create_soot_model(SOOT_MODELS[0].serialized_name);
        self.set_soot_model(soot_model);
        # Default PAH growth model
        PAH_growth_model = self._get_or_create_PAH_growth_model(PAH_GROWTH_MODELS[0].serialized_name);
        self.set_PAH_growth_model(PAH_growth_model);
        # Default surface reactions model
        surface_reactions_model = self._get_or_create_surface_reactions_model(SURFACE_REACTIONS_MODELS[0].serialized_name);
        self.set_surface_reactions_model(surface_reactions_model);

    # -----------------------------------------------------------------------------------------
    # Soot Model            
    def _get_or_create_soot_model(self, soot_model_name):
        if soot_model_name in self._soot_model_dict.keys():
            soot_model = self._soot_model_dict.get(soot_model_name);
            return soot_model;
        else:
            soot_model_registed_names = [model.serialized_name for model in SOOT_MODELS];
            if soot_model_name in soot_model_registed_names:
                SootModel = SOOT_MODELS[soot_model_registed_names.index(soot_model_name)]
                self._soot_model_dict[soot_model_name] = SootModel(self);
                return self._soot_model_dict[soot_model_name];
            else:
                raise ValueError(f"{soot_model_name} is not registered");

    @property
    def soot_model_type(self):
        return self.soot_model.serialized_name;

    @soot_model_type.setter
    def soot_model_type(self, model_type):
        current_soot_model = self._get_or_create_soot_model(model_type);
        self.set_soot_model(current_soot_model);
    
    
    # -----------------------------------------------------------------------------------------
    # PAH Growth Model       
    def _get_or_create_PAH_growth_model(self, PAH_growth_model_name):
        if PAH_growth_model_name in self._PAH_growth_model_dict.keys():
            PAH_growth_model = self._PAH_growth_model_dict.get(PAH_growth_model_name);
            return PAH_growth_model;
        else:
            PAH_growth_model_registed_names = [model.serialized_name for model in PAH_GROWTH_MODELS];
            if PAH_growth_model_name in PAH_growth_model_registed_names:
                PAHGrowthModel = PAH_GROWTH_MODELS[PAH_growth_model_registed_names.index(PAH_growth_model_name)]
                self._PAH_growth_model_dict[PAH_growth_model_name] = PAHGrowthModel(self);
                return self._PAH_growth_model_dict[PAH_growth_model_name];
            else:
                raise ValueError("The PAH growth model name is not registered")


    @property
    def PAH_growth_model_type(self):
        return self.PAH_growth_model.serialized_name;

    @PAH_growth_model_type.setter
    def PAH_growth_model_type(self, model_type):
        current_PAH_growth_model = self._get_or_create_PAH_growth_model(model_type);
        self.set_PAH_growth_model(current_PAH_growth_model);
    
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
    
    # -----------------------------------------------------------------------------------------
    # Surface Reactions Model    
    def _get_or_create_surface_reactions_model(self, surface_reactions_model_name):
        if surface_reactions_model_name in self._surface_reactions_model_dict.keys():
            surface_reactions_model = self._surface_reactions_model_dict.get(surface_reactions_model_name);
            return surface_reactions_model;
        else:
            surface_reactions_registered_names = [model.serialized_name for model in SURFACE_REACTIONS_MODELS];
            if surface_reactions_model_name in surface_reactions_registered_names:
                SurfaceReactionsModel = SURFACE_REACTIONS_MODELS[surface_reactions_registered_names.index(surface_reactions_model_name)];
                self._surface_reactions_model_dict[surface_reactions_model_name] = SurfaceReactionsModel(self);
                return self._surface_reactions_model_dict[surface_reactions_model_name];
            else:
                raise ValueError("The surface reactions model name is not registered")
                

    @property
    def surface_reactions_model_type(self):
        return self.surface_reactions_model.serialized_name;

    @surface_reactions_model_type.setter
    def surface_reactions_model_type(self, model_type):
        current_surface_reactions_model = self._get_or_create_surface_reactions_model(model_type);
        self.set_surface_reactions_model(current_surface_reactions_model);

    def __getattr__(self, name):
        if name in self.soot_model.soot_att:
            return self.soot_model.det_soot_att(name);
        else:
            raise ValueError(f"{name} is an attirbute of SootWrapper class")
