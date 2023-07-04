import numpy as np
import cantera as ct
from ..lib._omnisoot import CSootGas

class SootGas(CSootGas):
    key_species_list = ["C2H2", "H", "H2", "O2", "OH", "H2O", "CO"];
    MW_carbon = 12.011e-3 # [kg/mol]
    MW_hydrogen = 1.008e-3 # [kg/mol]
    def __init__(self, cantera_gas):
        super().__init__(cantera_gas);
        self.init_by_cantera(cantera_gas)

    # Molecular weights
    @property
    def molecular_weights(self):
        return self._MWs;
        
    @property
    def mean_molecular_weight(self):
        return self._mean_MW;
    
    # TPX
    @property
    def TPX(self):
        return self.cantera_gas.T, self.cantera_gas.P, self.cantera_gas.X;

    @TPX.setter
    def TPX(self, TPX):
        self.cantera_gas.TPX = TPX[0], TPX[1], TPX[2];
        self.update_thermo();
        self.update_transport();

    # TPY
    @property
    def TPY(self):
        return self.cantera_gas.T, self.cantera_gas.P, self.cantera_gas.Y;

    @TPY.setter
    def TPY(self, TPY):
        self.cantera_gas.TPY = TPY[0], TPY[1], TPY[2];
        self.update_thermo();
        self.update_transport();
    
    # TDX
    @property
    def TDX(self):
        return self.cantera_gas.T, self.cantera_gas.density, self.cantera_gas.X;

    @TDX.setter
    def TDX(self, TDX):
        self.cantera_gas.TDX = TDX[0], TDX[1], TDX[2];
        self.update_thermo();
        self.update_transport();

    # TDX
    @property
    def TDY(self):
        return self.cantera_gas.T, self.cantera_gas.density, self.cantera_gas.Y;

    @TDY.setter
    def TDY(self, TDY):
        self.cantera_gas.TDY = TDY[0], TDY[1], TDY[2];
        self.update_thermo();
        self.update_transport();
       
    # X
    @property
    def X(self):
        return self.cantera_gas.X;

    @X.setter
    def X(self, X):
        self.cantera_gas.X = X;
        self.update_thermo();
        self.update_transport();
   
    # X
    @property
    def Y(self):
        return self.cantera_gas.Y;

    @Y.setter
    def Y(self, Y):
        self.cantera_gas.Y = Y;
        self.update_thermo();
        self.update_transport();
   

    @property
    def density(self):
        return self.rho;

    # Transport 
    @property
    def viscosity(self):
        return self.mu;
  
    @property
    def mean_free_path(self):
        return self.lambda_gas;

    @property
    def state(self):
        return self.cantera_gas.state;
    

    def init_by_cantera(self, cantera_gas):
        self.update_MWs(cantera_gas.molecular_weights/1000.0);
        self.update_MW_carbon_array(self.get_no_element_array("C", cantera_gas)*self.MW_carbon);
        self.update_MW_hydrogen_array(self.get_no_element_array("H", cantera_gas)*self.MW_hydrogen);
        self.update_thermo_TPX(cantera_gas.T, cantera_gas.P, cantera_gas.X);
        self.update_transport();
        speices_indices_dict = self.build_speices_indices_dict(cantera_gas);
        self.set_species_indices(speices_indices_dict);
    
    def get_no_element_array(self, el, cantera_gas):
        gas = cantera_gas;
        no_element_array = np.zeros(gas.n_species);
        for i in range(gas.n_species):
            no_element_array[i] = gas.n_atoms(gas.species_names[i], el);
        return no_element_array;
    
    def build_speices_indices_dict(self, cantera_gas):
        speices_indices_dict = {};
        for species in self.key_species_list:
            if species in cantera_gas.species_names:
                speices_indices_dict[bytes(species, "ascii")] = cantera_gas.species_names.index(species);
            else:
                raise ValueError(f"{species} does not exist in cantera gas object!");
        return speices_indices_dict;

        
    def carbon_mass(self):
        gas = self.cantera_gas;
        return gas.elemental_mass_fraction('C')
    
    def elemental_mass_fraction(self, element_name):
        gas = self.cantera_gas;
        return gas.elemental_mass_fraction(element_name)
    
    def elemental_mole_fraction(self, element_name):
        gas = self.cantera_gas;
        return gas.elemental_mole_fraction(element_name)



