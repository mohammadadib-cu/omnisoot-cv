# Decelarations
import cantera as ct
import numpy as np
import scipy.integrate
from omnisoot.lib._omnisoot import CSootModel, CPFROde

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
    
    def check_mdot(self):
        if self.mdot is None:
            raise Exception("mdot of inlet is not set");
        elif self.mdot < 0.0:
            raise Exception("mdot of inlet is negative");

        return True;

    def check(self):
        if self.check_T() and self.check_mdot():
            return True;
        else:
            return False;

class CPFR(CPFROde):
    def __init__(self, phase):
        self.gas = phase;
        super().__init__(self.gas);

        self._res_time = 0;
        self._solved_length = 0;
        self._initailized = False;
        # Inlet
        self.inlet = Inlet(name = "reactor_inlet", phase = phase);

        # Solver Object
        self.solver = None;


        # Soot Properties
        self.turb_shear = False;
        self.wall_deposition = False;

        # Default PAH List
        self.PAH_species = ['A2', 'A3', 'A4', 'A2R5']


    def check_model(self):
        if len(self.PAH_species) == 0:
            raise ValueError("PAH species is not set!");
    

    def initialize(self):
        self._initailized = True;

        self.rho_u = self.inlet.mdot;
        self.gas.TPX = self.inlet.T, self.P, self.inlet.X;

        self.soot_enabled = True;

        # Initializing the soot model
        soot_model = CSootModel();
        soot_model.MWs = self.gas.molecular_weights /1000;
        
        self._set_PAH_species(soot_model);
        soot_model.species_names = self.gas.species_names

        # Enable/Disable differet stages of soot formation
        soot_model.carbonization_enabled = False;
        soot_model.inception_enabled = True;
        soot_model.growth_enabled = True;
        soot_model.coagulation_enabled = True;
        soot_model.oxidation_enabled = True;

        soot_model.steady_state_physical_dimerization = True;
        soot_model.steady_state_PAH_adsorption = True;

        # soot_model.use_empirical_alpha = True;
        # soot_model.use_unity_alpha = False;

        soot_model.modified_particle_addition = True;
        soot_model.update_gas(self.gas.T, self.gas.P, self.gas.density, self.gas.viscosity, self.gas.mean_molecular_weight, self.gas.X);
        soot_model.initialize();

        self.soot_model = soot_model;

        if self.soot_enabled:
            y0 = np.hstack((0.0, self.inlet.mdot / self.area, self.inlet.T, self.gas.Y, soot_model.N_agg, soot_model.N_pri, soot_model.C_tot, soot_model.H_tot));
        else:
            y0 = np.hstack((0.0, self.inlet.mdot / self.area, self.inlet.T, self.gas.Y));



        solver = scipy.integrate.ode(self);
        solver.set_integrator('vode', method='bdf', with_jacobian=True, max_step = 1e-3)
        solver.set_initial_value(y0, 0.0);
        self.solver = solver;


        # Setting initial conditions for the reactor
        self.gas.TPX = self.inlet.T, self.P, self.inlet.X;

    def advance(self, dz):
        if self.solver.successful():
            self.solver.integrate(self.solver.t + dz);
            self._solved_length += dz;

            if (self.soot_enabled):
                self.soot_model.update_solution(self.solver.y[self.soot_offset:])
    
    def _set_PAH_species(self, soot_model):
        PAH_indices = [self.gas.species_names.index(PAH) for PAH in self.PAH_species];
        PAH_n_C = [self.gas.n_atoms(PAH, 'C') for PAH in self.PAH_species];
        PAH_n_H = [self.gas.n_atoms(PAH, 'H') for PAH in self.PAH_species];

        soot_model.PAH_index_in_gas = PAH_indices;
        soot_model.PAH_Number_Carbon = PAH_n_C;
        soot_model.PAH_Number_Hydrogen = PAH_n_H;

    @property
    def PAH_list(self):
        return self.PAH_species;
    
    @PAH_list.setter
    def PAH_list(self, PAH_list):
        for PAH in PAH_list:
            if not PAH in self.gas.species_names:
                raise ValueError(f"{PAH} does not exist in the gas object");
        else:
            self.PAH_species = PAH_list;
    
    @property
    def res_time(self):
        return self.solver.y[self.Rt_offset]
    
    @property
    def T(self):
        return self.solver.y[self.T_offset]

    @property
    def gas_mdot(self):
        return self.solver.y[self.rho_u_vfv_offset]

    @property
    def solved_length(self):
        return self._solved_length;

    @property
    def d_p(self):
        return self.soot_model.d_p();

    @property
    def d_m(self):
        return self.soot_model.d_m();
    
    @property
    def d_g(self):
        return self.soot_model.d_g();

    @property
    def n_p(self):
        return self.soot_model.n_p();

    @property
    def N_agg(self):
        return self.soot_model.N_agg;
    
    @property
    def N_pri(self):
        return self.soot_model.N_pri;

    @property
    def A_tot(self):
        return self.soot_model.A_tot();
    
    @property
    def C_tot(self):
        return self.soot_model.C_tot;
    
    @property
    def H_tot(self):
        return self.soot_model.H_tot;

    @property
    def inception_rate(self):
        # return np.sum(self.soot_model.p_ROP);
        return self.soot_model.inception_N_agg();

    @property
    def HACA_growth(self):
        return self.soot_model.p_dC_tot_dt_HACA_alt;
    
    @property
    def PAH_growth(self):
        return self.soot_model.dC_tot_dt_PAH();

