# Third-party libraries
from scipy import integrate
import numpy as np
# Local imports
from omnisoot.lib._omnisoot import CPFRSoot, CConstUVSootHighCon, CConstUVSootLowCon, CPSRSoot
from omnisoot.apps.sootwrappers import SootWrapper
from omnisoot.apps.sootgas import SootGas

class ReactorAbstract:
    def __init__(self, soot_gas):
        super().__init__(soot_gas);

    def create_soot_wrapper(self):
        soot_wrapper = SootWrapper(self.soot_gas);
        self.set_soot_wrapper(soot_wrapper);


    def step(self):
        self.solver.step();
    
    def check_soot_array(self, soot_array):
        check = False;
        if isinstance(soot_array, np.ndarray):
            if soot_array.shape == (self.soot_wrapper.soot_model.n_eqns,):
                check = True;
        if not check:
            raise TypeError("Wrong inlet soot array");
        return check;

    @property
    def soot(self):
        return self.soot_wrapper;


class PlugFlowReactor(ReactorAbstract, CPFRSoot):
    def __init__(self, soot_gas):
        super().__init__(soot_gas);
        self.create_soot_wrapper();
        # Reactor pressure
        self.P = self.soot_gas.P;
        # Inlet
        self.inlet = Inlet(self);
        # Solver
        self.solver = None;
        # Max length
        self.max_length = 100;
        # First step
        self.first_step = None;
        # Max step
        self.max_step = 1e-3;
        # Tol
        self.rtol= 1e-5;
        self.atol= 1e-10;

    def start(self):
        if self.check_soot_array(self.inlet.soot):
            self.T = self.inlet.T;
            self.P = self.inlet.P;
            self.set_mdot(self.inlet.mdot);
            y0 = np.hstack((0.0, self.inlet.mdot, self.inlet.T, self.inlet.Y, self.inlet.soot));
            self.solver = integrate.LSODA(fun=self.derivatives, t0=0, y0=y0, t_bound=self.max_length,
                                        first_step=self.first_step, max_step = self.max_step,
                                        rtol=self.rtol, atol=self.atol);
            

    @property
    def total_carbon_flux(self) -> float:
        carbon_mass = self.soot_gas.elemental_mass_fraction('C') + self.soot_wrapper.soot_model.carbon_mass();
        return carbon_mass * self.mdot;
    
    @property
    def gas_carbon_flux(self) -> float:
        return self.soot_gas.elemental_mass_fraction('C') * self.mdot;

    @property
    def soot_carbon_flux(self) -> float:
        return self.soot_wrapper.soot_model.carbon_mass() * self.mdot;


    @property
    def total_hydrogen_flux(self) -> float:
        hydrogen_mass = self.soot_gas.elemental_mass_fraction('H') + self.soot_wrapper.soot_model.hydrogen_mass();
        return hydrogen_mass * self.mdot;
    
    @property
    def gas_hydrogen_flux(self) -> float:
        return self.soot_gas.elemental_mass_fraction('H') * self.mdot;

    @property
    def soot_hydrogen_flux(self) -> float:
        return self.soot_wrapper.soot_model.hydrogen_mass() * self.mdot;

    @property
    def gas_elemental_flux(self, element_name) -> float:
        return self.soot_gas.elemental_mass_fraction(element_name) * self.mdot;


class Inlet:
    def __init__(self, reactor, mdot = 0, soot = "zero_soot"):
        self.reactor = reactor;
        self._Y = self.reactor.soot_gas.X;
        self._X = self.reactor.soot_gas.Y;
        self.T = self.reactor.soot_gas.T;
        self.P = self.reactor.soot_gas.P;
        self.mdot = mdot;
        if soot == "zero_soot":
            self.soot = self.reactor.soot_wrapper.min_array;

    @property
    def X(self):
        return self._X;

    @X.setter
    def X(self, X):
        soot_gas = self.reactor.soot_gas;
        soot_gas.X = X;
        self._X = soot_gas.X
        self._Y = soot_gas.Y

    @property
    def Y(self):
        return self._Y;

    @Y.setter
    def Y(self, Y):
        soot_gas = self.reactor.soot_gas;
        soot_gas.Y = Y;
        self._Y = soot_gas.Y
        self._X = soot_gas.X

    @property
    def TPX(self):
        return self.T, self.P, self._X;

    @TPX.setter
    def TPX(self, TPX):
        soot_gas = self.reactor.soot_gas;
        soot_gas.TPX = TPX;
        self._Y = soot_gas.Y;
        self._X = soot_gas.X;

    @property
    def TPY(self):
        return self.T, self.P, self._Y;

    @TPY.setter
    def TPY(self, TPY):
        soot_gas = self.reactor.soot_gas;
        soot_gas.TPY = TPY;
        self._Y = soot_gas.Y
        self._X = soot_gas.X



class ConstantVolumeReactorMixin:
    def start(self):
        if self.check_soot_array(self.initial_soot):
            y0 = np.hstack((self.soot_gas.rho, self.soot_gas.T, self.soot_gas.Y_array, self.initial_soot));
            self.solver = integrate.LSODA(fun=self.derivatives, t0=0, y0=y0, t_bound=self.max_time,
                                        first_step=self.first_step, max_step = self.max_step,
                                        rtol=self.rtol, atol=self.atol);


    @property
    def total_carbon_mass(self) -> float:
        carbon_mass = self.soot_gas.elemental_mass_fraction('C') + self.soot_wrapper.soot_model.carbon_mass();
        return carbon_mass;

    @property
    def gas_carbon_mass(self) -> float:
        return self.soot_gas.elemental_mass_fraction('C');

    @property
    def soot_carbon_mass(self) -> float:
        return self.soot_wrapper.soot_model.carbon_mass();

    @property
    def total_hydrogen_mass(self) -> float:
        hydrogen_mass = self.soot_gas.elemental_mass_fraction('H') + self.soot_wrapper.soot_model.hydrogen_mass();
        return hydrogen_mass;

    @property
    def gas_hydrogen_mass(self) -> float:
        return self.soot_gas.elemental_mass_fraction('H');

    @property
    def soot_hydrogen_mass(self) -> float:
        return self.soot_wrapper.soot_model.hydrogen_mass();

    @property
    def gas_elemental_mass(self, element_name) -> float:
        return self.soot_gas.elemental_mass_fraction(element_name);

class ConstantVolumeReactorHighConcentration(ConstantVolumeReactorMixin, ReactorAbstract, CConstUVSootHighCon):
    def __init__(self, soot_gas):
        super().__init__(soot_gas);
        # super(ReactorAbstract, self).__init__();
        self.create_soot_wrapper();
        # Reactor pressure
        self.P = self.soot_gas.P;
        # Solver
        self.solver = None;
        # Max residence time
        self.max_time = 100;
        # First step
        self.first_step = None;
        # Max step
        self.max_step = 1e-3;
        # Tol
        self.rtol = 1e-7;
        self.atol = 1e-12;
        # Initial soot
        self.initial_soot = self.soot_wrapper.min_array;
        # Temperature solver
        self.temperature_solver_type = "energy_equation";



class ConstantVolumeReactorLowConcentration(ConstantVolumeReactorMixin, ReactorAbstract, CConstUVSootLowCon):
    def __init__(self, soot_gas):
        super().__init__(soot_gas);
        # super(ReactorAbstract, self).__init__();
        self.create_soot_wrapper();
        # Reactor pressure
        self.P = self.soot_gas.P;
        # Solver
        self.solver = None;
        # Max residence time
        self.max_time = 100;
        # First step
        self.first_step = None;
        # Max step
        self.max_step = 1e-3;
        # Tol
        self.rtol= 1e-7;
        self.atol= 1e-12;
        # Initial soot
        self.initial_soot = self.soot_wrapper.min_array;
        # Temperature solver
        self.temperature_solver_type = "energy_equation";

class ConstantVolumeReactor:
    def __new__(cls, soot_gas, high_concentration = False):
        if high_concentration:
            return ConstantVolumeReactorHighConcentration(soot_gas);
        else:
            return ConstantVolumeReactorLowConcentration(soot_gas);


class PerfectlyStirredReactor(ReactorAbstract, CPSRSoot):
    def __init__(self, soot_gas):
        super().__init__(soot_gas);
        self.create_soot_wrapper();
        # Reactor pressure
        self.P = self.soot_gas.P;
        # Solver
        self.solver = None;
        # Max residence time
        self.max_time = 100;
        # First step
        self.first_step = 1e-10;
        # Max step
        self.max_step = 1e-4;
        # Initial soot
        self.initial_soot = self.soot_wrapper.min_array;
        # Temperature solver
        self.temperature_solver_type = "isothermal";
        # External Heat
        self.Q_dot = 0;

        self.set_star();
        
    def set_modt(self):
        self.mdot = self.soot_gas.rho * self.reactor_volume / self.phy_restime;
    
    def check_params(self) -> bool:
        valid = True;
        if self.residence_time <= 0:
            valid = False;
            raise TypeError("Wrong residence time value!");
        
        if self.reactor_volume <= 0:
            valid = False;
            raise TypeError("Wrong reactor volume value!");            

        return valid;

    def start(self):
        if self.check_soot_array(self.initial_soot) and self.check_params():
            self.set_modt();
            y0 = np.hstack((self.soot_gas.rho, self.soot_gas.T, self.soot_gas.Y_array, self.initial_soot));
            self.solver = integrate.LSODA(fun=self.derivatives, t0=0, y0=y0, t_bound=self.max_time,
                                        first_step=self.first_step, max_step = self.max_step);

    @property
    def total_carbon_flux(self) -> float:
        carbon_mass = self.soot_gas.elemental_mass_fraction('C') + self.soot_wrapper.soot_model.carbon_mass();
        return carbon_mass * self.mdot * (1.0-self.soot_wrapper.soot_model.volume_fraction());
    
    @property
    def gas_carbon_flux(self) -> float:
        return self.soot_gas.elemental_mass_fraction('C') * self.mdot * (1.0-self.soot_wrapper.soot_model.volume_fraction());

    @property
    def soot_carbon_flux(self) -> float:
        return self.soot_wrapper.soot_gas.elemental_mass_fraction('C') * self.mdot * (1.0-self.soot_wrapper.soot_model.volume_fraction());


    @property
    def total_hydrogen_flux(self) -> float:
        hydrogen_mass = self.soot_gas.elemental_mass_fraction('H') + self.soot_wrapper.soot_model.hydrogen_mass();
        return hydrogen_mass * self.mdot * (1.0-self.soot_wrapper.soot_model.volume_fraction());
    
    @property
    def gas_hydrogen_flux(self) -> float:
        return self.soot_gas.elemental_mass_fraction('H') * self.mdot * (1.0-self.soot_wrapper.soot_model.volume_fraction());

    @property
    def soot_hydrogen_flux(self) -> float:
        return self.soot_wrapper.soot_gas.elemental_mass_fraction('H') * self.mdot * (1.0-self.soot_wrapper.soot_model.volume_fraction());

    @property
    def gas_elemental_flux(self, element_name) -> float:
        return self.soot_gas.elemental_mass_fraction(element_name) * self.mdot * (1.0-self.soot_wrapper.soot_model.volume_fraction());


    @property
    def residence_time(self) -> float:
        return self.phy_restime;

    @residence_time.setter
    def residence_time(self, residence_time):
        self.phy_restime = residence_time;

