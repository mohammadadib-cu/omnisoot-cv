from ..lib._omnisoot import h_mass_soot, u_mass_soot, cpv_mass_soot, h_mole_soot, u_mole_soot, cpv_mole_soot

class SootThermo:
    @staticmethod
    def u_mass_soot(T: float) -> float:
        return u_mass_soot(T);

    @staticmethod
    def u_mole_soot(T: float) -> float:
        return u_mole_soot(T);

    @staticmethod
    def h_mass_soot(T: float, P: float) -> float:
        return h_mass_soot(T);

    @staticmethod
    def h_mole_soot(T: float, P: float) -> float:
        return h_mole_soot(T);

    @staticmethod
    def cpv_mass_soot(T: float) -> float:
        return cpv_mass_soot(T);

    @staticmethod
    def cpv_mol_soot(T: float) -> float:
        return cpv_mole_soot(T);