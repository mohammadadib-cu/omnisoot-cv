from pathlib import Path

import cantera as ct
import numpy as np
from omnisoot import ConstantVolumeReactor, PressureReactor, PerfectlyStirredReactor, PlugFlowReactor, SootGas

from constants import (
    PRECURSOR_NAMES,
    INITIAL_COMPOSITION,
    NUMBER_OF_TEST_STEPS
)

def test_constant_volume_reactor(gas):
    soot_gas = SootGas(gas);
    SootGas.TPX = INITIAL_COMPOSITION, ct.one_atm, 2000;
    reactor = ConstantVolumeReactor(soot_gas);
    soot = reactor.soot;
    soot.set_precursor_names(PRECURSOR_NAMES);
    reactor.start();
    print("ConstantVolumeReactor started!")
    for i in range(NUMBER_OF_TEST_STEPS):
        reactor.step();
    print(f"ConstantVolumeReactor is integrated for {NUMBER_OF_TEST_STEPS} steps")

def test_pressure_reactor(gas):
    soot_gas = SootGas(gas);
    SootGas.TPX = INITIAL_COMPOSITION, ct.one_atm, 2000;
    reactor = PressureReactor(soot_gas);
    reactor.pressure_type = "constant";
    soot = reactor.soot;
    soot.set_precursor_names(PRECURSOR_NAMES);
    reactor.start();
    print("PressureReactor started!")
    for i in range(NUMBER_OF_TEST_STEPS):
        reactor.step();
    print(f"PressureReactor is integrated for {NUMBER_OF_TEST_STEPS} steps in the constant pressure mode");

    SootGas.TPX = INITIAL_COMPOSITION, ct.one_atm, 2000;
    reactor = PressureReactor(soot_gas);
    reactor.pressure_type = "imposed";
    pressure_profile = np.vstack(
        (np.array([0.0, 1e-3, 2e-3, 10e-3]),
        np.array([4*101350, 3.5*101350, 3*101350, 2.7*101350]))
    );
    reactor.set_fix_pressure_profile(pressure_profile);
    reactor.start();
    for i in range(NUMBER_OF_TEST_STEPS):
        reactor.step();
    print(f"PressureReactor is integrated for {NUMBER_OF_TEST_STEPS} steps in the imposed pressure mode");