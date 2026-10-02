from omnisoot import SootGas


def test_create_soot_gas(gas):
    soot_gas = SootGas(gas);
    assert not soot_gas is None;
    print("Soot gas is successfully loaded")
    