from collections import defaultdict

import cantera as ct
import numpy as np
from numpy.testing import assert_allclose

from omnisoot import ConstantVolumeReactor, SootGas, SootWrapper, constants

from constants import (
    SMALL,
    CONSERVED_ATOL,
    CONSERVED_RTOL,
    INITIAL_D_P,
)


def record_common_data(reactor, data_list, soot_headers, particle_dynamics_model_headers):
    # Reactor headers
    data_list["t"].append(reactor.restime);
    for header in soot_headers:
        data_list[header].append(getattr(reactor.soot, header));
    for header in particle_dynamics_model_headers:
        data_list[header].append(getattr(reactor.soot.particle_dynamics_model, header)());


def test_pure_coagulation_monodisperes(gas):
    """
    Test monodisperse particle dynamics under pure coagulation in the free-molecular regime.
    The source terms for inception, surface growth, and oxidation must be
    initialized to zero and remain zero throughout the simulation.
    The quantities d_p, N_pri, and C_tot must also remain constant.
    """
    print("test_pure_coagulation_monodisperes executed!")
    soot_gas = SootGas(gas);
    ## ---------------------------------------------------------------------------------
    ## Intial conditions of the reactor
    T = 1830; # K
    P = 101325; # Pa
    X = {'CH4':0.425, 'O2':0.435 , 'N2':0.14};
    ## ---------------------------------------------------------------------------------
    gas_density = 0.123;
    ## Reactor
    reactor = ConstantVolumeReactor(soot_gas);
    ### Initial conditions
    soot_gas.TPX = T, P, X;
    ## Simulation time
    max_res_time = 0.5;
    ## Getting the handle to soot
    soot = reactor.soot;
    ## Setting particle dynamics model
    soot.particle_dynamics_model_type = "Monodisperse";
    particle_dynamics_model = soot.particle_dynamics_model;
    ## Modifying the soot options to deactivate inception, PAH adsorption, surface growth and oxidation
    soot.inception_enabled = False;
    soot.PAH_adsorption_enabled = False;
    soot.surface_growth_enabled = False;
    soot.oxidation_enabled = False;
    soot.coagulation_enabled = True;
    ## Modifying the initial soot of the reactor
    reactor.initial_soot_type = "custom";
    ## Calulating the soot variables based on the initial number of concentration of spherical particles 2 nm in diameters;
    N_0_number = 2626056626709460000; # 1/m3 from t=0 of DEM simulations
    N_0 = N_0_number / constants.Av / gas_density; # mol/kg
    # Calculating and constructing the soot array
    N_agg_0 = N_0;
    N_pri_0 = N_0;
    C_tot_0 = N_0 * constants.incepient_soot_C;
    H_tot_0 = N_0 * constants.incepient_soot_H;
    soot_array = np.hstack((N_agg_0, N_pri_0, C_tot_0, H_tot_0));
    reactor.user_defined_initial_soot = soot_array;
    # The collision occurs in the free molecular regime, so the collision frequency interpolation type is set to "free_molecule"
    particle_dynamics_model.beta_interp_method_type = "free_molecule";
    ## Reactor integrator settings
    reactor.max_step = 8e-4;
    reactor.solver_type = "BDF";
    reactor.max_time = max_res_time;
    ## The energy equations is not solved, so the temperature is kept constants
    reactor.temperature_solver_type = "isothermal";
    ## Starting the reactor
    reactor.start();
    ## Recording the simulation data
    soot_headers = [        
        "N_agg",
        "N_pri",
        "C_tot", "H_tot", 
        "n_p", "d_p", "d_m", "d_g",
    ]
    particle_dynamics_model_headers = [
        "I_C_tot_growth", "I_H_tot_growth", 
        "I_C_tot_PAH", "I_H_tot_PAH", 
        "I_C_tot_oxidation", "I_H_tot_oxidation", 
        "I_N_pri_coal",
    ];
    # Checking initial values after starting the reactor
    assert_allclose(
        soot.d_p, 
        INITIAL_D_P,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'Initial d_p is not equal to {INITIAL_D_P}'
    );
    for header in particle_dynamics_model_headers:
        assert_allclose(
            getattr(reactor.soot.particle_dynamics_model, header)(),
            0.0,
            rtol=CONSERVED_RTOL,
            atol=CONSERVED_ATOL,
            err_msg=f'Initial {header} is not equal to zero'
        )
        

    data_list = defaultdict(list);
    ## Running the reactor
    step = 0;
    while reactor.restime < max_res_time:
        step += 1;
        reactor.step();
        if step % 5 == 0:
            record_common_data(reactor, data_list, soot_headers, particle_dynamics_model_headers)

    ## Check if time-history of properties are valid
    assert_allclose(
        np.array(data_list["d_p"]), 
        INITIAL_D_P,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'd_p is not conserved during a pure coagulation simulation.'
    );
    assert_allclose(
        np.array(data_list["N_pri"]), 
        N_pri_0,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'N_pri is not conserved during a pure coagulation simulation.'
    );
    assert_allclose(
        np.array(data_list["C_tot"]), 
        C_tot_0,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'C_tot is not conserved during a pure coagulation simulation.'
    );
    for header in particle_dynamics_model_headers:
        assert_allclose(
            getattr(reactor.soot.particle_dynamics_model, header)(),
            0.0,
            rtol=CONSERVED_RTOL,
            atol=CONSERVED_ATOL,
            err_msg=f'{header} must remain zero during a pure coagulation simulation.'
        );


def test_morphology_monodisperes(gas):
    """
    Test morphology calculations in the monodisperse model.

    Initialize the model for two cases: nearly spherical particles (n_p=1)
    and large agglomerates (n_p=10). Verify the following relationships:

        n_p = N_pri / N_agg
        d_p = (6/Pi * C_tot * MW_carbon / rho_soot / N_pri / Av) ** (1/3)
        d_m = d_p * n_p**0.45
        d_g = d_m / 1.29                     for n_p <= 1.5
        d_g = d_m / (n_p**-0.2 + 0.4)        for n_p > 1.5
    """
    print("test_morphology_monodisperes executed!")

    soot_gas = SootGas(gas);
    soot = SootWrapper(soot_gas);

    ### Case 1
    N_0 = 3.545e-5; # [mol/kg]
    N_agg_0 = N_0;
    N_pri_0 = N_0;
    C_tot_0 = N_0 * constants.incepient_soot_C;
    H_tot_0 = N_0 * constants.incepient_soot_H;
    soot_array = np.hstack((N_agg_0, N_pri_0, C_tot_0, H_tot_0));
    ## Setting particle dynamics model
    soot.particle_dynamics_model_type = "Monodisperse";
    particle_dynamics_model = soot.particle_dynamics_model;
    particle_dynamics_model.set_vars(soot_array);
    particle_dynamics_model.update();

    ## Testing n_p
    n_p_omnisoot = particle_dynamics_model.n_p();
    n_p_test = N_pri_0/N_agg_0;
    assert_allclose(
        n_p_omnisoot, 
        n_p_test,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'n_p calculation is not valid.'
    );

    ## Testing d_p
    d_p_omnisoot = particle_dynamics_model.d_p();
    d_p_test = (6.0/constants.Pi * C_tot_0 * constants.MW_carbon / constants.rho_soot / N_pri_0 / constants.Av) ** (1.0/3.0);
    assert_allclose(
        d_p_omnisoot, 
        d_p_test,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'd_p calculation is not valid.'
    );

    ## Testing d_m
    d_m_omnisoot = particle_dynamics_model.d_m();
    d_m_test = d_p_test * n_p_test**0.45;
    assert_allclose(
        d_m_omnisoot, 
        d_m_test,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'd_m calculation is not valid.'
    );

    ## Testing d_g
    d_g_omnisoot = particle_dynamics_model.d_g();
    d_g_test = d_m_test / 1.29
    assert_allclose(
        d_g_omnisoot, 
        d_g_test,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'd_g calculation is not valid.'
    );


    ### Case 2
    N_agg_0 = N_0/10.0;
    soot_array = np.hstack((N_agg_0, N_pri_0, C_tot_0, H_tot_0));
    particle_dynamics_model.set_vars(soot_array);
    particle_dynamics_model.update();


    ## Testing d_g
    d_g_omnisoot = particle_dynamics_model.d_g();
    n_p_test = N_pri_0/N_agg_0;
    d_p_test = (6.0/constants.Pi * C_tot_0 * constants.MW_carbon / constants.rho_soot / N_pri_0 / constants.Av) ** (1.0/3.0);
    d_m_test = d_p_test * n_p_test**0.45;
    d_g_test = d_m_test / (n_p_test**-0.2 + 0.4);
    assert_allclose(
        d_g_omnisoot, 
        d_g_test,
        rtol=CONSERVED_RTOL,
        atol=CONSERVED_ATOL,
        err_msg=f'd_g calculation is not valid for n_p=10.'
    );