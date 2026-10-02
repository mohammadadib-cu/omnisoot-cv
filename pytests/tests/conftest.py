import cantera as ct
import pytest

from constants import (
    MECHANISM_PATH,
    TEST_EXECUSION_ORDER,
)

@pytest.fixture(scope="session")
def gas():
    # Loading the mechanism is computationally expensive
    # It is expected to run once throughout the test session
    return ct.Solution(MECHANISM_PATH)



def deactive_pytest_collection_modifyitems(session, config, items):
    print("pytest_collection_modifyitems started")
    misc_tests = [item for item in items if not item.name in TEST_EXECUSION_ORDER];
    main_tests = [item for item in items if item.name in TEST_EXECUSION_ORDER];
    main_tests_ordered = [];
    for test_name in TEST_EXECUSION_ORDER:
        items_found = [item for item in main_tests if item.name == test_name];
        if len(items_found) == 1:
            main_tests_ordered.append(items_found[0])
    for item in (main_tests_ordered+misc_tests):
        print(item.path.name)
    print("pytest_collection_modifyitems ended")