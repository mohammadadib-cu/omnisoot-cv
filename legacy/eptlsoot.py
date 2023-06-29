from omnisoot.lib._omnisoot import CSootModel, CPFROde, CPFRSootOde, CFlameSolver, CFlameSolverOpt

class SootModel(CSootModel):
    def __init__(self):
        super().__init__()


class Reactor(CPFROde):
    def __init__(self, phase):
        super().__init__(phase)

class ReactorSoot(CPFRSootOde):
    def __init__(self, phase):
        super().__init__(phase)

class FlameSolver(CFlameSolver):
    def __init__(self, phase, flame, grid, soot_model):
        super().__init__(phase, flame, grid, soot_model)


class FlameSolverOpt(CFlameSolverOpt):
    def __init__(self, phase, flame, grid, soot_model):
        super().__init__(phase, flame, grid, soot_model)