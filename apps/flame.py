from ..lib._omnisoot import CTempFlameSolver, CTempFlameSolverOpt, CFVSolver, CFDSolver, CFDSolverTemp

class TempFlameSolver(CTempFlameSolver):
    def __init__(self, soot_gas, flame, grid, soot_model):
        super().__init__(soot_gas, flame, grid, soot_model)


# class TempFlameSolverVector(CTempFlameSolverVector):
#     def __init__(self, soot_gas, flame, grid, soot_model):
#         super().__init__(soot_gas, flame, grid, soot_model)


class TempFlameSolverOpt(CTempFlameSolverOpt):
    def __init__(self, soot_gas, flame, grid, soot):
        super().__init__(soot_gas, flame, grid, soot)   

class FVSolver(CFVSolver):
    def __init__(self, soot_gas, flame, grid, soot):
        super().__init__(soot_gas, flame, grid, soot)

class FDSolver(CFDSolver):
    def __init__(self, soot_gas, flame, grid, soot):
        super().__init__(soot_gas, flame, grid, soot)   

class FDSolverTemp(CFDSolverTemp):
    def __init__(self, soot_gas, flame, grid, soot):
        super().__init__(soot_gas, flame, grid, soot)   