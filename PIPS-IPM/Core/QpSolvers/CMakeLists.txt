add_library(ooqpmehrotra MehrotraSolver.C)
add_library(ooqpmehrotrastoch MehrotraStochSolver.C MehrotraSolver.C)
add_library(ooqpgondzio GondzioSolver.C)
add_library(ooqpgondziostoch GondzioStochSolver.C GondzioSolver.C GondzioStochLpSolver.C)



