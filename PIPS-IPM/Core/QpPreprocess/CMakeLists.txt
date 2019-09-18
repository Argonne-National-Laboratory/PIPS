add_library(ooqpqpscaler QpScaler.C)
add_library(ooqpequistochscaler EquiStochScaler.C GeoStochScaler.C QpScaler.C StochScaler.C)
add_library(ooqpqppresolver QpPresolver.C)
add_library(ooqpqpstochpresolver StochPostsolver.C QpPostsolver.C StochPresolver.C QpPresolver.C StochPresolverBase.C StochPresolverModelCleanup.C StochPresolverSingletonRows.C StochPresolverParallelRows.C StochPresolverSingletonColumns.C PresolveData.C StochPresolverBoundStrengthening.C)




