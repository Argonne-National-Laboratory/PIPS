add_library(ooqpqpscaler QpScaler.C)
add_library(ooqpequistochscaler EquiStochScaler.C GeoStochScaler.C QpScaler.C)
add_library(ooqpqppresolver QpPresolver.C)
#add_library(  StochPresolverSingletonColumns.C )
add_library(ooqpqpstochpresolver QpPostsolver.C StochPostsolver.C StochPresolver.C QpPresolver.C StochPresolverBase.C
 StochPresolverModelCleanup.C PresolveData.C StochPresolverParallelRows.C StochPresolverBoundStrengthening.C StochPresolverColumnFixation.C)




