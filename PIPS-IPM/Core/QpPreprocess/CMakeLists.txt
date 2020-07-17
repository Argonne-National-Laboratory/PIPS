add_library(ooqpqpscaler QpScaler.C)
add_library(ooqpequistochscaler EquiStochScaler.C GeoStochScaler.C QpScaler.C StochScaler.C)
add_library(ooqpqppresolver QpPresolver.C)
add_library(ooqppqpostsolver QpPostsolver.C )
#add_library(   StochPresolverSingletonColumns.C )
add_library(ooqpqpstochpresolver StochRowStorage.C StochColumnStorage.C QpPostsolver.C StochPostsolver.C StochPresolver.C QpPresolver.C StochPresolverBase.C
 StochPresolverModelCleanup.C PresolveData.C StochPresolverBoundStrengthening.C StochPresolverColumnFixation.C StochPresolverParallelRows.C
 StochPresolverSingletonRows.C StochPresolverSingletonColumns.C StochPresolverDualFixing.C )
