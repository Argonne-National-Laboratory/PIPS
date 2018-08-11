add_library(ooqpqpscaler QpScaler.C)
add_library(ooqpequistochscaler EquiStochScaler.C GeoStochScaler.C QpScaler.C)
add_library(ooqpqppresolver QpPresolver.C)
add_library(ooqpqpstochpresolver StochPresolver.C QpPresolver.C StochPresolverBase.C StochPresolverTinyEntries.C StochPresolverSingletonRows.C StochPresolverParallelRows.C StochPresolverSingletonColumns.C PresolveData.C StochPresolverBoundStrengthening.C)




