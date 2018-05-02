add_library(ooqpqpscaler QpScaler.C)
add_library(ooqpequistochscaler EquiStochScaler.C QpScaler.C)
add_library(ooqpqppresolver QpPresolver.C)
add_library(ooqpqpstochpresolver StochPresolver.C QpPresolver.C StochPresolverBase.C StochPresolverTinyEntries.C StochPresolverSingletonRows.C StochPresolverDuplicateRows.C PresolveData.C)




