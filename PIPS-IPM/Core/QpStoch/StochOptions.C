#include "StochOptions.h"

StochOptions::StochOptions()
{
   setDefaults();
}

void StochOptions::setDefaults()
{
   /* default int values */
   setIntParam("dummy", 1);

   /* default double values */
   setDoubleParam("dummy", 1.0);
}
