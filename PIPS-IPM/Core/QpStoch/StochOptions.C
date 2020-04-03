#include "StochOptions.h"

StochOptions::StochOptions()
{
   setDefaults();
}

void StochOptions::setDefaults()
{
   // TODO
   /* default bool values */
   bool_options["dummy"] = false;
   bool_options["postsolve"] = true;

   /* default int values */
   int_options["dummy"] = 1;

   /* default double values */
   double_options["dummy"] = 1.0;
}
