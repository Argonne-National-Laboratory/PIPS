/*
 * QpGenOptions.h
 *
 *  Created on: 01.07.2020
 *      Author: bzfkempk
 */

#include "QpGenOptions.h"
#include <limits>

namespace qpgen_options
{
   QpGenOptions::QpGenOptions()
   {
      QpGenOptions::setDefaults();
   }

   void QpGenOptions::setDefaults()
   {
      // TODO
      /* default bool values */
      bool_options["dummy"] = false;
      bool_options["OUTER_BICG_PRINT_STATISTICS"] = false;

      /* default int values */
      int_options["dummy"] = 1;
      int_options["OUTER_BICG_MAX_ITER"] = 75;
      int_options["OUTER_BICG_MAX_NORMR_DIVERGENCES"] = 4;
      int_options["OUTER_BICG_MAX_STAGNATIONS"] = 4;

      /* default double values */
      double_options["dummy"] = 1.0;
      double_options["OUTER_BICG_EPSILON"] = 1e-15;
   }
}
