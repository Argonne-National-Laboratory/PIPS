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
      setDefaults();
   }

   void QpGenOptions::setDefaults()
   {
      /// OUTER BIGCSTAB
      double_options["OUTER_BICG_EPSILON"] = 1e-15;

      bool_options["OUTER_BICG_PRINT_STATISTICS"] = false;

      int_options["OUTER_BICG_MAX_ITER"] = 75;
      int_options["OUTER_BICG_MAX_NORMR_DIVERGENCES"] = 4;
      int_options["OUTER_BICG_MAX_STAGNATIONS"] = 4;
   }
}
