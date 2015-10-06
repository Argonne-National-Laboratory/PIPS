/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef TWODIMFILTER_H
#define TWODIMFILTER_H

class FilterIPMOption;

class Filter2D  {
  public:
  double filter_ConNorm;
  double filter_Obj;
  Filter2D *nextfilter;
  
  Filter2D();
  Filter2D(double filterTheta_in, double filterPhi_in);
  ~Filter2D();

  //initialize the filter
  void Initialize(FilterIPMOption* FilterIPMOpt);

  //check if the given point is in the filter
  bool WithinFilter(double const curr_ConNorm, double const curr_Obj);
  
  //update the filter by given values
  void UpdateFilter(double const curr_ConNorm, double const curr_Obj);

} ;


#endif
