/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Tools/gb_rand.hxx

    This file is part of ConciBundle, a C/C++ library for convex optimization.

    ConicBundle is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ConicBundle is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

***************************************************************************** */



#ifndef CH_TOOLS__GB_RAND_HXX
#define CH_TOOLS__GB_RAND_HXX

#include <iostream>

/**  @file gb_rand.hxx
    @brief  Header declaring and (inline) implementing the random number generator CH_Tools::GB_rand
    @version 1.0
    @date 2005-03-01
    @author Christoph Helmberg

    The code is extracted/translated/adapted from an indirect 
    copy of a socalled Stanford graph-base library
*/

namespace CH_Tools {

/**@defgroup GB_rand Random Number Generator
*/
  //@{

  /** @brief device independent random number generator based on long int with seed

  */

class GB_rand
{
private:
  long A[56];   ///< storage for the next few precomputed random numbers 
  int ind;      ///< index of the next number to return

  /// compute (x-y) & 0x7fffffff
  long mod_diff(long x,long y){return ((x-y)&0x7fffffff);}

  /// compute next set of numbers
  long flip_cycle ()
  {
    register long *ii, *jj;
    for (ii = &A[1], jj = &A[32]; jj <= &A[55]; ii++, jj++)
      *ii = ((*ii-*jj) & 0x7fffffff);
    for (jj = &A[1]; ii <= &A[55]; ii++, jj++)
      *ii = ((*ii-*jj) & 0x7fffffff);
    ind = 54;
    return A[55];
  }
    
public:
  /// restart generator with seed
  void init(long seed=1)
  {
    register long i;
    register long prev = seed, next = 1;
    seed = prev = (prev & 0x7fffffff);
    A[55] = prev;
    for (i = 21; i; i = (i + 21) % 55) {
      A[i] = next;
      
      next = ((prev-next) & 0x7fffffff);
      if (seed & 1)
	seed = 0x40000000 + (seed >> 1);
      else
	seed >>= 1;
      next = ((next- seed) & 0x7fffffff);
      
      prev = A[i];
    }
    
    (void) flip_cycle ();
    (void) flip_cycle ();
    (void) flip_cycle ();
    (void) flip_cycle ();
    (void) flip_cycle ();
    
  }
  
  /// calls  init(seed)
  GB_rand(long seed=1)
  {ind=0;A[0]=-1;init(seed);}

  ///
  ~GB_rand(){}
    
  /// returns a random integer number "uniformly distributed" in {0,..,m-1} 
  long unif_long(long m)
  {
    const unsigned long two_to_the_31 = (unsigned long)0x80000000;
    register unsigned long t = two_to_the_31 - (two_to_the_31 % m);
    register long r;
    do {
      r = (A[ind]>=0?(A[ind--]):flip_cycle());
    } while (t <= (unsigned long) r);
    return r % m;
  }

  /// returns a random double number "uniformly distributed" in (0,1)
  double next()
  {
    long r=unif_long(0x40000000);
    return (double(r)+.5)/double(0x40000000);
  }
  
  /// save current configuration to out so as to continue identically after restore
  std::ostream& save(std::ostream& out) const
  {
    for(int i=0;i<56;i++) out<<A[i]<<"\n";
    out<<ind<<"\n"; return out;
  }

  /// restore settings as stored in save 
  std::istream& restore(std::istream& in)
  {
    for(int i=0;i<56;i++)
      in>>A[i];
    in>>ind; return in;
  }

};

  //@}

}

#endif

