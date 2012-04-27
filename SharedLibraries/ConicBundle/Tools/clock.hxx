/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Tools/clock.hxx

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



#ifndef CH_TOOLS__CLOCK_HXX
#define CH_TOOLS__CLOCK_HXX

/**  @file clock.hxx
    @brief Header declaring and (inline) implementing the classes CH_Tools::Microseconds and CH_Tools::Clock as well as some output routines for timing purposes
    @version 1.0
    @date 2006-06-19
    @author Christoph Helmberg

*/

#include <iostream>
#include <iomanip>
extern "C"{
#include <time.h>
}
#ifdef __unix
extern "C"{
#include <sys/resource.h>
//we need extern int getrusage(int who, struct rusage *r_usage);
}
#endif

/**@brief  Some convenient inline tools like a clock for timing, a random number generator, a heapsort template sorting indices 
*/	
namespace CH_Tools {

/**@defgroup Clock Clock (time measurement in Microseconds)
*/


  //@{

  /** @brief extra long integer number for expressing and computing time measurements in microseconds. 

  All operations assume that values and results of additions
  and subtractions are nonnegative, but this is not explicitly enforced.

  Internally Microseconds consits of two long values, seconds and microsecs,
  where microseconds is always <= 10^6 and any overflow is passed on to seconds. 
   */

  class Microseconds{
    bool infinity;     ///< if true, the value is regarded as infinity
    long seconds;      ///< counts time portion in seconds
    long microsecs;    ///< number in [0,10^6] counting the remaining microseconds without seconds

public:
  /** @name Constructors
   */
  //@{

    /// default constructor, value 0
    Microseconds(){seconds=0;microsecs=0;infinity=false;}
    /// constructor for setting value to 0 or infinity
    Microseconds(bool infty){seconds=0;microsecs=0;infinity=infty;}
    /// copy constructor
    Microseconds(const Microseconds& m)
    {seconds=m.seconds; microsecs=m.microsecs;infinity=m.infinity;}
    /// specify directly seconds and microseconds (in [0,10^6], no range check!)
    Microseconds(long secs, long msecs=0)
    {seconds=secs; microsecs=msecs;infinity=false;}
    
    /// specify directly seconds and microseconds (in [0,10^6], no range check!)
    Microseconds(int secs, int msecs=0)
    {seconds=secs; microsecs=msecs;infinity=false;}
    /// convert hours, minutes, secs, micros to Microseconds (no range check!)
    Microseconds(long hours,long minutes,long secs, long micros)
    {seconds=hours*3600+minutes*60+secs; microsecs=micros;infinity=false;}
    /// convert hours, minutes, secs, micros to Microseconds (no range check!)
    Microseconds(int hours,int minutes,int secs, int micros)
    {seconds=hours*3600+minutes*60+secs; microsecs=micros;infinity=false;}

    //@}

    
  /** @name Arithmetic operations and comparisons
   */
  //@{

    ///
    Microseconds& operator=(const Microseconds& m)
    {seconds=m.seconds; microsecs=m.microsecs; infinity=m.infinity;return *this;}
    
    ///
    Microseconds& operator+=(const Microseconds& m)
    {
      if (m.infinity){
	infinity=true;
	return *this;
      }
      microsecs+=m.microsecs;
      seconds+=m.seconds;
      while (microsecs>1000000){
	seconds++;
	microsecs-=1000000;
      }
      return *this;
    }
     
    ///
    Microseconds& operator-=(const Microseconds& m)
    {
     microsecs-=m.microsecs;
     seconds-=m.seconds;
     while (microsecs<0) {
         seconds--; microsecs+=1000000;
     }
     return *this;
    }

    /// 
    Microseconds operator-(const Microseconds& m) const
    {Microseconds s(*this); return s-=m;}
    ///
    Microseconds operator+(const Microseconds& m) const
    {Microseconds s(*this); return s+=m;}

    ///
    bool operator<(const Microseconds& m) const
    {
      if (infinity) return false;
      if (m.infinity) return true;
      return ((seconds<m.seconds)||((seconds==m.seconds)&&(microsecs<m.microsecs)));
    }

    ///
    bool operator>(const Microseconds& m) const
    {
      if (infinity&&(!m.infinity)) return true;
      if (m.infinity) return false;
      return ((seconds>m.seconds)||((seconds==m.seconds)&&(microsecs>m.microsecs)));}

    ///
    bool operator<=(const Microseconds& m) const
    {
      if (infinity&&(!m.infinity)) return false;
      if (m.infinity) return true;
      return ((seconds<m.seconds)||((seconds==m.seconds)&&(microsecs<=m.microsecs)));}

    ///
    bool operator>=(const Microseconds& m) const 
    {
      if (infinity) return true;
      if (m.infinity) return false;
      return ((seconds>m.seconds)||((seconds==m.seconds)&&(microsecs>=m.microsecs)));}

    ///
    bool operator==(const Microseconds& m) const
    {
     if ((infinity&&(!m.infinity))||(m.infinity)) return false;
     return ((seconds==m.seconds)&&(microsecs==m.microsecs));}

    //@}

  /** @name set, get, rounding, and conversions
   */
  //@{

    /// use true to regard value as infinity
    void set_infinity(bool infty){infinity=infty;}
    /// if true, value should be regarded as infinity
    bool get_infinity() const {return infinity;}

    /// convert to a double, where the size of one unit is one second(!)
    operator double() const 
    { return double(seconds)+double(microsecs)/1000000.; }

    /// convert and store the value of (*this) to hours, minutes, seconds 
    void hhmmss(long& hours,long& minutes,long& secs) const
    {
     long s=seconds;
     if (microsecs>=500000) s++;
     hours=s/3600;
     minutes=(s%3600)/60;
     secs=(s%60);
    }

    /// convert and store the value of (*this) to hours, minutes, seconds, hundredths 
    void hhmmssdd(long& hours,long& minutes,long& secs,long& hund) const
    {
     hund=(microsecs+5000)/10000;
     long s=seconds;
     if (hund==100) {s++;hund=0;}
     hours=s/3600;
     minutes=(s%3600)/60;
     secs=(s%60);
    }

    /// round the value to seconds
    long roundsecs() const
    {
     if (microsecs>=500000) return seconds+1;
     return seconds;
    }
    
    /// round the value to hundredths
    long roundhundredths() const
    {
     long hund=(microsecs+5000)/10000;
     hund+=seconds*100;
     return hund;
    }

    //@}

  /** @name Input and Output
   */
  //@{

    /// output in the format "seconds.microseconds" or "-1.000000" for infinity
    friend std::ostream& operator<<(std::ostream& out,const Microseconds& m)
    {if (m.infinity) return out<<"-1.000000";
     out<<m.seconds<<".";out.fill('0');out.width(6);out<<m.microsecs;
     out.fill(' ');return out;}


    /// input in the format "seconds.microseconds" , seconds<0 is regarded as infinity
    friend std::istream& operator>>(std::istream& in,Microseconds& m)
    {
     char c; m.infinity=false;
     in>>m.seconds>>c>>m.microsecs;
     if (m.seconds<0) {m.seconds=0; m.infinity=true;}
     return in;
    }     

  ///print Microseconds in the format "hh:mm:ss.dd" or "hh:mm:ss" 
    friend void print_time(std::ostream& out,     ///< output stream
			   const Microseconds& m, ///< time 
			   int secondsonly=0      ///< use =0 for "hh:mm:ss.dd", !=0 for "hh:mm:ss"
			   )
    {
      long hours,minutes,seconds,hunds;
      out.fill('0');
      if (secondsonly){
	m.hhmmss(hours,minutes,seconds);
	out.width(2);out<<hours;
	out<<":";out.width(2);out<<minutes;
	out<<":";out.width(2);out<<seconds;
      }
      else {
	m.hhmmssdd(hours,minutes,seconds,hunds);
	out.width(2);out<<hours;
	out<<":";out.width(2);out<<minutes;
	out<<":";out.width(2);out<<seconds;
	out<<".";out.width(2);out<<hunds;
      }
      out.fill(' ');
    }
    
    //@}

};



  /** @brief allows measuring time difference to its initialization time in Microseconds 

      the clock is initiliazed by calling start() (this also happens automatically 
      at construction) and returns, at each call to time(), the time passed
      since the last call to start() in Microseconds (optionally plus a given
      offset specified by set_offset() ). The routine elapsed_time() prints
      the current time difference together with current date and time 
      to an output strean. 
   */

class Clock
{
  private:
    Microseconds t_start;
    Microseconds offset;
#ifdef __unix
#else
    long l_start;
#endif

  public:
  /// read current time, all further time measurements will be in relation to this time
    void start(){
#ifdef __unix
        struct rusage r_usage;
        getrusage(RUSAGE_SELF,&r_usage);
        t_start=Microseconds(r_usage.ru_utime.tv_sec,r_usage.ru_utime.tv_usec);
        offset=Microseconds(0,0);
#else
        t_start=Microseconds(static_cast<int>(::time(0)),0);
#endif
    }

  /// calls start()
    Clock(){start();}
  /// nothing particular
    ~Clock(){}

  /// allows to specify an offset, that will furtheron be added to all time measurements 
    void set_offset(Microseconds offs){
        offset=offs;
    }

  /// return time elapsed since last call to start() in Microseconds (possibly adding an optional offset)
    Microseconds time() const{
#ifdef __unix
        struct rusage r_usage;
        getrusage(RUSAGE_SELF,&r_usage);
        Microseconds elapsed(r_usage.ru_utime.tv_sec,r_usage.ru_utime.tv_usec);
#else
        Microseconds elapsed(static_cast<int>(::time(0)),0);
#endif
        elapsed-=t_start;
        elapsed+=offset;
        return elapsed;
    }

  /// call time() and print the result in format "hh:mm:ss", togehter with current date and time, to out
    void elapsed_time(std::ostream& out) const
    {
        time_t current=::time(0);
        struct tm *loctime;                     
        Microseconds spent=time();
        out<<"elapsed time: ";
        print_time(out,spent);
        loctime=localtime(&current);
        out<<"   ----   "<<asctime(loctime)<<std::flush;
    }

  /// call cl.time() and print the result in format "hh:mm:ss" to out
    friend inline std::ostream& operator<<(std::ostream& out,const Clock& cl)
    {
        print_time(out,cl.time());
        return out;
    }
};

//@}

}

#endif

