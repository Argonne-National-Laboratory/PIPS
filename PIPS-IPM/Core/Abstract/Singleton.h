/*
 * Options.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_SINGLETON_H_
#define PIPS_IPM_CORE_ABSTRACT_SINGLETON_H_

/**
 * Abstract base class for options class.
 * Implements a inheritable to a singleton pattern.
 * The getInstance method must be specified in the BaseClasses,
 * thus for each base class there can only be exactly one instance.
 */




class Singleton
{
   /// delete copy constructor and assignment operator for all inheriting classes
   #ifdef PRE_CPP11
      Singleton(Singleton const&);
      void operator=(Singleton const&);
   #else
      Singleton(Singleton const&) = delete;
      void operator=(Singleton const&) = delete;
   #endif

   protected:
      Singleton() {};
      virtual ~Singleton() {};
};

#endif /* PIPS_IPM_CORE_ABSTRACT_SINGLETON_H_ */
