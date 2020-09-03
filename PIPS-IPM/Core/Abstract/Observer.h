/*
 * Observer.h
 *
 *  Created on: 30.06.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_OBSERVER_H_
#define PIPS_IPM_CORE_ABSTRACT_OBSERVER_H_

#include <list>
#include <string>
#include "pipsport.h"
/**
 * Abstract base classes for status based observer pattern - contains Subject and Observer
 * Observer can query Subject for ints, doubles, and bools defined via a std::string
 */

class Subject;

class Observer
{
   private:
      Subject *subj/* = nullptr*/;

   public:
      Observer();
      virtual ~Observer();

      bool hasSubject() const
      {
         return subj != nullptr;
      };

      void setSubject(Subject *subject);
      const Subject* getSubject() const;
      void removeSubject();

      virtual void notifyFromSubject() = 0;
};

class Subject
{
   private:
      std::list<Observer*> observers;

   public:
      virtual ~Subject() {};

      void registerObserver(Observer* observer);
      void unregisterObserver(Observer* observer);

      virtual int getIntValue(const std::string& s) const = 0;
      virtual double getDoubleValue(const std::string& s) const = 0;
      virtual bool getBoolValue(const std::string& s) const = 0;

      void notifyObservers();
};

#endif /* PIPS_IPM_CORE_ABSTRACT_OBSERVER_H_ */
