/*
 * Observer.h
 *
 *  Created on: 30.06.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_OBSERVER_H_
#define PIPS_IPM_CORE_ABSTRACT_OBSERVER_H_

#include <list>
/**
 * Abstract base classes for status based observer pattern - contains Subject and Observer
 * state is part of the subject and whenever it changes all observers get notified
 */

template<class T>
class Subject;

template<class T>
class Observer
{
   private:
      Subject<T> *subj = nullptr;

   public:
      virtual ~Observer() {};

      bool hasSubject() const
      {
         return subj != nullptr;
      };

      void setSubject(Subject<T> *subject);
      void removeSubject();

      virtual void update() = 0;

};

template<class T>
class Subject
{
   private:
      std::list<Observer<T>*> observers;
      T state;

   public:
      void registerObserver(Observer<T>* observer)
      {
         observers.push_back(observer);
      }

      void unregisterObserver(Observer<T>* observer)
      {
         observers.remove(observer);
      }

      const T& getState()
      {
         return state;
      }

      void setState(const T& s);
};

template<class T>
void Observer<T>::setSubject(Subject<T>* subject)
{
   if( subj != nullptr )
      subj->unregisterObserver(this);
   subj = subject;

   subject->registerObserver(this);
};

template<class T>
void Observer<T>::removeSubject()
{
   if( subj != nullptr )
      subj->unregisterObserver(this);
}

template<class T>
void Subject<T>::setState(const T& s)
{
      if( s != state )
      {
         state = s;

         for( Observer<T>* o : observers )
            o->update();
      }
}

#endif /* PIPS_IPM_CORE_ABSTRACT_OBSERVER_H_ */
