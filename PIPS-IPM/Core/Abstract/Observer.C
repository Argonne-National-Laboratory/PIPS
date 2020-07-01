/*
 * Observer.h
 *
 *  Created on: 30.06.2020
 *      Author: bzfkempk
 */
#include "Observer.h"

/// SUBJECT
void Subject::registerObserver(Observer* observer)
{
   observers.push_back(observer);
}

void Subject::unregisterObserver(Observer* observer)
{
   observers.remove(observer);
}

/// OBSERVER
const Subject* Observer::getSubject() const
{
   return subj;
}

Observer::~Observer()
{
   removeSubject();
}

void Observer::setSubject(Subject* subject)
{
   if( subj != nullptr )
      subj->unregisterObserver(this);
   subj = subject;

   subject->registerObserver(this);
};

void Observer::removeSubject()
{
   if( subj != nullptr )
      subj->unregisterObserver(this);
   subj = nullptr;
}

void Subject::notifyObservers()
{
   for( Observer* o : observers )
      o->notifyFromSubject();
}
