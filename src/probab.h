/*
 * probab.h
 *
 *  Created on: May 15, 2011
 *      Author: balwierz
 */
#include "TSC.h"
#include <list>

#ifndef PROBAB_H_
#define PROBAB_H_

double altprofprob(list<TSC*>::iterator it);
double profProbAllSamples(list<TSC*>::iterator it);
double priorprob(list<TSC*>::iterator it);
double priorprobEmpirical(TSC *segone, TSC *segtwo);

#endif /* PROBAB_H_ */
