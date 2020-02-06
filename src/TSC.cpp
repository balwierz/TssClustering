/*
 * TSC.cpp
 *
 *  Created on: May 14, 2011
 *      Author: balwierz
 */

#include "TSC.h"
#include <list>
#include <fstream>


using namespace std;

uint32_t TSC::namingCounter = 0;
string TSC::chromosome = string();
string TSC::strand = string();

void TSC::operator +=(const TSC& other)
{
	//cerr << "fusing (" << beg << "," << end << ") (" << other.beg << "," << other.end << ")\n";
	end = other.end;
	// find out the representative position:
	if(reprTot < other.reprTot)
	{
		reprTot = other.reprTot;
		repr = other.repr;
	}
	BOOST_FOREACH(const exprHashT::value_type &val, other.expr)
	{
		expr[val.first] += val.second; // tested: it works even if the key didn't exist. unkn as zero.
	}
}
//chr7:934,303-1,016,950
// prints in the BED format
ostream& operator<<(ostream &str, TSC* ptr)
{
	//str << "I" << endl;
	str << TSC::chromosome << '\t'						// 1
			<< (ptr->beg - 1) << "\t"
			<< ptr->end << '\t'
			<< ptr->name() << '\t'
			<< (ptr->end - ptr->beg + 1) << "\t" 		// 5
			<< TSC::strand << "\t"
			<< ptr->totalExpression() << '\t'
			<< ptr->maxExpression().norm << "\t"
			<< ptr->maxExpression().raw << '\t'
			<< ptr->surrDensity << '\t' 				// 10
			<< ptr->surrDensityOfMax << '\t'
			<< ptr->probLogRatio << '\t'
			<< ptr->surrDensity2 << '\n';				// 13

	return str;
}


void printClusters(list<TSC*> &result, const char* fileName, double maxTpmCutoff)
{
	ofstream out(fileName);
	for(list<TSC*>::iterator it = result.begin(); it != result.end(); it++)
	{
		if((*it)->maxExpression().norm > maxTpmCutoff)
			out << (*it);
	}
	out.close();
}

void printClustersBed(list<TSC*> &result, const char* fileName, double maxTpmCutoff)
{
	ofstream out(fileName);
	for(list<TSC*>::iterator it = result.begin(); it != result.end(); it++)
	{
		if((*it)->maxExpression().norm > maxTpmCutoff)
		{
			string name = string() + "TSC_v2.02_" + TSC::chromosome + "_" + TSC::strand + "_" + boost::lexical_cast<string>((*it)->beg);
			out << TSC::chromosome << '\t'
				<< ((*it)->beg - 1) << '\t'  //UCSC bed format
				<< (*it)->end << '\t'
				<< name << '\t'
				<< (*it)->totalExpression() << '\t'
				<< TSC::strand << '\n';
		}
	}
	out.close();
}
