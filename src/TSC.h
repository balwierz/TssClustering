/*
 * TSC.h
 *
 *  Created on: May 14, 2011
 *      Author: balwierz
 */

#ifndef TSC_H_
#define TSC_H_

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>
#include <ostream>
#include <boost/lexical_cast.hpp> // for ,,atoi''
#include <string>
#include <fstream>
#include <list>

using namespace std;

struct expression
{
	double raw;  // could be int, but then we have a conversion each time...
	double norm;
	expression(double r, double n) : raw(r), norm(n) {}
	expression(const expression& org) // maybe remove for speed ups
	{
//		/cerr << "Copying\n";
		raw = org.raw;
		norm = org.norm;
	}
	expression(): raw(0), norm(0)
	{
	}
	inline void operator += (const expression &other)
	{
		raw += other.raw;
		norm += other.norm;
	}
};
typedef int sampleIdT;
//typedef string sampleIdT;  //file names as keys
typedef boost::unordered_map<sampleIdT, expression> exprHashT; //key is a sample
struct TSC
{
	static uint32_t namingCounter;
	static string chromosome;
	static string strand;
	int repr;
	int beg;
	int end;
	double reprTot; // total expression of a representative position
	double probLogRatio;
	double surrDensity;
	double surrDensityOfMax;
	double surrDensity2;
	bool dropLeft;
	bool dropRight;
	int posInHeap;
	exprHashT expr; // hey is a sample
	TSC() : reprTot(0.0), probLogRatio(0.0), surrDensity(0.0), surrDensityOfMax(0.0), surrDensity2(0.0), dropLeft(false), dropRight(false)
	{
		//cerr << "constr:TSC\n";
	}
	void operator+=(const TSC&);
	string name()
	{
		string ret = "TSC_" + boost::lexical_cast<string>(namingCounter++);
		return ret;
	}
	double totalExpression() const
	{
		double ret = 0.0;
		BOOST_FOREACH(const exprHashT::value_type &sampleIexpr, expr )
		{
			ret += sampleIexpr.second.norm;
		}
		return ret;
	}
	uint32_t totalRaw() const // gives the RAW expression summed over samples
	{
		uint32_t ret = 0;
		BOOST_FOREACH(const exprHashT::value_type &sampleIexpr, expr )
		{
			ret += (uint32_t)sampleIexpr.second.raw; //this is double to int conversion!
		}
		return ret;
	}
//	uint32_t maxRaw() const // gives the RAW expression summed over samples
//	{
//		uint32_t ret = 0;
//		BOOST_FOREACH(const exprHashT::value_type &sampleIexpr, expr )
//		{
//			if(ret < sampleIexpr.second.raw)
//			ret = sampleIexpr.second.raw;
//		}
//		return ret;
//	}
	expression maxExpression() const
	{
		expression ret;
		BOOST_FOREACH(const exprHashT::value_type &sampleIexpr, expr )
		{
			if(ret.norm < sampleIexpr.second.norm)
				ret = sampleIexpr.second;
		}
		return ret;
	}
	int expressedInNSamples() const
	{
		return expr.size();
	}

};



ostream& operator<<(ostream &str, TSC* ptr);
void printClusters(list<TSC*> &result, const char* fileName, double  maxTpmCutoff = 0.0);
void printClustersBed(list<TSC*> &result, const char* fileName, double  maxTpmCutoff = 0.0);

#endif /* TSC_H_ */
