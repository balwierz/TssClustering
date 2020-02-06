/*
 * probab.cpp
 *
 *  Created on: May 15, 2011
 *      Author: balwierz
 */

#include "probab.h"
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>
#include <iostream>

extern int nSamples;
double sigmasq = 0.0169;  // noise size
double alpha = 3.1;    // gaussian prior size for TSS variance
double lenscale = 10; // coexpression prior scale depending on the distance
double normPseudocount = 0.5; // in tpms
double rawPseudocount = 0.5; // in tags
// the two above are not related with the scaling of depth.

boost::unordered_set<sampleIdT> samplesUnion;

// Calculates the likelihood log ratio:
// log( P(D|coexpressed) / P(D|independent )

double altprofprob(list<TSC*>::iterator it)
// we might remove the second one in case we only want to cluster i and i+1...
{
	TSC* segone = *it;
	it++;
	TSC* segtwo = *it;
	double L = 0;
	double Lindep = 0;

	double avwzz = 0;
	double avwz = 0;
	double avgamma = 0;
	double avw = 0;
	double avgammarho = 0;
	double avgammarhorho = 0;
	double avgammasigsig = 0;
	double avgammasig = 0;
	double avgammarhosig = 0;

	double avbetax = 0;
	double avbetay = 0;
	double avbetax_xx = 0;
	double avbetay_yy = 0;
	double avbetax_x = 0;
	double avbetay_y = 0;

	// go over all the samples in which we have expression of the first TSC
	// and check if we have also expression in the second
	// if true, then do the calculations.
	int pos = 0; // keeps the number of used samples
	//cout << "comparing " << segone->repr << " and " << segtwo->repr << endl;
	BOOST_FOREACH(const exprHashT::value_type &pair, segone->expr)
	{
		//string sample = pair.first;
		//cout << sample << '\t';
		if(segtwo->expr.find(pair.first) != segtwo->expr.end())
		{
			//cout << "M";
			double n = segone->expr[pair.first].raw;
			double m = segtwo->expr[pair.first].raw;
			double nnor = segone->expr[pair.first].norm;
			double mnor = segtwo->expr[pair.first].norm;
			// here can be some pseudocounts
			pos++;
			double x = log(nnor + normPseudocount);
			double y = log(mnor + normPseudocount);
			double z = x - y;
			double s = 0.5 * (x + y);
			double wx = 1 / (sigmasq + 1.0 / n);
			double wy = 1 / (sigmasq + 1.0 / m);
			double w = wx * wy * (1.0 + alpha / (wx + wy)) / (wx + wy + alpha);
			double rho = 0.5 * (wx - wy) / (wx + wy);
			double sig = s + rho * z;
			double gamma = alpha * (wx + wy) / (wx + wy + alpha);

			avwzz += w * z * z;
			avwz += w * z;
			avgammasigsig += gamma * sig * sig;
			avgammasig += gamma * sig;
			avgamma += gamma;
			avw += w;
			avgammarho += gamma * rho;
			avgammarhorho += gamma * rho * rho;
			avgammarhosig += gamma * rho * sig;

			L += 0.5 * log(wx * wy * alpha / (wx + wy + alpha));

			double betax = alpha * wx / (alpha + wx);
			double betay = alpha * wy / (alpha + wy);
			avbetax += betax;
			avbetay += betay;
			avbetax_xx += betax * x * x;
			avbetay_yy += betay * y * y;
			avbetax_x += betax * x;
			avbetay_y += betay * y;

			Lindep += 0.5 * log(betax * betay);
		}
	}
	/**normalize all the averages****/
	//cout << endl;
	if(pos > 0)
	{
		avwzz /= ((double) pos);
		avw /= ((double) pos);
		avwz /= ((double) pos);
		avgamma /= ((double) pos);
		avgammasig /= ((double) pos);
		avgammasigsig /= ((double) pos);
		avgammarho /= ((double) pos);
		avgammarhorho /= ((double) pos);
		avgammarhosig /= ((double) pos);

		avbetax /= ((double) pos);
		avbetay /= ((double) pos);
		avbetax_x /= ((double) pos);
		avbetay_y /= ((double) pos);
		avbetax_xx /= ((double) pos);
		avbetay_yy /= ((double) pos);

		/***contribution from the prefactor***/
		L -= 0.5 * log( avgamma * avw + avgamma * avgammarhorho - avgammarho * avgammarho);
		/**Note the 1/Nsamp, the 2pi's,  and the 1/(mumax-mumin) will cancel between dependent and independent model***/

		/***Expression in the exponent of the dependent model***/
		double Q = avwzz + avgammasigsig - avgammasig * avgammasig / avgamma;
		Q -= (avwz + avgammarhosig - avgammarho * avgammasig / avgamma) * (avwz + avgammarhosig - avgammarho * avgammasig / avgamma) / (avw + avgammarhorho - avgammarho * avgammarho / avgamma);

		L -= 0.5 * ((double) pos) *  Q;

		/***Now contribution from independent model***/
		Lindep -= 0.5 * ((double) pos) * (avbetax_xx - (avbetax_x * avbetax_x / avbetax));
		Lindep -= 0.5 * ((double) pos) * (avbetay_yy - (avbetay_y * avbetay_y / avbetay));
		Lindep -= 0.5 * log( avbetax * avbetay);

		/***Difference between the log-likelihoods***/
		L -= Lindep;
	}
	else
	{
		L = 0;
	}
	//cerr << "Returning " << L << endl;
	if(L != L) //(L == numeric_limits<double>::infinity())
	{
		cerr << "oops: we produced inf or NaN in altprofprob\n";
		cerr << "pos = " << pos << " Lindep = " << Lindep << " x = " << avbetax << " y " << avbetay
				<< " x_xx " << avbetax_xx << " y_yy " << avbetay_yy << " x_x " << avbetax_x << " y_y" << avbetay_y
				<< segone << '\n' << segtwo << endl;

		exit(1);
	}
	return L;
}


// this one is as the one above, but includes all the samples...

double profProbAllSamples(list<TSC*>::iterator it)
// we might remove the second one in case we only want to cluster i and i+1...
{
	TSC* segone = *it;
	it++;
	TSC* segtwo = *it;
	double L = 0;
	double Lindep = 0;

	double avwzz = 0;
	double avwz = 0;
	double avgamma = 0;
	double avw = 0;
	double avgammarho = 0;
	double avgammarhorho = 0;
	double avgammasigsig = 0;
	double avgammasig = 0;
	double avgammarhosig = 0;

	double avbetax = 0;
	double avbetay = 0;
	double avbetax_xx = 0;
	double avbetay_yy = 0;
	double avbetax_x = 0;
	double avbetay_y = 0;

	// go over all the samples in which we have expression of the first TSC
	// and check if we have also expression in the second
	// if true, then do the calculations.
	int pos = 0; // keeps the number of used samples

	// find out an union of the samples
	samplesUnion.clear();
	BOOST_FOREACH(const exprHashT::value_type &pair, segone->expr)
	{
		samplesUnion.emplace(pair.first);
	}
	BOOST_FOREACH(const exprHashT::value_type &pair, segtwo->expr)
	{
		samplesUnion.emplace(pair.first);
	}

	//
	BOOST_FOREACH(const sampleIdT sample, samplesUnion)
	{
		double n = rawPseudocount;
		double m = rawPseudocount;
		double nnor = 0.0;
		double mnor = 0.0;
		if(segone->expr.find(sample) != segone->expr.end())
		{
			n = segone->expr[sample].raw;
			nnor = segone->expr[sample].norm;
		}
		if(segtwo->expr.find(sample) != segtwo->expr.end())
		{
			m = segtwo->expr[sample].raw;
			mnor = segtwo->expr[sample].norm;
		}
		// here can be some pseudocounts
		pos++;
		double x = log(nnor + normPseudocount);
		double y = log(mnor + normPseudocount);
		double z = x - y;
		double s = 0.5 * (x + y);
		double wx = 1 / (sigmasq + 1.0 / n);
		double wy = 1 / (sigmasq + 1.0 / m);
		double w = wx * wy * (1.0 + alpha / (wx + wy)) / (wx + wy
				+ alpha);
		double rho = 0.5 * (wx - wy) / (wx + wy);
		double sig = s + rho * z;
		double gamma = alpha * (wx + wy) / (wx + wy + alpha);

		avwzz += w * z * z;
		avwz += w * z;
		avgammasigsig += gamma * sig * sig;
		avgammasig += gamma * sig;
		avgamma += gamma;
		avw += w;
		avgammarho += gamma * rho;
		avgammarhorho += gamma * rho * rho;
		avgammarhosig += gamma * rho * sig;

		L += 0.5 * log(wx * wy * alpha / (wx + wy + alpha));

		double betax = alpha * wx / (alpha + wx);
		double betay = alpha * wy / (alpha + wy);
		avbetax += betax;
		avbetay += betay;
		avbetax_xx += betax * x * x;
		avbetay_yy += betay * y * y;
		avbetax_x += betax * x;
		avbetay_y += betay * y;

		Lindep += 0.5 * log(betax * betay);
	}
	//**normalize all the averages***
	if(pos > 0)
	{
		avwzz /= ((double) pos);
		avw /= ((double) pos);
		avwz /= ((double) pos);
		avgamma /= ((double) pos);
		avgammasig /= ((double) pos);
		avgammasigsig /= ((double) pos);
		avgammarho /= ((double) pos);
		avgammarhorho /= ((double) pos);
		avgammarhosig /= ((double) pos);

		avbetax /= ((double) pos);
		avbetay /= ((double) pos);
		avbetax_x /= ((double) pos);
		avbetay_y /= ((double) pos);
		avbetax_xx /= ((double) pos);
		avbetay_yy /= ((double) pos);

		/***contribution from the prefactor***/
		L -= 0.5 * log( avgamma * avw + avgamma * avgammarhorho - avgammarho * avgammarho);
		/**Note the 1/Nsamp, the 2pi's,  and the 1/(mumax-mumin) will cancel between dependent and independent model***/

		/***Expression in the exponent of the dependent model***/
		double Q = avwzz + avgammasigsig - avgammasig * avgammasig / avgamma;
		Q -= (avwz + avgammarhosig - avgammarho * avgammasig / avgamma) * (avwz + avgammarhosig - avgammarho * avgammasig / avgamma) / (avw + avgammarhorho - avgammarho * avgammarho / avgamma);

		L -= 0.5 * ((double) pos) *  Q;

		/***Now contribution from independent model***/
		Lindep -= 0.5 * ((double) pos) * (avbetax_xx - (avbetax_x * avbetax_x / avbetax));
		Lindep -= 0.5 * ((double) pos) * (avbetay_yy - (avbetay_y * avbetay_y / avbetay));
		Lindep -= 0.5 * log( avbetax * avbetay);

		/***Difference between the log-likelihoods***/
		L -= Lindep;
	}
	else
	{
		cerr << "No samples used in likelihood calculation. Something went wrong." << endl;
		exit(1);
		L = 0;
	}
	//cerr << "Returning " << L << endl;
	if(L != L) //(L == numeric_limits<double>::infinity())
	{
		cerr << "oops: we produced inf or NaN in profProbAllSamples()\n";
		cerr << "pos = " << pos << " Lindep = " << Lindep << " x = " << avbetax << " y = " << avbetay
				<< " x_xx = " << avbetax_xx << " y_yy = " << avbetay_yy << " x_x = " << avbetax_x << " y_y = " << avbetay_y << '\n'
				<< segone << segtwo << endl;

		exit(1);
	}
	return L;
}



 /*

	for(i = 0; i < nSamples; ++i) // over samples
	{
		string sample = samples[i];
		n = segone->expr[sample].norm;
		m = segtwo->expr[sample].norm;

		if(n >= 1 || m >= 1 )
		{
			++pos;

			nnor = segone->normprofile[i];
			mnor = segtwo->normprofile[i];
			if(n < 1)
			{
				n = 1;
				nnor += 0.5 * 1000000.0 / (tottag[i] + 1.0);
			}
			else
			{
				++posx;
			}

			if(m < 1)
			{
				m = 1;
				mnor += 0.5 * 1000000.0 / (tottag[i] + 1.0);
			}
			else
			{
				++posy;
			}

			x = log(nnor);
			y = log(mnor);
			z = x - y;
			s = 0.5 * (x + y);
			wx = 1 / (sigmasq[i] + 1.0 / n);
			wy = 1 / (sigmasq[i] + 1.0 / m);
			w = wx * wy * (1.0 + alpha / (wx + wy)) / (wx + wy + alpha);
			rho = 0.5 * (wx - wy) / (wx + wy);
			sig = s + rho * z;
			gamma = alpha * (wx + wy) / (wx + wy + alpha);

			avwzz += w * z * z;
			avwz += w * z;
			avgammasigsig += gamma * sig * sig;
			avgammasig += gamma * sig;
			avgamma += gamma;
			avw += w;
			avgammarho += gamma * rho;
			avgammarhorho += gamma * rho * rho;
			avgammarhosig += gamma * rho * sig;

			L += 0.5 * log( wx * wy * alpha / (wx + wy + alpha));

			betax =  alpha * wx / (alpha + wx);
			betay =  alpha * wy / (alpha + wy);
			avbetax += betax;
			avbetay += betay;
			avbetax_xx += betax * x * x;
			avbetay_yy += betay * y * y;
			avbetax_x += betax * x;
			avbetay_y += betay * y;

			Lindep += 0.5 * log(betax * betay);
		}
	}
	//normalize all the averages
	if(pos > 0)
	{
		avwzz /= ((double) pos);
		avw /= ((double) pos);
		avwz /= ((double) pos);
		avgamma /= ((double) pos);
		avgammasig /= ((double) pos);
		avgammasigsig /= ((double) pos);
		avgammarho /= ((double) pos);
		avgammarhorho /= ((double) pos);
		avgammarhosig /= ((double) pos);

		avbetax /= ((double) pos);
		avbetay /= ((double) pos);
		avbetax_x /= ((double) pos);
		avbetay_y /= ((double) pos);
		avbetax_xx /= ((double) pos);
		avbetay_yy /= ((double) pos);

		//contribution from the prefactor
		L -= 0.5 * log( avgamma * avw + avgamma * avgammarhorho - avgammarho * avgammarho);
		//Note the 1/Nsamp, the 2pi's,  and the 1/(mumax-mumin) will cancel between dependent and independent model

		//Expression in the exponent of the dependent model
		Q = avwzz + avgammasigsig - avgammasig * avgammasig / avgamma;
		Q -= (avwz + avgammarhosig - avgammarho * avgammasig / avgamma) * (avwz + avgammarhosig - avgammarho * avgammasig / avgamma) / (avw + avgammarhorho - avgammarho * avgammarho / avgamma);

		L -= 0.5 * ((double) pos) *  Q;

		//Now contribution from independent model
		Lindep -= 0.5 * ((double) pos) * (avbetax_xx - (avbetax_x * avbetax_x / avbetax));
		Lindep -= 0.5 * ((double) pos) * (avbetay_yy - (avbetay_y * avbetay_y / avbetay));
		Lindep -= 0.5 * log( avbetax * avbetay);

		//Difference between the log-likelihoods
		L -= Lindep;
	}
	else
	{
		L = 0;
	}
	return L;
}
*/

/**prior probability for fusing segments depending on their distance***/
double priorprob(list<TSC*>::iterator it)
{
	TSC* segone = *it;
	it++;
	TSC* segtwo = *it;
	if ((void*)segone == (void*)segtwo)
	{
		cerr << "Comparing same objects!" << endl;
	}
	// the old way of calculating the distance:
	//int dist = segtwo->beg - segone->end ;

	// the new way of calculating the distance:
	int dist = segtwo->repr - segone->repr;

	//cout << "Distance " << dist << endl;
	/***test use instead distances between representative positions***/
	/*dist = segtwo->reppos-segone->reppos+1;*/
	double prior = -((double) dist) / 30.0;
	prior -= log(1.0 - exp(prior));
	if(isnan(prior))
	{
		cerr << "NAN prior\n" << endl;
		exit(1);
	}
	if(isinf(prior))
		{
			cerr << "INF prior, while dist = " << dist << endl;
			cerr << segone->beg << " " << segone->end << " " << segtwo->beg << " " << segtwo->end << endl;
			exit(1);
		}
	return prior;
}
// these formulas work
// 304*0.48 * exp(-x/30)  // fits the positive part roughly
// 304*0.75* exp(-x/70) - 0.3*304  // fits quite good apart from the 1

double priorprobEmpirical(TSC *segone, TSC *segtwo)
{
	int dist = segtwo->beg - segone->end ;
	if(dist < 1)
	{
		cerr << "Distance between TSCs is too small (=" << dist << ")\n";
		exit(1);
	}
	if(dist == 1)
	{
		return nSamples * 0.63;
	}
	else
	{
		return nSamples * 0.75 * exp(-(double)dist / 70.0) - 0.3 * nSamples;
	}
}
