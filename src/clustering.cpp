//============================================================================
// Name        : clustering.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : it shoule be run for each chr-strand combination separately
//					for speedups
//============================================================================
#include <fstream>
#include <iostream>
#include <list>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/unordered_map.hpp>
#include <boost/array.hpp>
#include "updateable_heap.h"
#include "TSC.h"
#include "probab.h"
#include <sys/time.h>
#include <boost/algorithm/string.hpp>


using namespace std;
namespace fs = boost::filesystem;
//typedef boost::unordered_set<std::string> setT;
typedef uint32_t tscKeyT;
typedef boost::unordered_map<tscKeyT, TSC> tscHashT; // key is a TSC ID. // in principle could be int
class transcriptMappingT
{
public:
	std::string name;
	vector<int> starts;
	vector<int> lenghts;
	std::string chr;
	std::string str;
	int beg;
};

//################### global variables #####################

string rootDir ("/import/bc2/home/nimwegen/GROUP/Promoterome/");
int nSamples = 0;
const int maxSamples = 1024;

timeval start,endTime;
boost::array<double,maxSamples> total;
//boost::array<string,maxSamples> samples;

//################### function declarations ################

//void printClusters(list<TSC*> &result, const char* fileName);
int cmpInt (const void * a, const void * b);
list<TSC*>::iterator listErase(list<TSC*> *l, const list<TSC*>::iterator &i );
double& listVal(list<TSC*>::iterator &i );
void readTranscripts(vector<transcriptMappingT> &transcriptMappings, const char* blatF);

/********** CLUSTERING PARAMETERS ***********/
//extern double sigmasq;  // noise size
//extern double alpha;    // gaussian prior size for TSS variance
//extern double lenscale; // coexpression prior scale depending on the distance
//extern double normPseudocount; // in tpms

int main(int argc, char* argv[]) // argv[1] - which chrStr to do
{

	//cerr << sizeof(int) << endl; // 4
	//##########################################################
	//               DATA READING                             ##
	//##########################################################

	string chrStr;
	string resultsDir;
	switch (argc)
	{
	case 4:
		TSC::chromosome = argv[1];
		TSC::strand = argv[2];
		chrStr = TSC::chromosome + TSC::strand;
		resultsDir = argv[3];
		break;
	case 3:
		chrStr = argv[1];
		TSC::chromosome = chrStr;
		TSC::strand = *(TSC::chromosome.end()-1);
		TSC::chromosome.erase(TSC::chromosome.end()-1);
		resultsDir = argv[2];
		break;
	case 2:
		TSC::chromosome = "chr19";
		TSC::strand = "-";
		chrStr = TSC::chromosome + TSC::strand;
		resultsDir = argv[1];
		break;
	}

	vector<transcriptMappingT> transcriptMappings;
	readTranscripts(transcriptMappings, "/import/bc2/home/nimwegen/GROUP/Mara/Annotations/hg19/refSeqAli_onlyNM_hg19.2011-07-29");

	string exprDirSplit = rootDir + resultsDir + "/NormalizedSplit/";
	string mappedDir = rootDir + resultsDir + "/Mapped/";
	string normalizedDir = rootDir + resultsDir + "/Normalized/";
	//fs::path normalizedPath( exprDirSplit + chrStr + "/");
	fs::path normalizedPath (normalizedDir);
	fs::directory_iterator end_iter;
	tscHashT allTSCs;
	int sampleCutoff = 0;
	try
	{
		for ( fs::directory_iterator dir_itr( normalizedPath ); dir_itr != end_iter; ++dir_itr ) //normalizedPath
		{
			if ( fs::is_regular_file( dir_itr->status() ) )
			{
				//if(sampleCutoff++ == 5)
				//	break;
				string file = dir_itr->path().filename().string();
				cerr << file << " " << flush;
				ifstream r ( (exprDirSplit + chrStr + "/" + file).c_str(), ifstream::in );
				if(! r.good())
				{
					cerr << "Could not open: " << (exprDirSplit + chrStr + "/" + file).c_str() << endl;
					return 1;
				}
				r >> skipws;

				while(1) // perhaps it should be .good(), looks like we are reading one more line...
				{
					tscKeyT pos = 0 ;
					double norm = -1;
					double raw = -1;
					char str;
					string chrom;
					r >> chrom;
					r >> str;
					r >> pos;
					r >> raw;
					r >> norm;
					if(r.eof())
						break;
					allTSCs[pos].expr[nSamples].norm = norm;  // we used to havee file == sample as a string key before
					allTSCs[pos].expr[nSamples].raw	 = raw;
					allTSCs[pos].beg = pos;
					allTSCs[pos].end = pos;
					allTSCs[pos].repr = pos;
					allTSCs[pos].reprTot += norm;
				}
				r.close();
				r.open( (mappedDir + file).c_str() );
				if(! r.good())
				{
					cerr << "Could not open: " << (mappedDir + file).c_str() << endl;
					return 1;
				}
				r >> total[nSamples];
				r.close();
				nSamples ++;
				//if(nSamples == 5)
				//{
				//	break;
				//}
			}
		}
	}
	catch (std::exception &e)
	{
		cerr << "Opening some directories or files went wrong. Here is the exception:\n" << e.what() << '\n';
		return 1;
	}
	cerr << "\n" << nSamples << "samples read" << endl;

	//###############################################################
	//              exon and promoter statistics                   ##
	//###############################################################

	BOOST_FOREACH(const transcriptMappingT &mapping, transcriptMappings)
	{
		// WARNING now the strands are wrong!!
		if(strcmp(mapping.chr.c_str(), TSC::chromosome.c_str()) == 0 && strcmp(mapping.str.c_str(), TSC::strand.c_str()) != 0)
		{
			// measure the promoter signal in -50 +50 bp
			double thisExpression = 0.0;
			for(int i=mapping.beg - 50; i<mapping.beg + 50; i++)
			{
				if(allTSCs.find(i) != allTSCs.end())
				{
					//cout << i << endl;
					thisExpression += allTSCs[i].totalExpression();
				}
			}
			thisExpression /= (nSamples);
			// now measure the exon signal:
			int summedExonLen = 0;
			double summedExonSignal = 0.0;
			//for(vector<int>::const_iterator it=mapping.starts.begin(); it != mapping.starts.end(); it++)
			for(unsigned int e=0; e < mapping.starts.size(); e++)
			{
				//cout << "Exon " << mapping.starts[e] << " of len " << mapping.lenghts[e] << endl;
				summedExonLen += mapping.lenghts[e];
				for(int i=0; i<mapping.lenghts[e]; i++)
				{
					//int pos = mapping.str.compare("+") ? (mapping.starts[e] - i) : (mapping.starts[e] + i); // compare !! 
					int pos = mapping.starts[e] + i; // true for both strands!
					if(allTSCs.find(pos) != allTSCs.end())
					{
						summedExonSignal += allTSCs[pos].totalExpression();
					}
				}
				
			}
			
			cout << mapping.name << "\t" << mapping.chr << "\t" << mapping.str << "\t" << mapping.beg << "\t";
			cout << thisExpression << "\t" << summedExonSignal << "\t" << summedExonSignal / (summedExonLen * nSamples)
				<< "\t" << summedExonSignal / ((summedExonLen - 30*mapping.starts.size()) * nSamples) << "\t" << summedExonLen << endl;
			//cout << mapping.chr << endl;
			//cout << mapping.str << endl;
			//sleep(1);
		}
	}
	cerr << "END\n";
	return 0;


	//###############################################################
	//        build a sorted list of the TSCs                      ##
	//###############################################################

	int nTSS = allTSCs.size();
	cerr << nTSS << " TSSs read" << endl;
	int *allPositions = new int[nTSS];
	unsigned int nUsedTSS=0;
	unsigned int nSingletons=0;
	unsigned int nSampleSpecific=0;
	BOOST_FOREACH(const tscHashT::value_type &elem, allTSCs) //loop over TSSs
	{


		//double totTpm = elem.second.totalExpression();
		//if(totTpm >= nSamples * 0.01) // if it passes the cutoff of expression
		//{
		//	allPositions[nUsedTSS++] = elem.first; // that is the position, which is a key.
		//}

		// RUN on chr21-:
		// 304 samples read
		// 2015631 TSSs read
		// Number of singletons (one tag in one sample) = 1519436
		// Number of TSS expressed in only one sample (possibly many tags) = 1531661
		// I will keep everything apart from singletons then!

		allPositions[nUsedTSS++] = elem.first;
//		if(elem.second.expressedInNSamples() == 1)
//		{
//			nSampleSpecific ++;
//			if(elem.second.totalRaw() == 1)
//				nSingletons ++;
//			else // many
//			{
//				allPositions[nUsedTSS++] = elem.first;
//			}
//		}
//		else
//		{
//			allPositions[nUsedTSS++] = elem.first;
//		}
	}
	cerr << "Number of singletons (one tag in one sample) = " << nSingletons << '\n'
			<< "Number of TSS expressed in only one sample (possibly many tags) = " << nSampleSpecific << endl;
	cerr << "Sorting positions" << endl;
	qsort(allPositions, nUsedTSS, sizeof(int), cmpInt);
	cerr << "Initializing" << endl;
	list<TSC*> workTSCs;
	for(unsigned int i = 0; i < nUsedTSS; ++i )
	{
		workTSCs.push_back(&allTSCs[allPositions[i]]); // I hope we are pushing an address of a real object, and not a copy
	}

	if(0) // turn on to do the filtering
	{
	//################################################################
	//                     SCANNING WITH A SLIDING WINDOW           ##
	//################################################################

	list<TSC*>::iterator itBeg = workTSCs.begin();
	list<TSC*>::iterator itEnd = workTSCs.begin();
	itEnd ++;
	int windowLen = 150;

	// produce stats for Erik: (or not)
	list<TSC*>::iterator itSummer; // from ,,summing''
	if(0)
	{
		ofstream movingLog((string("movingLog.") + chrStr).c_str());
		ofstream windowCounts((string("windowCounts.") + chrStr).c_str());
		// init
		int trueWinEnd = (*itBeg)->repr;
		int trueWinBeg = trueWinEnd - windowLen;

		while(1)
		{
			// the left end is the first one > left window position
			// the right one is the first one > right window position
			// if there is a TSS at the true window beg we DONT include it.
			// if there is a TSS at the true window end we DO include it

			//   (---------------------]     window
			//     ^(2)                   ^(3)  will result in
			//     (---------------------]

			// the current situation is NOT accounted yet. we move the window
			// and the jump is the number of positions that the current set of TSS inside
			// is valid

			// decide which one to move

			int jumpLeft = (*itBeg)->repr - trueWinBeg;
			int jumpRight;
			if (itEnd != workTSCs.end())
			{
				jumpRight = (*itEnd)->repr - trueWinEnd;
			}
			else
			{
				jumpRight = 999999999;
			}
			assert(jumpLeft > 0 && jumpRight > 0);
			int thisJump = min(jumpLeft, jumpRight);

			// current situation will span for thisJump, lets do the statistics now:
			// find the highest tpm;

			// if there is at least one TSS inside...
			if((*itBeg)->repr <= trueWinEnd)
			{
				double bestH = 0;
				for(int s=0; s<nSamples; s++)
				{
					double h = 0.0;
					itSummer = itBeg;
					itSummer ++;
					for(; itSummer != itEnd; itSummer++)
					{
						if((*itSummer)->expr.find(s) != (*itSummer)->expr.end())
							h += (*itSummer)->expr[s].norm;
					}
					// the end position:
					if((*itSummer)->expr.find(s) != (*itSummer)->expr.end())
						h += (*itSummer)->expr[s].norm;
					if(bestH < h)
					{
						bestH = h;
					}
				}
				// so we have bestH repeated thisJump times
				for(int i=0; i<thisJump; i++)
					windowCounts << bestH << endl;
			}

			// let's move the window and pointers now
			if(jumpLeft <= jumpRight)
			{
				// we move the left (and maybe right...)
				thisJump = jumpLeft;
				itBeg++;
				if (itBeg == workTSCs.end()) // here we possibly miss some windows in the end, does not matter in stats
					break;
			}
			if(jumpLeft >= jumpRight)
			{
				// we (also) move the right
				thisJump = jumpRight;
				itEnd++;
			}
			trueWinBeg += thisJump;
			trueWinEnd += thisJump;
			movingLog << (*itBeg)->repr << '\t' << (*itEnd)->repr << '\t'
					<< trueWinBeg << '\t' << trueWinEnd << '\t' << thisJump << endl;
		}
		movingLog.close();
		windowCounts.close();
		return 0;
	}
	//########## END OF ERIKS STUFF


	//######################################################################
	//				go with a small window -2 +2 and sum the total tpms
	//######################################################################
	// right now we don't have anything that does not have a high max tpm
	// and we have marked the regions which are above the 20 * max tpm curve
	// now we will add regions which have a large density of total tpms
	// regardless if they are below or under the curve
	windowLen = 7;
	itBeg = workTSCs.begin();
	itEnd = workTSCs.begin();
	itEnd ++;
	int trueWinEnd = (*itBeg)->repr;
	int trueWinBeg = trueWinEnd - windowLen;

	while(1)
	{
		// the left end is the first one > left window position
		// the right one is the first one > right window position
		// if there is a TSS at the true window beg we DONT include it.
		// if there is a TSS at the true window end we DO include it

		//   (---------------------]     window
		//     ^(2)                   ^(3)  will result in
		//     (---------------------]

		// the current situation is NOT accounted yet. we move the window
		// and the jump is the number of positions that the current set of TSS inside
		// is valid

		// decide which one to move

		int jumpLeft = (*itBeg)->repr - trueWinBeg;
		int jumpRight;
		if (itEnd != workTSCs.end())
		{
			jumpRight = (*itEnd)->repr - trueWinEnd;
		}
		else
		{
			jumpRight = 999999999;
			break; // FIXME ending should be fixed
		}
		assert(jumpLeft > 0 && jumpRight > 0);
		int thisJump = min(jumpLeft, jumpRight);

		// current situation will span for thisJump, lets do the statistics now:
		// find the highest tpm;

		// if there is at least one TSS inside...
		if((*itBeg)->repr <= trueWinEnd)
		{
			double h = 0.0;
			itSummer = itBeg;
			itSummer ++;
			for(; itSummer != itEnd; itSummer++)
			{
				h += (*itSummer)->totalExpression();
			}
			// the end position:
			if(itEnd != workTSCs.end())
				h += (*itSummer)->totalExpression();

			// so we have h repeated thisJump times
			// now if we hit the cutoff we paint the TSSs
			if(h > 0.0) // needs to be corrected for number of samples or total depth FIXME
			{
				itSummer = itBeg;
				itSummer ++;
				for(; itSummer != itEnd; itSummer++)
				{
					if(h > (*itSummer)->surrDensity2)
						(*itSummer)->surrDensity2 = h;
				}
				// the end position:
				if(itEnd != workTSCs.end())
					if(h > (*itSummer)->surrDensity2)
						(*itSummer)->surrDensity2 = h;
			}

		}

		// let's move the window and pointers now
		if(jumpLeft <= jumpRight)
		{
			// we move the left (and maybe right...)
			thisJump = jumpLeft;
			itBeg++;
			if (itBeg == workTSCs.end()) // here we possibly miss some windows in the end, does not matter in stats
				break;
		}
		if(jumpLeft >= jumpRight)
		{
			// we (also) move the right
			thisJump = jumpRight;
			itEnd++;
		}
		trueWinBeg += thisJump;
		trueWinEnd += thisJump;
	}

	printClusters(workTSCs, (string() + "tss.win." + chrStr).c_str()); // painting filtered




	///
	// find which windows have a TSS with max expression over cutoff 2.0

	itBeg = workTSCs.begin();
	itEnd = workTSCs.begin();
	int nGoodInside = 0; // (*itBeg)->maxExpression() > 2.0 ? 1 : 0;
	double totExprInside = 0.0;
	double sumMaxExprInside = 0.0;
	double exprCutoff = 0.02717 * nSamples; // 10 was used for human
	windowLen = 150;
	while(1)
	{
		while(itEnd != workTSCs.end() ? (itEnd.operator *()->repr - itBeg.operator *()->repr < windowLen ? 1 : 0) : 0)
		{
			expression tmp = (*itEnd)->maxExpression();
			if((tmp.norm > 2.0 && tmp.raw > 1) || (*itEnd)->surrDensity2 > exprCutoff)
			{
				nGoodInside++;
			}
			if(nGoodInside)
			{
				(*itEnd)->probLogRatio = 1; //just reusing the field for saving space
			}
			totExprInside += (*itEnd)->totalExpression() / nSamples;
			sumMaxExprInside += tmp.norm;
			(*itEnd)->surrDensity += totExprInside / windowLen / 2;
			(*itEnd)->surrDensityOfMax += sumMaxExprInside / windowLen / 2;
			itEnd++;
		}
		if(itEnd == workTSCs.end()) // special for the end
		{
			while(itBeg != workTSCs.end())
			{
				if(nGoodInside)
				{
					(*itBeg)->probLogRatio = 1;
				}
				expression tmp = (*itBeg)->maxExpression();
				if((tmp.norm > 2.0 && tmp.raw > 1) || (*itBeg)->surrDensity2 > exprCutoff)
				{
					nGoodInside--;
				}
				totExprInside -= (*itBeg)->totalExpression() / nSamples;
				sumMaxExprInside -= tmp.norm;
				(*itBeg)->surrDensity += totExprInside / windowLen / 2;
				(*itBeg)->surrDensityOfMax += sumMaxExprInside / windowLen / 2;
				itBeg++;
			}
			break;
		}
		while(itEnd.operator *()->repr - itBeg.operator *()->repr >= windowLen)
		{
			if(nGoodInside)
			{
				(*itBeg)->probLogRatio = 1;
			}
			expression tmp = (*itBeg)->maxExpression();
			if((tmp.norm > 2.0 && tmp.raw > 1) || (*itBeg)->surrDensity2 > exprCutoff)
			{
				nGoodInside--;
			}
			totExprInside -= (*itBeg)->totalExpression() / nSamples;
			sumMaxExprInside -= tmp.norm;
			(*itBeg)->surrDensity += totExprInside / windowLen / 2;
			(*itBeg)->surrDensityOfMax += sumMaxExprInside / windowLen / 2;
			itBeg++;
		}
	}

	printClusters(workTSCs, (string() + "tss.unfiltered." + chrStr).c_str());

	for(itBeg = workTSCs.begin(); itBeg != workTSCs.end(); )
	{
		if((*itBeg)->probLogRatio)
		{
			(*itBeg)->probLogRatio = 0;
			itBeg++;
		}
		else
		{
			itBeg = workTSCs.erase(itBeg);
			nUsedTSS --;
		}
	}


	//#####################################################################
	// now we have the densities and we use them for filtering the TSSs
	//#####################################################################
	double ratioCutoff = 20.0; // 0.054347826 * nSamples; // 20.0 for human  // ratio between current max tpm and the surrounding
	// cannot be that high. 20 seems also high
	double alwaysGoodCutoff = 25.0; // in tpm of max 	// I don't believe there is any exon paining max tpm > 30. per position

	itBeg = workTSCs.begin();
	itEnd = workTSCs.begin();
	while(1)
	{
		while(itEnd != workTSCs.end() ? (itEnd.operator *()->repr - itBeg.operator *()->repr < windowLen ? 1 : 0) : 0)
		{
			expression tmp = (*itEnd)->maxExpression();
			if((tmp.norm > ratioCutoff * (*itEnd)->surrDensityOfMax && tmp.raw > 1) || tmp.norm > alwaysGoodCutoff) // || (*itEnd)->surrDensity2 > 25.0)
			{
				nGoodInside++;
			}
			if(nGoodInside)
			{
				(*itEnd)->probLogRatio = 1; //just reusing the field for saving space
			}
			itEnd++;
		}
		if(itEnd == workTSCs.end()) // special for the end
		{
			while(itBeg != workTSCs.end())
			{
				if(nGoodInside)
				{
					(*itBeg)->probLogRatio = 1;
				}
				expression tmp = (*itBeg)->maxExpression();
				if((tmp.norm > ratioCutoff * (*itBeg)->surrDensityOfMax && tmp.raw > 1) || tmp.norm > alwaysGoodCutoff) // || (*itBeg)->surrDensity2 > 25.0)
				{
					nGoodInside--;
				}
				itBeg++;
			}
			break;
		}
		while(itEnd.operator *()->repr - itBeg.operator *()->repr >= windowLen)
		{
			if(nGoodInside)
			{
				(*itBeg)->probLogRatio = 1;
			}
			expression tmp = (*itBeg)->maxExpression();
			if((tmp.norm > ratioCutoff * (*itBeg)->surrDensityOfMax && tmp.raw > 1) || tmp.norm > alwaysGoodCutoff) // || (*itBeg)->surrDensity2 > 25.0)
			{
				nGoodInside--;
			}
			itBeg++;
		}
	}
	//##################################################################
	//   add to the good list (probLogRatio) the ones which           ##
	//   have an increase in surrDensityOfMax at the scale of 1000bp  ##
	//##################################################################

	windowLen = 1000;
	itBeg = workTSCs.begin();
	itEnd = workTSCs.begin();
	itEnd ++;
	trueWinEnd = (*itBeg)->repr;
	trueWinBeg = trueWinEnd - windowLen;

	while(1)
	{
		// the left end is the first one > left window position
		// the right one is the first one > right window position
		// if there is a TSS at the true window beg we DONT include it.
		// if there is a TSS at the true window end we DO include it

		//   (---------------------]     window
		//     ^(2)                   ^(3)  will result in
		//     (---------------------]

		// the current situation is NOT accounted yet. we move the window
		// and the jump is the number of positions that the current set of TSS inside
		// is valid

		// decide which one to move

		int jumpLeft = (*itBeg)->repr - trueWinBeg;
		int jumpRight;
		if (itEnd != workTSCs.end())
		{
			jumpRight = (*itEnd)->repr - trueWinEnd;
		}
		else
		{
			jumpRight = 999999999;
			break; // FIXME ending should be fixed
		}
		assert(jumpLeft > 0 && jumpRight > 0);
		int thisJump = min(jumpLeft, jumpRight);

		// current situation will span for thisJump, lets do the statistics now:
		// find the highest tpm;

		// if there is at least one TSS inside...
		if((*itBeg)->repr <= trueWinEnd)
		{
			double minVal = 1000000.0;
			itSummer = itBeg; // this time the summer does not sum but find the minimum
			itSummer ++;
			for(; (itSummer != itEnd) || (itSummer != itEnd && itEnd != workTSCs.end()); itSummer++)
			{
				if(minVal > (*itSummer)->surrDensity)
					minVal = (*itSummer)->surrDensity;
				if((*itSummer)->surrDensity > 20 * minVal)
					(*itSummer)->dropLeft = true;
			}
			// and the same opposite direction!
			itSummer = itEnd; // no need for --
			minVal = 10000000.0;
			for(; (itSummer != itBeg); itSummer--)
			{
				if(minVal > (*itSummer)->surrDensity)
					minVal = (*itSummer)->surrDensity;
				if((*itSummer)->surrDensity > 20 * minVal)
					(*itSummer)->dropRight = true;
			}
		}

		// let's move the window and pointers now
		if(jumpLeft <= jumpRight)
		{
			// we move the left (and maybe right...)
			thisJump = jumpLeft;
			itBeg++;
			if (itBeg == workTSCs.end()) // here we possibly miss some windows in the end, does not matter in stats
				break;
		}
		if(jumpLeft >= jumpRight)
		{
			// we (also) move the right
			thisJump = jumpRight;
			itEnd++;
		}
		trueWinBeg += thisJump;
		trueWinEnd += thisJump;
	}

	//###################################################################
	//       go over the list of TSSs and remove the non-good ones     ##
	//###################################################################
	//printClusters(workTSCs, "tss.paintigF"); // painting filtered

	for(itBeg = workTSCs.begin(); itBeg != workTSCs.end(); )
	{
		expression tmp = (*itBeg)->maxExpression();
		if((*itBeg)->probLogRatio || ((*itBeg)->dropLeft && (*itBeg)->dropRight))
		{
			(*itBeg)->probLogRatio = 0;
			itBeg++;
		}
		else
		{
			itBeg = workTSCs.erase(itBeg);
			nUsedTSS --;
		}
	}

	printClusters(workTSCs, (string() + "tss.final." + chrStr).c_str());
	}
	//################################################################
	//                BUILD THE HEAP                                ##
	//################################################################

	updateable_heap<heapElemT> clusteringQueue; // keeps probabilities to cluster
	list<TSC*>::iterator it = workTSCs.begin();
	cerr << nUsedTSS <<  " TSS after filtering" << endl;
	for(unsigned int i=0; i < nUsedTSS-1; i++)
	{
		double priorRatio =  priorprob(it);
		double likelihoodRatio = profProbAllSamples(it);
		double p = priorRatio + likelihoodRatio;
		(*it)->probLogRatio = p;
		clusteringQueue.push(it);
		it++;
	}
	// deal with the last element (never fuse)
	(*it)->probLogRatio = -999999999.0;
	clusteringQueue.push(it);


	//##################################################
	//############### CLUSTERING LOOP ##################
	//##################################################

	cerr << "Entering the clustering loop" << endl;
	// get the clustering curve...
	int nTSC = nUsedTSS;
	ofstream fusingLog((string() +"fusions." + chrStr).c_str());
	gettimeofday(&endTime, NULL);
	while(nTSC > 2)
	{
		//cerr << "Clusters left " << nTSC << endl;

		if(!(nTSC % 1000))
		{
			start = endTime;
			gettimeofday(&endTime, NULL);
			cerr << "Clusters left " << nTSC << "\t" <<
					(endTime.tv_sec - start.tv_sec) *1000 +
					(endTime.tv_usec - start.tv_usec)/1000.0 + 0.5 <<
					endl;
		}
		//if(workTSCs.size() != clusteringQueue.size()) // size might be linear!!
		//{
		//	cerr << "discrepancy! " << workTSCs.size() << " " << (clusteringQueue.size()) << endl;
		//	exit(1);
		//}
		list<TSC*>::iterator itBestTSC = clusteringQueue.peek().left;
		if(isinf((*itBestTSC)->probLogRatio))
		{
			cerr << "queue size " << clusteringQueue.size() << endl;
			cerr << "tsc queue " << workTSCs.size() << endl;
			exit(1);
		}
		if((*itBestTSC)->probLogRatio < 0)
		{
			break;
			// we just finished clustering
		}
		list<TSC*>::iterator itNextTSC = itBestTSC;
		itNextTSC++;
		if(itNextTSC == workTSCs.end())
		{
			cerr << "ERROR: We are trying to fuse the last TSC with non-existing right TSC." << endl;
			exit(1);
		}
		//cerr << "Fusing:\tp=" << (*itBestTSC)->probLogRatio << flush << "\t" << (*itBestTSC)->beg << flush << "\t" << (*itBestTSC)->end << endl;
		//cerr << "And   :\tp=" << (*itNextTSC)->probLogRatio << flush << "\t" << (*itNextTSC)->beg << flush << "\t" << (*itNextTSC)->end << endl;
		fusingLog << (*itBestTSC)->probLogRatio << '\t' <<  ((*itNextTSC)->beg - (*itBestTSC)->end) << '\n';
		**itBestTSC += **(itNextTSC); // fusion

		clusteringQueue.erase((*itNextTSC)->posInHeap);
		//itNextTSC = workTSCs.erase(itNextTSC);  // this always exists, doesn't invalidate the iterator
		itNextTSC = listErase(&workTSCs, itNextTSC);

		// update the neighbouring fusion probabilities
		if(itBestTSC != workTSCs.begin())
		{
			//cerr << "Left pair update" << endl;
			list<TSC*>::iterator itPrevTSC = itBestTSC;
			itPrevTSC--;
			double priorRatio =  priorprob(itPrevTSC);
			double likelihoodRatio = profProbAllSamples(itPrevTSC);
			listVal(itPrevTSC) = priorRatio + likelihoodRatio;
			//(*itPrevTSC)->probLogRatio = priorRatio + likelihoodRatio;
			// cerr << "We updated " << (*itPrevTSC)->posInHeap << endl;
			clusteringQueue.updated((*itPrevTSC)->posInHeap);
		}
		if(itNextTSC != workTSCs.end())
		{
			// cerr << "Right pair update" << endl;
			double priorRatio =  priorprob(itBestTSC);
			double likelihoodRatio = profProbAllSamples(itBestTSC);
			//(*itBestTSC)->probLogRatio = priorRatio + likelihoodRatio;
			listVal(itBestTSC) = priorRatio + likelihoodRatio;
			//cerr << "Updating " << (*itBestTSC)->posInHeap << endl;
			clusteringQueue.updated((*itBestTSC)->posInHeap);
		}
		else // we just removed the last element. current last is Best.
		{
			//cerr << "Dealing with the last element" << endl;
			listVal(itBestTSC) =  -999999999;
			//(*itBestTSC)->probLogRatio = -999999999;
			clusteringQueue.updated((*itBestTSC)->posInHeap);
		}
		nTSC --;
	}
	//printClusters(workTSCs, (string() + "tsc.final." + chrStr).c_str(), 2.0);
	printClustersBed(workTSCs, (string() + "TSC." + chrStr).c_str(), 2.0);
	delete[] allPositions;
	return 0;
}

//########################################################
//                 FUNCTIONS                            ##
//########################################################

inline list<TSC*>::iterator listErase(list<TSC*> *l, const list<TSC*>::iterator &i )
{
	return l->erase(i);
}

inline double& listVal(list<TSC*>::iterator &i )
{
	return (*i)->probLogRatio;
}


int cmpInt (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

void readTranscripts(vector<transcriptMappingT> &transcriptMappings, const char* blatF)
{
	ifstream blat(blatF);
	int foo, matches, mismatches, repmatches, e, f, g, h, len, qsize, qstart, qend, chrlen, beg, end, blockCount;
	char str;
	char ref[32];
	char chr[32];
	char fooStr[10000];
	std::vector<std::string> starts;
	std::vector<std::string> lenghts;

	while(blat.good())
	{
		blat >> foo >> matches >> mismatches >> repmatches >> e >> f >> g >> h
			>> len >> str >> ref >> qsize >> qstart >> qend >> chr >> chrlen >> beg >> end >> blockCount >> fooStr;
		//cout << "A" << blat.good();
		// now we have block sizes in foo string
		boost::split(lenghts, fooStr, boost::is_any_of(","));
		blat >> fooStr;
		//cout << "B" << blat.good();
		// now starts in the query
		blat >> fooStr;
		//cout << "C" << blat.good();
		// now starts in the chromosome
		boost::split(starts, fooStr, boost::is_any_of(","));
		transcriptMappingT thisMapping;
		thisMapping.name = ref;
		thisMapping.chr = chr;
		thisMapping.str = str;
		if(qstart > 5)
			continue;
		thisMapping.beg = (str == '+') ? beg : end;
		BOOST_FOREACH(const std::string &len, lenghts)
		{
			if(len.length())
			{
				thisMapping.lenghts.push_back(atoi(len.c_str()));
			}
		}
		BOOST_FOREACH(const std::string &start, starts)
		{
			if(start.length())
			{
				thisMapping.starts.push_back(atoi(start.c_str()));
			}
		}
		transcriptMappings.push_back(thisMapping);
		
	}
	blat.close();
	cerr << transcriptMappings.size() << " mappings read" << endl;

}

