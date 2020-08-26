#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include "SampleDecoder.h"
#include "MTRand.h"
#include "BRKGA.h"
#include "bossa_timer.h"
#include "ArgPack.hpp"
#include <fstream>

using namespace std;

vector<vector<double> > points;

int main(int argc, char* argv[]) {

	ArgPack single_ap(argc, argv);

	unsigned n = 0;			// size of chromosomes
	const unsigned p = ArgPack::ap().population;		// size of population
	const double pe = ArgPack::ap().populationElite;		// fraction of population to be the elite-set
	const double pm = ArgPack::ap().populationMutants;		// fraction of population to be replaced by mutants
	const double rhoe = ArgPack::ap().rhoe;	// probability that offspring inherit an allele from elite parent
	const unsigned K = ArgPack::ap().K;		// number of independent populations
	const unsigned MAXT = ArgPack::ap().threads;	// number of threads for parallel decoding

	const double cutoff_time = ArgPack::ap().time;

	// Reading instance
	string s;
	string file = ArgPack::ap().inputFile;

	ifstream f(file);

	int n_points, dim;
	f >> n_points >> dim;

	points = vector<vector<double> > (n_points);
	for(int i=0; i<n_points; ++i){
		points[i] = vector<double> (dim);
	}

	for(int i=0; i<n_points ;++i) {
		for(int d=0; d < dim ; ++d) {
			f >> points[i][d];
		}
	}

	SampleDecoder decoder;			// initialize the decoder
	const long unsigned rngSeed = ArgPack::ap().rngSeed;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator

	n = n_points;

	// initialize the BRKGA-based heuristic
	BossaTimer timer;

	double bestValue = -1;
	double timerToBest;
	bool verbose = ArgPack::ap().verbose;
	int k_max = sqrt(n_points); // 2 <= k <= sqrt(n)
	int number_pop = k_max - 1;
	vector<BRKGA<SampleDecoder, MTRand>*> NPopulations(number_pop);
	vector<SampleDecoder*> vec_decoders(number_pop);
	vector<double> best_values(number_pop);

	// Random select a number between [0,number_permutation) 
	vector<int> indexs;

	int best_population = 0;
	double sum_best= 0.0;

	const unsigned X_INTVL =  ArgPack::ap().exchangeBest;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = ArgPack::ap().exchangeTop;	// exchange top 2 best
	const unsigned MAX_GENS = ArgPack::ap().generations;	// run for 1000 gens
	

	if(verbose){
		cout << "Initializing populations\n";
	}

	timer.start(); 
	double best_temp;
	for(int i=0; i < number_pop; ++i) {
		vec_decoders[i] = new SampleDecoder();
		vec_decoders[i]->set_k(i+2);
		NPopulations[i] = new BRKGA<SampleDecoder, MTRand> (n, p, pe, pm, rhoe, *vec_decoders[i], rng, K, MAXT);
		best_temp = (-1)*NPopulations[i]->getBestFitness();
		if(best_temp > bestValue) {
			timerToBest = timer.getTime();
			best_population = i;
			bestValue = best_temp;
		}
	}

	// Tournament
	vector<int> aux;
	unsigned generation = 1;		// current generation
	unsigned tournament = 1;

	cout << "Cut off " << cutoff_time << endl; 
	do {
		indexs.clear();
		indexs.reserve(number_pop);
		for( int i = 0 ; i < number_pop ; ++i ){
			indexs.push_back(i);
		}
		random_shuffle(indexs.begin(), indexs.end());
		cout << "Tournament " << tournament << endl;
		while (indexs.size() > 1 and timer.getTime() < cutoff_time){
			
			cout << "[" << generation << "] Sequence:";
			for( int i = 0 ; i < indexs.size() ; ++i ){
				cout << " " << indexs[i];
			}
			cout << endl;

			swap(aux,indexs);
			indexs.clear();

			int size = aux.size();
			for (int i = 0; i < size/2; ++i) {
				// Battle: 2*i vs 2*i + 1
				NPopulations[aux[2*i]]->evolve();
				NPopulations[aux[2*i+1]]->evolve();
				// TODO: Olhar se o operator está correto
				// Caso queira, se der empate escolher aleatoriamente um dos dois
				if( NPopulations[aux[2*i]]->getBestFitness() < NPopulations[aux[2*i + 1]]->getBestFitness() ){
					indexs.push_back(aux[2*i]);
				} else {
					indexs.push_back(aux[2*i+1]);
				}
			}
			if ( size & 1 ){
				NPopulations[aux[size-1]]->evolve();
				indexs.push_back(aux[size-1]);
			}
			generation++;
		}
		generation++;
		tournament++;
		if(timer.getTime() < cutoff_time) NPopulations[indexs[0]]->evolve();
	} while(generation < MAX_GENS and timer.getTime() < cutoff_time);

	cout << "Best solution " << bestValue << " from pop " << best_population << " k = " << best_population +2 << endl;
	cout << "Total time = " << timer.getTime() << endl;
	cout << "Time to Best ttb = " << timerToBest << endl;

	for(int i=0;i<number_pop; ++i) {
		delete NPopulations[i];
		delete vec_decoders[i];
	}


	return 0;
}
