#include "spWf.h"

#include "constants.h"
#include "TH1.h"

#include <iostream>
#include <cstdlib>

#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

spWf::spWf(){
	fill_arrays();
}

spWf::~spWf(){}

void spWf::fill_arrays(){


	std::ifstream f_corr("../include/2H_momDist.txt");

	double testInt = 0.;	

	char line[256]; 
	for( int i = 0 ; i < 29 ; i ++) f_corr.getline(line, 256); 
	double k, rhokS, rhokD, rhoK;
	int ctr = 0;
	if (f_corr.is_open() == true){ 
		while (!f_corr.eof()) { 
			f_corr >> k >> rhokS >> rhokD >> rhoK; 
			double bin = k/0.05;
			int BIN = round(bin);
			density[BIN] = rhoK; 	// fm^3
			kPts[BIN] = k;		// fm^-1
		} 
	} 
	else{ 
		cout << "\tFailed to find wavefunction...\n";
	}  

	splDensity = new TSpline3("splDensity",kPts,density,201,"b2e2",0,0);
}

double spWf::getDensity( double k ){
	// Assumes k is in GeV
	double k_inFermi = k / (GeVfm);
	return splDensity->Eval(k_inFermi) / ( pow(GeVfm,3) );
}

