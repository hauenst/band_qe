#ifndef __SP_WF_H__
#define __SP_WF_H__

#include "TSpline.h"
#include "TGraph.h"

class spWf
{
	public:
		spWf();
		~spWf();

		double getDensity(double k);

	private:
		double density[201];
		double kPts[201];
		TSpline3 * splDensity;
		void fill_arrays();

};

#endif
