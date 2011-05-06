#ifndef __eicModel_h
#define __eicModel_h

#include "eicInput.h"
#include <math.h>
#include <cstdlib>

// EIC electron beam

class eicModel{
    public:
	eicModel(eicInput *inp = 0);
	~eicModel();

	int    GetModel() { return model; }
	double GetRadLen() {return radlen; }
	double GetLength() {return length; }
	double GetRadius() {return radius; }
	void   SetModel(int m) {  model = m; }
	void   SetRadLen(double r) { radlen = r; }

    private:
	int model;
	double radlen;
	double length;
	double radius;
};
#endif//__eicModel_h
