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
	double GetLx() {return len_x; }
	double GetLy() {return len_y; }
	void   SetModel(int m) {  model = m; }
	void   SetRadLen(double r) { radlen = r; }

    private:
	int model;
	double radlen;
	double length;
	double len_x;
	double len_y;
};
#endif//__eicModel_h
