#ifndef __eicModel_h
#define __eicModel_h

#include "eicInput.h"

// EIC electron beam

class eicModel{
    public:
	eicModel(eicInput *inp = 0);
	~eicModel();

	int GetModel() { return model; }
	void   SetModel(int m) {  model = m; }

    private:
	int model;
};
#endif//__eicModel_h
