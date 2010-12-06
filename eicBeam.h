#ifndef __eicBeam_h
#define __eicBeam_h

#include "eicInput.h"

// EIC electron beam

class eicBeam{
    public:
	eicBeam(eicInput *inp = 0);
	~eicBeam();

	double GetEnergy() { return fEnergy; }
	void   SetEnergy(double e) {  fEnergy = e; }
	double GetLumin() { return fLumin; }
	void   SetLumin(double l) {  fLumin = l; }
    private:
	double fLumin;
	double fEnergy;
};
#endif//__eicBeam_h
