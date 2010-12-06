#ifndef __eicIon_h
#define __eicIon_h

#include "eicInput.h"

// EIC Ion beam

class eicIon{
    public:
	eicIon(eicInput *inp = 0);
	~eicIon();

	double GetEnergy() { return fEnergy; }
	void   SetEnergy(double e) {  fEnergy = e; }
	double GetLumin() { return fLumin; }
	void   SetLumin(double l) {  fLumin = l; }

	double GetMass() { return fMass; }
	int    GetZ() { return fZ; }
	int    GetN() { return fN; }

	void   SetMass( double m ){ fMass = m; }
	void   SetZ( int z ){ fZ = z; }
	void   SetN( int n ){ fN = n; }

	enum emcType {kNone};

	emcType GetEMCtype(){ return fEMC; }
	void SetEMCtype(emcType t){ fEMC = t; }

    private:
	double fLumin;
	double fEnergy;
	double fMass;
	int    fZ;
	int    fN;

	emcType fEMC;
};
#endif//__eicIon_h
