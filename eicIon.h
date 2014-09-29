#ifndef __eicIon_h
#define __eicIon_h

#include "eicInput.h"

// EIC Ion beam

class eicIon{
    public:
	eicIon(eicInput *inp = 0);
	~eicIon();
	
	double GetMom() { return fMom; }
	double GetTheta() { return fTheta; }
	double GetPhi() { return fPhi; }	
	double GetMass() { return fMass; }
	double GetEnergy() { return fEnergy; }	
	int    GetZ() { return fZ; }
	int    GetN() { return fN; }
	
	void   SetMom(double p) {  fMom = p; }	
	void   SetTheta( double theta ){ fTheta = theta; }	
	void   SetPhi( double phi ){ fPhi = phi; }	
	void   SetMass( double m ){ fMass = m; }
	void   SetEnergy(double e) {  fEnergy = e; }			
	void   SetZ( int z ){ fZ = z; }
	void   SetN( int n ){ fN = n; }

	enum emcType {kNone};

	emcType GetEMCtype(){ return fEMC; }
	void SetEMCtype(emcType t){ fEMC = t; }

    private:
	double fMom;
	double fTheta;
	double fPhi;
	double fMass;	
	double fEnergy;
	int    fZ;
	int    fN;

	emcType fEMC;
};
#endif//__eicIon_h
