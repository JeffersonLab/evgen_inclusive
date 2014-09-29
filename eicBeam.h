#ifndef __eicBeam_h
#define __eicBeam_h

#include "eicInput.h"

// EIC electron beam

class eicBeam{
    public:
	eicBeam(eicInput *inp = 0);
	~eicBeam();

	double GetMom() { return fMom; }
	double GetTheta() { return fTheta; }
	double GetPhi() { return fPhi; }	
	double GetMass() { return fMass; }
	double GetEnergy() { return fEnergy; }	
	double GetLumin() { return fLumin; }	
	
	void   SetMom(double p) {  fMom = p; }	
	void   SetTheta( double theta ){ fTheta = theta; }	
	void   SetPhi( double phi ){ fPhi = phi; }	
	void   SetMass( double m ){ fMass = m; }
	void   SetEnergy(double e) {  fEnergy = e; }		
	void   SetLumin(double l) {  fLumin = l; }
	
    private:
	double fMom;
	double fTheta;
	double fPhi;
	double fMass;
	double fEnergy;	
	double fLumin;	
};
#endif//__eicBeam_h
