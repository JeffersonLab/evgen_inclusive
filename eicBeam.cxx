#include "eicBeam.h"

#include <math.h>

eicBeam::eicBeam(eicInput *inp):fMom(0),fTheta(0),fPhi(0),fMass(0),fEnergy(0),fLumin(0){
    if( inp ){
	fMom = inp->Get_eMom();
	fTheta = inp->Get_eTheta();
	fPhi = inp->Get_ePhi();
	fMass = inp->Get_eMass();
	fEnergy = sqrt(fMom*fMom+fMass*fMass);
	
	fLumin  = inp->GetLumin();	
	
	printf("Luminosity =  %e Hz/m^2\n", fLumin);
	printf("e Mom %f GeV, e Theta %f rad, e Phi %f rad, e mass %f GeV\n",fMom, fTheta,fPhi,fMass);		
    }

    return;
}

eicBeam::~eicBeam(){
    return;
}
