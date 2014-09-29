#include "eicIon.h"

#include <math.h>

eicIon::eicIon(eicInput *inp):
    fMom(0), fTheta(0), fPhi(0),fMass(0),fEnergy(0),fZ(0), fN(0)
{
    if( inp ){      
	fMom = inp->Get_ionMom();
	fTheta = inp->Get_ionTheta();
	fPhi = inp->Get_ionPhi();	
	fMass   = inp->Get_ionMass();
	fEnergy = sqrt(fMom*fMom+fMass*fMass);
	
	fZ      = inp->Get_ionZ();
	fN      = inp->Get_ionN();
	
	printf("ion Mom %f GeV, ion Theta %f rad, ion Phi %f rad, ion mass %f GeV\n",fMom, fTheta,fPhi,fMass);	
	printf("ion Z %i, ion N %i\n", fZ, fN);		
    }
    return;
}

eicIon::~eicIon(){
    return;
}
