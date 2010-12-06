#include "eicBeam.h"

eicBeam::eicBeam(eicInput *inp):fLumin(0),fEnergy(0){
    if( inp ){
	fLumin  = inp->GetLumin();
	fEnergy = inp->Get_eEnergy();

	printf("Luminosity =  %e Hz/m^2; E_beam = %f GeV\n", fLumin, fEnergy);
    }

    return;
}

eicBeam::~eicBeam(){
    return;
}
