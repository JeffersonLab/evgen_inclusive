#include "eicIon.h"

eicIon::eicIon(eicInput *inp):
    fLumin(0), fEnergy(0), fMass(0),  fZ(0), fN(0)
{
    if( inp ){
	fEnergy = inp->Get_ionEnergy();
	fMass   = inp->Get_ionMass();
	fZ      = inp->Get_ionZ();
	fN      = inp->Get_ionN();
//	printf("ion mass = %f\n", fMass );
    }
    return;
}

eicIon::~eicIon(){
    return;
}
