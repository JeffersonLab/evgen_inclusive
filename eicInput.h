#ifndef __eicInput_h
#define __eicInput_h

#include "TObject.h"

struct inputdata {
    int    nevt;
    int    nprnt;
    double lumin;
    double runtime;
    double e_energy;
    double ion_energy;
    double ion_mass;
    int    ion_Z;
    int    ion_N;
    double dens;
    double length;
    int	   model;  // 1= electron ; 2=pi+ ; 3=pi- ; 4=pi0 
    int    out_fmt; // 1=root 2=root+LUND 3=root+SOLLUND

};

class eicInput {
    public:
	 eicInput(const char *file);
	~eicInput();

	int    GetNevt(){ return fData.nevt; }
	int    GetNprnt(){ return fData.nprnt; }
	double GetLumin(){ return fData.lumin; }
	double GetRunTime(){ return fData.runtime; }
	double Get_eEnergy(){ return fData.e_energy; }
	double Get_ionEnergy(){ return fData.ion_energy; }
	double Get_ionMass(){ return fData.ion_mass; }
	int    Get_ionZ(){ return fData.ion_Z; }
	int    Get_ionN(){ return fData.ion_N; }
	double Get_dens(){ return fData.dens; }
	double Get_length(){ return fData.length; }
	int    Get_model(){ return fData.model; }
	int    Get_fmt(){ return fData.out_fmt; }
	


    private:
	inputdata fData;
};
#endif//__eicInput_h
