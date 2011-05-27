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
    double lenx;
    double leny;

    double off_x;
    double off_y;
    double off_z;

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
	double Get_lenx(){ return fData.lenx; }
	double Get_leny(){ return fData.leny; }
	int    Get_model(){ return fData.model; }
	int    Get_fmt(){ return fData.out_fmt; }

	double GetTgtXoff(){ return fData.off_x; }
	double GetTgtYoff(){ return fData.off_y; }
	double GetTgtZoff(){ return fData.off_z; }
	


    private:
	inputdata fData;
};
#endif//__eicInput_h
