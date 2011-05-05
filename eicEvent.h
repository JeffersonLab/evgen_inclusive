#ifndef __eicEvent_h
#define __eicEvent_h

#include "TObject.h"
#include "TLorentzVector.h"

struct eventdata {
	double weight;
	double Wpweight;
	double Wmweight;
	double theta;
	double jetth;
	double phi;
	double ef;
	double x;
	double W;
	double y;
	double crs;
	double f1;
	double f2;
	double f1gz;
	double f3gz;
	double g1gz;
	double g5gz;
	double f1wp;
	double f3wp;
	double g1wp;
	double g5wp;
	double f1wm;
	double f3wm;
	double g1wm;
	double g5wm;
	double Q2;
	double eta_gZ;
	double rate;
	double Wprate;
	double Wmrate;
	double pdf[6];
	double Dpdf[5];

	double Abeam;
	double Y;

	double A_L;
	double A_L_g1;
	double A_L_g5;
	double Y_L;
	double eta_L;
	
	double A_T;
	double Y_T;
	double eta_T;

	double AWp;
	double AWm;

	double AWp_g1;
	double AWp_g5;
	double AWm_g1;
	double AWm_g5;

	double jetE;
	double jetp;

        int particle_id;
        int charge;
        int Z_ion;
        int N_ion;
        double mom;
        double mass;

        TLorentzVector pi0_g1;
        TLorentzVector pi0_g2;
  
};

class eicEvent {
    public:
	eicEvent();
	~eicEvent();

	void SetEventData( eventdata d ){ fData = d; }
	eventdata GetEventData(){ return fData; }

    private:
	eventdata fData;
	
    public:
//	ClassDef(eicEvent,1);
};
#endif//__eicEvent_h
