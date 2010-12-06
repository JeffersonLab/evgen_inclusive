#ifndef __eicEvent_h
#define __eicEvent_h

#include "TObject.h"

struct eventdata {
	double weight;
	double theta;
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
	double Q2;
	double eta_gZ;
	double rate;
	double pdf[6];
	double Dpdf[5];

	double Abeam;
	double Y;

	double A_L;
	double Y_L;
	double eta_L;
	
	double A_T;
	double Y_T;
	double eta_T;
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
