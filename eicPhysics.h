#ifndef __eicPhysics_h
#define __eicPhysics_h

#include "cteqpdf.h"

#include "eicBeam.h"
#include "eicIon.h"
#include "eicEvent.h"


#include "TRandom2.h"

#define MASS_P 0.938


class eicPhysics {
    public:
	eicPhysics();
	~eicPhysics();

	void MakeEvent( eicBeam *, eicIon *, eicEvent * );

	enum nucl {kProton, kNeutron};
    private:
	 cteq_pdf_t *fPDF;

	 TRandom2 *fRandom;

	 double GetTheta( double min = 0, double max = 3.14159 );
	 double GetPhi(   double min = 0, double max = 3.14159*2 );
	 double GetEf( double, double );
	 double GetEfmax( double, double );

	 double F2(double x, double Q2, nucl n );
	 double F1(double x, double Q2, nucl n );

	 double F1gZ(double x, double Q2, nucl n );
	 double F3gZ(double x, double Q2, nucl n );
	 double g1gZ(double x, double Q2, nucl n );
	 double g5gZ(double x, double Q2, nucl n );

	 void ReadPolTable();
	 double getDeltaq( double x, double Q2, int quark );

	 int fNx;
	 int fNQ2;
	 double fDeltaq[5][50][50];
	 double fDeltaqx[50][50];
	 double fDeltaqQ2[50][50];
};
#endif//__eicPhysics_h
