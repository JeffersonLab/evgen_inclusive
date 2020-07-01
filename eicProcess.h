#ifndef __eicProcess_h
#define __eicProcess_h

#include "eicBeam.h"
#include "eicIon.h"
#include "eicPhysics.h"
#include "eicInput.h"
#include "eicOutput.h"
#include "eicEvent.h"
#include "eicModel.h"
#include "TString.h"


class eicProcess {
    public:
         eicProcess(const char *file,char *file2, double nevent, int whichmodel);
	~eicProcess();

	void Run();

    private:

	eicInput    *finp;
	eicOutput   *fout;
	eicBeam     *fbeam;
	eicIon      *fion;
	eicPhysics  *fphy;
	eicEvent    *fevt;
	eicModel    *fmodel;
	
	int nevt;
	int model;
};
#endif//__eicProcess_h
