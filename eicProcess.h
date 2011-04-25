#ifndef __eicProcess_h
#define __eicProcess_h

#include "eicBeam.h"
#include "eicIon.h"
#include "eicPhysics.h"
#include "eicInput.h"
#include "eicOutput.h"
#include "eicEvent.h"
#include "eicModel.h"


class eicProcess {
    public:
         eicProcess(const char *file,char *file2);
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

};
#endif//__eicProcess_h
