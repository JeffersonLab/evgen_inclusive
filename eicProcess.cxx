#include "eicProcess.h"

eicProcess::eicProcess(const char *file, char *file2){
    finp   = new eicInput(file);

    fphy  = new eicPhysics();
    fbeam = new eicBeam(finp);
    fion  = new eicIon(finp);

    fout  = new eicOutput(finp,file2);

    fevt  = new eicEvent();
    return;
}

void eicProcess::Run(){
    int nevt = finp->GetNevt();

    printf("nevt = %d\n", nevt);

    int evt;
    int nprnt = finp->GetNprnt();

    for( evt = 0; evt < nevt; evt++ ){
	if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	fphy->MakeEvent(fbeam, fion, fevt);
	fout->Write(fevt);
    }

    fout->Close();

    return;
}

eicProcess::~eicProcess(){
    return;
}
