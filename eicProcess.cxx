#include "eicProcess.h"

eicProcess::eicProcess(const char *file, char *file2){
    finp   = new eicInput(file);

    fphy  = new eicPhysics();
    fbeam = new eicBeam(finp);
    fion  = new eicIon(finp);

    fout  = new eicOutput(finp,file2);

    fevt  = new eicEvent();
    fmodel = new eicModel(finp);
    return;
}

void eicProcess::Run(){
    int nevt = finp->GetNevt();

    printf("nevt = %d\n", nevt);

    int evt;
    int nprnt = finp->GetNprnt();
    int fmt = finp->Get_fmt();
    int model = fmodel->GetModel();

    printf("Model = %d \n", model);
    printf("Format = %d \n", fmt);

    if (model == 1) { // electron DIS generator
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    else { // pion generator
      printf("Pion cross section \n");
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent2(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    fout->Close();

    if (fmt == 2 ) { // also LUND format 
      printf("Creating LUND output \n") ;
      fout->MakeFileLUND();


    }

    if (fmt == 3 ) { // also SOLLUND format
      printf("Creating SOLLUND output \n") ;
      fout->MakeFileSOLLUND();

    } 
    

    return;
}

eicProcess::~eicProcess(){
    return;
}
