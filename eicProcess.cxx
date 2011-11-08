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
    TString fmtst[3];
    TString modelst[4];
    modelst[0] = "Electron DIS";
    modelst[1] = "Pi+";
    modelst[2] = "Pi-";
    modelst[3] = "Pi0";
    fmtst[0] = "ROOT";
    fmtst[1] = "ROOT and LUND";
    fmtst[2] = "ROOT and SOLLUND";
    int fmt = finp->Get_fmt();
    int model = fmodel->GetModel();

    printf("Model = %s \n", modelst[(model-1)].Data() );
    printf("Format = %s \n", fmtst[(fmt-1)].Data() );

    if (model == 1) { // electron DIS generator
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    else if (model == 5) { // electron ELASTIC generator
      printf("Elastic cross section \n");
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent3(fbeam, fion, fevt, fmodel);
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
