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
    TString modelst[10];
    modelst[0] = "electron DIS";
    modelst[1] = "pi+";
    modelst[2] = "pi-";
    modelst[3] = "pi0";
    modelst[4] = "electron elastic";        
    modelst[5] = "K+";
    modelst[6] = "K-";
    modelst[7] = "Ks";
    modelst[8] = "Kl";
    modelst[9] = "p";   
    TString fmtst[3];    
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
