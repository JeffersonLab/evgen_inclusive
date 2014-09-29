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
    TString modelst[23];
    modelst[0] = "Pi+";
    modelst[1] = "Pi-";
    modelst[2] = "Pi0";
    modelst[3] = "K+";
    modelst[4] = "K-";
    modelst[5] = "Ks";
    modelst[6] = "Kl";
    modelst[7] = "p";
    modelst[8] = "p-bar";    

    modelst[10] = "Pi+ Uniform";
    modelst[11] = "Pi- Uniform";
    modelst[12] = "Pi0 Uniform";
    modelst[13] = "K+ Uniform";
    modelst[14] = "K- Uniform";
    modelst[15] = "Ks Uniform";
    modelst[16] = "Kl Uniform";
    modelst[17] = "p Uniform";
    modelst[18] = "p-bar Uniform";   
    
    modelst[20] = "Electron DIS";    
    modelst[21] = "Elastic";
    modelst[22] = "Moller";

    fmtst[0] = "ROOT";
    fmtst[1] = "ROOT and LUND";
    fmtst[2] = "ROOT and SOLLUND";
    int fmt = finp->Get_fmt();
    int model = fmodel->GetModel();

    printf("Model = %s \n", modelst[model].Data() );
    printf("Format = %s \n", fmtst[fmt].Data() );

    if (model == 20) { // electron DIS generator
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    else if (model == 21) { // electron ELASTIC generator
      printf("Elastic cross section \n");
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent3(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    else if (model == 22) { // electron Moller generator
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent4(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    else if(0<=model && model <=8 ) { // standard hadron generator
      printf("hadron cross section with weight factor\n");
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent2(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    else if(10<=model && model <=18 ) { // Uniform hadron generator
      printf("hadron cross section without weight factor\n");	      
      for( evt = 0; evt < nevt; evt++ ){
	  if( (evt%nprnt) == 0 ){printf("Event %10d \n", evt);}
	  fphy->MakeEvent5(fbeam, fion, fevt, fmodel);
	  fout->Write(fevt);
      }
    }
    else {
      printf("Unsupported model, check input file please.\n");	
      exit(-1);
    }
    fout->Close();

    if (fmt == 1 ) { // also LUND format 
      printf("Creating LUND output \n") ;
      fout->MakeFileLUND();
    }

    if (fmt == 2 ) { // also SOLLUND format
      printf("Creating SOLLUND output \n") ;
      fout->MakeFileSOLLUND();
    } 
    

    return;
}

eicProcess::~eicProcess(){
    return;
}
