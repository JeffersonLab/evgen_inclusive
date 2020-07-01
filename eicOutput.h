#ifndef __eicOutput_h
#define __eicOutput_h

#include "eicInput.h"
#include "eicEvent.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std; 


class eicOutput {
    public:
        eicOutput(eicInput *inp = 0, char *fileout="output.root", double nevent=-1,int whichmodel=-1);
	~eicOutput();

	void Write( eicEvent *ev );
	void Close();
	void MakeFileLUND();
    private:
	void InitTree();

	char fOutName[255];

	TFile *fOutfile;
	TTree *fTree;

	double fRunTime;
	double fNevt;
	
	int fModel;

	eventdata fData;

	eicInput *fInput;


};
#endif//__eicOutput_h
