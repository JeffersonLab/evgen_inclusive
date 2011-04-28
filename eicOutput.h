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
        eicOutput(eicInput *inp = 0, char *fileout="output.root");
	~eicOutput();

	void Write( eicEvent *ev );
	void Close();
	void MakeFileLUND();
	void MakeFileSOLLUND();
    private:
	void InitTree();

	char fOutName[255];

	TFile *fOutfile;
	TTree *fTree;

	double fRunTime;
	double fNevt;

	eventdata fData;


};
#endif//__eicOutput_h
