#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TVector3.h"
#include "snprintf.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <math.h>

using namespace std;

void Print_Usage();
void Parse_Args(int *argc, char **argv);

extern TSystem *gSystem;


//////////////////////////////////////////////////////////////
// Options Globals
//////////////////////////////////////////////////////////////
char  outfilename[128];
char  infilename[100];
double theta_min;
double theta_max;
double phi_min;
double phi_max;
int lund;
int geantino;


int main (int argc, char **argv) {

  if(argc < 5){
    Print_Usage();
    exit(1);
  }

  theta_min = -1000;
  theta_max = -1000;
  phi_min = -1000;
  phi_max = -1000;
  lund = -1000;
  geantino = -1000;

  Parse_Args( &argc, argv );
 
  if( outfilename[0] == 0 ) {
    cout << "Specify Output File Name" << endl;
    Print_Usage();
    exit(1);
  }
  
  if( infilename[0] == 0 ) {
    cout << "Specify Input File Name" << endl;
    Print_Usage();
    exit(1);
  }  

  if( theta_min == -1000 ) {
    theta_min = 0;
  }  
  cout << "Minimum of theta set to " << theta_min << endl; 

  if( theta_max == -1000 ) {
    theta_max = TMath::Pi();
  }  
  cout << "Maximum of theta set to " << theta_max << endl; 

  if( phi_min == -1000 ) {
    phi_min = -TMath::Pi();
  }  
  cout << "Minimum of phi set to " << phi_min << endl; 

  if( phi_max == -1000 ) {
    phi_max = TMath::Pi();
  } 
  cout << "Maximum of phi set to " << phi_max << endl; 
  if( lund == -1000 || lund == 0) {
    cout << "Using root format as input file" << endl;
    lund = 0;
  } 
  else cout << "Using Lund format as input file" << endl;
  if( geantino == -1000 || geantino == 0) {
    geantino = 0;
  } 
  else cout << "Using geantino/charged geantino in output" << endl;

  
  if (phi_min >= phi_max) {
    cout << "phi minimum >= phi maximum. Please recheck your values" << endl;
    Print_Usage();
    exit(1);
  }  
    
  if (theta_min >= theta_max) {
    cout << "theta minimum >= theta maximum. Please recheck your values" << endl;
    Print_Usage();
    exit(1);
  }

  double weight;
  double theta;
  double phi;
  double Ef;
  double x;
  double W;
  double y;
  double crs;
  double f1;
  double f2;
  double f1gz;
  double f3gz;
  double g1gz;
  double g5gz;
  double Q2;
  double eta_gZ;
  double rate;
  double pdf[6];
  double Dpdf[5];
  
  double Abeam;
  double Y;
  
  double A_L;
  double Y_L;
  double eta_L;
  
  double A_T;
  double Y_T;
  double eta_T;

  int particle_id;
  int charge;
  double mass;


  // Variables used for determine the lund output
  
  double nu; // energy of the Gamma*
  double MASS_p = 0.938; // mass proton in GeV 
  double px,py,pz ;
  double vx,vy,vz ;

  if (geantino == 0) {
    particle_id = 11; // Setting to have an electron or a geantino in the output
    charge = -1;
    mass = 0.00051;
  } 
  else if (geantino == 1) {
    particle_id = 0;
    charge = 0;
    mass = 0.0;
  }
  else if (geantino == 2) {
    particle_id = 0;
    charge = 1;
    mass = 0.0;
  }  
  if (lund == 0) {

    TChain input_chain("T");
    input_chain.Add(infilename);

    input_chain.SetBranchAddress("weight",&weight);
    input_chain.SetBranchAddress("theta",&theta);
    input_chain.SetBranchAddress("phi",&phi);
    input_chain.SetBranchAddress("Ef",&Ef);
    input_chain.SetBranchAddress("x",&x);
    input_chain.SetBranchAddress("W",&W);
    input_chain.SetBranchAddress("y",&y);
    input_chain.SetBranchAddress("crs",&crs);
    input_chain.SetBranchAddress("F1",&f1);
    input_chain.SetBranchAddress("F2",&f2);
    input_chain.SetBranchAddress("F1gZ",&f1gz);
    input_chain.SetBranchAddress("F3gZ",&f3gz);
    input_chain.SetBranchAddress("g1gZ",&g1gz);
    input_chain.SetBranchAddress("g5gZ",&g5gz);
    input_chain.SetBranchAddress("Q2",&Q2);
    input_chain.SetBranchAddress("eta_gZ",&eta_gZ);
    input_chain.SetBranchAddress("rate",&rate);
    input_chain.SetBranchAddress("pdf",&pdf);
    input_chain.SetBranchAddress("Dpdf",&Dpdf);
  
    input_chain.SetBranchAddress("Abeam",&Abeam);
    input_chain.SetBranchAddress("Y",&Y);
  
    input_chain.SetBranchAddress("A_L",&A_L);
    input_chain.SetBranchAddress("Y_L",&Y_L);
    input_chain.SetBranchAddress("eta_L",&eta_L);
  
    input_chain.SetBranchAddress("A_T",&A_T);
    input_chain.SetBranchAddress("Y_T",&Y_T);
    input_chain.SetBranchAddress("eta_T",&eta_T);  

    Int_t nentries = (Int_t)input_chain.GetEntries();
  
  
    ofstream OUT (outfilename);
    for (int i=0; i<nentries ; i++) {
      input_chain.GetEntry(i);
      if(i % 100000 == 0 ){
	printf("Analyzed %09d events of total %09d \n",i,nentries);
      }

      if (x>0)   nu = Q2 / (2 * MASS_p * x) ;
      else nu = 0;

      px = Ef * sin(theta) * cos(phi);
      py = Ef * sin(theta)* sin(phi);
      pz = Ef * cos(theta);

      vx = 0.0;
      vy = 0.0;
      vz = 0.0;

  
      if (phi >=phi_min && phi <=phi_max && theta >= theta_min && theta <= theta_max) {  
	
#ifdef USE_SOLLUND
	OUT << "1" << " \t " << "2"  << " \t " << "1"  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
#else
	OUT << "1" << " \t " << weight <<  " \t " << "2"  << " \t " << "1"  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
#endif
	OUT << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ef << " \t " << mass << " \t " << vx  << " \t " << vy << " \t " << vz << endl; 

      }
    }

    OUT.close();
  
  }
  
  else {  // Filtering events from Lund format
    ifstream IN (infilename);
    ofstream OUT (outfilename);
    char tmp[10];
    TVector3 v_p3;

    if (IN.is_open())
      {
	while ( IN.good() )
	  {
	    IN >> tmp  >> tmp   >> tmp   >> tmp   >> tmp   >> x  >> y   >> W   >> Q2   >> nu ;
	    IN  >> tmp  >> tmp  >> tmp  >> tmp  >> tmp  >> tmp  >> px  >> py  >> pz  >> Ef  >> tmp  >> vx   >> vy  >> vz;  
	    v_p3.SetXYZ(px,py,pz);
	    theta = v_p3.Theta();
	    phi = v_p3.Phi();
	    if (phi >=phi_min && phi <=phi_max && theta >= theta_min && theta <= theta_max) {
	      OUT << "1" << " \t " << "2"  << " \t " << "1"  << " \t " << "0"  << " \t " << "0" << " \t "  << x << " \t " << y  << " \t " << W  << " \t " << Q2  << " \t " << nu << endl;
	      OUT << " \t " << "1" << " \t " << charge << " \t " << "1" << " \t " << particle_id << " \t " << "0" << " \t " << "0" << " \t " << px << " \t " << py << " \t " << pz << " \t " << Ef << " \t " << mass << " \t " << vx  << " \t " << vy << " \t " << vz << endl;
	    }
	  }
	IN.close();
	OUT.close();
      }



  }
}




//***********************************************************************************
//  Parse arguments 
//***********************************************************************************
void Parse_Args(int *argc, char **argv){
  // Read and parse all the options, leaving only input files in the
  // argv array.
  // This version, all options flags are globals.
#define REMOVE_ONE {(*argc)--;for(j=i;j<(*argc);j++)argv[j]=argv[j+1];i--;}
#define I_PLUS_PLUS if((i+1)<(*argc)){i++;}else{break;}
  int  i,j;


  for( i=1; i<(*argc); i++ )  
    {
      if( argv[i][0] == '-' ) {
	if(strcmp(argv[i],"-o")==0 || strcmp(argv[i],"-outfile")==0)
	  {
	    I_PLUS_PLUS;
	    strcpy(outfilename,argv[i]);
	    REMOVE_ONE;
	  }
	else if(strcmp(argv[i],"-i")==0 || strcmp(argv[i],"-infile")==0)
	  {
	    I_PLUS_PLUS;
	    strcpy(infilename,argv[i]);
	    REMOVE_ONE;
	  }
	else if(strcmp(argv[i],"-th_min")==0)
          {
            I_PLUS_PLUS;
            sscanf(argv[i],"%lf",&theta_min);
            REMOVE_ONE;
          }
	else if(strcmp(argv[i],"-th_max")==0)
          {
            I_PLUS_PLUS;
            sscanf(argv[i],"%lf",&theta_max);
            REMOVE_ONE;
          }
	else if(strcmp(argv[i],"-ph_min")==0)
          {
            I_PLUS_PLUS;
            sscanf(argv[i],"%lf",&phi_min);
            REMOVE_ONE;
          }
	else if(strcmp(argv[i],"-ph_max")==0)
          {
            I_PLUS_PLUS;
            sscanf(argv[i],"%lf",&phi_max);
            REMOVE_ONE;
          }
	else if(strcmp(argv[i],"-l")==0)
          {
            I_PLUS_PLUS;
            sscanf(argv[i],"%d",&lund);
            REMOVE_ONE;
          }
	else if(strcmp(argv[i],"-g")==0)
          {
            I_PLUS_PLUS;
            sscanf(argv[i],"%d",&geantino);
            REMOVE_ONE;
          }
	else if(strcmp(argv[i],"-help")==0||strcmp(argv[i],"-h")==0)  
	  {
	    Print_Usage();
	    exit(1);
	  }
	else
	  {
	    fprintf(stderr,"\nI did not understand the option : %s\n",argv[i]);
	    Print_Usage();
	    exit(1);
	  }
	/* KILL the option from list */
	REMOVE_ONE;
      }
    }
}

//***********************************************************************************
//  Print Usage 
//***********************************************************************************

void Print_Usage() {
  cout << " root2lund : This Program translate the output from eicRate in LUND format and Filter the scattered electron angles\n";  
  cout << " Usage: root2lund -o outputfile -i inputfile [-th_min theta_min] [-th_max theta_max] [-ph_min phi_min] [-ph_max phi_max] [-l 1-0] [-g 1-0]\n";  
  cout << "     -o outputfile     : output file name (example pluto.lund)  \n";  
  cout << "     -i inputfile      : input file name (example pluto.root)  \n";  
  cout << "     -h help, print this message \n";
  cout << "   [optional] \n";
  cout << "     -th_min theta_min : specify theta minimum for scattered electron \n";
  cout << "     -th_max theta_max : specify theta maximum for scattered electron \n"; 
  cout << "     -ph_min phi_min   : specify phi minimum for scattered electron \n";
  cout << "     -ph_max phi_max   : specify phi maximum for scattered electron \n";  
  cout << "     -l 1-0            : specify format input file ( lund-> 1 (just angle filtering); root->0 (default) )\n";
  cout << "     -g 2-1-0          : specify if particle in output is a geantino ( geantino-> 1 ; charged geantino-> 2 ; normal->0 (default) )\n";
  cout << " --- \n\n";
}
