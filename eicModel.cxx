#include "eicModel.h"

eicModel::eicModel(eicInput *inp):model(0){
  double length,dens,X0=0.;
  int Z,N;
    if( inp ){
	model  = inp->Get_model();


	Z = inp->Get_ionZ();
	N = inp->Get_ionN();
	dens = inp->Get_dens();
	length = inp->Get_length();

	if (dens == 0.) {
	  if (Z==1 && N==0) {// hydrogen target
            printf("Density Hydrogen target \n");
	    X0 = 890.4;
	  }
	  else if (Z==1 && N==1) {// deuterium target
	    printf("Density Deuterium target \n");
	    X0 = 745.4;
	  }
	  else {
	    printf("Target not in this database, update density in the input file!\n"); exit(1);
	  }
	}
	else {
	  X0 = 716.4 * double(Z+N) / ( double(Z) * double(Z+1) * log(287/pow(double(Z),0.5))) / dens;
	}
	if (length == 0.) {
	  if (Z==1 && N==0) {// hydrogen target
            printf("Length Hydrogen target \n");
	    length = 40.;
	  }
	  else if (Z==1 && N==1) {// deuterium target
	    printf("Length Deuterium target \n");
	    length = 40.;
	  }
	  else {
	    printf("Target not in this database, update length in the input file!\n"); exit(1);
	  }
	}    

	if (length> 0 && X0> 0) {
	  radlen = length / X0;
	  printf("Radiation length = %e \n",radlen);
	}
	else {
	  printf("Radiation length less than 0, update for radiation lenght in the input file!\n"); exit(1);
	}
    }

    return;
}

eicModel::~eicModel(){
    return;
}
