#include "eicModel.h"

eicModel::eicModel(eicInput *inp, int whichmodel):model(0){
  double dens,X0=0.;
  int Z,N;
    if( inp ){
	model  = whichmodel;


	Z = inp->Get_ionZ();
	N = inp->Get_ionN();
	dens = inp->Get_dens();
	length = inp->Get_length();
	len_x = inp->Get_lenx();
	len_y = inp->Get_leny();

	offset.SetXYZ( inp->GetTgtXoff(),
		       inp->GetTgtYoff(),
		       inp->GetTgtZoff() );

// 	if (dens == 0.) {
// 	  if (Z==1 && N==0) {// hydrogen target
//             printf("Density Hydrogen target \n");
// 	    X0 = 890.4;
// 	  }
// 	  else if (Z==1 && N==1) {// deuterium target
// 	    printf("Density Deuterium target \n");
// 	    X0 = 745.4;
// 	  }
// 	  else {
// 	    printf("Target not in this database, update density in the input file!\n"); exit(1);
// 	  }
// 	}
// 	else {
// 	  X0 = 716.4 * double(Z+N) / ( double(Z) * double(Z+1) * log(287/pow(double(Z),0.5))) / dens;
// 	}
// 	if (length == 0.) {
// 	  if (Z==1 && N==0) {// hydrogen target
//             printf("Length Hydrogen target \n");
// 	    length = 40.;
// 	  }
// 	  else if (Z==1 && N==1) {// deuterium target
// 	    printf("Length Deuterium target \n");
// 	    length = 40.;
// 	  }
// 	  else {
// 	    printf("Target not in this database, update length in the input file!\n"); exit(1);
// 	  }
// 	} 
// 	if (len_x < 0. || len_y < 0. ) {
// 	  printf("Raster lx or ly not in this database, update radius in the input file!\n"); exit(1);
// 	}
// 
// 	if (length> 0 && X0> 0) {
// 	  radlen = length / X0;
// 	  printf("Radiation length = %e \n",radlen);
// 	}
// 	else {
// 	  printf("Radiation length less than 0, update for radiation lenght in the input file!\n"); exit(1);
// 	}
	
	if (len_x == 0. && len_y == 0. && dens== 0. && length==0. ) {
	  radlen=0;
	  printf("target density,length,ras_x and ras_y are all 0, it must be a collider!\n"); 	  
	}
	else if (dens > 0 && length > 0){
 	  X0 = 716.4 * double(Z+N) / ( double(Z) * double(Z+1) * log(287/pow(double(Z),0.5))) / dens;	  
	  radlen = length / X0;
	  printf("target has radlen %f\n",radlen);
	}
	else {
	  printf("Something wrong with the input!\n");
	  exit(1);
	}
	
    }

    length = length / 100. ; // conversion distances in m
    len_x = len_x / 100. ; 
    len_y = len_y / 100. ; 
    offset *= 1e-2; // conversion distances in cm to m

    return;
}

eicModel::~eicModel(){
    return;
}
