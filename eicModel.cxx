#include "eicModel.h"

eicModel::eicModel(eicInput *inp):model(0){
    if( inp ){
	model  = inp->Get_model();


	printf("Model = %d \n", model);
    }

    return;
}

eicModel::~eicModel(){
    return;
}
