#include "eicModel.h"

eicModel::eicModel(eicInput *inp):model(0){
    if( inp ){
	model  = inp->Get_model();
	radlen = inp->Get_radlen();
    }

    return;
}

eicModel::~eicModel(){
    return;
}
