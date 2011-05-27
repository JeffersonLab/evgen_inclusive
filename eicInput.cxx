#include "eicInput.h"
#include <stdio.h>
#include <cstdlib>

eicInput::eicInput(const char *file){
    printf("Reading %s for input\n", file);

    FILE *f = fopen(file, "r");

    if( !f ){ printf("%s cannot be opened\n", file); exit(1); }

    char dummy[255];

    fscanf(f, "%s%d", dummy, &fData.nevt);
    fscanf(f, "%s%d", dummy, &fData.nprnt);
    fscanf(f, "%s%lf%s", dummy, &fData.lumin, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.runtime, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.e_energy, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.ion_energy, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.ion_mass, dummy);
    fscanf(f, "%s%d", dummy, &fData.ion_Z);
    fscanf(f, "%s%d", dummy, &fData.ion_N);
    fscanf(f, "%s%lf%s", dummy, &fData.dens, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.length, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.lenx, dummy);
    fscanf(f, "%s%lf%s", dummy, &fData.leny, dummy);

    fscanf(f, "%s%lf%lf%lf%s", dummy, &fData.off_x, &fData.off_y, 
	    		       &fData.off_z, dummy);
    fscanf(f, "%s%d", dummy, &fData.model);
    fscanf(f, "%s%d", dummy, &fData.out_fmt);


    fData.lumin *= 1e4; // Convert cm^-2 to m^-2
    fData.runtime *= 3600; // Convert hours to s

    fclose(f);
}

eicInput::~eicInput(){
}
