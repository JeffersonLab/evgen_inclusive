#include "string.h"
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cstring>
char input_gen_file[50];
char output_gen_file[50];
double nevent=-1;
int whichmodel=-1;
void Print_Usage();
void Parse_Args(int *argc, char **argv);
