#include "eicRate.h"
#include "eicProcess.h"


int main(int argc, char **argv){
    char default_file[255] = "input.dat";
    char defaultname[255] = "output.root";
    Parse_Args( &argc, argv );
   
 
    eicProcess *proc;
    if (strlen(input_gen_file) == 0)  {
      strcpy(input_gen_file,default_file);
      printf ("No input file specified. Using default input file %s \n",default_file);  
      Print_Usage();
    }
    if (strlen(output_gen_file) == 0)  {
      strcpy(output_gen_file,defaultname);
      printf ("No output file specified. Using default output file %s \n",defaultname); 
      Print_Usage();
    }
 
 
    proc = new eicProcess(input_gen_file,output_gen_file);
 
    proc->Run();

    return 0;
}

void Parse_Args(int *argc, char **argv) {
#define REMOVE_ONE {(*argc)--;for(j=i;j<(*argc);j++)argv[j]=argv[j+1];i--;}
#define I_PLUS_PLUS if((i+1)<(*argc)){i++;}else{break;}
  int  i,j;
  for( i=1; i<(*argc); i++ )  
    {
      if( argv[i][0] == '-' ) {
	if(strcmp(argv[i],"-i")==0 || strcmp(argv[i],"-input")==0)
	  {
	    I_PLUS_PLUS;
	    strcpy(input_gen_file,argv[i]);
	    REMOVE_ONE;
	  }
	else if(strcmp(argv[i],"-o")==0||strcmp(argv[i],"-output")==0)  
	  {
	    I_PLUS_PLUS;
	    strcpy(output_gen_file,argv[i]);
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


void Print_Usage() {
  printf (" eicRate : This Program Generates events \n");  
  printf (" Usage: eicRate -i inputfile -o outputfile \n");  
  printf ("     -i (or -input) inputfile    input file name with correct flags (example input.dat)  \n");  
  printf ("     -o (or -output) outputfile output file name with correct flags (example output.root)  \n");  
  printf ("     -h (or -help) print this message \n");
  printf (" --- \n\n");
}
