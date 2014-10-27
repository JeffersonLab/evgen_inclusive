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
    if (nevent == -1)  {
      printf ("I don't know how many events to run\n"); 
      Print_Usage();
      exit(1);      
    } 
    if (whichmodel == -1)  {
      printf ("I don't know which model to run\n"); 
      Print_Usage();
      exit(1);      
    }     
 
    proc = new eicProcess(input_gen_file,output_gen_file,nevent,whichmodel);
 
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
	else if(strcmp(argv[i],"-n")==0||strcmp(argv[i],"-nevent")==0)  
	  {
	    I_PLUS_PLUS;
	    nevent = atof(argv[i]);
	    REMOVE_ONE;
	  }	  
	else if(strcmp(argv[i],"-m")==0||strcmp(argv[i],"-model")==0)  
	  {
	    I_PLUS_PLUS;
	    whichmodel = atoi(argv[i]);
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
  printf (" Usage: eicRate -i inputfile -o outputfile -n nevent -m model\n");  
  printf ("     -i (or -input) inputfile    input file name with correct flags (example input.dat)  \n");  
  printf ("     -o (or -output) outputfile output file name with correct flags (example output.root)  \n");
  printf ("     -n (or -nevent) number of events to run (example 1000)  \n");
  printf ("     -m (or -model) which model to run (example 0) as listed below\n");
  printf ("		0  pi+		(uniformly distributed with weight)\n");
  printf ("		1  pi-		(uniformly distributed with weight)\n");
  printf ("		2  pi0		(uniformly distributed with weight)\n");
  printf ("		3  K+		(uniformly distributed with weight)\n");
  printf ("		4  K-		(uniformly distributed with weight)\n");
  printf ("		5  Ks		(uniformly distributed with weight)\n");
  printf ("		6  Kl		(uniformly distributed with weight)\n");
  printf ("		7  p		(uniformly distributed with weight)\n");
  printf ("		8  pbar		(uniformly distributed with weight)\n");
  printf ("		10  pi+		(distributed according to crossection)\n");
  printf ("		11  pi-		(distributed according to crossection)\n");
  printf ("		12  pi0		(distributed according to crossection)\n");
  printf ("		13  K+		(distributed according to crossection)\n");
  printf ("		14  K-		(distributed according to crossection)\n");
  printf ("		15  Ks		(distributed according to crossection)\n");
  printf ("		16  Kl		(distributed according to crossection)\n");
  printf ("		17  p		(distributed according to crossection)\n");
  printf ("		18  pbar	(distributed according to crossection)\n");
  printf ("		20  e- DIS      (uniformly distributed with weight)\n");
  printf ("		21  e- Elastic  (uniformly distributed with weight)\n");
  printf ("		22  e- Moller 	(uniformly distributed with weight)\n");  
  printf ("     -h (or -help) print this message \n");
  printf (" --- \n\n");
}
