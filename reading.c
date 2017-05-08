/****************************************************************************************************
NAME: conf2dump
FUNCTION: Reads the input file with parameters
INPUT: Parameters file
RETURN: 0
****************************************************************************************************/

int conf2dump( char filename[] )
{
  int nread;
  char cmd[1000];
  /*
  sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | gawk -F\"=\" '{print $2}' > %s.dump", 
	   filename, filename );
  */
  sprintf( cmd, "grep -v \"#\" %s | grep -v \"^$\" | awk -F\"=\" '{print $2}' > %s.dump", 
	   filename, filename );
  nread = system( cmd );
  return 0;
}


/****************************************************************************************************
NAME: read_parameters
FUNCTION: Reads the parameters
INPUT: Parameters file
RETURN: 0
****************************************************************************************************/
int read_parameters( char filename[] )
{
  int nread;
  char cmd[1000], filenamedump[1000];
  FILE *file;
  
  /*+++++ Loading the file +++++*/
  file = fopen( filename, "r" );
  if( file==NULL )
    {
      printf( "  * The file '%s' doesn't exist!\n", filename );
      return 1;
    }
  fclose(file);
  
  /*+++++ Converting to plain text +++++*/
  conf2dump( filename );
  sprintf( filenamedump, "%s.dump", filename );
  file = fopen( filenamedump, "r" );
  
  /*+++++ Parameters for binary data +++++*/
  nread = fscanf(file, "%d", &GV.NCELLS);
  nread = fscanf(file, "%s", GV.FILENAMEREC1);
  nread = fscanf(file, "%s", GV.FILENAMEREC2);  
  nread = fscanf(file, "%d", &GV.NRays);

  printf("Number of cells %10d\n", GV.NCELLS);
  printf("Number of rays %10d\n", GV.NRays);  
  printf("Data file at %s\n", GV.FILENAMEREC1);
  printf("Data file at %s\n", GV.FILENAMEREC2);

  
  fclose( file );
  
  printf( "  * The file '%s' has been loaded!\n", filename );
  
  sprintf( cmd, "rm -rf %s.dump", filename );
  nread = system( cmd );
  
  return 0;
}



/**************************************************************************************************** 
NAME: read_binary
FUNCTION: Reads the binary data file
INPUT: None
RETURN: 0 
****************************************************************************************************/
int read_binary(void)
{
  int i, nread;
  double pos_aux[3], dummy;
  FILE *inFile=NULL;

  /****** Reconstruction 1 ******/
  inFile = fopen(GV.FILENAMEREC1, "r");

  printf("Reading simulation parameters for reconstruction 1\n");

  
  /*+++++ Saving Simulation parameters +++++*/
  nread = fread(&GV.BoxSize,  sizeof(double), 1, inFile);  //Box Size
  nread = fread(&GV.Omega_M0, sizeof(double), 1, inFile);  //Matter density parameter
  nread = fread(&GV.Omega_L0, sizeof(double), 1, inFile);  //Cosmological constant density parameter
  nread = fread(&GV.z_RS,     sizeof(double), 1, inFile);  //Redshift
  nread = fread(&GV.H0,       sizeof(double), 1, inFile);  //Hubble parameter
  nread = fread(&GV.NCELLS,   sizeof(int),    1, inFile);  //Number or cells

  GV.a_SF = 1.0 / (1.0 + GV.z_RS);
  GV.NTOTALCELLS = GV.NCELLS * GV.NCELLS * GV.NCELLS;

  /*+++++ Memory allocation +++++*/
  gp1 = (struct grid *) calloc((size_t) GV.NTOTALCELLS, sizeof(struct grid));
  printf("Memory allocated!\n");
  printf("--------------------------------------------------\n");
  
  
  printf("-----------------------------------------------\n");
  printf("Cosmological parameters:\n");
  printf("OmegaM0=%lf OmegaL0=%lf redshift=%lf HubbleParam=%lf\n",
	 GV.Omega_M0,
	 GV.Omega_L0,
	 GV.z_RS,
	 GV.H0);
  printf("-----------------------------------------------\n");
  
  printf("Simulation parameters:\n");
  printf("L=%lf\n",
	 GV.BoxSize);
  printf("-----------------------------------------------\n");

      
  for(i=0; i<GV.NTOTALCELLS; i++)
    {
      /*..... File app2 .....*/
        nread = fread(&gp1[i].potDot, sizeof(double), 1, inFile);
	
	if(i%10000000==0)
	  {
	    printf("Ready for i=%d with PotDot=%lf\n", 
		   i, gp1[i].potDot);
	  }//if
    }//for m	      

  printf("Data read for reconstruction 1!\n");
  fclose(inFile);


  /****** Reconstruction 2 ******/
  inFile = fopen(GV.FILENAMEREC2, "r");

  printf("Reading reconstruction 2 parameters\n");

  
  /*+++++ Saving Simulation parameters +++++*/
  nread = fread(&dummy, sizeof(double), 1, inFile);  //Box Size
  nread = fread(&dummy, sizeof(double), 1, inFile);  //Matter density parameter
  nread = fread(&dummy, sizeof(double), 1, inFile);  //Cosmological constant density parameter
  nread = fread(&dummy, sizeof(double), 1, inFile);  //Redshift
  nread = fread(&dummy, sizeof(double), 1, inFile);  //Hubble parameter
  nread = fread(&dummy, sizeof(int),    1, inFile);  //Number or cells


  /*+++++ Memory allocation +++++*/
  gp2 = (struct grid *) calloc((size_t) GV.NTOTALCELLS, sizeof(struct grid));
  printf("Memory allocated!\n");
  printf("--------------------------------------------------\n");
  
      
  for(i=0; i<GV.NTOTALCELLS; i++)
    {
      /*..... File app2 .....*/
      nread = fread(&gp2[i].potDot, sizeof(double), 1, inFile);
      
      if(i%10000000==0)
	{
	  printf("Ready for i=%d with PotDot=%lf\n", 
		 i, gp2[i].potDot);
	}//if
    }//for m	      
  
  printf("Data read for reconstruction 2!\n");
  fclose(inFile);


  return 0;
}//read_binary


/**************************************************************************************************** 
NAME: write_binary
FUNCTION: writes file in binary format
INPUT: None
RETURN: 0 
****************************************************************************************************/
int write_binary()
{
  FILE *outFile=NULL;
  int m, i, j, k;
  
  /*+++++ Reconstruction 1 +++++*/

#ifdef COMPLETERAY_REC1
  outFile = fopen("./../../Processed_data/ISW_radial_app2_Rec1_complete.bin", "w");
#endif 

#ifdef UNTIL90MPC_REC1
  outFile = fopen("./../../Processed_data/ISW_radial_app2_Rec1_until90Mpc.bin", "w");
#endif


  fwrite(&GV.BoxSize,  sizeof(double), 1, outFile);  // Box Size                                 
  fwrite(&GV.Omega_M0, sizeof(double), 1, outFile);  // Matter density parameter
  fwrite(&GV.Omega_L0, sizeof(double), 1, outFile);  // Cosmological constant density parameter
  fwrite(&GV.z_RS,     sizeof(double), 1, outFile);  // Redshift
  fwrite(&GV.H0,       sizeof(double), 1, outFile);  // Hubble parameter
  fwrite(&GV.NCELLS,   sizeof(int),    1, outFile);  // Number of cells per axis
  fwrite(&GV.NRays,    sizeof(int),    1, outFile);  // Number of radial rays
  
  for(m=0; m<GV.NRays; m++)
    {
      fwrite(&(ray_r1[m].rad),      sizeof(double), 1, outFile);
      fwrite(&(ray_r1[m].theta),    sizeof(double), 1, outFile);
      fwrite(&(ray_r1[m].phi),      sizeof(double), 1, outFile);
      fwrite(&(ray_r1[m].ISW_temp), sizeof(double), 1, outFile);
    }//for m                                                                                            
  
  fclose(outFile);
  
  
  /*+++++ Reconstruction 2 +++++*/
#ifdef COMPLETERAY
  outFile = fopen("./../../Processed_data/ISW_radial_app2_Rec2_complete.bin", "w");
#endif  
  
#ifdef UNTIL90MPC
  outFile = fopen("./../../Processed_data/ISW_radial_app2_Rec2_until90Mpc.bin", "w");
#endif  
  
#ifdef FROM90MPC
  outFile = fopen("./../../Processed_data/ISW_radial_app2_Rec2_from90Mpc.bin", "w");
#endif  
  
  fwrite(&GV.BoxSize,  sizeof(double), 1, outFile);  // Box Size                   
  fwrite(&GV.Omega_M0, sizeof(double), 1, outFile);  // Matter density parameter   
  fwrite(&GV.Omega_L0, sizeof(double), 1, outFile);  // Cosmological constant density parameter
  fwrite(&GV.z_RS,     sizeof(double), 1, outFile);  // Redshift                                     
  fwrite(&GV.H0,       sizeof(double), 1, outFile);  // Hubble parameter                      
  fwrite(&GV.NCELLS,   sizeof(int),    1, outFile);  // Number of cells per axis
  fwrite(&GV.NRays,    sizeof(int),    1, outFile);  // Number of radial rays
  
  for(m=0; m<GV.NRays; m++)
    {
      fwrite(&(ray_r2[m].rad),      sizeof(double), 1, outFile);
      fwrite(&(ray_r2[m].theta),    sizeof(double), 1, outFile);
      fwrite(&(ray_r2[m].phi),      sizeof(double), 1, outFile);
      fwrite(&(ray_r2[m].ISW_temp), sizeof(double), 1, outFile);
    }//for m            

  fclose(outFile);

  return 0;
}//write_binary
