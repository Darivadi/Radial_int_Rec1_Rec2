/****************************************************************************************************
NAME: rand_rays_coordinates
FUNCTION: Generates random spherical coordinates for each ray
INPUT: 
RETURN: 0
****************************************************************************************************/
int rand_rays_coordinates(void)
{
  int m;
  const gsl_rng_type * T; /*Define el tipo de generador de números aleatorios. No hay que liberarlo*/
  gsl_rng * r; /*Análogo al w. Puntero que contiene la info sobre cual generador se va a usar,cantidad de memoria a usar, etc.*/
  long seed;

  /*+++++ Initializing random generation of numbers +++++*/
  gsl_rng_env_setup();//Inicializa las rutinas de generación                                                    
  T = gsl_rng_default;/*Inicialización de T con esta variable de GSL que es la default*/
  r = gsl_rng_alloc(T);/*Alocación de memoria*/
  
  seed = time(NULL)*getpid();
  
  gsl_rng_set(r, seed);/*Recibe puntero de inicialización de generación y  un entero largo como semilla*/
  
  for(m=0; m<GV.NRays; m++)
    {     
      /*+++++ ray for reconstruction 1 +++++*/
#if defined(COMPLETERAY_REC1) || defined(FROM90MPC_REC1)
      ray_r1[m].rad   = 0.5*GV.BoxSize;
#endif

#ifdef UNTIL90MPC_REC1 
      ray_r1[m].rad   = 90.0;
#endif

      ray_r1[m].theta = acos(1.0 - 2.0 * gsl_rng_uniform (r)); //Polar
      ray_r1[m].phi   = 2.0*M_PI * gsl_rng_uniform (r); //Azimuth

      /*+++++ ray for reconstruction 2 +++++*/
#if defined(COMPLETERAY) || defined(FROM90MPC)
      ray_r2[m].rad   = 0.5*GV.BoxSize;
#endif

#ifdef UNTIL90MPC
      ray_r2[m].rad   = 90.0;
#endif

      ray_r2[m].theta = ray_r1[m].theta; //Polar
      ray_r2[m].phi   = ray_r1[m].phi; //Azimuth

      
      /*----- Initial position's vector -----*/
      
      ray_r1[m].vec_ini[X] = 0.0;
      ray_r1[m].vec_ini[Y] = 0.0;
      ray_r1[m].vec_ini[Z] = 0.0;

      ray_r2[m].vec_ini[X] = 0.0;
      ray_r2[m].vec_ini[Y] = 0.0;
      ray_r2[m].vec_ini[Z] = 0.0;

      /*----- Vector's Final position -----*/
      ray_r1[m].vec_end[X] = ray_r1[m].rad * sin(ray_r1[m].theta) * cos(ray_r1[m].phi);
      ray_r1[m].vec_end[Y] = ray_r1[m].rad * sin(ray_r1[m].theta) * sin(ray_r1[m].phi);
      ray_r1[m].vec_end[Z] = ray_r1[m].rad * cos(ray_r1[m].theta);

      ray_r2[m].vec_end[X] = ray_r2[m].rad * sin(ray_r2[m].theta) * cos(ray_r2[m].phi);
      ray_r2[m].vec_end[Y] = ray_r2[m].rad * sin(ray_r2[m].theta) * sin(ray_r2[m].phi);
      ray_r2[m].vec_end[Z] = ray_r2[m].rad * cos(ray_r2[m].theta);

      
      if(m<2 || m==GV.NRays-1)
	{
	  printf("Ray m=%d generated at (theta, phi) = (%10.5lf, %10.5lf)\n", 
		 m, ray_r1[m].theta, ray_r1[m].phi);
	  printf("This corresponds to (x,y,z) = (%10.5lf, %10.5lf, %10.5lf)\n", 
		 ray_r1[m].vec_end[X], ray_r1[m].vec_end[Y], ray_r1[m].vec_end[Z]);
	}//if
      
    }//for m
  
  gsl_rng_free(r);

  printf("Rays generated\n");
  printf("--------------------------------------------------------------------------------------\n");
 
  return 0;
}//rand_rays_coordinates


/****************************************************************************************************
NAME: intersect_integ(void)
FUNCTION: Finds the intersections of the rays with the grid and performs the integration as 
INPUT: 
RETURN: 0
****************************************************************************************************/
int intersect_integ(struct radial_rays *ray, struct grid *gp)
{
  int i, j, k, n, m, p, q, r, s;
  double x0, y0, z0, xf, yf, zf, rad_max;
  double tMax_x, tMax_y, tMax_z, tMin_all; 
  double* dist_trav=NULL, *PotDot=NULL;

  printf("Let's find intersections\n");
  printf("--------------------------------------------------------------------------------------\n");

  dist_trav = (double *) calloc((size_t) (GV.NCELLS), sizeof(double) );
  PotDot    = (double *) calloc((size_t) (GV.NCELLS), sizeof(double) );

  printf("Memory allocated for dist_trav and PotDot arrays\n");
  printf("--------------------------------------------------------------------------------------\n");
  
  
  for(m=0; m<GV.NRays; m++)
    {      
      rad_max = 0.0;
      x0 = y0 = z0 = xf = yf = zf = 0.0;      
              
      for(p=0; p<GV.NCELLS; p++)
	{
	  dist_trav[p] = 0.0;
	  PotDot[p]    = 0.0;
	}//for p            

      p = 1;
      
      if(m%100000==0)
	printf("Entering to the algorithm of parameter t of line for m=%d with p=%d\n", m, p);

      do 
	{
	  x0 = xf;
	  y0 = yf;
	  z0 = zf;	 
	  
	  /*+++++ Computing parameter t in equation \vec[u] + t * \vec[v] = vec_end +++++*/
	  if( (ray[m].vec_end[X] - ray[m].vec_ini[X]) > 0.0 || (ray[m].vec_end[X] - ray[m].vec_ini[X]) < 0.0 )
	    tMax_x = fabs( (q*GV.CellSize - ray[m].vec_ini[X]) / (ray[m].vec_end[X] - ray[m].vec_ini[X]) );
	  else
	    tMax_x = ray[m].rad;
	  
	  if( (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) > 0.0 || (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) < 0.0 )
	    tMax_y = fabs( (r*GV.CellSize - ray[m].vec_ini[Y]) / (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) );
	  else
	    tMax_y = ray[m].rad;
	  
	  if( (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) > 0.0 || (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) < 0.0 )
	    tMax_z = fabs( (s*GV.CellSize - ray[m].vec_ini[Z]) / (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) );
	  else
	    tMax_z = ray[m].rad;

	  
	  if( m%100000==0 && p<2)
	    {
	      printf("-------------------------------------------------------------------------------\n");
	      printf("tx, ty, tz = %10.5lf, %10.5lf, %10.5lf\n", tMax_x, tMax_y, tMax_z);
	    }//if
	  
	  /*+++++ Finding the minimum of the 3 t parameters +++++*/
	  if (tMax_x < tMax_y)
	    {
	      tMin_all = tMax_x;
	      
	      if (tMax_z < tMin_all)
		{
		  tMin_all = tMax_z;	      
		}//if 1
	    }//if
	  else
	    {
	      tMin_all = tMax_y;
	      
	      if (tMax_z < tMin_all)
		{
		  tMin_all = tMax_z;
		}//if 1
	    }//else
	  
	  /*----- Depending on the minimum, the iterator per axis changes -----*/
          if(fabs(tMin_all - tMax_x) <= GV.ZERO)
            q += 1;
          else if(fabs(tMin_all - tMax_y) <= GV.ZERO)
            r += 1;
          else if(fabs(tMin_all -tMax_z) <= GV.ZERO)
            s += 1;	  	  
	  
	  if(m%100000==0 && p<2)
	    {	      
	      printf("For m=%d, p=%d, tMin_all=%lf\n", m, p, tMin_all);
	    }//if
	  
	  xf = ray[m].vec_ini[X] + tMin_all * (ray[m].vec_end[X] - ray[m].vec_ini[X]);
	  yf = ray[m].vec_ini[Y] + tMin_all * (ray[m].vec_end[Y] - ray[m].vec_ini[Y]);
	  zf = ray[m].vec_ini[Z] + tMin_all * (ray[m].vec_end[Z] - ray[m].vec_ini[Z]);

	  
	  if(fabs(xf) > fabs(ray[m].vec_end[X]) )
	    {
	      //printf("xf=%lf > ray[m].vec_end[X]=%lf for p=%d\n", xf, ray[m].vec_end[X], p);
	      xf = ray[m].vec_end[X];
	      //break;
	    }//if

	  if(fabs(yf) > fabs(ray[m].vec_end[Y]) )
	    {
	      //printf("yf=%lf > ray[m].vec_end[Y]=%lf for p=%d\n", yf, ray[m].vec_end[Y], p);
	      yf = ray[m].vec_end[Y];
	      //break;
	    }

	  if(fabs(zf) > fabs(ray[m].vec_end[Z]) )
	    {
	      //printf("zf=%lf > ray[m].vec_end[Z]=%lf for p=%d\n", zf, ray[m].vec_end[Z], p);
	      zf = ray[m].vec_end[Z];
	      //break;
	    }
	  
	  	  
	  /*+++++ Computing index of the corresponding cell +++++*/
	  i = floor( ( (xf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  j = floor( ( (yf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  k = floor( ( (zf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  n = INDEX_C_ORDER(i,j,k);

	  if(i >= GV.NCELLS)
	    {
	      printf("Warning! i=%d >= GV.NCELLS=%d with p=%d\n", i, GV.NCELLS, p);
	      //return 1;
	    }

	  if(j >= GV.NCELLS)
	    {
	      printf("Warning! j=%d >= GV.NCELLS=%d with p=%d\n", j, GV.NCELLS, p);
	      //return 1;
	    }

	  if(k>= GV.NCELLS)
	    {
	      printf("Warning! k=%d >= GV.NCELLS=%d with p=%d\n", k, GV.NCELLS, p);
	      //return 1;
	    }

	  if(n >= GV.NTOTALCELLS)
	    {
	      printf("Warning! n=%d >= GV.NTOTALCELLS=%d with p=%d\n", i, GV.NTOTALCELLS, p);
	      //return 1;
	    }
	  
	  /*+++++ Assigning distance and PotDot to compute Integral  +++++*/
	  rad_max = sqrt( POW2(xf) + POW2(yf) + POW2(zf) );

#ifdef FROM90MPC
	  if(rad_max < 90.0)
	    {
	      dist_trav[p-1] = 0.0;
	      PotDot[p-1]    = 0.0;
	      p++;
	      continue;
	    }//if
#endif

	  dist_trav[p-1] = sqrt( POW2(xf - x0) + POW2(yf - y0) + POW2(zf - z0) );	  
	  PotDot[p-1] = gp[n].potDot;
	  
	  if(m%100000==0 && p<2)
	    {	      
	      printf("Intersection at (x,y,z)=(%10.5lf,%10.5lf,%10.5lf) at cell n=%d\n", xf, yf, zf, n);
	      printf("Distance traveled: %10.5lf and PotDot = %10.5lf\n", dist_trav[p-1], PotDot[p-1]);
	    }//if	  	  
	  	  	  

	  if(p>GV.NCELLS)
	    {
	      printf("Oops! p=%d > NCells=%d for m=%d and rad_max=%lf\n", 
		     p, GV.NCELLS, m, rad_max);
	      return 1;
	    }//if
	  p++;	  	  	  

	} //while( rad_max <= (ray[m].rad - 1.0e-10) );
      while( fabs(ray[m].rad - rad_max) >  1.0e-10 );


      
      if(m%100000==0)
	{
	  printf("For the last cell:\n");
	  printf("rad_max=%10.5lf\n", rad_max);
	  printf("Intersection at (x,y,z)=(%10.5lf,%10.5lf,%10.5lf) at cell n=%d\n", xf, yf, zf, n);
	  printf("Distance traveled: %10.5lf and PotDot = %10.5lf\n", dist_trav[p-1], PotDot[p-1]);
	}//if
      

      if(m%100000==0)
	printf("Intersections ready for m=%d with p=%d and rad_max=%lf\n", m, p, rad_max);
      
      /*+++++ Computing integral +++++*/
      
      if(m%100000==0)
	printf("Computing integral for m=%d\n", m);
      
      ray[m].ISW_temp = 0.0;
	
      for(p=0; p<(GV.NCELLS); p++)
	{
	  ray[m].ISW_temp += dist_trav[p] * PotDot[p] * GV.a_SF;
	}//for p

      
      if(m%100000==0)
	{
	  printf("For ray m=%d, ISW temperature is %10.5lf\n", m, ( 2.0*GV.CMB_T0/(POW3(GV.c_SL)) ) * ray[m].ISW_temp);
	  printf("Integral ready for m=%d\n", m);
	}//if
      
    }//for m

  free(dist_trav);
  free(PotDot);
  
  return 0;
}//intersect_integ


/****************************************************************************************************
NAME: intersect_integ(void)
FUNCTION: Finds the intersections of the rays with the grid and performs the integration as 
INPUT: 
RETURN: 0
****************************************************************************************************/
int intersect_integ_Rec1(struct radial_rays *ray, struct grid *gp)
{
  int i, j, k, n, m, p, q, r, s;
  double x0, y0, z0, xf, yf, zf, rad_max;
  double tMax_x, tMax_y, tMax_z, tMin_all; 
  double* dist_trav=NULL, *PotDot=NULL;

  printf("Let's find intersections\n");
  printf("--------------------------------------------------------------------------------------\n");

  dist_trav = (double *) calloc((size_t) (GV.NCELLS), sizeof(double) );
  PotDot    = (double *) calloc((size_t) (GV.NCELLS), sizeof(double) );

  printf("Memory allocated for dist_trav and PotDot arrays\n");
  printf("--------------------------------------------------------------------------------------\n");
  
  
  for(m=0; m<GV.NRays; m++)
    {      
      rad_max = 0.0;
      x0 = y0 = z0 = xf = yf = zf = 0.0;      
              
      for(p=0; p<GV.NCELLS; p++)
	{
	  dist_trav[p] = 0.0;
	  PotDot[p]    = 0.0;
	}//for p            

      p = 1;
      
      if(m%100000==0)
	printf("Entering to the algorithm of parameter t of line for m=%d with p=%d\n", m, p);

      do 
	{
	  x0 = xf;
	  y0 = yf;
	  z0 = zf;	 
	  
	  /*+++++ Computing parameter t in equation \vec[u] + t * \vec[v] = vec_end +++++*/
	  if( (ray[m].vec_end[X] - ray[m].vec_ini[X]) > 0.0 || (ray[m].vec_end[X] - ray[m].vec_ini[X]) < 0.0 )
	    tMax_x = fabs( (q*GV.CellSize - ray[m].vec_ini[X]) / (ray[m].vec_end[X] - ray[m].vec_ini[X]) );
	  else
	    tMax_x = ray[m].rad;
	  
	  if( (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) > 0.0 || (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) < 0.0 )
	    tMax_y = fabs( (r*GV.CellSize - ray[m].vec_ini[Y]) / (ray[m].vec_end[Y] - ray[m].vec_ini[Y]) );
	  else
	    tMax_y = ray[m].rad;
	  
	  if( (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) > 0.0 || (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) < 0.0 )
	    tMax_z = fabs( (s*GV.CellSize - ray[m].vec_ini[Z]) / (ray[m].vec_end[Z] - ray[m].vec_ini[Z]) );
	  else
	    tMax_z = ray[m].rad;

	  
	  if( m%100000==0 && p<2)
	    {
	      printf("-------------------------------------------------------------------------------\n");
	      printf("tx, ty, tz = %10.5lf, %10.5lf, %10.5lf\n", tMax_x, tMax_y, tMax_z);
	    }//if
	  
	  /*+++++ Finding the minimum of the 3 t parameters +++++*/
	  if (tMax_x < tMax_y)
	    {
	      tMin_all = tMax_x;
	      
	      if (tMax_z < tMin_all)
		{
		  tMin_all = tMax_z;	      
		}//if 1
	    }//if
	  else
	    {
	      tMin_all = tMax_y;
	      
	      if (tMax_z < tMin_all)
		{
		  tMin_all = tMax_z;
		}//if 1
	    }//else

	  if(fabs(tMin_all - tMax_x) <= GV.ZERO)
            q += 1;
          else if(fabs(tMin_all - tMax_y) <= GV.ZERO)
            r += 1;
          else if(fabs(tMin_all -tMax_z) <= GV.ZERO)
            s += 1;	  	  

	  
	  if(m%100000==0 && p<2)
	    {	      
	      printf("For m=%d, p=%d, tMin_all=%lf\n", m, p, tMin_all);
	    }//if
	  
	  xf = ray[m].vec_ini[X] + tMin_all * (ray[m].vec_end[X] - ray[m].vec_ini[X]);
	  yf = ray[m].vec_ini[Y] + tMin_all * (ray[m].vec_end[Y] - ray[m].vec_ini[Y]);
	  zf = ray[m].vec_ini[Z] + tMin_all * (ray[m].vec_end[Z] - ray[m].vec_ini[Z]);

	  
	  if(fabs(xf) > fabs(ray[m].vec_end[X]) )
	    {
	      //printf("xf=%lf > ray[m].vec_end[X]=%lf for p=%d\n", xf, ray[m].vec_end[X], p);
	      xf = ray[m].vec_end[X];
	      //break;
	    }//if

	  if(fabs(yf) > fabs(ray[m].vec_end[Y]) )
	    {
	      //printf("yf=%lf > ray[m].vec_end[Y]=%lf for p=%d\n", yf, ray[m].vec_end[Y], p);
	      yf = ray[m].vec_end[Y];
	      //break;
	    }

	  if(fabs(zf) > fabs(ray[m].vec_end[Z]) )
	    {
	      //printf("zf=%lf > ray[m].vec_end[Z]=%lf for p=%d\n", zf, ray[m].vec_end[Z], p);
	      zf = ray[m].vec_end[Z];
	      //break;
	    }
	  
	  	  
	  /*+++++ Computing index of the corresponding cell +++++*/
	  i = floor( ( (xf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  j = floor( ( (yf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  k = floor( ( (zf + 0.5*GV.BoxSize) / GV.BoxSize) * GV.NCELLS );
	  n = INDEX_C_ORDER(i,j,k);

	  if(i >= GV.NCELLS)
	    {
	      printf("Warning! i=%d >= GV.NCELLS=%d with p=%d\n", i, GV.NCELLS, p);
	      //return 1;
	    }

	  if(j >= GV.NCELLS)
	    {
	      printf("Warning! j=%d >= GV.NCELLS=%d with p=%d\n", j, GV.NCELLS, p);
	      //return 1;
	    }

	  if(k>= GV.NCELLS)
	    {
	      printf("Warning! k=%d >= GV.NCELLS=%d with p=%d\n", k, GV.NCELLS, p);
	      //return 1;
	    }

	  if(n >= GV.NTOTALCELLS)
	    {
	      printf("Warning! n=%d >= GV.NTOTALCELLS=%d with p=%d\n", i, GV.NTOTALCELLS, p);
	      //return 1;
	    }
	  
	  /*+++++ Assigning distance and PotDot to compute Integral  +++++*/
	  rad_max = sqrt( POW2(xf) + POW2(yf) + POW2(zf) );

#ifdef FROM90MPC_REC1
	  if(rad_max < 90.0)
	    {
	      dist_trav[p-1] = 0.0;
	      PotDot[p-1]    = 0.0;
	      p++;
	      continue;
	    }//if
#endif

	  dist_trav[p-1] = sqrt( POW2(xf - x0) + POW2(yf - y0) + POW2(zf - z0) );	  
	  PotDot[p-1] = gp[n].potDot;
	  
	  if(m%100000==0 && p<2)
	    {	      
	      printf("Intersection at (x,y,z)=(%10.5lf,%10.5lf,%10.5lf) at cell n=%d\n", xf, yf, zf, n);
	      printf("Distance traveled: %10.5lf and PotDot = %10.5lf\n", dist_trav[p-1], PotDot[p-1]);
	    }//if	  	  
	  	  	  

	  if(p>GV.NCELLS)
	    {
	      printf("Oops! p=%d > HalfNCells=%d\n and m=%d with rad_max=%lf", 
		     p, GV.NCELLS, m, rad_max);
	      return 1;
	    }//if
	  p++;	  	  	  

	} //while( rad_max <= (ray[m].rad - 1.0e-10) );
      while( fabs(ray[m].rad - rad_max) >  1.0e-10 );


      
      if(m%100000==0)
	{
	  printf("For the last cell:\n");
	  printf("rad_max=%10.5lf\n", rad_max);
	  printf("Intersection at (x,y,z)=(%10.5lf,%10.5lf,%10.5lf) at cell n=%d\n", xf, yf, zf, n);
	  printf("Distance traveled: %10.5lf and PotDot = %10.5lf\n", dist_trav[p-1], PotDot[p-1]);
	}//if
      

      if(m%100000==0)
	printf("Intersections ready for m=%d with p=%d and rad_max=%lf\n", m, p, rad_max);
      
      /*+++++ Computing integral +++++*/
      
      if(m%100000==0)
	printf("Computing integral for m=%d\n", m);
      
      ray[m].ISW_temp = 0.0;
	
      for(p=0; p<(GV.NCELLS); p++)
	{
	  ray[m].ISW_temp += dist_trav[p] * PotDot[p] * GV.a_SF;
	}//for p

      
      if(m%100000==0)
	{
	  printf("For ray m=%d, ISW temperature is %10.5lf\n", m, ( 2.0*GV.CMB_T0/(POW3(GV.c_SL)) ) * ray[m].ISW_temp);
	  printf("Integral ready for m=%d\n", m);
	}//if
      
    }//for m

  free(dist_trav);
  free(PotDot);
  
  return 0;
}//intersect_integ
