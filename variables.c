/*********************************************************************
Index conventions:
i = x axis
j = y axis
k = z axis
l = particle
m = grid cell  
*********************************************************************/


/********************************************************************
                                   UNITS
********************************************************************
Frame's time = 1
Redshift = 0
Flagsfr = 0
Flagsfed = 0
Flagcool = 0
numfiles = 1
Box Size = 1
Matter density parameter \Omega_{m,0} = 0.258
Dark energy density parameter \Omega_{\Lambda,0} = 0.742
Hubble's parameter h = 0.72
Hubble's constant H_0 = 1
Particle's mass = 3.41454 (* 10^{10} M_{Sun}/h)
Total number of particles = 134217728
Mass unit = 1 * 10^{10} M_{Sun}/h
Lenght unit = 1 Mpc/h
Gravitational constant in the interal units G = 43.0071
*******************************************************************/


/*****************************************************************
                            STRUCTURES
******************************************************************/

struct grid
{
  double potDot;   //Potential's time derivative (exact solution)
}*gp1=NULL, *gp2=NULL; //grid



struct radial_rays
{
  /*+++++ Spherical coordinates +++++*/    
  double rad;   // Radius
  double theta; // Polar angle
  double phi;   // Azimuth angle

  /*+++++ Cartesian Coordinates +++++*/
  double vec_ini[3]; // Initial positions (all are 0.0)
  double vec_end[3]; // Final positions (maximum values are BoxSize/2)
  
  /*+++++ ISW temperature fluctuation +++++*/
  double ISW_temp;

}*ray_r1, *ray_r2, *ray_s; //radial_rays


struct GlobalVariables
{
  char FILENAMEREC1[1000]; //Path of the data file
  char FILENAMEREC2[1000]; //Path of the data file

  /*+++ Grid constants +++*/
  double BoxSize;      // Size of the simulation box in one axis (all must be the same)
  int NCELLS;       // Number of cells in one axis
  int NTOTALCELLS;  // Total number of cell
  

  double CellSize; // Size of the cell
  double ZERO;     // Zero for the computer
  int NRays;       // Number of radial rays to perform integration


  /*+++ Cosmological Parameters +++*/
  double H0;      //= 1.0 Hubble's constant in the inner units
  double z_RS;    // = 0.0 Redshift of the simulation
  double a_SF;    // Scale factor's time derivative
  double Hz;      //Hubble's parameter a_dot/a
  double Omega_M0; //= 0.258 Density parameter of matter
  double Omega_L0; //= Density parameter of cosmological constant
  double MeanDen; // MeanDens;= 7.160809 Units *1E10 M_Sun/h
  double c_SL; // Speed of light 300000 km/s
  double CMB_T0; //Mean temperature of CMB in K
}GV;//globalVariables


/***************************************************************
                       DEFINITIONS
 ***************************************************************/

#define X 0
#define Y 1
#define Z 2
#define POW2(x) ((x)*(x))
#define POW3(x) ((x)*(x)*(x))
#define INDEX_C_ORDER(i,j,k) (k)+GV.NCELLS*((j)+GV.NCELLS*(i)) //Index in C-order
#define INDEX_2D_XLOS(j,k) (k)+GV.NCELLS*(j)//Index in C-order in 2D with x as LOS
#define INDEX_2D_YLOS(i,k) (k)+GV.NCELLS*(GV.NCELLS*(i)) //Index in C-order in 2D with y as LOS
#define INDEX_2D_ZLOS(i,j) GV.NCELLS*((j)+GV.NCELLS*(i)) //Index in C-order in 2D with z as LOS

