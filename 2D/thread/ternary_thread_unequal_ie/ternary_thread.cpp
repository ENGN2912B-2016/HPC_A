/* High Performance computing project for course ENGN 2912B
 * Code written by Arjun, Aravind and Sayan
 * 
 *
 *Input file format (each value has to be a double):
 *X (interaction parameter, default - 3.5)
 *kA kB kC (gradient energy coefficient, default - 4.0)
 *mAA mBB mAB  (mobilities, range: -1 to 1, default: 1.0)
 *cA cB cC (initial conc, range: 0 to 1 (excluding 0), sum must be 1, no default)
 *delta_x delta_y delta_t (scale, default: 1.0 1.0 0.5)
 *
 *Please note that the dafault value would be assigned only if the value is 0 or garbage (with an warning ofcourse), the input files always has to be a 5 lines file. If that is not followed, then the program throws and error and exits.
 *
 *To compile type $g++ ternary.cpp -lfftw3 -lm -std=c++11 -g
 */ 


#define _USE_MATH_DEFINES
#include <cmath> //for log and M_PI constants
#include <iostream> //for console output
#include <fftw3.h> //for discreet fourier transforms
#include <fstream> //for reading input and writing output
#include <chrono> //for the time seed
#include <random> //for the initial random function generator
#include <thread>
using namespace std;

class ternary
{
  fftw_complex *cA, *cB, *cC, *gA, *gB, *cAtemp, *cBtemp; //composition and free energy variables
  fftw_plan planFcA, planBcA, planFcB, planBcB, planFgA, planFgB, planBgA, planBgB, planFcAtemp, planFcBtemp, planBcAtemp, planBcBtemp; //fourier transform plans

  int n_x, n_y; //grid size
  double kx, ky, kx2, ky2, k2, k4; //reciprocal space vectors
  double X, kC, kAA, kBB, mAA, mBB, mAB; //physical constants
  double delta_x, delta_y, delta_kx, delta_ky; //spatial step
  double delta_t; //time step

public:

  //ternary();//default constructor
  ternary(double X1, double kA, double kB, double kC1, double icA, double icB, double icC, double mAA, double mBB, double mAB, double dx, double dy, double dt);//parameterized constructor

  //eint writeData(fstream *f); //method to write to file

  void simulate(); //method which runs the simulation
};

ternary::ternary(double X1, double kA, double kB, double kC1, double icA, double icB, double icC, double moAA, double moBB, double moAB, double dx, double dy, double dt)
  {
    X = X1;
    kAA = kA + kC1;
    kBB = kB + kC1;
    kC = kC1;
    n_x = n_y = 512;
    delta_x = dx;
    delta_y = dy;
    delta_t = dt;
    kx = ky = 0.0;
    kx2 = 0.0;
    ky2 = 0.0;
    k2 = 0.0;
    k4 = 0.0;
    mAA = moAA;
    mBB = moBB;
    mAB = moAB;

    //setting up the random number engine to introduce the initial perturbation in concentration
    //chrono::high_resolution_clock myclock;
    //myclock::duration d = myclock::now();
    //unsigned seed = d.count();
    mt19937_64 generator(time(0));
    uniform_real_distribution<double> distribution(0.0,1.0);

    fftw_init_threads();
    fftw_plan_with_nthreads(12);

    //declaring memory for the concentration and free energy variable
    cA = (fftw_complex*)fftw_malloc(n_x*n_y*sizeof(fftw_complex));
    cB = (fftw_complex*)fftw_malloc(n_x*n_y*sizeof(fftw_complex));
    cC = (fftw_complex*)fftw_malloc(n_x*n_y*sizeof(fftw_complex));
    gA = (fftw_complex*)fftw_malloc(n_x*n_y*sizeof(fftw_complex));
    gB = (fftw_complex*)fftw_malloc(n_x*n_y*sizeof(fftw_complex));
    cAtemp = (fftw_complex*)fftw_malloc(n_x*n_y*sizeof(fftw_complex));
    cBtemp = (fftw_complex*)fftw_malloc(n_x*n_y*sizeof(fftw_complex));

    //setting up the forward plans
    planFcA = fftw_plan_dft_2d(n_x,n_y,cA,cA,FFTW_FORWARD,FFTW_ESTIMATE);
    planFcB = fftw_plan_dft_2d(n_x,n_y,cB,cB,FFTW_FORWARD,FFTW_ESTIMATE);
    planFgA = fftw_plan_dft_2d(n_x,n_y,gA,gA,FFTW_FORWARD,FFTW_ESTIMATE);
    planFgB = fftw_plan_dft_2d(n_x,n_y,gB,gB,FFTW_FORWARD,FFTW_ESTIMATE);
    planFcAtemp = fftw_plan_dft_2d(n_x,n_y,cAtemp,cAtemp,FFTW_FORWARD,FFTW_ESTIMATE);
    planFcBtemp = fftw_plan_dft_2d(n_x,n_y,cBtemp,cBtemp,FFTW_FORWARD,FFTW_ESTIMATE);

    //setting up the inverse fourier transform plans
    planBcA = fftw_plan_dft_2d(n_x,n_y,cA,cA,FFTW_BACKWARD,FFTW_ESTIMATE);
    planBcB = fftw_plan_dft_2d(n_x,n_y,cB,cB,FFTW_BACKWARD,FFTW_ESTIMATE);
    planBgA = fftw_plan_dft_2d(n_x,n_y,gA,gA,FFTW_BACKWARD,FFTW_ESTIMATE);
    planBgB = fftw_plan_dft_2d(n_x,n_y,gB,gB,FFTW_BACKWARD,FFTW_ESTIMATE);
    planBcAtemp = fftw_plan_dft_2d(n_x,n_y,cAtemp,cAtemp,FFTW_FORWARD,FFTW_ESTIMATE);
    planBcBtemp = fftw_plan_dft_2d(n_x,n_y,cBtemp,cBtemp,FFTW_FORWARD,FFTW_ESTIMATE);

    for(int i1 = 0; i1 < n_x; ++i1)
      {
	for(int i2 = 0; i2 < n_y; ++i2)
	  {
	    //genrating the first disturbed microstructure around the initial composition
	    double u = distribution(generator);
	    cA[i2+n_y*i1][0] = icA + (icA - u)*5e-3;
	    cA[i2+n_y*i1][1] = 0.0;
	    u = distribution(generator);
	    cB[i2+n_y*i1][0] = icB + (icB - u)*5e-3;
	    cB[i2+n_y*i1][1] = 0.0;
	    
	    cC[i2+n_y*i1][0] = icC + (icC - u)*5e-3;
	    cC[i2+n_y*i1][1] = 0.0;

	    cAtemp[i2+n_y*i1][0] = 0.0;
	    cBtemp[i2+n_y*i1][0] = 0.0;
	    cAtemp[i2+n_y*i1][1] = 0.0;
	    cBtemp[i2+n_y*i1][1] = 0.0;
	  }
      }

    delta_kx = (2.0*M_PI)/(n_x*delta_x);
    delta_ky = (2.0*M_PI)/(n_y*delta_y);
  }

 void ternary::simulate()
 {
   double coe1, coe2, coe3, coe4; //coeeficient of the equation
   double det; //determinant, required to solve using cramers rule
   int half_nx = (int)(n_x/2); //required for periodic boundary condition
   int half_ny = (int)(n_y/2);
   ofstream out;
   string s;
   int i1,i2;

   //running for all time steps
   for(int INDEX = 0; INDEX <= 100000; ++INDEX)
     {
       for(i1 = 0; i1 < n_x; ++i1)
	 {
	   for(i2 = 0; i2 < n_y; ++i2)
	     {
	       //calculating the free energies
	       gA[i2+n_y*i1][0] = X*(cC[i2+n_y*i1][0] - cA[i2+n_y*i1][0]) + log(cA[i2+n_y*i1][0]) - log(cC[i2+n_y*i1][0]);
	       gB[i2+n_y*i1][0] = X*(cC[i2+n_y*i1][0] - cB[i2+n_y*i1][0]) + log(cB[i2+n_y*i1][0]) - log(cC[i2+n_y*i1][0]);

	       gA[i2+n_y*i1][1] = 0.0;
	       gB[i2+n_y*i1][1] = 0.0;
	     }
	 }

       //execute forward fourier transform
       fftw_execute_dft(planFcA,cA,cA);
       fftw_execute_dft(planFgA,gA,gA);
       fftw_execute_dft(planFcB,cB,cB);
       fftw_execute_dft(planFgB,gB,gB);
       fftw_execute_dft(planFcAtemp,cAtemp,cAtemp);
       fftw_execute_dft(planFcBtemp,cAtemp,cAtemp);

       for(i1 = 0; i1 < n_x; i1++)
	 {
	   //incorporate periodic boundary condition along xe
	   if(i1 < half_nx)
	     kx = i1*delta_kx;
	   else
	     kx = (i1 - n_x)*delta_kx;

	   kx2 = kx*kx;
	   for(i2 = 0; i2 < n_y; i2++)
	     {
	       //another similar boundary condition along y
	       if(i2 < half_ny)
		 ky = i2*delta_ky;
	       else
		 ky = (i2 - n_y)*delta_ky;

	       ky2 = ky*ky;
	       k2 = kx2 + ky2;
	       k4 = k2*k2;

	       coe1 = 1 + 2*k4*delta_t*(mAA*kAA + mAB*kC);
	       coe2 = 2*k4*delta_t*(mAA*kC + mAB*kBB);
	       coe3 = 2*k4*delta_t*(mAB*kAA + mBB*kC);
	       coe4 = 1 + 2*k4*delta_t*(mAB*kC + mBB*kBB);
	       det = coe1*coe4 - coe2*coe3;

	       //solving cA and cB simultaneously
	       cAtemp[i2 + n_y*i1][0] = (coe4*(cA[i2 + n_y*i1][0] - k2*delta_t*(mAA*gA[i2 + n_y*i1][0] + mAB*gB[i2 + n_y*i1][0])) - coe2*(cB[i2 + n_y*i1][0] - k2*delta_t*(mAB*gA[i2 + n_y*i1][0] + mBB*gB[i2 + n_y*i1][0])))/det;

	       cBtemp[i2 + n_y*i1][0] = (coe1*(cB[i2 + n_y*i1][0] - k2*delta_t*(mAB*gA[i2 + n_y*i1][0] + mBB*gB[i2 + n_y*i1][0])) - coe3*(cA[i2 + n_y*i1][0] - k2*delta_t*(mAA*gA[i2 + n_y*i1][0] + mAB*gB[i2 + n_y*i1][0])))/det;

	       cA[i2 + n_y*i1][0] = cAtemp[i2 + n_y*i1][0];
	       cB[i2 + n_y*i1][0] = cBtemp[i2 + n_y*i1][0];

	     }
	 }

       //convert back from reciprocal space -> real space by inverse fourier transform

       fftw_execute_dft(planBcA,cA,cA);
       fftw_execute_dft(planBcB,cB,cB);
       fftw_execute_dft(planBgA,gA,gA);
       fftw_execute_dft(planBgB,gB,gB);
       fftw_execute_dft(planBcAtemp,cAtemp,cAtemp);
       fftw_execute_dft(planBcBtemp,cBtemp,cBtemp);

       for(i1 = 0; i1 < n_x; i1++)
	 {
	   for(i2 = 0; i2 < n_y; i2++)
	     {
	       //normalizing the values
	       cA[i2 + n_y*i1][0] = cA[i2 + n_y*i1][0]/(n_x*n_y);
	       cB[i2 + n_y*i1][0] = cB[i2 + n_y*i1][0]/(n_x*n_y);

	       cA[i2 + n_y*i1][1] = 0.0;
	       cB[i2 + n_y*i1][1] = 0.0;

	       //calculate cC
	       cC[i2 + n_y*i1][0] = 1.0 - cA[i2 + n_y*i1][0] - cB[i2 + n_y*i1][0];
	       cC[i2 + n_y*i1][1] = 0.0;

	       //cout << cA[i2+n_y*i1][0] << " " << cB[i2+n_y*i1][0] << " " << cC[i2+n_y*i1][0] << "\n";
	     }
	 }

       //prints data to fiile for every 1000 timestep
       if(INDEX%1000 == 0)
	 {
	   s = "./output/timestep"+to_string(INDEX)+".dat";
	   out.open(s);
	   for(i1 = 0; i1 < n_x; i1++)
	     {
	       for(i2 = 0; i2 < n_y; i2++)
		 {
		   double r = 10*cA[i2+n_y*i1][0] + 30*cB[i2+n_y*i1][0] + 20*cC[i2+n_y*i1][0];
		   out <<  to_string(r) << " ";
		 }
	       out << "\n";
	     }
	   out.close();
	 }
       
     }
       //release the variables
       fftw_free(cA);
       fftw_free(cB);
       fftw_free(cC);
       fftw_free(gA);
       fftw_free(gB);
       fftw_free(cAtemp);
       fftw_free(cBtemp);

       //destroy the plans
       fftw_destroy_plan(planFcA);
       fftw_destroy_plan(planFcB);
       fftw_destroy_plan(planFgA);
       fftw_destroy_plan(planFgB);
       fftw_destroy_plan(planBcA);
       fftw_destroy_plan(planBcB);
       fftw_destroy_plan(planBgA);
       fftw_destroy_plan(planBgB);
 }

int main(int argc, char **argv)
{
  ternary *t = new ternary(3.5, atof(argv[1]), atof(argv[2]), atof(argv[3]), 0.25, 0.25, 0.5, 1.0, 1.0,-0.5, 1.0, 1.0, 0.05);

  t->simulate();
  return 1;
}
