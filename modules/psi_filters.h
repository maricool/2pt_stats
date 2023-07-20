#ifndef PSI_FILTERS_H
#define PSI_FILTERS_H
#include "function_cosebis.h"
//BesselJ0_Zeros and UFilter_roots are saved here
#include "psi_filters_integrand_zeros.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iterator>
#include <string>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include "matrix.h"
#include "Integrate.h"

//Only goes to lmax=5e5 above that everything goes to zero.
// We can switch to determin integration limit to go higher

class psi_filters : public function_cosebis
{
	public:
		//Constructors
		psi_filters();
		~psi_filters();

		// Initialize the global parameters for the functions
		psi_filters(double minTH, double maxTH, int Bins, int Nmax, double LMIN, 
			double LMAX, int LBINS, int Nmin=1,
			string FolderName="../psi/WFilters/",
			string WFileName="W_psi");

		void initialise(double minTH, double maxTH, int Bins, int Nmax,double LMIN, 
			double LMAX, int LBINS, int Nmin,string FolderName,string WFileName);

		//^^^^ calls the next two functions to set the parameters
		void setValues(double minTH, double maxTH, int Bins, double LMIN, double LMAX, int LBINS);
		// Sets the maximum mode being analysed.
		void initialize_mode( int m_max, int m_min=1);
		void setNames(string FolderName1,string WFileName1);
		//This calculates W_gg(l) for a given l
		double get( double l); 
		//void loadWZeros(string filename="../psi/zeros_j0.txt");
		double valueFuncW(double l, double x);
		//Determines the minima and maxima of a function so then can compute the guassian integration in-between pts.
		vector<double> determine_integration_limits(double last_x_value);
		vector<double> find_zeros_of_W_integrand(double last_x_value, vector<double> UFilter_roots_for_mode);
		void roots_zeros();
		//This function tests the root finder call
		void test_find_zeros_of_W_integrand(int m, double elle, int ttbins);
		int sum_of_digits(int max_digit);

		//Define the integrants used within the class
		number integrant(number l);
		// number integrant_Q(number l);
		// number integrant_Wgm(number l);
		// number integrant_Wgg(number l);


		number U(int mode, double theta);
		number U1(number theta);
		number Un(int mode, number theta);
		double radians_to_arcmin(double t_in);

		number analyticalQ(double theta, int n);

		void set(int order, string WnFileName="W_nn_file_");


		function_cosebis Wn_table;

	private:

		/* Global parameters of the class 

		+ mode = COSEBI's mode currently being analysed.
		+ Nmodes = Maximum mode analysed.
		+ thetabins = No. of thetabins (linspace) analysed by ATHENA.
		+ lbins = No. of (logspace) lbins in the integral over power spectra to obtain the FAB stats.

		+ thetaMin,thetaMax = Theta range of statistics.
		+ thetaBar = (thetaMin+thetaMax)/2.
		+ deltaTheta = thetaMax - thetaMin.
		+ lmin,lmax = ell-limits of Cl integral.
		+ lmode = current ell of which the statistic is being analysed.

		*/

		int mode, Nmodes, thetaBins, lbins; 
		double thetaMin, thetaMax, thetaBar, deltaTheta, lmin, lmax, lmode;
		vector<double> all_roots;

		number legendre_coeffs(number sum_idx);
		number factorial(number num_in);

		string FolderName,WFileName;
		string integration_type;
		//vector<number> BesselJ0_Zeros;


};

#endif