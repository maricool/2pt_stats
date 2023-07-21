#ifndef PSI_FILTERS_H
#define PSI_FILTERS_H
#include "function_cosebis.h"
//BesselJ0_Zeros and UFilter_roots are saved here
#include "psi_filters_integrand_zeros.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include "matrix.h"
#include "Integrate.h"

// Wgm is not correct probably, Q needs to be fixed

class psi_filters : public function_cosebis
{
	public:
		//Constructors
		psi_filters();
		~psi_filters();

		// Initialize the global parameters for the functions
		psi_filters(number thetamin, number thetamax, int nMaximum,
			number LMIN, number LMAX, int LBINS,
			string FolderName=COSEBIS_DIR "/psi_Wn/",
			string WFileName= "W_psi");

		void initialise(number thetamin, number thetamax, int nMaximum,
			number LMIN, number LMAX, int LBINS,string FolderName,string WFileName);

		//^^^^ calls the next two functions to set the parameters
		void setValues(number thetamin1, number thetamax1, number LMIN, number LMAX, int LBINS);
		// Sets the maximum mode being analysed.
		void initialize_mode(int nMaximum);
		void setNames(string FolderName1,string WFileName1);
		void set_integration_type(string integration_type1);
		//This calculates W_gg(l) for a given l
		number get( number l); 
		//void loadWZeros(string filename="../psi/zeros_j0.txt");
		number valueFuncW(number l, number x);
		//Determines the minima and maxima of a function so then can compute the guassian integration in-between pts.
		vector<number> determine_integration_limits(number last_x_value);
		vector<number> find_zeros_of_W_integrand(number last_x_value, vector<number> UFilter_roots_for_mode);
		void roots_zeros();
		//This function tests the root finder call
		void test_find_zeros_of_W_integrand(int m, number elle, int ttbins);
		int sum_of_digits(int max_digit);

		//Define the integrants used within the class
		number integrant(number l);
		// number integrant_Q(number l);
		// number integrant_Wgm(number l);
		// number integrant_Wgg(number l);


		number U(int mode, number theta);
		number U1(number theta);
		number Un(int mode, number theta);
		number radians_to_arcmin(number t_in);

		number analyticalQ(number theta, int n);

		void set(int order, string WnFileName="W_nn_file_");


		function_cosebis Wn_table;

	private:

		int mode, nMaximum, lbins; 
		number thetamin, thetamax, thetaBar, deltaTheta, lmin, lmax, lmode;
		vector<number> all_roots;

		number legendre_coeffs(number sum_idx);
		number factorial(number num_in);

		string FolderName,WFileName;
		string integration_type;

};

#endif