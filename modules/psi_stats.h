#ifndef PSI_STATS_H
#define PSI_STATS_H

#include "psi_filters.h"
#include "Integrate.h"


const number high=500.;

class psi_stats: public function_cosebis
{
	public:


		psi_stats();
		~psi_stats();
		///constructor
		psi_stats(number theta_min_arc, number theta_max_arc, number LMIN, number LMAX, 
			int LBINS, 
			string FolderName=COSEBIS_DIR "/psi_Wn/",
			string WFileName= "W_psi");
		//initialisation function, calls some of the functions below
		void initialise(number theta_min_arc, number theta_max_arc, number LMIN, number LMAX, 
			int LBINS, string WFolderName1,string WFileName1);

		void setIntegrationType(string integration_type1);

		///initialises Wn_vec
		void setWFilters(int nMaximum);

		///integrant has options for gg, gm and covariance
		number integrant(number l);
		///this needs ot be set just once for each weight function
		vector<number> determine_integration_limits(number last_x_value);
		/// value of psi for a given n.
		number value_psi(int n);
		/// calculates Psi for gg or gm
		vector<number> calPsi(string integration_type1);
		/// reads input Cl and make interpolation table
		void setCl(vector<number> log_ell,vector<number> Cl);

		///////////////////////////////////////////////////////////////////
		// Covariance section
		//read in power spectra for galaxy clustering: cl_gg, galaxy-galaxy lensing: cl_gm, cosmic shear: cl_mm
		void input_for_setCl(vector<number> log_ell, vector<vector<number> > cl_gg, vector<vector<number> > cl_gm,vector<vector<number> > cl_mm);
		
		/// set noise for covariance
		void setNoise(vector<number> sigma_e_vec1,vector<number> nBar_shear_vec,vector<number> nBar_gals_vec);
		
		vector<vector<matrix* > > determine_integration_limitsCov();
		number valueCov(int n1,int m1,vector<vector<matrix* > > integ_limits_cov);
		matrix calCov(number Area_GG,number Area_GM,number Area_cross,bool Cross=true);
		
		number delta(int i1, int j1);
		int calP(int nBins,int fbin,int sbin);
		int calP_cross(int nBins_GG,int nBins_MM,int fbin,int sbin);


	private:

		number thetamin,thetamax,thetaBar,deltaTheta,lmin,lmax,theta_eval,lthresh;
		int lbins,thetaBins,nMaximum,nW,mW,integ_step,rp;
		bool find_min_max,WnSet;
		int rp1,rp2,rp3,rp4;
		int nPair_MM, nPair_GG, nPair_GM;
		int nBins_MM, nBins_GG;
		number delta1noise,delta2noise,delta3noise,delta4noise;
		string cov_case;
		string integration_type;

		//Create a vector of class objects
		vector<psi_filters> Wn_vec;
		function_cosebis Cl;

		//these are for covaraince
		vector<function_cosebis> Cl_GG_vec;
		vector<function_cosebis> Cl_GM_vec;
		vector<function_cosebis> Cl_MM_vec;

		vector<number> noise_mm_vec,noise_gg_vec;
		vector<vector<number> > integ_limits_vec_vec;
		string WFileName,WFolderName;

};

#endif

