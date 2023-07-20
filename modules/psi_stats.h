#ifndef PSI_STATS_H
#define PSI_STATS_H

#include "psi_filters.h"
// #include <gsl/gsl_sf_bessel.h>
#include "Integrate.h"


const number high=500.;

class psi_stats: public function_cosebis
{
	public:


		psi_stats();
		~psi_stats();
		///some inputs for the constructor
		psi_stats(number minTH, number maxTH, int nThetaBins, number LMIN, number LMAX, 
			int LBINS, int Nstep_Integ, number lthresh_set, int nMinimum,
			string FolderName=COSEBIS_DIR "/psi_Wn/",
			string WFileName= "W_psi");

		void initialise(number minTH, number maxTH, int nThetaBins, number LMIN, number LMAX, 
			int LBINS, int Nstep_Integ, number lthresh_set, int nMinimum,
			string WFolderName,string WFileName);

		//CosmoSIS Real Space input analysis
		//void input_for_setRealSpace(vector<number> theta, vector<vector<number> > w_CosmoSIS, vector<vector<number> > ggl_CosmoSIS);
		//void setPower_realspace(vector<number> log_theta, vector<vector<number> > w_CosmoSIS, vector<vector<number> > ggl_CosmoSIS);
		//vector<number> InterpolateRealSpace(vector<number> log_theta, vector<vector<number> > real_space_corr, vector<number> interp_theta);

		///initializes Wn_vec
		void setWFilters(int nMaximum, int nMinimum);
		//void set_theta_eval(number t_eval);
		//void set_theta_eval_rad(number t_eval);
		void set_mode(int mode_eval);
		number arcmin_to_radians(number t_arcmin); //Converts theta from arcmin --> radians
		vector<number> radians_to_arcmin_vec(vector<number> t_rad_in);
		vector<number> arcmin_to_radians_vec(vector<number> t_arcmin);
		number radians_to_arcmin(number t_rad_in);
		number integrant(number l);
		vector<number> determine_integration_limits(number last_x_value);
		number value_psi(int n, string integration_type1);
		vector<number> openFile(string FilePath,string FileNames);

		void setCl_GG(vector<number> log_ell, vector<number> Cl_GG);
		void input_for_setCl(vector<number> ell, vector<vector<number> > cl_gg, vector<vector<number> > cl_gm,vector<vector<number> > cl_mm);
		void setPower(vector<number> log_ell,vector<vector<number> > Cl_GG, vector<vector<number> > Cl_GM,vector<vector<number> > InputPower_matter);
		void setCl(vector<number> log_ell,vector<number> Cl);

		void setNoise(vector<number> sigma_e_vec1,vector<number> nBar_shear_vec,vector<number> nBar_gals_vec);
		
		vector<vector<matrix* > > determine_integration_limitsCov();
		number valueCov(int n1,int m1,vector<vector<matrix* > > integ_limits_cov);
		matrix calCov(number Area_GG,number Area_GM,number Area_cross,bool Cross=true);
		
		number delta(int i1, int j1);
		int calP(int nBins,int fbin,int sbin);
		int calP_cross(int nBins_GG,int nBins_MM,int fbin,int sbin);


	private:

		number thetaMin,thetaMax,thetaBar,deltaTheta,lmin,lmax,theta_eval,lthresh;
		int lbins,thetaBins,nMaximum,nMinimum,nW,mW,integ_step,rp;
		bool find_min_max,WnSet;
		int rp1,rp2,rp3,rp4;
		int nPair_MM, nPair_GG, nPair_GM;
		int nBins_MM, nBins_GG;
		number delta1noise,delta2noise,delta3noise,delta4noise;
		string cov_case;
		string integration_type;

		//Create a vector of class objects
		vector<psi_filters> Wn_vec;
		vector<function_cosebis> Cl_GG_vec;
		vector<function_cosebis> Cl_GM_vec;
		vector<function_cosebis> Cl_MM_vec;
		vector<number> integ_limits;
		vector<number> noise_mm_vec,noise_gg_vec;
		vector<vector<number> > integ_limits_vec_vec;
		string WFileName,WFolderName;
		function_cosebis Cl_func;


	// void test_WFilters_root_finder(int nMaximum);

};

#endif

