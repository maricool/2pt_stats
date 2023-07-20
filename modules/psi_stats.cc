#include "psi_stats.h"


psi_stats::psi_stats(){}
psi_stats::~psi_stats(){}

psi_stats::psi_stats(number minTH, number maxTH, int nThetaBins, number LMIN, number LMAX, 
	int LBINS, int Nstep_Integ, number lthresh_set,
    int nMinimum,string WFolderName, string WFileName)
{ 
	initialise(minTH, maxTH, nThetaBins, LMIN, LMAX, 
		LBINS, Nstep_Integ, lthresh_set,nMinimum,WFolderName,WFileName);
}

void psi_stats::initialise(number minTH, number maxTH, int nThetaBins, number LMIN, number LMAX, 
			int LBINS, int Nstep_Integ, number lthresh_set, int nMinimum1,
			string WFolderName1,string WFileName1)
{
	thetaMin = minTH;
	thetaMax = maxTH;
	thetaBar = (thetaMax + thetaMin) / 2.;
	deltaTheta = thetaMax - thetaMin;
	thetaBins = nThetaBins;

	//The number of integration steps for the Cl -> psi integral
	//integ_step = Nstep_Integ;

	//Set paramters for the l-range which the power-spectra is integrated 
	// (i.e an approximation from l=0 --> l=infinity in reality, 
	// but because of the limited support of the filter functions this is not necessary)
	lmin = LMIN;
	lmax = LMAX;
	lbins = LBINS;


	//lthresh = lthresh_set;

	nMinimum = nMinimum1;

	WFolderName = WFolderName1;
	WFileName = WFileName1;


	//Set-up the boolean to find min/max here
	find_min_max = true;
	//Wn are not set yet
	WnSet =false; 
}


number psi_stats::arcmin_to_radians(number theta)
{
	return 2.*pi*theta/(60.*360.);
}

vector<number> psi_stats::arcmin_to_radians_vec(vector<number> theta_in_vec)
{
	vector<number> theta_radians_vec; 
	theta_radians_vec.clear();

	for (int element=0; element<int(theta_in_vec.size()); element++)
		theta_radians_vec.push_back(2.*pi*theta_in_vec[element]/(60.*360.));	
	
	return theta_radians_vec;

}

vector<number> psi_stats::radians_to_arcmin_vec(vector<number> t_rad_in)
{

	vector<number> t_arcmin_vec; t_arcmin_vec.clear();

	for (int element=0; element<int(t_rad_in.size()); element++)
		t_arcmin_vec.push_back(t_rad_in[element]*(60.*360.)/(2.*pi));	
	
	return t_arcmin_vec;

}

number psi_stats::radians_to_arcmin(number t_rad_in)
{
 	return t_rad_in*(60.*360.)/(2.*pi);
}

void psi_stats::set_mode(int mode_eval)
{
	nW=mode_eval-1;
}


//This method initializes the Wns
void psi_stats::setWFilters(int nMaximum1, int nMinimum)
{
	//Set-up of global parameter describing the max no. of modes analysed
	if((!WnSet) || (nMaximum!=nMaximum1))
	{
		nMaximum = nMaximum1;
		//Initiate W-Class
		psi_filters Wn(thetaMin, thetaMax, thetaBins, nMaximum, lmin, 
			lmax, lbins, nMinimum,WFolderName,WFileName);

		//Clear W_gg(n,l) vector
		Wn_vec.clear();
		
		///these need to be in separate for loops, otherwise a segmentation fault happens!
		for(int n=nMinimum; n<=nMaximum; n++)
			Wn_vec.push_back(Wn);
		for(int n=0; n<=(nMaximum-nMinimum); n++)
			Wn_vec[n].set(n+nMinimum);

		WnSet=true;
	}
}


// Need to change this into one function:

number psi_stats::integrant(number l)
{
	if(integration_type=="gg")
	{
		return Wn_vec[nW].value(l) * l * Cl_GG_vec[rp].value(l);
	}
	else if (integration_type=="gm")
	{
		return Wn_vec[nW].value(l) * l * Cl_GM_vec[rp].value(l);
	}
	else if (integration_type=="cov")
	{
		number integ = 0.0;
		if(cov_case=="clustering")
		{
			number power1=Cl_GG_vec[rp1].value(l)+delta1noise;
			number power2=Cl_GG_vec[rp2].value(l)+delta2noise;
			number power3=Cl_GG_vec[rp3].value(l)+delta3noise;
			number power4=Cl_GG_vec[rp4].value(l)+delta4noise;
			number powers=power1*power2+power3*power4;
			integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
			return integ;
		}
		else if(cov_case=="ggl")
		{
			number power1=Cl_GG_vec[rp1].value(l)+delta1noise;
			number power2=Cl_MM_vec[rp2].value(l)+delta2noise;
			number power3=Cl_GM_vec[rp3].value(l);
			number power4=Cl_GM_vec[rp4].value(l);
			number powers=power1*power2+power3*power4;
			integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
			return integ;
		}
		else if(cov_case=="cross")
		{
			number power1=Cl_GG_vec[rp1].value(l)+delta1noise;
			number power2=Cl_GM_vec[rp2].value(l);
			number power3=Cl_GM_vec[rp3].value(l);
			number power4=Cl_GG_vec[rp4].value(l)+delta4noise;
			number powers=power1*power2+power3*power4;
			integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
			return integ;
		}
		else
		{
			clog<<"not a recognised covariance case, exiting now ..."<<endl;
			exit(1);
		}	
	}
	else
	{
		clog<<"not a recognised integration case, exiting now ..."<<endl;
		exit(1);
	}
}


// // integrants start here:
// number psi_stats::integrant_Cl_gg(number l)
// {
// 	//here nW = current mode being analysed - 1 because counter starts at 0 and mode counter starts at 1
// 	//here function will do the interpolation for l values that are not in the Wn_vec[nW] table.
// 	//clog << l << " " << Wn_vec[nW].value(l) << endl;
// 	number integ = Wn_vec[nW].value(l) * l * Cl_func.value(l);
// 	return integ;
// }


// number psi_stats::integrant_Cl_gm(number l)
// {
// 	//here nW = current mode being analysed - 1 because counter starts at 0 and mode counter starts at 1
// 	//here function will do the interpolation for l values that are not in the Wn_vec[nW] table. 
// 	number integ=Wn_vec[nW].value(l) * l * Cl_func.value(l);

// 	return integ;
// }

// number psi_stats::integrant_w_theta(number l)
// {
// 	number integ= ( gsl_sf_bessel_Jn(0,(l*theta_eval)) * l * Cl_GG_vec[0].value(l) ) / (2.*pi);

// 	return integ;
// }

// number psi_stats::integrant_gamma_T(number l)
// {

// 	number integ= ( gsl_sf_bessel_Jn(2,(l*theta_eval)) * l * Cl_GM_vec[0].value(l) ) / (2.*pi);

// 	return integ;
// }

// ///not done
// number psi_stats::integrant_cov(number l)
// {
// 	number integ;
// 	if(cov_case=="clustering")
// 	{
// 		number power1=Cl_GG_vec[rp1].value(l)+delta1noise;
// 		number power2=Cl_GG_vec[rp2].value(l)+delta2noise;
// 		number power3=Cl_GG_vec[rp3].value(l)+delta3noise;
// 		number power4=Cl_GG_vec[rp4].value(l)+delta4noise;
// 		number powers=power1*power2+power3*power4;
// 		integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
// 		return integ;
// 	}
// 	else if(cov_case=="ggl")
// 	{
// 		number power1=Cl_GG_vec[rp1].value(l)+delta1noise;
// 		number power2=Cl_MM_vec[rp2].value(l)+delta2noise;
// 		number power3=Cl_GM_vec[rp3].value(l);
// 		number power4=Cl_GM_vec[rp4].value(l);
// 		number powers=power1*power2+power3*power4;
// 		integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
// 		return integ;
// 	}
// 	else if(cov_case=="cross")
// 	{
// 		number power1=Cl_GG_vec[rp1].value(l)+delta1noise;
// 		number power2=Cl_GM_vec[rp2].value(l);
// 		number power3=Cl_GM_vec[rp3].value(l);
// 		number power4=Cl_GG_vec[rp4].value(l)+delta4noise;
// 		number powers=power1*power2+power3*power4;
// 		integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
// 		return integ;
// 	}
// 	else
// 	{
// 		clog<<"not a recognised covariance case, exiting now ..."<<endl;
// 		exit(1);
// 	}	
// }
// integrants end here


//find extrema for W_psi to integrate over
vector<number> psi_stats::determine_integration_limits(number last_x_value)
{
	//clog << "\n \n. LAST l value is .... -> l = " << last_x_value << "\n \n" << endl;
	//The value where the input function begins
	number first_x_value = lmin;

    //number of bins for the table of the oscillating function, 
    //you can play with this to see what effect it has on the function.
	const int Nbins = 1000000;

	// make table of the function values on a very fine grid
    vector<number> table_y(Nbins);
    vector<number> table_x(Nbins);

    // free old list, very important with vectors to do this. Easy to forget and mess up everything.
    vector<number> integ_limits_Cl;

	lthresh = 2.*pi/thetaMax/20.;
	number first_logx = log(lthresh);
	number last_logx = log(last_x_value);
	number dx = (last_logx - first_logx)/(Nbins-1.);
	//clog<<"first_logx="<<first_logx<<" last_logx="<<last_logx<<" dx="<<dx<<endl;
	//clog<<"exp(first_logx)="<<exp(first_logx)<<endl;

	//clog << "going through the table, type is:"<<Type<< endl;

	for(int i=0;i<Nbins;i++)
	{
		//clog<<"i="<<i<<endl;
        //here I'm making a log table, but you can do linear depending on what your function is.
        //table_x[i]=first_x_value+(last_x_value-first_x_value)/(Nbins-1.)*i;
        //log-space x-values
        table_x[i]=exp(first_logx + i*dx);
        //clog<<"table_x["<<i<<"]="<<table_x[i]<<endl;
        // table_y[i]=integrand_Q(mode,table_x[i]);
        table_y[i]=Wn_vec[nW].value(table_x[i]);
        
	}

	// go through list and pick minima/maxima (sort-of; does not need to be awfully exact)
	integ_limits_Cl.push_back(first_x_value);
	integ_limits_Cl.push_back(lthresh);


	for(int i=1;i<Nbins-1;i++)
	{
		if(  ((table_y[i-1]<table_y[i]) && (table_y[i+1]<table_y[i]))
		   || ((table_y[i-1]>table_y[i]) && (table_y[i+1]>table_y[i])))
			integ_limits_Cl.push_back(table_x[i]);
	}

	integ_limits_Cl.push_back(last_x_value);
	//clog<<"found the integ limits"<<endl;
    return integ_limits_Cl;
}


// number psi_stats::valueFunc_w_theta()
// {
// 	number result=0.;
// 	int steps=100;

// 	integration_type = "w_theta";

// 	vector<number> integ_limits_w;
// 	integ_limits_w.clear();
// 	integ_limits_w = determine_integration_limits(lmax);

// 	for(unsigned int i=0; (i+1)<integ_limits_w.size(); i++)
// 	{
// 		number res=gaussianIntegrate_gsl(*this,integ_limits_w[i],integ_limits_w[i+1],steps);
// 		result+=res;
// 	}

// 	return result;
// }

// number psi_stats::valueFunc_gam_T_theta()
// {
// 	number result=0.;
// 	int steps=100;

// 	integration_type = "gamT";

// 	vector<number> integ_limits_gamT;
// 	integ_limits_gamT.clear();
// 	integ_limits_gamT = determine_integration_limits(lmax);

// 	for(unsigned int i=0; (i+1)<integ_limits_gamT.size(); i++)
// 	{

// 		number res=gaussianIntegrate_gsl(*this,integ_limits_gamT[i],integ_limits_gamT[i+1],steps);

// 		result+=res;
// 	}

// 	return result;
// }


number psi_stats::value_psi(int n, string integration_type1)
{

	nW=n-nMinimum;

	number result=0.;
	int steps;

	integration_type = integration_type1;

	for(unsigned int i=0; (i+1)<integ_limits_vec_vec[nW].size(); i++)
	{
		if (i==0)
			steps=100;
		else
			steps=10;

		number res=gaussianIntegrate_gsl(*this,integ_limits_vec_vec[nW][i],integ_limits_vec_vec[nW][i+1],steps);
		result+=res;
	}

	return result / (2.*pi);
}


vector<number> psi_stats::openFile(string FilePath,string FileNames)
{
	string FileName=string(FilePath)+string(FileNames);
	string str;
	number temp;
	vector<number>  Ksi_vecvec;
	Ksi_vecvec.clear();
	ifstream KsiFile((FileName).c_str());
 	clog<<"reading, "<<FileName<<endl;
	if(KsiFile.fail())
	{
		clog<<"error occured during opening: "<<FileName<<endl;
		exit(1);
	}
	///lets change this
	string line;
	while(getline(KsiFile,line))
	{
		stringstream stream(line);
		str=stream.peek();
		if (str=="#")
		{

		}
		else
		{

			stream>>temp;
			// clog<<temp<<endl;
			Ksi_vecvec.push_back(temp);

		}
	}
	KsiFile.close();
	return Ksi_vecvec;
}



void psi_stats::input_for_setCl(vector<number> ell, vector<vector<number> > cl_gg, vector<vector<number> > cl_gm,vector<vector<number> > cl_mm)
{

	//clog << "In input for Cl" << endl;
	setPower(ell,cl_gg,cl_gm,cl_mm);
	//Set a global boolean parameter which tells you wether you have computed the min/max already at this point
	//If find_min_max=true then you must locate the min/max for each mode
	// if (find_min_max)
	// {
	// 	vector<number> integ_lim;
	// 	//Ok here now compute the min-max for all modes! -> nmax has already been set
	// 	//clog << "\n \n \n &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n \n \n " << endl;
	// 	//clog << "I am now computing the integration limits" << endl;
	// 	//clog << "nMaximum = " << nMaximum << endl;

	// 	//Loop of the modes
	// 	for (int modeCount=0; modeCount<=nMaximum-nMinimum; modeCount++)
	// 	{

	// 		//This pointer points to the object for which the member function is called. Static member functions do not have a this pointer.
	// 		nW = modeCount;
	// 		rp = 0;
	// 		integ_lim.clear();

	// 		//auto start = chrono::high_resolution_clock::now();

	// 		integ_lim = determine_integration_limits(lmax,"FAB_gg");
	// 		// cout.precision(10);
	// 		// if (modeCount==0){for (int j=0;j<integ_limits_Cl.size();j++){cout << integ_limits_Cl[j] << endl;}}

	// 		//clog << "integ_limits_Cl.size() = " << integ_lim.size() << endl;
	// 		//auto finish = chrono::high_resolution_clock::now();
	// 		//chrono::duration<number> elapsed = finish - start;
	// 		//clog << "Elapsed time to compute integ limits for mode " << nW+1 << " : " << elapsed.count() << " s\n";
	// 		integ_limits_vec_vec.push_back(integ_lim);

	// 	}

	// 	find_min_max = false;
	// }
	//clog<<"input_for_setCl finished"<<endl;
}


// void psi_stats::input_for_setRealSpace(vector<number> theta, vector<vector<number> > w_CosmoSIS, 
// 	vector<vector<number> > ggl_CosmoSIS)
// {
// 	clog << "In input for real space" << endl;
// 	setPower_realspace(theta,w_CosmoSIS,ggl_CosmoSIS);
// }

// vector<number> psi_stats::InterpolateRealSpace(vector<number> log_theta, vector<vector<number> > real_space_corr, 
// 	vector<number> interp_theta)
// {

// 	//Clear the function COSEBIs object
// 	real_space_vec.clear();


// 	for(int r=0; r<real_space_corr.size(); r++)
// 		real_space_vec.push_back(function_cosebis());

// 	for(int r=0; r<real_space_corr.size(); r++)
// 	{
// 		//Load values of the real space correlation function into a function COSEBIs table in ARCMINUTES
// 		clog << "Size of Cosmosis_theta_input = " << log_theta.size() << "\n" << endl;
// 		clog << "Size of real_space_corr[r] = " << real_space_corr[r].size() << "\n" << endl;
// 		real_space_vec[r].loadWithValues(log_theta,real_space_corr[r],false);
// 		real_space_vec[r].extrapolationOn();
// 	}

// 	//Define the interpolated vector of values of the real space correlation vector
// 	vector<number> real_space_vec_interp;
// 	real_space_vec_interp.clear();

// 	clog << "Cosmosis_theta_input[0] = " << log_theta[0] << "\n" << "real_space_corr[0] =  " << real_space_corr[0][0] << "\n" << "real_space_vec[0].value(interp_theta[0]) = " << real_space_vec[0].value(interp_theta[0])  << "\n" << "interp_theta[0]" << interp_theta[0] << endl;

// 	//Loop over the interpolated theta values (Note these interp_theta values are ALSO arcminute!!)
// 	for (int mbin=0; mbin<interp_theta.size(); mbin++)
// 	{
// 		real_space_vec_interp.push_back(real_space_vec[0].value(interp_theta[mbin]));
// 	}	

// 	return real_space_vec_interp;

// }

// void psi_stats::setPower_realspace(vector<number> log_theta, vector<vector<number> > w_CosmoSIS, 
// 	vector<vector<number> > ggl_CosmoSIS)
// {

// 	ggl_CosmoSIS_vec.clear();
// 	w_CosmoSIS_vec.clear();

// 	for(int r=0; r<ggl_CosmoSIS.size(); r++)
// 	{
// 		ggl_CosmoSIS_vec.push_back(function_cosebis());
// 		w_CosmoSIS_vec.push_back(function_cosebis());
// 	}

// 	for(int r=0; r<ggl_CosmoSIS.size(); r++)
// 	{
// 		ggl_CosmoSIS_vec[r].loadWithValues(log_theta,ggl_CosmoSIS[r],false);
// 		ggl_CosmoSIS_vec[r].extrapolationOn();

// 		w_CosmoSIS_vec[r].loadWithValues(log_theta,w_CosmoSIS[r],false);
// 		w_CosmoSIS_vec[r].extrapolationOn();
// 	}
// }


void psi_stats::setPower(vector<number> log_ell,vector<vector<number> > InputPower_clustering,
 vector<vector<number> > InputPower_ggl,vector<vector<number> > InputPower_matter)
{
	
	Cl_GG_vec.clear();
	Cl_GM_vec.clear();
	Cl_MM_vec.clear();
	nPair_GG=InputPower_clustering.size();
	nPair_GM=InputPower_ggl.size();
	nPair_MM=InputPower_matter.size();
	nBins_GG=(sqrt(8*nPair_GG+1)-1)/2;
	nBins_MM=(sqrt(8*nPair_MM+1)-1)/2;

	clog<<"nBins_GG="<<nBins_GG<<" nBins_MM="<<nBins_MM<<endl;
	clog<<"nPair_GM="<<nPair_GM<<endl;
	//these have to be separate otherwise does not work
	for(int r=0; r<nPair_GM; r++)
		Cl_GM_vec.push_back(function_cosebis());
	for(int r=0; r<nPair_GM; r++)
	{
		Cl_GM_vec[r].loadWithValues(log_ell,InputPower_ggl[r],true);
		Cl_GM_vec[r].extrapolationOn();
	}

	for(int r=0; r<nPair_GG; r++)
		Cl_GG_vec.push_back(function_cosebis());
	for(int r=0; r<nPair_GG; r++)
	{
		Cl_GG_vec[r].loadWithValues(log_ell,InputPower_clustering[r],true);
		Cl_GG_vec[r].extrapolationOn();
	}

	for(int r=0; r<nPair_MM; r++)
		Cl_MM_vec.push_back(function_cosebis());
	for(int r=0; r<nPair_MM; r++)
	{
		Cl_MM_vec[r].loadWithValues(log_ell,InputPower_matter[r],true);
		Cl_MM_vec[r].extrapolationOn();
	}
}


void psi_stats::setCl(vector<number> log_ell,vector<number> Cl)
{
	Cl_func.loadWithValues(log_ell,Cl,true);
	Cl_func.extrapolationOn();

	integration_type = "psi_gg";
	if (find_min_max)
	{
		vector<number> integ_lim;
		//Ok here now compute the min-max for all modes! -> nmax has already been set
		//clog << "\n \n \n &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& \n \n \n " << endl;
		//clog << "I am now computing the integration limits" << endl;
		//clog << "nMaximum = " << nMaximum << endl;

		//Loop of the modes
		for (int modeCount=0; modeCount<=nMaximum-nMinimum; modeCount++)
		{

			//This pointer points to the object for which the member function is called. Static member functions do not have a this pointer.
			nW = modeCount;
			integ_lim.clear();

			//auto start = chrono::high_resolution_clock::now();

			integ_lim = determine_integration_limits(lmax);
			// cout.precision(10);
			// if (modeCount==0){for (int j=0;j<integ_limits_Cl.size();j++){cout << integ_limits_Cl[j] << endl;}}

			//clog << "integ_limits_Cl.size() = " << integ_lim.size() << endl;
			//auto finish = chrono::high_resolution_clock::now();
			//chrono::duration<number> elapsed = finish - start;
			//clog << "Elapsed time to compute integ limits for mode " << nW+1 << " : " << elapsed.count() << " s\n";
			integ_limits_vec_vec.push_back(integ_lim);

		}
		find_min_max = false;
	}
}


// number psi_stats::Trap(number a, number b, int steps)
// {
//   number s = 0;
//   number h = (b-a)/steps;

//   for (int i = 0; i < steps; ++i)
//     s += funcTrap(a + h*i, h);

//   return h*s;

// }

// number psi_stats::funcTrap(number x, number h)
// {
//     return (integrant_w_theta(x) + integrant_w_theta(x+h))/2;
// }


// void psi_stats::input_for_setMocks(int thetamin_arc, int thetamax_arc, int NBINS,int LOSmin,int LOSmax, 
// 	string mock_w_path, string mock_ggl_path)
// {

// 	char U_gg_filename[256];
// 	char U_gm_filename[256];

// 	number thetaMinArc = radians_to_arcmin(thetaMin);
// 	number thetaMaxArc = radians_to_arcmin(thetaMax);


//     snprintf(U_gg_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GG_min_%.1f_max_%.1f_LOS_%d.txt", thetaMinArc, thetaMaxArc ,LOSmin);
//     printf(U_gg_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GG_min_%.1f_max_%.1f_LOS_%d.txt", thetaMinArc, thetaMaxArc,LOSmin);
//     clog << " " << endl;
//     snprintf(U_gm_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GM_min_%.1f_max_%.1f_LOS_%d.txt", thetaMinArc, thetaMaxArc,LOSmin);
//     printf(U_gm_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GM_min_%.1f_max_%.1f_LOS_%d.txt", thetaMinArc, thetaMaxArc,LOSmin);

//  	//Initialize all variables
// 	string fileName_ggl,fileName_w,fileName_theta,U_filter_name,Q_filter_name;
// 	vector<number> ggl_mock,w_mock,psi_GM_vec,psi_GG_vec,U_mock_vec,Q_mock_vec;
// 	vector<vector<number> > psi_GG_mock,psi_GM_mock,U_mock,Q_mock;
// 	number psi_GG,psi_GM;

// 	int LOS_counter=0; //The number of exsisting LOS
// 	int LOS_success = 682; //402


// 	if (LOSmax-LOSmin==0)
// 	{
// 		LOS_success = 1;
// 	}


// 	//(Column, row)
// 	matrix tableUGMmock(LOS_success,nMaximum);
// 	matrix tableUGGmock(LOS_success,nMaximum);



// 	//Open theta file
// 	fileName_theta = string("theta_")+toString(thetamin_arc)+string(".0_")+toString(thetamax_arc)+string(".0_")+toString(NBINS)+string("_bins_LOS_74.txt");
// 	vector<number> theta_mock = arcmin_to_radians_vec(openFile(mock_ggl_path,fileName_theta));

// 	//Initialize the class psi_filters
// 	psi_filters WnL_mock(thetaMin, thetaMax, NBINS, nMaximum, lmin, lmax, lbins);
// 	clog<<"?????????????????????????????????????????????"<<endl;
// 	clog << "theta_mock.size() " << theta_mock.size() << endl;
// 	clog<<"?????????????????????????????????????????????"<<endl;

// 	clog<<"?????????????????????????????????????????????"<<endl;
// 	clog << "theta_mock[0] " << theta_mock[0] << "minTH = " << thetaMin << endl;
// 	clog<<"?????????????????????????????????????????????"<<endl;
// 	clog << "theta_mock[NBINS-1] " << theta_mock[NBINS-1] << endl;
// 	clog<<"?????????????????????????????????????????????"<<endl;

// 	//Set the maximum number of modes to be analysed
// 	WnL_mock.initialize_mode(nMaximum);

// 	//Calculate the step-size of integration in theta-space
// 	//CHECK that this is true....is it log-space or is it linear space.....check all
// 	number dt = theta_mock[2] - theta_mock[1];
// 	clog << "dt = " << dt << endl;


// 	matrix U_filter_mat(nMaximum+1,theta_mock.size());
// 	matrix Q_filter_mat(nMaximum+1,theta_mock.size());

// 	for (int i=0; i<theta_mock.size(); i++)
// 	{
// 		U_filter_mat.load(0,i,theta_mock[i]);
// 		Q_filter_mat.load(0,i,theta_mock[i]);
// 	}


// 	//Loop over each individual LOS from the mocks
// 	for (int LOS=LOSmin; LOS<=LOSmax; LOS++)
// 	{


// 		clog << "Doing analysis of LOS " << LOS << endl;

// 		//Declare the ggl and clustering file names for the mocks
// 		// char U_gg_filename[256];
// 		// char U_gm_filename[256];
// 	 //    snprintf(fileName_ggl, 256, "/disk09/ifrisw/psi_Stats/mock_outputs/ggl_LOS_%d_bins_%d_min_%.1f_max_%.1f",LOS, NBINS, thetamin_arc, thetamax_arc);
// 	 //    printf(fileName_ggl, 256, "/disk09/ifrisw/psi_Stats/mock_outputs/ggl_LOS_%d_bins_%d_min_%.1f_max_%.1f",LOS, NBINS, thetamin_arc, thetamax_arc);
// 	 //    snprintf(fileName_w, 256, "/disk09/ifrisw/psi_Stats/mock_outputs/w_theta_%.1f_%.1f_%d_bins_LOS_%d_.txt",thetamin_arc, thetamax_arc, LOS, NBINS);
// 	 //    printf(fileName_w, 256, "/disk09/ifrisw/psi_Stats/mock_outputs/w_theta_%.1f_%.1f_%d_bins_LOS_%d_.txt",thetamin_arc, thetamax_arc, LOS, NBINS);

// 		fileName_ggl = string("ggl_LOS_")+toString(LOS)+string("_bins_")+toString(NBINS)+
// 						string("_min_")+toString(thetamin_arc)+string(".0_max_")+toString(thetamax_arc)+string(".0");
// 		fileName_w = string("w_theta_") + toString(thetamin_arc) + string(".0_")+toString(thetamax_arc)
// 						+string(".0_")+toString(NBINS)+string("_bins_LOS_")+toString(LOS)+string(".txt");
	
// 		ifstream testLOSggl((string(mock_ggl_path)+string(fileName_ggl)).c_str());
// 		ifstream testLOSw((string(mock_w_path)+string(fileName_w)).c_str());

// 		//Check if this LOS exsists
// 		if(testLOSggl.fail())
// 		{
// 			//If it doesn't exsist log the LOS ID
// 			clog << "LOS " << LOS << " does not exist." << endl;
// 			clog << string(mock_ggl_path) << string(fileName_ggl) << endl;
// 		}

// 		else if(testLOSw.fail())
// 		{
// 			//If it doesn't exsist log the LOS ID
// 			clog << "LOS " << LOS << " does not exist." << endl;
// 			clog << string(mock_w_path) << string(fileName_w) << endl;
// 		}

// 		else
// 		{

// 			//If it does exsist open the file and save into a vector
// 			clog << "Doing LOS " << LOS << endl;
// 			w_mock = openFile(mock_w_path,fileName_w);		
// 			ggl_mock = openFile(mock_ggl_path,fileName_ggl);

// 			//Loop through modes.
// 			for (int modeIDX=1; modeIDX<=nMaximum; modeIDX++)
// 			{
// 				psi_GG=0.;
// 				psi_GM=0.;
// 				//Loop through each of the Theta Bins.
// 				for (int thetaIDX=0; thetaIDX<NBINS; thetaIDX++)
// 				{

// 					U_filter_mat.load(modeIDX,thetaIDX,WnL_mock.U(modeIDX,theta_mock[thetaIDX]));
// 					Q_filter_mat.load(modeIDX,thetaIDX,WnL_mock.analyticalQ(theta_mock[thetaIDX],modeIDX));

// 					//If theta bin is first OR last theta bin then down-weight
// 					if ((thetaIDX==0) || (thetaIDX==NBINS-1))
// 					{

// 						psi_GG  += (dt/2.)*theta_mock[thetaIDX]*WnL_mock.U(modeIDX,theta_mock[thetaIDX])*w_mock[thetaIDX];
// 						psi_GM  += (dt/2.)*theta_mock[thetaIDX]*WnL_mock.analyticalQ(theta_mock[thetaIDX],modeIDX)*ggl_mock[thetaIDX];

// 					}
// 					else
// 					{
// 						psi_GG  += dt*theta_mock[thetaIDX]*WnL_mock.U(modeIDX,theta_mock[thetaIDX])*w_mock[thetaIDX];
// 						psi_GM  += dt*theta_mock[thetaIDX]*WnL_mock.analyticalQ(theta_mock[thetaIDX],modeIDX)*ggl_mock[thetaIDX];
// 					}
// 				}
// 				tableUGGmock.load(LOS_counter,modeIDX-1,psi_GG);
// 				tableUGMmock.load(LOS_counter,modeIDX-1,psi_GM);
// 			}
// 			LOS_counter += 1;
// 		}
// 		//After each LOS important to clear the vectors
// 		w_mock.clear();
// 		ggl_mock.clear();
// 	}

// 	clog << LOS_counter <<endl;

// 	// char U_gg_filename[256];
// 	// char U_gm_filename[256];
//  //    snprintf(U_gg_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GG_min_%.1f_max_%.1f_LOS_%d.txt", thetamin_arc, thetamax_arc ,LOSmin);
//  //    printf(U_gg_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GG_min_%.1f_max_%.1f_LOS_%d.txt", thetamin_arc, thetamax_arc ,LOSmin);
//  //    snprintf(U_gm_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GM_min_%.1f_max_%.1f_LOS_%d.txt", thetamin_arc, thetamax_arc ,LOSmin);
//  //    printf(U_gm_filename, 256, "/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GM_min_%.1f_max_%.1f_LOS_%d.txt", thetamin_arc, thetamax_arc ,LOSmin);

// 	// U_gg_filename = string("/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GG_min_")+toString(thetamin_arc)+string("_max_")+toString(thetamax_arc)+string("_LOS_")+toString(LOS)+string(".txt");
// 	// U_gm_filename = string("/disk09/ifrisw/psi_Stats/mock_psi_outputs/psi_GM_min_")+toString(thetamin_arc)+string("_max_")+toString(thetamax_arc)+string("_LOS_")+toString(LOS)+string(".txt");
// 	// U_filter_name = string("/disk09/ifrisw/psi_Stats/mock_psi_outputs/U_filter_mock_")+toString(thetamin_arc)+string("_max_")+toString(thetamax_arc)+string(".txt");
// 	// Q_filter_name = string("/disk09/ifrisw/psi_Stats/mock_psi_outputs/Q_filter_mock_")+toString(thetamin_arc)+string("_max_")+toString(thetamax_arc)+string(".txt");

// 	tableUGGmock.printOut(U_gg_filename);
// 	tableUGMmock.printOut(U_gm_filename);
// 	// U_filter_mat.printOut(U_filter_name.c_str());
// 	// Q_filter_mat.printOut(Q_filter_name.c_str());

// }


// vector<number> psi_stats::psi_real_space_calculator(vector<number> w_in,vector<number> ggl_in,
// 			 vector<number>theta_in)
// {

// 	//Initialize the class psi_filters
// 	psi_filters WnL_in(thetaMin, thetaMax, theta_in.size(), nMaximum, lmin, lmax, lbins);
// 	//Set the maximum number of modes to be analysed
// 	WnL_in.initialize_mode(nMaximum);

// 	clog << "\n" << "\n" << "\n" << "\n psi_real_space_calculator" << "\n" << "\n" << "\n" << endl;

// 	number dt = theta_in[3]-theta_in[2];
// 	number psi_GG,psi_GM,thetaEvaluated;

// 	vector<number> psi_GM_vec,psi_GG_vec;
// 	psi_GM_vec.clear();psi_GG_vec.clear();

// 	//Loop through modes.
// 	for (int modeIDX=1; modeIDX<=nMaximum; modeIDX++)
// 	{
// 		psi_GG=0.;psi_GM=0.;
		
// 		//Loop through each of the Theta Bins.
// 		for (int thetaIDX=0; thetaIDX<theta_in.size(); thetaIDX++)
// 		{
// 			thetaEvaluated = theta_in[thetaIDX];
// 			//If theta bin is first OR last theta bin then down-weight
// 			if ((thetaIDX==0) || (thetaIDX==theta_in.size()-1))
// 			{
// 				// psi_GG  += (dt/2.)*thetaEvaluated*WnL_mock.U(modeIDX,thetaEvaluated)*w_CosmoSIS_vec[0].value(thetaEvaluated);
// 				// psi_GM  += (dt/2.)*thetaEvaluated*WnL_mock.analyticalQ(thetaEvaluated,modeIDX)*ggl_CosmoSIS_vec[0].value(thetaEvaluated);
// 				psi_GG  += (dt/2.)*thetaEvaluated*WnL_in.U(modeIDX,thetaEvaluated)*w_in[thetaIDX];
// 				psi_GM  += (dt/2.)*thetaEvaluated*WnL_in.analyticalQ(thetaEvaluated,modeIDX)*ggl_in[thetaIDX];
// 				// clog << thetaEvaluated << " " << WnL_mock.U(modeIDX,thetaEvaluated) << endl;
// 			}
// 			else
// 			{

// 				// psi_GG  += dt*thetaEvaluated*WnL_mock.U(modeIDX,thetaEvaluated)*w_CosmoSIS_vec[0].value(thetaEvaluated);
// 				// psi_GM  += dt*thetaEvaluated*WnL_mock.analyticalQ(thetaEvaluated,modeIDX)*ggl_CosmoSIS_vec[0].value(thetaEvaluated);
// 				psi_GG  += dt*thetaEvaluated*WnL_in.U(modeIDX,thetaEvaluated)*w_in[thetaIDX];
// 				psi_GM  += dt*thetaEvaluated*WnL_in.analyticalQ(thetaEvaluated,modeIDX)*ggl_in[thetaIDX];
// 			}
// 		}

// 		psi_GM_vec.push_back(psi_GM);
// 		psi_GG_vec.push_back(psi_GG);
// 		clog << psi_GG << " " << psi_GM << endl;
// 	}


// 	vector<number> psi;
// 	psi.clear();

// 	for (int ct=0; ct<psi_GG_vec.size();ct++)
// 		psi.push_back(psi_GG_vec[ct]);
// 	for (int ct=0; ct<psi_GM_vec.size();ct++)
// 		psi.push_back(psi_GM_vec[ct]);

// 	return(psi);
// }



// vector<number> psi_stats::mean_mocks(vector<vector<number> > input_vec_vec, int Ncols, int Nrows)
// {
// 	number element; vector<number> mean_vec;
// 	for(int rowIDX=0; rowIDX<Nrows; rowIDX++)
// 	{
// 		element=0.;
// 		for(int colIDX=0; colIDX<Ncols; colIDX++)
// 		{
// 			element += input_vec_vec[colIDX][rowIDX];
// 			clog << rowIDX << " " << colIDX << endl;
// 		}
// 		mean_vec.push_back(element/Ncols);
// 	}
// 	return mean_vec;
// }


//vec[a][b] --> b=mode:1,2,3, a=LOS increasing


// vector<matrix> psi_stats::returnApertureStatisticsIntegrant(vector <number> ell_vec, number theta, 
// 	string ap_type, int multi, int redshift)
// {

// 	//Initialise vectors
// 	vector<matrix> IntegAndU_mat; IntegAndU_mat.clear();
// 	matrix integ_mat((multi+1), ell_vec.size());
// 	matrix U_mat(2, ell_vec.size());


// 	//Loop over the ell modes - This U-filter function is the I[l0] in the literature and is the same for all aperture stats
// 	for(unsigned int i=0; i<ell_vec.size();i++)
// 	{
// 		//theta_l = theta/2*l
// 		number thetaell=theta/2.*ell_vec[i];
// 		//U_t = 24*J4(theta*l)/(theta_l)^2
// 		number U_t=24.*gsl_sf_bessel_Jnu(4,thetaell)/(thetaell*thetaell);
// 		//Load values into table ---|ell|U_t|---
// 		U_mat.load(0,i,ell_vec[i]);
// 		U_mat.load(1,i,U_t);
// 	}
	

// 	int redshiftPair=redshift;

// 	for(unsigned int i=0; i<ell_vec.size();i++)
// 	{
// 		integ_mat.load(0,i,ell_vec[i]);

// 		//Do I want to calculate both NAP and MNAP or just one ;)
// 		if (multi==2)
// 		{

// 			number integ_NAP=U_mat.get(1,i)*U_mat.get(1,i)*ell_vec[i]*Cl_GG_vec[redshiftPair].value(ell_vec[i]);
// 			number integ_MNAP=U_mat.get(1,i)*U_mat.get(1,i)*ell_vec[i]*Cl_GM_vec[redshiftPair].value(ell_vec[i]);

// 			integ_mat.load(1,i,integ_NAP);
// 			integ_mat.load(2,i,integ_MNAP);
// 			clog << ell_vec[i] << " " << integ_NAP << " " << integ_MNAP << " " << endl;

// 		}

// 		else if (multi==1)
// 		{

// 			if (ap_type=="NAP")
// 			{
// 				number integ=U_mat.get(1,i)*U_mat.get(1,i)*ell_vec[i]*Cl_GG_vec[redshiftPair].value(ell_vec[i]);
// 				integ_mat.load(1,i,integ);
// 			}
// 			else if (ap_type=="MAPNAP")
// 			{
// 				number integ=U_mat.get(1,i)*U_mat.get(1,i)*ell_vec[i]*Cl_GM_vec[redshiftPair].value(ell_vec[i]);
// 				integ_mat.load(1,i,integ);
// 			} 

// 		}

// 	}

// 	IntegAndU_mat.push_back(integ_mat);
// 	// IntegAndU_mat.push_back(U_mat);
// 	return IntegAndU_mat;
// }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//covariance functions:



void psi_stats::setNoise(vector<number> sigma_e_vec,vector<number> nBar_shear_vec,vector<number> nBar_gals_vec)
{
	clog<<"setting noise in psi_stats"<<endl;
	noise_gg_vec.clear();
	noise_mm_vec.clear();

	for(int bin=0; bin<nBar_gals_vec.size(); bin++)
	{
		noise_gg_vec.push_back(1./nBar_gals_vec[bin]);
		clog<<"noise_gg_vec["<<bin<<"]="<<noise_gg_vec[bin]<<endl;
	}

	for(int bin=0; bin<sigma_e_vec.size(); bin++)
	{
		noise_mm_vec.push_back(sigma_e_vec[bin]*sigma_e_vec[bin]/(2.*nBar_shear_vec[bin]));
		clog<<"noise_mm_vec["<<bin<<"]="<<noise_mm_vec[bin]<<endl;
	}
}

//NOTE: the value of LHIGH is important for the off-diagonals. 
// For better precision use a bigger high
// void psi_stats::determine_integration_limitsCov()
// {
// /* Idea: find a possibly complete list of consecutive local 
// minima/maxima of oscillating integrant and integrate between them
// */
// 	const int Nbins = 1000000;
// 	number LHIGH=high*20./thetaMax;
// 	//clog<<"LHIGH="<<LHIGH<<endl;
// 	number LLOW=lmin;
// 	// free old list
// 	integ_limits.clear();
// 	// make table of integrant values (Wn's only) on a very fine grid
// 	matrix table(2,Nbins);
// 	number lthresh=pi/thetaMax/2.;
// 	for(int i=0;i<Nbins;i++)
// 	{
// 		table.load(0,i,exp(log(lthresh)+log(LHIGH/lthresh)/(Nbins-1.)*i));
// 		table.load(1,i,Wn_vec[nW].value(table.get(0,i))*Wn_vec[mW].value(table.get(0,i)));
// 	}
// // go through list and pick minima/maxima (sort-of; does not need to be awfully exact)
// 	integ_limits.push_back(LLOW);
// 	integ_limits.push_back(lthresh);
// 	for(int i=1;i<Nbins-1;i++)//2->1
// 		if ((table.get(1,i-1)<table.get(1,i)&& table.get(1,i+1)<table.get(1,i))
// 		|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
// 		      integ_limits.push_back(table.get(0,i));
// 	integ_limits.push_back(LHIGH);
// }

// number psi_stats::valueCov(int n1,int m1)
// {
// 	nW=n1-1;
// 	mW=m1-1;
// 	determine_integration_limitsCov();
// 	number result= 0.;	
// 	for(unsigned int i=0;(i+1)<integ_limits.size();i++)
// 	{
// 		number res=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],40,"covariance");
// 		result+=res;
// 	}
// 	return result/2./pi;
// }


vector<vector<matrix* > > psi_stats::determine_integration_limitsCov()
{
	// Idea: find a possibly complete list of consecutive local 
	//minima/maxima of oscillating integrant and integrate between them

	int Nbins = 1000000;
	number LHIGH=high*20./thetaMax;
	//clog<<"LHIGH="<<LHIGH<<endl;
	number LLOW=lmin;
	number lthresh=pi/thetaMax/2.;
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(lthresh)+log(LHIGH/lthresh)/(Nbins-1.)*i));
	}

	//integ_limits_cov.clear();
	vector<vector<matrix* > > integ_limits_cov;
	vector <matrix*> integ_mat_vec;
	for(int n=0; n<(nMaximum-nMinimum+1); n++)
	{
		integ_mat_vec.clear();
		for(int m=0; m<(nMaximum-nMinimum+1); m++)
		{
			vector <number> integ_lim_vec;
			for(int i=0;i<Nbins;i++)
			{
				table.load(1,i,Wn_vec[n].value(table.get(0,i))*Wn_vec[m].value(table.get(0,i)));
			}
		// go through list and pick minima/maxima (sort-of; does not need to be awfully exact)
			integ_lim_vec.push_back(LLOW);
			integ_lim_vec.push_back(lthresh);
			for(int i=1;i<Nbins-1;i++)//2->1
				if ((table.get(1,i-1)<table.get(1,i)&& table.get(1,i+1)<table.get(1,i))
				|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
				      integ_lim_vec.push_back(table.get(0,i));
			integ_lim_vec.push_back(LHIGH);

			int size=integ_lim_vec.size();
			clog<<"n="<<n<<" m="<<m<<" integ_lims="<<size<<endl;
			integ_mat_vec.push_back(new matrix(size));
			for(int i=0; i<integ_lim_vec.size(); i++)
			{
				integ_mat_vec[m]->load(i,integ_lim_vec[i]);
			}
		}
		integ_limits_cov.push_back(integ_mat_vec);
	}
	return integ_limits_cov;
}





number psi_stats::valueCov(int n1,int m1, vector<vector<matrix* > > integ_limits_cov)
{
	nW=n1-1;
	mW=m1-1;
	number result= 0.;	
	integration_type = "cov";
	for(unsigned int i=0;(i+1)<integ_limits_cov[nW][mW]->size();i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits_cov[nW][mW]->get(i),integ_limits_cov[nW][mW]->get(i+1),10);
		result+=res;
	}
	return result/2./pi;
}

// number psi_stats::valueCov(int n1,int m1)
// {
// 	return 0.;
// }



matrix psi_stats::calCov(number Area_GG,number Area_GM,number Area_cross,bool Cross)
{
	clog<<"calculating the covariance in psi_stats"<<endl;
	matrix CMT(nMaximum*(nPair_GG+nPair_GM),nMaximum*(nPair_GG+nPair_GM));
	setWFilters(nMaximum, nMinimum);
	//for testing set it all to zero 
	//CMT.zero();

	clog<<"power spectra are set, nBins_GG="<<nBins_GG<<endl;
	cov_case="clustering";
	vector<vector<matrix* > >integ_limits_cov=determine_integration_limitsCov();
	clog<<"cov_case="<<cov_case<<endl;
	for(int bin1=0; bin1<nBins_GG; bin1++)
	{
		for(int bin2=bin1; bin2<nBins_GG; bin2++)
		{
			for(int bin3=bin1; bin3<nBins_GG; bin3++)
			{
				for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins_GG;bin4++)
				{
					clog<<"GG: bin1="<<bin1<<" bin2="<<bin2<<" bin3="<<bin3<<" bin4="<<bin4<<endl;
					//clog<<"noise_vec="<<noise_vec[bin1]<<"  "<<noise_vec[bin2]<<endl;
					int n,m=0;
					rp1=calP(nBins_GG,bin1,bin3);
					rp2=calP(nBins_GG,bin2,bin4);
					rp3=calP(nBins_GG,bin1,bin4);
					rp4=calP(nBins_GG,bin2,bin3);
					//this should depend on the type of correlation
					delta1noise=delta(bin1,bin3)*noise_gg_vec[bin1];
					delta2noise=delta(bin2,bin4)*noise_gg_vec[bin2];
					delta3noise=delta(bin1,bin4)*noise_gg_vec[bin1];
					delta4noise=delta(bin2,bin3)*noise_gg_vec[bin2];
					//the pair considered for the Eparam_vec
					int p1=calP(nBins_GG,bin1,bin2);
					int p2=calP(nBins_GG,bin3,bin4);
					for(int i=(nMaximum-nMinimum+1)*p1,n=nMinimum;i<nMaximum*(p1+1);i++,n++)
					{
						for(int j=(nMaximum-nMinimum+1)*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
						{
							//clog<<"n="<<n<<" m="<<m<<endl;
							CMT.load(i,j,valueCov(n,m,integ_limits_cov)/Area_GG);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
					CMT.printOut((string("C_GG.ascii")).c_str(),10);
				}
			}
		}
	}
	//CMT.printOut((string("C_GG.ascii")).c_str(),4);

	//CMT.zero();
	cov_case="ggl";
	clog<<"cov_case="<<cov_case<<endl;
	for(int bin1=0; bin1<nBins_GG; bin1++)
	{
		for(int bin2=0; bin2<nBins_MM; bin2++)
		{
			for(int bin3=0; bin3<nBins_GG; bin3++)
			{
				for(int bin4=0; bin4<nBins_MM;bin4++)
				{
					clog<<"GM: bin1="<<bin1<<" bin2="<<bin2<<" bin3="<<bin3<<" bin4="<<bin4<<endl;
					int n,m=0;
					rp1=calP(nBins_GG,bin1,bin3);
					rp2=calP(nBins_MM,bin2,bin4);
					rp3=calP_cross(nBins_GG,nBins_MM,bin4,bin1);
					rp4=calP_cross(nBins_GG,nBins_MM,bin2,bin3);
					//clog<<"rp1="<<rp1<<" rp2="<<rp2<<" rp3="<<rp3<<" rp4="<<rp4<<endl;
					//this should depend on the type of correlation
					delta1noise=delta(bin1,bin3)*noise_gg_vec[bin1];
					delta2noise=delta(bin2,bin4)*noise_mm_vec[bin2];
					delta3noise=0.;
					delta4noise=0.;
					//the pair considered for the Eparam_vec
					int p1=calP_cross(nBins_GG,nBins_MM,bin2,bin1);
					int p2=calP_cross(nBins_GG,nBins_MM,bin4,bin3);
					for(int i=(nMaximum-nMinimum+1)*p1,n=nMinimum;i<nMaximum*(p1+1);i++,n++)
					{
						for(int j=(nMaximum-nMinimum+1)*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
						{
							//clog<<"n="<<n<<" m="<<m<<endl;
							int i_ggl=i+nMaximum*nPair_GG;
							int j_ggl=j+nMaximum*nPair_GG;
							CMT.load(i_ggl,j_ggl,valueCov(n,m,integ_limits_cov)/Area_GM);
							CMT.load(j_ggl,i_ggl,CMT.get(i_ggl,j_ggl));
							// CMT.load(i_ggl+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i_ggl,j_ggl));
							// CMT.load(j_ggl+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i_ggl,j_ggl));
						}
					}
				}
			}
		}
	}
	//CMT.printOut((string("C_GM.ascii")).c_str(),4);

	//CMT.zero();
	if(Cross)
	{
		cov_case="cross";
		clog<<"cov_case="<<cov_case<<endl;
		//clog<<"nPair_GG="<<nPair_GG<<endl;
		for(int bin1=0; bin1<nBins_GG; bin1++)
		{
			for(int bin2=bin1; bin2<nBins_GG; bin2++)
			{
				for(int bin3=0; bin3<nBins_GG; bin3++)
				{
					for(int bin4=0; bin4<nBins_MM;bin4++)
					{
						clog<<"Cross: bin1="<<bin1<<" bin2="<<bin2<<" bin3="<<bin3<<" bin4="<<bin4<<endl;
						int n,m=0;
						rp1=calP(nBins_GG,bin1,bin3);
						rp2=calP_cross(nBins_GG,nBins_MM,bin4,bin2);
						rp3=calP_cross(nBins_GG,nBins_MM,bin4,bin1);
						rp4=calP(nBins_GG,bin2,bin3);
						//this should depend on the type of correlation
						delta1noise=delta(bin1,bin3)*noise_gg_vec[bin1];
						delta2noise=delta(bin2,bin4);
						delta3noise=0.;
						delta4noise=delta(bin2,bin3)*noise_gg_vec[bin2];
						//the pair considered for the Eparam_vec
						int p1=calP(nBins_GG,bin1,bin2);
						int p2=calP_cross(nBins_GG,nBins_MM,bin4,bin3);
						for(int i=(nMaximum-nMinimum+1)*p1,n=nMinimum;i<nMaximum*(p1+1);i++,n++)
						{
							for(int j=(nMaximum-nMinimum+1)*p2,m=nMinimum; j<nMaximum*(p2+1); j++,m++)
							{
								//clog<<"n="<<n<<" m="<<m<<endl;
								//clog<<"j+nPair_GG*nMaximum="<<j+nPair_GG*nMaximum<<endl;
								CMT.load(i,j+nPair_GG*nMaximum,valueCov(n,m,integ_limits_cov)/Area_cross);
								CMT.load(j+nPair_GG*nMaximum,i,CMT.get(i,j+nPair_GG*nMaximum));
							}
						}
					}
				}
			}
		}
		//CMT.printOut((string("C_X.ascii")).c_str(),4);
	}
	clog<<"calculated covariance"<<endl;
	return CMT;
}


number psi_stats::delta(int i1, int j1)
{
	return i1==j1? 1.: 0.;
}


int psi_stats::calP(int nBins,int fbin,int sbin)
{
	if(fbin>sbin)
	{
		int swap=fbin;
		fbin=sbin;
		sbin=swap;
	}
	int p=fbin*nBins;

	for(int i=0; i<fbin; i++)
	  	p-=i;
	return p+sbin-fbin;
}

int psi_stats::calP_cross(int nBins_GG,int nBins_MM,int mbin,int gbin)
{
	return gbin*nBins_MM+mbin;
}




	