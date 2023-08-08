#include "psi_stats.h"


psi_stats::psi_stats(){}
psi_stats::~psi_stats(){}

psi_stats::psi_stats(number theta_min_arc, number theta_max_arc, number LMIN, number LMAX, 
	int LBINS, string WFolderName, string WFileName)
{ 
	initialise(theta_min_arc, theta_max_arc, LMIN, LMAX, 
		LBINS,WFolderName,WFileName);
}

void psi_stats::initialise(number theta_min_arc, number theta_max_arc, number LMIN, number LMAX, 
	int LBINS, string WFolderName1,string WFileName1)
{
	thetamin = theta_min_arc*arcmin;
	thetamax = theta_max_arc*arcmin;

	//Set paramters for the l-range which the power-spectra is integrated 
	lmin = LMIN;
	lmax = LMAX;
	lbins = LBINS;

	WFolderName = WFolderName1;
	WFileName = WFileName1;

	//Set-up the boolean to find min/max here
	find_min_max = true;
	//Wn are not set yet
	WnSet =false; 
	setIntegrationType("dummy");
}

void psi_stats::setIntegrationType(string integration_type1)
{
	integration_type = integration_type1;
}

//This method initializes the Wns
void psi_stats::setWFilters(int nMaximum1)
{
	//Set-up of global parameter describing the max no. of modes analysed
	if((!WnSet) || (nMaximum!=nMaximum1))
	{
		nMaximum = nMaximum1;

		//Initiate psi_filters object
		psi_filters Wn(thetamin, thetamax, nMaximum, lmin, 
			lmax, lbins, WFolderName, WFileName);

		Wn_vec.clear();
		
		///these need to be in separate for loops, otherwise a segmentation fault happens!
		for(int n=0; n<nMaximum; n++)
			Wn_vec.push_back(Wn);
		for(int n=0; n<nMaximum; n++)
			Wn_vec[n].set(n+1);

		WnSet=true;
	}

	if (find_min_max)
	{
		integ_limits_vec_vec.clear();

		//Loop of the modes
		for (int modeCount=0; modeCount<nMaximum; modeCount++)
		{
			nW = modeCount;
			vector<number> integ_lim_vec = determine_integration_limits(lmax);
			integ_limits_vec_vec.push_back(integ_lim_vec);
		}
		find_min_max = false;
	}
}


// Need to change this into one function:

number psi_stats::integrant(number l)
{
	if(integration_type=="gg")
	{
		number integ = Wn_vec[nW].value(l) * l * Cl.value(l);
		return integ;
	}
	else if (integration_type=="gm")
	{
		return Wn_vec[nW].value(l) * l * Cl.value(l);
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


//find extrema for W_psi to integrate over
vector<number> psi_stats::determine_integration_limits(number last_x_value)
{

	vector<number> integ_limits;

    integ_limits.push_back(lmin);

    lthresh = 2.*pi/thetamax/20.;
	integ_limits.push_back(lthresh);

	number first_logx = log(lthresh);
	number last_logx = log(last_x_value);

    //number of bins for the table of the oscillating function, 
    //you can play with this to see what effect it has on the function.
	const int Nbins = 1000000;

	// make table of the function values on a very fine grid
    vector<number> table_y(Nbins);
    vector<number> table_x(Nbins);
	
	number dx = (last_logx - first_logx)/(Nbins-1.);

	for(int i=0;i<Nbins;i++)
	{
        //log-space x-values
        table_x[i]=exp(first_logx + i*dx);
        table_y[i]=Wn_vec[nW].value(table_x[i]);
	}

	for(int i=1;i<Nbins-1;i++)
	{
		if(   ( (table_y[i-1]<table_y[i]) && (table_y[i+1]<table_y[i]) )
		   || ( (table_y[i-1]>table_y[i]) && (table_y[i+1]>table_y[i]) ) )
			integ_limits.push_back(table_x[i]);
	}

	integ_limits.push_back(last_x_value);
    return integ_limits;
}

void psi_stats::setCl(vector<number> log_ell,vector<number> Cl_in)
{
	Cl.loadWithValues(log_ell,Cl_in,true);
	Cl.extrapolationOn();
}

number psi_stats::value_psi(int n)
{

	nW=n-1;
	number result=0.;
	int steps;

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


vector<number> psi_stats::calPsi(string integration_type1)
{
	integration_type = integration_type1;

	vector<number> psi_vec;

	for(int n = 1; n<=nMaximum; n++)
		psi_vec.push_back(value_psi(n));

	return psi_vec;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// COVARIANCE section


void psi_stats::input_for_setCl(vector<number> log_ell, vector<vector<number> > cl_gg, 
	vector<vector<number> > cl_gm,vector<vector<number> > cl_mm)
{
	Cl_GG_vec.clear();
	Cl_GM_vec.clear();
	Cl_MM_vec.clear();
	nPair_GG=cl_gg.size();
	nPair_GM=cl_gm.size();
	nPair_MM=cl_mm.size();
	nBins_GG=(sqrt(8*nPair_GG+1)-1)/2;
	nBins_MM=(sqrt(8*nPair_MM+1)-1)/2;

	//these have to be separate otherwise does not work
	for(int r=0; r<nPair_GM; r++)
		Cl_GM_vec.push_back(function_cosebis());

	for(int r=0; r<nPair_GM; r++)
	{
		Cl_GM_vec[r].loadWithValues(log_ell,cl_gm[r],true);
		Cl_GM_vec[r].extrapolationOn();
	}

	for(int r=0; r<nPair_GG; r++)
		Cl_GG_vec.push_back(function_cosebis());

	for(int r=0; r<nPair_GG; r++)
	{
		Cl_GG_vec[r].loadWithValues(log_ell,cl_gg[r],true);
		Cl_GG_vec[r].extrapolationOn();
	}

	for(int r=0; r<nPair_MM; r++)
		Cl_MM_vec.push_back(function_cosebis());

	for(int r=0; r<nPair_MM; r++)
	{
		Cl_MM_vec[r].loadWithValues(log_ell,cl_mm[r],true);
		Cl_MM_vec[r].extrapolationOn();
	}
}


void psi_stats::setNoise(vector<number> sigma_e_vec,vector<number> nBar_shear_vec,vector<number> nBar_gals_vec)
{
	clog<<"setting noise in psi_stats"<<endl;
	noise_gg_vec.clear();
	noise_mm_vec.clear();

	for(int bin=0; bin<nBar_gals_vec.size(); bin++)
	{
		noise_gg_vec.push_back(1./nBar_gals_vec[bin]);
		//clog<<"bin"<<bin<<" noise gg is "<<noise_gg_vec[bin]<<endl;
	}

	for(int bin=0; bin<sigma_e_vec.size(); bin++)
	{
		noise_mm_vec.push_back(sigma_e_vec[bin]*sigma_e_vec[bin]/(2.*nBar_shear_vec[bin]));
		//clog<<"bin"<<bin<<" noise mm is "<<noise_mm_vec[bin]<<endl;
	}
}


vector<vector<matrix* > > psi_stats::determine_integration_limitsCov()
{
	// Idea: find a possibly complete list of consecutive local 
	//minima/maxima of oscillating integrant and integrate between them

	int Nbins = 1000000;
	number LHIGH=high*20./thetamax;
	number LLOW=lmin;
	number lthresh=pi/thetamax/2.;

	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
		table.load(0,i,exp(log(lthresh)+log(LHIGH/lthresh)/(Nbins-1.)*i));

	vector<vector<matrix* > > integ_limits_cov;
	vector <matrix*> integ_mat_vec;
	for(int n=0; n<nMaximum; n++)
	{
		integ_mat_vec.clear();
		for(int m=0; m<nMaximum; m++)
		{
			vector <number> integ_lim_vec;
			for(int i=0;i<Nbins;i++)
				table.load(1,i,Wn_vec[n].value(table.get(0,i))*Wn_vec[m].value(table.get(0,i)));

			integ_lim_vec.push_back(LLOW);
			integ_lim_vec.push_back(lthresh);
			for(int i=1;i<Nbins-1;i++)//2->1
				if ((table.get(1,i-1)<table.get(1,i)&& table.get(1,i+1)<table.get(1,i))
				|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
				      integ_lim_vec.push_back(table.get(0,i));
			integ_lim_vec.push_back(LHIGH);

			int size=integ_lim_vec.size();
			integ_mat_vec.push_back(new matrix(size));
			for(int i=0; i<integ_lim_vec.size(); i++)
				integ_mat_vec[m]->load(i,integ_lim_vec[i]);
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


matrix psi_stats::calCov(number Area_GG,number Area_GM,number Area_cross,bool Cross)
{
	clog<<"calculating the covariance in psi_stats"<<endl;
	matrix CMT(nMaximum*(nPair_GG+nPair_GM),nMaximum*(nPair_GG+nPair_GM));
	setWFilters(nMaximum);

	//clog<<"power spectra are set, nBins_GG="<<nBins_GG<<endl;
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
					for(int i=nMaximum*p1,n=1;i<nMaximum*(p1+1);i++,n++)
					{
						for(int j=nMaximum*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
						{
							CMT.load(i,j,valueCov(n,m,integ_limits_cov)/Area_GG);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
				}
			}
		}
	}

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
					int n,m=0;
					rp1=calP(nBins_GG,bin1,bin3);
					rp2=calP(nBins_MM,bin2,bin4);
					rp3=calP_cross(nBins_GG,nBins_MM,bin4,bin1);
					rp4=calP_cross(nBins_GG,nBins_MM,bin2,bin3);
					//this should depend on the type of correlation
					delta1noise=delta(bin1,bin3)*noise_gg_vec[bin1];
					delta2noise=delta(bin2,bin4)*noise_mm_vec[bin2];
					delta3noise=0.;
					delta4noise=0.;
					//the pair considered for the Eparam_vec
					int p1=calP_cross(nBins_GG,nBins_MM,bin2,bin1);
					int p2=calP_cross(nBins_GG,nBins_MM,bin4,bin3);
					for(int i=nMaximum*p1,n=1;i<nMaximum*(p1+1);i++,n++)
					{
						for(int j=nMaximum*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
						{
							
							int i_ggl=i+nMaximum*nPair_GG;
							int j_ggl=j+nMaximum*nPair_GG;
							CMT.load(i_ggl,j_ggl,valueCov(n,m,integ_limits_cov)/Area_GM);
							CMT.load(j_ggl,i_ggl,CMT.get(i_ggl,j_ggl));
						}
					}
				}
			}
		}
	}

	if(Cross)
	{
		cov_case="cross";
		clog<<"cov_case="<<cov_case<<endl;
		for(int bin1=0; bin1<nBins_GG; bin1++)
		{
			for(int bin2=bin1; bin2<nBins_GG; bin2++)
			{
				for(int bin3=0; bin3<nBins_GG; bin3++)
				{
					for(int bin4=0; bin4<nBins_MM;bin4++)
					{
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
						for(int i=nMaximum*p1,n=1;i<nMaximum*(p1+1);i++,n++)
						{
							for(int j=nMaximum*p2,m=1; j<nMaximum*(p2+1); j++,m++)
							{
								CMT.load(i,j+nPair_GG*nMaximum,valueCov(n,m,integ_limits_cov)/Area_cross);
								CMT.load(j+nPair_GG*nMaximum,i,CMT.get(i,j+nPair_GG*nMaximum));
							}
						}
					}
				}
			}
		}
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

	