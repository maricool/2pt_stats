#include "psi_filters.h"


/* This class creates the filter functions which are used in conjunction with the 2D-power spectra to determine the FAB-statistics */


psi_filters::psi_filters(){}
psi_filters::~psi_filters(){}

//------------------------------------------------------------------------//
//------------------------------Setup-------------------------------------//
//------------------------------------------------------------------------//

psi_filters::psi_filters(number minTH, number maxTH, int Bins, int nMaximum,
	number LMIN, number LMAX, int LBINS,string FolderName,string WFileName)
{
	initialise(minTH, maxTH, Bins, nMaximum,LMIN, LMAX, LBINS,FolderName,WFileName);
}

void psi_filters::initialise(number minTH, number maxTH, int Bins, int nMaximum,
	number LMIN, number LMAX, int LBINS,string FolderName,string WFileName)
{
	clog << "Initializing psi_filters" << endl;
	setValues(minTH, maxTH, Bins, LMIN, LMAX, LBINS);
	initialize_mode(nMaximum);
	setNames(FolderName,WFileName);
}


void psi_filters::setValues(number minTH, number maxTH, int Bins, number LMIN, 
	number LMAX, int LBINS)
{ 

	//NB. All the angles within the C-codes are in radians
	thetaMin   = minTH;
	thetaMax   = maxTH;
	thetaBar   = (thetaMax + thetaMin) / 2.;
	deltaTheta = thetaMax - thetaMin;
	thetaBins  = Bins;

	//Set paramters for the l-range which the power-spectra is integrated over (i.e an approximation from l=0 --> l=infinity in reality, but because of filter functions this is not necessary)
	lmin = LMIN;
	lmax = LMAX;
	lbins = LBINS;
	clog<<"lmin="<<lmin<<" lmax="<<lmax<<" lbins="<<lbins<<endl;
}

void psi_filters::initialize_mode(int nMaximum1)
{ 
	nMaximum=nMaximum1; 
	mode=1;
}

void psi_filters::setNames(string FolderName1,string WFileName1)
{
	FolderName=FolderName1;
	WFileName=WFileName1;
}

//------------------------------------------------------------------------//
//---------------------------Get Filters----------------------------------//
//------------------------------------------------------------------------//

number psi_filters::get( number l )
{
	//I have left this corrType part in as it gives the user an option to go via the Q 
	//(ie._gm integration path although not required)
	integration_type="W_gg";
	return valueFuncW(l,thetaMax);
}

// void psi_filters::loadWZeros(string filename)
// {
// 	matrix zeros;
// 	zeros.readFromASCII_marika(filename.c_str());
// 	int rows= zeros.rows;
// 	clog<<"loaded "<<rows <<" zeros"<<endl;
// 	//vector<number> BesselJ0_Zeros(cols);
// 	BesselJ0_Zeros.clear();
// 	for(int i=0; i<rows; i++)
// 	{
// 		BesselJ0_Zeros.push_back(zeros.get(i));
// 	}
// }

//This sets up the root vector
void psi_filters::roots_zeros()
{
	all_roots.clear();
	int root_loop = 0;
	number root=0;

	//int sizeBesselJ0=sizeof(BesselJ0_Zeros);
	clog<<"sizeBesselJ0_zeros="<<sizeBesselJ0_zeros<<endl;
	//Loop through all the J0 Bessel function roots
	int max=sizeBesselJ0_zeros-1;
	clog<<"theta_max*lmax="<<thetaMax*lmax<<" BesselJ0_Zeros["<<max<<"]="<<BesselJ0_Zeros[max]<<endl;
	while ((root < thetaMax*lmax)&&(root_loop<sizeBesselJ0_zeros))
	{
		root = BesselJ0_Zeros[root_loop];
		//Check if root inside range of interest
		//First root will probably always be 2.404825557695773 as always very small value
		if (root >= thetaMin*lmin)
		{
			//clog << "no. = " << root_loop << ", root = " << root << ", max_val = " << thetaMax*lmax << endl; 
			all_roots.push_back(root);
		}
		root_loop++;
	}

	clog <<"Number of roots found = " << all_roots.size() << endl;
	//clog <<"BesselJ0_Zeros[99999]="<<BesselJ0_Zeros[99999]<<endl;
}


/*
void psi_filters::StepFinder()
{
	//StepN is a vector of consequitive zeros of integrant for each n

	clog<<"StepFinder begins"<<endl;
	int i0=0;
	int iR=0;
	number arg=0.;
	StepN.clear();
	while((arg<thetaMax*lmax)&&(i0<100000)&&(iR<100000))
	{
		if(iR>n)
		{
			arg=BesselJ0_Zeros[i0];
			i0++;
		}
		else
		{
			if(BesselJ0_Zeros[i0]<rootTheta[n-1][iR])
			{
				arg=BesselJ0_Zeros[i0];
				i0++;
			}
			else
			{
				arg=rootTheta[n-1][iR];
				iR++;
				clog<<"iR="<<iR<<endl;
			}
		}

		StepN.push_back(arg);

	}
	clog<<"StepFinder ended"<<endl;
	clog<<"i0="<<i0<<"iR="<<iR<<endl;

}*/


number psi_filters::valueFuncW(number l, number x)
{	
	lmode = l;
	//clog << "l = " << lmode << endl;
	int accuracyG = 40;
	number xmin = (lmode*thetaMin);
	number xmax = (lmode*thetaMax);
	number resultG=0.;
	int i=0;

	//clog<<"lmode="<<lmode<<endl;
	if (all_roots[0]>xmax)
	{
		resultG=gaussianIntegrate_gsl(*this, xmin, xmax,accuracyG);
 		//clog<< "# in if (StepN[0]>xmax) resultG="<<resultG<<endl;
	}
	else
	{
		for(i=0; xmin>(all_roots[i]); i++);

		number result1=gaussianIntegrate_gsl(*this, xmin, all_roots[i],accuracyG);
		resultG+=result1;
		//clog<< "# in else (StepN[0]<xmax) resultG="<<resultG<<endl;
		while(all_roots[i+1]<xmax)
		{
			number result1=gaussianIntegrate_gsl(*this, all_roots[i],all_roots[i+1],accuracyG);
			resultG+=result1;
			i++;
		}
		resultG+=gaussianIntegrate_gsl(*this, all_roots[i],xmax,accuracyG);
		//clog<< "# end of while resultG="<<resultG<<endl;
	}
	return resultG/lmode/lmode;
}


//------------------------------------------------------------------------//
//---------------------------Integrants-----------------------------------//
//------------------------------------------------------------------------//

number psi_filters::integrant(number x)
{
	// Note Wgm and Wgg should be exactly the same.
	if (integration_type=="Wgg")
		return ( x * gsl_sf_bessel_Jn(0,x) * U(mode,x/lmode) );	
	else if(integration_type=="Wgm")
		return ( x * gsl_sf_bessel_Jn(2,lmode*x) * analyticalQ(x,mode) );
	else if(integration_type=="Q")
		return ( x * U(mode,x) );
	else
	{
		clog<<"In psi_filters, not a recognised intergration type:"<<integration_type<<"exiting now ..."<<endl;
		exit(1);
	}
}

//------------------------------------------------------------------------//
//----------------------------Calulate U----------------------------------//
//------------------------------------------------------------------------//

number psi_filters::U(int n, number theta)
{
	//This if statement does the same job as the heaviside function.
	if ((theta > thetaMax)||(theta < thetaMin))
		return 0.;
	if (n==1)
		return 12. * deltaTheta * (theta - deltaTheta) / ( pow(deltaTheta, 3) * sqrt( pow(deltaTheta,2) + 24. * pow (thetaBar,2)));
	else
	{

		number prefactor = 1./(pow(deltaTheta,2))*sqrt(( 2. * n + 1.) / 2.);

		number x = 2.*(theta - thetaBar)/deltaTheta;

		//This condition has been put in purely because for some reason at the 15th decimal place 
		//stuff was getting put in and breaching the gsl domain for legendre polynomials
		//if (x<-1){x=-1.0;}

		number legendre = gsl_sf_legendre_Pl(n, x);
		return prefactor*legendre;

	}
}

//------------------------------------------------------------------------//
//----------------------------Calulate Q----------------------------------//
//------------------------------------------------------------------------//

number psi_filters::analyticalQ(number theta, int n)
{
	mode = n;
	number sum=0; 
	number numerator, denominator, factor, integration_term, integration_term_low,Q1;
	int M = floor(n/2);

	number thetaTerm = theta - thetaBar ;
	number thetaTermmin = thetaMin - thetaBar;

	if (n == 1)
	{

		number A_theta = pow(theta,2)*( (4.* theta * pow(thetaBar,2) - 6. * deltaTheta)/deltaTheta - 0.5) ;
		number A_thetaMin = pow(thetaMin,2)*( (4.* thetaMin * pow(thetaBar,2) - 6. * deltaTheta)/deltaTheta - 0.5) ;

		number prefactor = 2./(pow(theta,2) * deltaTheta * sqrt(2*deltaTheta+ 24.* pow(theta,2)));
		number Q1 = prefactor * (A_theta - A_thetaMin) - U(n,theta);
		return Q1;
	}
	else
	{
		for (int m=0; m<=M; m++)
		{
			numerator = pow(-1,m) * factorial(2*n-2*m) * pow((2./deltaTheta),(n-2*m));

			denominator = pow(2,n) * factorial(m) * factorial(n-m) * factorial(n-2*m);

			factor =  numerator / denominator;

			number nm = (n - 2*m + 1);
			integration_term = (1./nm) * ( theta * pow(thetaTerm,nm) - (pow(thetaTerm,(nm+1))/(nm+1)));

			integration_term_low = (1./nm) * ( thetaMin * pow(thetaTermmin,nm) - (pow(thetaTermmin,(nm+1))/(nm+1)));

			sum += (integration_term - integration_term_low) * factor ;
		}
		return  sqrt((2*n + 1)/2.)/pow((theta * deltaTheta), 2) * sum - U(n,theta);
	}
}

number psi_filters::factorial(number num_in)
{
	number factorial_out=1.;
	for (int i=1; i<=num_in; i++)
		factorial_out *= i;
	return factorial_out;
}

number psi_filters::radians_to_arcmin(number theta)
{
	return (60.*360.)*theta/(2.*pi);
}


//------------------------------------------------------------------------//
//--------------------------Set the table---------------------------------//
//------------------------------------------------------------------------//


void psi_filters::set(int n, string WnFileName)
{
	/*
	This function first checks if the WFilters file exists.
	If it does, then function cosebis just reads in file.
	*/
	mode = n;

	clog << "\n \n In  psi_filters::set, setting up the W filters" << endl;
	clog << "Theta min (arcmin): " << radians_to_arcmin(thetaMin) << ", Theta max (arcmin): " <<  radians_to_arcmin(thetaMax) << endl;
	clog << "lmin=" << lmin << ", lmax=" << lmax <<  ", lbins=" << lbins << endl;

	// is there a table on disk?
	string myname =FolderName+string("/")+string(WFileName)+toString(mode)+string("-")
		 +toString(radians_to_arcmin(thetaMin),2)+string("-")
		 +toString(radians_to_arcmin(thetaMax),2)+string("-lmin-")
		 +toString(lmin,1)+string("-lmax-")
		 +toString(lmax,1)+string("-lbins-")+toString(lbins);
	//this sets the name of the table to be saved or loaded from disk, it is a function in function_cosebis
	setName(myname.c_str(),function_cosebis::NONAMECOUNTER);

	// checks if the table already exists, otherwise run the stepFinder in preparation  for the integration
	ifstream fhandler((myname+string(".table")).c_str());
	if (fhandler.fail())
	{
		if(!CheckFolderExist(FolderName))
		{
			clog<<"making folder for W_psi:"<<FolderName<<endl;
			mkdir((FolderName).c_str(), 0777);
		}
		clog<<"writing table:"<<myname<<endl;
		//loadWZeros();
		roots_zeros();
		//need to set integration limits here
	}
	//close the file
	fhandler.close();

	///loads the Wn_table either from file or if the file is not available it makes it and saves to file, function_cosebis calls the get function to make the table
	//clog<<"lmin="<<lmin<<" lmax="<<lmax<<" lbins="<<lbins<<endl;
	loadTable(lmin,lmax,lbins,true);
	clog << "Done mode  " << mode << endl;
	///this sets off extrapolation for Wn_table 
	extrapolationOff();
}

