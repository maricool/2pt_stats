///for cpp use the .hh and for c the .h version
//This deals with the inputs and outputs

#include "datablock/datablock.hh"
#include "datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to E/B - COSEBIs
*/
#include "COSEBIs.h"


extern "C" 
{
	const string shear_cl = SHEAR_CL_SECTION;
	const number MinPowerCOSEBIs=0.1;
	const number MaxPowerCOSEBIs=1e6;
	const int PowerTableNumberCOSEBIs=200;

	typedef struct COSEBIs_config {
		string sectionName;
		string input_section_name;
		string output_section_name;
		int n_max;
		number theta_min;
		number theta_max;
		bool calNoiseCov;
		string Cov_En_name;
		int nBins; //this is only needed if the noise only covariance is to be estimated from input nPair
		COSEBIs *cosebis;
		number sigma_m;
		matrix sigma_m_cov;
		bool sigma_m_cov_read;
		bool cal_nonGaussian_cov;//if true calculates the non-Gaussian covariance
		matrix InputCl_Cov_mat; //needs and input Cl matrix with ordering: bin_1_1: ell1,ell2,... bin_2_1: ell1, ell2,...
		vector<number> Input_ell_vec; // The input ell_vec that corresponds to the given input Cl matrix
	} COSEBIs_config;

	///define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as cosebis. COSEBIs is a class that produces En/Bn and 
	//their covariance, etc. 

	int get_option(cosmosis::DataBlock * options, const string &name, string &parameter)
	{
	
		auto status = options->get_val(OPTION_SECTION, name, parameter);
		if (status!=DBS_SUCCESS) 
		{
			parameter = "";
			cerr<< "Could not find or understand parameter in cosebis section: " 
				<< name << std::endl; 
			return 1;
		}
		return 0;
	}


  	void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
  	{

  		COSEBIs_config * config = new COSEBIs_config;

  		string sectionName=OPTION_SECTION;
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		clog<<endl<<endl;
		clog<<"*********in COSEBIs covariance interface setup*********"<<endl;
		clog<<endl;

		
		status=options->get_val<string>(sectionName, string("output_section_name"), config->output_section_name);
		if (status) 
		{
			clog<<"Could not load out_section_name to COSEBIs, ";
			clog<<"setting to default: cosebis"<<endl;
			config->output_section_name=string("cosebis");
		}
		else
			clog<<"Got the value of output_section_name:"<<config->output_section_name<<endl;

		//get input section name, default= shear_cl
		status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
		if (status) 
		{
			clog<<"Could not load input_section_name to COSEBIs,";
			clog<<" setting to default: shear_cl"<<endl;
			config->input_section_name=string("shear_cl");
		}
		else
			clog<<"Got the value of input_section_name:"<<config->input_section_name<<endl;
	

  		status=options->get_val<number>(sectionName, string("theta_min"), 1. , config->theta_min);
  		if (status) 
		{
			clog<<"Could not load theta_min to COSEBIs"<<endl;
			clog<<"Setting it to the default value of "<<1.<<endl;
		}
		else
			clog<<"Got the value of theta_min="<<config->theta_min<<endl;

	    status=options->get_val<number>(sectionName, string("theta_max"),100., config->theta_max);

	    if (status) 
	    {
			clog<<"Could not load theta_max to COSEBIs"<<endl;
			clog<<"Setting it to the default value of "<<100.<<endl;
	    }
		else
			clog<<"Got the value of theta_max="<<config->theta_max<<endl;


	    status=options->get_val<int>(sectionName, string("n_max"),10, config->n_max);
	    if (status) 
	    {
			clog<<"Could not load n_max to COSEBIs"<<endl;
			clog<<"setting it to the default value of "<<10<<endl;
	    }
		else
			clog<<"Got the value of n_max="<<config->n_max<<endl;

		//get Wn, Tn and output Tn folder names
		string WnFolderName,TnFolderName,OutputTnFolderName,WnFileName;
		WnFolderName= COSEBIS_DIR "WnLog/";
		WnFileName   = "WnLog";
		TnFolderName= COSEBIS_DIR "TLogsRootsAndNorms/";
		OutputTnFolderName=COSEBIS_DIR "TpnLog/";

		// int precision = 20;
		// int Nlbins    = 1000000;

		int precision = 10;
		int Nlbins    = 100000;

		status=options->get_val<string>(sectionName, string("Wn_Output_FolderName"), WnFolderName);
		if(status)
		{
			clog<<"Could not find WnLog folder name in Wn_Output_FolderName, ";
			clog<<"setting to default: "<<WnFolderName<<endl;
		}
		else
			clog<<"WnLog folder name is:"<<WnFolderName<<endl;

		status=options->get_val<string>(sectionName, string("Wn_file_name"), WnFileName);
		if(status)
		{
			clog<<"Could not find Wn file name in Wn_file_name, ";
			clog<<"setting to default: "<<WnFileName<<endl;
		}
		else
			clog<<"Wn file name is:"<<WnFileName<<endl;

		status=options->get_val<int>(sectionName, string("table_precision"), precision);
		if(status)
		{
			clog<<"Could not find the number of digits in table_precision, ";
			clog<<"setting to default: "<<precision<<endl;
		}
		else
			clog<<"table_precision is:"<<precision<<endl;

		status=options->get_val<int>(sectionName, string("number_of_Wn_l_bins"), Nlbins);
		if(status)
		{
			clog<<"Could not find number_of_Wn_l_bins, ";
			clog<<"setting to default: "<<Nlbins<<endl;
		}
		else
			clog<<"number_of_Wn_l_bins is:"<<Nlbins<<endl;

		status=options->get_val<string>(sectionName, string("Roots_n_Norms_FolderName"), TnFolderName);
		if(status)
		{
			clog<<"Could not find Root and Norms folder name in Roots_n_Norms_FolderName,"; 
			clog<<" setting to default: "<<TnFolderName<<endl;
		}
		else
			clog<<"Root and Norms folder name is:"<<TnFolderName<<endl;


		status=options->get_val<string>(sectionName, string("Tn_Output_FolderName"), OutputTnFolderName);
		if(status)
		{
			clog<<"Could not find T_pm folder name in Tn_Output_FolderName, ";
			clog<<"setting to default: "<<OutputTnFolderName<<endl;
	  	}
		else
			clog<<"T_pm folder name is:"<<OutputTnFolderName<<endl;

		//   initialize COSEBIs
		COSEBIs *cosebis = new COSEBIs();
		cosebis->initialize(config->n_max,config->theta_min,config->theta_max,1 //npair set to one for now, will be set seperately in execute to the correct value
				,WnFolderName,TnFolderName,OutputTnFolderName,WnFileName,precision);

		cosebis->setWns(config->n_max,Nlbins);

		config->Cov_En_name="EnCovarianceFromTheory";
		status=options->get_val<string>(sectionName, string("cov_name"), config->Cov_En_name);
		clog<<"cov name start is set to: "<<config->Cov_En_name<<endl;

		//
		string sigma_e_file;
		vector<number> sigma_e;
		// first read sigma_e if not given exit
		status=options->get_val<string>(sectionName, string("sigma_e_file"),sigma_e_file);
		if(status)
		{
			clog<<"Didn't find sigma_e_file for covariance, going to look for sigma_e values"<<endl;
		
			status=options->get_val<std::vector<number> >(sectionName, string("sigma_e"),sigma_e);
			if(status)
			{
				clog<<"Didn't find sigma_e values for covariance"<<endl;
				exit(1);
		  	}
			else
			{
				clog<<"Found "<<sigma_e.size()<<" sigma_e values"<<endl;
				for(int i=0; i<sigma_e.size(); i++)
				{
					clog<<i<<":"<<sigma_e[i]<<endl;
					sigma_e[i]*=sqrt(2.);
				}
			}
		}
		else
		{
			clog<<"found sigma_e_file:"<<sigma_e_file<<endl;
			matrix sigma_e_mat;
			sigma_e_mat.readFromASCII_marika((sigma_e_file).c_str());
			for(int i=0; i<sigma_e_mat.size(); i++)
			{
				clog<<i<<":"<<sigma_e_mat.get(i)<<endl;
				sigma_e.push_back(sqrt(2.)*sigma_e_mat.get(i));
			}
		}

		vector<number> ngal_effective;
		string ngal_file;
		status=options->get_val<string>(sectionName, string("ngal_file"),ngal_file);
		if(status)
		{
			clog<<"Didn't find ngal_file for covariance, going to look for ngal_effective values"<<endl;
		
			status=options->get_val<std::vector<number> >(sectionName, string("ngal_effective"),ngal_effective);
			if(status)
			{
				clog<<"Didn't find ngal_effective values for covariance"<<endl;
				exit(1);
		  	}
			else
			{
				clog<<"Found "<<ngal_effective.size()<<" ngal_effective values"<<endl;
				for(int i=0; i<ngal_effective.size(); i++)
				{
					clog<<i<<":"<<ngal_effective[i]<<endl;
					ngal_effective[i]*=1./arcmin/arcmin;
				}
			}
		}
		else
		{
			clog<<"found ngal_file:"<<ngal_file<<endl;
			matrix ngal_mat;
			ngal_mat.readFromASCII_marika((ngal_file).c_str());
			for(int i=0; i<ngal_mat.size(); i++)
			{
				clog<<i<<":"<<ngal_mat.get(i)<<endl;
				ngal_effective.push_back(ngal_mat.get(i)/arcmin/arcmin);
			}
		}

		number Area;
		status=options->get_val<number>(sectionName, string("Area"),Area);
		if(status)
		{
			clog<<"Didn't find Area values for covariance"<<endl;
			exit(1);
	  	}
		else
		{
			clog<<"Found Area="<<Area<<endl;
			Area*=pow(pi/180.,2);//change to radians
		}

		cosebis->setNoise(Area,sigma_e,ngal_effective);

		string inputCl_cov;
		config->cal_nonGaussian_cov=false;//if true calculates the non-Gaussian covariance
		status=options->get_val<string>(sectionName, string("input_nonGaussian_Cl_cov"),inputCl_cov);
		if(status)
		{
			clog<<"Did not find an input nonGaussian term"<<endl;
		}
		else
		{
			clog<<"got the input nonGaussian Cl covariance name:"<<inputCl_cov<<endl;
			config->InputCl_Cov_mat.readFromASCII_marika(inputCl_cov.c_str());

			string input_ell;
			status=options->get_val<string>(sectionName, string("input_nonGaussian_Cl_ell_vec"),input_ell);
			if(status)
			{
				clog<<"please also give a file with a vector of corresponding ell values for the non-Gaussian covariance"<<endl;
				clog<<"Going to skip calculating non-Gaussian covariance"<<endl;
			}
			else
			{
				config->cal_nonGaussian_cov=true;
				clog<<"Got the ell_vec file. Going to calculate the non-Gaussian terms for COSEBIs Cov"<<endl;
				matrix Input_ell_mat;
				Input_ell_mat.readFromASCII_marika(input_ell.c_str());
				for(int i=0; i<Input_ell_mat.rows; i++)
				{
					config->Input_ell_vec.push_back(Input_ell_mat.get(i));
				}
				clog<<"initialised Input_ell_vec with "<<config->Input_ell_vec.size()<<" ell bins"<<endl;
			}
		}

		config->sigma_m=0.;
		status=options->get_val<number>(sectionName, string("sigma_m"),config->sigma_m);
		if(status)
		{
			clog<<"Sigma_m is not set."<<endl;
		}
		else
		{
			clog<<"Got the value of sigma_m="<<config->sigma_m<<endl;
		}

		string sigma_m_cov_file;
		config->sigma_m_cov_read=false;
		status=options->get_val<string>(sectionName, string("sigma_m_cov_file"),sigma_m_cov_file);
		if(status)
		{
			clog<<"No covariance file given for sigma_m."<<endl;
		}
		else
		{
			config->sigma_m_cov_read=true;
			clog<<"got the value of sigma_m_cov_file="<<sigma_m_cov_file<<endl;
			config->sigma_m_cov.readFromASCII_marika((sigma_m_cov_file).c_str());
		}

		status=options->get_val(sectionName, string("nBins"), config->nBins);
	    if (status) 
	    {
			clog<<"Could not load nBins to COSEBIs."<<endl;
	    }
		else
			clog<<"Got the value of nBins="<<config->nBins<<endl;

		string input_nPair_files_prefix;
		status=options->get_val(sectionName, string("input_nPair_files_prefix"),input_nPair_files_prefix);
		if(status)
		{
			config->calNoiseCov=false;
			clog<<"No input_nPair_files_preffix was given."<<endl;
		}
		else
		{
			config->calNoiseCov=true;
			clog<<input_nPair_files_prefix<<" is the input nPair prefix."<<endl;
			string input_nPair_files_suffix="";
			status=options->get_val(sectionName, string("input_nPair_files_suffix"),input_nPair_files_suffix);
			int Athena_input=1;
			bool Athena= true;
			status=options->get_val(sectionName, string("Athena_input"),Athena_input);
			if(Athena_input)
			{
				clog<<"Assuming input is athena, so going to devide the noise contribution for autocorrelated z-bins by 2"<<endl;
				Athena=true;
			}
			else
			{
				clog<<"Assuming input is not athena. Maybe it is treecorr. No correction to the autocorrelated z-bins is done."<<endl;
				Athena=false;
			}
			int nCol_nPair=8;//default for Athena
			int nCol_theta=1;
			status=options->get_val(sectionName, string("nPairs_column"),nCol_nPair);
			status=options->get_val(sectionName, string("theta_column"),nCol_theta);
			if(config->nBins)
			{
				vector<string> FileName_vec;
				for(int bin1=0; bin1<config->nBins; bin1++)
				{
					for(int bin2=bin1; bin2<config->nBins; bin2++)
					{
						string FileName=input_nPair_files_prefix+
							+("_nBins_")+toString(config->nBins)+string("_Bin")
							+toString(bin1+1)+string("_Bin")+toString(bin2+1)+input_nPair_files_suffix;
						FileName_vec.push_back(FileName);
					}
				}
				//-1 because c++ starts from 0
				cosebis->readNpairs(FileName_vec,nCol_theta-1,nCol_nPair-1,Athena);
			}
		}
		config->cosebis=cosebis;

		clog<<endl<<endl;
		clog<<"*********end of COSEBIs covariance interface setup*********"<<endl;
		clog<<endl;

  		return (void *) config;
  		// config is sent to execute 
	}

	DATABLOCK_STATUS execute(cosmosis::DataBlock *block, void *config_in) 
	{
		// Config is whatever you returned from setup above
		// Block is the collection of parameters and calculations for
		// this set of cosmological parameters
		
		COSEBIs_config *config= (COSEBIs_config*) config_in;
		
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		//get cl from cosmosis
		int num_z_bin_A;
		int num_z_bin_B;
		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			if(status)
			{
				clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
				return status;
			}
			num_z_bin_B = num_z_bin_A;
		}
		else
		{
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
			if(status)
			{
				clog<<"looked for nbin_b it is not set"<<endl;
				return status;
			}
		}

		// get ell from cosmosis and turn into ln(ell)
		vector<number> ell,logell;
		status = block->get_val(config->input_section_name, string("ell"), ell);
		int nell=ell.size();
		for(int i=0; i<nell; i++)
			logell.push_back(log(ell[i]));

		
		if (status) 
		{
			clog<<"Could not load ell in C_ell to COSEBIs"<<endl;
			return status;
		}
		

		// Put Cl a vector of vectors to be sent to cosebis
		vector <vector<number> > InputPower_vec_vec;
		int nPairs=0;
		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=1; j_bin<=num_z_bin_B; j_bin++) 
			{
				// read in C(l)
				vector<number> C_ell;
				string name_in=string("bin_")+toString(j_bin)+string("_")+toString(i_bin);
				bool has_val = block->has_val(config->input_section_name, name_in);
				if (has_val) 
				{
					status = block->get_val<vector<number> >(config->input_section_name, name_in, C_ell);
					InputPower_vec_vec.push_back(C_ell);
					nPairs++;
				}
			}
		}

		//put cl and nPairs into cosebis
		clog<<"nPairs="<<nPairs<<endl;
		config->cosebis->setZbins(nPairs);
		config->cosebis->setPower(logell,InputPower_vec_vec);

		// Calculate the En vector, this is saved later and also used for getting the sigma_m covariance
		matrix En_mat;
		En_mat=config->cosebis->calEn();

		//set these to zero
		matrix Cov_sigm,Cov_nG_th;
		Cov_sigm.zero(0);
		Cov_nG_th.zero(0);

		if(config->sigma_m)
		{
			Cov_sigm=config->cosebis->calCovForSigma_m(config->sigma_m);
			Cov_sigm.printOut((config->Cov_En_name+string("_sigma_m_")+toString(config->sigma_m,4)+(".ascii")).c_str(),20);
		}

		if(config->sigma_m_cov_read)
		{
			Cov_sigm=config->cosebis->calCovForSigma_m_from_m_cov(config->sigma_m_cov);
			Cov_sigm.printOut((config->Cov_En_name+string("_sigma_m_from_m_cov_")+toString(config->sigma_m,4)+(".ascii")).c_str(),20);
		}

		if(config->cal_nonGaussian_cov)
		{
			Cov_nG_th=config->cosebis->calCovFromInputPowerCov(config->InputCl_Cov_mat, config->Input_ell_vec);
			Cov_nG_th.printOut((config->Cov_En_name+string("_non_Gaussian")+(".ascii")).c_str(),20);
		}

		matrix Cov_Bn,CovNoise,Cov_En,Cov_cosmicVar,Cov_mixed;

		// calculate and save theory Bn covariance. This assumes simple geometry.
		Cov_Bn=config->cosebis->calBCov();
		Cov_Bn.printOut((config->Cov_En_name+string("_TheoryBn.ascii")).c_str(),20);

		// This section calculates Bn covariance using the input number of galaxy pairs. This is also used as the noise term for En covariance.
		if(config->calNoiseCov)
		{
			if((int((config->nBins)*(config->nBins+1))/2)==nPairs)
			{
				CovNoise=config->cosebis->calNoiseCov_fromInputNpair();
				CovNoise.printOut((config->Cov_En_name+string("_NoiseOnly.ascii")).c_str(),20);				
			}
			else
			{
				clog<<"!!!!WARNING!!!! the number of bins for the noise only case does not match the input number of redshift bins"<<endl;
			}
		}

		// Calculate the En theory Gaussian covariance assuming simple geometry.
		Cov_En=config->cosebis->calCov();
		Cov_En.printOut((config->Cov_En_name+string("_TheoryEn.ascii")).c_str(),20);

		// Calculates the Cosmic variance term for En covariance
		config->cosebis->setNoiseToZero();
		Cov_cosmicVar=config->cosebis->calCov();
		Cov_cosmicVar.printOut((config->Cov_En_name+string("_TheoryCosmicVar.ascii")).c_str(),20);

		// Calculates the Mixed term (assuming simple geometry) for En covariance
		Cov_mixed=Cov_En-Cov_Bn-Cov_cosmicVar;
		Cov_mixed.printOut((config->Cov_En_name+string("_TheoryMixed.ascii")).c_str(),20);

		// Calculates the Gaussian covariance including the correct geometry effectes for the noise-only term.
		if(config->calNoiseCov)
		{
			if((int((config->nBins)*(config->nBins+1))/2)==nPairs)
			{
				Cov_En=Cov_cosmicVar+Cov_mixed+CovNoise;
				Cov_En.printOut((config->Cov_En_name+string("_NoiseJustForNoise.ascii")).c_str(),20);
			}
			else
				clog<<"!!!!WARNING!!!! the number of bins for the noise only case does not match the input number of redshift bins"<<endl;
		}

		// Adds the non-Gaussian term to the covariance calculated in the previous step
		if(config->cal_nonGaussian_cov)
			Cov_En=Cov_En+Cov_nG_th;

		// Adds the sigma_m covaraince to the covariance calculated in the previous step
		if( (config->sigma_m) || (config->sigma_m_cov_read) )
			Cov_En=Cov_En+Cov_sigm;

		// prints out the main result
		Cov_En.printOut((config->Cov_En_name+string(".ascii")).c_str(),20);

		int n_max=config->n_max;

		// puts the En for the parameters used to calculate the covariance matrix in block
		vector<number> En_vec(n_max);
		vector<int> n_vals(n_max);
		int p1=0;
		string name_En;
		for(int i_bin=0; i_bin<num_z_bin_A; i_bin++) 
		{
			for (int j_bin=i_bin; j_bin<num_z_bin_B; j_bin++) 
			{
				int m=0;
				for(int n1=n_max*p1,m=0 ;n1<n_max*(p1+1) ;n1++,m++)
					En_vec[m]=En_mat.get(n1);
				name_En=string("bin_")+toString(j_bin+1)+string("_")+toString(i_bin+1);
				status = block->put_val<vector<double> >(config->output_section_name, name_En, En_vec);
				p1++;
			}
		}

		status = block->put_val<int>(config->output_section_name, string("n_mode"), n_max);
		status = block->put_val<double>(config->output_section_name, string("theta_min"), config->theta_min);
		status = block->put_val<double>(config->output_section_name, string("theta_max"), config->theta_max);
		status = block->put_val<double>(config->output_section_name, string("nbin_a"), num_z_bin_A);
		status = block->put_val<double>(config->output_section_name, string("nbin_b"), num_z_bin_B);
		
		for(int n=0;n<n_max;n++)
			n_vals[n]=n+1;

	    status = block->put_val<vector<int> >(config->output_section_name, string("n"), n_vals);

	    return status;
	}
}// end of extern C

