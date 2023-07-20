#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/section_names.h"
//#include <iostream>
//#include <cmath>
//#include <stdlib.h>
//#include <vector>
//#include <string>
//#include <stdio.h>
//#include <fstream>
//#include <sstream>
//#include <string.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
//#include <chrono>
#include "psi_stats.h"
#include "matrix.h"


/*
TODO:
1) extend the WFilter integral to higher ell
2) test derivatives
*/

using namespace std;

extern "C" {


typedef struct psi_config 
{
    int nmode;                           // this is the maximum number of modes analysed
    number theta_min;                    // theta_min this is the minimum scale analysed
    number theta_max;                    // theta_max this is the maximum scale analysed
    int ell_bins;                        // no. of ell bins which the W-filters are calculated at
    number ell_min;                      // minimum ell-scale used in the Cl-integration
    number ell_max;                      // maximum ell-scale used in the Cl-integration
    psi_stats *psi;                  // Pointer to the object of the class psi_stats
    int corr_type_key;                   // If 0 then Cl_gg+Cl_gm, if 1 then clustering + ggl, if 2 then mock analysis
    int mock_theta_bins;                 // No. of mock theta bins
    int Integ_nsteps;                    // No. of integration steps for Cl->psi Integral
    number set_lthresh;                  // If want to set l-thresh
    int nMin;
    string input_section_name;
    string output_section_name;
    string type;
} 
psi_config;


///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///


// DATABLOCK_STATUS save_matrix(cosmosis::DataBlock * block,
//                              std::vector<std::vector<double> > &matrix)
// {
//     // Get the size of the matrix.
//     // This will fail if the first row of the matrix is empty
//     // Should check for that and return an error here
//     // (exercise!)
//     size_t ni = matrix.size();
//     size_t nj = matrix[0].size();

//     // We have to flatten matrix first.
//     // This is the space for that
//     std::vector<double> F; 

//     // Now loop through each row of the matrix
//     // and append it to the flattened vector.
//     for (int i=0; i<ni; i++)
//     {
//         // Get row
//         std::vector<double> &row = matrix[i];
//         // Copy row to flattened vector
//         F.insert(F.end(), row.begin(), row.end());
//     }

//     // record the dimensions of the matrix
//     std::vector<size_t> dims = {ni, nj};

//     // Make an ndarray object out of our flattened array and the dimensions
//     cosmosis::ndarray<double> M(F, dims);
    
//     // Save the ndarray to block
//     DATABLOCK_STATUS status = block->put_val("data_vector", "psi_inverse_covariance", M);

//     return status;
// }


///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///


void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
{

    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    clog<<"In setup of psi_module"<<endl;
    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

    psi_config * config = new psi_config;

    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

    //Read in the filename of the inverse covariance matrix
    //string cov_filename;
    string WFolderName,WFileName;
    double thetamin,thetamax,ellmin,ellmax,set_lthresh,step_size;
    int Nmodes,ellbins,tins,tcorr,mbins,Nstep_Integ,lbins_for_file,set_nMin;


    //read options from the ini file, there are default values for each of these
    //checkout: https://bitbucket.org/joezuntz/cosmosis/wiki/modules_cpp
    //options->get_val(OPTION_SECTION, "inverse_covariance_filename", cov_filename);
    options->get_val(OPTION_SECTION, "theta_min_arcmin",1., thetamin);
    options->get_val(OPTION_SECTION, "theta_max_arcmin",100., thetamax);
    options->get_val(OPTION_SECTION, "n_modes",10, Nmodes);
    //options->get_val(OPTION_SECTION, "l_bins",1000, ellbins); //not used currently
    options->get_val(OPTION_SECTION, "l_min",1.0, ellmin);
    options->get_val(OPTION_SECTION, "l_max",5e5, ellmax);
    //options->get_val(OPTION_SECTION, "Theta_Bins",10000000, tins); //not used currently
    options->get_val(OPTION_SECTION, "Type_of_correlation_input",1, tcorr); //for now this is the only option that works
    //options->get_val(OPTION_SECTION, "Mock_Theta_bins",10000, mbins); //not used now
    options->get_val(OPTION_SECTION, "lthresh",10., set_lthresh); //we should fix this to a value
    options->get_val(OPTION_SECTION, "n_min",1, set_nMin);
    options->get_val(OPTION_SECTION, "step_size",0.02, step_size);
    options->get_val<string>(OPTION_SECTION, string("W_output_folder_name"),"../psi/WFilters/", WFolderName);
    options->get_val<string>(OPTION_SECTION, string("W_file_name"),"W_Psi", WFileName);


    options->get_val<string>(OPTION_SECTION, string("type"),"gg", config->type);
    // set if statements here so that it automatically sets the input and output sections names based on type

    options->get_val<string>(OPTION_SECTION, string("input_section_name"),"galaxy_galaxy_cl", config->input_section_name);
    options->get_val<string>(OPTION_SECTION, string("output_section_name"),"galaxy_galaxy_psi", config->output_section_name);


    // options->get_val<string>(OPTION_SECTION, string("input_section_name_ggl"),"galaxy_shear_cl", config->input_section_name_ggl);
    // options->get_val<string>(OPTION_SECTION, string("input_section_name_clustering"),"galaxy_cl", config->input_section_name_clustering);
    // options->get_val<string>(OPTION_SECTION, string("input_section_name_matter"),"shear_cl", config->input_section_name_matter);
    // options->get_val<string>(OPTION_SECTION, string("output_section_name_ggl"),"galaxy_shear_psi", config->output_section_name_ggl);
    // options->get_val<string>(OPTION_SECTION, string("output_section_name_clustering"),"galaxy_psi", config->output_section_name_clustering);


    // make it such that you have to run it separately for ggl and clustering the same as band powers.

    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    clog<<"The Setup parameters are:"<<endl;
    clog<<endl;

    clog << "theta_min_arcmin = " << thetamin << endl;
    clog << "theta_max_arcmin = " << thetamax << endl;
    clog << "nmodes = " << Nmodes << endl;
    clog << "l_bins = " << ellbins << endl;
    clog << "l_min = " << ellmin << endl;
    clog << "l_max = " << ellmax << endl;
    //clog << "Theta_Bins = " << tins << endl;
    clog << "Type_of_correlation_input = " << tcorr << endl;
    //clog << "Mock_Theta_bins = " << mbins << endl;
    //clog << "Using the covariance file at location: " << cov_filename << endl;
    clog << "Integration_Steps = " << Nstep_Integ << endl;
    clog << "lthresh = " << set_lthresh << endl;
    clog << "nMin = " << set_nMin << endl;
    clog << "Step-size of derivative = " << step_size << endl;
    clog << "tcorr = " << tcorr << endl;


    if (tcorr==0)
        clog << "You are using Real-space inputs to compute psi" << endl;
    else if (tcorr==1)
        clog << "You are using Fourier-space inputs to compute psi" << endl;
    else
    {
        clog << "Not a recognised value, either use 0 for read space or 1 for Fourier space analysis, exiting now ..." << endl;
        exit(1);
    }

    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    clog<<endl<<endl;
    

    psi_stats FABtemp;
    double thetamin_rad,thetamax_rad;
    thetamin_rad = FABtemp.arcmin_to_radians(thetamin);
    thetamax_rad = FABtemp.arcmin_to_radians(thetamax);

    //clog << "\n" << "\n" << "\n" << "\n Initiate psi_stats" << "\n" << "\n" << "\n" << endl;
    psi_stats *psi = new psi_stats(thetamin_rad, thetamax_rad, tins, ellmin, ellmax, ellbins,
                                    Nstep_Integ, set_lthresh, set_nMin,WFolderName,WFileName);
    
    //clog << "\n" << "\n" << "\n" << "\n testWFilters" << "\n" << "\n" << "\n" << endl;
    //clog << "\n" << "\n" << "\n" << "\n setWFilters" << "\n" << "\n" << "\n" << endl;
    psi->setWFilters(Nmodes,set_nMin);
    //clog << "\n" << "\n" << "\n" << "\n Write to config" << "\n" << "\n" << "\n" << endl;

    config->psi=psi;
    config->nmode = Nmodes;
    config->theta_min = thetamin_rad;
    config->theta_max = thetamax_rad;
    config->ell_bins = ellbins;
    config->ell_min = ellmin;
    config->ell_max = ellmax;
    //config->theta_bins = tins;
    config->corr_type_key = tcorr;
    //config->mock_theta_bins = mbins;
    //config->inv_cov = cov_vecvec;
    config->Integ_nsteps = Nstep_Integ;
    config->set_lthresh = set_lthresh;
    config->nMin = set_nMin;

return (void *) config;
}


///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///


  
DATABLOCK_STATUS execute(cosmosis::DataBlock * block, void * config_in)
{
    // clog<<endl<<endl;
    // clog<<"?????????????????????????????????????????????"<<endl;
    // clog<<"-------------In Psi execute------------------"<<endl;
    // clog<<"?????????????????????????????????????????????"<<endl;
    // clog<<endl;

    ///------------------------------------------------------------------------------------------///
    ///------------------------------------------------------------------------------------------///

    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

    psi_config * config = (psi_config*) config_in;

    

    int nMin=config->nMin;      // Initial mode number of analysis
    vector<double> psi;         // Container for the psi vectors
    double val; 


    int corr_type = config->corr_type_key;

    // fix this
    ///at the start we said if it is ==0 the cl
    ///------------------------------------------------------------------------------------------///
    ///------------------------------------------------------------------------------------------///
    //If corr_type=1, then code goes from Cl->Psi-statistics
    ///------------------------------------------------------------------------------------------///
    ///------------------------------------------------------------------------------------------///

    if (corr_type==1)
    {
        //first do ggl
        //get the number of redshift bins from cosmosis
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

        vector<number> ell,logell;
        //get ell vector
        status = block->get_val(config->input_section_name, string("ell"), ell);
        if (status) 
        {
            clog<<"Could not load ell in C_ell to psi"<<endl;
            return status;
        }

        int nell=ell.size();
        //make logell to send to psi_statss
        for(int i=0; i<nell; i++)
        {
            logell.push_back(log(ell[i]));
        }


        vector<int> n_vals;
        for(int n=0; n<config->nmode-nMin+1; n++)
            n_vals.push_back(n+nMin);

        for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
        {
            for (int j_bin=1; j_bin<=num_z_bin_B; j_bin++) 
            {
                //check which bins are available
                vector<number> psi_vec;
                string name_in=string("bin_")+toString(i_bin)+string("_")+toString(j_bin);
                string index_name_in=string("index_")+name_in;
                //check if C(l) exists for this z-bin combination
                bool has_val = block->has_val(config->input_section_name, name_in);
                if (has_val) 
                {
                    //clog<<name_in<<endl;
                    for (int mode_count=nMin; mode_count<=config->nmode; mode_count++)
                    {
                        //clog<<"mode_count="<<mode_count<<endl;
                        vector<number> C_ell;
                        status = block->get_val<vector<number> >(config->input_section_name, name_in, C_ell);
                        config->psi->setCl(logell,C_ell);
                        val = config->psi->value_psi(mode_count,config->type);
                        psi_vec.push_back(val);
                    }
                    block->put_val(config->output_section_name, name_in, psi_vec);
                    block->put_val(config->output_section_name, index_name_in, n_vals);
                }
            }
        }
        block->put_val(config->output_section_name, "nbin_a", num_z_bin_A); 
        block->put_val(config->output_section_name, "nbin_b", num_z_bin_B);
        block->put_val(config->output_section_name, "type", config->type);
        block->put_val(config->output_section_name, "input_section_name", config->input_section_name);
    }
    return status;
}

} // end of extern C
