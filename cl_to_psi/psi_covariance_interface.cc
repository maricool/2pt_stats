#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/section_names.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <chrono>
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
    double theta_min;                    // theta_min this is the minimum scale analysed
    double theta_max;                    // theta_max this is the maximum scale analysed
    int ell_bins;                        // no. of ell bins which the W-filters are calculated at
    double ell_min;                      // minimum ell-scale used in the Cl-integration
    double ell_max;                      // maximum ell-scale used in the Cl-integration
    psi_stats *FABmath;                    // Pointer to the object of the class psi_stats
    vector<vector<double> > Derivs;      // theta-bins at which the estimators were computed at
    int corr_type_key;                   // If 0 then Cl_gg+Cl_gm, if 1 then clustering + ggl, if 2 then mock analysis
    int mock_theta_bins;                 // No. of mock theta bins
    int Integ_nsteps;                    // No. of integration steps for Cl->psi Integral
    double set_lthresh;                  // If want to set l-thresh
    int nMin;
    string input_section_name_ggl;
    string input_section_name_clustering;
    string output_section_name_ggl;
    string output_section_name_clustering;
    string input_section_name_matter;
    number Area_GG;
    number Area_GM;
    number Area_cross;
    bool calCov;
    string cov_file_name;
    bool Cross_cov;
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
    options->get_val<string>(OPTION_SECTION, string("input_section_name_ggl"),"galaxy_shear_cl", config->input_section_name_ggl);
    options->get_val<string>(OPTION_SECTION, string("input_section_name_clustering"),"galaxy_cl", config->input_section_name_clustering);
    options->get_val<string>(OPTION_SECTION, string("input_section_name_matter"),"shear_cl", config->input_section_name_matter);
    options->get_val<string>(OPTION_SECTION, string("output_section_name_ggl"),"galaxy_shear_psi", config->output_section_name_ggl);
    options->get_val<string>(OPTION_SECTION, string("output_section_name_clustering"),"galaxy_psi", config->output_section_name_clustering);


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
    psi_stats *FABmath = new psi_stats(thetamin_rad, thetamax_rad, tins, ellmin, ellmax, ellbins,
                                    Nstep_Integ, set_lthresh, set_nMin,WFolderName,WFileName);
    
    //clog << "\n" << "\n" << "\n" << "\n testWFilters" << "\n" << "\n" << "\n" << endl;
    //clog << "\n" << "\n" << "\n" << "\n setWFilters" << "\n" << "\n" << "\n" << endl;
    FABmath->setWFilters(Nmodes,set_nMin);
    //clog << "\n" << "\n" << "\n" << "\n Write to config" << "\n" << "\n" << "\n" << endl;

    int calculateCov;
    status=options->get_val<int>(OPTION_SECTION, string("calculateCov"), calculateCov);
    if(status)
    {
        clog<<"calculate cov not set: "<<calculateCov<<endl;
        config->calCov=false;
    }
    else
    {
        clog<<"calculateCov set to "<<calculateCov<<endl;
        if(calculateCov)
        {
            config->calCov=true;
            status=options->get_val<string>(OPTION_SECTION, string("cov_file_name"),config->cov_file_name);
            if(status)
            {
                clog<<"Could not find cov_file_name, setting to the default values of cov.ascii"<<endl;
            }
            else
            {
                clog<<"Found cov_file_name="<<config->cov_file_name<<endl;
            }

            int Cross_cov=1;
            status=options->get_val(OPTION_SECTION, string("Cross_cov"),Cross_cov);
            if(status)
            {
                clog<<"Could not find Cross_cov, setting to the default values of 1. Going to calculate cross covariance."<<endl;
                config->Cross_cov=true;
            }
            else
            {
                clog<<"Found Cross_cov="<<Cross_cov<<endl;
                if(Cross_cov)
                    config->Cross_cov=true;
                else
                    config->Cross_cov=false;
            }

            vector<number> sigma_e,ngal_shear,ngal_position;

            status=options->get_val<std::vector<number> >(OPTION_SECTION, string("sigma_e"),sigma_e);
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
            status=options->get_val<std::vector<number> >(OPTION_SECTION, string("ngal_shear"),ngal_shear);
            if(status)
            {
                clog<<"Didn't find ngal_shear values for covariance"<<endl;
                exit(1);
            }
            else
            {
                clog<<"Found "<<ngal_shear.size()<<" ngal_shear values"<<endl;
                for(int i=0; i<ngal_shear.size(); i++)
                {
                    clog<<i<<":"<<ngal_shear[i]<<endl;
                    ngal_shear[i]*=1./arcmin/arcmin;
                }
            }
            status=options->get_val<std::vector<number> >(OPTION_SECTION, string("ngal_position"),ngal_position);
            if(status)
            {
                clog<<"Didn't find ngal_position values for covariance"<<endl;
                exit(1);
            }
            else
            {
                clog<<"Found "<<ngal_position.size()<<" ngal_position values"<<endl;
                for(int i=0; i<ngal_position.size(); i++)
                {
                    clog<<i<<":"<<ngal_position[i]<<endl;
                    ngal_position[i]*=1./arcmin/arcmin;
                }
            }
            FABmath->setNoise(sigma_e,ngal_shear,ngal_position);


            number Area;
            status=options->get_val<number>(OPTION_SECTION, string("Area"),Area);
            if(status)
            {
                clog<<"Didn't find the Area value for covariance, will look for Area_GG, Area_GM and Area_cross"<<endl;
                status=options->get_val<number>(OPTION_SECTION, string("Area_GG"),config->Area_GG);
                if(status)
                {
                    clog<<"Didn't find the Area_GG value for covariance, exiting now"<<endl;
                    exit(1);
                }
                else
                {
                    clog<<"found Area_GG="<<config->Area_GG<<endl;
                    config->Area_GG*=pow(pi/180.,2);
                }
                status=options->get_val<number>(OPTION_SECTION, string("Area_GM"),config->Area_GM);
                if(status)
                {
                    clog<<"Didn't find the Area_GM value for covariance, exiting now"<<endl;
                    exit(1);
                }
                else
                {
                    clog<<"found Area_GM="<<config->Area_GM<<endl;
                    config->Area_GM*=pow(pi/180.,2);
                }
                status=options->get_val<number>(OPTION_SECTION, string("Area_cross"),config->Area_cross);
                if(status)
                {
                    clog<<"Didn't find the Area_cross value for covariance, exiting now"<<endl;
                    exit(1);
                }
                else
                {
                    clog<<"found Area_cross="<<config->Area_cross<<endl;
                    config->Area_cross*=pow(pi/180.,2);
                }

            }
            else
            {
                clog<<"Found Area="<<Area<<endl;
                Area*=pow(pi/180.,2);//change to radians
                config->Area_GG=Area;
                config->Area_GM=Area;
                config->Area_cross=Area;
            }
            
            
        }
        else
        {
            config->calCov=false;
        }
    }

    config->FABmath=FABmath;
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
        int num_z_bin_A_ggl;
        int num_z_bin_B_ggl;
        status = block->get_val(config->input_section_name_ggl, string("nbin_a"), num_z_bin_A_ggl);
        if(status)
        {
            status = block->get_val(config->input_section_name_ggl, string("nbin"), num_z_bin_A_ggl);
            if(status)
            {
                clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
                return status;
            }
            num_z_bin_B_ggl = num_z_bin_A_ggl;
        }
        else
        {
            status = block->get_val(config->input_section_name_ggl, string("nbin_b"), num_z_bin_B_ggl);
            if(status)
            {
                clog<<"looked for nbin_b it is not set"<<endl;
                return status;
            }
        }

        vector<number> ell,logell;
        //get ell vector
        status = block->get_val(config->input_section_name_ggl, string("ell"), ell);
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

        //clustering
        int num_z_bin_A_clustering;
        int num_z_bin_B_clustering;
        status = block->get_val(config->input_section_name_clustering, string("nbin_a"), num_z_bin_A_clustering);
        if(status)
        {
            status = block->get_val(config->input_section_name_clustering, string("nbin"), num_z_bin_A_clustering);
            if(status)
            {
                clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
                return status;
            }
            num_z_bin_B_clustering = num_z_bin_A_clustering;
        }
        else
        {
            status = block->get_val(config->input_section_name_clustering, string("nbin_b"), num_z_bin_B_clustering);
            if(status)
            {
                clog<<"looked for nbin_b it is not set"<<endl;
                return status;
            }
        }


        if(config->calCov)
        {
            //clog<<"Reading ggl"<<endl;
            vector <vector<number> > InputPower_ggl_vec_vec;
            int nPairs_ggl=0;
            for (int i_bin=1; i_bin<=num_z_bin_A_ggl; i_bin++) 
            {
                for (int j_bin=1; j_bin<=num_z_bin_B_ggl; j_bin++) 
                {
                    // read in C(l)
                    vector<number> C_ell;
                    string name_in=string("bin_")+toString(i_bin)+string("_")+toString(j_bin);
                    string index_name_in=string("index_")+name_in;
                    //check if C(l) exists for this z-bin combination
                    bool has_val = block->has_val(config->input_section_name_ggl, name_in);
                    if (has_val) 
                    {
                        //clog<<name_in<<endl;
                        //if C(l) exists then read in
                        status = block->get_val<vector<number> >(config->input_section_name_ggl, name_in, C_ell);
                        InputPower_ggl_vec_vec.push_back(C_ell);
                        ///put Cl in a vector to send to COSEBIs.
                        nPairs_ggl++;
                    }
                }
            }
            //clog<<"Reading Clustering"<<endl;
            vector <vector<number> > InputPower_clustering_vec_vec;
            int nPairs_clustering=0;
            for (int i_bin=1; i_bin<=num_z_bin_A_clustering; i_bin++) 
            {
                for (int j_bin=1; j_bin<=num_z_bin_B_clustering; j_bin++) 
                {
                    // read in C(l)
                    vector<number> C_ell;
                    string name_in=string("bin_")+toString(j_bin)+string("_")+toString(i_bin);
                    string index_name_in=string("index_")+name_in;
                    //check if C(l) exists for this z-bin combination
                    bool has_val = block->has_val(config->input_section_name_clustering, name_in);
                    if (has_val) 
                    {
                        //clog<<name_in<<endl;
                        //if C(l) exists then read in
                        status = block->get_val<vector<number> >(config->input_section_name_clustering, name_in, C_ell);
                        InputPower_clustering_vec_vec.push_back(C_ell);
                        ///put Cl in a vector to send to COSEBIs.
                        nPairs_clustering++;
                    }
                }
            }
           // clog<<"clustering read, nPairs_clustering="<<nPairs_clustering<<endl;

            //now matter, needed for covariance
            int num_z_bin_A_matter;
            int num_z_bin_B_matter;
            status = block->get_val(config->input_section_name_matter, string("nbin_a"), num_z_bin_A_matter);
            if(status)
            {
                status = block->get_val(config->input_section_name_matter, string("nbin"), num_z_bin_A_matter);
                if(status)
                {
                    clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
                    return status;
                }
                num_z_bin_B_matter = num_z_bin_A_matter;
            }
            else
            {
                status = block->get_val(config->input_section_name_matter, string("nbin_b"), num_z_bin_B_matter);
                if(status)
                {
                    clog<<"looked for nbin_b it is not set"<<endl;
                    return status;
                }
            }

    //  shear_cl ---> matter power
            //clog<<"Reading matter"<<endl;
            vector <vector<number> > InputPower_matter_vec_vec;
            int nPairs_matter=0;
            for (int i_bin=1; i_bin<=num_z_bin_B_matter; i_bin++) 
            {
                for (int j_bin=1; j_bin<=num_z_bin_B_matter; j_bin++) 
                {
                    // read in C(l)
                    vector<number> C_ell;
                    string name_in=string("bin_")+toString(j_bin)+string("_")+toString(i_bin);
                    string index_name_in=string("index_")+name_in;
                    //check if C(l) exists for this z-bin combination
                    bool has_val = block->has_val(config->input_section_name_matter, name_in);
                    if (has_val) 
                    {
                        //clog<<name_in<<endl;
                        //if C(l) exists then read in
                        status = block->get_val<vector<number> >(config->input_section_name_matter, name_in, C_ell);
                        InputPower_matter_vec_vec.push_back(C_ell);
                        ///put Cl in a vector to send to COSEBIs.
                        nPairs_matter++;
                    }
                }
            }
    //        clog<<"matter read, nPairs_matter="<<nPairs_matter<<endl;

            //set the input power spectrum for psi
            //clog<<"put powers into FABmath"<<endl;
            config->FABmath->input_for_setCl(logell,InputPower_clustering_vec_vec,InputPower_ggl_vec_vec
                ,InputPower_matter_vec_vec);

            //have to set a check here that matter power is also set and all the needed combinations are there
            matrix Cov_mat=config->FABmath->calCov(config->Area_GG,config->Area_GM, config->Area_cross, config->Cross_cov);
            Cov_mat.printOut((config->cov_file_name).c_str(),20);
            // matrix icov=Cov_mat.inverse();
            // icov.printOut((config->cov_file_name+string(".inverse")).c_str(),20);
            // clog<<"covariance done"<<endl;
        }

        //clog<<"set the cls in psi"<<endl;
        vector<number> psi_all,n_vals;

        for(int n=0; n<config->nmode-nMin+1; n++)
            n_vals.push_back(n+nMin);


        // now read the clustering cl --> galaxy_cl
        int nPairs=0;
        for (int i_bin=1; i_bin<=num_z_bin_A_clustering; i_bin++) 
        {
            for (int j_bin=1; j_bin<=num_z_bin_B_clustering; j_bin++) 
            {
                //check which bins are available
                vector<number> psi;
                string name_in=string("bin_")+toString(i_bin)+string("_")+toString(j_bin);
                string index_name_in=string("index_")+name_in;
                //check if C(l) exists for this z-bin combination
                bool has_val = block->has_val(config->input_section_name_clustering, name_in);
                if (has_val) 
                {
                    //clog<<name_in<<endl;
                    for (int mode_count=nMin; mode_count<=config->nmode; mode_count++)
                    {
                        //clog<<"mode_count="<<mode_count<<endl;
                        vector<number> C_ell;
                        status = block->get_val<vector<number> >(config->input_section_name_clustering, name_in, C_ell);
                        config->FABmath->setCl(logell,C_ell);
                        val = config->FABmath->valueFunc_Cl(mode_count,"FAB_gg",nPairs);
                        //clog<<"psi_gg="<<val<<endl;
                        psi.push_back(val);
                        psi_all.push_back(val);
                    }
                    block->put_val(config->output_section_name_clustering, name_in, psi);
                    block->put_val(config->output_section_name_clustering, index_name_in, n_vals);
                    nPairs++;
                }
            }
        }
        //clog<<"done clustering"<<endl;
        block->put_val(config->output_section_name_clustering, "nbin_a", num_z_bin_A_clustering); 
        block->put_val(config->output_section_name_clustering, "nbin_b", num_z_bin_B_clustering);

// 
        nPairs=0;
        for (int i_bin=1; i_bin<=num_z_bin_A_ggl; i_bin++) 
        {
            for (int j_bin=1; j_bin<=num_z_bin_B_ggl; j_bin++) 
            {
                // read in C(l)
                vector<number> psi;
                string name_in=string("bin_")+toString(i_bin)+string("_")+toString(j_bin);
                string index_name_in=string("index_")+name_in;
                //check if C(l) exists for this z-bin combination
                bool has_val = block->has_val(config->input_section_name_ggl, name_in);
                if (has_val) 
                {
                    for (int mode_count=nMin; mode_count<=config->nmode; mode_count++)
                    {
                        vector<number> C_ell;
                        status = block->get_val<vector<number> >(config->input_section_name_ggl, name_in, C_ell);
                        config->FABmath->setCl(logell,C_ell);
                        val = config->FABmath->valueFunc_Cl(mode_count,"FAB_gm",nPairs);
                        //clog<<"psi_gm="<<val<<endl;
                        psi.push_back(val);
                        psi_all.push_back(val);
                    }
                    block->put_val(config->output_section_name_ggl, name_in, psi);
                    block->put_val(config->output_section_name_ggl, index_name_in, n_vals);
                    nPairs++;
                }
            }
        }
        block->put_val(config->output_section_name_ggl, "nbin_a", num_z_bin_A_ggl); 
        block->put_val(config->output_section_name_ggl, "nbin_b", num_z_bin_B_ggl);     
        block->put_val("data_vector", "psi_theory", psi_all);     
    }
    return status;
}

} // end of extern C
