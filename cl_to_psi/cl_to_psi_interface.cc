#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/section_names.h"
#include "psi_stats.h"
#include "matrix.h"


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
    psi_stats *psi;                      // Pointer to the object of the class psi_stats
    int Integ_nsteps;                    // No. of integration steps for Cl->psi Integral
    number set_lthresh;                  // If want to set l-thresh
    string input_section_name;           // Input section name where Cl resides.
    string output_section_name;          // Section name for outputs
    string type;                         // type of integration: either gg or gm
} 
psi_config;


int get_option(cosmosis::DataBlock * options, const string &name, string &parameter)
{

    auto status = options->get_val(OPTION_SECTION, name, parameter);
    if (status!=DBS_SUCCESS) 
    {
        parameter = "";
        cerr<< "Could not find or understand parameter in psi section: " 
            << name << std::endl; 
        return 1;
    }
    return 0;
}



void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
{

    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    clog<<"In setup of cl_to_psi_interface"<<endl;
    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

    psi_config * config = new psi_config;

    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;


    string WFolderName,WFileName;
    number thetamin,thetamax,ellmin,ellmax,set_lthresh,step_size;
    int n_modes,ellbins,tins,tcorr,mbins,Nstep_Integ,lbins_for_file;



    status = options->get_val<string>(OPTION_SECTION, string("type"),"gg", config->type);
    if (status) 
    {
        clog<<"Could not load type to psi,";
        clog<<"setting to default: gg"<<endl;
    }
    else
        clog<<"got the value of type:"<<config->type<<endl;


    status = options->get_val<string>(OPTION_SECTION, string("input_section_name"), config->input_section_name);
    if (status) 
    {
        if( (config->type) == "gg")
        {
            config->input_section_name = "galaxy_galaxy_cl";
        }
        else if ( (config->type) == "gm")
        {
             config->input_section_name = "galaxy_matter_cl";
        }
        else
        {
            clog<<"WARNINING!!!!!!!! not a recognised type ="<<config->type<<" use gg or gm or set input_section_name manually"<<endl;
            clog<<"setting config->type to "<<"gg"<<endl;
            config->type = "gg";
        }

        clog<<"Could not load input_section_name to psi, used type to set its value to "<<config->input_section_name<<endl;
    }
    else
        clog<<"input_section_name = "<<config->input_section_name<<endl;

    status = options->get_val<string>(OPTION_SECTION, string("output_section_name"), config->output_section_name);
    if (status) 
    {
        if ( (config->type) == "gg")
        {
            config->output_section_name = "galaxy_galaxy_psi";
        }
        else if ( (config->type) == "gm")
        {
             config->output_section_name = "galaxy_matter_psi";
        }
        else
        {
            clog<<"WARNINING!!!!!!!! not a recognised type ="<<config->type<<" use gg or gm or set output_section_name manually"<<endl;
            clog<<"setting config->type to "<<"gg"<<endl;
            config->type = "gg";
        }

        clog<<"Could not load output_section_name to psi, used type to set its value to "<<config->output_section_name<<endl;
    }
    else
        clog<<"output_section_name = "<<config->output_section_name<<endl;



    status = options->get_val(OPTION_SECTION, "theta_min",1., thetamin);
    if (status) 
    {
        clog<<"Could not load theta_min to psi,";
        clog<<"setting to default: "<<thetamin<<endl;
    }
    else
        clog<<"theta_min (arcmin) = "<<thetamin<<endl;

    status = options->get_val(OPTION_SECTION, "theta_max",100., thetamax);
    if (status) 
    {
        clog<<"Could not load theta_max to psi,";
        clog<<"setting to default: "<<thetamin<<endl;
    }
    else
        clog<<"theta_max (arcmin) = "<<thetamax<<endl;


    status = options->get_val(OPTION_SECTION, "n_modes",10, n_modes);
    if (status) 
    {
        clog<<"Could not load n_modes to psi,";
        clog<<"setting to default: "<<thetamin<<endl;
    }
    else
        clog<<"n_modes="<<n_modes<<endl;


    status = options->get_val(OPTION_SECTION, "l_min",1.0, ellmin);
    if (status) 
    {
        clog<<"Could not load l_min to psi,";
        clog<<"setting to default: "<<ellmin<<endl;
    }
    else
        clog<<"l_min="<<ellmin<<endl;

    status = options->get_val(OPTION_SECTION, "l_max",1.0, ellmax);
    if (status) 
    {
        clog<<"Could not load l_max to psi,";
        clog<<"setting to default: "<<ellmax<<endl;
    }
    else
        clog<<"l_max="<<ellmax<<endl;


    string WnFolderName;
    WnFolderName = COSEBIS_DIR  "Wn_psi/";


    options->get_val(OPTION_SECTION, "lthresh",10., set_lthresh); //we should fix this to a value
    options->get_val(OPTION_SECTION, "step_size",0.02, step_size);
    options->get_val<string>(OPTION_SECTION, string("W_output_folder_name"), WFolderName);
    options->get_val<string>(OPTION_SECTION, string("W_file_name"),"W_psi", WFileName);

    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    clog<<"The Setup parameters are:"<<endl;
    clog<<endl;

    clog << "theta_min = " << thetamin << endl;
    clog << "theta_max = " << thetamax << endl;
    clog << "nmodes = " << n_modes << endl;
    clog << "l_bins = " << ellbins << endl;
    clog << "l_min = " << ellmin << endl;
    clog << "l_max = " << ellmax << endl;
    clog << "Integration_Steps = " << Nstep_Integ << endl;
    clog << "lthresh = " << set_lthresh << endl;

    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    clog<<endl<<endl;

    psi_stats *psi = new psi_stats(thetamin, thetamax, tins, ellmin, ellmax, ellbins,
                                    Nstep_Integ, set_lthresh,WFolderName,WFileName);
    
    psi->setWFilters(n_modes);

    config->psi=psi;
    config->nmode = n_modes;
    config->theta_min = thetamin;
    config->theta_max = thetamax;
    config->ell_bins = ellbins;
    config->ell_min = ellmin;
    config->ell_max = ellmax;
    config->Integ_nsteps = Nstep_Integ;
    config->set_lthresh = set_lthresh;

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

    vector<double> psi;         // Container for the psi vectors
    double val; 


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
    for(int n=0; n<config->nmode+1; n++)
        n_vals.push_back(n+1);

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
                for (int mode_count=0; mode_count<=config->nmode; mode_count++)
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
    return status;
}

} // end of extern C
