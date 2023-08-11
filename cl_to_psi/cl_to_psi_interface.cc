#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/section_names.h"
#include "psi_stats.h"


using namespace std;

extern "C" {

typedef struct psi_config 
{
    int nMaximum;                        // this is the maximum number of modes analysed
    number theta_min;                    // theta_min this is the minimum scale analysed
    number theta_max;                    // theta_max this is the maximum scale analysed
    psi_stats *psi;                      // Pointer to the object of the class psi_stats
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

    psi_config * config = new psi_config;

    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;


    string WFolderName,WFileName;
    number ellmin,ellmax;
    int ellbins;



    status = options->get_val<string>(OPTION_SECTION, string("type"),"gg", config->type);
    if (status) 
    {
        clog<<"Could not load type to psi,";
        clog<<"setting to default: gg"<<endl;
    }
    else
        clog<<"type = "<<config->type<<endl;


    status = options->get_val<string>(OPTION_SECTION, string("input_section_name"), config->input_section_name);
    if (status) 
    {
        if( (config->type) == "gg")
        {
            config->input_section_name = "galaxy_cl";
        }
        else if ( (config->type) == "gm")
        {
             config->input_section_name = "galaxy_shear_cl";
        }
        else
        {
            clog<<"WARNINING!!!!!!!! not a recognised type ="<<config->type<<" use gg or gm or set input_section_name manually"<<endl;
            clog<<"setting config->type to "<<"gg"<<endl;
            config->type = "gg";
            config->input_section_name = "galaxy_cl";
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
            config->output_section_name = "galaxy_psi";
        }
        else if ( (config->type) == "gm")
        {
             config->output_section_name = "galaxy_shear_psi";
        }
        else
        {
            clog<<"WARNINING!!!!!!!! not a recognised type ="<<config->type<<" use gg or gm or set output_section_name manually"<<endl;
            clog<<"setting config->type to "<<"gg"<<endl;
            config->type = "gg";
            config->output_section_name = "galaxy_psi";
        }

        clog<<"Could not load output_section_name to psi, used type to set its value to "<<config->output_section_name<<endl;
    }
    else
        clog<<"output_section_name = "<<config->output_section_name<<endl;



    status = options->get_val(OPTION_SECTION, "theta_min",1., config->theta_min);
    if (status) 
    {
        clog<<"Could not load theta_min to psi,";
        clog<<"setting to default: "<<config->theta_min<<endl;
    }
    else
        clog<<"theta_min (arcmin) = "<<config->theta_min<<endl;

    status = options->get_val(OPTION_SECTION, "theta_max",100., config->theta_max);
    if (status) 
    {
        clog<<"Could not load theta_max to psi,";
        clog<<"setting to default: "<<config->theta_max<<endl;
    }
    else
        clog<<"theta_max (arcmin) = "<<config->theta_max<<endl;


    status = options->get_val(OPTION_SECTION, "n_max",10, config->nMaximum);
    if (status) 
    {
        clog<<"Could not load n_max to psi,";
        clog<<"setting to default: "<<config->nMaximum<<endl;
    }
    else
        clog<<"n_max="<<config->nMaximum<<endl;


    status = options->get_val(OPTION_SECTION, "l_min",1.0, ellmin);
    if (status) 
    {
        clog<<"Could not load l_min to psi,";
        clog<<"setting to default: "<<ellmin<<endl;
    }
    else
        clog<<"l_min="<<ellmin<<endl;

    status = options->get_val(OPTION_SECTION, "l_max",1e6, ellmax);
    if (status) 
    {
        clog<<"Could not load l_max to psi,";
        clog<<"setting to default: "<<ellmax<<endl;
    }
    else
        clog<<"l_max="<<ellmax<<endl;

    status = options->get_val(OPTION_SECTION, "l_bins",1000, ellbins);
    if (status) 
    {
        clog<<"Could not load l_bins to psi,";
        clog<<"setting to default: "<<ellbins<<endl;
    }
    else
        clog<<"l_bins="<<ellbins<<endl;


    string WnFolderName;
    WnFolderName = COSEBIS_DIR  "Wn_psi/";

    options->get_val<string>(OPTION_SECTION, string("W_output_folder_name"), WFolderName);
    options->get_val<string>(OPTION_SECTION, string("W_file_name"),"W_psi", WFileName);

    clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

    psi_stats *psi = new psi_stats(config->theta_min, config->theta_max, ellmin, ellmax, ellbins,WFolderName,WFileName);
    
    psi->setWFilters(config->nMaximum);

    config->psi=psi;

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
    enable_gsl_error_handling();

    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
    psi_config * config = (psi_config*) config_in;

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


    //get ell vector
    vector<number> ell,logell;
    status = block->get_val(config->input_section_name, string("ell"), ell);
    if (status) 
    {
        clog<<"Could not load ell in C_ell to psi"<<endl;
        return status;
    }

    int nell=ell.size();
    //make logell to send to psi_stats
    for(int i=0; i<nell; i++)
    {
        logell.push_back(log(ell[i]));
    }


    vector<int> n_vals;
    for(int n=0; n<config->nMaximum; n++)
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
                vector<number> C_ell;
                status = block->get_val< vector<number> >(config->input_section_name, name_in, C_ell);
                config->psi->setCl(logell,C_ell);
                vector<number> psi_vec = config->psi->calPsi(config->type);

                block->put_val(config->output_section_name, name_in, psi_vec);
                block->put_val(config->output_section_name, index_name_in, n_vals);
            }
        }
    }
    block->put_val(config->output_section_name, "nbin_a", num_z_bin_A); 
    block->put_val(config->output_section_name, "nbin_b", num_z_bin_B);
    block->put_val(config->output_section_name, "theta_min", config->theta_min);
    block->put_val(config->output_section_name, "theta_max", config->theta_max);
    block->put_val(config->output_section_name, "type", config->type);
    block->put_val(config->output_section_name, "input_section_name", config->input_section_name);

    disable_gsl_error_handling();
    return status;
}

} // end of extern C
