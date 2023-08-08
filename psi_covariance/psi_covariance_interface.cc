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
    string output_section_name;          // Section name for outputs
    string input_section_name_gg;
    string input_section_name_gm;
    string input_section_name_mm;
    number Area_GG;
    number Area_GM;
    number Area_cross;
    vector<number> ngal_shear;
    vector<number> ngal_position;
    vector<number> sigma_e;
    string cov_file_name;
    //string type;
    bool Cross_cov;
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

///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///
///------------------------------------------------------------------------------------------///


void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
{
    psi_config * config = new psi_config;

    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

    //Read in the filename of the inverse covariance matrix
    string cov_filename;
    string WFolderName,WFileName;
    number ellmin,ellmax;
    int ellbins;

    status = options->get_val<string>(OPTION_SECTION, string("input_section_name_gg"), config->input_section_name_gg);
    if (status) 
    {

        config->input_section_name_gg = "galaxy_cl";
    }
    else
        clog<<"input_section_name_gg = "<<config->input_section_name_gg<<endl;

    status = options->get_val<string>(OPTION_SECTION, string("input_section_name_gm"), config->input_section_name_gm);
    if (status) 
    {

        config->input_section_name_gm = "galaxy_shear_cl";
    }
    else
        clog<<"input_section_name_gm = "<<config->input_section_name_gm<<endl;

    status = options->get_val<string>(OPTION_SECTION, string("input_section_name_mm"), config->input_section_name_mm);
    if (status) 
    {

        config->input_section_name_mm = "shear_cl";
    }
    else
        clog<<"input_section_name_mm = "<<config->input_section_name_mm<<endl;

    status = options->get_val<string>(OPTION_SECTION, string("output_section_name"),"psi_covariance", config->output_section_name);
    if (status) 
    {

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

    status=options->get_val<string>(OPTION_SECTION, string("cov_file_name"),config->cov_file_name);
    if(status)
    {
        config->cov_file_name = COSEBIS_DIR "/psi_cov.ascii";
        clog<<"Could not find cov_file_name, setting to the default values ="<< config->cov_file_name<<endl;
    }
    else
    {
        clog<<"cov_file_name="<<config->cov_file_name<<endl;
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
        clog<<"Setting Cross_cov to ";
        if(Cross_cov)
        {
            config->Cross_cov=true;
            clog<<"true"<<endl;
        }
        else
        {
            config->Cross_cov=false;
            clog<<"false"<<endl;
        }
    }

    status=options->get_val< vector<number> >(OPTION_SECTION, string("sigma_e"),config->sigma_e);
    if(status)
    {
        clog<<"Didn't find sigma_e values for covariance"<<endl;
        exit(1);
    }
    else
    {
        clog<<"Found "<<config->sigma_e.size()<<" sigma_e values"<<endl;
        for(int i=0; i<config->sigma_e.size(); i++)
        {
            clog<<i<<":"<<config->sigma_e[i]<<endl;
            config->sigma_e[i]*=sqrt(2.);
        }
    }

    status=options->get_val<std::vector<number> >(OPTION_SECTION, string("ngal_shear"),config->ngal_shear);
    if(status)
    {
        clog<<"Didn't find ngal_shear values for covariance"<<endl;
        exit(1);
    }
    else
    {
        clog<<"Found "<<config->ngal_shear.size()<<" ngal_shear values"<<endl;
        for(int i=0; i<config->ngal_shear.size(); i++)
        {
            clog<<i<<":"<<config->ngal_shear[i]<<endl;
            config->ngal_shear[i]*=1./arcmin/arcmin;
        }
    }

    status=options->get_val<std::vector<number> >(OPTION_SECTION, string("ngal_position"),config->ngal_position);
    if(status)
    {
        clog<<"Didn't find ngal_position values for covariance"<<endl;
        exit(1);
    }
    else
    {
        clog<<"Found "<<config->ngal_position.size()<<" ngal_position values"<<endl;
        for(int i=0; i<config->ngal_position.size(); i++)
        {
            clog<<i<<":"<<config->ngal_position[i]<<endl;
            config->ngal_position[i]*=1./arcmin/arcmin;
        }
    }

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
            clog<<"Area_GG="<<config->Area_GG<<" degress squared"<<endl;
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
            clog<<"found Area_GM="<<config->Area_GM<<" degress squared"<<endl;
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
            clog<<"found Area_cross="<<config->Area_cross<<" degress squared"<<endl;
            config->Area_cross*=pow(pi/180.,2);
        }

    }
    else
    {
        clog<<"Found Area="<<Area<<" degress squared"<<endl;
        Area*=pow(pi/180.,2);//change to radians
        config->Area_GG=Area;
        config->Area_GM=Area;
        config->Area_cross=Area;
    }


    ///// setup psi for covariace calculation
    psi_stats *psi = new psi_stats(config->theta_min, config->theta_max, ellmin, ellmax, ellbins,WFolderName,WFileName);
    
    psi->setWFilters(config->nMaximum);
    psi->setNoise(config->sigma_e,config->ngal_shear,config->ngal_position);

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

    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

    psi_config * config = (psi_config*) config_in;


    ////////////////////////////////////////////////////////////////////////////////////////
    int num_z_bin_A_gm;
    int num_z_bin_B_gm;
    status = block->get_val(config->input_section_name_gm, string("nbin_a"), num_z_bin_A_gm);
    if(status)
    {
        status = block->get_val(config->input_section_name_gm, string("nbin"), num_z_bin_A_gm);
        if(status)
        {
            clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
            return status;
        }
        num_z_bin_B_gm = num_z_bin_A_gm;
    }
    else
    {
        status = block->get_val(config->input_section_name_gm, string("nbin_b"), num_z_bin_B_gm);
        if(status)
        {
            clog<<"looked for nbin_b it is not set"<<endl;
            return status;
        }
    }

    // we are assuming that the ell bining is the same for all three.
    vector<number> ell,logell;
    //get ell vector
    status = block->get_val(config->input_section_name_gm, string("ell"), ell);
    if (status) 
    {
        clog<<"Could not load ell in C_ell to psi"<<endl;
        return status;
    }

    int nell=ell.size();
    //make logell to send to psi
    for(int i=0; i<nell; i++)
    {
        logell.push_back(log(ell[i]));
    }

    // reading galaxy-shear Cls
    vector <vector<number> > InputPower_gm_vec_vec;
    for (int i_bin=1; i_bin<=num_z_bin_A_gm; i_bin++) 
    {
        for (int j_bin=1; j_bin<=num_z_bin_B_gm; j_bin++) 
        {
            // read in C(l)
            vector<number> C_ell;
            string name_in=string("bin_")+toString(i_bin)+string("_")+toString(j_bin);
            //check if C(l) exists for this z-bin combination
            bool has_val = block->has_val(config->input_section_name_gm, name_in);
            if (has_val) 
            {
                status = block->get_val<vector<number> >(config->input_section_name_gm, name_in, C_ell);
                InputPower_gm_vec_vec.push_back(C_ell);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    //clustering
    int num_z_bin_A_gg;
    int num_z_bin_B_gg;
    status = block->get_val(config->input_section_name_gg, string("nbin_a"), num_z_bin_A_gg);
    if(status)
    {
        status = block->get_val(config->input_section_name_gg, string("nbin"), num_z_bin_A_gg);
        if(status)
        {
            clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
            return status;
        }
        num_z_bin_B_gg = num_z_bin_A_gg;
    }
    else
    {
        status = block->get_val(config->input_section_name_gg, string("nbin_b"), num_z_bin_B_gg);
        if(status)
        {
            clog<<"looked for nbin_b it is not set"<<endl;
            return status;
        }
    }
            
    // read in galaxy-galaxy Cl
    vector <vector<number> > InputPower_gg_vec_vec;
    for (int i_bin=1; i_bin<=num_z_bin_A_gg; i_bin++) 
    {
        for (int j_bin=1; j_bin<=num_z_bin_B_gg; j_bin++) 
        {
            vector<number> C_ell;
            string name_in=string("bin_")+toString(j_bin)+string("_")+toString(i_bin);

            bool has_val = block->has_val(config->input_section_name_gg, name_in);
            if (has_val) 
            {
                status = block->get_val<vector<number> >(config->input_section_name_gg, name_in, C_ell);
                InputPower_gg_vec_vec.push_back(C_ell);
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    //read in matter-matter Cl
    int num_z_bin_A_mm;
    int num_z_bin_B_mm;
    status = block->get_val(config->input_section_name_mm, string("nbin_a"), num_z_bin_A_mm);
    if(status)
    {
        status = block->get_val(config->input_section_name_mm, string("nbin"), num_z_bin_A_mm);
        if(status)
        {
            clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
            return status;
        }
        num_z_bin_B_mm = num_z_bin_A_mm;
    }
    else
    {
        status = block->get_val(config->input_section_name_mm, string("nbin_b"), num_z_bin_B_mm);
        if(status)
        {
            clog<<"looked for nbin_b it is not set"<<endl;
            return status;
        }
    }

    vector <vector<number> > InputPower_mm_vec_vec;
    // int nPairs_matter=0;
    for (int i_bin=1; i_bin<=num_z_bin_B_mm; i_bin++) 
    {
        for (int j_bin=1; j_bin<=num_z_bin_B_mm; j_bin++) 
        {
            // read in C(l)
            vector<number> C_ell;
            string name_in=string("bin_")+toString(j_bin)+string("_")+toString(i_bin);
            bool has_val = block->has_val(config->input_section_name_mm, name_in);
            if (has_val) 
            {
                status = block->get_val<vector<number> >(config->input_section_name_mm, name_in, C_ell);
                InputPower_mm_vec_vec.push_back(C_ell);
                // nPairs_matter++;
            }
        }
    }

    // we are assuming that the ell bining is the same for all three.
    config->psi->input_for_setCl(logell,InputPower_gg_vec_vec,InputPower_gm_vec_vec
        ,InputPower_mm_vec_vec);

    //have to set a check here that matter power is also set and all the needed combinations are there
    matrix Cov_mat=config->psi->calCov(config->Area_GG,config->Area_GM, config->Area_cross, config->Cross_cov);
    Cov_mat.printOut((config->cov_file_name).c_str(),20);

    status = block->put_val(config->output_section_name, "cov_file_name", config->cov_file_name);
    status = block->put_val(config->output_section_name, "nbin_a_gg", num_z_bin_A_gg); 
    status = block->put_val(config->output_section_name, "nbin_b_gg", num_z_bin_B_gg);

    status = block->put_val(config->output_section_name, "nbin_a_gm", num_z_bin_A_gm); 
    status = block->put_val(config->output_section_name, "nbin_b_gm", num_z_bin_B_gm);

    status = block->put_val(config->output_section_name, "nbin_a_mm", num_z_bin_A_mm); 
    status = block->put_val(config->output_section_name, "nbin_b_mm", num_z_bin_B_mm);

    status = block->put_val(config->output_section_name, "Cross_cov", config->Cross_cov);
    status = block->put_val(config->output_section_name, "Area_gg (deg^2)", config->Area_GG/pow(pi/180.,2));
    status = block->put_val(config->output_section_name, "Area_gm (deg^2)", config->Area_GM/pow(pi/180.,2));
    status = block->put_val(config->output_section_name, "Area_cross (deg^2)", config->Area_cross/pow(pi/180.,2));

    status = block->put_val(config->output_section_name, "input_section_name_gg", config->input_section_name_gg);
    status = block->put_val(config->output_section_name, "input_section_name_gm", config->input_section_name_gm);
    status = block->put_val(config->output_section_name, "input_section_name_mm", config->input_section_name_mm);

    status = block->put_val(config->output_section_name, "theta_min", config->theta_min);
    status = block->put_val(config->output_section_name, "theta_max", config->theta_max);
    status = block->put_val(config->output_section_name, "nMaximum", config->nMaximum);

    status = block->put_val(config->output_section_name, "sigma_e", config->sigma_e);
    status = block->put_val(config->output_section_name, "ngal_shear", config->ngal_shear);
    status = block->put_val(config->output_section_name, "ngal_position", config->ngal_position);

    return status;
}

} // end of extern C
