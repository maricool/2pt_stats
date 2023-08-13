include ${COSMOSIS_SRC_DIR}/config/compilers.mk

COLIB     = modules/libcolib.a

ALL_LIBS = band_powers/band_powers.so \
           bin_xi/bin_xi.so \
           cl_to_cosebis/cl_to_cosebis.so \
           cosebis_covariance/cosebis_covariance.so \
           cl_to_psi/cl_to_psi.so \
           psi_covariance/psi_covariance.so

# ALL_LIBS = cl_to_psi/cl_to_psi.so psi_covariance/psi_covariance.so


all: $(ALL_LIBS)
     
$(COLIB): modules/*.cc modules/*.h
	cd modules && $(MAKE)


band_powers/band_powers.so: $(COLIB) band_powers/band_powers_interface.cc
	cd band_powers && $(MAKE)

bin_xi/bin_xi.so: $(COLIB) bin_xi/bin_xi_interface.cc
	cd bin_xi && $(MAKE)

cl_to_cosebis/cl_to_cosebis.so: $(COLIB) cl_to_cosebis/cl_to_cosebis_interface.cc
	cd cl_to_cosebis && $(MAKE)

cl_to_psi/cl_to_psi.so: $(COLIB) cl_to_psi/cl_to_psi_interface.cc
	cd cl_to_psi && $(MAKE)

psi_covariance/psi_covariance.so: $(COLIB) psi_covariance/psi_covariance_interface.cc
	cd psi_covariance && $(MAKE)

# xi_to_cosebis/xi_to_cosebis.so: $(COLIB) xi_to_cosebis/xi_to_cosebis_interface.cc
# 	cd xi_to_cosebis && $(MAKE)

cosebis_covariance/cosebis_covariance.so: $(COLIB) cosebis_covariance/cosebis_covariance_interface.cc
	cd cosebis_covariance && $(MAKE)

clean:
	rm -rf */*.so */*.so.dSYM
	cd modules && $(MAKE) clean

	
	
