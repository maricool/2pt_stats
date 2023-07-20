include ${COSMOSIS_SRC_DIR}/config/compilers.mk

COLIB     = modules/libcolib.a

ALL_LIBS = bandpowers/bandpowers_interface.so \
           bin_xi/bin_xi_interface.so \
           cl_to_cosebis/cl_to_cosebis.so \
           xi_to_cosebis/xi_to_cosebis.so \
           cosebis_covariance/cosebis_covariance.so


all: $(ALL_LIBS)
     
$(COLIB): modules/*.cc modules/*.h
	cd modules && $(MAKE)


bandpowers/bandpowers_interface.so: $(COLIB) bandpowers/bandpowers_interface.cc
	cd bandpowers && $(MAKE)

bin_xi/bin_xi_interface.so: $(COLIB) bin_xi/bin_xi_interface.cc
	cd bin_xi && $(MAKE)

cl_to_cosebis/cl_to_cosebis.so: $(COLIB) cl_to_cosebis/cl_to_cosebis_interface.cc
	cd cl_to_cosebis && $(MAKE)

# xi_to_cosebis/xi_to_cosebis.so: $(COLIB) xi_to_cosebis/xi_to_cosebis_interface.cc
# 	cd xi_to_cosebis && $(MAKE)

cosebis_covariance/cosebis_covariance.so: $(COLIB) cosebis_covariance/cosebis_covariance_interface.cc
	cd cosebis_covariance && $(MAKE)

clean:
	rm -rf */*.so */*.so.dSYM
	cd modules && $(MAKE) clean

	
	
