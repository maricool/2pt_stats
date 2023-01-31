include ${COSMOSIS_SRC_DIR}/config/compilers.mk

COLIB     = modules/libcolib.a

ALL_LIBS = bin_cl/bin_cl_interface.so \
           bin_xi/bin_xi_interface.so \
           cl_to_cosebis/cl_to_cosebis.so \
           xi_to_cosebis/xi_to_cosebis.so \
           cosebis_covariance/cosebis_covariance.so

all: $(ALL_LIBS)
     
$(COLIB): modules/*.cc modules/*.h
	cd modules && $(MAKE)


bin_cl/bin_cl_interface.so: $(COLIB) bin_cl/bin_cl_interface.cc
	cd bin_cl && $(MAKE)

bin_xi/bin_xi_interface.so: $(COLIB) bin_xi/bin_xi_interface.cc
	cd bin_xi && $(MAKE)

cl_to_cosebis/cl_to_cosebis.so: $(COLIB) cl_to_cosebis/cl_to_cosebis_interface.cc
	cd cl_to_cosebis && $(MAKE)

xi_to_cosebis/xi_to_cosebis.so: $(COLIB) xi_to_cosebis/xi_to_cosebis_interface.cc
	cd xi_to_cosebis && $(MAKE)

cosebis_covariance/cosebis_covariance.so: $(COLIB) cosebis_covariance/cosebis_covariance_interface.cc
	cd cosebis_covariance && $(MAKE)

clean:
	rm -rf */*.so */*.so.dSYM
	cd modules && $(MAKE) clean



# $(submodules): %: %/%_interface.cc
# 	cd $@ && $(MAKE)

# #-----------------------------------
# #--- rules
# #-----------------------------------
# all: $(MODULE)/$(COLIB) libcosebis.so libcosebis_cov.so libbandpower.so libxipm_binned.so libcosebis_2pcfs.so 

# # tells how to make my library (another makefile)

# $(MODULE)/$(COLIB)::
# 	cd $(MODULE) && $(MAKE)



# *.so: COSEBIs_2PCFs_interface.cc
# 	$(CXX) $(CXXFLAGS) COSEBIs_2PCFs_interface.cc -shared -o libcosebis_2pcfs.so $(LDFLAGS) $(USER_LDFLAGS)

# libcosebis.so: COSEBIs_interface.cc
# 	$(CXX) $(CXXFLAGS) COSEBIs_interface.cc -shared -o libcosebis.so $(LDFLAGS) $(USER_LDFLAGS)

# libcosebis_cov.so: COSEBIs_covariance_interface.cc
# 	$(CXX) $(CXXFLAGS) COSEBIs_covariance_interface.cc -shared -o libcosebis_cov.so $(LDFLAGS) $(USER_LDFLAGS)

# libbandpower.so: BandPower_interface.cc
# 	$(CXX) $(CXXFLAGS) BandPower_interface.cc -shared -o libbandpower.so $(LDFLAGS) $(USER_LDFLAGS)

# libxipm_binned.so: xipm_binned.cc
# 	$(CXX) $(CXXFLAGS) xipm_binned.cc -shared -o libxipm_binned.so $(LDFLAGS) $(USER_LDFLAGS)


# # Clean up
# .PHONY: clean
# clean:
# 	cd $(MODULE) && $(MAKE) clean
# 	rm -f libcosebis.so
# 	rm -f libcosebis_cov.so
# 	rm -f libcosebis_2pcfs.so
# 	rm -f libbandpower.so
# 	rm -f libxipm_binned.so
# 	test -n "./" && rm -rf ./libcosebis_2pcfs.so.dSYM/
# 	test -n "./" && rm -rf ./libcosebis.so.dSYM/
# 	test -n "./" && rm -rf ./libcosebis_cov.so.dSYM/
# 	test -n "./" && rm -rf ./libbandpower.so.dSYM/
# 	test -n "./" && rm -rf ./libxipm_binned.so.dSYM/

	
	
