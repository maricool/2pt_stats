include ${COSMOSIS_SRC_DIR}/config/compilers.mk
#-----------------------------------
#--- objects to build
#-----------------------------------
USER_LDFLAGS += -L${GSL_LIB}  -lgsl -lgslcblas
USER_CXXFLAGS += -I${GSL_INC}

OBJS     = matrix.o function_cosebis.o ticker.o tostring.o tolower.o timestamp.o \
		   globaldef.o Integrate.o COSEBIs.o WnLog.o TpnLogRadian.o \
		   BandPower_g.o BandPower_W.o BandPower.o

LIBNAME  = libcolib.a

HEADER   = $(OBJS:.o=.h) 

#-----------------------------------
#--- libraries
#-----------------------------------

# make objects (components) of library
%.o : %.cc 
	$(CXX) $(USER_LDFLAGS) $(CXXFLAGS) ${USER_CXXFLAGS} -c $<
#%.o: %.cc
#	$(CXX) -o $@ -c $(CXXFLAGS) $*.cc

# make library
all: $(HEADER) $(OBJS) 
	rm -f $(LIBNAME)
	ar rv $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

.PHONY: clean
clean:
	rm -f $(LIBNAME)
	rm -f *.o

