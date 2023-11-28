EXE_NAME=run/shallow_water
FC=gfortran
ALL_OBJS = `find . -name "*.o"`
NCINCLUDES=-I/util/academic/netcdf/4.3.3.1/include 
NETCDF=-L/util/academic/netcdf/4.3.3.1/lib -lnetcdff -lnetcdf

export NCINCLUDES
export FC
export NETCDF

OBJS   = shallow_water.o initialize.o

all: $(OBJS) net_cdf dyn_cor
	$(FC) $(FFLAGS) -o $(EXE_NAME) $(ALL_OBJS) $(NETCDF)

dyn_cor:
	cd dyncor; $(MAKE) all

net_cdf:
	cd netcdf; $(MAKE) all

clean: 
	rm $(EXE_NAME)
	rm $(ALL_OBJS)
