NCSRCS = ncinit.f ncwritehist.f ncreadhist.f

all: $(NCSRCS)
	$(FC) $(NCINCLUDES) -c $(NCSRCS) 

clean:
	rm -rf *.o
