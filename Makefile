include Makefile.arch

#------------------------------------------------------------------------------

EICRATES  := $(wildcard *.$(SrcSuf)) eicRateDict.$(SrcSuf)
EICRATEH  := $(wildcard *.h)
EICRATEO  := $(EICRATES:.$(SrcSuf)=.$(ObjSuf)) 
EICRATE    = eicRate$(ExeSuf)

OBJS          = $(EICRATEO)

PROGRAMS      = $(EICRATE)

FC=gfortran
FFLAGS= -O2

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(EICRATE):  $(EICRATEO) piajp.o
		@echo "this: $@, $(EICRATEO)"
		$(LD) $(LDFLAGS) $^ $(LIBS) -ggdb -lgfortran -O2 -lEG -lGed -lHtml -lThread -lGeom -L. -lcteqpdf $(OutPutOpt)$@
		@echo "$@ done"

eicRate.$(SrcSuf): 
		@echo "Generating dictionary $@..."
		@rootcint -f $@ -c $^

piajp.o : piajp.f
		$(FC) -c $(FFLAGS) $< -o $@

clean:
		@rm -f $(OBJS) *.res core piajp.o

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.geom *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

###

.$(SrcSuf).$(ObjSuf): 
		$(CXX) $(CXXFLAGS) -c -g $<

eicRateDict.$(SrcSuf): LinkDef.h
		@echo "Generating dictionary $@..."
		$(ROOTCINT) -f $@ -c $^

