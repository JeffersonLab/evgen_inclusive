include Makefile.arch

#------------------------------------------------------------------------------

EICRATES  := $(wildcard *.$(SrcSuf)) eicRateDict.$(SrcSuf)
EICRATEH  := $(wildcard *.h)
EICRATEO  := $(EICRATES:.$(SrcSuf)=.$(ObjSuf)) 
EICRATE    = eicRate$(ExeSuf)

OBJS          = $(EICRATEO)

PROGRAMS      = $(EICRATE)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

$(EICRATE):  $(EICRATEO)
		@echo "this: $@, $(EICRATEO)"
		$(LD) $(LDFLAGS) $^ $(GLIBS) -O2 -lEG -lRGL -lGed -lHtml -lThread -lGeom -lcteqpdf $(OutPutOpt)$@
		cp eicRate output/
		@echo "$@ done"

eicRate.$(SrcSuf): 
		@echo "Generating dictionary $@..."
		@rootcint -f $@ -c $^

clean:
		@rm -f $(OBJS) *.res core

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.geom *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

###

.$(SrcSuf).$(ObjSuf): 
		$(CXX) $(CXXFLAGS) -c $<

eicRateDict.$(SrcSuf): LinkDef.h
		@echo "Generating dictionary $@..."
		$(ROOTCINT) -f $@ -c $^

