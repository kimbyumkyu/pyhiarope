PROGRAM       = pythia_rope

version       = JTKT
CXX           = g++
#CXXFLAGS      = -O -Wall -g -Wno-deprecated -bind_at_load -D$(version)
CXXFLAGS      = -O -Wall -g -Wno-deprecated -D$(version)  #-ggdb
LD            = g++
LDFLAGS       = -O 
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)
#CXXFLAGS += $(shell $(FASTJET)/fastjet-config --cxxflags )
#LDFLAGS += $(shell $(FASTJET)/fastjet-config --libs --plugins ) 
LIBDIRARCH      = lib
LDFLAGS += -L$(PYTHIA8)/$(LIBDIRARCH) -lpythia8 -ldl
INCS    += -I$(PYTHIA8)/include
CXXFLAGS  += $(INCS) 

HDRSDICT = 
           
HDRS	+= $(HDRSDICT)   nanoDict.h


SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS) $(PROGRAM).C
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -lEG -lPhysics -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)
		@echo "done"

%.cxx:


clean:
		rm -f $(OBJS) core *Dict* $(PROGRAM).o *.d $(PROGRAM) $(PROGRAM).sl

cl:  clean $(PROGRAM)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -c -D$(version) $(HDRSDICT)
