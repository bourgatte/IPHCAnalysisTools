include ../make.inc
#################################################
### Defining Compiling variables              ###
#################################################



CXX           = g++ -Wno-deprecated -Wall -std=c++11
LD            = g++ -Wno-deprecated -Wall -std=c++11
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

override CXXFLAGS += -I$(ROOTSYS)/include -I./ -c
LDFLAGS       =  -fPIC
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint  -lPostscript -lMatrix -lPhysics -lMinuit2 -lGui -LObj -lThread -rdynamic -Wl,--rpath $(ROOTSYS)/lib

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

SRCS    = $(wildcard *.cxx)
HDRS    = $(SRCS:.cxx=.h)

OBJS          = $(SRCS:.cxx=.o) 

LIB     = lib/libDataFormats.so

OBJS          = MyDict.o

SHAREDLIB     = lib/libDataFormats.so

$(SHAREDLIB): $(OBJS)
	@echo "======================================================="
	@echo "Linking SharedLib: $(SHAREDLIB) ..."
	@echo $(LD) $(LDFLAGS) -shared $(OBJS) $(LIBS) -o $(SHAREDLIB)
	@$(LD) $(LDFLAGS) -shared $(OBJS) $(LIBS) -o $(SHAREDLIB)
	@echo "Linking SharedLib: $(SHAREDLIB) Complete"
	@echo "======================================================="

##########################################################
###                                                    ###
##########################################################

ROOTCINT: $(HDRS)  Linkdef.h
	 rootcint -f MyDict.cxx -c  $^

$(OBJS): MyDict.cxx 
	$(LD) -shared -rdynamic  -o $@ `root-config --ldflags` $(CXXFLAGS) $(LDFLAGS) $(ROOTLIBS) $^

.PHONY: clean install 

install: ROOTCINT $(OBJS) $(SHAREDLIB)

clean:
	@rm MyDict.*
	@rm lib/*.so
