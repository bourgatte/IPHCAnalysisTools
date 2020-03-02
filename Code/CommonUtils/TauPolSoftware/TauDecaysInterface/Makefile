#################################################
### Defining Compiling variables              ###
#################################################



CXX           = g++ -Wno-deprecated -Wall  -std=c++11 #-I$(ROOTSYS)/include -pthread -std=c++17 -m64 -I/opt/sbg/cms/safe1/cms/bourgatte/Tau/Temp/include -I../../ -c
LD            = g++ -Wno-deprecated -Wall  -std=c++11 #-I$(ROOTSYS)/include -pthread -std=c++17 -m64 -I/opt/sbg/cms/safe1/cms/bourgatte/Tau/Temp/include -I../../ -c
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################
USERINC    = 
USERLIBS   = 
SHAREDLIBFLAGS =


override CXXFLAGS += -I$(ROOTSYS)/include $(shell root-config --cflags) -I../../ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib $(shell root-config --glibs)

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

SF_SRCS    = $(wildcard *.cc)
SF_HDRS    = $(SF_SRCS:.cc=.h)

HDRS          =		interface/a1Helper.h \
		interface/rhoHelper.h \
		interface/TauPolInterface.h \
		interface/fonction_a1.h \
		interface/PolarimetricA1.h \
		interface/SCalculator.h


SRCS          = 	src/a1Helper.cc \
	        src/rhoHelper.cc \
		src/TauPolInterface.cc \
		src/fonction_a1.cc \
		src/PolarimetricA1.cc \
		src/SCalculator.cc \
		bin/example.cc 

OBJS1          = $(SRCS:.cc=.o) 
OBJS2         = $(SRCS:.cc=.o) 


PROGRAM       = run

SHAREDLIB     = libTDCILib.so
$(SHAREDLIB): $(OBJS2)
	@echo "======================================================="
	@echo "Print out environment:  "
	@echo "OBJS:   $(OBJS2)  "
	@echo "Linking SharedLib: $(SHAREDLIB) ..."
	@echo $(LD)  -fPIC $(LIBS)  -c  $(SRCS) -o $(OBJS2)
	@g++  -shared -o $(SHAREDLIB) $(OBJS2)
	@echo "Linking SharedLib: $(SHAREDLIB) Complete"
	@echo $(LD) $(LDFLAGS) src/*.o bin/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) src/*.o bin/*.o $(LIBS) -o $(PROGRAM)
	@mv libTDCILib.so bin
	@mv run bin
	@echo "======================================================="



$(PROGRAM): $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) src/*.o bin/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) src/*.o bin/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"


##########################################################
###clean  and PHONY rules                              ###
##########################################################

$(OBJS2): %.o: %.cc
	$(CXX) $(CXXFLAGS) $(DEFS) -fPIC $< -o $@




.PHONY: clean install 

install: $(SHAREDLIB)
all:  $(SHAREDLIB) $(PROGRAM)
clean:
	@rm src/*.o bin/*.o
	@rm libUserLib.so
	@rm run.exe
