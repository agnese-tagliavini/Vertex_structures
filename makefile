# MAKEFILE for fRG_O2 Project

#--------------------------------------General------------------------------------------
CXX 	:= g++ # This is the main compiler
#CXX    := icpc # This is the main compiler
SRCDIR  := src
HEADDIR := include
BUILDDIR := build
TARGET  := bin/run


#--------------------------------------Sources and header files------------------------------------------
SRCEXT  := cpp
HEADEXT := h
#SOURCES := $(shell find $(SRCDIR) -name '*.$(SRCEXT)') 
#HEADERS := $(shell find $(HEADDIR) -name '*.$(HEADEXT)') 
SOURCES := $(shell ls -1 $(SRCDIR)/*.$(SRCEXT)) 
HEADERS := $(shell ls -1 $(HEADDIR)/*.$(HEADEXT)) 
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

#--------------------------------------Settings for fRG flow------------------------------------------

CFLAGS += -D NO_MOMENTA # 					Set flag for calculating impurity models only (no k dependence)
#CFLAGS += -D FORCED_ZEROS # 					Use forced zeros for the vertex
CFLAGS += -D COMPRESS #						Compress HDF5 datafile 

# -- Specify hybridization function in case of impurity setup (default = Wide Band limit)
#CFLAGS += -D QMC_BATH
#CFLAGS += -D ED_BATH

#----Choose if input file to read
CFLAGS += -D READIN

# -- Choose cuttoff scheme
CFLAGS += -D INT_FLOW #			Interaction Flow
#CFLAGS += -D EBERL_FLOW #		Eberlein Flow
#CFLAGS += -D OMEGA_FLOW #		Omega flow according to Salmhofer/Gierig ( arxiv:1208.6131 )
#CFLAGS += -D RES_FLOW #		Reservoir cutoff scheme

#------Choose which method for the inversion of the Bethe-Salpeter Equations
CFLAGS += -D METHOD2 #Agnese's method
#CFLAGS += -D METHOD1 #Stefan's method
#CFLAGS += -D METHOD3 # Brute force inversion

# -- Higher order corrections
#CFLAGS += -D KATANIN #			Use katanin scheme
#CFLAGS += -D TWOLOOP #			Use twoloop corrections

#--------------------------------------Compiler settings------------------------------------------

CFLAGS += -std=c++11 -D OPENMP -fopenmp #			General compiler flags
CFLAGS += -D BOOST_DISABLE_ASSERTS #				Disable boost assert checks
DBFLAGS := -O2 -g #						Compiler flags for debugging
PROFFLAGS := -O2 -g #						Compiler flags for profiling
RUNFLAGS := -O3 #						Compiler flags for quick compile and run
OPTIMIZEFLAGS := -flto -march=native -O3 #			GCC Compiler flags for optimal speed
#OPTIMIZEFLAGS := -O3 -fp-model fast=2 -xHost # -no-prec-div #	Intel Compiler flags for optimal speed
H5LIB := $(shell h5c++ -show|cut -d ' ' -f2-) #			HDF5 libraries
LIB := -lgsl -lgslcblas -lm $(H5LIB) -fopenmp -lmpi #		Library flags
INC := -I include #						Additional include paths
INC += -I /usr/include/eigen3 		#Additional include paths
INC += -I /usr/include/hdf5/serial/ #Additional include paths

#--------------------------------------Targets ------------------------------------------
$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@mkdir -p bin
	@mkdir -p dat
	@echo " $(CXX) $^ -o $(TARGET) $(LIB)"; $(CXX) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

run: 	CFLAGS += $(RUNFLAGS)
run: 	$(TARGET)

debug: 	CFLAGS += $(DBFLAGS)
debug: 	$(TARGET)

prof: 	CFLAGS += $(PROFFLAGS)
prof: 	$(TARGET)

optimize: CFLAGS += $(OPTIMIZEFLAGS)
optimize: $(SOURCES) $(HEADERS)
	@mkdir -p bin
	@mkdir -p dat
	$(CXX) $(CFLAGS) $(INC) $(SOURCES) -o $(TARGET) $(LIB)
	#rm *.o

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean

