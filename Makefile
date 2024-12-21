############################################################
# The Target Binary Program
TARGET     := Twinkle

SRCDIR     := src
USRDIR     ?= usr
INCDIR     :=
BUILDDIR   := obj
TARGETDIR  := bin
SRCEXT     := cpp
DEPEXT     := d
OBJEXT     := o

ARCH       ?= CUDA
DEBUG      ?= 0
MPI        ?= 0

# Flags, Libraries and Includes
CFLAGS_CMD :=
INCDEP     :=
ifeq ($(ARCH),CUDA)
	CC     := nvcc
	ifeq ($(DEBUG),1)
	CFLAGS := -std=c++17 -rdc=true -arch=sm_70 -O0 -g\
              -include ./src/utils/debug_macros.h \
	          -G -D __GPU_DEBUG__ -diag-suppress 3
	else   #DEBUG
	CFLAGS := -std=c++17 -rdc=true -arch=sm_70 -O3 \
              --use_fast_math $(CFLAGS_CMD)
    endif  #DEBUG
    INC    := -x cu
	LIB    := $(CFLAGS)
else ifeq ($(ARCH),HIPCPU)
	CC     := g++-10
	ifeq ($(DEBUG),1) 
	CFLAGS := -std=c++17 -O0 -g -fmax-errors=2 \
              -D __CPU_DEBUG__ \
              -include ./src/utils/debug_macros.h
	else   #DEBUG
	CFLAGS := -std=c++17 -O3 $(CFLAGS_CMD)
	endif  #DEBUG
	LIB    := $(CFLAGS) -lpthread -ltbb
	INC    := -include hip/hip_runtime.h
else ifeq ($(ARCH),HIP)
	CC     := hipcc
	ifeq ($(DEBUG),1) 
	CFLAGS := -std=c++17 -O0 -g -include hip/hip_runtime.h \
	        -fgpu-rdc -fPIC
	else   #DEBUG
	CFLAGS := -std=c++17 -O2 -include hip/hip_runtime.h \
	        -fgpu-rdc -fPIC -ffast-math \
			-ffp-contract=fast
	endif  #DEBUG
	LIB    := $(CFLAGS) --hip-link 
	INC    :=
else ifeq ($(ARCH),MUSA)
	CC     := mcc
	ifeq ($(DEBUG),1)
    CFLAGS := -std=c++17 -O2 -g $(CFLAGS_CMD)
	else   #DEBUG
	CFLAGS := -std=c++17 -Ofast -ffast-math $(CFLAGS_CMD) \
	         # --offload-arch=mp_21
    endif  #DEBUG
    INC    := -x musa
	LIB    := $(CFLAGS) -lmusa -lmusart -lpthread
endif

ifeq ($(MPI),1) 
CFLAGS     += -D__MPI__
LIB        += -lmpi
endif

ifeq ($(PRECISION),0)
CFLAGS     += -DPRECISION=0
else ifeq ($(PRECISION),2)
CFLAGS     += -DPRECISION=2
endif

############################################################
# Automatic generating objs and deps

SOURCES := $(shell find -L $(SRCDIR) -type f \
             -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,\
            $(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
USRSRC  := $(shell find -L $(USRDIR) -type f \
             -name *.$(SRCEXT))
OBJECTS += $(patsubst $(USRDIR)/%,$(BUILDDIR)/$(USRDIR)/%,\
            $(USRSRC:.$(SRCEXT)=.$(OBJEXT)))

.PHONY : all dirs clean remake print

#Defauilt Make
all: dirs $(TARGET)

# print : ; $(info $$(SOURCES) is [${SOURCES}])
print : ; $(info $$(OBJECTS) is [${OBJECTS}])

# Remake
remake: clean all

# Make the Directories
dirs:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(BUILDDIR)/$(USRDIR)

#C lean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)/* $(TARGETDIR)/*

# Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

# Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

# Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
	@$(CC) $(CFLAGS) $(INCDEP) -M $(SRCDIR)/$*.$(SRCEXT) >\
        $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) \
           $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' \
         < $(BUILDDIR)/$*.$(DEPEXT).tmp \
         > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' \
          < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | \
          sed -e 's/^ *//' -e 's/$$/:/' >> \
          $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

$(BUILDDIR)/$(USRDIR)/%.$(OBJEXT): $(USRDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
	@$(CC) $(CFLAGS) $(INCDEP) -M $(USRDIR)/$*.$(SRCEXT) >\
        $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT) \
           $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$(USRDIR)/$*.$(OBJEXT):|' \
         < $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT).tmp \
         > $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' \
        < $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT).tmp | fmt -1 |\
          sed -e 's/^ *//' -e 's/$$/:/' >> \
          $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$(USRDIR)/$*.$(DEPEXT).tmp
