# Very basic makefile
SILENT      ?= @

CXXFLAGS = -fPIC -g  -O2 -W -Wall -Wextra -Wshadow 
SOMAKER = g++
SOFLAGS = -shared 
SOSUFFIX = so 

ROOT_INCLUDE_FLAGS := $(shell root-config --cflags) 
ROOT_LIB_FLAGS := $(shell root-config --libs) -lRooFit -lRooFitCore -lMinuit -lFoam
ROOTCINT := $(shell root-config --bindir)/rootcint
ROOT := $(shell root-config --bindir)/root

# --- Directories -------------------------------------------------------------------
SRCDIR := ../src
BINDIR := ../bin
BUILDDIR := ../build
INCDIR := ../include
LIBDIR := ../lib
LINKDEFIN := $(wildcard $(INCDIR)/*LinkDef.h.in)

LIBNAME := nEXOSensitivity

# Make dependencies
MAKEDEPEND = rm -f $@; \
	$(CXX) $(CCFLAGS) -MM $(INCLUDEFLAGS) $< > $@.$$$$;                           \
        cp $@.$$$$ $@.temp.$$$$; sed -e 's/\#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
        -e '/^$$/ d' -e 's/$$/ :/' < $@.$$$$ >> $@.temp.$$$$;                         \
	sed 's,\($*\)\.o[ :]*,$(BUILDDIR)/\1.o $@ : ,g' < $@.temp.$$$$ > $@ ;         \
        rm -f $@*.$$$$;


EXECUTABLES = 
SHLIB := $(LIBDIR)/lib$(LIBNAME).$(SOSUFFIX)
AUTOLINKDEF := $(BUILDDIR)/LinkDef_$(LIBNAME).h
PYTHON := python
ROOTSOSUFFIX := so
ROOTDICTCCNAME := $(BUILDDIR)/$(LIBNAME)DICT.C
ROOTDICTHHNAME := $(BUILDDIR)/$(LIBNAME)DICT.h
ROOTDICTOBJECT := $(ROOTDICTCCNAME:.C=.o)
INCLUDEFLAGS = $(ROOT_INCLUDE_FLAGS) -I$(INCDIR)
LIBFLAGS = $(ROOT_LIB_FLAGS) 
DICTINCLUDEFLAGS= $(filter-out -stdlib=libc++ -std=c++11, $(INCLUDEFLAGS))

SOURCES := $(wildcard $(SRCDIR)/*.cxx)
OBJECTS := $(patsubst $(SRCDIR)/%.cxx, $(BUILDDIR)/%.o, $(SOURCES))
DEPCCFILES := $(patsubst $(SRCDIR)/%.cxx, $(BUILDDIR)/%.d, $(SOURCES))
#OBJECTS := $(BUILDDIR)/$(SOURCES:.cxx=.o)
ROOTDICTHEADERS := $(wildcard $(INCDIR)/*.hh)
BUILDDIRCONF := $(BUILDDIR)/.buildconf

.PHONY: all depend lib clean static shared

all: exe lib 

exe: $(EXECUTABLES) lib

lib: shared 

shared: $(SHLIB)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPCCFILES) 
endif

$(BUILDDIRCONF):
	$(SILENT)if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
	$(SILENT)touch $@

$(DEPCCFILES): $(BUILDDIR)/%.d: $(SRCDIR)/%.cxx $(BUILDDIRCONF)
	$(SILENT)echo "Update dependencies........ $(<F)"
	$(SILENT)$(MAKEDEPEND)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cxx 
	$(SILENT)echo "Compiling file............. $(<F)"
	$(SILENT)$(CXX) -c $(CXXFLAGS) $(INCLUDEFLAGS) $< -o $@

%.o: %.C
	$(SILENT)$(CXX) -c $(CXXFLAGS) $(INCLUDEFLAGS) $< -o $@

$(EXECUTABLES): $(OBJECTS) $(ROOTDICTOBJECT) 
	$(SILENT)$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

$(ROOTDICTCCNAME): $(ROOTDICTHEADERS) $(AUTOLINKDEF) 
	$(SILENT)echo "Rootifying files........... $(@F)"
	$(SILENT)rm -f $(ROOTDICTCCNAME) $(ROOTDICTHHNAME) 
	$(SILENT)$(ROOTCINT) $@ -c -p $(CXXFLAGS) $(ROOTDICTINCLUDE) $(DICTINCLUDEFLAGS) $^ 
 
$(SHLIB): $(OBJECTS) $(ROOTDICTOBJECT)
	$(SILENT)echo "Building library........... $(@F)"
	$(SILENT)if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR); fi
	$(SILENT)$(SOMAKER) $(SOFLAGS) -o $@ $^ $(LIBFLAGS)

$(AUTOLINKDEF): $(LINKDEFIN) $(ROOTDICTHEADERS) 
	$(SILENT)echo "Auto-making Linkdef........ $(@F)"
	$(SILENT)$(PYTHON) makeLinkDef.py\
        --output=$@  \
        --basebuild=$(INCDIR)    \
        --input_linkdef=$(LINKDEFIN) $(ROOTDICTHEADERS)

clean:
	$(SILENT)rm -f $(SHLIB) $(ROOTDICTCCNAME) $(ROOTDICTHHNAME) $(BUILDDIR)/*.o *~ $(BUILDDIR)/*.d
	$(SILENT)rm -f $(AUTOLINKDEF) 
	$(SILENT)rm -f $(EXECUTABLES)


