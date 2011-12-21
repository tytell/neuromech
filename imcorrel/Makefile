MATLAB = /usr/local/bin/matlab
MATLABROOT = /Applications/MATLAB_R2011b.app
CC = gcc
CXX = g++
LD = gcc

CPPFLAGS1 = -fno-common -no-cpp-precomp -arch x86_64 -isysroot \
	/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.5  -fexceptions \
	-DMATLAB_MEX_FILE
CINCLUDEDIR = -I$(MATLABROOT)/extern/include -I$(MATLABROOT)/simulink/include 
#use for 32bit back compatibility
#CFLAGS32 = -DMX_COMPAT_32
CDEBUGFLAGS = -g
COPTIMFLAGS = -O2 -DNDEBUG

CPPFLAGS = $(CPPFLAGS1) $(CINCLUDEDIR) $(CDEBUGFLAGS) 

CPPFILES = imcorrel.cpp image.cpp error.cpp
OFILES = imcorrel.o image.o error.o

LDFLAGS1 = -Wl,-twolevel_namespace -undefined error -arch x86_64 \
	-Wl,-syslibroot,/Developer/SDKs/MacOSX10.6.sdk -mmacosx-version-min=10.5 \
	-bundle \
	-Wl,-exported_symbols_list\
	-Wl,/Applications/MATLAB_R2011b.app/extern/lib/maci64/mexFunction.map
LDDEBUGFLAGS = -g
LDOPTIMFLAGS = -O
LDEXTENSION = .mexmaci64
LIBS = -L/Applications/MATLAB_R2011b.app/bin/maci64 -lmx -lmex -lmat -lstdc++

LDFLAGS = $(LDFLAGS1) $(LDDEBUGFLAGS)

all: imcorrel$(LDEXTENSION)

imcorrel$(LDEXTENSION): $(OFILES)
	$(LD) $(LDFLAGS) $(LDDEBUGFLAGS)  $(OFILES) $(LIBS) -o $@

mexcpp$(LDEXTENSION): mexcpp.o
	$(LD) $(LDFLAGS) $(LDDEBUGFLAGS)  mexcpp.o $(LIBS) -o $@
