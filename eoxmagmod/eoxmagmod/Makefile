CC=gcc
LD=gcc
RM=rm -fv 

INC_WMM=-I/usr/include/geomag -I./include
LIB_WMM=-lgeomag -lm
INC_PYTHON=$(shell python-config --includes )
LIB_PYTHON=$(shell python-config --libs )

INCS=$(INC_PYTHON) $(INC_WMM)
LIBS=$(LIB_PYTHON) $(LIB_WMM) 

#CCFLAGS= -c -g -DDEBUG=1 -Wall -fPIC $(INCS)
CCFLAGS= -c -O2 -Wall -fPIC $(INCS)
#CCFLAGS= -c -O2 -g -DDEBUG=1 -Wall -fPIC $(INCS)
#CCFLAGS= -c -O2 -Wall -fPIC $(INC)

HDR=pywmm_aux.h pywmm_geomag.h pywmm_sphargrd.h\
    pywmm_legendre.h pywmm_lonsincos.h pywmm_relradpow.h\
    pywmm_vrot_common.h pywmm_vrot_sph2cart.h pywmm_vrot_sph2geod.h\
    pywmm_cconv.h pywmm_coord.h pywmm_exc.h
OBJ=pywmm.o
DST=_pywmm.so

# --------------------------------------------

all: pywmm

%.o: %.c
	$(CC) $< $(CCFLAGS) -o $@

pywmm: $(OBJ)
	$(LD) -shared -o $(DST) $(OBJ) $(LIBS) 

clean:
	$(RM) *.o

purge: clean 
	$(RM) $(DST)

pywmm.o: pywmm.c $(HDR)

