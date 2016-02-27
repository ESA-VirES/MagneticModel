CC=gcc
LD=gcc
RM=rm -fv

INC_WMM=
LIB_WMM=-lgeomag -lm
INC_QD=-I./include
LIB_QD=-lqdipole -lm
INC_PYTHON=$(shell python-config --includes )
LIB_PYTHON=$(shell python-config --libs )

INCS=$(INC_PYTHON) $(INC_WMM) $(INC_QD)
LIBS_WMM=$(LIB_PYTHON) $(LIB_WMM)
LIBS_QD=$(LIB_PYTHON) $(LIB_QD)

#CCFLAGS= -c -g -DDEBUG=1 -Wall -fPIC $(INCS)
CCFLAGS= -c -O2 -Wall -Wno-unused-function -fPIC $(INCS)
#CCFLAGS= -c -O2 -g -DDEBUG=1 -Wall -fPIC $(INCS)
#CCFLAGS= -c -O2 -Wall -fPIC $(INC)

HDR_WMM=pywmm_aux.h pywmm_sheval.h pywmm_sphargrd.h pywmm_spharpot.h\
    pywmm_legendre.h pywmm_lonsincos.h pywmm_relradpow.h\
    pywmm_vrot_common.h pywmm_vrot_sph2cart.h pywmm_vrot_sph2geod.h\
    pywmm_cconv.h pywmm_coord.h pywmm_exc.h
OBJ_WMM=pywmm.o
DST_WMM=_pywmm.so

HDR_QD=pydq_eval_apex.h pyqd_exc.h pywmm_aux.h
OBJ_QD=pyqd.o
DST_QD=_pyqd.so

# --------------------------------------------

all: pywmm pyqd

%.o: %.c
	$(CC) $< $(CCFLAGS) -o $@

pywmm: $(OBJ_WMM)
	$(LD) -shared -o $(DST_WMM) $(OBJ_WMM) $(LIBS_WMM)

pyqd: $(OBJ_QD)
	$(LD) -shared -o $(DST_QD) $(OBJ_QD) $(LIBS_QD)

clean:
	$(RM) *.o

purge: clean
	$(RM) $(DST_WMM)

pywmm.o: pywmm.c $(HDR_WMM)
pyqd.o: pyqd.c $(HDR_QD)

