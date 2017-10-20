CC=gcc
CPP=g++
LD=gcc
RM=rm -fv

INC_WMM=
LIB_WMM=-lgeomag -lm
INC_QD=-I./include
LIB_QD=-lqdipole -lm
INC_SP=-I./include
LIB_SP=-lm
INC_PYTHON=$(shell python-config --includes )
LIB_PYTHON=$(shell python-config --libs )

INCS=$(INC_PYTHON) $(INC_WMM) $(INC_QD)
LIBS_WMM=$(LIB_PYTHON) $(LIB_WMM)
LIBS_QD=$(LIB_PYTHON) $(LIB_QD)
LIBS_QD=$(LIB_PYTHON) $(LIB_QD)

#CCFLAGS= -c -g -DDEBUG=1 -Wall -fPIC $(INCS)
CCFLAGS= -c -O2 -Wall -Wno-unused-function -fPIC $(INCS)
#CCFLAGS= -c -O2 -g -DDEBUG=1 -Wall -fPIC $(INCS)
#CCFLAGS= -c -O2 -Wall -fPIC $(INC)

HDR_WMM=pywmm_aux.h pywmm_sheval.h pywmm_sphargrd.h pywmm_spharpot.h\
    pywmm_legendre.h pywmm_lonsincos.h pywmm_relradpow.h\
    pywmm_vrot_common.h pywmm_vrot_sph2cart.h pywmm_vrot_sph2geod.h\
    pywmm_vrot_cart2sph.h pywmm_cconv.h pywmm_coord.h pywmm_exc.h
OBJ_WMM=pywmm.o
DST_WMM=_pywmm.so

HDR_QD=pywmm_aux.h pyqd_exc.h\
    pyqd_eval_qdlatlon.h pyqd_eval_mlt.h pyqd_eval_subsol.h pyqd_eval_apex.h
OBJ_QD=pyqd.o
DST_QD=_pyqd.so

HDR_SP=pywmm_aux.h pysunpos.h pysunpos_original.h
OBJ_SP=pysunpos.o
DST_SP=_pysunpos.so

# --------------------------------------------

all: pywmm pyqd pysunpos

%.o: %.c
	$(CC) $< $(CCFLAGS) -o $@

%.o: %.cpp
	$(CC) $< $(CCFLAGS) -o $@

pywmm: $(OBJ_WMM)
	$(LD) -shared -o $(DST_WMM) $(OBJ_WMM) $(LIBS_WMM)

pyqd: $(OBJ_QD)
	$(LD) -shared -o $(DST_QD) $(OBJ_QD) $(LIBS_QD)

pysunpos: $(OBJ_SP)
	$(LD) -shared -o $(DST_SP) $(OBJ_SP) $(LIBS_SP)


clean:
	$(RM) *.o

purge: clean
	$(RM) $(DST_WMM)

pywmm.o: pywmm.c $(HDR_WMM)
pyqd.o: pyqd.c $(HDR_QD)
pysunpos.o: pysunpos.c $(HDR_SP)
