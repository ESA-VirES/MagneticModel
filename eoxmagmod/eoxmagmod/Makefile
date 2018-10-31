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

HDR_WMM=pymm_aux.h pymm_sheval.h pymm_sphargrd.h pymm_spharpot.h\
    pymm_legendre.h pymm_lonsincos.h pymm_relradpow.h\
    pymm_vrot_common.h pymm_vrot_sph2cart.h pymm_vrot_sph2geod.h\
    pymm_vrot_cart2sph.h pymm_cconv.h pymm_coord.h pymm_exc.h
OBJ_WMM=pymm.o
DST_WMM=_pymm.so

HDR_QD=pymm_aux.h pyqd_exc.h\
    pyqd_eval_qdlatlon.h pyqd_eval_mlt.h pyqd_eval_subsol.h pyqd_eval_apex.h
OBJ_QD=pyqd.o
DST_QD=_pyqd.so

HDR_SP=pymm_aux.h pysunpos.h pysunpos_original.h
OBJ_SP=pysunpos.o
DST_SP=_pysunpos.so

HDR_TC=pymm_aux.h pytimeconv_mjd2000_to_decimal_year.h pytimeconv_decimal_year_to_mjd2000.h
OBJ_TC=pytimeconv.o
DST_TC=_pytimeconv.so

# --------------------------------------------

all: pymm pyqd pysunpos pytimeconv

%.o: %.c
	$(CC) $< $(CCFLAGS) -o $@

%.o: %.cpp
	$(CC) $< $(CCFLAGS) -o $@

pymm: $(OBJ_WMM)
	$(LD) -shared -o $(DST_WMM) $(OBJ_WMM) $(LIBS_WMM)

pyqd: $(OBJ_QD)
	$(LD) -shared -o $(DST_QD) $(OBJ_QD) $(LIBS_QD)

pysunpos: $(OBJ_SP)
	$(LD) -shared -o $(DST_SP) $(OBJ_SP) $(LIBS_SP)

pytimeconv: $(OBJ_TC)
	$(LD) -shared -o $(DST_TC) $(OBJ_TC) $(LIBS_TC)

clean:
	$(RM) *.o

purge: clean
	$(RM) $(DST_WMM)

pymm.o: pymm.c $(HDR_WMM)
pyqd.o: pyqd.c $(HDR_QD)
pysunpos.o: pysunpos.c $(HDR_SP)
pytimeconv.o: pytimeconv.c $(HDR_TC)
