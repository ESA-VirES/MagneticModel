CC?=gcc
RM=rm -fv

INC_MM=
LIB_MM=-lm
INC_QD=-I./include
LIB_QD=-lqdipole -lm
INC_SP=-I./include
LIB_SP=-lm
PYTHON_CONFIG ?=python3-config
INC_NUMPY=-I$(shell python -c 'import numpy; print(numpy.get_include())')
INC_PYTHON=$(shell $(PYTHON_CONFIG) --includes )
LIB_PYTHON=$(shell $(PYTHON_CONFIG) --ldflags )

INCS=$(INC_NUMPY) $(INC_PYTHON) $(INC_MM) $(INC_QD) $(INCLUDE)
LIBS_MM=$(LIB_PYTHON) $(LIB_MM)
LIBS_QD=$(LIB_PYTHON) $(LIB_QD)
LIBS_QD=$(LIB_PYTHON) $(LIB_QD)

#CCFLAGS= -c -g -DDEBUG=1 -Wall -fPIC $(INCS)
CCFLAGS= -c -O2 -Wall -Wno-unused-function -fPIC $(INCS)
#CCFLAGS= -c -O2 -g -DDEBUG=1 -Wall -fPIC $(INCS)
#CCFLAGS= -c -O2 -Wall -fPIC $(INC)

HDR_MM=pymm_aux.h pymm_shevaltemp.h pymm_sheval2dfs.h pymm_sheval.h \
    pymm_sphar_common.h pymm_sphargrd.h pymm_spharpot.h \
    pymm_legendre.h pymm_loncossin.h pymm_relradpow.h\
    pymm_vrot_common.h pymm_vrot_sph2cart.h pymm_vrot_sph2geod.h\
    pymm_vrot_cart2sph.h pymm_cconv.h pymm_coord.h pymm_bisect.h\
    pymm_interp.h pymm_fourier2d.h
OBJ_MM=pymm.o
DST_MM=_pymm.so

HDR_QD=pymm_aux.h pyqd_eval_qdlatlon.h pyqd_eval_mlt.h pyqd_eval_subsol.h
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

pymm: $(OBJ_MM)
	$(CC) -shared -o $(DST_MM) $(OBJ_MM) $(LIBS_MM)

pyqd: $(OBJ_QD)
	$(CC) -shared -o $(DST_QD) $(OBJ_QD) $(LIBS_QD)

pysunpos: $(OBJ_SP)
	$(CC) -shared -o $(DST_SP) $(OBJ_SP) $(LIBS_SP)

pytimeconv: $(OBJ_TC)
	$(CC) -shared -o $(DST_TC) $(OBJ_TC) $(LIBS_TC)

clean:
	$(RM) *.o

purge: clean
	$(RM) $(DST_MM)

pymm.o: pymm.c $(HDR_MM)
pyqd.o: pyqd.c $(HDR_QD)
pysunpos.o: pysunpos.c $(HDR_SP)
pytimeconv.o: pytimeconv.c $(HDR_TC)
