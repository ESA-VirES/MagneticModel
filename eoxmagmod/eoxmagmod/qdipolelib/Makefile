F90=gfortran
CC=gcc
LD=gfortran
RM=rm -fv
CFLAGS=-fPIC -O2 -g
LFLAGS=-shared -Wl,-rpath=$(TARGET_LIB_DIR)
TARGET_LIB_DIR=/usr/lib64/qdipole
INSTAL_LIB_DIR=/usr/lib64/qdipole
INSTAL_INC_DIR=/usr/include/qdipole
INSTAL_DOC_DIR=/usr/share/doc/qdipole
LDD_CONF_FILE=/etc/ld.so.conf.d/qdipole.conf
TARGET=libqdipole.so
HDR=cqdipole.h
OBJ=cqdipole.o qdipole.o apex.o apexsh.o makeapexsh.o
DOC=README
MOD=alfbasismodule.mod apxshmodule.mod

%.o: %.c $(HDR)
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.f90 $(HDR)
	$(F90) -c -o $@ $< $(CFLAGS)

all: build

build: $(TARGET)

build_test:
	$(CC) -o qdipole_test qdipole_test.c -lqdipole -L$(INSTAL_LIB_DIR)

clean_test:
	rm -f qdipole_test 

install: $(TARGET) $(HDR) $(DOC)
	mkdir -p $(INSTAL_LIB_DIR)
	cp -v $(TARGET) $(INSTAL_LIB_DIR) 
	chmod 0755 $(INSTAL_LIB_DIR)/$(TARGET)
	mkdir -p $(INSTAL_INC_DIR)
	ldconfig 
	-for F in $(HDR) ;\
	do \
	   cp -v $$F $(INSTAL_INC_DIR) && \
	   chmod 0644 $(INSTAL_INC_DIR)/$$F ; \
	done
	mkdir -p $(INSTAL_DOC_DIR)
	-for F in $(DOC) ;\
	do \
	   cp -v $$F $(INSTAL_DOC_DIR) && \
	   chmod 0644 $(INSTAL_DOC_DIR)/$$F ; \
	done
	echo $(INSTAL_LIB_DIR) > $(LDD_CONF_FILE)
	ldconfig 

clean:
	$(RM) $(OBJ) $(TARGET) $(MOD)

uninstall:
	-for F in $(HDR) ;\
	do \
	   rm -f $(INSTAL_INC_DIR)/$$F ; \
	done
	rmdir $(INSTAL_INC_DIR)
	-for F in $(DOC) ;\
	do \
	   rm -f $(INSTAL_DOC_DIR)/$$F ; \
	done
	rmdir $(INSTAL_DOC_DIR)
	rm -f $(INSTAL_LIB_DIR)/$(TARGET)
	rmdir $(INSTAL_LIB_DIR)
	rm -f $(INSTAL_LIB_DIR)
	ldconfig 

$(TARGET): $(OBJ) 
	$(LD) -o $(TARGET) $(OBJ) $(LFLAGS)
