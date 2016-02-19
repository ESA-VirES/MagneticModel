CC=gcc
RM=rm -fv
CFLAGS=-fPIC -g -I/usr/include/geomag
LFLAGS=-lgeomag -lm
INSTAL_EXE_DIR=/usr/bin
TARGET=wmm_file.exe wmm_grid.exe wmm_point.exe
HDR=
OBJ=wmm_file.o wmm_grid.o wmm_point.o

%.exe: %.o $(HDR)
	$(CC) -o $@ $< $(LFLAGS)

%.o: %.c $(HDR)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build 

build: $(TARGET)

install: $(TARGET) 
	mkdir -p $(INSTAL_EXE_DIR)
	-for F in $(TARGET) ;\
	do \
	   cp -v $$F $(INSTAL_EXE_DIR) && \
	   chmod 0755 $(INSTAL_EXE_DIR)/$$F ; \
	done

clean:
	$(RM) $(OBJ) $(TARGET)

