CC=gcc
RM=rm -fv
CFLAGS=-fPIC -O2 -g
LFLAGS=-shared -Wl,-rpath=$(TARGET_LIB_DIR)
TARGET_LIB_DIR=/usr/lib64
INSTAL_LIB_DIR=/usr/lib64
INSTAL_INC_DIR=/usr/include/geomag
INSTAL_DOC_DIR=/usr/share/doc/geomag
TARGET=libgeomag.so
HDR=GeomagnetismHeader.h EGM9615.h
OBJ=GeomagnetismLibrary.o
DOC=Changelog.txt WMM_README.txt README

%.o: %.c $(HDR)
	$(CC) -c -o $@ $< $(CFLAGS)

all: build install

build: $(TARGET)

install: $(TARGET) $(HDR) $(DOC)
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

clean:
	$(RM) $(OBJ) $(TARGET)

$(TARGET): $(OBJ) 
	$(CC) -o $(TARGET) $(OBJ) $(LFLAGS)
