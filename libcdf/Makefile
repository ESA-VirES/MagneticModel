#
# install CDF library from sources
#

VERSION = 39_2
DIRNAME = ./cdf$(VERSION)-dist-main/
FILENAME = ./cdf$(VERSION)-dist-cdf.tar.gz
SOURCE_URL = https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf$(VERSION)/linux/cdf$(VERSION)-dist-cdf.tar.gz
MD5CHECKSUM = 9d1abe5559640ede524426c2213af303
INSTALLDIR ?= /usr/local/cdf

all: build

clean:
	-rm -fR $(DIRNAME)

install: build test
	make -C $(DIRNAME) OS=linux ENV=gnu INSTALLDIR=$(INSTALLDIR) install -f Makefile-classic

test: build
	make -C $(DIRNAME) test -f Makefile-classic

build: $(DIRNAME)
	make -C $(DIRNAME) OS=linux ENV=gnu CURSES=no all -f Makefile-classic

$(DIRNAME): $(FILENAME) md5check
	tar -xzf $(FILENAME)

md5check: $(FILENAME)
	echo $(MD5CHECKSUM) $(FILENAME) | md5sum -c

download: $(FILENAME) md5check

$(FILENAME):
	curl $(SOURCE_URL) > $(FILENAME)
