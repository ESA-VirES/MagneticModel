#
# install CDF library from sources
#

VERSION = 39_2
DIRNAME = ./cdf$(VERSION)-dist-main/
FILENAME = ./cdf$(VERSION)-dist-cdf.tar.gz
SOURCE_URL = https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf$(VERSION)/linux/cdf$(VERSION)-dist-cdf.tar.gz
MD5CHECKSUM = b0708e2120933c6eac83bab5025154f6
INSTALLDIR ?= /usr/local/cdf

all: build

clean:
	-rm -fR $(DIRNAME)

install: build test
	touch $(DIRNAME)/README_cdf_tools.txt # 3.9.2 error - missing file
	make -C $(DIRNAME) INSTALLDIR=$(INSTALLDIR) install -f Makefile-classic

test: build
	make -C $(DIRNAME) test -f Makefile-classic

build: $(DIRNAME)
	make -C $(DIRNAME) OS=linux ENV=gnu CURSES=no all -f Makefile-classic

$(DIRNAME): $(FILENAME) md5check
	tar -xzf $(FILENAME)

md5check: $(FILENAME)
	echo $(MD5CHECKSUM) $(FILENAME) | md5sum -c

$(FILENAME):
	curl $(SOURCE_URL) > $(FILENAME)
