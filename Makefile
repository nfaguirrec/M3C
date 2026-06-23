MAKEFLAGS = -s
FC = ifort
SCIFT_HOME = /home/aguirre/Develop/scift
INTELC_VERSION = `icpc -dumpversion`
MACHINE = `uname -m`
VERSION = `cat VERSION`

all: build_src build_utils

build_src:
	make -C src FC="$(FC)" SCIFT_HOME="$(SCIFT_HOME)"

build_utils:
	make -C utils/libmsym/

clean:
	make -C src clean
	make -C utils/libmsym/ clean

binary:
	mkdir M3C-v$(VERSION)
	mkdir M3C-v$(VERSION)/doc
	cp README.md M3C-v$(VERSION)/
	cp M3Cvars.sh M3C-v$(VERSION)/
	cp doc/tutorial/tutorial-*.pdf M3C-v$(VERSION)/doc/
#	cp -r examples M3C-v$(VERSION)
	mkdir M3C-v$(VERSION)/utils
	find utils/ -maxdepth 1 -type f -exec cp {} M3C-v$(VERSION)/utils \;
	mkdir M3C-v$(VERSION)/bin
	find src -type f -executable -exec cp {} M3C-v$(VERSION)/bin \;
	find $(SCIFT_HOME)/examples/ -name "molecule.*" -type f -executable -exec cp {} M3C-v$(VERSION)/bin \;
	sleep 10s
	tar cvfz M3C-v$(VERSION)-intelc-$(INTELC_VERSION)-$(MACHINE).tar.gz M3C-v$(VERSION)
	rm -rf M3C-v$(VERSION)
