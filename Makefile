INTELC_VERSION = `icpc -dumpversion`
MACHINE = `uname -m`
VERSION = `cat VERSION`
SCIFT_HOME = /home/nestor/Develop/scift

all:
	cd src; make; cd ..
# 	cd examples; make; cd ..

clean:
	cd src; make clean
	
# doc: doxyfile src/libscift.a
# 	doxygen doxyfile

binary:
	mkdir M3C-v${VERSION}
	mkdir M3C-v${VERSION}/doc
	cp README M3C-v${VERSION}/
	cp NEWS M3C-v${VERSION}/
	cp M3Cvars.sh M3C-v${VERSION}/
	cp doc/tutorial/M3C-tutorial.pdf M3C-v${VERSION}/doc/
#	cp -r examples M3C-v${VERSION}
	cp -r utils M3C-v${VERSION}
	mkdir M3C-v${VERSION}/bin
	find src -type f -executable -exec cp {} M3C-v${VERSION}/bin \;
	find ${SCIFT_HOME}/examples/ -name "molecule.*" -type f -executable -exec cp {} M3C-v${VERSION}/bin \;
	sleep 10s
	tar cvfz M3C-v${VERSION}-intelc-${INTELC_VERSION}-${MACHINE}.tar.gz M3C-v${VERSION}
	rm -rf M3C-v${VERSION}
        
