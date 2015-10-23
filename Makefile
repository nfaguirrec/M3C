INTELC_VERSION = `icpc -dumpversion`
MACHINE = `uname -m`
VERSION = "1.2"
SCIFT_HOME = /home/nestor/Develop/scift

# all:
# 	cd src; make; cd ..
# 	cd examples; make; cd ..

# clean:
# 	cd src; make clean
	
# install: src/libscift.a
# 	echo -n 'Installing scift in ${HOME}/Software/scift ... '
# 	mkdir -p ${HOME}/Software/scift/finclude/
# 	cp src/*.mod ${HOME}/Software/scift/finclude/
# 	cp src/lib*.a ${HOME}/Software/scift/
# 	echo 'OK'
# 
# uninstall: ${HOME}/Software/scift/libscift.a
# 	echo -n 'Uninstalling scift from ${HOME}/Software/scift ... '
# 	rm -rf ${HOME}/Software/scift
# 	echo 'OK'

# doc: doxyfile src/libscift.a
# 	doxygen doxyfile

binary:
	mkdir M3C-v${VERSION}
	mkdir M3C-v${VERSION}/doc
	cp doc/tutorial/M3C-tutorial.pdf M3C-v${VERSION}/doc/
#	cp -r examples M3C-v${VERSION}
	cp -r utils M3C-v${VERSION}
	mkdir M3C-v${VERSION}/bin
	find src -type f -executable -exec cp {} M3C-v${VERSION}/bin \;
	cp ${SCIFT_HOME}/examples/molecule.compare M3C-v${VERSION}/bin
	cp ${SCIFT_HOME}/examples/molecule.fv M3C-v${VERSION}/bin
	cp ${SCIFT_HOME}/examples/molecule.inertia M3C-v${VERSION}/bin
	cp ${SCIFT_HOME}/examples/molecule.mass M3C-v${VERSION}/bin
	cp ${SCIFT_HOME}/examples/molecule.random M3C-v${VERSION}/bin
	cp ${SCIFT_HOME}/examples/molecule.radius M3C-v${VERSION}/bin
	cp ${SCIFT_HOME}/examples/molecule.minMult M3C-v${VERSION}/bin
	cp ${SCIFT_HOME}/examples/molecule.rotate M3C-v${VERSION}/bin
	sleep 10s
	tar cvfz M3C-v${VERSION}-intelc-${INTELC_VERSION}-${MACHINE}.tar.gz M3C-v${VERSION}
	rm -rf M3C-v${VERSION}
        
