INTELC_VERSION = `icpc -dumpversion`
MACHINE = `uname -m`
VERSION = "1.1"
#SCIFT_HOME = /media/14eacb85-2460-4c08-90dc-dacc9adabbef/nestor/Dropbox/UAM/scift
SCIFT_HOME = /home/nestor/Dropbox/UAM/scift

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
	mkdir M3C
	mkdir M3C/doc
	cp doc/tutorial/M3C-tutorial.pdf M3C/doc/
#	cp -r examples M3C
	cp -r utils M3C
	mkdir M3C/bin
	find src -type f -executable -exec cp {} M3C/bin \;
	cp ${SCIFT_HOME}/examples/molecule.compare M3C/bin
	cp ${SCIFT_HOME}/examples/molecule.fv M3C/bin
	cp ${SCIFT_HOME}/examples/molecule.inertia M3C/bin
	cp ${SCIFT_HOME}/examples/molecule.mass M3C/bin
	cp ${SCIFT_HOME}/examples/molecule.random M3C/bin
	cp ${SCIFT_HOME}/examples/molecule.radius M3C/bin
	sleep 10s
	tar cvfz M3C-${VERSION}-intelc-${INTELC_VERSION}-${MACHINE}.tar.gz M3C
	rm -rf M3C
        
