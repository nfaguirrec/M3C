MAKEFLAGS = -s
FC = ifort
GFCFLAGS = -w -I. -I/home/aguirre/Software/miniconda3/envs/python3.10/include -I/usr/include -cpp -O3 -g -fPIC -fdec-structure -ffree-line-length-none -I$(SCIFT_HOME)/src
GLDFLAGS = -L$(SCIFT_HOME)/src -L/home/aguirre/Software/miniconda3/envs/python3.10/lib -L/usr/lib -Wl,-rpath,/home/aguirre/Software/miniconda3/envs/python3.10/lib -lscift -lfftw3 -lopenblas -llapack
IFCFLAGS = -w -I. -O3 -g -fPIC -fpp -mkl -I$(SCIFT_HOME)/src
ILDFLAGS = -L$(SCIFT_HOME)/src -lscift
TLIB = libM3C.a

ifeq ($(findstring gfortran,$(FC)),gfortran)
    FCFLAGS = $(GFCFLAGS)
    LDFLAGS = $(GLDFLAGS)
else
    FCFLAGS = $(IFCFLAGS)
    LDFLAGS = $(ILDFLAGS)
endif

all:  AzizSlamanPotential.o Fragment.o FragmentsDB.o FragmentsListBase.o FragmentsList.o GOptionsM3C.o MarkovChain.o ModelPotential.o NNLS.o Reactor.o  $(TLIB) M3CBR M3C M3CfitBR molecule.ZPE 

AzizSlamanPotential.o azizslamanpotential_.mod: AzizSlamanPotential.f90 
	fmake
	echo -n 'Building AzizSlamanPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c AzizSlamanPotential.f90

Fragment.o fragment_.mod: Fragment.f90 goptionsm3c_.mod
	fmake
	echo -n 'Building Fragment.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c Fragment.f90

FragmentsDB.o fragmentsdb_.mod: FragmentsDB.f90 fragment_.mod azizslamanpotential_.mod goptionsm3c_.mod modelpotential_.mod
	fmake
	echo -n 'Building FragmentsDB.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c FragmentsDB.f90

FragmentsListBase.o fragmentslistbase_.mod: FragmentsListBase.f90 fragment_.mod azizslamanpotential_.mod goptionsm3c_.mod modelpotential_.mod fragmentsdb_.mod
	fmake
	echo -n 'Building FragmentsListBase.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c FragmentsListBase.f90

FragmentsList.o fragmentslist_.mod: FragmentsList.f90 fragment_.mod azizslamanpotential_.mod goptionsm3c_.mod modelpotential_.mod fragmentsdb_.mod fragmentslistbase_.mod
	fmake
	echo -n 'Building FragmentsList.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c FragmentsList.f90

GOptionsM3C.o goptionsm3c_.mod: GOptionsM3C.f90 
	fmake
	echo -n 'Building GOptionsM3C.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c GOptionsM3C.f90

M3CBR: M3CBR.f90 NNLS.o nnls_.mod
	fmake
	echo -n 'Building M3CBR.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) M3CBR.f90 NNLS.o -o M3CBR $(LDFLAGS)

M3C: M3C.f90 Fragment.o AzizSlamanPotential.o GOptionsM3C.o ModelPotential.o FragmentsDB.o FragmentsListBase.o FragmentsList.o MarkovChain.o Reactor.o fragment_.mod azizslamanpotential_.mod goptionsm3c_.mod modelpotential_.mod fragmentsdb_.mod fragmentslistbase_.mod fragmentslist_.mod markovchain_.mod reactor_.mod
	fmake
	echo -n 'Building M3C.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) M3C.f90 Fragment.o AzizSlamanPotential.o GOptionsM3C.o ModelPotential.o FragmentsDB.o FragmentsListBase.o FragmentsList.o MarkovChain.o Reactor.o -o M3C $(LDFLAGS)

M3CfitBR: M3CfitBR.f90 NNLS.o nnls_.mod
	fmake
	echo -n 'Building M3CfitBR.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) M3CfitBR.f90 NNLS.o -o M3CfitBR $(LDFLAGS)

MarkovChain.o markovchain_.mod: MarkovChain.f90 fragment_.mod azizslamanpotential_.mod goptionsm3c_.mod modelpotential_.mod fragmentsdb_.mod fragmentslistbase_.mod fragmentslist_.mod reactor_.mod
	fmake
	echo -n 'Building MarkovChain.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c MarkovChain.f90

ModelPotential.o modelpotential_.mod: ModelPotential.f90 goptionsm3c_.mod azizslamanpotential_.mod
	fmake
	echo -n 'Building ModelPotential.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c ModelPotential.f90

molecule.ZPE: molecule.ZPE.f90 Fragment.o GOptionsM3C.o fragment_.mod goptionsm3c_.mod
	fmake
	echo -n 'Building molecule.ZPE.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) molecule.ZPE.f90 Fragment.o GOptionsM3C.o -o molecule.ZPE $(LDFLAGS)

NNLS.o nnls_.mod: NNLS.f90 
	fmake
	echo -n 'Building NNLS.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c NNLS.f90

Reactor.o reactor_.mod: Reactor.f90 fragment_.mod azizslamanpotential_.mod goptionsm3c_.mod modelpotential_.mod fragmentsdb_.mod fragmentslistbase_.mod fragmentslist_.mod
	fmake
	echo -n 'Building Reactor.f90 '
	time -f "(%E)" $(FC) $(FCFLAGS) -c Reactor.f90

$(TLIB): AzizSlamanPotential.o Fragment.o FragmentsDB.o FragmentsListBase.o FragmentsList.o GOptionsM3C.o MarkovChain.o ModelPotential.o NNLS.o Reactor.o 
	echo 'Building library $(TLIB)'
	ar rvs $(TLIB) AzizSlamanPotential.o Fragment.o FragmentsDB.o FragmentsListBase.o FragmentsList.o GOptionsM3C.o MarkovChain.o ModelPotential.o NNLS.o Reactor.o  > /dev/null 2> /dev/null

clean:
	rm -f *.o *.mod *~ .deps M3CBR M3C M3CfitBR molecule.ZPE  2> /dev/null

distclean:
	rm -f *.o *.mod *~ .deps Makefile M3CBR M3C M3CfitBR molecule.ZPE  2> /dev/null

Makefile: .deps
	fmake

.deps:
	fmake
