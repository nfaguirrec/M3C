!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                   !!
!! This file is part of M3C project                                                  !!
!!                                                                                   !!
!! Copyright (c) 2019-2020 by authors                                                !!
!! Authors:                                                                          !!
!!                         * Néstor F. Aguirre (2020-2020)                           !!
!!                           nfaguirrec@gmail.com                                    !!
!!                                                                                   !!
!! Copyright (c) 2013-2016 Departamento de Química                                   !!
!!                         Universidad Autónoma de Madrid                            !!
!!                         All rights reserved.                                      !!
!!                                                                                   !!
!!                         * Néstor F. Aguirre (2013-2016)                           !!
!!                           nestor.aguirre@uam.es                                   !!
!!                         * Sergio Díaz-Tendero (2013-2016)                         !!
!!                           sergio.diaztendero@uam.es                               !!
!!                         * M. Paul-Antoine Hervieux (2013-2015)                    !!
!!                           Paul-Antoine.Hervieux@ipcms.unistra.fr                  !!
!!                         * Manuel Alcamí (2013-2016)                               !!
!!                           manuel.alcami@uam.es                                    !!
!!                         * Fernando Martín (2013-2016)                             !!
!!                           fernando.martin@uam.es                                  !!
!!                                                                                   !!
!!  Redistribution and use in source and binary forms, with or without               !!
!!  modification, are permitted provided that the following conditions are met:      !!
!!                                                                                   !!
!!  1. Redistributions of source code must retain the above copyright notice, this   !!
!!     list of conditions and the following disclaimer.                              !!
!!  2. Redistributions in binary form must reproduce the above copyright notice,     !!
!!     this list of conditions and the following disclaimer in the documentation     !!
!!     and/or other materials provided with the distribution.                        !!
!!  3. Neither the name of the copyright holders nor the names of its contributors   !!
!!     may be used to endorse or promote products derived from this software         !!
!!     without specific prior written permission.                                    !!
!!                                                                                   !!
!!  The copyright holders provide no reassurances that the source code provided      !!
!!  does not infringe any patent, copyright, or any other intellectual property      !!
!!  rights of third parties.  The copyright holders disclaim any liability to any    !!
!!  recipient for claims brought against recipient by any third party for            !!
!!  infringement of that parties intellectual property rights.                       !!
!!                                                                                   !!
!!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  !!
!!  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    !!
!!  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           !!
!!  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR  !!
!!  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   !!
!!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     !!
!!  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      !!
!!  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       !!
!!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    !!
!!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     !!
!!                                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program M3CBR
	use GOptions_
	use Math_
	use String_
	use IOStream_
	use Grid_
	use RNFunction_
	use NIntegrator_
	use Matrix_
	use RVector_
	use CommandLineParser_
	use BlocksIFileParser_
	
	use NNLS_
	implicit none
	
	integer :: i, j, k
	character(100), allocatable :: sTokens(:)
	character(1000) :: cArrBuffer
	real(8) :: rBuffer
	type(String) :: sBuffer
	type(Grid) :: bufferGrid
	logical :: lBuffer
	
	real(8) :: BR, rms, error
	real(8) :: averE, stdevE, skewE
	
	type(CommandLineParser) :: programOptions
	type(BlocksIFileParser) :: iParser
	
	type(String), allocatable :: basisTable(:)
	
	type(String) :: errorType
	type(String), allocatable :: expBRTable(:)
	type(String), allocatable :: expBRKey(:)
	real(8), allocatable :: expBRerrors(:)
	type(Matrix) :: R !< Experimental branching ratios matrix
	
	type(String) :: inputFileName
	type(String) :: dataFileName
	type(String) :: keysDataFileName
	type(String) :: energyDistFileName
	type(String) :: BRFileName
	integer :: nChannels
	type(IFStream) :: keysFile
	type(String), allocatable :: key(:)
	integer, allocatable :: idData2Exp(:)
	logical, allocatable :: isMappedExp2Data(:)
	
	integer :: n, l !< Global parameters which are used in basisFunction
	integer, allocatable :: NL(:,:) !< Global parameters which are used in basisFunction
	integer :: basisSize
	type(Grid) :: energyGrid !< Energy grid which is loaded from dataFileName
	type(RNFunction), allocatable :: P(:)  !< Channel probabilities which are loaded from dataFileName
	type(RNFunction), allocatable :: B(:)  !< Numerical representation for each basis function
	type(NIntegrator) :: integrator
	
	type(Matrix) :: A !< Integrals matrix
	type(Matrix) :: C !< Fitting coefficients matrix
	
	type(RNFunction) :: f !< Final energy distribution function

	!-------------------------------------------------------------------
	! Loading the command line options
	!-------------------------------------------------------------------
	call programOptions.init()
	
	if( command_argument_count() /= 6 ) then
		call usage()
		stop
	end if
	
	write(*,*) ""
	
	inputFileName = programOptions.getString( "-i" )
	write(*,"(A40,A)") "Input file = ", trim(inputFileName.fstr)
	
	dataFileName = programOptions.getString( "-d" )
	write(*,"(A40,A)") "Input data file = ", trim(dataFileName.fstr)
	
	keysDataFileName = programOptions.getString( "-k" )
	write(*,"(A40,A)") "Input keys data file = ", trim(keysDataFileName.fstr)
	
	write(*,*) ""
	
	!-------------------------------------------------------------------
	! Loading the data from the input file
	!-------------------------------------------------------------------
	call iParser.init( inputFileName.fstr )
	
	if( .not. iParser.isThereBlock( "BRANCHING_RATIOS" ) ) then
		write(*,*) "### ERROR ### M3CBR: BRANCHING_RATIOS block is required"
		stop
	end if
	
	energyDistFileName = iParser.getString( "BRANCHING_RATIOS:eDistOFile" )
	write(*,"(A40,A)"), "Energy distribution file = ", trim(energyDistFileName.fstr)
	
	BRFileName = iParser.getString( "BRANCHING_RATIOS:BRfile" )
	write(*,"(A40,A)"), "Branching ratios file = ", trim(BRFileName.fstr)
	
	!-------------------------------------------------------------------
	! Loading the basis set from the input file
	!-------------------------------------------------------------------
	call iParser.getBlock( "BRANCHING_RATIOS", basisTable )
	
	basisSize = size(basisTable)
	
	call C.init( basisSize, 1 )
	
	write(*,*) ""
	write(*,*) "Basis set"
	write(*,*) "---------"
	write(*,"(A)") " energy grid from "//trim(dataFileName.fstr)
	write(*,*) ""
	write(*,"(A15,I10)") "basis size = ", basisSize
	write(*,*) ""
	
	call energyGrid.fromFile( dataFileName.fstr, column=1 )

	allocate( NL(basisSize,2) )
	allocate( B(basisSize) )
	
	do i=1,basisSize
		call basisTable(i).split( sTokens, " " )
		call C.set( i, 1, FString_toReal( sTokens(1) ) )
		
		sBuffer = sTokens(2)
		call sBuffer.split( sTokens, "(,)" )
		
		n = FString_toInteger( sTokens(2) )
		l = FString_toInteger( sTokens(3) )
		
		NL(i,1) = n
		NL(i,2) = l
		
		call B(i).fromFunction( energyGrid, basisFunction )
		
		write(*,"(F15.5,2I5)") C.get(i,1), n, l
	end do
	write(*,*) ""
	
	!-------------------------------------------------------------------
	! Loading the keys
	!-------------------------------------------------------------------
	call keysFile.init( keysDataFileName.fstr )
	
	allocate( key(keysFile.numberOfLines) )
	
	do i=1,keysFile.numberOfLines
		key(i) = keysFile.readLine()
	end do
	
	write(*,*) "Available keys in "//trim(keysDataFileName.fstr)
	write(*,*) "----------------------------------------"
	do i=1,keysFile.numberOfLines
		write(*,"(I5,A50)") i, key(i).fstr
	end do
	write(*,*) ""
	
	if( iParser.isThereBlock( "EXPERIMENTAL_BRANCHING_RATIOS" ) ) then
		
		!-------------------------------------------------------------------
		! Loading the experimental values
		!-------------------------------------------------------------------
		errorType = iParser.getString( "EXPERIMENTAL_BRANCHING_RATIOS:error", "absolute" )
		
		call iParser.getBlock( "EXPERIMENTAL_BRANCHING_RATIOS", expBRTable )
		
		nChannels = size(expBRTable)
		
		call R.init( nChannels, 1 )
		allocate( expBRKey(nChannels) )
		allocate( expBRerrors(nChannels) )
		
		do i=1,nChannels
			read(expBRTable(i).fstr,*) cArrBuffer, rBuffer, expBRerrors(i)
			
			if( errorType == "relative" ) then
				expBRerrors(i) = rBuffer*expBRerrors(i)/100.0_8
			else if( errorType == "absolute" ) then
				expBRerrors(i) = expBRerrors(i)
			else
				write(*,*) "### ERROR ### unknown value for ''error'' parameter."
				write(*,*) "              available options: relative or absolute"
				stop
			end if
			
			expBRKey(i) = trim(cArrBuffer)
			call R.set( i, 1, rBuffer )
		end do
		
		deallocate(expBRTable)
		
		write(*,*) ""
		write(*,*) "Experimental branching ratios"
		write(*,*) "-----------------------------"
		do i=1,nChannels
			write(*,"(I5,A50,F15.6,A10,F15.6)") i, expBRKey(i).fstr, R.get(i,1), " +/- ", expBRerrors(i)
		end do
		
		write(*,*) ""
		write(*,"(A)") "Mapping channels"
		write(*,"(A)") "----------------"
		write(*,*) ""
		
		allocate( idData2Exp(keysFile.numberOfLines) )
		allocate( isMappedExp2Data(nChannels) )
		
		isMappedExp2Data = .false.
		do i=1,size(key)
			write(*,"(A)",advance="no") "Mapping "//trim(key(i).fstr)//" ... "
			
			lBuffer = .false.
			do j=1,size(expBRKey)
				if( trim(key(i).fstr) == trim(expBRKey(j).fstr) ) then
					write(*,"(A)") "OK ("//trim(FString_fromInteger(i))//":"//trim(FString_fromInteger(j))//")"
					idData2Exp(j) = i
					isMappedExp2Data(j) = .true.
					lBuffer = .true.
					exit
				end if
			end do
			
			if( .not. lBuffer ) then
				write(*,"(A)") "Failed"
	! 			write(*,*) ""
	! 			write(*,*) "### ERROR ### M3CfitBR: Key data inconsistency from both files"
	! 			write(*,*) "                "//trim(key(i).fstr)//" cannot be mapped"
	! 			stop
			end if
		end do
		write(*,*) ""
		
		do i=1,size(isMappedExp2Data)
			if( isMappedExp2Data(i) == .false. ) then
				write(*,"(A)") "@@@ WARNING @@@ "//trim(expBRKey(i).fstr)//" has not been mapped. Selecting zero probability in all range of energy"
			end if
		end do
		write(*,*) ""

	end if
	
	!-------------------------------------------------------------------
	! Loading the channel probabilities
	!-------------------------------------------------------------------
	
	write(*,*) "----------------------------------------"
	write(*,*) "Loading probabilities from "//trim(dataFileName.fstr)
	write(*,*) "----------------------------------------"
	
	call energyGrid.fromFile( dataFileName.fstr, column=1 )
	
	nChannels = keysFile.numberOfLines
	allocate( P(nChannels) )
	do i=1,size(P)
		write(*,"(A)",advance="no") "Reading channel "//trim(FString_fromInteger(i))//":"//trim(FString_fromInteger(2*i))//" ... "
		
		if( iParser.isThereBlock( "EXPERIMENTAL_BRANCHING_RATIOS" ) ) then
			if( isMappedExp2Data(i) .and. idData2Exp(i) /= 0 ) then
				call bufferGrid.fromFile( dataFileName.fstr, column=2*idData2Exp(i) )
			else
				bufferGrid = energyGrid
				bufferGrid.data = 0.0_8
			end if
		else
			call bufferGrid.fromFile( dataFileName.fstr, column=2*i )
		end if
		
		call P(i).fromGridArray( energyGrid, bufferGrid.data )
		write(*,"(A)") "OK"
	end do
	
	!-------------------------------------------------------------------
	! Final report
	!-------------------------------------------------------------------
	write(*,*) ""
	write(*,*) "----------------------------------------"
	write(*,*) "Final branching ratios ("//trim(BRFileName.fstr)//")"
	write(*,*) "----------------------------------------"
	
	open(11,file=BRFileName.fstr)
	
	write(*,"(35X,2A15,A30,A15)") "   ", "      ", "    Experimental     ", "         "
	write(*,"(35X,5A15)")         "key", "  BR  ", "  BR  ", "  error  ", "  delta  "
	write(*,"(35X,5A15)")         "---", "------", "------", "---------", "---------"

	write(11,"(A)") "# Branching ratios fitted with M3CBR"
	write(11,"(A1,27X,7X,2A15)") "#", "key", "BR"
	write(11,"(A1,27X,7X,2A15)") "#", "---", "---"
	
	call A.init( nChannels, basisSize )
	
	rms = 0.0_8
	do i=1,nChannels
		
		BR = 0.0_8
		do k=1,basisSize
			f = B(k)*P(i)
			call integrator.init( f, NIntegrator_BOOLE )
			BR = BR + C.get(k,1)*integrator.evaluate()
		end do
		
		if( iParser.isThereBlock( "EXPERIMENTAL_BRANCHING_RATIOS" ) ) then
			
			if( isMappedExp2Data(i) .and. idData2Exp(i) /= 0 ) then
				error = BR - R.get(i,1)
				
				write(*,"(A50,4F15.5)") key(idData2Exp(i)).fstr, BR, R.get(i,1), expBRerrors(i), error
				write(11,"(A50,4F15.5)") key(idData2Exp(i)).fstr, BR, R.get(i,1), expBRerrors(i), error
! 			else
! 				write(*,"(A50,4F15.5,A)") key(idData2Exp(i)).fstr, BR, 0.0, 0.0, 0.0, " *"
! 				write(11,"(A50,4F15.5,A)") key(idData2Exp(i)).fstr, BR, 0.0, 0.0, 0.0, " *"
			end if
			
			rms = rms + error**2
		else
			write(*,"(A50,F15.5)") key(i).fstr, BR, error
			write(11,"(A50,F15.5)") key(i).fstr, BR, error
		end if
		
	end do
	rms = sqrt(rms/nChannels)
	
	close(11)
	
	! Maxima:
	! B(E):=E**l*exp(-E/n)/(n**(l+1)*gamma(l+1));
	! integrate(E*B(E),E,0,inf); #---> (l+1)*n
	! integrate((E-mu)**2*B(E),E,0,inf); #---> (l+1)*(l+2)*n**2-2*(l+1)*n*mu+mu**2
	! integrate(((E-mu)/sigma)**3*B(E),E,0,inf); #--> ( (l+1)*(l+2)*(l+3)*n**3-3*(l+1)*(l+2)*n**2*mu+3*(l+1)*n*mu**2-mu**3 )/sigma**3
	
	averE = 0.0_8
	do k=1,basisSize
		n = NL(k,1)
		l = NL(k,2)
		
		averE = averE + (C.get(k,1)/100.0_8)*real((l+1)*n,8)
	end do
	
	stdevE = 0.0_8
	do k=1,basisSize
		n = NL(k,1)
		l = NL(k,2)
		
		stdevE = stdevE + (C.get(k,1)/100.0_8)*( real((l+1)*(l+2)*n**2,8) - 2.0_8*real((l+1)*n,8)*averE + averE**2 )
	end do
	stdevE = sqrt(stdevE)
	
	skewE = 0.0_8
	do k=1,basisSize
		n = NL(k,1)
		l = NL(k,2)
		
		skewE = skewE + (C.get(k,1)/100.0_8)*( real((l+1)*(l+2)*(l+3)*n**3,8) - 3.0_8*real((l+1)*(l+2)*n**2,8)*averE + 3.0_8*real((l+1)*n,8)*averE**2-averE**3 )/stdevE**3
	end do
	
	call f.fromFunction( energyGrid, energyFunction )
	call f.save( energyDistFileName.fstr )
	call integrator.init( f, NIntegrator_BOOLE )
	write(*,*) ""
	write(*,"(A15,F10.5)") "rms = ", rms
	write(*,"(A15,F10.5)") "Integral = ", integrator.evaluate()
	write(*,"(A15,F10.5)") "<E> = ", averE
	write(*,"(A15,F10.5)") "stdev(E) = ", stdevE
	write(*,"(A15,F10.5)") "skew(E) = ", skewE
	
	if( integrator.evaluate() < 98.0_8 ) then
		write(*,*) ""
		write(*,*) " ### ERROR ### The energy distribution function is not correctly normalized"
		write(*,*) ""
		stop
	end if
	
	contains
	
	!>
	!! 
	!!
	subroutine usage()
		write(*,"(A)") "NAME"
		write(*,"(A)") "        M3CBR - Program to calculate the branching ratios an energy distribution from M3C program results"
		write(*,"(A)") ""
		write(*,"(A)") "SYNOPSYS"
		write(*,"(A)") "        M3CBR -i File -d File -k File"
		write(*,"(A)") ""
		write(*,"(A)") "DESCRIPTION"
		write(*,"(A)") "        M3CBR -i File -d File -k File"
		write(*,"(A)") ""
		write(*,"(A)") "OPTIONS"
		write(*,"(A)") "        M3CBR accepts the following options."
		write(*,"(A)") ""
		write(*,"(A)") "        -i File"
		write(*,"(A)") "                M3C iput file containing a section like this:"
		write(*,"(A)") ""
		write(*,"(A)") "                BEGIN BRANCHING_RATIOS"
		write(*,"(A)") "                		method = NNLS"
		write(*,"(A)") "                		basis = 5,5"
		write(*,"(A)") "                		eDistOFile = edist.out"
		write(*,"(A)") "                		BRfile = fitBR.out"
		write(*,"(A)") ""
		write(*,"(A)") "                	#---------------------------------------"
		write(*,"(A)") "                	#  Channel                    BR   error"
		write(*,"(A)") "                	#---------------------------------------"
		write(*,"(A)") "                	   C_5                      15.0     1.0"
		write(*,"(A)") "                	   C_1+C_4                   8.0     1.0"
		write(*,"(A)") "                	   C_2+C_3                  59.0     3.0"
		write(*,"(A)") "                	   C_1+C_1+C_3               8.0     1.0"
		write(*,"(A)") "                	   C_1+C_2+C_2               7.0     1.0"
		write(*,"(A)") "                	   C_1+C_1+C_1+C_2           2.5     0.6"
		write(*,"(A)") "                	   C_1+C_1+C_1+C_1+C_1      0.27     0.2"
		write(*,"(A)") "                END BRANCHING_RATIOS"
		write(*,"(A)") "                "
		write(*,"(A)") "        -d File"
		write(*,"(A)") "                Channel or species energy curves in column format"
		write(*,"(A)") "                <Energy>    <prob1>  <error1>   <prob2>  <error2>  ..."
		write(*,"(A)") "                "
		write(*,"(A)") "        -k File"
		write(*,"(A)") "                Key file in agreement with -d option"
		write(*,"(A)") "                "
		write(*,"(A)") "AUTHORS"
		write(*,"(A)") "                "
		write(*,"(A)") "EXAMPLES"
		write(*,"(A)") "                "
		write(*,"(A)") "ACKNOWLEDGEMENTS"
	end subroutine usage
	
	!>
	!! Representa la forma analítica de las funciones de base.
	!! Los valores de n y l son compartidos con el programa principal
	!!
	function basisFunction( E ) result( output )
		real(8), intent(in) :: E
		real(8) :: output
		
		output = (1.0_8/gamma(real(l+1,8)))*(1.0_8/real(n,8))**(l+1)*E**l*exp(-E/real(n,8))
	end function basisFunction
	
	!>
	!! Representa la forma analítica de la función de depósito de energía.
	!! Los valores de C, k, n y l son compartidos con el programa principal
	!!
	function energyFunction( E ) result( output )
		real(8), intent(in) :: E
		real(8) :: output
		
		output = 0.0_8
		
		do k=1,basisSize
			n = NL(k,1)
			l = NL(k,2)
			
			output = output + C.get(k,1)*basisFunction( E )
		end do
	end function energyFunction
	
end program M3CBR
