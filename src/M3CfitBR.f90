program M3CfitBR
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
	character(1000) :: cArrBuffer
	real(8) :: rBuffer, ssum
	character(10), allocatable :: sTokens(:)
	type(String) :: sBuffer
	type(Grid) :: bufferGrid
	logical :: lBuffer
	
	type(CommandLineParser) :: programOptions
	type(BlocksIFileParser) :: iParser
	
	integer :: method ! 0=NNLS and 1=NNWLS
	type(String) :: inputFileName
	type(String) :: dataFileName
	type(String) :: keysDataFileName
	type(String) :: energyDistFileName
	type(String) :: BRFileName

	type(String) :: errorType
	type(String), allocatable :: expBRTable(:)
	type(String), allocatable :: expBRKey(:)
	real(8), allocatable :: expBRerrors(:)
	integer :: nChannels
	type(IFStream) :: keysFile
	type(String), allocatable :: key(:)
	integer, allocatable :: idData2Exp(:)
	logical, allocatable :: isMappedExp2Data(:)
	
	integer :: n, l !< Global parameters which are used in basisFunction
	integer :: Nmax, Lmax
	integer :: basisSize
	type(Grid) :: energyGrid !< Energy grid which is loaded from dataFileName
	type(RNFunction), allocatable :: P(:)  !< Channel probabilities which are loaded from dataFileName
	type(RNFunction), allocatable :: B(:)  !< Numerical representation for each basis function
	type(NIntegrator) :: integrator
	
	type(Matrix) :: A !< Integrals matrix
	type(Matrix) :: C !< Fitting coefficients matrix
	type(Matrix) :: R !< Experimental branching ratios matrix
	
	type(Matrix) :: W  !< weight matrix
	type(Matrix) :: A_ !< Necessary types for NNLS method
	type(Matrix) :: R_
	integer :: info
	real(8), allocatable :: work(:)
	integer, allocatable :: indx(:)
	
	type(RNFunction) :: f !< Final energy distribution function
	real(8) :: sumExp, sumFit, sumError, norm
	
	call system( "rm -rf edist.out fitBR.out" )
	
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
	!    method: NNLS or NNWLS
	!    basis to use
	!    experimental branching ratios and its errors
	!    output files for energy distribution and BR fit results
	!-------------------------------------------------------------------
	call iParser.init( inputFileName.fstr )
	
	if( .not. iParser.isThereBlock( "FIT_BRANCHING_RATIOS" ) ) then
		write(*,*) "### ERROR ### M3CfitBR: FIT_BRANCHING_RATIOS block is required"
		stop
	end if
	
	if( .not. iParser.isThereBlock( "EXPERIMENTAL_BRANCHING_RATIOS" ) ) then
		write(*,*) "### ERROR ### M3CfitBR: EXPERIMENTAL_BRANCHING_RATIOS block is required"
		stop
	end if
	
	sBuffer = iParser.getString( "FIT_BRANCHING_RATIOS:method" )
	if( trim(sBuffer.fstr) == "NNLS" ) then
		method = 0
	else if( trim(sBuffer.fstr) == "NNWLS" ) then
		method = 1
	else
		write(*,*) "### ERROR ### M3CfitBR: Method "//trim(sBuffer.fstr)//" is not available"
		stop
	end if
	write(*,"(A40,A)"), "method = ", trim(sBuffer.fstr)
	
	sBuffer = iParser.getString( "FIT_BRANCHING_RATIOS:basis" )
	
	call sBuffer.split( sTokens, "," )
	
	if( size(sTokens) /= 2 ) then
		write(*,*) "### ERROR ### M3CfitBR: Error in basis set definition"
		stop
	end if
	
	Nmax = FString_toInteger( sTokens(1) )
	Lmax = FString_toInteger( sTokens(2) )
	write(*,"(A40,2I7)"), "basis = ", Nmax, Lmax
		
	energyDistFileName = iParser.getString( "FIT_BRANCHING_RATIOS:eDistfile" )
	write(*,"(A40,A)"), "Energy distribution file = ", trim(energyDistFileName.fstr)
	
	BRFileName = iParser.getString( "FIT_BRANCHING_RATIOS:BRfile" )
	write(*,"(A40,A)"), "Branching ratios file = ", trim(BRFileName.fstr)
	
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
		write(*,"(A50,F15.6,A10,F15.6)") expBRKey(i).fstr, R.get(i,1), " +/- ", expBRerrors(i)
	end do
	
	!-------------------------------------------------------------------
	! Loading the keys
	! El fichero contiene las columnas en el formato:
	! <energy>   <P1>  <S1>   <P2>  <S2>    ...   <PN>  <SN>
	!             ^            ^                   ^
	!           col=2        col=4               col=2*n
	!
	! conde Pi son las probabilidades y Si las desviaciones estandar
	!-------------------------------------------------------------------
	call keysFile.init( keysDataFileName.fstr )
	
	allocate( idData2Exp(keysFile.numberOfLines) )
	allocate( key(keysFile.numberOfLines) )
	
	do i=1,keysFile.numberOfLines
		key(i) = keysFile.readLine()
	end do
	
! 	if( size(key) /= size(expBRKey) ) then
! 		write(*,*) ""
! 		write(*,*) "### ERROR ### M3CfitBR: Key data inconsistency from both files"
! 		write(*,*) ""
! 		
! 		write(*,*) "----------------------------------------"
! 		write(*,*) "Available keys in "//trim(keysDataFileName.fstr)
! 		write(*,*) "----------------------------------------"
! 		do i=1,keysFile.numberOfLines
! 			write(*,"(I5,A50)") i, key(i).fstr
! 		end do
! 		
! 		write(*,*) "----------------------------------------"
! 		write(*,*) "Available keys in "//trim(inputFileName.fstr)
! 		write(*,*) "----------------------------------------"
! 		do i=1,size(expBRKey)
! 			write(*,"(I5,A50)") i, expBRKey(i).fstr
! 		end do
! 		
! 		stop
! 	end if
	
	write(*,*) ""
	write(*,"(A)") "Mapping channels"
	write(*,"(A)") "----------------"
	write(*,*) ""
	
	!----------------------------------------------------------------------
	! Mapea los canales que estan en el calculo, respecto a los del input.
	! Si hay uno que encuentra que hay uno en el calculo, pero no en el
	! fichero de entrada, lo reporta como Failed, pero continua con los
	! que tiene
	!----------------------------------------------------------------------
	allocate( isMappedExp2Data(nChannels) )
	isMappedExp2Data = .false.
	do i=1,size(key)
		write(*,"(A50,A)",advance="no") trim(key(i).fstr), " ... "
		
		lBuffer = .false.
		do j=1,size(expBRKey)
			if( trim(key(i).fstr) == trim(expBRKey(j).fstr) ) then
				idData2Exp(j) = i
				isMappedExp2Data(j) = .true.
				lBuffer = .true.
				
				write(*,"(A10,A)") "OK       (", trim(FString_fromInteger(j))//":"//trim(FString_fromInteger(2*i))//")"
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
			write(*,*) "### ERROR ### M3CfitBR: "//trim(expBRKey(i).fstr)//" has not been mapped."
			write(*,*) "                        Remove it from the input file"
			write(*,*) ""
			stop
		end if
	end do
	
	!-------------------------------------------------------------------
	! Loading the channel probabilities
	!-------------------------------------------------------------------
	call energyGrid.fromFile( dataFileName.fstr, column=1 )
	
	allocate( P(nChannels) )
	do i=1,nChannels
		write(*,"(A)",advance="no") "Reading channel "//trim(FString_fromInteger(i))//":"//trim(FString_fromInteger(2*i))//" ... "
				
		call bufferGrid.fromFile( dataFileName.fstr, column=2*idData2Exp(i) )
		call P(i).fromGridArray( energyGrid, bufferGrid.data )
		write(*,"(A)") "OK"
	end do
	
	!-------------------------------------------------------------------
	! Building the numerical functions for basis
	!-------------------------------------------------------------------
	k=1
	do n=1,Nmax
		do l=1,Lmax
			k = k+1
		end do
	end do
	basisSize = k-1
	
	write(*,*) ""
	write(*,"(A15,I10)") "basis size = ", basisSize
	write(*,*) ""
	
	allocate( B(basisSize) )
	
	k=1
	do n=1,Nmax
		do l=1,Lmax
			call B(k).fromFunction( energyGrid, basisFunction )
			k=k+1
		end do
	end do
	
	!-------------------------------------------------------------------
	! Building the integrals matrix
	!-------------------------------------------------------------------
	call A.init( nChannels, basisSize )
	
	do i=1,nChannels
		do k=1,basisSize
			f = B(k)*P(i)
			call integrator.init( f, NIntegrator_BOOLE )
			call A.set(i, k, integrator.evaluate() )
		end do
	end do
	
	!-------------------------------------------------------------------
	! Building the matrices for NNLS method
	!-------------------------------------------------------------------
	call C.init( basisSize+1, 1 )
	
	allocate( indx(basisSize) )
	allocate( work(basisSize) )
	
	call A_.init( A.nRows+1, A.nCols )
	A_.data(1:A.nRows,1:A.nCols) = A.data(1:A.nRows,1:A.nCols)
	A_.data(A.nRows+1,:) = 1.0_8
	
	call R_.init( R.nRows+1, R.nCols )
	R_.data(1:R.nRows,1:R.nCols) = R.data(1:R.nRows,1:R.nCols)
	R_.data(R.nRows+1,1:R.nCols) = 100.0_8
	
	call W.init( A.nRows+1, A.nRows+1, 0.0_8 )
	
	if( method == 0 ) then !< method NNLS
		do i=1,nChannels
			call W.set( i, i, 1.0_8 )
		end do
	else if( method == 1 ) then !< method NNWLS
		do i=1,nChannels
			call W.set( i, i, 1.0_8/expBRerrors(i)**2 )
		end do
	end if
	
	call W.set( A.nRows+1, A.nRows+1, 1d8 )
	
	A_ = W*A_
	R_ = W*R_
		
	call NNLS( A_.data, nChannels+1, basisSize, R_.data(:,1), C.data(:,1), rBuffer, work, indx, info )
	
	if( info == 1 ) then
		write(*,"(A)") " >>> Successful fitting <<<"
	else
		write(*,"(A)") " ### ERROR ### Failed to converge fitting"
		stop
	end if
	write(*,*) ""
	
	!-------------------------------------------------------------------
	! Final report
	!-------------------------------------------------------------------
	write(*,*) "Final energy distibution funtion"
	write(*,*) "--------------------------------"
	write(*,*) ""
	
	ssum = 0.0_8
	k=1
	do n=1,Nmax
		do l=1,Lmax
			if( C.data(k,1) > 1e-6 ) then
				write(*,"(I7,F15.6,A)") k, C.data(k,1), "  B("//trim(FString_fromInteger(n))//","//trim(FString_fromInteger(l))//")"
				ssum = ssum + C.data(k,1)
			end if
			
			k = k+1
		end do
	end do
	
	write(*,"(7X,A15)") "------------"
	write(*,"(7X,F15.6)") ssum
	write(*,*) ""

	write(*,*) "Final branching ratios"
	write(*,*) "----------------------"
	write(*,*) ""
	
	open(11,file=BRFileName.fstr)
	
	write(*,"(35X,4A15)") "key", "fitted", "exact", "error"
	write(*,"(35X,4A15)") "---", "------", "-----", "-----"
	
	write(11,"(A)") "# Branching ratios fitted with M3CfitBR"
	write(11,"(A1,27X,7X,4A15)") "#", "key", "fit.BR", "exp.BR", "exp.BR.error"
	write(11,"(A1,27X,7X,4A15)") "#", "---", "------", "------", "------------"
	
	norm = 0.0_8
	do i=1,nChannels
		norm = norm + sum(A.data(i,:)*C.data(:,1))
	end do
	C.data = 100.0_8*C.data/norm
	
	ssum = 0.0_8
	sumExp = 0.0_8
	sumFit = 0.0_8
	sumError = 0.0_8
	do i=1,nChannels
		write(*,"(A50,3F15.5)") expBRKey(i).fstr, sum(A.data(i,:)*C.data(1:basisSize,1)), R.get(i,1), ( sum(A.data(i,:)*C.data(1:basisSize,1))-R.get(i,1) )
		
		write(11,"(A50,3F15.5)") expBRKey(i).fstr, sum(A.data(i,:)*C.data(1:basisSize,1)), R.get(i,1), expBRerrors(i)
		
		sumFit = sumFit + sum(A.data(i,:)*C.data(:,1))
		sumExp = sumExp + R.get(i,1)
		sumError = sumError + ( sum(A.data(i,:)*C.data(1:basisSize,1))-R.get(i,1) )
		
		if( method == 0 ) then !< method NNLS
			ssum = ssum + ( sum(A.data(i,:)*C.data(1:basisSize,1))-R.get(i,1) )**2
		else if( method == 1 ) then !< method NNWLS
			ssum = ssum + ( sum(A.data(i,:)*C.data(1:basisSize,1))-R.get(i,1) )**2/expBRerrors(i)**2
		end if
	end do
	write(*,"(35X,4A15)") "", "------", "-----", "-----"
	write(*,"(A50,3F15.5)") "", sumFit, sumExp, sumError
	
	write(*,*) ""
	write(*,"(A15,F10.5)") "rms = ", sqrt(ssum)
	close(11)
	
	call f.fromFunction( energyGrid, energyFunction )
	call f.save( energyDistFileName.fstr )
	call integrator.init( f, NIntegrator_BOOLE )
	write(*,"(A15,F10.5)") "Integral = ", integrator.evaluate()
	
	if( integrator.evaluate() < 98.0_8 ) then
		write(*,*) ""
		write(*,*) " ### ERROR ### The energy distribution function is not correctly normalized"
		write(*,*) ""
		stop
	end if
	
	contains
	
	!>
	!! Representa la forma analítica de las funciones de base.
	!! Los valores de n y l son compartidos con el programa principal
	!!
	subroutine usage()
		write(*,"(A)") "NAME"
		write(*,"(A)") "        M3CfitBR - Program to calculate the branching ratios an energy distribution from M3C program results"
		write(*,"(A)") ""
		write(*,"(A)") "SYNOPSYS"
		write(*,"(A)") "        M3CfitBR -i File -d File -k File"
		write(*,"(A)") ""
		write(*,"(A)") "DESCRIPTION"
		write(*,"(A)") "        M3CfitBR -i File -d File -k File"
		write(*,"(A)") ""
		write(*,"(A)") "OPTIONS"
		write(*,"(A)") "        M3CfitBR accepts the following options."
		write(*,"(A)") ""
		write(*,"(A)") "        -i File"
		write(*,"(A)") "                M3C iput file containing a section like this:"
		write(*,"(A)") ""
		write(*,"(A)") "                BEGIN FIT_BRANCHING_RATIOS"
		write(*,"(A)") "                		method = NNLS"
		write(*,"(A)") "                		basis = 5,5"
		write(*,"(A)") "                		eDistfile = edist.out"
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
		write(*,"(A)") "                END FIT_BRANCHING_RATIOS"
		write(*,"(A)") "                "
		write(*,"(A)") "        -d File"
		write(*,"(A)") "                The energy distribution function will be saved in this file"
		write(*,"(A)") "                "
		write(*,"(A)") "        -k File"
		write(*,"(A)") "                The branching ratios will be saved in this file"
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
		
		k=1
		do n=1,Nmax
			do l=1,Lmax
				output = output + C.get(k,1)*basisFunction( E )
				
				k = k+1
			end do
		end do
	end function energyFunction
	
end program M3CfitBR
