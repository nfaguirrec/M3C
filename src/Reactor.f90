!>
!! @brief
!!
module Reactor_
	use Math_
	use String_
	use IOStream_
	use Timer_
	use UnitsConverter_
	use IVector_
	use RandomUtils_
	use StringList_
	use AtomicElementsDB_
	use BlocksIFileParser_
	use RealHistogram_
	use StringIntegerMap_
	use StringRealMap_
	use StringRealPair_
	use StringHistogram_
	use GOptions_
	
	use Fragment_
	use FragmentsList_
	use FragmentsDB_
	
	implicit none
	private
	
	enum, BIND(c)
		enumerator :: NULL_REACTOR
		enumerator :: STRUCTURE_REACTOR
		enumerator :: TRANSLATIONAL_REACTOR
		enumerator :: ROTATIONAL_REACTOR
		enumerator :: VIBRATIONAL_REACTOR
	end enum
	
	character(3), public :: REACTOR_LABEL(4) = [ 'S', 'T', 'R', 'V' ]
	
	type, public :: Reactor
		type(FragmentsList) :: reactives
		type(FragmentsList) :: products
		logical :: state
		
		character(3), private :: name
		integer, private :: type
		integer, allocatable :: dNFrag(:)
		
		contains
			procedure :: init ! To overwride
			procedure :: initReactor
			generic :: assignment(=) => copyReactor
			
			procedure :: copyReactor
			final :: destroyReactor
			procedure :: showReactorHeader
			procedure :: setType
			procedure :: run
			procedure, private :: changeComposition
			
			procedure :: execute
			procedure :: executeMinFragmentationEnergy
	end type Reactor
	
	integer, private :: internalReactivesComposition(AtomicElementsDB_nElems) = 0 !< For convenience in changeAllComposition procedure
	integer, private :: internalMassNumber = 0
	integer, private :: internalCharge = 0
	integer, private :: internalNTrials = 0
	
	contains
	
	!>
	!! @brief Constructor
	!!
	subroutine init( this, reactives, excitationEnergy )
		class(Reactor) :: this
		type(FragmentsList), intent(in) :: reactives
		real(8), intent(in) :: excitationEnergy
		
		call this.initReactor( reactives, excitationEnergy )
		call this.reactives.initialGuessFragmentsList()
	end subroutine init
	
	!>
	!! @brief Constructor
	!!
	subroutine initReactor( this, reactives, excitationEnergy )
		class(Reactor) :: this
		type(FragmentsList), intent(in) :: reactives
		real(8), intent(in) :: excitationEnergy
		
		real(8) :: rBuffer1, rBuffer2
		
		this.name = "R"
		this.state = .true.
		
		this.reactives = reactives
		
		! Este es el cero de energía
		call this.reactives.setReactorEnergy( excitationEnergy )
		
		if( GOptions_printLevel >= 2 ) then
			call GOptions_section( "REACTOR INITIALIZATION" )
			
			call GOptions_valueReport( "Eelec", reactives.electronicEnergy()/eV, "eV", indent=1 )
			call GOptions_valueReport( "Eexcit", excitationEnergy/eV, "eV", indent=1 )
			call GOptions_valueReport( "Reacts", trim(this.reactives.label()), indent=1 )
		end if
		
! 		call this.reactives.buildInitialConfiguration()
		
		! La energía total es igual a la energía de excitación
! 		call this.reactives.setTranslationalEnergy( this.reactives.reactorEnergy() - this.reactives.internalEnergy() )
		
		if( GOptions_printLevel >= 2 ) then
			call GOptions_paragraph( "Energy balance", indent=1 )
			call GOptions_valueReport( "Reactor energy", this.reactives.reactorEnergy()/eV, "eV", indent=1 )
			call GOptions_valueReport( "Int. energy", this.reactives.internalEnergy()/eV, "eV", indent=1 )
			call GOptions_valueReport( "T energy", this.reactives.kineticEnergy()/eV, "eV", indent=1 )
			write(*,*) ""
		end if
		
	end subroutine initReactor
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copyReactor( this, other )
		class(Reactor), intent(out) :: this
		class(Reactor), intent(in) :: other
		
		call GOptions_error( &
			"This method is not implemented", &
			"Reactor.copyReactor()" &
		)
	end subroutine copyReactor
	
	!>
	!! @brief Destructor
	!!
	subroutine destroyReactor( this )
		type(Reactor) :: this
		
		if( allocated(this.dNFrag) ) deallocate( this.dNFrag )
		
	end subroutine destroyReactor
	
	!>
	!! @brief Set the reactor type by using strings like these V, R(10), T(10), S:0, S:1:-1
	!!
	subroutine setType( this, strId )
		class(Reactor) :: this
		character(*), intent(in) :: strId
		
		integer :: i
		character(10), allocatable :: tokens(:)
		
		this.name = trim(strId)
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! En caso que el strId sea de un reactor de estructura
		! se obtiene el cambio en el número de fragmentos
		call FString_split( strId, tokens, ":" )
		
		if( size(tokens) > 1 ) then
			
			if( allocated(this.dNFrag) ) deallocate(this.dNFrag)
			allocate(this.dNFrag(size(tokens)-1))
			
			do i=2,size(tokens)
				this.dNFrag(i-1) = FString_toInteger( tokens(i) )
			end do
		end if
		
		do i=1,size(REACTOR_LABEL)
			if( trim(REACTOR_LABEL(i)) == trim( tokens(1) ) ) then
				this.type = i
				return
			end if
		end do
		
		call GOptions_error( &
			"Unknown reactor id="//trim(FString_fromInteger(i))//" or strId="//trim(strId), &
			"Reactor.setType()" &
		)
	end subroutine setType
	
	!>
	!! @brief Shows the reactor header
	!!
	subroutine showReactorHeader( this )
		class(Reactor), intent(in) :: this
		
		if( GOptions_printLevel >= 2 ) then
			call GOptions_section( "REACTOR STARTS RUNNING", indent=1 )
			
			select case( this.type )
				case( STRUCTURE_REACTOR )
					call GOptions_valueReport( "type", "STRUCTURE", indent=1 )
					call GOptions_valueReport( "dNFrag", this.dNFrag, indent=1 )
				case( TRANSLATIONAL_REACTOR )
					call GOptions_valueReport( "type", "TRANSLATIONAL", indent=1 )
				case( ROTATIONAL_REACTOR )
					call GOptions_valueReport( "type", "ROTATIONAL", indent=1 )
				case( VIBRATIONAL_REACTOR )
					call GOptions_valueReport( "type", "VIBRATIONAL", indent=1 )
			end select
			
			write(*,*) ""
		end if
	end subroutine showReactorHeader
	
	!>
	!! @brief Change the composition of the system
	!!
	subroutine changeComposition( this, dNfrag )
		class(Reactor), intent(inout) :: this
		integer, intent(in) :: dNfrag(:)
		
		integer :: n
		
		n = RandomUtils_uniform( [ 1, size(this.dNFrag) ] )
		
		if( GOptions_printLevel >= 2 ) then
			call GOptions_paragraph("CHANGE COMPOSITION")
			
			write(*,"(A10,I20,5X,A)") "dN", this.dNFrag(n), "used for reactor"
		end if
		
		call changeAllComposition( this.reactives, this.products, this.dNFrag(n) )
		
		return
	end subroutine changeComposition
	
	!>
	!! This is only necessary for changeAllComposition
	!!
	function reactorConstraint( multisetPositions, current ) result( output )
		integer, allocatable, intent(in) :: multisetPositions(:)
		integer, intent(in) :: current
		logical :: output
		
		integer :: i
		real(8) :: totalMass
		integer :: totalCharge
		integer :: productsComposition( AtomicElementsDB_nElems )
		
		totalMass = 0.0_8
		totalCharge = 0
		productsComposition = 0
		do i=1,current
				totalMass = totalMass + FragmentsDB_instance.clusters( multisetPositions(i) ).mass()
				totalCharge = totalCharge + FragmentsDB_instance.clusters( multisetPositions(i) ).charge
				productsComposition = productsComposition + FragmentsDB_instance.clusters( multisetPositions(i) ).composition
		end do
		
		if( current == size(multisetPositions) ) then
			!----------------------------------------------
			! ¿ El multiset encontrado es correcto ?
			!----------------------------------------------
			output = .false.
			
			! @warning Acá estoy asumiendo que los reactivos tienen actualizada la formula y así su composición
			if( all( productsComposition == internalReactivesComposition ) .and. ( internalCharge == totalCharge ) ) then
				output = .true.
			end if
			
			internalNTrials = internalNTrials + 1  ! Por el momento este valor es solo informativo
		else
			!----------------------------------------------
			! ¿ El camino recorrido es incorrecto ?
			!----------------------------------------------
			output = .true.
			
			if( all( productsComposition <= internalReactivesComposition ) .and. totalCharge <= internalCharge ) then
				output = .false.
			end if
		end if
	end function reactorConstraint
	
	!>
	!! @brief Change the composition of the system
	!!
	subroutine changeAllComposition( reactives, products, dNfrag )
		type(FragmentsList), intent(in) :: reactives
		type(FragmentsList), intent(inout) :: products
		integer, intent(in) :: dNfrag
		
		integer, allocatable :: ids(:)
		integer, allocatable :: channelInfo(:) ! Ids para camino de reacción aleatorio
		integer :: nChannels ! number of channels
		integer :: nProducts ! number of products in one channel
		integer :: i
		
		nProducts = reactives.nMolecules() + dNfrag
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Si el cluster no se puede fragmentar mas, mantenga los reactivos
		if( nProducts > reactives.nAtoms() ) then
			if( GOptions_printLevel >= 2 ) then
				call GOptions_info( &
					"The fragmentation limit has been reached", &
					"Reactor.changeAllComposition()", &
					"The reactives composition is kept." &
				)
			end if
			
			products = reactives
			return
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Si los clusters no se pueden fusionar mas, mantenga los reactivos
		if( nProducts < 1 ) then
			if( GOptions_printLevel >= 2 ) then
				call GOptions_info( &
					"The fussion limit has been reached", &
					"Reactor.changeAllComposition()", &
					"The reactives composition is kept." &
				)
			end if
			
			products = reactives
			return
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se reserva la memoria necesaria para almacenar
		! todos los canales
		nChannels = Math_multisetNumber( FragmentsDB_instance.nMolecules(), nProducts ) ! El tamaño es multiset( N_db, N_prod )
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Si solo hay un canal es porque los reactivos
		! son los mismos productos
! 		if( nChannels <= 1 ) then
! 			products = reactives
! 			return
! 		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se calcula la masa y carga total las cuales se utilizarán
		! en las retricciones de conservación
		allocate( ids(FragmentsDB_instance.nMolecules()) )
		
		internalMassNumber = 0
		internalCharge = 0
		internalReactivesComposition = 0
		do i=1,reactives.nMolecules()
			internalMassNumber = internalMassNumber + reactives.clusters(i).massNumber()
			internalCharge = internalCharge + reactives.clusters(i).charge
			internalReactivesComposition = internalReactivesComposition + reactives.clusters(i).composition
		end do
		
		internalNTrials = 0
		
		do i=1,FragmentsDB_instance.nMolecules()
			ids(i) = i
		end do
		
		call RandomUtils_randomMultiset( ids, nProducts, channelInfo, reactorConstraint )
		deallocate(ids)
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se muestran lo valores importantes
		if( GOptions_printLevel >= 2 ) then
			write(*,"(A10,I20,5X,A)") "massNumber", internalMassNumber, "used for reactor"
			write(*,"(A10,I20,5X,A)") "charge", internalCharge, "used for reactor"
			write(*,"(A10,I20,5X,A)") "nTrials", internalNTrials, "used for reactor"
		end if
		
		call products.init( nProducts )
		
		do i=1,nProducts
			call products.set( i, FragmentsDB_instance.clusters( channelInfo(i) ) )
		end do
		
		if( GOptions_printLevel >= 2 ) then
			write(*,*) ""
			write(*,"(5X,A)")      "Choosen channel: "
			write(*,"(5X,A,5X,A)") "                 ", trim(reactives.label())//" --> "//trim(products.label())
			write(*,*) ""
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se libera la memoria solicitada
		deallocate( channelInfo )
		
	end subroutine changeAllComposition

	!>
	!! @brief Override run method
	!!
	subroutine run( this )
		class(Reactor), intent(inout) :: this
		
		real(8) :: min_dLnW
		
		integer :: n, nMin
		real(8) :: rBuffer
		type(String) :: sBuffer
		
		this.state = .true.
		
		call this.showReactorHeader()
		
		select case( this.type )
			case( STRUCTURE_REACTOR )
			
				! Se actualiza la composición
				call this.changeComposition( this.dNFrag )
				
				! Se le asocia la energía del reactor para asegurar que calcula un peso Wt es adecuado para los productos
				call this.products.setReactorEnergy( this.reactives.reactorEnergy() )
				
				! Para que fuerce los centros aleatorios en la siguiente iteración
				this.products.forceRandomCenters = .true.
	! 			call this.products.initialGuessFragmentsList()
				
				! Los productos utilizan parte de la energía
				call this.products.changeVibrationalEnergy()
				call this.products.changeOrientations()
! 				call this.products.changeAngularMomenta()
				call this.products.changeGeometry()
				
			case( VIBRATIONAL_REACTOR )
				
				! La composición es igual antes y después
				this.products = this.reactives
				
				! Los productos utilizan parte de la energía cinética en vibracional
				call this.products.changeVibrationalEnergy()
			
			case( TRANSLATIONAL_REACTOR )
				! La composición es igual antes y después
				this.products = this.reactives
				
				! Los productos utilizan parte de la energía
				call this.products.changeGeometry()
				
			case( ROTATIONAL_REACTOR )
				! La composición es igual antes y después
				this.products = this.reactives
				
				! Los productos utilizan parte de la energía
				call this.products.changeOrientations()
! 				call this.products.changeAngularMomenta()
		end select
		
		if( GOptions_printLevel >= 2 ) then
			write(*,"(A)") "#--------------------------------------------------------------------------------------------------------------"
			write(*,"(A1,3X,6A15,5X,A)") "#", "LnWe", "LnWv", "LnWn", "LnWr", "LnWt", "LnW", "formula"
			write(*,"(A1,3X,6A15,5X,A)") "#", "arb.", "arb.", "arb.", "arb.", "arb.", "arb.", ""
			write(*,"(A)") "#--------------------------------------------------------------------------------------------------------------"
			
			write(*,"(4X,A)") "Reactives"
			sBuffer = this.reactives.weightHistoryLine()
			write(*,"(A)") trim(sBuffer.fstr)
			write(*,"(A)") ""
			write(*,"(4X,A)") "Products"
			
			sBuffer = this.products.weightHistoryLine()
			write(*,"(A)") trim(sBuffer.fstr)

			write(*,"(A)") ""
			write(*,"(A)") "#--------------------------------------------------------------------------------------------------------------"
			write(*,"(A1,3X,5A15,5X,A)") "#", "trans", "intermol", "vib", "rot", "tot", "formula"
			write(*,"(A1,3X,5A15,5X,A)") "#", "eV", "eV", "eV", "eV", ""
			write(*,"(A)") "#--------------------------------------------------------------------------------------------------------------"

			write(*,"(4X,A)") "Reactives"
			
			sBuffer = this.reactives.energyHistoryLine()
			write(*,"(A)") trim(sBuffer.fstr)
			write(*,"(A)") ""
			write(*,"(4X,A)") "Products"
			
			sBuffer = this.products.energyHistoryLine()
			write(*,"(A)") trim(sBuffer.fstr)
			write(*,"(A)") "#--------------------------------------------------------------------------------------------------------------"
			write(*,"(A)") ""
		end if
		
		if( this.products.translationalEnergy() < 0.0_8 ) then
			if( GOptions_printLevel >= 2 ) then
				write(*,*) ""
				write(*,*) "### Warning ### The kinetic energy is negative"
				write(*,"(3X,A,F15.5,A)") "kEnergy = ", this.products.kineticEnergy()/eV, "  eV"
				write(*,*) "products <= reactives"
				write(*,*) ""
			end if
			
			this.products = this.reactives
			this.state = .false.
		end if
		
! 		if( clusDB.isForbidden( this.reactives, this.productos ) ) then
! 			if( GOptions_printLevel >= 2 ) then
! 				write(*,*) ""
! 				write(*,*) "### Warning ### This reaction is forbidden"
! 				write(*,*) "products <= reactives"
! 				write(*,*) ""
! 			end if
! 			
! 			this.products = this.reactives
! 			this.state = .false.
! 		end if
		
	end subroutine run
	
	!>
	!! @brief
	!!
	subroutine execute( this, iParser )
		class(Reactor) :: this
		type(BlocksIFileParser), intent(in) :: iParser
		
		type(FragmentsList) :: reactives
		character(20), allocatable :: reactiveTokens(:)
		
		type(String) :: strReactives
		type(String) :: sBuffer
		real(8) :: rBuffer, p
		integer :: i, iBuffer
		
! 		if( size(FragmentsDB_instance.clusters) >= 6 ) then
! 			write(*,"(A)") "### ERROR ###: Reactor.execute(). The fragments database limit of 6 have been exceeded."
! 			write(*,"(A)") "               Please contact us:"
! 			write(*,"(A)") "                      Dr. Nestor F. Aguirre ( nestor.aguirre@uam.es )"
! 			write(*,"(A)") "                      Dr. Sergio Díaz-Tendero ( sergio.diaztendero@uam.es )"
! 			write(*,"(A)") "                      Prof. M. Paul-Antoine Hervieux ( Paul-Antoine.Hervieux@ipcms.unistra.fr )"
! 			write(*,"(A)") "                      Prof. Fernando Martín ( fernando.martin@uam.es )"
! 			write(*,"(A)") "                      Prof. Manuel Alcamí ( manuel.alcami@uam.es )"
! 			stop
! 		end if
		
		sBuffer = iParser.getString( "REACTOR:reactives" )
		rBuffer = iParser.getReal( "REACTOR:excitationEnergy", def=10.0_8 )*eV
		
		call sBuffer.split( reactiveTokens, ":" )
		
		if( reactiveTokens(1) == "file" ) then
			call reactives.loadXYZ( reactiveTokens(2) )
			reactives.forceInitializing = .false.
			
			call this.initReactor( reactives, rBuffer )
			write(*,"(A40,F15.5,A)") "excitationEnergy = ", rBuffer/eV, " eV"
		else
			strReactives = FragmentsDB_instance.extendFragmentsListName( sBuffer.fstr )
			call strReactives.split( reactiveTokens, "+" )
			
			call reactives.init( size(reactiveTokens) )
			do i=1,size(reactiveTokens)
				iBuffer = FragmentsDB_instance.getIdFromName( reactiveTokens(i) )
				call reactives.set( i, FragmentsDB_instance.clusters(iBuffer) )
			end do
			
			call this.init( reactives, rBuffer )
			write(*,"(A40,F15.5,A)") "excitationEnergy = ", rBuffer/eV, " eV"
		end if
		
! 		call reactives.initialGuessFragmentsList()
		
		sBuffer = iParser.getString( "REACTOR:type", def="V" )
		call this.setType( trim(adjustl(sBuffer.fstr)) )
		write(*,"(A40,A)") "type = ", sBuffer.fstr
		
		call this.run()
		
		sBuffer = iParser.getString( "REACTOR:geomReactivesFile", def="#@NONE@#" )
		if( trim(sBuffer.fstr) /= "#@NONE@#" ) then
			call this.reactives.save( sBuffer.fstr )
		end if
		
		sBuffer = iParser.getString( "REACTOR:geomProductsFile", def="#@NONE@#" )
		if( trim(sBuffer.fstr) /= "#@NONE@#" ) then
			call this.products.save( sBuffer.fstr )
		end if
		
		if( this.state ) then
			! @todo Creo que es bueno dejar este Pi como una función de la clase reactor
			rBuffer = this.products.LnW()-this.reactives.LnW()
			
			write(*,*) ""
			write(*,*) " log(Wn+1/Wn) = ", rBuffer
			write(*,*) ""
			
			if( rBuffer > 0.0_8 ) then
				write(*,"(1X,A)") "@@@ ACCEPTED @@@"
			else
				p = log( RandomUtils_uniform( [0.0_8, 1.0_8] ) )
				
				if( p <= rBuffer ) then
					write(*,"(1X,A)") "@@@ ACCEPTED @@@ p="//trim(adjustl(FString_fromReal(p,"(F10.3)")))&
						//",Pi="//trim(adjustl(FString_fromReal(rBuffer,"(F10.3)")))//"  (p<Pi)"
				else
					write(*,"(1X,A)") "@@@ REJECTED @@@ p="//trim(adjustl(FString_fromReal(p,"(F10.3)")))&
						//",Pi="//trim(adjustl(FString_fromReal(rBuffer,"(F10.3)")))//"  (p>Pi)"
				end if
			end if
		end if
		
		stop
	end subroutine execute
	
	!>
	!! @brief
	!!
	subroutine executeMinFragmentationEnergy( this, iParser )
		class(Reactor) :: this
		type(BlocksIFileParser), intent(in) :: iParser
		
		type(FragmentsList) :: reactives
		character(20), allocatable :: reactiveTokens(:)
		
		type(String) :: strReactives
		type(String) :: sBuffer
		real(8) :: rBuffer, p
		integer :: i, k, iBuffer, dN, id
		logical :: lBuffer
		integer :: nSteps
		logical :: detailed
		
		type(StringHistogram) :: fragmentsHistogram
		class(StringRealMapIterator), pointer :: iter
		type(StringRealPair) :: pair
		real(8) :: energy, minValue, minNegativeValue
		type(String) :: labelMinEnergy, labelMinNegativeEnergy
		logical :: warningNegativeEnergy
		
		if( .not. iParser.isThereBlock( "FRAGMENTS_DATABASE" ) ) then
			return
		end if
		
		nSteps = iParser.getInteger( "FRAGMENTS_DATABASE"//":maxVibNSteps", def=10000 )
		write(*,"(A40,I10)") "maxVibNSteps = ", nSteps
		write(*,*) ""
		
		detailed = iParser.getLogical( "FRAGMENTS_DATABASE"//":maxVibDetailed", def=.false. )
		
		call reactives.init( 1 )
		
		write(*,"(A30,A30,A10)") "reactive", "channel", "energy"
		write(*,"(A30,A30,A10)") "", "", "eV"
		write(*,"(A30,A30,A10)") "--------", "-------", "------"
		
		do id=1,size(FragmentsDB_instance.clusters)
			if( FragmentsDB_instance.clusters(id).nAtoms() == 1 ) cycle
			
			call reactives.set( 1, FragmentsDB_instance.clusters(id) )
			call FragmentsDB_instance.setEnergyReference( FragmentsDB_instance.clusters(id).electronicEnergy )
			
			call this.initReactor( reactives, 100.0_8*eV )
			
			warningNegativeEnergy = .false.
			minValue = Math_INF
			minNegativeValue = Math_INF
			do dN=1,reactives.nAtoms()
				call this.setType( "S:"//trim(FString_fromInteger(dN)) )
				
				call fragmentsHistogram.init()
				
				do k=1,nSteps
					call this.run()
					
					if( this.products.nMolecules() == dN+1 ) then
						call fragmentsHistogram.add( FString_toString( trim(this.products.label()) ) )
					end if
				end do
				
				call fragmentsHistogram.build()
				
				k=1
				call fragmentsHistogram.densityBegin( iter )
				do while( associated(iter) )
					pair = fragmentsHistogram.pair( iter )
					energy = ( FragmentsDB_instance.getEelecFromName(pair.first.fstr)-FragmentsDB_instance.getEelecFromName(reactives.label()) )/eV
					
					if( detailed ) &
						write(*,"(A30,F10.4)") trim(adjustl(pair.first.fstr)), energy
					
					if( energy < minValue ) then
						if( energy < 0.0_8 ) then
							warningNegativeEnergy = .true.
							
							if( energy < minNegativeValue ) then
								minNegativeValue = energy
								labelMinNegativeEnergy = pair.first.fstr
							end if
						else
							minValue = energy
							labelMinEnergy = pair.first.fstr
						end if
					end if
					
					iter => iter.next
					k = k+1
				end do
				
				
				call fragmentsHistogram.clear()
			end do
			
			if( detailed ) &
				write(*,*) ""
				
			if( .not. warningNegativeEnergy ) then
				write(*,"(A30,A30,F10.5)") trim(adjustl(reactives.label())), trim(adjustl(labelMinEnergy.fstr)), minValue
			else
				write(*,"(A30,A30,F10.5,A,A30,F10.5,A)") trim(adjustl(reactives.label())), trim(adjustl(labelMinEnergy.fstr)), minValue, &
										   "   ==> ( ", trim(adjustl(labelMinNegativeEnergy.fstr)), minNegativeValue, " )"
			end if
			
			if( detailed ) then
				write(*,*) ""
				write(*,*) ""
			end if
		end do
		
		write(*,*) ""
		
	end subroutine executeMinFragmentationEnergy

end module Reactor_
