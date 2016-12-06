!>
!! @brief
!!
module FragmentsDB_
	use UnitsConverter_
	use Math_
	use IOStream_
	use String_
	use BlocksIFileParser_
	use RandomUtils_
	use GOptions_
	use AtomicElementsDB_
	use StringIntegerMap_
	
	use Fragment_
	use ModelPotential_
	use GOptionsM3C_
	
	implicit none
	private
	
	public :: &
		FragmentsDB_test
		
	type, public :: FragmentsDB
		type(Fragment), allocatable :: clusters(:)
		type(ModelPotential), allocatable :: potentials(:,:)
		type(ModelPotential) :: atomicPotentials( AtomicElementsDB_nElems, AtomicElementsDB_nElems )
		
		type(String), allocatable :: forbiddenReactions(:)
		logical :: useForbiddenReactionsDetails
		
! 		"C2(d1)+C(t1)-->C3(s1)" --> 5 --> transitionState(5)
! 		|----------str2id_TS----------|
		type(Fragment), allocatable :: transitionState(:)
		type(StringIntegerMap) :: str2id_TS
		logical, allocatable :: involvedInTS(:) ! One for each group. For speed purposes only.
		
		real(8), private :: energyReference_
		logical, private :: useAtomicPotentials
		
		contains
			generic :: init => initDefault, setFragmentsTable
			
			procedure :: initDefault
			procedure :: fromInputFile
			procedure, private :: setFragmentsTable
			procedure, private :: setForbiddenReactionsTable
			procedure, private :: setTransitionStatesTable
			procedure, private :: setPotentialTable
			procedure, private :: setAtomicPotentialTable
			procedure, private :: checkPotentials
			final :: destroyFragmentsDB
			procedure :: nMolecules
			procedure :: getIdFromName
			procedure :: getEelecFromName
			procedure :: extendFragmentsListName
			procedure :: potential
			procedure :: isThereAModel
			procedure :: isForbidden
			procedure :: energyReference
			procedure :: setEnergyReference
	end type FragmentsDB
	
	type(FragmentsDB), public :: FragmentsDB_instance
	
	contains
	
	!>
	!! @brief Constructor
	!!
	subroutine initDefault( this )
		class(FragmentsDB) :: this 
		
	end subroutine initDefault
	
	!>
	!! @brief
	!!
	subroutine fromInputFile( this, iParser )
		class(FragmentsDB) :: this
		type(BlocksIFileParser), intent(in) :: iParser
		
		type(String), allocatable :: fragmentsTable(:)
		type(String), allocatable :: forbiddenReactionsTable(:)
		type(String), allocatable :: potentialTable(:)
		
		!--------------------------------------------------------------------
		! Loading fragments database
		!--------------------------------------------------------------------
		if( iParser.isThereBlock( "FRAGMENTS_DATABASE" ) ) then
			call iParser.getBlock( "FRAGMENTS_DATABASE", fragmentsTable )
			
			call this.setFragmentsTable( &
				fragmentsTable, &
				iParser.getString( "FRAGMENTS_DATABASE:reference", def="@@NONE@@" ), &
				store=iParser.getString( "FRAGMENTS_DATABASE:store", def="." ) &
			)
			
			deallocate(fragmentsTable)
		else
			write(*,*) "### ERROR ### M3C: FRAGMENTS_DATABASE block is required"
			stop
		end if
		
		!--------------------------------------------------------------------
		! Loading the forbidden reactions
		!--------------------------------------------------------------------
		if( iParser.isThereBlock( "FORBIDDEN_REACTIONS" ) ) then
			call iParser.getBlock( "FORBIDDEN_REACTIONS", forbiddenReactionsTable )
			
			call this.setForbiddenReactionsTable( &
				forbiddenReactionsTable, &
				details=iParser.getLogical( "FORBIDDEN_REACTIONS:details", def=.false. ) &
			)
			
			deallocate(forbiddenReactionsTable)
		end if
		
		!--------------------------------------------------------------------
		! Loading transition states
		!--------------------------------------------------------------------
		if( iParser.isThereBlock( "TRANSITION_STATES_DATABASE" ) ) then
			call iParser.getBlock( "TRANSITION_STATES_DATABASE", fragmentsTable )
			
			call this.setTransitionStatesTable( &
				fragmentsTable, &
				store=iParser.getString( "TRANSITION_STATES_DATABASE:store", def="." ) &
			)
			
			deallocate(fragmentsTable)
		end if

		!--------------------------------------------------------------------
		! Loading potentials
		!--------------------------------------------------------------------
		this.useAtomicPotentials = .false.
		
		if( iParser.isThereBlock( "POTENTIAL_TABLE" ) ) then
			call iParser.getBlock( "POTENTIAL_TABLE", potentialTable )
			
			call this.setPotentialTable( potentialTable, iParser.getLogical( "POTENTIAL_TABLE:coulombContribution", def=.true. ) )
			deallocate(potentialTable)
			
			if( iParser.getLogical( "POTENTIAL_TABLE:check", def=.false. ) ) then
				call this.checkPotentials( &
					iParser.getReal( "POTENTIAL_TABLE:check.rMin", def=0.0_8 )*angs, &
					iParser.getReal( "POTENTIAL_TABLE:check.rMax", def=10.0_8 )*angs, &
					iParser.getReal( "POTENTIAL_TABLE:check.rStep", def=0.1_8 )*angs, &
					iParser.getString( "POTENTIAL_TABLE:check.outputFile", def="#@NONE@#" ) &
				)
			end if
			
			return ! Hay preferencia de los potentiales moleculares que atomicos
		end if
		
		if( iParser.isThereBlock( "ATOMIC_POTENTIAL_TABLE" ) ) then
			call iParser.getBlock( "ATOMIC_POTENTIAL_TABLE", potentialTable )
			
			call this.setAtomicPotentialTable( potentialTable, iParser.getLogical( "ATOMIC_POTENTIAL_TABLE:coulombContribution", def=.true. ) )
			deallocate(potentialTable)
			
			if( iParser.getLogical( "ATOMIC_POTENTIAL_TABLE:check", def=.false. ) ) then
				call this.checkPotentials( &
					iParser.getReal( "ATOMIC_POTENTIAL_TABLE:check.rMin", def=0.0_8 )*angs, &
					iParser.getReal( "ATOMIC_POTENTIAL_TABLE:check.rMax", def=10.0_8 )*angs, &
					iParser.getReal( "ATOMIC_POTENTIAL_TABLE:check.rStep", def=0.1_8 )*angs, &
					iParser.getString( "ATOMIC_POTENTIAL_TABLE:check.outputFile", def="#@NONE@#" ) &
				)
			end if
			
			return ! Hay preferencia de los potentiales moleculares que atomicos
		end if
		
	end subroutine fromInputFile
	
	!>
	!! @brief
	!!
	subroutine setFragmentsTable( this, fragmentsTable, strReference, store )
		class(FragmentsDB) :: this
		type(String), allocatable, intent(in) :: fragmentsTable(:)
		type(String), intent(in) :: strReference
		type(String), optional, intent(in) :: store
		
		integer :: i, n
		character(100), allocatable :: tokens(:)
		real(8) :: rBuffer
		
		if( allocated(this.clusters) ) deallocate( this.clusters )
		allocate( this.clusters(size(fragmentsTable)) )
		
		call GOptions_section( "FRAGMENTS DATABASE INITIALIZATION", indent=1 )
		
		do i=1,size(fragmentsTable)
			call this.clusters(i).fromMassTableRow( fragmentsTable(i).fstr, id=i, store=store.fstr )
			
			call FString_split( fragmentsTable(i).fstr, tokens, " " )
			
			!------------------------------------------
			! Choosing the maximum vibrational energy
			!------------------------------------------
			if( size(tokens) >= 8 .and. this.clusters(i).nAtoms() /= 1 ) then
				if( FString_isNumeric( tokens(8) ) ) then
					this.clusters(i).maxEvib = FString_toReal( tokens(8) )*eV
				else
					this.clusters(i).maxEvib = this.getEelecFromName( tokens(8) ) &
							 - this.clusters(i).electronicEnergy
				end if
			else
				this.clusters(i).maxEvib = 0.0_8
			end if
			
			write(IO_STDOUT,"(4X,A22,F15.7,A)")  "           maxEvib = ", this.clusters(i).maxEvib/eV, "   eV"
			
			deallocate( tokens )
		end do
		
		if( strReference /= "@@NONE@@" ) then
			this.energyReference_ = this.getEelecFromName( strReference.fstr )
		else
			this.energyReference_ = 0.0_8
			do i=1,size(fragmentsTable)
				if( this.clusters(i).electronicEnergy < this.energyReference_ ) then
					this.energyReference_ = this.clusters(i).electronicEnergy
				end if
			end do
		end if
		
		write(IO_STDOUT,*) ""
		write(IO_STDOUT,*) ""
		write(IO_STDOUT,"(A,F15.7,A)")  "     Electronic energy reference = ", this.energyReference_/eV, "   eV"
		write(IO_STDOUT,*) ""
		
		call GOptions_section( "END FRAGMENTS DATABASE INITIALIZATION", indent=1 )
		
	end subroutine setFragmentsTable
	
	!>
	!! @brief
	!!
	subroutine setForbiddenReactionsTable( this, forbiddenReactionsTable, details )
		class(FragmentsDB) :: this
		type(String), allocatable, intent(in) :: forbiddenReactionsTable(:)
		logical, optional, intent(in) :: details
		
		integer :: i, j, k, n
		logical :: duplicated
		character(100), allocatable :: tokens(:)
		real(8) :: rBuffer
		
		this.useForbiddenReactionsDetails = .false.
		if( present(details) ) this.useForbiddenReactionsDetails = details
		
		call GOptions_section( "FORBIDDEN REACTIONS INITIALIZATION", indent=1 )
		
		write(IO_STDOUT,"(4X,A22,L15)")  "           details = ", details
		
		if( allocated(this.forbiddenReactions) ) deallocate( this.forbiddenReactions )
		
		call forbiddenReactionsTable(1).split( tokens, "<-->" )
		if( size(tokens) == 2 ) then
			allocate( this.forbiddenReactions(2*size(forbiddenReactionsTable)) )
			write(IO_STDOUT,"(4X,A22,L15)")  "           directed = ", .false.
		else if( size(tokens) == 1 ) then
			allocate( this.forbiddenReactions(size(forbiddenReactionsTable)) )
			write(IO_STDOUT,"(4X,A22,L15)")  "           directed = ", .true.
		else
			call GOptions_error( &
					"Wrong format in FORBIDDEN_REACTIONS block", &
					"SMoleculeDB.setForbiddenReactionsTable()", &
					trim(forbiddenReactionsTable(1).fstr) &
				)
		end if
		
		write(IO_STDOUT,*) ""
		
		j=1
		do i=1,size(forbiddenReactionsTable)
			call forbiddenReactionsTable(i).split( tokens, "<-->" )
			
			duplicated = .false.
			do k=1,j-1
				if( this.forbiddenReactions(k) == trim(tokens(1))//"-->"//trim(tokens(2)) ) then
					duplicated = .true.
				end if
			end do
			
			if( .not. duplicated ) then
				if ( size(tokens) == 2 ) then
					
					this.forbiddenReactions(j) = trim(tokens(1))//"-->"//trim(tokens(2))
					write(IO_STDOUT,"(4X,I6,3X,A)")  j, trim(this.forbiddenReactions(j).fstr)
					
					! Esto evita que se almacenen doblemente las isomerizaciones
					if( trim(tokens(1)) /= trim(tokens(2)) ) then
						this.forbiddenReactions(j+1) = trim(tokens(2))//"-->"//trim(tokens(1))
						write(IO_STDOUT,"(4X,I6,5X,A)")  j+1, trim(this.forbiddenReactions(j+1).fstr)
						
						j = j + 2
					else
						j = j + 1
					end if
					
				else if( .not. duplicated .and. size(tokens) == 1 ) then
				
					this.forbiddenReactions(j) = trim(adjustl(forbiddenReactionsTable(i).fstr))
					write(IO_STDOUT,"(4X,I6,3X,A)")  j, trim(this.forbiddenReactions(j).fstr)
					
					j = j + 1
				else
					call GOptions_error( &
							"Wrong format in FORBIDDEN_REACTIONS block", &
							"SMoleculeDB.setForbiddenReactionsTable()", &
							trim(forbiddenReactionsTable(i).fstr) &
						)
				end if
			else
				call GOptions_error( &
						"Duplicated reaction in FORBIDDEN_REACTIONS block", &
						"SMoleculeDB.setForbiddenReactionsTable()", &
						trim(forbiddenReactionsTable(i).fstr) &
					)
			end if
		end do
		
		call GOptions_section( "END FORBIDDEN REACTIONS INITIALIZATION", indent=1 )
		
		if( allocated(tokens) ) deallocate( tokens )
		
	end subroutine setForbiddenReactionsTable
	
	!>
	!! @brief
	!!
	subroutine setTransitionStatesTable( this, transitionStatesTable, store )
		class(FragmentsDB) :: this
		type(String), allocatable, intent(in) :: transitionStatesTable(:)
		type(String), optional, intent(in) :: store
		
		integer :: i, j, k, n
		character(100), allocatable :: tokens(:), tokens2(:), tokens3(:)
		real(8) :: rBuffer
		
		if( allocated(this.transitionState) ) deallocate( this.transitionState )
		allocate( this.transitionState(size(transitionStatesTable)) )
		
		call this.str2id_TS.init()
		
		if( allocated(this.involvedInTS) ) deallocate( this.involvedInTS )
		allocate( this.involvedInTS(size(this.clusters)) )
		
		call GOptions_section( "TRANSITION_STATES DATABASE INITIALIZATION", indent=1 )
		
		do i=1,size(transitionStatesTable)
			call this.transitionState(i).fromMassTableRow( transitionStatesTable(i).fstr, id=i, store=store.fstr )
			
			call FString_split( transitionStatesTable(i).fstr, tokens, " " )
			
			!------------------------------------------
			! Choosing the maximum vibrational energy
			!------------------------------------------
			if( size(tokens) >= 8 .and. this.transitionState(i).nAtoms() /= 1 ) then
				if( FString_isNumeric( tokens(8) ) ) then
					this.transitionState(i).maxEvib = FString_toReal( tokens(8) )*eV
				else
					this.transitionState(i).maxEvib = this.getEelecFromName( tokens(8) ) &
							 - this.transitionState(i).electronicEnergy
				end if
			else
				this.transitionState(i).maxEvib = 0.0_8
			end if
			
			write(IO_STDOUT,"(4X,A22,F15.7,A)")  "           maxEvib = ", this.transitionState(i).maxEvib/eV, "   eV"
			
			!------------------------------------------
			! Choosing the reactives and products
			!------------------------------------------
			if( size(tokens) >= 9 .and. this.transitionState(i).nAtoms() /= 1 ) then
				call FString_split( trim(adjustl(tokens(9))), tokens2, "-->" )
				
				write(*,*) "Reactivos"
				call FString_split( trim(adjustl(tokens2(1))), tokens3, "+" )
				do j=1,size(tokens3)
					write(*,*) "token = ", trim(tokens3(j))
					do k=1,size(this.clusters)
						if( trim(tokens3(j)) == this.clusters(k).label() ) then
							write(IO_STDOUT,*) k, this.clusters(k).id
							this.involvedInTS(k) = .true.
						end if
					end do
				end do
				deallocate( tokens3 )
				
				write(*,*) "Productos"
				call FString_split( trim(adjustl(tokens2(2))), tokens3, "+" )
				do j=1,size(tokens3)
					write(*,*) "token = ", trim(tokens3(j))
					do k=1,size(this.clusters)
						if( trim(tokens3(j)) == this.clusters(k).label() ) then
							write(IO_STDOUT,*) k, this.clusters(k).id
							this.involvedInTS(k) = .true.
						end if
					end do
				end do
				deallocate( tokens3 )
				
				deallocate( tokens2 )
				
				write(*,*) "Agregado: ", FString_toString(trim(adjustl(tokens(9)))//"-->"//trim(adjustl(tokens(10)))), ">>", i
				call this.str2id_TS.insert( FString_toString(trim(adjustl(tokens(9)))//"-->"//trim(adjustl(tokens(10)))), i )
			else
				call GOptions_error( &
						"Bad number of atoms in transition state (N=0)", &
						"SMoleculeDB.setTransitionStatesTable()", &
						trim(transitionStatesTable(i).fstr) &
					)
			end if
			
			deallocate( tokens )
		end do
		
		write(IO_STDOUT,*) ""
		
		call GOptions_section( "END TRANSITION_STATES DATABASE INITIALIZATION", indent=1 )
		
	end subroutine setTransitionStatesTable
	
	!>
	!! @brief
	!!
	subroutine checkPotentials( this, rMin, rMax, rStep, oFileName )
		class(FragmentsDB) :: this
		real(8), intent(in) :: rMin, rMax, rStep
		type(String), intent(in) :: oFileName
		
		type(OFStream) :: oFile
		integer :: i, j
		real(8):: r
		integer :: maxMass
		
		
		write(*,*) ""
		write(*,*) "-------------------------------"
		write(*,*) " CHECKING POTENTIAL"
		write(*,*) ""
		write(*,"(A30,F10.5,A)") " rMin = ", rMin/angs, " A"
		write(*,"(A30,F10.5,A)") " rMax = ", rMax/angs, " A"
		write(*,"(A30,F10.5,A)") "rStep = ", rStep/angs, " A"
		if( trim(oFileName.fstr) /= "#@NONE@#" ) then
			write(*,"(A30,A)") "outputFile = ", oFileName.fstr
			call oFile.init( oFileName.fstr )
		end if
		write(*,"(A)") ""
		write(*,"(A)") ""
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Escribe la identidad de las moleculas
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		write(oFile.unit,"(A)") "# Molecules identity"
		write(oFile.unit,"(A1,A10,A20,5X,A)") "#", "covR", "Ee-E0", "label"
		write(oFile.unit,"(A1,A10,A20,5X,A)") "#", "----", "-----", "-----"
		maxMass = 0
		do i=1,this.nMolecules()
			if( this.clusters(i).nAtoms() > 0 ) then
				write(oFile.unit,"(1X,F10.2,F20.8,5X,A)") this.clusters(i).radius( type=GOptionsM3C_radiusType )/angs, &
					( this.clusters(i).electronicEnergy - this.energyReference() )/eV, &
					trim(this.clusters(i).label())
					
			end if
			
			if( this.clusters(i).massNumber() > maxMass ) then
				maxMass = this.clusters(i).massNumber()
			end if
		end do
		
		write(oFile.unit,"(A)") ""
		write(oFile.unit,"(A)") ""
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Escribe la identidad de los limites disociativos
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		write(oFile.unit,"(A)") "# Dissociative limits"
		write(oFile.unit,"(A1,A10,A20,5X,A)") "#", "Rmax", "V", "label"
		write(oFile.unit,"(A1,A10,A20,5X,A)") "#", "----", "-----", "-----"
		do i=1,this.nMolecules()
			do j=i,this.nMolecules()
				if( this.clusters(i).massNumber()+this.clusters(j).massNumber() < maxMass ) then
					if( this.potentials(i,j).getId() /= 0 ) then
						write(oFile.unit,"(F10.2,F20.8,5X,A)") rMax/angs, &
							( this.clusters(i).electronicEnergy &
							+ this.clusters(j).electronicEnergy &
							+ this.potential( i, j, rMax ) &
							- this.energyReference() &
							)/eV, &
							trim(this.clusters(i).label())//"+"//trim(this.clusters(j).label())
					end if
				end if
			end do
		end do
		
		write(oFile.unit,"(A)") ""
		write(oFile.unit,"(A)") ""
		
		write(oFile.unit,"(A)") "# Dissociative limits"
		write(oFile.unit,"(A1,2A20,5X,A)") "#", "r", "V(r)", "label"
		write(oFile.unit,"(A1,2A20,5X,A)") "#", "---", "-----", "-----"
		do i=1,this.nMolecules()
			do j=i,this.nMolecules()
! 					do j=i+1,this.nMolecules()
				if( this.potentials(i,j).getId() > 1 ) then  ! Esto evita que imprima NONE y HARDSPHERE
				
					write(*,"(A)") "#  " &
						//this.clusters(i).name &
						//"+" &
						//this.clusters(j).name
					
					if( trim(oFileName.fstr) /= "#@NONE@#" ) then
						write(oFile.unit,"(A)") "#  " &
							//this.clusters(i).name &
							//"+" &
							//this.clusters(j).name
					end if

					r = rMin
					do while( r <= rMax )
						write(*,"(2F20.7)") r/angs, ( &
							+ this.clusters(i).electronicEnergy &
							+ this.clusters(j).electronicEnergy &
							+ this.potential( i, j, r ) &
							- this.energyReference() &
							)/eV
						
						if( trim(oFileName.fstr) /= "#@NONE@#" ) then
							write(oFile.unit,"(2F20.7)") r/angs, ( &
								+ this.clusters(i).electronicEnergy &
								+ this.clusters(j).electronicEnergy &
								+ this.potential( i, j, r ) &
								- this.energyReference() &
								)/eV
						end if
						
						r = r + rStep
					end do
					write(*,"(A)") ""
					write(*,"(A)") ""
					
					if( trim(oFileName.fstr) /= "#@NONE@#" ) then
						write(oFile.unit,"(A)") ""
						write(oFile.unit,"(A)") ""
					end if
					
				end if
			end do
		end do
		
! 				do i=1,this.nMolecules()
! 					do j=1,this.nMolecules()
! 						do k=1,this.nMolecules()
! 							if( this.potentials(i,j,k).getId() /= 0 ) then
! 							
! 								write(*,"(A)") "#  " &
! 									//this.clusters(i).name &
! 									//" + " &
! 									//this.clusters(j).name &
! 									//" -> " &
! 									//this.clusters(k).name
! 								
! 								if( trim(oFileName.fstr) /= "#@NONE@#" ) then
! 									write(oFile.unit,"(A)") "#  " &
! 										//this.clusters(i).name &
! 										//" + " &
! 										//this.clusters(j).name &
! 										//" -> " &
! 										//this.clusters(k).name
! 								end if
! 
! 								r = rMin
! 								do while( r <= rMax )
! 									write(*,"(2F20.7)") r/angs, this.potential( i, j, k, r )/eV
! 									
! 									if( trim(oFileName.fstr) /= "#@NONE@#" ) then
! 										write(oFile.unit,"(2F20.7)") r/angs, this.potential( i, j, k, r )/eV
! 									end if
! 									
! 									r = r + rStep
! 								end do
! 								write(*,"(A)") ""
! 								write(*,"(A)") ""
! 								
! 								if( trim(oFileName.fstr) /= "#@NONE@#" ) then
! 									write(oFile.unit,"(A)") ""
! 									write(oFile.unit,"(A)") ""
! 								end if
! 								
! 							end if
! 						end do
! 					end do
! 				end do
		
		if( trim(oFileName.fstr) /= "#@NONE@#" ) then
			call oFile.close()
		end if
		
		stop
	end subroutine checkPotentials
	
	!>
	!! @brief
	!!
	subroutine setPotentialTable( this, potentialTable, coulombContribution )
		class(FragmentsDB) :: this
		type(String), allocatable, intent(in) :: potentialTable(:)
		logical, intent(in) :: coulombContribution
		
		integer :: i, j, n, idR1, idR2
		character(100), allocatable :: cols(:)
		character(100), allocatable :: tokens(:)
		logical :: firstTime
		integer :: maxMass
		
		if( allocated(this.potentials) ) deallocate( this.potentials )
		allocate( this.potentials(this.nMolecules(),this.nMolecules()) )
		
		write(*,"(A)") ""
		write(*,"(<GOptions_indentLength*1>X,A)") "-----------------------------------------------------------------------"
		write(*,"(<GOptions_indentLength*1>X,A)") " MODEL POTENTIALS"
		write(*,"(<GOptions_indentLength*1>X,A)") "-----------------------------------------------------------------------"
		write(*,"(A)") ""
		write(*,"(<GOptions_indentLength*1>X,A10,A10,11X,A15,A10)") "R1", "R2", "potential", "params"
		write(*,"(<GOptions_indentLength*1>X,A10,A10,11X,A15,A10)") "----", "----", "---------", "------"
		
		firstTime = .true.
		do i=1,size(potentialTable)
		
			if( firstTime ) then
				write(*,*) ""
				write(*,"(A)") "USER"
				firstTime = .false.
			end if
		
			call FString_split( potentialTable(i).fstr, cols, " " )
			
			if( size(cols) >= 2 ) then
				call FString_split( cols(1), tokens, "+" )
				
				idR1 = this.getIdFromName( tokens(1) )
				idR2 = this.getIdFromName( tokens(2) )
				
				write(*,"(<GOptions_indentLength*1>X,A10,A10,5X,2I3)", advance="no") &
						this.clusters(idR1).name, this.clusters(idR2).name, idR1, idR2
						
				if( this.potentials(idR1,idR2).getId() == 0 ) then
					call this.potentials(idR1,idR2).init( cols(2) )
					call this.potentials(idR2,idR1).init( cols(2) )
					
					write(*,"(A15,10F10.3)") &
						trim(MODEL_NAME(this.potentials(idR1,idR2).getId())), this.potentials(idR1,idR2).potParams
					exit
				end if
			else
				call GOptions_error( &
						"Intermolecular potential table row incomplete, this row should have had 3 columns", &
						"SMoleculeDB.setPotentialTable()", &
						trim(potentialTable(i).fstr) &
					)
			end if
			
		end do
		
		maxMass = 0
		do i=1,this.nMolecules()
			if( this.clusters(i).massNumber() > maxMass ) then
				maxMass = this.clusters(i).massNumber()
			end if
		end do
		
		if( coulombContribution ) then
			firstTime = .true.
			do i=1,this.nMolecules()
				do j=i,this.nMolecules()
					if( this.potentials(i,j).getId() == 0 .and. this.clusters(i).charge > 0 .and. this.clusters(j).charge > 0 &
						.and. this.clusters(i).massNumber() + this.clusters(j).massNumber() <= maxMass ) then
						
						if( firstTime ) then
							write(*,*) ""
							write(*,"(A)") "COULOMB"
							firstTime = .false.
						end if
						
						write(*,"(<GOptions_indentLength*1>X,A10,A10,5X,2I3)", advance="no") &
								this.clusters(i).name, this.clusters(j).name, i, j
								
						call this.potentials(i,j).init( "COULOMB(1.0)" )
						call this.potentials(j,i).init( "COULOMB(1.0)" )
						
						write(*,"(A15,10F10.3)") &
							trim(MODEL_NAME(this.potentials(i,j).getId())), this.potentials(i,j).potParams
					end if
				end do
			end do
		end if
		
		firstTime = .true.
		do i=1,this.nMolecules()
			do j=i,this.nMolecules()
				if( this.potentials(i,j).getId() == 0 &
					.and. this.clusters(i).massNumber() + this.clusters(j).massNumber() <= maxMass ) then
					
					if( firstTime ) then
						write(*,*) ""
						write(*,"(A)") "DEFAULT"
						firstTime = .false.
					end if
				
					write(*,"(<GOptions_indentLength*1>X,A10,A10,5X,2I3)", advance="no") &
							this.clusters(i).name, this.clusters(j).name, i, j
							
					call this.potentials(i,j).init( "HARDSPHERE()" )
					call this.potentials(j,i).init( "HARDSPHERE()" )
					
					write(*,"(A15,10F10.3)") &
						trim(MODEL_NAME(this.potentials(i,j).getId())), this.potentials(i,j).potParams
				end if
			end do
		end do
		
		write(*,*) ""
		do i=1,this.nMolecules()
			write(*,*) i, "--"//trim(adjustl(this.clusters(i).name))//"--"
		end do
		write(*,*) ""
		
		write(*,"(<GOptions_indentLength*1>X,A)") "-----------------------------------------------------------------------"
		
		if( allocated(cols) ) deallocate( cols )
		if( allocated(tokens) ) deallocate( tokens )
	end subroutine setPotentialTable
	
	!>
	!! @brief
	!!
	subroutine setAtomicPotentialTable( this, potentialTable, coulombContribution )
		class(FragmentsDB) :: this
		type(String), allocatable, intent(in) :: potentialTable(:)
		logical, intent(in) :: coulombContribution
		
		integer :: i, j, n, idR1, idR2
		character(100), allocatable :: cols(:)
		character(100), allocatable :: tokens(:)
		logical :: firstTime
		integer :: maxMass
		integer :: composition( AtomicElementsDB_nElems )
		
		write(*,"(A)") ""
		write(*,"(<GOptions_indentLength*1>X,A)") "-----------------------------------------------------------------------"
		write(*,"(<GOptions_indentLength*1>X,A)") " ATOMIC MODEL POTENTIALS"
		write(*,"(<GOptions_indentLength*1>X,A)") "-----------------------------------------------------------------------"
		write(*,"(A)") ""
		write(*,"(<GOptions_indentLength*1>X,A10,A10,11X,A15,A10)") "R1", "R2", "potential", "params"
		write(*,"(<GOptions_indentLength*1>X,A10,A10,11X,A15,A10)") "----", "----", "---------", "------"
		
		firstTime = .true.
		do i=1,size(potentialTable)
		
			if( firstTime ) then
				write(*,*) ""
				write(*,"(A)") "USER"
				firstTime = .false.
			end if
		
			call FString_split( potentialTable(i).fstr, cols, " " )
			
			if( size(cols) >= 2 ) then
				call FString_split( cols(1), tokens, "+" )
				
				idR1 = AtomicElementsDB_instance.atomicNumber( tokens(1) )
				idR2 = AtomicElementsDB_instance.atomicNumber( tokens(2) )
				
				write(*,"(<GOptions_indentLength*1>X,A10,A10,5X,2I3)", advance="no") &
						trim(tokens(1)), trim(tokens(2)), idR1, idR2
						
				if( this.atomicPotentials(idR1,idR2).getId() == 0 ) then
					call this.atomicPotentials(idR1,idR2).init( cols(2) )
					call this.atomicPotentials(idR2,idR1).init( cols(2) )
					
					write(*,"(A15,10F10.3)") &
						trim(MODEL_NAME(this.atomicPotentials(idR1,idR2).getId())), this.atomicPotentials(idR1,idR2).potParams
					exit
				end if
			else
				call GOptions_error( &
						"Intermolecular potential table row incomplete, this row should have had 3 columns", &
						"SMoleculeDB.setPotentialTable()", &
						trim(potentialTable(i).fstr) &
					)
			end if
			
		end do
		
! 		maxMass = 0
! 		do i=1,this.nMolecules()
! 			if( this.clusters(i).massNumber() > maxMass ) then
! 				maxMass = this.clusters(i).massNumber()
! 			end if
! 		end do
		
		composition = 0
		do i=1,this.nMolecules()
			composition = composition + this.clusters(i).composition
		end do
		composition = merge( 1, 0, composition>0 )
		
		firstTime = .true.
		do i=1,size(composition)
			do j=i,size(composition)
				if( composition(i) /= 0 .and. composition(j) /= 0 .and. this.atomicPotentials(i,j).getId() == 0 ) then
! 					.and. this.clusters(i).massNumber() + this.clusters(j).massNumber() < maxMass ) then
					
					if( firstTime ) then
						write(*,*) ""
						write(*,"(A)") "DEFAULT"
						firstTime = .false.
					end if
				
					write(*,"(<GOptions_indentLength*1>X,A10,A10,5X,2I3)", advance="no") &
							AtomicElementsDB_instance.symbol(i), &
							AtomicElementsDB_instance.symbol(j), i, j
							
					call this.atomicPotentials(i,j).init( "HARDSPHERE()" )
					call this.atomicPotentials(j,i).init( "HARDSPHERE()" )
					
					write(*,"(A15,10F10.3)") &
						trim(MODEL_NAME(this.atomicPotentials(i,j).getId())), this.atomicPotentials(i,j).potParams
				end if
			end do
		end do
		
		write(*,*) ""
		do i=1,size(composition)
			if( composition(i) /= 0 ) then
				write(*,*) i, "--"//trim(AtomicElementsDB_instance.symbol(i))//"--"
			end if
		end do
		write(*,*) ""
		
		write(*,"(<GOptions_indentLength*1>X,A)") "-----------------------------------------------------------------------"
		
		if( allocated(cols) ) deallocate( cols )
		if( allocated(tokens) ) deallocate( tokens )
	end subroutine setAtomicPotentialTable
	
	!>
	!! @brief Destructor
	!!
	subroutine destroyFragmentsDB( this )
		type(FragmentsDB) :: this
		
		if( allocated(this.clusters) ) deallocate(this.clusters)
		if( allocated(this.potentials) ) deallocate(this.potentials)
		if( allocated(this.forbiddenReactions) ) deallocate(this.forbiddenReactions)
	end subroutine destroyFragmentsDB
	
	!>
	!! @brief Returns the number of clusters into the database
	!!
	pure function nMolecules( this ) result( output )
		class(FragmentsDB), intent(in) :: this
		integer :: output
		
		output = size( this.clusters )
	end function nMolecules
	
	!>
	!! @brief
	!!
	function getIdFromName( this, name ) result( output )
		class(FragmentsDB), intent(in) :: this
		character(*), intent(in) :: name
		integer :: output
		
		integer :: i
		
		output = -1
		do i=1,size( this.clusters )
			if( trim(this.clusters(i).name) == trim(adjustl(name)) ) then
				output = i
				return
			end if
		end do
		
		call GOptions_error( &
				"The cluster --"//trim(adjustl(name))//"-- doesn't exist in the database", &
				"FragmentsDB.getIdFromName()" &
			)
	end function getIdFromName
	
	!>
	!! @brief
	!!
	function getEelecFromName( this, name ) result( output )
		class(FragmentsDB), intent(in) :: this
		character(*), intent(in) :: name
		real(8) :: output
		
		integer :: i, id, factor
		character(20), allocatable :: reactiveTokens(:)
		character(20), allocatable :: subTokens(:)
		
		call FString_split( name, reactiveTokens, "+" )
		
		output = 0.0_8
		do i=1,size(reactiveTokens)
			call FString_split( reactiveTokens(i), subTokens, "*" )
			
			
			if( size(subTokens) < 2 ) then
				id = this.getIdFromName( subTokens(1) )
				factor = 1
			else
				id = this.getIdFromName( subTokens(2) )
				factor = FString_toInteger( subTokens(1) )
			end if
			
			output = output + factor*this.clusters(id).electronicEnergy
			
			deallocate( subTokens )
		end do
		
		if( allocated(reactiveTokens) ) deallocate( reactiveTokens )
	end function getEelecFromName
	
	!>
	!! @brief Expande las formulas que utilizan algunos operadores como *
	!!
	!! Ejemplo:
	!!       input = sHe+10*C10s+2B
	!!      output = sHe+C10s+C10s+C10s+C10s+C10s+C10s+C10s+C10s+C10s+C10s+2B
	!!
	!!       input = sHe+3*C10s+2*2B
	!!      output = sHe+C10s+C10s+C10s+2B+2B
	!!
	function extendFragmentsListName( this, name ) result( output )
		class(FragmentsDB), intent(in) :: this
		character(*), intent(in) :: name
		type(String) :: output
		
		integer :: i, j, factor
		character(20), allocatable :: tokens(:)
		character(20), allocatable :: subTokens(:)
		
		call FString_split( name, tokens, "+" )
		
		output = ""
		do i=1,size(tokens)
			call FString_split( tokens(i), subTokens, "*" )
			
			if( size(subTokens) < 2 ) then
				if( i /= 1 ) then
					output = output + "+" + trim(subTokens(1))
				else
					output = output + trim(subTokens(1))
				end if
			else
				factor = FString_toInteger( subTokens(1) )
				
				do j=1,factor
					if( i /= 1 .or. j /= 1 ) then
						output = output + "+" + trim(subTokens(2))
					else
						output = output + trim(subTokens(2))
					end if
				end do
			end if
			
			deallocate( subTokens )
		end do
		
		deallocate( tokens )
	end function extendFragmentsListName
		
	!>
	!! @brief Returns the intermolecular potential
	!!
	function potential( this, idR1, idR2, r ) result ( Vij )
		class(FragmentsDB), intent(in) :: this
		integer, intent(in) :: idR1
		integer, intent(in) :: idR2
		real(8), intent(in) :: r
		real(8) :: Vij
		
		real(8), allocatable :: params(:)
		integer :: i, j, atomId1, atomId2
		
		if( .not. this.useAtomicPotentials ) then
			if( .not. allocated(this.potentials) ) then
				Vij = 0.0_8
				return
			end if
			
			if( idR1 > size(this.potentials,dim=1) .or. idR2 > size(this.potentials,dim=1) ) then
				call GOptions_error( &
					"Unknown reactive(s) ("//trim(FString_fromInteger(idR1))//","//trim(FString_fromInteger(idR1))//")", &
					"FragmentsDB.potential()" &
				)
			end if
			
			if( this.potentials(idR1,idR2).getId() == 0 ) then
				Vij = 0.0_8
				return
			end if
			
			allocate( params(3) )
			
			params(1) = this.clusters(idR1).radius( type=GOptionsM3C_radiusType ) &
					+ this.clusters(idR2).radius( type=GOptionsM3C_radiusType ) ! Ra+Rb = Re
			params(2) = real(this.clusters(idR1).charge,8) ! qa = q1
			params(3) = real(this.clusters(idR2).charge,8) ! qb = q2
			
			Vij = this.potentials(idR1,idR2).evaluate( r, params )
			
			deallocate( params )
		else
			allocate( params(3) )
			
			Vij = 0.0_8
			do i=1,this.clusters(idR1).nAtoms()
				do j=1,this.clusters(idR2).nAtoms()
					
					atomId1 = AtomicElementsDB_instance.atomicNumber( trim(this.clusters(idR1).atoms(i).symbol) )
					atomId2 = AtomicElementsDB_instance.atomicNumber( trim(this.clusters(idR2).atoms(j).symbol) )
					
					if( this.atomicPotentials(atomId1,atomId2).getId() == 0 ) then
						Vij = 0.0_8
						exit
					end if
					
					params(1) = this.clusters(idR1).atoms(i).radius( type=GOptionsM3C_radiusType ) &
							+ this.clusters(idR2).atoms(j).radius( type=GOptionsM3C_radiusType ) ! Ra+Rb = Re
! 					params(2) = real(this.clusters(atomId1).charge,8) ! qa = q1
! 					params(3) = real(this.clusters(atomId2).charge,8) ! qb = q2
					params(2) = 0.0_8 ! qa = q1
					params(3) = 0.0_8 ! qb = q2
					
					Vij = Vij + this.atomicPotentials(atomId1,atomId2).evaluate( r, params )
				end do
			end do
			
			deallocate( params )
		end if
	end function potential
	
	!>
	!! @brief Returns true if the there is a potential model associated to the clusters
	!!        which ids are idR1 and idR2
	!!
	function isThereAModel( this, idR1, idR2 ) result ( output )
		class(FragmentsDB), intent(in) :: this
		integer, intent(in) :: idR1
		integer, intent(in) :: idR2
		logical :: output
		
		output = .false.
		
		if( .not. allocated(this.potentials) ) return
		if( this.potentials(idR1, idR2).getId() /= 0 ) output = .true.
	end function isThereAModel
	
	!>
	!! @brief Returns true if the reaction associated with the label is forbidden
	!!        according with the user
	!!
	function isForbidden( this, label ) result ( output )
		class(FragmentsDB), intent(in) :: this
		type(String), intent(in) :: label
		logical :: output
		
		integer :: i
		
		output = .false.
		do i=1,size(this.forbiddenReactions)
			if( label == this.forbiddenReactions(i) ) then
				output = .true.
				return
			end if
		end do
	end function isForbidden
	
	!>
	!! @brief
	!!
	function energyReference( this ) result ( output )
		class(FragmentsDB), intent(in) :: this
		real(8) :: output
		
		output = this.energyReference_
	end function energyReference
	
	!>
	!! @brief
	!!
	subroutine setEnergyReference( this, energyReference )
		class(FragmentsDB) :: this
		real(8) :: energyReference
		
		this.energyReference_ = energyReference
	end subroutine setEnergyReference
	
	!>
	!! @brief Test method
	!!
	subroutine FragmentsDB_test()
		integer :: i, j
		type(FragmentsDB) :: db
		type(String), allocatable :: fragmentsTable(:)
		type(String), allocatable :: potential(:)
		character(:), allocatable :: forbidden
		character(100), allocatable :: tokens(:), items(:)
		real(8) :: r, rBuffer
		type(String) :: strBuffer
		
! 		write(*,*) "-------------------------------------------"
! 		write(*,*) " Checking for extendFragmentsListName Method"
! 		write(*,*) "-------------------------------------------"
! 		write(*,*) "  input = 10*sHe+C10s+2B"
! 		strBuffer = this.extendFragmentsListName( "10*sHe+C10s+2B" )
! 		write(*,*) " output = "//trim(strBuffer.fstr)
! 		write(*,*) "  input = sHe+10*C10s+2B"
! 		strBuffer = this.extendFragmentsListName( "sHe+10*C10s+2B" )
! 		write(*,*) " output = "//trim(strBuffer.fstr)
! 		write(*,*) "  input = sHe+3*C10s+2*2B"
! 		strBuffer = this.extendFragmentsListName( "sHe+3*C10s+2*2B" )
! 		write(*,*) " output = "//trim(strBuffer.fstr)
		
! 		allocate( fragmentsTable(16) )
! 		
! 		!                     Label    Z  G    M  L  SYM         geomFile            Eelec       ZPE  Chann.Vib  Jmax
! 		fragmentsTable(1)  = "       sC1    0  F    1  2    0   C1.xyz             -1026.581828  0.000000"
! 		fragmentsTable(2)  = "       tC1    0  F    3  1    0   C1.xyz             -1028.031662  0.000000"
! 		fragmentsTable(3)  = "      slC2    0  T    1  0    2   C2S-linear.xyz     -2062.282893  0.232200      2*tC1   300"
! 		fragmentsTable(4)  = "      tlC2    0  T    3  1    2   C2T-linear.xyz     -2061.703744  0.209800      2*tC1   300"
! 		fragmentsTable(5)  = "      slC3    0  T    1  0    2   C3S-linear.xyz     -3097.388207  0.049062   tC1,slC2   400"
! 		fragmentsTable(6)  = "      tlC3    0  T    3  1    2   C3T-linear.xyz     -3095.315621  0.047823   tC1,slC2   400"
! 		fragmentsTable(7)  = "      scC3    0  F    1  0    2   C3S-cyclic.xyz     -3096.460007  0.176400   tC1,slC2   300"
! 		fragmentsTable(8)  = "      tcC3    0  F    3  0    6   C3T-cyclic.xyz     -3096.641729  0.161110   tC1,slC2   300"
! 		fragmentsTable(9)  = "      slC4    0  T    1  0    2   C4S-linear.xyz     -4129.754238  0.067399   tC1,slC3  1000"
! 		fragmentsTable(10) = "      tlC4    0  T    3  0    2   C4T-linear.xyz     -4129.134441  0.068026   tC1,slC3  1000"
! 		fragmentsTable(11) = "      scC4    0  F    1  0    4   C4S-cyclic.xyz     -4130.294751  0.099871   tC1,slC3   300"
! 		fragmentsTable(12) = "      tcC4    0  F    3  0    4   C4T-cyclic.xyz     -4129.378872  0.089866   tC1,slC3   300"
! 		fragmentsTable(13) = "      slC5    0  T    1  0    2   C5S-linear.xyz     -5165.261895  0.066190  slC2,slC3  1000"
! 		fragmentsTable(14) = "      tlC5    0  T    3  1    2   C5T-linear.xyz     -5162.952781  0.078498  slC2,slC3  1000"
! 		fragmentsTable(15) = "      scC5    0  F    1  0    1   C5S-cyclic.xyz     -5160.722344  0.137100  slC2,slC3   300"
! 		fragmentsTable(16) = "      tcC5    0  F    3  0    1   C5T-cyclic.xyz     -5162.267819  0.093690  slC2,slC3   300"
! 		
! 		call db.init( fragmentsTable )
! 		
! 		write(*,*) "getEelecFromName('slC3,tcC4,sC1') = ", (-3097.388207+0.049062-4129.378872+0.089866-1026.581828), " eV"
! 		write(*,*) db.getEelecFromName( "slC3,tcC4,sC1" )/eV, " eV"
! 		
! 		write(*,*) "getEelecFromName('slC3,2*tcC4,3*sC1') = ", (-3097.388207+0.049062+2.0*(-4129.378872+0.089866)+3.0*(-1026.581828)), " eV"
! 		write(*,*) db.getEelecFromName( "slC3,2*tcC4,3*sC1" )/eV, " eV"
! 		
! 		forbidden = "  slC2,slC4,slC4  -->  tlC2,slC3,tlC2,tlC2  "
! 		call FString_split( forbidden, tokens, "->" )
! 		
! 		do i=1,size(tokens)
! 			call FString_split( tokens(i), items, "," )
! 			
! 			do j=1,size(items)
! 				write(*,*) i, trim(adjustl(items(j))), db.getIdFromName( items(j) )
! 			end do
! 			write(*,*) ""
! 		end do
! 		
! 		allocate( potential(3) )
! 		potential(1) = "  tC1+tC1->tlC2    3.0   MORSE(1.5,1.0,1.0,-1.5)  "
! 		potential(2) = "  slC2+slC3->scC5  1.5   MORSE(10.5,1.0,5.0,-10.5)"
! 		potential(3) = "  slC4+tC1->slC5   1.5   MORSE(0.5,1.0,3.0,-0.5)  "
! 		
! 		call db.setPotentialTable( potential )
		
! 		write(*,*) "------------------------------------"
! 		write(*,*) " Checking reactivity"
! 		write(*,*) "------------------------------------"
! 		
! 		call clist.init(2)
! 		call clist.set( 1, db.clusters(2) )
! 		call clist.set( 2, db.clusters(2) )
! 		
! 		call clist.show( formatted=.true. )
! 		clist.kineticEnergy = 3.0_8*eV
! 		rBuffer = clist.LnW()

! 		r = 0.0_8
! 		do while( r <= 10.0_8 )
! 			write(*,*) r, db.potential( db.getIdFromName( "tC1" ), db.getIdFromName( "tC1" ), r*angs )/eV
! 			r = r + 0.2
! 		end do
		
! 		do i=1,size(clusters)
! 			write(*,*) i, clusters(i).formula()
! 		end do
		
	end subroutine FragmentsDB_test
	
end module FragmentsDB_
