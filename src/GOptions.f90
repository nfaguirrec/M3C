!>
!! @brief
!!
module GOptions_
	use UnitsConverter_
	use IOStream_
	use Timer_
	
	implicit none
	public
	
	real(8) :: GOptions_zero = 1e-12
	
	real(8) :: GOptions_systemRadius = 10.0_8*angs
	real(8) :: GOptions_randomWalkStepRadius = 2.0_8*angs
	real(8) :: GOptions_overlappingRadius = 0.0_8*angs
	logical :: GOptions_useWeightedWalkStep = .false.
	logical :: GOptions_useRandomWalkers = .true.
	logical :: GOptions_useZPECorrection = .false.
	logical :: GOptions_useLCorrection = .false.
	logical :: GOptions_useLDOSContrib = .false.
	logical :: GOptions_useLWeightContrib = .false.
	real(8) :: GOptions_gammaLCorrection = 1.0_8
	logical :: GOptions_useLReference = .false.
	
	! 1) MMMC
	! 2) MMMC, Reactor
	! 3) MMMC, Reactor, RigidMoleculeList
	! 4) MMMC, Reactor, RigidMoleculeList, RigidMolecule
	logical :: GOptions_printLevel = 1
	
	! 1) NOTHING
	! 2) INFO
	! 3) INFO + WARNING
	logical :: GOptions_debugLevel = 1
	
	integer :: GOptions_indentLength = 5
	
	type(Timer) :: GOptions_timer
	
    interface GOptions_valueReport
            module procedure GOptions_iValueReport
            module procedure GOptions_rValueReport
            module procedure GOptions_sValueReport
            module procedure GOptions_rArrValueReport
            module procedure GOptions_iArrValueReport
    end interface GOptions_valueReport
	
	contains
	
	!>
	!! @brief
	!!
	subroutine GOptions_section( message, indent )
		character(*), intent(in) :: message
		integer, intent(in), optional :: indent
		
		integer :: effIndent
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(A)") ""
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") "+"//repeat("-",len_trim(message)+2)//"+"
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") "| "//trim(message)//" |"
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") "+"//repeat("-",len_trim(message)+2)//"+"
		write(STDOUT,"(A)") ""
	end subroutine GOptions_section
	
	!>
	!! @brief
	!!
	subroutine GOptions_subsection( message, indent )
		character(*), intent(in) :: message
		integer, intent(in), optional :: indent
		
		integer :: effIndent
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(A)") ""
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") repeat("-",len_trim(message)+2)
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") " "//trim(message)//" "
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") repeat("-",len_trim(message)+2)
		write(STDOUT,"(A)") ""
	end subroutine GOptions_subsection

	!>
	!! @brief
	!!
	subroutine GOptions_paragraph( message, indent )
		character(*), intent(in) :: message
		integer, intent(in), optional :: indent
		
		integer :: effIndent
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(A)") ""
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") " "//trim(message)//" "
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A)") repeat("-",len_trim(message)+2)
		write(STDOUT,"(A)") ""
	end subroutine GOptions_paragraph
	
	!>
	!! @brief
	!!
	subroutine GOptions_error( message, where, addinfo )
		character(*), intent(in) :: message
		character(*), intent(in), optional :: where
		character(*), intent(in), optional :: addinfo
		
		character(1000) :: effWhere
		character(1000) :: effAddinfo
		
		effWhere = ""
		if( present(where) ) effWhere = where

		effAddinfo = ""
		if( present(addinfo) ) effAddinfo = addinfo
		
		write(STDOUT,"(A)") ""
		write(STDERR,"(A)") "### ERROR ### "//trim(effWhere)//": "//trim(message)
		write(STDERR,"(A)") "              "//trim(effAddinfo)
		
		stop
	end subroutine GOptions_error
	
	!>
	!! @brief
	!!
	subroutine GOptions_warning( message, where, addinfo )
		character(*), intent(in) :: message
		character(*), intent(in), optional :: where
		character(*), intent(in), optional :: addinfo
		
		character(1000) :: effWhere
		character(1000) :: effAddinfo
		
		effWhere = ""
		if( present(where) ) effWhere = where

		effAddinfo = ""
		if( present(addinfo) ) effAddinfo = addinfo
		
		if( GOptions_debugLevel >= 3 ) then
			write(STDOUT,"(A)") ""
			write(STDERR,"(A)") "!!! WARNING ¡¡¡ "//trim(effWhere)//": "//trim(message)
			write(STDERR,"(A)") "                "//trim(effAddinfo)
		end if
	end subroutine GOptions_warning
	
	!>
	!! @brief
	!!
	subroutine GOptions_info( message, where, addinfo )
		character(*), intent(in) :: message
		character(*), intent(in), optional :: where
		character(*), intent(in), optional :: addinfo
		
		character(1000) :: effWhere
		character(1000) :: effAddinfo
		
		effWhere = ""
		if( present(where) ) effWhere = where
		
		effAddinfo = ""
		if( present(addinfo) ) effAddinfo = addinfo
		
		if( GOptions_debugLevel >= 2 ) then
			write(STDOUT,"(A)") ""
			write(STDOUT,"(A)") "%%% INFO %%% "//trim(effWhere)//": "//trim(message)
			write(STDOUT,"(A)") "             "//trim(effAddinfo)
		end if
	end subroutine GOptions_info
	
	!>
	!! @brief
	!!
	subroutine GOptions_rValueReport( varName, value, message, units, indent )
		character(*), intent(in) :: varName
		real(8), intent(in) :: value
		character(*), intent(in), optional :: message
		character(*), intent(in), optional :: units
		integer, intent(in), optional :: indent
		
		character(1000) :: effMessage
		character(100) :: effUnits
		integer :: effIndent
		
		effMessage = ""
		if( present(message) ) effMessage = message
		
		effUnits = ""
		if( present(units) ) effUnits = units
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A10,F20.5,A8,5X,A)") trim(varName), value, trim(effUnits), trim(effMessage)
	end subroutine GOptions_rValueReport
	
	!>
	!! @brief
	!!
	subroutine GOptions_iValueReport( varName, value, message, units, indent )
		character(*), intent(in) :: varName
		integer, intent(in) :: value
		character(*), intent(in), optional :: message
		character(*), intent(in), optional :: units
		integer, intent(in), optional :: indent
		
		character(1000) :: effMessage
		character(100) :: effUnits
		integer :: effIndent
		
		effMessage = ""
		if( present(message) ) effMessage = message
		
		effUnits = ""
		if( present(units) ) effUnits = units
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A10,I14,6X,A8,5X,A)") trim(varName), value, trim(effUnits), trim(effMessage)
	end subroutine GOptions_iValueReport

	!>
	!! @brief
	!!
	subroutine GOptions_sValueReport( varName, value, message, indent )
		character(*), intent(in) :: varName
		character(*), intent(in) :: value
		character(*), intent(in), optional :: message
		integer, intent(in), optional :: indent
		
		character(1000) :: effMessage
		integer :: effIndent
		
		effMessage = ""
		if( present(message) ) effMessage = message
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A10,A20,A8,5X,A)") trim(varName), trim(value), trim(effMessage)
	end subroutine GOptions_sValueReport
	
	!>
	!! @brief
	!!
	subroutine GOptions_rArrValueReport( varName, values, message, units, indent )
		character(*), intent(in) :: varName
		real(8), intent(in) :: values(:)
		character(*), intent(in), optional :: message
		character(*), intent(in), optional :: units
		integer, intent(in), optional :: indent
		
		character(1000) :: effMessage
		character(100) :: effUnits
		integer :: effIndent
		
		effMessage = ""
		if( present(message) ) effMessage = message
		
		effUnits = ""
		if( present(units) ) effUnits = units
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(<GOptions_indentLength*effIndent>X,A10,<size(values)>F15.5,A8,5X,A)") trim(varName), values, trim(effUnits), trim(effMessage)
	end subroutine GOptions_rArrValueReport
	
	!>
	!! @brief
	!!
	subroutine GOptions_iArrValueReport( varName, values, message, units, indent )
		character(*), intent(in) :: varName
		integer, intent(in) :: values(:)
		character(*), intent(in), optional :: message
		character(*), intent(in), optional :: units
		integer, intent(in), optional :: indent
		
		character(1000) :: effMessage
		character(100) :: effUnits
		integer :: effIndent
		
		effMessage = ""
		if( present(message) ) effMessage = message
		
		effUnits = ""
		if( present(units) ) effUnits = units
		
		effIndent = 0
		if( present(indent) ) effIndent = indent
		
		write(STDOUT,"(A10,<size(values)>I15,A8,5X,A)") trim(varName), values, trim(effUnits), trim(effMessage)
	end subroutine GOptions_iArrValueReport
	
end module GOptions_
