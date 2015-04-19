!    Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009  J. Williams, I. Danshita, R. Mishmash, D. Schirmer, M. L. Wall
!    This file is part of OpenSourceTEBD.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
MODULE propagation_module
!
! Purpose: Module to propagate TEBD form wavefunctions
! for OpenSourceTEBD v1.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   9/28/09   M. L. Wall	v2.0 release
!
USE system_parameters
USE TEBDtools_module
USE Hamiltonian_tools_module
USE bose_hubbard_module
USE fermi_hubbard_module
USE spinS_module
USE rotation_module
USE heisenberg_module
USE local_operations_module	
USE observables_module
USE io_module
IMPLICIT NONE

	
CONTAINS

SUBROUTINE TrotterStep(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2)
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(INOUT) :: totalTruncerr


IF(trotterOrder==2) THEN
	IF(BoundaryCond=='O') THEN
		CALL TrotterStep2ndOrder(Udt, Gammas, Lambdas, totalTruncerr)
	ELSE IF(BoundaryCond=='P') THEN
		CALL TrotterStep2ndOrderPBC(Udt, Gammas, Lambdas, totalTruncerr)
	ELSE
		STOP "BoundCond speicification not recognized in TrotterStep!"
	END IF
ELSE IF(trotterOrder==5) THEN
	IF(BoundaryCond=='O') THEN
		CALL TrotterStep5thOrder(Udt, Gammas, Lambdas, totalTruncerr)
	ELSE IF(BoundaryCond=='P') THEN
		CALL TrotterStep5thOrderPBC(Udt, Gammas, Lambdas, totalTruncerr)
	ELSE
		STOP "BoundCond speicification not recognized in TrotterStep!"
	END IF
	
ELSE
	STOP "trotterOrder specifier not recognized in TrotterStep!  Use 2 or 5"
END IF

END SUBROUTINE TrotterStep


SUBROUTINE TrotterStep2ndOrder(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2)
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
	totalTruncerr = 0.0_rKind
!!! Operate Exp(-i Hodd dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Heven dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd dt/2)			
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp )
		totalTruncerr = totalTruncerr + trunctemp
	END DO
END SUBROUTINE TrotterStep2ndOrder

SUBROUTINE TrotterStep2ndOrderPBC(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion using swapping routines
! to deform the chain so that the edge propagator may be efficiently applied
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
REAL(KIND=rKind) :: trunctemp
INTEGER :: i


	totalTruncerr = 0.0_rKind
!!! Operate Exp(-i Hodd dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
IF(MOD(systemSize,2)==1) THEN
!!! Operate Exp(-i Heven dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge dt)
		CALL TwoSiteOp(1, Udt(systemSize)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

ELSE

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge dt)
		CALL TwoSiteOp(1, Udt(systemSize)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
END IF

!!! Operate Exp(-i Hodd dt/2)			
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp )
		totalTruncerr = totalTruncerr + trunctemp
	END DO
END SUBROUTINE TrotterStep2ndOrderPBC

SUBROUTINE TrotterStep5thOrder(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2)
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
	totalTruncerr = 0.0_rKind
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Heven (1-2*theta) dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
END SUBROUTINE TrotterStep5thOrder

SUBROUTINE TrotterStep5thOrderPBC(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2)
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
	totalTruncerr = 0.0_rKind

	!Odd numbers of sites
	IF(MOD(systemSize,2)==1) THEN

!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(3*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)*theta dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(3*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+3)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta*(1-2*theta) dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+4)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)**2 dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+3)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+4)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta*(1-2*theta) dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+3)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(3*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)*theta dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOp(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(3*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Even numbers of sites
ELSE

!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta dt)
		CALL TwoSiteOp(1, Udt(2*(systemSize-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta) dt)
		CALL TwoSiteOp(1, Udt(2*(systemSize-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-2*theta) dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta dt)
		CALL TwoSiteOp(1, Udt(2*(systemSize-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL Swapping(i,trunctemp,Gammas,Lambdas)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOp(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

END IF
END SUBROUTINE TrotterStep5thOrderPBC


SUBROUTINE CanonicalFormAll(Gammas,Lambdas)
!
!Purpose: Make all Bipartite splittings canonical. 
!Used to reorthogonalize the Schmidt basis after an ITP timestep.
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: k
	DO k=1,(systemSize-1)
		CALL CanonicalForm(Lambdas(k)%v,Gammas(k)%t,Lambdas(k+1)%v,Gammas(k+1)%t,Lambdas(k+2)%v)
	END DO
	DO k=(systemSize-1),1,(-1)
		CALL CanonicalForm(Lambdas(k)%v,Gammas(k)%t,Lambdas(k+1)%v,Gammas(k+1)%t,Lambdas(k+2)%v)
	END DO
END SUBROUTINE CanonicalFormAll

SUBROUTINE ImagTimeProp(H, GammasOuter, LambdasOuter, chiIn, intDegFree)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!and is the internal degree of freedom component whose number we wish to see (should be the ground state)
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: chiIn
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: GammasOuter(:), GammasInner(:)
TYPE(vector), POINTER :: LambdasOuter(:), LambdasInner(:)
TYPE(matrix), POINTER :: Uitp(:)
INTEGER, INTENT(IN),OPTIONAL :: intDegFree
COMPLEX(KIND=rKind) :: eye, dt
REAL(KIND=rKind) :: LrgstLamsBfr(systemSize+1), LrgstLamsAft(systemSize+1), GapBfrAft(systemSize+1) 
REAL(KIND=rKind) :: numberInner, cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l, chi, chiInc

!!! Define the trotter time steps for ITP. 
	eye=CMPLX(0.0,1.0,KIND=rKind)
	dt=-eye*dtITP
!!! Construct the Imaginary time propagator
	CALL AllocateProp(Uitp)
	CALL ConstructPropagators(H, Uitp, dt)
!!! This 'if' statement is for preventing the chi increment to be 0 when chiIn=chiMax.		
	IF(chiIn==chiMax) THEN
		chiInc=1
	ELSE
		chiInc=chiMax-chiIn
	END IF

!!! We conduct the ITP for both chi=chiMin and chiMax for efficiency.
!!! The ITP for chi=chiMin usually creates a good initial state for the ITP for chi=chiMax. 			
	DO i=chiIn,chiMax,chiInc
		chi=i

IF(print_switch) THEN
	PRINT *, 'chi', chi
END IF

!!! The convergence criterion in the ITP is tighter for chi=chiMax than for chi=chiMin.
	IF((chiIn==chiMax).OR.(chi==chiMax)) THEN
		convCr=convCriterion2
	ELSE
		convCr=convCriterion1
	END IF
	
!!! Initialize the parameters for SVD.
CALL SVDInit(chi)

!!! Initialize MPS for ITP.			
		CALL AllocateGamLam(GammasInner, LambdasInner, chi)
		CALL CopyGamLam(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
		CALL DeallocateGamLam(GammasOuter, LambdasOuter)

!Store the largest lambda of each splitting initially
		DO l=1,systemSize+1
			LrgstLamsBfr(l)=LambdasInner(l)%v(1)
		END DO

!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
!!! or until the convergence criterion is satisfied.
		DO j=1,maxITPsteps
!!! This 'if' statement is for judging the convergence.
!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
			IF(MOD(j,stepsForJudge)==0) THEN
			
!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN
					CALL TotalNumber(numberInner, GammasInner, LambdasInner, intDegFree)
!Internal degree(s) of freedom absent			
ELSE
					CALL TotalNumber(numberInner, GammasInner, LambdasInner)
END IF
					CALL TotalEnergy(energy,H, GammasInner, LambdasInner)

!Store the largest lambda of each splitting after stepsForJudge time steps
				DO l=1,systemSize+1
					LrgstLamsAft(l)=LambdasInner(l)%v(1)
!Find the percent differences of each lambda
					GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
				END DO
!The lambda with the largest percent difference after stepsForJudge time steps determines convergence
				testTol=MAXVAL(GapBfrAft)
!Store the location of the lambda with the largest percent difference after stepsForJudge time steps 
				lrgstPoint=MAXLOC(GapBfrAft)
				PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(lrgstPoint(1))%v(1), &
						'found at position', lrgstPoint(1)
				PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr
IF(PRESENT(intDegFree)) THEN
				PRINT *, 'number in the',intDegFree,'mode', numberInner, 'Energy',energy
ELSE
				PRINT *, 'number', numberInner, 'Energy',energy
END IF
				IF(testTol < convCr) EXIT
!Store the new largest lambda of each splitting
				DO l=1,systemSize+1
					LrgstLamsBfr(l)=LrgstLamsAft(l)
				END DO
			END IF

!Time step
			CALL TrotterStep(Uitp, GammasInner, LambdasInner, totalTruncerr)
!Reorthogonalize
			CALL CanonicalFormAll(GammasInner,LambdasInner)
			
		END DO

!!! Deallocate the parameters for SVD.
CALL SVDFinish()

!!! Reset GammasOuter and LambdasOuter.		
		CALL AllocateGamLam(GammasOuter, LambdasOuter, chi)
		CALL CopyGamLam(GammasOuter, LambdasOuter, GammasInner, LambdasInner)
!Clean up
		CALL DeallocateGamLam(GammasInner, LambdasInner)
	END DO
		
!	chiIn=chi	
	CALL DeallocateProp(Uitp) 

END SUBROUTINE ImagTimeProp

SUBROUTINE ImagTimePropSpin(H, GammasOuter, LambdasOuter, chiIn)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: chiIn
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: GammasOuter(:), GammasInner(:)
TYPE(vector), POINTER :: LambdasOuter(:), LambdasInner(:)
TYPE(matrix), POINTER :: Uitp(:)
COMPLEX(KIND=rKind) :: eye, dt, syInner
REAL(KIND=rKind) :: LrgstLamsBfr(systemSize+1), LrgstLamsAft(systemSize+1), GapBfrAft(systemSize+1) 
REAL(KIND=rKind) :: sxInner,szInner, cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l, chi, chiInc

!!! Define the trotter time steps for ITP. 
	eye=CMPLX(0.0,1.0,KIND=rKind)
	dt=-eye*dtITP
!!! Construct the Imaginary time propagator
	CALL AllocateProp(Uitp)
	CALL ConstructPropagators(H, Uitp, dt)
!!! This 'if' statement is for preventing the chi increment to be 0 when chiIn=chiMax.		
	IF(chiIn==chiMax) THEN
		chiInc=1
	ELSE
		chiInc=chiMax-chiIn
	END IF

!!! We conduct the ITP for both chi=chiMin and chiMax for efficiency.
!!! The ITP for chi=chiMin usually creates a good initial state for the ITP for chi=chiMax. 			
	DO i=chiIn,chiMax,chiInc
		chi=i

IF(print_switch) THEN
	PRINT *, 'chi', chi
END IF

!!! The convergence criterion in the ITP is tighter for chi=chiMax than for chi=chiMin.
	IF((chiIn==chiMax).OR.(chi==chiMax)) THEN
		convCr=convCriterion2
	ELSE
		convCr=convCriterion1
	END IF
	
!!! Initialize the parameters for SVD.
CALL SVDInit(chi)

!!! Initialize MPS for ITP.			
		CALL AllocateGamLam(GammasInner, LambdasInner, chi)
		CALL CopyGamLam(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
		CALL DeallocateGamLam(GammasOuter, LambdasOuter)

!Store the largest lambda of each splitting initially
		DO l=1,systemSize+1
			LrgstLamsBfr(l)=LambdasInner(l)%v(1)
		END DO

!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
!!! or until the convergence criterion is satisfied.
		DO j=1,maxITPsteps
!!! This 'if' statement is for judging the convergence.
!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
			IF(MOD(j,stepsForJudge)==0) THEN
			
					CALL TotalOneSite(sxInner,Sx_opS%mr, GammasInner, LambdasInner)
					CALL TotalOneSite(syInner,Sy_opS%m, GammasInner, LambdasInner)
					CALL TotalOneSite(szInner,Sz_opS%mr, GammasInner, LambdasInner)
					CALL TotalEnergy(energy,H, GammasInner, LambdasInner)

!Store the largest lambda of each splitting after stepsForJudge time steps
				DO l=1,systemSize+1
					LrgstLamsAft(l)=LambdasInner(l)%v(1)
!Find the percent differences of each lambda
					GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
				END DO
!The lambda with the largest percent difference after stepsForJudge time steps determines convergence
				testTol=MAXVAL(GapBfrAft)
!Store the location of the lambda with the largest percent difference after stepsForJudge time steps 
				lrgstPoint=MAXLOC(GapBfrAft)
				PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(lrgstPoint(1))%v(1), &
						'found at position', lrgstPoint(1)
				PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr
				PRINT *, 'Sx', SxInner,'Sy', SyInner,'Sz', SzInner, 'Energy',energy

				IF(testTol < convCr) EXIT
!Store the new largest lambda of each splitting
				DO l=1,systemSize+1
					LrgstLamsBfr(l)=LrgstLamsAft(l)
				END DO
			END IF

!Time step
			CALL TrotterStep(Uitp, GammasInner, LambdasInner, totalTruncerr)
!Reorthogonalize
			CALL CanonicalFormAll(GammasInner,LambdasInner)
			
		END DO

!!! Deallocate the parameters for SVD.
CALL SVDFinish()

!!! Reset GammasOuter and LambdasOuter.		
		CALL AllocateGamLam(GammasOuter, LambdasOuter, chi)
		CALL CopyGamLam(GammasOuter, LambdasOuter, GammasInner, LambdasInner)
!Clean up
		CALL DeallocateGamLam(GammasInner, LambdasInner)
	END DO
		
!	chiIn=chi	
	CALL DeallocateProp(Uitp) 

END SUBROUTINE ImagTimePropSpin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Number conserving method starts !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TrotterStepNC(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr, intDegFree)
!
!Purpose: Propagate a TEBD form wavefunction via the Trotter expansion keeping number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
totalTruncerr=0.0_rKind
!Internal degree(s) of freedom present	
IF(PRESENT(intDegfree)) THEN

IF(trotterOrder==2) THEN
	IF(BoundaryCond=='O') THEN
		CALL TrotterStep2ndOrderNC(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr, intDegFree)
	ELSE IF(BoundaryCond=='P') THEN
		CALL TrotterStep2ndOrderPBCNC(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr, intDegFree)
	ELSE
		STOP "BoundCond speicification not recognized in TrotterStep!"
	END IF
ELSE IF(trotterOrder==5) THEN
	IF(BoundaryCond=='O') THEN
		CALL TrotterStep5thOrderNC(Udt,Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr, intDegFree)
	ELSE IF(BoundaryCond=='P') THEN
		CALL TrotterStep5thOrderPBCNC(Udt, Gammas, Lambdas,  LabelLeft, LabelRight, totalTruncerr, intDegFree)
	ELSE
		STOP "BoundCond speicification not recognized in TrotterStep!"
	END IF
ELSE
	STOP "trotterOrder specifier not recognized in TrotterStep!  Use 2 or 5"
END IF

ELSE
IF(trotterOrder==2) THEN
	IF(BoundaryCond=='O') THEN
		CALL TrotterStep2ndOrderNC(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr)
	ELSE IF(BoundaryCond=='P') THEN
		CALL TrotterStep2ndOrderPBCNC(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr)
	ELSE
		STOP "BoundCond speicification not recognized in TrotterStep!"
	END IF
ELSE IF(trotterOrder==5) THEN
	IF(BoundaryCond=='O') THEN
		CALL TrotterStep5thOrderNC(Udt,Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr)
	ELSE IF(BoundaryCond=='P') THEN
		CALL TrotterStep5thOrderPBCNC(Udt, Gammas, Lambdas,  LabelLeft, LabelRight, totalTruncerr)
	ELSE
		STOP "BoundCond speicification not recognized in TrotterStep!"
	END IF

ELSE
	STOP "trotterOrder specifier not recognized in TrotterStep!  Use 2 or 5"
END IF
END IF


END SUBROUTINE TrotterStepNC

SUBROUTINE TrotterStep2ndOrderNC(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr, intDegFree)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2) keeping number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is an array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
totalTruncerr=0.0_rKind
!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN

!!! Operate Exp(-i Hodd dt/2)	
		DO i=1,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
!!! Operate Exp(-i Heven dt)
		DO i=2,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
!!! Operate Exp(-i Hodd dt/2)			
		DO i=1,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
!Internal degree(s) of freedom absent			
ELSE

!!! Operate Exp(-i Hodd dt/2)	
		DO i=1,(systemSize-1),2
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
!!! Operate Exp(-i Heven dt)
		DO i=2,(systemSize-1),2
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
!!! Operate Exp(-i Hodd dt/2)			
		DO i=1,(systemSize-1),2
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO

END IF
END SUBROUTINE TrotterStep2ndOrderNC

SUBROUTINE TrotterStep2ndOrderPBCNC(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr, intDegFree)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2) keeping number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is an array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
totalTruncerr=0.0_rKind
!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN

!!! Operate Exp(-i Hodd dt/2)	
		DO i=1,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
IF(MOD(systemSize,2)==1) THEN
!!! Operate Exp(-i Heven dt/2)
		DO i=2,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge dt)
		CALL TwoSiteOpNC(i, Udt(systemSize)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven dt/2)
		DO i=2,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO

ELSE

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge dt)
		CALL TwoSiteOpNC(i, Udt(systemSize)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven dt)
		DO i=2,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
END IF

!!! Operate Exp(-i Hodd dt/2)			
		DO i=1,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
		END DO

!Internal degree(s) of freedom absent			
ELSE
!!! Operate Exp(-i Hodd dt/2)	
		DO i=1,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
IF(MOD(systemSize,2)==1) THEN
!!! Operate Exp(-i Heven dt/2)
		DO i=2,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
!		PRINT *, 'trunctemp from swapping', i,'=',trunctemp
	END DO

!!! Operate Exp(-i Hedge dt)
		CALL TwoSiteOpNC(i, Udt(systemSize)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven dt/2)
		DO i=2,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO

ELSE

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge dt)
		CALL TwoSiteOpNC(i, Udt(systemSize)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven dt)
		DO i=2,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO
END IF

!!! Operate Exp(-i Hodd dt/2)			
		DO i=1,(systemSize-1),2

			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
		END DO

END IF
END SUBROUTINE TrotterStep2ndOrderPBCNC



SUBROUTINE TrotterStep5thOrderNC(Udt, Gammas, Lambdas,  LabelLeft, LabelRight, totalTruncerr, intDegFree)
!
!Purpose: Propagate a TEBD form wavefunction via the 5th order Trotter expansion
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
	totalTruncerr = 0.0_rKind
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m,  Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Heven (1-2*theta) dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m,  Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m,  Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m,  Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
END SUBROUTINE TrotterStep5thOrderNC


SUBROUTINE TrotterStep5thOrderPBCNC(Udt, Gammas, Lambdas,  LabelLeft, LabelRight, totalTruncerr, intDegFree)
!
!Purpose: Propagate a TEBD form wavefunction via the 5th order Trotter expansion
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
REAL(KIND=rKind) :: trunctemp
INTEGER :: i
	totalTruncerr = 0.0_rKind

	!Odd numbers of sites
	IF(MOD(systemSize,2)==1) THEN

IF(PRESENT(intdegFree)) THEN
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)*theta dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+3)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta*(1-2*theta) dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+4)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)**2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+3)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+4)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta*(1-2*theta) dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+3)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+2)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)*theta dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
ELSE
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)*theta dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+3)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta*(1-2*theta) dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+4)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)**2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+3)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+4)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta*(1-2*theta) dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta*(1-2theta) dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+3)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta)*theta dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-theta)*theta dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta^2 dt)
		CALL TwoSiteOpNC(1, Udt(3*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Heven theta^2 dt/2)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2+3*(i-2)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO


!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(3*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
END IF
!Even numbers of sites
ELSE

IF(PRESENT(intDegFree)) THEN
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta dt)
		CALL TwoSiteOpNC(1, Udt(2*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta) dt)
		CALL TwoSiteOpNC(1, Udt(2*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-2*theta) dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta dt)
		CALL TwoSiteOpNC(1, Udt(2*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
ELSE
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta dt)
		CALL TwoSiteOpNC(1, Udt(2*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge (1-2*theta) dt)
		CALL TwoSiteOpNC(1, Udt(2*(systemSize-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven (1-2*theta) dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m,Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd (1-theta) dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+2)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!Swap L to 2
	DO i=systemSize-1,2,(-1)
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Hedge theta dt)
		CALL TwoSiteOpNC(1, Udt(2*(systemSize-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp

!Swap back
	DO i=2,systemSize-1,1
		CALL SwappingNC(i, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
			totalTruncerr = totalTruncerr + trunctemp
	END DO

!!! Operate Exp(-i Heven theta dt)
	DO i=2,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!!! Operate Exp(-i Hodd theta dt/2)	
	DO i=1,(systemSize-1),2
		CALL TwoSiteOpNC(i, Udt(2*(i-1)+1)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
END IF
END IF
END SUBROUTINE TrotterStep5thOrderPBCNC


SUBROUTINE CanonicalFormAllNC(Gammas, Lambdas, LabelLeft, LabelRight,intDegFree)
!
!Purpose: Make all bipartite splittings canonical keeping number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind) :: idenMat(localSize*localSize,localSize*localSize)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: i, k
REAL(KIND=rKind) :: trunctemp

		idenMat=CMPLX(0.0,KIND=rKind)
		DO k=1,localSize*localSize
			idenMat(k,k)=CMPLX(1.0,KIND=rKind)
		END DO

!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN
		
		DO i=1,systemSize-1,1
			CALL TwoSiteOpNC(i, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		END DO
		DO i=systemSize-1,1,-1
			CALL TwoSiteOpNC(i, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
		END DO
!Internal degree(s) of freedom absent			
ELSE
		DO i=1,systemSize-1,1
			CALL TwoSiteOpNC(i, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		END DO
		DO i=systemSize-1,1,-1
			CALL TwoSiteOpNC(i, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
		END DO
END IF

END SUBROUTINE CanonicalFormAllNC

	
SUBROUTINE ImagTimePropNC(H, GammasOuter, LambdasOuter, LabelLeftOuter, LabelRightOuter, chiIn, intDegFree)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state keeping number conservaton
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!and is the internal degree of freedom component whose number we wish to see (should be the ground state)
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: chiIn
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: GammasOuter(:)
TYPE(vector), POINTER :: LambdasOuter(:)
TYPE(vectorInt), POINTER :: LabelLeftOuter(:), LabelRightOuter(:)
TYPE(tensor), POINTER :: GammasInner(:)
TYPE(vector), POINTER :: LambdasInner(:)
TYPE(vectorInt), POINTER :: LabelLeftInner(:), LabelRightInner(:)
TYPE(matrix), POINTER :: Uitp(:)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
COMPLEX(KIND=rKind) :: eye, dt
REAL(KIND=rKind) :: LrgstLamsBfr(systemSize+1), LrgstLamsAft(systemSize+1), GapBfrAft(systemSize+1) 
REAL(KIND=rKind) :: numberInner, cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l, chi, chiInc

!!! Define the time steps for ITP. 
		eye=CMPLX(0.0,1.0,KIND=rKind)
		dt=-eye*dtITP
!!! Construct the propagation operator for ITP.

		CALL AllocateProp(Uitp)
		CALL ConstructPropagators(H, Uitp, dt)

!!! This 'if' statement is for preventing The increment of Chi to be 0 when chiIn=chiMax.		
		IF(chiIn==chiMax) THEN
			chiInc=1
		ELSE
			chiInc=chiMax-chiIn
		END IF

!!! We conduct the ITP for both chi=chiMin and chiMax for efficiency.
!!! The ITP for chi=chiMin usually creates a good initial state for the ITP for chi=chiMax. 			
		DO i=chiIn,chiMax,chiInc

		chi=i
		PRINT *, chi
			
!!! The convergence criterion in the ITP is tighter for chi=chiMax than for chi=chiMin.
			IF((chiIn==chiMax).OR.(chi==chiMax)) THEN
				convCr=convCriterion2
			ELSE
				convCr=convCriterion1
			END IF
!!! Initialize MPS for ITP.			
			CALL AllocateGamLam(GammasInner, LambdasInner, chi)
			CALL AllocateLabel(LabelLeftInner, LabelRightInner, chi)
			CALL CopyGamLam(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
			CALL CopyLabel(LabelLeftInner, LabelRightInner, LabelLeftOuter, LabelRightOuter)
			CALL DeallocateGamLam(GammasOuter, LambdasOuter)
			CALL DeallocateLabel(LabelLeftOuter, LabelRightOuter)

!Store the largest lambda of each splitting initially
			DO l=1,systemSize+1
				LrgstLamsBfr(l)=LambdasInner(l)%v(1)
			END DO
		IF(print_switch) THEN	
			PRINT *, 'maxitpsteps', maxITPsteps
		END IF	
!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
!!! or until the convergence criterion is satisfied.
			DO j=1,maxITPsteps

!!! This 'if' statement is for judging the convergence.
!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
			IF(MOD(j,stepsForJudge)==0) THEN
			
!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN
					CALL TotalNumber(numberInner, GammasInner, LambdasInner, intDegFree)
!Internal degree(s) of freedom absent			
ELSE
					CALL TotalNumber(numberInner, GammasInner, LambdasInner)
END IF
					CALL TotalEnergy(energy,H, GammasInner, LambdasInner)

!Store the largest lambda of each splitting after stepsForJudge time steps
				DO l=1,systemSize+1
					LrgstLamsAft(l)=LambdasInner(l)%v(1)
!Find the percent differences of each lambda
					GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
				END DO
!The lambda with the largest percent difference after stepsForJudge time steps determines convergence
				testTol=MAXVAL(GapBfrAft)
!Store the location of the lambda with the largest percent difference after stepsForJudge time steps 
				lrgstPoint=MAXLOC(GapBfrAft)
				PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(lrgstPoint(1))%v(1), &
						'found at position', lrgstPoint(1)
				PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr
IF(PRESENT(intDegFree)) THEN
				PRINT *, 'number in the',intDegFree,'mode', numberInner, 'Energy',energy
ELSE
				PRINT *, 'number', numberInner, 'Energy',energy
END IF
				IF(testTol < convCr) EXIT
!Store the new largest lambda of each splitting
				DO l=1,systemSize+1
					LrgstLamsBfr(l)=LrgstLamsAft(l)
				END DO
			END IF
!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN
!Time step				
				CALL TrotterStepNC(Uitp, GammasInner, LambdasInner, LabelLeftInner, LabelRightInner, totalTruncerr, intDegFree)
!Reorthogonalize
				CALL CanonicalFormAllNC(GammasInner, LambdasInner, LabelLeftInner, LabelRightInner,intDegFree)
!Internal degree(s) of freedom absent			
ELSE
!Time step		
				CALL TrotterStepNC(Uitp, GammasInner, LambdasInner, LabelLeftInner, LabelRightInner, totalTruncerr)
!Reorthogonalize
				CALL CanonicalFormAllNC(GammasInner, LambdasInner, LabelLeftInner, LabelRightInner)
END IF


			END DO
!!! Reset GammasOuter and LambdasOuter.		
			CALL AllocateGamLam(GammasOuter, LambdasOuter, chi)
			CALL AllocateLabel(LabelLeftOuter, LabelRightOuter, chi)
			CALL CopyGamLam(GammasOuter, LambdasOuter, GammasInner, LambdasInner)
			CALL CopyLabel(LabelLeftOuter, LabelRightOuter, LabelLeftInner, LabelRightInner)
			CALL DeallocateGamLam(GammasInner, LambdasInner)
			CALL DeallocateLabel(LabelLeftInner, LabelRightInner)
		END DO

		CALL DeallocateProp(Uitp)
END SUBROUTINE ImagTimePropNC


END MODULE propagation_module
