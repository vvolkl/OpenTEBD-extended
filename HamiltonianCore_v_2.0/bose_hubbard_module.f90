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
MODULE bose_hubbard_module
!
! Purpose: Module to construct spinless Bose-Hubbard (BH) model operators/Hamiltonian
! for OpenSourceTEBD v1.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   9/28/09   M. L. Wall	v2.0 release
!
USE system_parameters
USE TEBDtools_module
USE io_module
USE Hamiltonian_tools_module

	
IMPLICIT NONE

INTERFACE HamiltonianBoseHubbard
MODULE PROCEDURE HamiltonianBoseHubbardScalar, HamiltonianBoseHubbardUVector,&
				HamiltonianBoseHubbardJVector, HamiltonianBoseHubbardJUVector
END INTERFACE


CONTAINS

INTEGER FUNCTION boseHubbardLocalDim()
!
!Purpose: Calculate local dimension for spinless Bose-Hubbard system.
!
IMPLICIT NONE
	boseHubbardLocalDim = maxFilling+1
END FUNCTION boseHubbardLocalDim

SUBROUTINE CreateFieldOps()
!
!Purpose: Allocate and define operators to define and characterize the Bose-Hubbard Hamiltonian
!
IMPLICIT NONE
INTEGER :: i, j, d, dsq, n, nprime
COMPLEX(KIND=rKind) :: pie, eye
	d=localSize
	dsq=d*d
	ALLOCATE(a_op%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for SpinlessOps'
			END IF
	ALLOCATE(t_op%mr(dsq,dsq), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate t for SpinlessOps'
			END IF
	ALLOCATE(one_op%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate one for SpinlessOps'
			END IF
	ALLOCATE(PBphase_op%m(d, d), STAT=statInt) ! Allocate for Pegg-Barnett phase operator.
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate PBPhase for SpinlessOps'
			END IF
	a_op%mr=CMPLX(0.0,KIND=rKind)
	one_op%mr=CMPLX(0.0,KIND=rKind)
		
	DO j=1,(d-1)
!Define the unit operator
		one_op%mr(j,j)=1.0_rKind
!Define the destruction operator
		a_op%mr(j,j+1)=REAL(SQRT(j*1.0_rKind),KIND=rKind)
	END DO
	one_op%mr(d,d)=1.0_rKind
!Define the tunneling operator	
	t_op%mr = tensorProd(TRANSPOSE(a_op%mr),a_op%mr) + TRANSPOSE(tensorProd(TRANSPOSE(a_op%mr),a_op%mr))

	! Create matrix representation of Pegg-Barnett phase operator.
	PBphase_op%m = CMPLX(0.0, KIND=rKind)
	eye = CMPLX(0.0, 1.0, KIND=rKind)
	pie=3.1415926535897_rKind
	DO n = 1, d
		PBphase_op%m(n, n) = CMPLX((d-1)*pie,KIND=rKind)/CMPLX(d,KIND=rKind)
	END DO
	
	DO n = 1, d-1
		DO nprime = d, n+1, -1
		PBphase_op%m(n, nprime) = 2.0_rKind*pie/d * 1.0_rKind/( EXP(eye*(n-nprime)*2.0_rKind*pie/d) - 1.0_rKind )
		END DO
	END DO
	
	DO n = 2, d
		DO nprime = 1, n-1
		PBphase_op%m(n, nprime) = 2.0_rKind*pie/d * 1.0_rKind/( EXP(eye*(n-nprime)*2.0_rKind*pie/d) - 1.0_rKind )
		END DO
	END DO
	
	IF(print_switch) THEN		
	PRINT *, 'Bose-Hubbard operators created!'
	END IF
END SUBROUTINE CreateFieldOps
	
SUBROUTINE DestroyFieldOps()
!
!Purpose: Deallocate the Bose-Hubbard Hamiltonian operators.
!
IMPLICIT NONE
	DEALLOCATE(a_op%mr,one_op%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SpinlessOps'
			END IF
	DEALLOCATE(t_op%mr,PBphase_op%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SpinlessOps'
			END IF
		IF(print_switch) THEN
		PRINT *, 'Bose-Hubbard operators destroyed!'
		END IF
END SUBROUTINE DestroyFieldOps

SUBROUTINE HamiltonianBoseHubbardScalar(H, jTunn, U0, mu0,  V0, extPot)
!
!Purpose: Construct the (extended) Bose-Hubbard Hamiltonian (see manual) in TEBD form
!Here jTunn is the tunneling energy, U0 is the on-site repulsion, mu0 is the chemical potential, 
!V0 is the nearest neighbor repulsion, and extPot is an optional vector argument specifing an
!arbitrary site-dependent trap.
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) ::  jTunn, U0, mu0, V0
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: extPot(:)
INTEGER :: i

	DO i=1,(systemSize-1)
H(i)%m=0.0_rKind
!One site operations
!"Left" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Left" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Two site operations
!Tunneling part
H(i)%m=H(i)%m-jTunn*t_op%mr

!Nearest-neighbor part
H(i)%m=H(i)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
	END DO

!Fox box boundary conditions first and last site get additional left and right
!one site operations, respectively

!"Left" part of U/2(n(n-1)) for first site
H(1)%m=H(1)%m+0.25_rKind*U0*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Left" part of mu*n for first site
H(1)%m=H(1)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)

!"Right" part of U/2(n(n-1)) for last site
H(systemSize-1)%m=H(systemSize-1)%m+0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Right" part of mu*n for last site
H(systemSize-1)%m=H(systemSize-1)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Application of an arbitrary external potential
IF(PRESENT(extPot)) THEN

DO i=1,(systemSize-1)
H(i)%m=H(i)%m+0.5_rKind*extPot(i)*(tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr))+0.5_rKind*extPot(i+1)*(tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))
END DO
H(systemSize-1)%m=H(systemSize-1)%m+0.5_rKind*extPot(systemSize)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
H(1)%m=H(1)%m+0.5_rKind*extPot(1)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
END IF

IF(BoundaryCond=='P') THEN
H(systemSize)%m=0.0_rKind
H(systemSize)%m=H(systemSize)%m-jTunn*t_op%mr

!Nearest-neighbor part
H(systemSize)%m=H(systemSize)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
END IF

END SUBROUTINE HamiltonianBoseHubbardScalar

SUBROUTINE HamiltonianBoseHubbardUVector(H, jTunn, U0, mu0,  V0, extPot)
!
!Purpose: Construct the (extended) Bose-Hubbard Hamiltonian (see manual) in TEBD form
!Here jTunn is the tunneling energy, U0 is the on-site repulsion, mu0 is the chemical potential, 
!V0 is the nearest neighbor repulsion, and extPot is an optional vector argument specifing an
!arbitrary site-dependent trap.
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) ::  jTunn, U0(:), mu0, V0
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: extPot(:)
INTEGER :: i

	DO i=1,(systemSize-1)
H(i)%m=0.0_rKind
!One site operations
!"Left" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0(i)*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0(i)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0(i+1)*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0(i+1)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Left" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Two site operations
!Tunneling part
H(i)%m=H(i)%m-jTunn*t_op%mr

!Nearest-neighbor part
H(i)%m=H(i)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
	END DO

!Fox box boundary conditions first and last site get additional left and right
!one site operations, respectively

!"Left" part of U/2(n(n-1)) for first site
H(1)%m=H(1)%m+0.25_rKind*U0(1)*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0(1)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Left" part of mu*n for first site
H(1)%m=H(1)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)

!"Right" part of U/2(n(n-1)) for last site
H(systemSize-1)%m=H(systemSize-1)%m+0.25_rKind*U0(systemSize-1)*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0(systemSize-1)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Right" part of mu*n for last site
H(systemSize-1)%m=H(systemSize-1)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Application of an arbitrary external potential
IF(PRESENT(extPot)) THEN

DO i=1,(systemSize-1)
H(i)%m=H(i)%m+0.5_rKind*extPot(i)*(tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr))+0.5_rKind*extPot(i+1)*(tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))
END DO
H(systemSize-1)%m=H(systemSize-1)%m+0.5_rKind*extPot(systemSize)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
H(1)%m=H(1)%m+0.5_rKind*extPot(1)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
END IF

IF(BoundaryCond=='P') THEN
H(systemSize)%m=0.0_rKind
H(systemSize)%m=H(systemSize)%m-jTunn*t_op%mr

!Nearest-neighbor part
H(systemSize)%m=H(systemSize)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
END IF

END SUBROUTINE HamiltonianBoseHubbardUVector

SUBROUTINE HamiltonianBoseHubbardJVector(H, jTunn, U0, mu0,  V0, extPot)
!
!Purpose: Construct the (extended) Bose-Hubbard Hamiltonian (see manual) in TEBD form
!Here jTunn is the tunneling energy, U0 is the on-site repulsion, mu0 is the chemical potential, 
!V0 is the nearest neighbor repulsion, and extPot is an optional vector argument specifing an
!arbitrary site-dependent trap.
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) ::  jTunn(:), U0, mu0, V0
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: extPot(:)
INTEGER :: i

	DO i=1,(systemSize-1)
H(i)%m=0.0_rKind
!One site operations
!"Left" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Left" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Two site operations
!Tunneling part
H(i)%m=H(i)%m-jTunn(i)*t_op%mr

!Nearest-neighbor part
H(i)%m=H(i)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
	END DO

!Fox box boundary conditions first and last site get additional left and right
!one site operations, respectively

!"Left" part of U/2(n(n-1)) for first site
H(1)%m=H(1)%m+0.25_rKind*U0*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Left" part of mu*n for first site
H(1)%m=H(1)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)

!"Right" part of U/2(n(n-1)) for last site
H(systemSize-1)%m=H(systemSize-1)%m+0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Right" part of mu*n for last site
H(systemSize-1)%m=H(systemSize-1)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Application of an arbitrary external potential
IF(PRESENT(extPot)) THEN

DO i=1,(systemSize-1)
H(i)%m=H(i)%m+0.5_rKind*extPot(i)*(tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr))+0.5_rKind*extPot(i+1)*(tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))
END DO
H(systemSize-1)%m=H(systemSize-1)%m+0.5_rKind*extPot(systemSize)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
H(1)%m=H(1)%m+0.5_rKind*extPot(1)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
END IF

IF(BoundaryCond=='P') THEN
H(systemSize)%m=0.0_rKind
H(systemSize)%m=H(systemSize)%m-jTunn(systemSize)*t_op%mr

!Nearest-neighbor part
H(systemSize)%m=H(systemSize)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
END IF

END SUBROUTINE HamiltonianBoseHubbardJVector

SUBROUTINE HamiltonianBoseHubbardJUVector(H, jTunn, U0, mu0,  V0, extPot)
!
!Purpose: Construct the (extended) Bose-Hubbard Hamiltonian (see manual) in TEBD form
!Here jTunn is the tunneling energy, U0 is the on-site repulsion, mu0 is the chemical potential, 
!V0 is the nearest neighbor repulsion, and extPot is an optional vector argument specifing an
!arbitrary site-dependent trap.
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) ::  jTunn(:), U0(:), mu0, V0
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: extPot(:)
INTEGER :: i


DO i=1,(systemSize-1)
H(i)%m=0.0_rKind
!One site operations
!"Left" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0(i)*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0(i)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0(i+1)*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0(i+1)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Left" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Two site operations
!Tunneling part
H(i)%m=H(i)%m-jTunn(i)*t_op%mr

!Nearest-neighbor part
H(i)%m=H(i)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
	END DO

!Fox box boundary conditions first and last site get additional left and right
!one site operations, respectively

!"Left" part of U/2(n(n-1)) for first site
H(1)%m=H(1)%m+0.25_rKind*U0(1)*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0(1)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Left" part of mu*n for first site
H(1)%m=H(1)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)

!"Right" part of U/2(n(n-1)) for last site
H(systemSize-1)%m=H(systemSize-1)%m+0.25_rKind*U0(systemSize-1)*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0(systemSize-1)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Right" part of mu*n for last site
H(systemSize-1)%m=H(systemSize-1)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Application of an arbitrary external potential
IF(PRESENT(extPot)) THEN

DO i=1,(systemSize-1)
H(i)%m=H(i)%m+0.5_rKind*extPot(i)*(tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr))+0.5_rKind*extPot(i+1)*(tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))
END DO
H(systemSize-1)%m=H(systemSize-1)%m+0.5_rKind*extPot(systemSize)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
H(1)%m=H(1)%m+0.5_rKind*extPot(1)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
END IF

IF(BoundaryCond=='P') THEN
H(systemSize)%m=0.0_rKind
H(systemSize)%m=H(systemSize)%m-jTunn(i)*t_op%mr

!Nearest-neighbor part
H(systemSize)%m=H(systemSize)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
END IF


END SUBROUTINE HamiltonianBoseHubbardJUVector

SUBROUTINE SetupBHName(baseName,diRectory)
!
!Purpose: Begin a file name in the directory diRectory that contains all Bose-Hubbard parameters
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: diRectory

CALL createFileName(baseName,diRectory)
CALL appendBaseName(baseName,'BH_')
CALL appendBaseName(baseName,'L',systemSize)
IF(ncswitch) THEN
	CALL appendBaseName(baseName,'N',totNum)
ELSE
	CALL appendBaseName(baseName,'mu',2,mu0)
END IF
CALL appendBaseName(baseName,'Chi',chiMax)
CALL appendBaseName(baseName,'jTunn',2,jTunn)
CALL appendBaseName(baseName,'U0',2,U0)
baseName=TRIM(baseName)//TRIM(ADJUSTL(BoundaryCond))//'BC'
IF(V0.ne.0.0) THEN
CALL appendBaseName(baseName,'V0',2,V0)
END IF
IF(ALLOCATED(extPot)) THEN
CALL appendBaseName(baseName,'extPot')
END IF

END SUBROUTINE SetupBHName

END MODULE bose_hubbard_module
