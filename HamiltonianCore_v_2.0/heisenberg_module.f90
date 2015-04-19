!    Copyright (C) 2009  M. L. Wall
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
MODULE heisenberg_module
!
! Purpose: Module to construct heisenberg spin chain operators/Hamiltonian
!for OpenSourceTEBD v1.0
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

!Heisenberg Parameters
REAL(KIND=rKind) :: Jx, Jy, Jz!Spin Couplings
REAL(KIND=rKind) :: Jz1, Jz2
REAL(KIND=rKind) :: magH !magnetic field


CONTAINS

INTEGER FUNCTION heisenbergLocalDim()
!
!Purpose: Calculate local dimension for heisenberg spin chain.
!
IMPLICIT NONE
IF(spin.lt.0.5) THEN
STOP "Heisenberg spin must be at least 1/2!"
END IF
	heisenbergLocalDim = spinSize
END FUNCTION heisenbergLocalDim

SUBROUTINE CreateHeisenbergOps()
!
!Purpose: Create the operators needed to define
!and characterize the heisenberg spin chain
!
IMPLICIT NONE

INTEGER :: i, j, k, d, dsq
REAL(KIND=rKind) :: mi, mj
COMPLEX(KIND=rKind) :: preFactor, eye

		d=localSize
		dsq = d*d
		eye = CMPLX(0.0, 1.0, KIND=rKind)

!Otherwise just allocate the operator
		ALLOCATE(one_op%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate one for HeisenbergOps!'
			END IF
		ALLOCATE(Sx_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sx for HeisenbergOps!'
			END IF
		ALLOCATE(Sy_opS%m(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sy for HeisenbergOps!'
			END IF
		ALLOCATE(Sz_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sz for HeisenbergOps!'
			END IF
		ALLOCATE(Ssq_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Ssq for HeisenbergOps!'
			END IF

!Define the unity operator
		one_op%mr=0.0_rKind
		DO j=1,d,1
			one_op%mr(j,j)=1.0_rKind
		END DO
!Define the spin operators	
Sx_opS%mr=0.0_rKind
Sy_opS%m=0.0_rKind
Sz_opS%mr=0.0_rKind

	DO i=1,spinSize
	Sz_opS%mr(i,i)=spin-(i-1)*1.0_rKind
	IF(i.ne.spinSize) THEN
		Sx_opS%mr(i,i+1)=SQRT((spin+(spin-(i-1)*1.0_rKind))*(spin+1.0_rKind-(spin-(i-1)*1.0_rKind)))
		Sy_opS%m(i,i+1)=SQRT((spin+(spin-(i-1)*1.0_rKind))*(spin+1.0_rKind-(spin-(i-1)*1.0_rKind)))
	END IF
	END DO

Sx_opS%mr=0.5_rKind*(Sx_opS%mr+TRANSPOSE(Sx_opS%mr))
Sy_opS%m=eye*0.5_rKind*(TRANSPOSE(Sy_opS%m)-Sy_opS%m)

	Ssq_opS%mr=MATMUL(Sx_opS%mr,Sx_opS%mr)+REAL(MATMUL(Sy_opS%m,Sy_opS%m)) &
					+MATMUL(Sz_opS%mr,Sz_opS%mr)


		IF(print_switch) THEN
		PRINT *, "Heisenberg operators created!"
		END IF

END SUBROUTINE CreateHeisenbergOps

SUBROUTINE DestroyHeisenbergOps()
!
!Purpose: Deallocate heisenberg operators.
!
IMPLICIT NONE
INTEGER ::  k
		
!Deallocate scalar operators
		DEALLOCATE(Sx_opS%mr, Sy_opS%m, Sz_opS%mr, Ssq_opS%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate HeisenbergOps!'
			END IF
		DEALLOCATE(one_op%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate HeisenbergOps!'
			END IF

		IF(print_switch) THEN
		PRINT *, "Heisenberg operators destroyed!"
		END IF
	END SUBROUTINE DestroyHeisenbergOps


SUBROUTINE HamiltonianHeisenberg(H , Jx, Jy, Jz, magField)
!
!Purpose: Create TEBD form of heisenberg chain, where Jx, Jy, and Jz are the couplings
!Between spins in the x, y, and z directions and magField is an optional scalar argument
!specifying the presence of a uniform magnetic field.  
!
! See manual for more detail
!
IMPLICIT NONE
		TYPE(matrix), POINTER :: H(:)
		REAL(KIND=rKind), INTENT(IN) :: Jx, Jy, Jz
		REAL(KIND=rKind), INTENT(IN), OPTIONAL :: magField
		INTEGER :: i
		DO i = 1, (systemSize-1)
			H(i)%m = Jx*TensorProd(Sx_opS%mr,Sx_opS%mr)+Jy*TensorProd(Sy_opS%m,Sy_opS%m)+Jz*TensorProd(Sz_opS%mr,Sz_opS%mr)
		IF(PRESENT(magField)) THEN
		H(i)%m=H(i)%m+magField*HamiOneSite(Sz_opS)
		END IF
		END DO

!Boundary corrections
		IF(PRESENT(magField)) THEN
		H(1)%m=H(1)%m+magField*HamiLeft(Sz_opS)
		H(systemSize-1)%m=H(systemSize-1)%m+magField*HamiRight(Sz_opS)
		END IF

		IF(BoundaryCond=='P') THEN
			H(systemSize)%m = Jx*TensorProd(Sx_opS%mr,Sx_opS%mr)+Jy*TensorProd(Sy_opS%m,Sy_opS%m)+Jz*TensorProd(Sz_opS%mr,Sz_opS%mr)
		END IF


END SUBROUTINE HamiltonianHeisenberg

SUBROUTINE perturbedChain(H, JL, JR, Delta, magField)
IMPLICIT NONE
		TYPE(matrix), POINTER :: H(:)
		REAL(KIND=rKind), INTENT(IN) :: JL, JR, Delta
		REAL(KIND=rKind), INTENT(IN), OPTIONAL :: magField
		INTEGER :: i
		DO i = 1, (systemSize-1)
		IF(i.eq.FLOOR((systemSize-1)/2.)) THEN
			H(i)%m = JL*TensorProd(Sx_opS%mr,Sx_opS%mr)+JL*TensorProd(Sy_opS%m,Sy_opS%m)+Delta*TensorProd(Sz_opS%mr,Sz_opS%mr)
                ELSE IF (i.eq.FLOOR((systemSize+1)/2.)) THEN
                        H(i)%m = JR*TensorProd(Sx_opS%mr,Sx_opS%mr)+JR*TensorProd(Sy_opS%m,Sy_opS%m)+Delta*TensorProd(Sz_opS%mr,Sz_opS%mr)
		ELSE
			H(i)%m = TensorProd(Sx_opS%mr,Sx_opS%mr)+TensorProd(Sy_opS%m,Sy_opS%m)+Delta*TensorProd(Sz_opS%mr,Sz_opS%mr)
		END IF
		IF(PRESENT(magField)) THEN
		H(i)%m=H(i)%m+magField*HamiOneSite(Sz_opS)
		END IF
		END DO

!Boundary corrections
		IF(PRESENT(magField)) THEN
		H(1)%m=H(1)%m+magField*HamiLeft(Sz_opS)
		H(systemSize-1)%m=H(systemSize-1)%m+magField*HamiRight(Sz_opS)
		END IF

		IF(BoundaryCond=='P') THEN
			H(systemSize)%m = TensorProd(Sx_opS%mr,Sx_opS%mr)+TensorProd(Sy_opS%m,Sy_opS%m)+Jz2*TensorProd(Sz_opS%mr,Sz_opS%mr)
		END IF
END SUBROUTINE perturbedChain

SUBROUTINE J1J2Chain(H , Jz1, Jz2, magField)
!
!Purpose: Create TEBD form of a non-uniform heisenberg chain made of a junction of 2 half-chains
!linked together in the middle of the system. Here Jx and Jy are set to 1
!and Jz1 (resp. Jz2) is the anisotropy parameter (coupling in the z-direction) on the left (resp. right) 
!hand-side of the interface. 
!Between spins in the x, y, and z directions and magField is an optional scalar argument
!specifying the presence of a uniform magnetic field.  
!
! G.P. modified
!
IMPLICIT NONE
		TYPE(matrix), POINTER :: H(:)
		REAL(KIND=rKind), INTENT(IN) :: Jz1, Jz2
		REAL(KIND=rKind), INTENT(IN), OPTIONAL :: magField
		INTEGER :: i
		DO i = 1, (systemSize-1)
		IF(i.le.FLOOR((systemSize-1)/2.)) THEN
			H(i)%m = TensorProd(Sx_opS%mr,Sx_opS%mr)+TensorProd(Sy_opS%m,Sy_opS%m)+Jz1*TensorProd(Sz_opS%mr,Sz_opS%mr)
		ELSE
			H(i)%m = TensorProd(Sx_opS%mr,Sx_opS%mr)+TensorProd(Sy_opS%m,Sy_opS%m)+Jz2*TensorProd(Sz_opS%mr,Sz_opS%mr)
		END IF
		IF(PRESENT(magField)) THEN
		H(i)%m=H(i)%m+magField*HamiOneSite(Sz_opS)
		END IF
		END DO

!Boundary corrections
		IF(PRESENT(magField)) THEN
		H(1)%m=H(1)%m+magField*HamiLeft(Sz_opS)
		H(systemSize-1)%m=H(systemSize-1)%m+magField*HamiRight(Sz_opS)
		END IF

		IF(BoundaryCond=='P') THEN
			H(systemSize)%m = TensorProd(Sx_opS%mr,Sx_opS%mr)+TensorProd(Sy_opS%m,Sy_opS%m)+Jz2*TensorProd(Sz_opS%mr,Sz_opS%mr)
		END IF


END SUBROUTINE J1J2Chain

SUBROUTINE SetupHeisenbergName(baseName,diRectory)
!
!Purpose: Begin a file name in the directory diRectory that contains all Heisenberg parameters
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: diRectory

CALL createFileName(baseName,diRectory)
CALL appendBaseName(baseName,'Heisen_')
CALL appendBaseName(baseName,'L',systemSize)
CALL appendBaseName(baseName,'Chi',chiMax)
CALL appendBaseName(baseName,'spin',1,spin)
CALL appendBaseName(baseName,'Jx',2,Jx)
CALL appendBaseName(baseName,'Jy',2,Jy)
CALL appendBaseName(baseName,'Jz',2,Jz)
IF(magH.ne.0.0) THEN
CALL appendBaseName(baseName,'magH',2,magH)
END IF
baseName=TRIM(baseName)//TRIM(ADJUSTL(BoundaryCond))//'BC'

END SUBROUTINE SetupHeisenbergName



END MODULE heisenberg_module

