!    Copyright (C)  2008, 2009  M. L. Wall
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
MODULE spinS_module
!
! Purpose: Module to construct spin-S Bose-Hubbard (BH) model operators/Hamiltonian
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

!Spin-One parameters
!In addition to the Bose-Hubbard parameters above we have
REAL(KIND=rKind) :: U2 !Spin dependent interaction
REAL(KIND=rKind) :: VB !Zeeman interaction


CONTAINS

INTEGER FUNCTION spinSLocalDim()
!
!Purpose: Calculate local dimension for spin-S system.
!
IMPLICIT NONE
! (N+2s+1)!/(N!(2s+1)!)
	spinSLocalDim = BinomialCoef(maxFilling+spinSize,spinSize)
END FUNCTION spinSLocalDim


SUBROUTINE CreateSpinSops()
!
!Purpose: Create the operators needed to define
!and characterize the spin-s BH model
!
IMPLICIT NONE

INTEGER :: i, j, k, d, dsq
REAL(KIND=rKind) :: mi, mj
TYPE(matrix) :: stateList,newState
COMPLEX(KIND=rKind) :: preFactor, eye

		d=localSize
		dsq = d*d
		eye = CMPLX(0.0, 1.0, KIND=rKind)


!For the spin indexed operators we must:
!1. Allocate the number of lists
		ALLOCATE(a_opS(spinSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for SpinSOps'
			END IF

!2. Allocate each list
		DO i=1,spinSize,1
		ALLOCATE(a_opS(i)%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for SpinSOps'
			END IF
		END DO


!Otherwise just allocate the operator
		ALLOCATE(one_op%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate one for SpinSOps'
			END IF
		ALLOCATE(Sx_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sx for SpinSOps'
			END IF
		ALLOCATE(Sy_opS%m(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sy for SpinSOps'
			END IF
		ALLOCATE(Sz_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sz for SpinSOps'
			END IF
		ALLOCATE(VB_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate VB for SpinSOps'
			END IF
		ALLOCATE(Ssq_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Ssq for SpinSOps'
			END IF
		ALLOCATE(ntot_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ntot for SpinSOps'
			END IF
		ALLOCATE(ttot_opS%mr(dsq, dsq), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ttot for SpinSOps'
			END IF

!State enumeration stuff
		ALLOCATE(stateList%m(d, spinSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate stateList for SpinSOps'
			END IF
		ALLOCATE(newState%m(1, spinSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate newState for SpinSOps'
			END IF
		ALLOCATE(Conserv%vi(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Conserv for SpinSOps'
			END IF

!Define the unity operator
		one_op%mr=0.0_rKind
		DO j=1,d,1
			one_op%mr(j,j)=1.0_rKind
		END DO

!Create and index the Fock space
		CALL onsiteStateListIdof(stateList%m, spinSize)

!Define the destruction operator
		DO k=1,spinSize,1
		a_opS(k)%mr = 0.0_rKind
			DO j = 1, d
			newState%m(1, :) = stateList%m(j, :)
			preFactor = SQRT(newState%m(1, k))
			newState%m(1, k) = newState%m(1, k) - 1.0_rKind
				DO i = 1, d		
					a_opS(k)%mr(i, j) = preFactor*kronDelta(stateList%m(i, :), &
								newState%m(1, :), spinSize)
				END DO			
			END DO
		END DO
!Deallocate the state enumeration stuff
		DEALLOCATE(stateList%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate stateList for SpinSOps'
			END IF
		DEALLOCATE(newState%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate newState for SpinSOps'
			END IF


ntot_opS%mr=0.0_rKind
!Define the total number operator
	DO k=1,spinSize,1
	ntot_opS%mr=ntot_opS%mr+MATMUL(TRANSPOSE(a_opS(k)%mr), a_opS(k)%mr)	
	END DO

Sx_opS%mr=0.0_rKind
Sy_opS%m=0.0_rKind
Sz_opS%mr=0.0_rKind
VB_opS%mr=0.0_rKind

!Define the Spin operators using
!\hat{S}_{\nu}=sum_{k,k'}\hat{a}^{\dagger}_{k} [F]^{\nu}_{kk'} \hat{a}_{k'}
!Here [F]^{\nu} is the 2s+1 matrix representation (pauli matrix)
!of the spin operator in the \nu^{th} Cartesian direction
		DO i=1,spinSize,1
!S_+ elements of [F]^{\nu}
			IF(i.ne.spinSize) THEN
			mi=(2*spin+1-i)*(i)
			Sx_opS%mr=Sx_opS%mr+0.5_rKind*SQRT(mi)*MATMUL(Transpose(a_opS(i)%mr),a_opS(i+1)%mr)
			Sy_opS%m=Sy_opS%m-eye*0.5_rKind*SQRT(mi)*MATMUL(Transpose(a_opS(i)%mr),a_opS(i+1)%mr)
			END IF
			IF(i.ne.1) THEN
!S_- elements of [F]^{\nu}
			mj=(2*spin-i+2)*(1-(2-i))
			Sx_opS%mr=Sx_opS%mr+0.5_rKind*SQRT(mj)*MATMUL(Transpose(a_opS(i)%mr),a_opS(i-1)%mr)
			Sy_opS%m=Sy_opS%m+eye*0.5_rKind*SQRT(mj)*MATMUL(Transpose(a_opS(i)%mr),a_opS(i-1)%mr)
			END IF
			Sz_opS%mr=Sz_opS%mr+(spin+1-i)*MATMUL(Transpose(a_opS(i)%mr),a_opS(i)%mr)
!Define the Zeeman operator
			VB_opS%mr=VB_opS%mr+Abs(spin+1-i)*Abs(spin+1-i)*MATMUL(Transpose(a_opS(i)%mr),a_opS(i)%mr)
		END DO
!Define the total spin operator
		Ssq_opS%mr=MATMUL(Sx_opS%mr,Sx_opS%mr)+REAL(MATMUL(Sy_opS%m,Sy_opS%m)) &
						+MATMUL(Sz_opS%mr,Sz_opS%mr)


!Define the tunneling operator
	ttot_opS%mr=0.0
	DO k=1,spinSize,1
		ttot_opS%mr = ttot_opS%mr+TensorProd(Transpose(a_opS(k)%mr),a_opS(k)%mr) &
		+ Transpose(TensorProd(Transpose(a_opS(k)%mr),a_opS(k)%mr))
	END DO


END SUBROUTINE CreateSpinSops


SUBROUTINE DestroySpinSops()
!
!Purpose: Deallocate spin-S operators.
!
IMPLICIT NONE
INTEGER ::  k
		
!Deallocate scalar operators
		DEALLOCATE(Sx_opS%mr, Sy_opS%m, Sz_opS%mr, Ssq_opS%mr, ntot_opS%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SpinSOps'
			END IF
		DEALLOCATE(ttot_opS%mr,one_op%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SpinSOps'
			END IF
		DEALLOCATE(Conserv%vi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SpinSOps'
			END IF
		DEALLOCATE(VB_opS%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate VB for SpinSOps'
			END IF

!Deallocate vector operators
		DO k=1,spinSize,1
		DEALLOCATE(a_opS(k)%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SpinSOps'
			END IF
		END DO
		DEALLOCATE(a_opS, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SpinSOps'
			END IF
		IF(print_switch) THEN
		PRINT *, "Spin-S operators destroyed!"
		END IF
	END SUBROUTINE DestroySpinSops


SUBROUTINE HamiltonianSpinOne(H , J, U0, U2, mu0, VB)
!
!Purpose: Create TEBD form of spin-1 Bose-Hubbard Hamiltonian with quadratic Zeeman 
! strength VB as described in the manual.  
!
! J, U0, U2 are the usual spin-1 Bose-Hubbard parameters and mu0 is 
! the chemical potential.
! This subroutine needs to be modified if an arbitary external potenial 
! on the spin components is desired.
!
! See manual for more detail
!
IMPLICIT NONE
		TYPE(matrix), POINTER :: H(:)
		REAL(KIND=rKind), INTENT(IN) :: J, U0, U2, mu0, VB
		INTEGER :: i
		DO i = 1, (systemSize-1)
			H(i)%m = ((-1.0_rKind*mu0 - U0/2.0_rKind - U2)*TensorProd(ntot_opS%mr,one_op%mr) +&
				 U0/2.0_rKind*TensorProd(MATMUL(ntot_opS%mr,ntot_opS%mr),one_op%mr))/2.0_rKind &
				+ ((-1.0_rKind*mu0 - U0/2.0_rKind - U2)*TensorProd(one_op%mr,ntot_opS%mr) +&
				 U0/2.0_rKind*TensorProd(one_op%mr,MATMUL(ntot_opS%mr,ntot_opS%mr)))/2.0_rKind &
				+ (U2/2.0_rKind*TensorProd(Ssq_opS%mr,one_op%mr))/2.0_rKind + &
				((U2/2.0_rKind*TensorProd(one_op%mr,Ssq_opS%mr))/2.0_rKind)+ &
				(- J*ttot_opS%mr+VB*TensorProd(VB_opS%mr,one_op%mr)/2.0_rKind)+&
					VB*TensorProd(one_op%mr,VB_opS%mr)/2.0_rKind
			
		END DO

		H(1)%m = H(1)%m + ((-1.0_rKind*mu0 - U0/2.0_rKind - U2)*TensorProd(ntot_opS%mr,one_op%mr) +&
			 U0/2.0_rKind*TensorProd(MATMUL(ntot_opS%mr,ntot_opS%mr),one_op%mr))/2.0_rKind &
			+ (U2/2.0_rKind*TensorProd(Ssq_opS%mr,one_op%mr))/2.0_rKind +&
			VB*TensorProd(VB_opS%mr,one_op%mr)/2.0_rKind
		H(systemSize-1)%m = H(systemSize-1)%m + ((-1.0_rKind*mu0 - U0/2.0_rKind - U2)*TensorProd(one_op%mr,ntot_opS%mr) +&
			 U0/2.0_rKind*TensorProd(one_op%mr,MATMUL(ntot_opS%mr,ntot_opS%mr)))/2.0_rKind &
			+ (U2/2.0_rKind*TensorProd(one_op%mr,Ssq_opS%mr))/2.0_rKind +&
					VB*TensorProd(one_op%mr,VB_opS%mr)/2.0_rKind

END SUBROUTINE HamiltonianSpinOne

SUBROUTINE HamiltonianSpinS(H , J, gS, mu0, VB)
!
!Purpose: Create TEBD form of spin-S Bose-Hubbard Hamiltonian with quadratic Zeeman 
! strength VB as described in the manual.  
!
! J, and the gS(:) are the usual spin-S Bose-Hubbard parameters and mu0 is 
! the chemical potential.
! This subroutine needs to be modified if an arbitary external potenial 
! on the spin components is desired.
!
! See manual for more detail
!
IMPLICIT NONE
		TYPE(matrix), POINTER :: H(:)
		REAL(KIND=rKind), INTENT(IN) :: J, gS(:), mu0, VB
		TYPE(matrixReal) :: O_opS
		INTEGER :: i, alpha1, alpha2, sTot, alpha
		
ALLOCATE(O_opS%mr(localSize, localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate O_opS for HamiltonianSpinS'
			END IF		
		
!Spin-independent part of Hamiltonian		
		DO i = 1, (systemSize-1)
			H(i)%m = (-0.5_rKind*mu0 )*TensorProd(ntot_opS%mr,one_op%mr) +&
			(-0.5_rKind*mu0 )*TensorProd(one_op%mr,ntot_opS%mr) &
				- J*ttot_opS%mr+VB*TensorProd(VB_opS%mr,one_op%mr)/2.0_rKind+&
					VB*TensorProd(one_op%mr,VB_opS%mr)/2.0_rKind
			
		END DO

		H(1)%m = H(1)%m + (-0.5_rKind*mu0 )*TensorProd(ntot_opS%mr,one_op%mr) +&
			VB*TensorProd(VB_opS%mr,one_op%mr)/2.0_rKind
		H(systemSize-1)%m = H(systemSize-1)%m + (-0.5_rKind*mu0 )*TensorProd(one_op%mr,ntot_opS%mr) +&
					VB*TensorProd(one_op%mr,VB_opS%mr)/2.0_rKind


!Construct the O operators to form spin-dependent part of the Hamiltonian
DO sTot=0, spinSize-1
	DO alpha=-sTot,sTot

	O_opS%mr=0.0_rKind
	DO alpha1=1,spinSize
		DO alpha2=1,spinSize
		O_opS%mr=O_opS%mr+Clebsch(FLOOR(2.0*spin),2*FLOOR(alpha1*1.0_rKind-spin-1.0_rKind),FLOOR(2.0*spin),2*FLOOR(alpha2*1.0_rKind-spin-1.0_rKind),2*sTot, 2*alpha)&
				*MATMUL(a_opS(alpha1)%mr,a_opS(alpha2)%mr)*SQRT(gS(sTot+1))
		END DO
	END DO

		DO i = 1, (systemSize-1)
		H(i)%m = H(i)%m+0.25_rKind*TensorProd(MATMUL(TRANSPOSE(O_opS%mr),O_opS%mr),one_op%mr)+&
				0.25_rKind*TensorProd(one_op%mr,MATMUL(TRANSPOSE(O_opS%mr),O_opS%mr))
		END DO
		
	H(1)%m = H(1)%m+0.25_rKind*TensorProd(MATMUL(TRANSPOSE(O_opS%mr),O_opS%mr),one_op%mr)
	H(systemSize-1)%m = H(systemSize-1)%m+0.25_rKind*TensorProd(one_op%mr,MATMUL(TRANSPOSE(O_opS%mr),O_opS%mr))



	END DO
END DO


DEALLOCATE(O_opS%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate O_opS'
			END IF		


END SUBROUTINE HamiltonianSpinS


SUBROUTINE SetupBHSpinOneName(baseName,diRectory)
!
!Purpose: Begin a file name in the directory diRectory that contains all spin One Bose-Hubbard parameters
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: diRectory

CALL createFileName(baseName,diRectory)
CALL appendBaseName(baseName,'BHSpinOne')
CALL appendBaseName(baseName,'L',systemSize)
IF(ncswitch) THEN
	CALL appendBaseName(baseName,'N',totNum)
ELSE
	CALL appendBaseName(baseName,'mu',2,mu0)
END IF
CALL appendBaseName(baseName,'Chi',chiMax)
CALL appendBaseName(baseName,'jTunn',2,jTunn)
CALL appendBaseName(baseName,'U0',2,U0)
CALL appendBaseName(baseName,'U2',2,U2)
CALL appendBaseName(baseName,'VB',2,VB)
baseName=TRIM(baseName)//TRIM(ADJUSTL(BoundaryCond))//'BC'

END SUBROUTINE SetupBHSpinOneName

END MODULE spinS_module
