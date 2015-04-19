!    Copyright (C)  2009  M. L. Wall
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
MODULE fermi_hubbard_module
!
! Purpose: Module to construct spin-S Fermi-Hubbard (BH) model operators/Hamiltonian
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

CONTAINS

INTEGER FUNCTION spinSFermiLocalDim()
!
!Purpose: Calculate local dimension for Fermionic spin-S system.
!
IMPLICIT NONE
INTEGER :: i
IF(maxFilling.gt.spinSize) THEN
PRINT *, 'Warning: maxfilling is larger than allowed for a spin-s fermionic system!'
PRINT *,  'It is being set to the largest possible value=',spinSize
maxFilling=spinSize
END IF

spinSFermiLocalDim=0
DO i=0,maxFilling
	spinSFermiLocalDim=spinSFermiLocalDim+FLOOR(BinomialCoef(spinSize,i))
END DO

END FUNCTION spinSFermiLocalDim

SUBROUTINE CreateFermiSops()
!
!Purpose: Create the operators needed to define
!and characterize the spin-s Fermion model
!
IMPLICIT NONE

INTEGER :: i, j, k,l, d, dsq,fPhase, fermiSwitch
REAL(KIND=rKind) :: mi, mj
TYPE(matrix) :: stateList,newState
TYPE(matrixReal) :: Basis
COMPLEX(KIND=rKind) :: preFactor, eye

		d=localSize
		dsq = d*d
		eye = CMPLX(0.0, 1.0, KIND=rKind)


!For the spin indexed operators we must:
!1. Allocate the number of lists
		ALLOCATE(a_opS(spinSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for FermiSops'
			END IF

!2. Allocate each list
		DO i=1,spinSize,1
		ALLOCATE(a_opS(i)%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for FermiSops'
			END IF
		END DO

!Otherwise just allocate the operator
		ALLOCATE(one_op%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate one for FermiSops'
			END IF
		ALLOCATE(ntot_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ntot for FermiSops'
			END IF
		ALLOCATE(Sx_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sx for FermiSops'
			END IF
		ALLOCATE(Sy_opS%m(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sy for FermiSops'
			END IF
		ALLOCATE(Sz_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Sz for FermiSops'
			END IF
		ALLOCATE(Ssq_opS%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Ssq for FermiSops'
			END IF
		ALLOCATE(ttot_opS%mr(dsq, dsq), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ttot for FermiSops'
			END IF
		ALLOCATE(fermiPhase_op%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for FermiSops'
			END IF


!State enumeration stuff
		ALLOCATE(stateList%m(d, spinSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate stateList for FermiSops'
			END IF
		ALLOCATE(newState%m(1, spinSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate newState for FermiSops'
			END IF
		ALLOCATE(Conserv%vi(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Conserv for FermiSops'
			END IF

!Define the unity operator
		one_op%mr=0.0_rKind
		DO j=1,d,1
			one_op%mr(j,j)=1.0_rKind
		END DO

fermiSwitch=1
!Create and index the Fock space
		CALL onsiteStateListIdof(stateList%m, spinSize, fermiSwitch)

!Define "bosonic" (phaseless) operators to create Spin operators in the Schwinger representation
		DO k=1,spinSize,1
		a_opS(k)%mr = 0.0_rKind
			DO j = 1, d
			newState%m(1, :) = stateList%m(j, :)
			newState%m(1, k) = newState%m(1, k) - 1.0_rKind
				DO i = 1, d		
					a_opS(k)%mr(i, j) = kronDelta(stateList%m(i, :), &
								newState%m(1, :), spinSize)
				END DO			
			END DO
		END DO

Sx_opS%mr=0.0_rKind
Sy_opS%m=0.0_rKind
Sz_opS%mr=0.0_rKind

!Define spin operators
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

		END DO
!Define the total spin operator
		Ssq_opS%mr=MATMUL(Sx_opS%mr,Sx_opS%mr)+REAL(MATMUL(Sy_opS%m,Sy_opS%m)) &
						+MATMUL(Sz_opS%mr,Sz_opS%mr)

!Define the destruction operator with the fermi phase
	DO k=1,spinSize,1
	a_opS(k)%mr = 0.0_rKind
		DO j = 1, d
		newState%m(1, :) = stateList%m(j, :)
		newState%m(1, k) = newState%m(1, k) - 1.0_rKind
			DO i = 1, d
			fPhase=0
				DO l=1,k-1
				fPhase=fPhase+FLOOR(REAL(stateList%m(j, l)))
				END DO		
		a_opS(k)%mr(i, j) = ((-1.0_rKind)**fPhase)*kronDelta(stateList%m(i, :), &
							newState%m(1, :), spinSize)
			END DO			
		END DO
	END DO

!Define the unity operator
		one_op%mr=0.0_rKind
		DO j=1,d,1
			one_op%mr(j,j)=1.0_rKind
		END DO

!Deallocate the state enumeration stuff
		DEALLOCATE(stateList%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate stateList for FermiSops'
			END IF
		DEALLOCATE(newState%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate newState for FermiSops'
			END IF

ntot_opS%mr=0.0_rKind
!Define the total number operator
	DO k=1,spinSize,1
	ntot_opS%mr=ntot_opS%mr+MATMUL(TRANSPOSE(a_opS(k)%mr), a_opS(k)%mr)	
	END DO

!Construct the fermi phase operator to define phase for n-site observables
	fermiPhase_op%mr=0.0_rKind
	DO i=1,localSize
	fermiPhase_op%mr(i,i)=(-1.0_rKIND)**(ntot_opS%mr(i,i))
	END DO

!Store the on-site quantum numbers
		Conserv%vi=0
	DO i=1,localSize,1
		Conserv%vi(i)=FLOOR(ntot_opS%mr(i,i))
	END DO

!Define the tunneling operator-note the fermi phase
	ttot_opS%mr=0.0
	DO k=1,spinSize,1
	ttot_opS%mr=ttot_opS%mr+MATMUL(TRANSPOSE(TensorProd(a_opS(k)%mr,one_op%mr)),TensorProd(fermiPhase_op%mr,a_opS(k)%mr))&
						+MATMUL(TRANSPOSE(TensorProd(fermiPhase_op%mr,a_opS(k)%mr)),TensorProd(a_opS(k)%mr,one_op%mr))
	
	END DO
		
	IF(print_switch) THEN		
	PRINT *, 'Spin-s Fermi operators created!'
	END IF
END SUBROUTINE CreateFermiSops


SUBROUTINE DestroyFermiSops()
!
!Purpose: Deallocate spin-S operators.
!
IMPLICIT NONE
INTEGER ::  k
		
!Deallocate scalar operators
		DEALLOCATE(Sx_opS%mr, Sy_opS%m, Sz_opS%mr, Ssq_opS%mr, ntot_opS%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate FermiSops'
			END IF
		DEALLOCATE(ttot_opS%mr,one_op%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate FermiSops'
			END IF
		DEALLOCATE(Conserv%vi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate FermiSops'
			END IF
		DEALLOCATE(fermiPhase_op%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate FermiSops'
			END IF
!Deallocate vector operators
		DO k=1,spinSize,1
		DEALLOCATE(a_opS(k)%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate FermiSops'
			END IF
		END DO
		DEALLOCATE(a_opS, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate FermiSops'
			END IF
		IF(print_switch) THEN
		PRINT *, "Spin-s Fermi operators destroyed!"
		END IF
END SUBROUTINE DestroyFermiSops

SUBROUTINE HamiltonianHubbard(H ,J,U0, mu0,V0, extPot)
!
!Purpose: Create TEBD form of Hubbard model.
!
! See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) :: J, mu0, U0
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: extPot(:), V0
INTEGER :: i,k 

IF(spin.gt.0.5_rKind) THEN
PRINT *, 'Warning!  Hubbard Hamiltonian has interactions for only the lowest two components!'
ELSE IF(spin==0.0) THEN
PRINT *, 'Spinless fermions have no on-site interactions!  U0 is being ignored!'
END IF

!Spin-independent part of Hamiltonian		
DO i = 1, (systemSize-1)

	H(i)%m = (-0.5_rKind*mu0 )*TensorProd(ntot_opS%mr,one_op%mr) +&
	(-0.5_rKind*mu0 )*TensorProd(one_op%mr,ntot_opS%mr) &
		- J*ttot_opS%mr
	IF(spin.ge.0.5_rKind) THEN
	H(i)%m =H(i)%m	+0.5_rKind*U0*TensorProd(&
	MATMUL(MATMUL(TRANSPOSE(a_opS(1)%mr), a_opS(1)%mr)-0.5_rKind*one_op%mr,MATMUL(TRANSPOSE(a_opS(2)%mr), a_opS(2)%mr)-0.5_rKind*one_op%mr)&
	,one_op%mr)+0.5_rKind*U0*TensorProd(one_op%mr,&
	MATMUL(MATMUL(TRANSPOSE(a_opS(1)%mr), a_opS(1)%mr)-0.5_rKind*one_op%mr,MATMUL(TRANSPOSE(a_opS(2)%mr), a_opS(2)%mr)-0.5_rKind*one_op%mr))
	END IF		
END DO

H(1)%m = H(1)%m + (-0.5_rKind*mu0 )*TensorProd(ntot_opS%mr,one_op%mr)
	IF(spin.ge.0.5_rKind) THEN
	H(1)%m = H(1)%m +0.5_rKind*U0*TensorProd(&
	MATMUL(MATMUL(TRANSPOSE(a_opS(1)%mr), a_opS(1)%mr)-0.5_rKind*one_op%mr,MATMUL(TRANSPOSE(a_opS(2)%mr), a_opS(2)%mr)-0.5_rKind*one_op%mr)&
	,one_op%mr)
	END IF
H(systemSize-1)%m = H(systemSize-1)%m + (-0.5_rKind*mu0 )*TensorProd(one_op%mr,ntot_opS%mr) 
	IF(spin.ge.0.5_rKind) THEN
	H(systemSize-1)%m = H(systemSize-1)%m+0.5_rKind*U0*TensorProd(one_op%mr,&
	MATMUL(MATMUL(TRANSPOSE(a_opS(1)%mr), a_opS(1)%mr)-0.5_rKind*one_op%mr,MATMUL(TRANSPOSE(a_opS(2)%mr), a_opS(2)%mr)-0.5_rKind*one_op%mr))
	END IF
	
IF(PRESENT(extPot)) THEN
DO k=1,spinSize
DO i=1,(systemSize-1)
H(i)%m=H(i)%m+0.5_rKind*extPot(i)*tensorProd(MATMUL(TRANSPOSE(a_opS(k)%mr),a_opS(k)%mr),one_op%mr)+0.5_rKind*extPot(i+1)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_opS(k)%mr),a_opS(k)%mr))
END DO
H(systemSize-1)%m=H(systemSize-1)%m+0.5_rKind*extPot(systemSize)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_opS(k)%mr),a_opS(k)%mr))
H(1)%m=H(1)%m+0.5_rKind*extPot(1)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_opS(k)%mr),a_opS(k)%mr))
END DO
END IF

IF(PRESENT(V0)) THEN
	DO k=1,spinSize
		DO i=1,(systemSize-1)
			H(i)%m=H(i)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_opS(k)%mr),a_opS(k)%mr),MATMUL(TRANSPOSE(a_opS(k)%mr),a_opS(k)%mr))
		END DO
	END DO
END IF


IF(BoundaryCond=='P') THEN
	STOP "Periodic Boundary conditions not supported for Fermionic systems!"
END IF

END SUBROUTINE HamiltonianHubbard

SUBROUTINE SetupFHName(baseName,diRectory)
!
!Purpose: Begin a file name in the directory diRectory that contains all Fermi Hubbard parameters
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: diRectory

CALL createFileName(baseName,diRectory)
CALL appendBaseName(baseName,'FH_')
CALL appendBaseName(baseName,'L',systemSize)
IF(ncswitch) THEN
	CALL appendBaseName(baseName,'N',totNum)
ELSE
	CALL appendBaseName(baseName,'mu',2,mu0)
END IF
CALL appendBaseName(baseName,'Chi',chiMax)
CALL appendBaseName(baseName,'jTunn',2,jTunn)
CALL appendBaseName(baseName,'U0',2,U0)
IF(V0.ne.0.0) THEN
CALL appendBaseName(baseName,'V0',2,V0)
END IF
IF(ALLOCATED(extPot)) THEN
CALL appendBaseName(baseName,'extPot')
END IF

END SUBROUTINE SetupFHName

END MODULE fermi_hubbard_module
