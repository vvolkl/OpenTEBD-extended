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
MODULE rotation_module
!
! Purpose: Module to construct Molecular Hubbard operators/Hamiltonian
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

!Rotation system parameters
INTEGER, PARAMETER :: rotLevel=2 !highest rotational level kept
INTEGER :: rotSize=(1+rotLevel)*(1+rotLevel) !shorthand for size of M-degenerate subspace
INTEGER, PARAMETER :: Jcut=25 !Stark dressing cutoff
INTEGER, PARAMETER :: qDC=0 !spherical polarization of DC field 
INTEGER, PARAMETER :: qAC=0 !spherical polarization of AC field 
REAL(KIND=rKind), PARAMETER :: rotConst=1.0_rKind !B, the rotational constant
REAL(KIND=rKind), PARAMETER :: dip=(0.566_rKind)/2.54_rKind !permanent dipole, term in parentheses is value in Debye
REAL(KIND=rKind) :: eDC=0.0_rKind/(dip*rotConst) !Strength of DC field
REAL(KIND=rKind) :: eAC=5.0_rKind/(dip*rotConst) !10 times strength of AC field
REAL(KIND=rKind) :: omega=0.0_rKind !frequency of AC field (defined later in code)
REAL(KIND=rKind) :: detuning=0.0_rkind !detuning from resonance
REAL(KIND=rKind), PARAMETER :: alphaBar=237.0_rKind !Average polarizability
REAL(KIND=rKind), PARAMETER :: deltaAlpha=165.8_rKind !Polarizability anisotropy
REAL(KIND=rKind), PARAMETER :: Udipdip=10.0_rKind !dipole-dipole energy scale
REAL(KIND=rKind), PARAMETER :: ERecoil=10.0_rKind !Recoil energy
REAL(KIND=rKind) :: LattHeight !Ersatz height of lattice


! *** MHH operators ***
TYPE(matrixReal) :: EDC_opR !DC Field operator
TYPE(matrixReal) :: EAC_opR !AC field operator
TYPE(matrixReal) :: ttot_opR !Total tunneling operator
TYPE(matrixReal) :: dipdip_opR !dipole-dipole operator

!Stark dressing allocations
INTEGER :: LWORK
REAL(KIND=rKind), ALLOCATABLE :: W(:), Work(:)

CONTAINS

SUBROUTINE AllocateDsyev(dimen)
!
!Purpose: Allocate the necessary variables for the LAPACK diagonalization
!routine Dsyev
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: dimen
	LWORK=3*dimen+2
	ALLOCATE(W(dimen), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate W for Dsyev'
			END IF
	ALLOCATE(Work(LWORK), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Work for Dsyev'
			END IF	
END SUBROUTINE AllocateDsyev

SUBROUTINE DeAllocateDsyev()
!
!Purpose: Dellocate the necessary variables for the LAPACK diagonalization
!routine Dsyev
!
IMPLICIT NONE
	DEALLOCATE(W, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate W for Dsyev'
			END IF
	DEALLOCATE(Work, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Work for Dsyev'
			END IF	
END SUBROUTINE DeAllocateDsyev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!M=0 code begins!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


INTEGER FUNCTION rotationLocalDimMzero()
!
!Purpose: Calculate physical site local dimension for a system
!with at most maxFilling particles on site, each allowed to be in at most
!the rotLevel^th rotational level with M=0
!
IMPLICIT NONE
	rotationLocalDimMzero = BinomialCoef(maxFilling+1+rotLevel,maxFilling)
END FUNCTION rotationLocalDimMzero


SUBROUTINE DiagDsyev(dimen, Mat,OutMat,OutEig,dipoles)
!
!Purpose: Diagonalize the Rotational+DC Hamiltonian, return energies, eigenvectors, and dipoles
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: dimen
REAL(KIND=rKind), INTENT(OUT) :: dipoles(:,:)
REAL(KIND=rKind), INTENT(OUT) :: OutEig(:),OutMat(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: Mat(:,:)
REAL(KIND=rKind)  :: A(dimen,dimen) 
CHARACTER(1) :: JOBZ, UPLO
INTEGER :: INFO,i,j,k,l, tuhcount(11)

!Store Matrix in a dummy array
A=Mat

!Save only the dipole portion of the field free states
DO i=1,dimen,1
	DO j=1,dimen,1
		IF(i==j) THEN
		Mat(i,j)=0.0_rKind
		END IF
		
		IF(Mat(i,j).ne.0.0) THEN
		Mat(i,j)=-Mat(i,j)/eDC
		END IF

	END DO
END DO


		JOBZ='V'
		UPLO='U'
	CALL DSYEV(JOBZ,UPLO,dimen,A,dimen,W,WORK,LWORK,INFO)

!Add in the zero-particle state
		OutMat=0.0_rKind
		OutMat(1,1)=1.0_rKind
		OutEig(1)=0.0_rKind

!Write the eigenvalues/vectors into output
		DO i=1,localSize-1,1	
			OutEig(i+1)=W(i)
			DO j=1,localSize-1,1
			OutMat(i+1,j+1)=A(i,j)
			END DO
		END DO	


!Write the induced dipoles into output
	Mat=MATMUL(TRANSPOSE(A),MATMUL(Mat,A))
	DO i=1,localSize,1
	dipoles(i,1)=0.0_rKind
	dipoles(1,i)=0.0_rKind
	END DO
	DO i=1,localSize-1
		DO j=1,localSize-1,1
			IF(i==j) THEN
			dipoles(i+1,j+1)=Mat(i,j)
			ELSE
			dipoles(i+1,j+1)=ABS(Mat(i,j))
			END IF
		END DO
	END DO

END SUBROUTINE DiagDsyev		


SUBROUTINE CreateRotationopsMzero()
!
!Purpose: Create the operators needed to define and characterize the Molecular-Hubbard Hamiltonian, M=0
!
IMPLICIT NONE
INTEGER :: i, j, k,l, d, dsq, j1, j2, m1, m2, ecount, passSize
REAL(KIND=rKind) :: mi, mj, mis,mjs, ell, norm(localSize)
REAL(KIND=rKind) ::  TunnJ(1+rotLevel),TunnJalt(1+rotLevel)
TYPE(matrix) :: stateList,newState 
COMPLEX(KIND=rKind) :: preFactor, eye
TYPE(matrixReal) :: dipoles
!Stark dressing allocations
TYPE(vector) :: EigAlt
TYPE(matrixReal) :: Ealt,EigVecAlt

passSize=1+rotLevel
		d=localSize
		dsq = d*d
		eye = CMPLX(0.0, 1.0, KIND=rKind)

IF(maxfilling>=2) THEN
STOP "Filling>=2 not supported in rotational code!"
END IF

!For the rotationally indexed operators we must:
!1. Allocate the number of lists
		ALLOCATE(a_opS(1+rotLevel), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for RotationOps'
			END IF

!2. Allocate each list
		DO i=1,1+rotLevel,1
		ALLOCATE(a_opS(i)%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for RotationOps'
			END IF
		END DO

!Otherwise just allocate the matrix
		ALLOCATE(EDC_opR%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EDC for RotationOps'
			END IF
		ALLOCATE(EAC_opR%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EAC for RotationOps'
			END IF
		ALLOCATE(one_op%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for RotationOps'
			END IF
		ALLOCATE(dipoles%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate dipoles for RotationOps'
			END IF

		ALLOCATE(ttot_opR%mr(dsq, dsq), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ttot for RotationOps'
			END IF
		ALLOCATE(dipdip_opR%mr(dsq, dsq), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate dipdip for RotationOps'
			END IF

!State enumeration stuff
		ALLOCATE(stateList%m(d, 1+rotLevel), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate statList for RotationOps'
			END IF
		ALLOCATE(newState%m(1, 1+rotLevel), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate newState for RotationOps'
			END IF
!Dressing stuff
		ALLOCATE(Ealt%mr(1+Jcut, 1+Jcut), STAT=statInt)
			IF(statInt/=0) THEN
			WRITE(*,*) 'Failed to allocate Ealt for RotationOps'
			END IF
		ALLOCATE(EigVecalt%mr(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EigVecalt for RotationOps'
			END IF
		ALLOCATE(EigAlt%v(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EigAlta for RotationOps'
			END IF
!Number Conserving stuff
		ALLOCATE(Conserv%vi(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Conserv for RotationOps'
			END IF

!Create state list
		CALL onsiteStateListIdof(stateList%m,passSize)

!Define the unit operator
		one_op%mr=0.0_rKind
		DO j=1,d,1
			one_op%mr(j,j)=1.0_rKind
		END DO

!Define the destruction operator
		DO k=1,1+rotLevel,1
		a_opS(k)%mr = 0.0_rKind
			DO j = 1, d
			newState%m(1, :) = stateList%m(j, :)
			preFactor = SQRT(newState%m(1, k))
			newState%m(1, k) = newState%m(1, k) - 1.0_rKind
				DO i = 1, d		
					a_opS(k)%mr(i, j) = preFactor*kronDelta(stateList%m(i, :), &
								newState%m(1, :), 1+rotLevel)
				END DO			
			END DO

		END DO
		DEALLOCATE(stateList%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate statelist for RotationOps'
			END IF
		DEALLOCATE(newState%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate newState for RotationOps'
			END IF

!Initialize the Clebsch Calculators
	CALL SetupLogFac()

!Define DC field operator on the extended (highest rotational state=Jcut) manifold
	Ealt%mr=0.0_rKind
		DO i=0,Jcut
			DO k=0,Jcut
		IF(k.ge.i) THEN						    
		Ealt%mr(i+1,k+1)=Ealt%mr(i+1,k+1)-dip*(eDC)*Clebsch(2*k,0,2,2*qDC,2*i,0)*Clebsch(2*k,0,2,0,2*i,0)*SQRT(((2*k+1)*1.0_rKind)/((2*i+1)*1.0_rKind))	
		ELSE
		Ealt%mr(i+1,k+1)=Ealt%mr(k+1,i+1)
		END IF

	
!Add in the rotational energies
		IF(k==i) THEN
		Ealt%mr(i+1,k+1)=rotConst*i*(i+1)
		END IF
			END DO
		END DO


	dipoles%mr=0.0_rKind
	CALL AllocateDsyev(SIZE(Ealt%mr,1))

!Diagonalize the DC+rotational Hamiltonian, store the eigenvalues/eigenvectors
!of the lowest rotLevel+1 dressed levels (to provide appropriate AC dressing
!of the lowest 2) as well as the induced dipoles
	CALL DiagDsyev(SIZE(Ealt%mr,1), Ealt%mr,EigVecAlt%mr,EigAlt%v, dipoles%mr)

	CALL DeAllocateDsyev()

!We now have our DC dressed energies, we work in the basis where this is diagonal
	EDC_opR%mr=0.0_rKind
	DO i=1,1+rotLevel,1
		EDC_opR%mr=EDC_opR%mr+EigAlt%v(i+1)*MATMUL(TRANSPOSE(a_opS(i)%mr), a_opS(i)%mr)
	END DO




!Renormalize the reduced eigenspace (convergence plots in the NJP paper)
	DO i=1,localSize,1
	norm(i)=0.0_rKind
		DO j=1,localSize,1
		norm(i)=norm(i)+EigVecAlt%mr(j,i)*EigVecAlt%mr(j,i)
		END DO
	END DO

	DO i=1,localSize,1
		DO j=1,localSize,1	
		EigVecAlt%mr(j,i)=EigVecAlt%mr(j,i)/SQRT(norm(i))
		END DO
	END DO


!Define the AC perturbation Hamiltonian of 0.1 unit field strength using
!the induced dipoles from the DC dressing
		DO i=1,d,1
			DO j=1,d,1
				IF(j.ge.i) THEN	
				EAC_opR%mr(i,j)=-0.1*dipoles%mr(i,j)
				ELSE
				EAC_opR%mr(i,j)=EAC_opR%mr(j,i)
				END IF	
			END DO
		END DO	


!Define the dipole-dipole operator
	dipdip_opR%mr=0.0_rKind

	DO j1=2,d,1
		DO m1=2,d,1
			DO j2=2,d,1
				DO m2=2,d,1

		dipdip_opR%mr=dipdip_opR%mr+0.5_rKind*Udipdip*dipoles%mr(j1,m1)*dipoles%mr(j2,m2) &
		*TensorProd(MATMUL(TRANSPOSE(a_opS(j1-1)%mr),a_opS(m1-1)%mr),MATMUL(TRANSPOSE(a_opS(j2-1)%mr),a_opS(m2-1)%mr))

				END DO
			END DO
		END DO
	END DO	

!Define the field-free tunneling energy
	DO k=1,1+rotLevel,1
	TunnJ(k)=Erecoil*((LattHeight*(1.0_rKind+2.0_rKind*(deltaAlpha/alphaBar)*(REAL(k*(k-1),KIND=rKind)/(REAL((2*k-3)*(2*k+1),KIND=rKind)))))**(1.051))&
			*Exp(-2.121_rKind*SQRT(LattHeight*(1.0_rKind+2.0_rKind*(deltaAlpha/alphaBar)*(REAL(k*(k-1),KIND=rKind)/(REAL((2*k-3)*(2*k+1),KIND=rKind))))))
	IF(print_switch) THEN
	PRINT *, TunnJ(k)
	END IF
	END DO

!Print the effective tunneling strength
	DO i=1,1+rotLevel,1
	TunnJalt(i)=0.0_rKind
		DO j=1,1+rotLevel,1
	TunnJalt(i)=TunnJalt(i)+EigVecAlt%mr(j+1,i+1)*EigVecAlt%mr(j+1,i+1)*TunnJ(j)
		END DO
	IF(print_switch) THEN
		PRINT *,'T in the ',i-1,'th dressed level=', TunnJalt(i), 'Resonant dipole interaction in the ',i-1,'th level=',Udipdip*dipoles%mr(i+1,i+1)*dipoles%mr(i+1,i+1)
		PRINT *,'t/u in the ',i-1,'th dressed level=', TunnJalt(i)/(Udipdip*dipoles%mr(i+1,i+1)*dipoles%mr(i+1,i+1))
	END IF
	END DO

	DO i=1,3,1
	IF(print_switch) THEN
		PRINT *,'T in the ',i-1,'th dressed level=', TunnJalt(i), 'Resonant dipole interaction in the ',i-1,'th level=',Udipdip*dipoles%mr(i+1,i+1)*dipoles%mr(i+1,i+1)
		PRINT *,'t/u in the ',i-1,'th dressed level=', TunnJalt(i)/(Udipdip*dipoles%mr(i+1,i+1)*dipoles%mr(i+1,i+1))
	END IF
	END DO

	IF(print_switch) THEN
	PRINT *, 't/U ratios=', (TunnJalt(1)/(Udipdip*dipoles%mr(1+1,1+1)*dipoles%mr(1+1,1+1)))/(TunnJalt(2)/(Udipdip*dipoles%mr(2+1,2+1)*dipoles%mr(2+1,2+1)))
	END IF


!Transform tunneling operator to the eigenbasis of the DC dressed Hamiltonian
	ttot_opR%mr=0.0_rKind
	DO k=1,1+rotLevel,1
	ttot_opR%mr = ttot_opR%mr+TunnJ(k)*TensorProd(Transpose(MATMUL(TRANSPOSE(EigVecAlt%mr),MATMUL(a_opS(k)%mr,EigVecAlt%mr))),MATMUL(TRANSPOSE(EigVecAlt%mr),MATMUL(a_opS(k)%mr,EigVecAlt%mr))) &
		+ TunnJ(k)*Transpose(TensorProd(Transpose(MATMUL(TRANSPOSE(EigVecAlt%mr),MATMUL(a_opS(k)%mr,EigVecAlt%mr))),MATMUL(TRANSPOSE(EigVecAlt%mr),MATMUL(a_opS(k)%mr,EigVecAlt%mr))))
	END DO


!Define the number conservation vector
		Conserv%vi=0
		DO i=2,localSize,1
		Conserv%vi(i)=1
		END DO
	
		DEALLOCATE(Ealt%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Ealt for RotationOps'
			END IF
		DEALLOCATE(EigVecalt%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate EigVecalt for RotationOps'
			END IF
		DEALLOCATE(EigAlt%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate EigAlt for RotationOps'
			END IF
		DEALLOCATE(dipoles%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate dipoles for RotationOps'
			END IF

	IF(print_switch) THEN
	PRINT *, "M=0 Rotation-R operators created!"	
	END IF
END SUBROUTINE CreateRotationopsMzero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!M=0 code ends!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!M!=0 code begins!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION rotationLocalDim()
!
!Purpose: Calculate physical site local dimension for a system
!with at most maxFilling particles on site, each allowed to be in at most
!the rotLevel^th rotational level
!
	rotationLocalDim = BinomialCoef(maxFilling+rotSize,maxFilling)
END FUNCTION rotationLocalDim

SUBROUTINE CreateRotationops()
!
!Purpose: Create the operators needed to define and characterize the Molecular-Hubbard Hamiltonian, M=0
!
IMPLICIT NONE
INTEGER :: i, j, k,l, d, dsq, j1, j2, m1, m2, ecount, counter1, counter2, passSize
REAL(KIND=rKind) :: mi, mj, mis,mjs, ell, norm(localSize)
REAL(KIND=rKind) ::  TunnJ(rotSize),TunnJalt(rotSize)
TYPE(matrix) :: stateList,newState 
COMPLEX(KIND=rKind) :: preFactor, eye
TYPE(matrixReal) :: dipoles
!Stark dressing allocations
TYPE(vector) :: EigAlt
TYPE(matrixReal) :: Ealt,EigVecAlt

passSize=rotSize

		d=localSize
		dsq = d*d
		eye = CMPLX(0.0, 1.0, KIND=rKind)

IF(maxfilling>=2) THEN
STOP "Filling>=2 not supported in rotational code!"
END IF

!For the rotationally indexed operators we must:
!1. Allocate the number of lists
		ALLOCATE(a_opS(rotSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for RotationOps'
			END IF

!2. Allocate each list
		DO i=1,rotSize,1
		ALLOCATE(a_opS(i)%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for RotationOps'
			END IF
		END DO

!Otherwise just allocate the matrix
		ALLOCATE(EDC_opR%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EDC for RotationOps'
			END IF
		ALLOCATE(EAC_opR%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EAC for RotationOps'
			END IF
		ALLOCATE(one_op%mr(d, d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate a for RotationOps'
			END IF
		ALLOCATE(dipoles%mr(d,d), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate dipoles for RotationOps'
			END IF

		ALLOCATE(ttot_opR%mr(dsq, dsq), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ttot for RotationOps'
			END IF
		ALLOCATE(dipdip_opR%mr(dsq, dsq), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate dipdip for RotationOps'
			END IF

!State enumeration stuff
		ALLOCATE(stateList%m(d, rotSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate statList for RotationOps'
			END IF
		ALLOCATE(newState%m(1,rotSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate newState for RotationOps'
			END IF
!Dressing stuff
		ALLOCATE(Ealt%mr(Floor(BinomialCoef(maxFilling+(1+Jcut)*(1+Jcut),maxFilling)),Floor(BinomialCoef(maxFilling+(1+Jcut)*(1+Jcut),maxFilling))), STAT=statInt)
			IF(statInt/=0) THEN
			WRITE(*,*) 'Failed to allocate Ealt for RotationOps'
			END IF
		ALLOCATE(EigVecalt%mr(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EigVecalt for RotationOps'
			END IF
		ALLOCATE(EigAlt%v(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate EigAlta for RotationOps'
			END IF
!Number Conserving stuff
		ALLOCATE(Conserv%vi(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Conserv for RotationOps'
			END IF




!Create state list
		CALL onsiteStateListIdof(stateList%m, passSize)

!Define the unit operator
		one_op%mr=0.0_rKind
		DO j=1,d,1
			one_op%mr(j,j)=1.0_rKind
		END DO

!Define the destruction operator
		DO k=1,rotSize,1
		a_opS(k)%mr = 0.0_rKind
			DO j = 1, d

			newState%m(1, :) = stateList%m(j, :)
			preFactor = SQRT(newState%m(1, k))
			newState%m(1, k) = newState%m(1, k) - 1.0_rKind
				DO i = 1, d	

					a_opS(k)%mr(i, j) = preFactor*kronDelta(stateList%m(i, :), &
								newState%m(1, :), rotSize)

				END DO			
			END DO
		END DO


		DEALLOCATE(stateList%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate statelist for RotationOps'
			END IF
		DEALLOCATE(newState%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate newState for RotationOps'
			END IF

!Initialize the Clebsch Calculators
	CALL SetupLogFac()

!Define DC+rotational field operator on the extended (highest rotational state=Jcut) manifold
	counter1=1
	counter2=1
	Ealt%mr=0.0_rKind
	DO i=0,Jcut
		DO j=i,-i,-1
		counter1=counter1+1
			DO k=0,Jcut
				DO l=k,-k,-1
				counter2=k*k+k-l+2
					IF(counter2.ge.counter1) THEN						    
					Ealt%mr(counter1,counter2)=Ealt%mr(counter1,counter2)-dip*eDC*Clebsch(2*k,2*l,2,2*qDC,2*i,2*j)*Clebsch(2*k,0,2,0,2*i,0)*SQRT(((2*k+1)*1.0_rKind)/((2*i+1)*1.0_rKind))		
					ELSE
					Ealt%mr(counter1,counter2)=Ealt%mr(counter2,counter1)
					END IF	
					!Add in the rotational energies
					IF(counter2==counter1) THEN
					Ealt%mr(counter1,counter2)=rotConst*i*(i+1)
					END IF
			END DO
		END DO
	END DO
END DO	
		


	dipoles%mr=0.0_rKind

	CALL AllocateDsyev(SIZE(Ealt%mr,1))

!Diagonalize the DC+rotational Hamiltonian, store the eigenvalues/eigenvectors
!of the lowest rotLevel+1 dressed levels (to provide appropriate AC dressing
!of the lowest 2) as well as the induced dipoles
	CALL DiagDsyev(SIZE(Ealt%mr,1), Ealt%mr,EigVecAlt%mr,EigAlt%v, dipoles%mr)

	CALL DeAllocateDsyev()

!We now have our DC dressed energies, we work in the basis where this is diagonal
	EDC_opR%mr=0.0_rKind
	DO i=1,rotSize,1
		EDC_opR%mr=EDC_opR%mr+EigAlt%v(i+1)*MATMUL(TRANSPOSE(a_opS(i)%mr), a_opS(i)%mr)
	END DO


!Renormalize the reduced eigenspace (convergence plots in the documentation)
	DO i=1,localSize,1
	norm(i)=0.0_rKind
		DO j=1,localSize,1
		norm(i)=norm(i)+EigVecAlt%mr(j,i)*EigVecAlt%mr(j,i)
		END DO
	END DO


	DO i=1,localSize,1
		DO j=1,localSize,1	
		EigVecAlt%mr(j,i)=EigVecAlt%mr(j,i)/SQRT(norm(i))
		END DO
	END DO

!Define the AC perturbation Hamiltonian of 0.1 unit field strength using
!the induced dipoles from the DC dressing
		DO i=1,d,1
			DO j=1,d,1
				IF(j.ge.i) THEN	
				EAC_opR%mr(i,j)=-0.1*dipoles%mr(i,j)
				ELSE
				EAC_opR%mr(i,j)=EAC_opR%mr(j,i)
				END IF	
			END DO
		END DO	




!!!!!!!!!!!!!!!!!!!!!!DONE TO HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
		DEALLOCATE(Ealt%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Ealt for RotationOps'
			END IF
		DEALLOCATE(EigVecalt%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate EigVecalt for RotationOps'
			END IF
		DEALLOCATE(EigAlt%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate EigAlt for RotationOps'
			END IF
		DEALLOCATE(dipoles%mr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate dipoles for RotationOps'
			END IF

	IF(print_switch) THEN
	PRINT *, "Rotation-R operators created!"	
	END IF
END SUBROUTINE CreateRotationops

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!M!=0 code ends!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE HamiltonianRotationTI(H)
!
!Purpose: Construct the TEBD form of Rotational Hamiltonian: Time-independent version for ITP
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
INTEGER :: i

DO i = 1, (systemSize-1)
	H(i)%m =-ttot_opR%mr+dipdip_opR%mr+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr)
END DO

	H(1)%m = H(1)%m+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)
	H(systemSize-1)%m = H(systemSize-1)%m+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr)
END SUBROUTINE HamiltonianRotationTI


SUBROUTINE HamiltonianRotationTD(H,time)
!
!Purpose: Construct the TEBD form of Rotational Hamiltonian: Time-dependent version for RTP
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) ::   time
REAL (KIND=rKind) :: ACfac
INTEGER :: i
	ACfac=2.0_rKind*eAC*sin(omega*time)

	DO i = 1, (systemSize-1)
		H(i)%m =-ttot_opR%mr+dipdip_opR%mr+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr) &
		+0.5_rKind*TensorProd(ACfac*EAC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(one_op%mr,ACfac*EAC_opR%mr)
	END DO

	H(1)%m = H(1)%m+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(ACfac*EAC_opR%mr,one_op%mr)
	H(systemSize-1)%m = H(systemSize-1)%m+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr)+0.5_rKind*TensorProd(one_op%mr,ACfac*EAC_opR%mr)
END SUBROUTINE HamiltonianRotationTD



SUBROUTINE DestroyRotationops()
!
!Purpose: Deallocate the operators needed to define and characterize the Molecular-Hubbard Hamiltonian
!
IMPLICIT NONE
INTEGER ::  k
		
!Deallocate scalar operators
DEALLOCATE(ttot_opR%mr,one_op%mr, STAT=statInt)
	IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate RotationOps'
	END IF
!Deallocate vector operators
DO k=1,1+rotLevel,1
	DEALLOCATE(a_opS(k)%mr, STAT=statInt)
		IF(statInt.ne.0) THEN
		PRINT *, 'Failed to deallocate RotationOps'
		END IF
END DO

DEALLOCATE(a_opS, STAT=statInt)
	IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate RotationOps'
	END IF
DEALLOCATE(Conserv%vi, STAT=statInt)
	IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate RotationOps'
	END IF
DEALLOCATE(EDC_opR%mr, STAT=statInt)
	IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate RotationOps'
	END IF
DEALLOCATE(EAC_opR%mr, STAT=statInt)
	IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate RotationOps'
	END IF
DEALLOCATE(dipdip_opR%mr, STAT=statInt)
	IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate RotationOps'
	END IF
	IF(print_switch) THEN
	PRINT *, "Rotation-R operators destroyed!"
	END IF
END SUBROUTINE DestroyRotationops

SUBROUTINE SetupRotName(baseName,diRectory)
!
!Purpose: Begin a file name in the directory diRectory that contains all MHH parameters
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: diRectory

CALL createFileName(baseName,diRectory)
CALL appendBaseName(baseName,'Rot')
IF(ncswitch) THEN
	CALL appendBaseName(baseName,'N',totNum)
ELSE
	CALL appendBaseName(baseName,'mu',2,mu0)
END IF
CALL appendBaseName(baseName,'L',systemSize)
CALL appendBaseName(baseName,'EDC',2,eDC)
CALL appendBaseName(baseName,'LH',2,LattHeight)
CALL appendBaseName(baseName,'Chi',chiMax)

END SUBROUTINE SetupRotName

END MODULE rotation_module
