!    Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009  J. Williams, I. Danshita, R. Mishmash, D. Schirmer,B. Schneider,  M. L. Wall
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
MODULE TEBDtools_module
!
! Purpose: Module Containing derived types/allocations/matrix manipulations
! for OpenSourceTEBD v1.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   9/28/09   M. L. Wall	v2.0 release
!
USE system_parameters
IMPLICIT NONE

! *** DERIVED TYPES ***

TYPE vector
	REAL(KIND=rKind), POINTER :: v(:)
END TYPE vector

TYPE vectorComplex
	COMPLEX(KIND=rKind), POINTER :: vc(:)
END TYPE vectorComplex

TYPE matrix
	COMPLEX(KIND=rKind), POINTER :: m(:,:)
END TYPE matrix

TYPE tensor
	COMPLEX(KIND=rKind), POINTER :: t(:,:,:)
END TYPE tensor
	
TYPE vectorInt
	INTEGER, POINTER :: vi(:)
END TYPE vectorInt
	
TYPE matrixInt
	INTEGER, POINTER :: mi(:,:)
END TYPE matrixInt
	
TYPE matrixReal
	REAL(KIND=rKind), POINTER :: mr(:,:)
END TYPE matrixReal

!Type for local measures
TYPE mlocal
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	COMPLEX(KIND=rKIND), POINTER :: value(:)
END TYPE mlocal

!Type for site-averaged measures
TYPE mavg
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	COMPLEX(KIND=rKIND) :: value
END TYPE mavg

!Type for correlation functions
TYPE mcorr
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	COMPLEX(KIND=rKIND), POINTER :: value(:,:)
END TYPE mcorr

!Type for correlation functions with fermi phases
TYPE mcorrf
	COMPLEX(KIND=rKIND), POINTER :: Op(:,:)
	COMPLEX(KIND=rKIND), POINTER :: value(:,:)
END TYPE mcorrf

TYPE entropy
	REAL(KIND=rKIND) :: qme
	REAL(KIND=rKIND), POINTER :: vN(:)
	REAL(KIND=rKIND), POINTER :: chain(:)
	REAL(KIND=rKIND), POINTER :: tbvN(:,:)
END TYPE entropy

!*** INTERFACES
!Matrix exponential of real/complex matrices
INTERFACE matrix_exponential
MODULE PROCEDURE matrix_exponential_r,&
				 matrix_exponential_c
END INTERFACE  matrix_exponential

INTERFACE ConstructPropagators
	MODULE PROCEDURE OLDConstructPropagators, NEWConstructPropagators
END INTERFACE ConstructPropagators

!Function/Subroutine tensor product of rr/cc/rc/cr matrices
INTERFACE tensorProd
MODULE PROCEDURE tensorProd_r,&
				 tensorProd_c,&
				 tensorProd_rc,&
				 tensorProd_cr
END INTERFACE  tensorProd

!kronecker delta of real/complex vectors
INTERFACE kronDelta
MODULE PROCEDURE kronDelta_r,&
				 kronDelta_c
END INTERFACE  kronDelta

!Trace of A*B AB=rr/cc/rc/cr
INTERFACE TraceMatmul
MODULE PROCEDURE TraceMatmul_rf,&
				 TraceMatmul_cf,&
				 TraceMatmul_rcf,&
				 TraceMatmul_crf
END INTERFACE  TraceMatmul
	
CONTAINS


!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE matrix_exponential!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE matrix_exponential_r(A,Exp_A,tau,n)
!
!Purpose: If matrix_exponential is called with argument types A=real, Exp_A=real, tau=Real
!         Then compute Exp_A=EXP(-tau*A)
!
!Based on routines by schneider, b. i.(nsf)
!
IMPLICIT NONE
INTEGER                                :: n, i, j, k, info
REAL(KIND=rKind), DIMENSION(:,:)                 :: A
REAL(KIND=rKind), DIMENSION(:,:)                 :: Exp_A
REAL(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
REAL(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
REAL(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Temp
REAL(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Scratch
REAL(KIND=rKind)                                 :: tau
REAL(KIND=rKind)								 :: expeig
CHARACTER (LEN=80)                     :: title
! Allocate storage for diagonalization routine.
ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Scratch(10*n), Temp(n,n) , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate matrix_exp variables'
			END IF 
Eigen_Vectors = A
!Call LAPACK routine to diagonalize double precision real symmetric matrix
CALL DSYEV('v','l',n,Eigen_Vectors,n,Eigen_Values,           &
              Scratch,10*n,info)
			  
! Form the matrix with exponentials of the eigenvalues on the diagonal
! Then similarity transform back into the original basis
DO i=1,n
  	DO j=1,n
	Exp_A(i,j)=0.0_rKind
		DO k=1,n
  	  	expeig=exp(-tau*Eigen_Values(k))
        Exp_A(i,j) = Exp_A(i,j) + Eigen_Vectors(i,k)*expeig*Eigen_Vectors(j,k)
		END DO
    END DO
END DO

! Deallocate the unneeded storage
DEALLOCATE ( Eigen_Values, Eigen_Vectors, Scratch, Temp  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate matrix_exp variables'
			END IF   
END SUBROUTINE matrix_exponential_r


SUBROUTINE matrix_exponential_c(A,Exp_A,t,n)
!
!Purpose: If matrix_exponential is called with argument types A=complex, Exp_A=complex, t=complex
!         Then compute Exp_A=EXP(-i*t*A)
!
!Based on routines by schneider, b. i.(nsf)
!
IMPLICIT NONE
INTEGER                                    :: n, i, j, k, info
COMPLEX(KIND=rKind), DIMENSION(:,:)                 :: A
COMPLEX(KIND=rKind), DIMENSION(:,:)                 :: Exp_A
COMPLEX(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Eigen_Vectors
REAL(KIND=rKind),     DIMENSION(:),   ALLOCATABLE    :: Eigen_Values
COMPLEX(KIND=rKind), DIMENSION(:),   ALLOCATABLE    :: Workv
COMPLEX(KIND=rKind), DIMENSION(:,:), ALLOCATABLE    :: Temp
REAL(KIND=rKind),     DIMENSION(:),   ALLOCATABLE    :: Rworkv
COMPLEX(KIND=rKind)                                 :: t
CHARACTER (LEN=80)                         :: title
COMPLEX(KIND=rKind)                                 :: eye=(0.0_rKind,1.0_rKind)
COMPLEX(KIND=rKind)								 :: expeig
! Allocate some storage for diagonalization routine.
ALLOCATE ( Eigen_Values(n), Eigen_Vectors(n,n), Workv(10*n), Rworkv(10*n), Temp(n,n)  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate matrix_exp variables'
			END IF   
Eigen_Vectors = A
!Call LAPACK routine to diagonalize double precision hermitian matrix
CALL ZHEEV('v','l',n,Eigen_Vectors,n,Eigen_Values,              &
              Workv,10*n,Rworkv,info)

! Form the matrix with exponentials of the eigenvalues on the diagonal
! Then similarity transform back into the original basis
DO i=1,n
	DO j=1,n
	Exp_A(i,j)=CMPLX(0.0,KIND=rKind)
		DO k=1,n
  	  	expeig=exp(-eye*t*Eigen_Values(k))
        Exp_A(i,j) = Exp_A(i,j) + Eigen_Vectors(i,k)*expeig*CONJG(Eigen_Vectors(j,k))
		END DO
    END DO
END DO
  ! Deallocate the unneeded storage
DEALLOCATE ( Eigen_Values, Eigen_Vectors, Workv, Rworkv, Temp  , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate matrix_exp variables'
			END IF     
END SUBROUTINE matrix_exponential_c


!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE matrix_exponential!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE tensorProd!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION tensorProd_r(A,B)
!
!Purpose: Return the tensor product of real matrices A and B
!    
IMPLICIT NONE     
REAL(KIND=rKind), INTENT(IN) :: A(:,:), B(:,:)
REAL(KIND=rKind) :: tensorProd_r(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_r((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_r

FUNCTION tensorProd_c(A,B)
!
!Purpose: Return the tensor product of complex matrices A and B
!  
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: A(:,:), B(:,:)
COMPLEX(KIND=rKind) :: tensorProd_c(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_c((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_c

FUNCTION tensorProd_rc(A,B)
!
!Purpose: Return the tensor product of real matrix A and complex matrix B
!  
IMPLICIT NONE
REAL(KIND=rkind), INTENT(IN) :: A(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: B(:,:)
COMPLEX(KIND=rKind) :: tensorProd_rc(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_rc((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_rc

FUNCTION tensorProd_cr(A,B)
!
!Purpose: Return the tensor product of complex matrix A and real matrix B
!  
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: A(:,:)
REAL(KIND=rkind), INTENT(IN) :: B(:,:)
COMPLEX(KIND=rKind) :: tensorProd_cr(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))
INTEGER i,j,k,l, dA1, dA2, dB1, dB2
dA1 = SIZE(A,1)
dA2 = SIZE(A,2)
dB1 = SIZE(B,1)
dB2 = SIZE(B,2)
	DO i=1,dA1
		DO j=1,dB1
			DO k=1,dA2
				DO l=1,dB2
					tensorProd_cr((i-1)*dB1+j,(k-1)*dB2+l)=A(i,k)*B(j,l)
				END DO
			END DO
		END DO
	END DO
END FUNCTION tensorProd_cr

!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE tensorProd!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE kronDelta !!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION kronDelta_r(vec1, vec2, dim)
!
!Purpose: Kronecker delta function defined for two real vectors vec1 and vec2.
!
IMPLICIT NONE	
REAL(KIND=rKind), INTENT(IN) :: vec1(:), vec2(:)
INTEGER, INTENT(IN) :: dim
INTEGER :: dim1, dim2, i, j
INTEGER :: booles
	dim1 = SIZE(vec1)
	dim2 = SIZE(vec2)
		IF (dim1 /= dim .OR. dim2 /= dim) THEN
			STOP "Dimensions of input vectors in function kronDelta must be the same."
		END IF
		DO i = 1, dim
			IF (vec1(i) == vec2(i)) THEN
				booles = 1
			ELSE
				booles = 0
				EXIT
			END IF
		END DO
	kronDelta_r = booles
END FUNCTION kronDelta_r	

INTEGER FUNCTION kronDelta_c(vec1, vec2, dim)
!
!Purpose: Kronecker delta function defined for two complex vectors vec1 and vec2.
!
IMPLICIT NONE	
COMPLEX(KIND=rKind), INTENT(IN) :: vec1(:), vec2(:)
INTEGER, INTENT(IN) :: dim
INTEGER :: dim1, dim2, i, j
INTEGER :: booles
	dim1 = SIZE(vec1)
	dim2 = SIZE(vec2)
		IF (dim1 /= dim .OR. dim2 /= dim) THEN
			STOP "Dimensions of input vectors in function kronDelta must be the same."
		END IF
		DO i = 1, dim
			IF (vec1(i) == vec2(i)) THEN
				booles = 1
			ELSE
				booles = 0
				EXIT
			END IF
		END DO
	kronDelta_c = booles
END FUNCTION kronDelta_c	
!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE kronDelta !!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE TraceMatmul !!!!!!!!!!!!!!!!!!!!!!!!


REAL(KIND=rKind) FUNCTION TraceMatmul_rf(A,B)
!
!Purpose: Function to calculate the trace of real matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: A(:,:), B(:,:)
INTEGER i, j
	TraceMatmul_rf=0.0_rKind
	DO i=1,SIZE(A,1)
		DO j=1,SIZE(A,2)
			TraceMatmul_rf = TraceMatmul_rf + A(i,j)*B(j,i)
		END DO
	END DO

END FUNCTION TraceMatmul_rf

COMPLEX(KIND=rKind) FUNCTION TraceMatmul_cf(A,B)
!
!Purpose: Function to calculate the trace of complex matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: A(:,:), B(:,:)
INTEGER i, j
	TraceMatmul_cf=CMPLX(0.0,KIND=rKind)
	DO i=1,SIZE(A,1)
		DO j=1,SIZE(A,2)
			TraceMatmul_cf = TraceMatmul_cf + A(i,j)*B(j,i)
		END DO
	END DO

END FUNCTION TraceMatmul_cf

COMPLEX(KIND=rKind) FUNCTION TraceMatmul_rcf(A,B)
!
!Purpose: Function to calculate the trace of real/complex matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: A(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: B(:,:)
INTEGER i, j
	TraceMatmul_rcf=CMPLX(0.0,KIND=rKind)
	DO i=1,SIZE(A,1)
		DO j=1,SIZE(A,2)
			TraceMatmul_rcf = TraceMatmul_rcf + A(i,j)*B(j,i)
		END DO
	END DO

END FUNCTION TraceMatmul_rcf

COMPLEX(KIND=rKind) FUNCTION TraceMatmul_crf(A,B)
!
!Purpose: Function to calculate the trace of complex/real matrices A*B
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: A(:,:)
REAL(KIND=rKind), INTENT(IN) :: B(:,:)
INTEGER i, j
	TraceMatmul_crf=CMPLX(0.0,KIND=rKind)
	DO i=1,SIZE(A,1)
		DO j=1,SIZE(A,2)
			TraceMatmul_crf = TraceMatmul_crf + A(i,j)*B(j,i)
		END DO
	END DO
END FUNCTION TraceMatmul_crf

!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE TraceMatmul !!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION Factorial(n)
!
!Purpose: Return the factorial of an integer n
!  
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER :: k
Factorial = 1
	DO k = 2, n
	Factorial = Factorial * k
	END DO
END FUNCTION Factorial


REAL(KIND=rKind) FUNCTION BinomialCoef(n,m)
!
!Purpose: Return the Binomial Coefficient _nC_m
!
IMPLICIT NONE  
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN) :: m
INTEGER :: k
BinomialCoef=1.0_rKind
	DO k=1,m,1
	BinomialCoef=BinomialCoef*(n-k+1)*1.0_rKind/(k*1.0_rKind)		
	END DO
END FUNCTION BinomialCoef

SUBROUTINE AllocateGamLam(Gammas, Lambdas, chi)
!
!Purpose: Allocate gammas and Lambdas based on a value of chi
!
IMPLICIT NONE  
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER, INTENT(IN) :: chi
INTEGER :: i	
!Gammas live on sites-there are systemSize of them
ALLOCATE(Gammas(systemSize))
!Lambdas live on links, the extra 2 assist in computing two-site observables
ALLOCATE(Lambdas(systemSize+1) , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate lambdas'
			END IF 
	DO i=1,systemSize
	ALLOCATE(Gammas(i)%t(chi,localSize,chi), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate gammas'
			END IF 
	Gammas(i)%t=CMPLX(0.0,KIND=rKind)
	ALLOCATE(Lambdas(i)%v(chi), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate lambdas'
			END IF 
	Lambdas(i)%v=0.0_rKind
	END DO
ALLOCATE(Lambdas(systemSize+1)%v(chi), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate lambdas'
			END IF 
Lambdas(systemSize+1)%v=0.0_rKind


END SUBROUTINE AllocateGamLam

SUBROUTINE CopyGamLam(GammasCopy, LambdasCopy, GammasOrig, LambdasOrig)
!
!Purpose: Copy Gammas and Lambdas from Orig to Copy
!
IMPLICIT NONE
TYPE(tensor), POINTER :: GammasCopy(:), GammasOrig(:)
TYPE(vector), POINTER :: LambdasCopy(:), LambdasOrig(:)
INTEGER :: i, alpha, j, beta, chimin
!If one of the Gammas is larger than the other, only sum up to the 
!maximum indices of the smaller

 
	chimin=MIN(SIZE(GammasCopy(1)%t,3),SIZE(GammasOrig(1)%t,3))
	DO i=1,(systemSize+1)
	LambdasCopy(i)%v=0.0_rKind
		DO alpha=1,chimin
		LambdasCopy(i)%v(alpha)=LambdasOrig(i)%v(alpha)
		END DO
	END DO
	
	DO i=1,systemSize
	GammasCopy(i)%t=0.0_rKind
		DO alpha=1,chimin
			DO beta=1,chimin
				DO j=1,localSize
					GammasCopy(i)%t(alpha,j,beta)=GammasOrig(i)%t(alpha,j,beta)
				END DO
			END DO
		END DO
	END DO
END SUBROUTINE CopyGamLam

SUBROUTINE DeallocateGamLam(Gammas, Lambdas)
!
!Purpose: Deallocate gammas and Lambdas
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER :: i
!Deallocate each site/link object
	DO	i=1,(systemSize)
		DEALLOCATE(Gammas(i)%t, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Gammas'
			END IF 
		DEALLOCATE(Lambdas(i)%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Lambdas'
			END IF 
	END DO
!Deallocate the list of objects
	DEALLOCATE(Lambdas(systemSize+1)%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Lambdas'
			END IF 
	DEALLOCATE(Gammas, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Gammas'
			END IF 
	DEALLOCATE(Lambdas, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Lambdas'
			END IF 
END SUBROUTINE DeallocateGamLam

SUBROUTINE AllocateLabel(LabelLeft, LabelRight, chi)
!
!Purpose: Allocate labels needed for conserving a single Abelian symmetry
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: chi
INTEGER :: i
!LableLeft(l)%vi(alpha) gives the conserved quantity associated with the alpha^th 
!left Schmidt vector for a bipartite splitting at link l, likewise for LabelRight
ALLOCATE(LabelLeft(systemSize+1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate labelLeft'
			END IF 
ALLOCATE(LabelRight(systemSize+1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate labelRight'
			END IF 
	DO i=1, (systemSize+1)
	ALLOCATE(LabelLeft(i)%vi(chi), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate labelLeft'
			END IF 
	ALLOCATE(LabelRight(i)%vi(chi), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate labelRight'
			END IF 
!Label=10000 lets us know that we haven't altered the Label (to avoid confusion from those 
!that are identically zero later)
	LabelLeft(i)%vi=10000
	LabelRight(i)%vi=10000
	END DO
END SUBROUTINE AllocateLabel	

SUBROUTINE CopyLabel(LabLCopy, LabRCopy, LabLOrig, LabROrig)
!
!Purpose: Copy Labels from Orig to Copy
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabLCopy(:), LabRCopy(:), LabLOrig(:), LabROrig(:)
INTEGER :: i, alpha, chimin
!If one of the Labels is larger than the other, only sum up to the 
!maximum indices of the smaller 	
chimin=MIN(SIZE(LabLCopy(1)%vi,1),SIZE(LabLOrig(1)%vi,1))
	DO i=1,(systemSize+1)
	LabLCopy(i)%vi=10000
	LabRCopy(i)%vi=10000
		DO alpha=1,chimin
		LabLCopy(i)%vi(alpha)=LabLOrig(i)%vi(alpha)
		LabRCopy(i)%vi(alpha)=LabROrig(i)%vi(alpha)
		END DO
	END DO
END SUBROUTINE CopyLabel

SUBROUTINE DeallocateLabel(LabelLeft, LabelRight)
!
!Purpose: Deallocate labels needed for conserving a single Abelian symmetry
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: i
	DO i=1, (systemSize+1)
		DEALLOCATE(LabelLeft(i)%vi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate LabelLeft'
			END IF 
		DEALLOCATE(LabelRight(i)%vi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate LabelRight'
			END IF 
	END DO
	DEALLOCATE(LabelLeft, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate LabelLeft'
			END IF 
	DEALLOCATE(LabelRight, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate LabelRight'
			END IF 
END SUBROUTINE DeallocateLabel

SUBROUTINE AllocateIndexLR(indL, indR, BlockSize)
!
!Purpose: Allocate indices for the Block diagonal Theta
!		  Used in Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks
NumOfBlocks=SIZE(BlockSize,1)
ALLOCATE(indL(NumOfBlocks), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate indexLeft'
			END IF 
ALLOCATE(indR(NumOfBlocks), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate indexRight'
			END IF 
	DO k=1,NumOfBlocks
	ALLOCATE(indL(k)%mi(BlockSize(k,2),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate indexLeft'
			END IF 
	ALLOCATE(indR(k)%mi(BlockSize(k,3),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate indexRight'
			END IF 
	indL(k)%mi=0
	indR(k)%mi=0
	END DO
END SUBROUTINE AllocateIndexLR

SUBROUTINE DeallocateIndexLR(indL, indR)
!
!Purpose: Deallocate indicies for the Block diagonal Theta
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER :: k, NumOfBlocks
	NumOfBlocks=SIZE(indL,1)
	DO k=1,NumOfBlocks
		DEALLOCATE(indL(k)%mi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate indL'
			END IF 
		DEALLOCATE(indR(k)%mi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate indR'
			END IF 
	END DO
	DEALLOCATE(indL, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate indL'
			END IF 
	DEALLOCATE(indR, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate indR'
			END IF 
END SUBROUTINE DeallocateIndexLR

SUBROUTINE AllocateBlockTheta(BlockTheta, BlockSize)
!
!Purpose: Allocate Block diagonal Theta based on the number and size of each block
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: BlockTheta(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks
NumOfBlocks=SIZE(BlockSize,1)
ALLOCATE(BlockTheta(NumOfBlocks))
	DO k=1,NumOfBlocks
		!If the block has both left and right Schmidt vectors (i.e. is of nonzero dimension) 
		!then allocate to full size
		IF((BlockSize(k,2).ne.0).AND.(BlockSize(k,3).ne.0)) THEN
		ALLOCATE(BlockTheta(k)%m(BlockSize(k,2),BlockSize(k,3)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockTheta'
			END IF 
		ELSE
		!ELSE allocate a 1 by 1 block
		ALLOCATE(BlockTheta(k)%m(1,1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockTheta'
			END IF 
		END IF
	BlockTheta(k)%m=CMPLX(0.0,KIND=rKind)
	END DO
END SUBROUTINE AllocateBlockTheta

SUBROUTINE DeallocateBlockTheta(BlockTheta)
!
!Purpose: Deallocate Block diagonal Theta based on the number and size of each block
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: BlockTheta(:)
INTEGER :: k, NumOfBlocks
	NumOfBlocks=SIZE(BlockTheta,1)
	DO k=1,NumOfBlocks
		DEALLOCATE(BlockTheta(k)%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate BlockTheta'
			END IF 
	END DO
	DEALLOCATE(BlockTheta, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate BlockTheta'
			END IF 
END SUBROUTINE DeallocateBlockTheta

SUBROUTINE AllocateUSV(US, SS, VS, BlockSize)
!
!Purpose: Allocate stuff for the singular value decomposition LAPACK routine
!		  Specific to Block diagonal Theta 
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: US(:), VS(:)
TYPE(vector), POINTER :: SS(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks
	NumOfBlocks=SIZE(BlockSize,1)
!Allocate a list of SVD variables
!The index corresponds to a "block" of the block diagonal theta
	ALLOCATE(US(NumOfBlocks))
	ALLOCATE(SS(NumOfBlocks))
	ALLOCATE(VS(NumOfBlocks))
		DO k=1,NumOfBlocks
		!Allocate U, S, and V only if the block is greater than 1X1 in size
		IF((BlockSize(k,2).ne.0).AND.(BlockSize(k,3).ne.0)) THEN
			ALLOCATE(US(k)%m(BlockSize(k,2),BlockSize(k,2)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate US'
			END IF 
			ALLOCATE(SS(k)%v(MIN(BlockSize(k,2),BlockSize(k,3))), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SS'
			END IF 
			ALLOCATE(VS(k)%m(BlockSize(k,3),BlockSize(k,3)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate VS'
			END IF 
		ELSE
			ALLOCATE(US(k)%m(1,1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate US'
			END IF 
			ALLOCATE(SS(k)%v(1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SS'
			END IF 
			ALLOCATE(VS(k)%m(1,1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate VS'
			END IF 
		END IF
		US(k)%m=CMPLX(0.0,KIND=rKind)
			SS(k)%v=0.0_rKind
			VS(k)%m=CMPLX(0.0,KIND=rKind)
		END DO
END SUBROUTINE AllocateUSV

SUBROUTINE DeallocateUSV(US, SS, VS)
!
!Purpose: Allocate stuff for the singular value decomposition LAPACK routine
!		  Specific to Block diagonal Theta 
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: US(:), VS(:)
TYPE(vector), POINTER :: SS(:)
INTEGER :: k, NumOfBlocks
	NumOfBlocks=SIZE(SS,1)
	DO k=1,NumOfBlocks
		DEALLOCATE(US(k)%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate US'
			END IF 
		DEALLOCATE(SS(k)%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SS'
			END IF 
		DEALLOCATE(VS(k)%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate VS'
			END IF 
	END DO
	DEALLOCATE(US, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate US'
			END IF 
	DEALLOCATE(SS, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SS'
			END IF 
	DEALLOCATE(VS, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate VS'
			END IF 
END SUBROUTINE DeallocateUSV
	
SUBROUTINE AllocateSSflat(ssflat, BlockSize)
!
!Purpose: Allocate vector to hold the singular values from all blocks
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(vector) :: ssflat
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks, SizeOfSS
	NumOfBlocks=SIZE(BlockSize,1)
	SizeOfSS=0
		DO k=1, NumOfBlocks
		!There are min(m,n) singular values of an mXn matrix
		!Count the singular values from each block
		SizeOfSS=SizeOfSS+MIN(BlockSize(k,2),BlockSize(k,3))
		END DO
	ALLOCATE(ssflat%v(SizeOfSS), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ssflat'
			END IF 
	ssflat%v=0.0_rKind
END SUBROUTINE AllocateSSflat
	
SUBROUTINE AllocateOps(Ops,numops,opsize)
!
!Purpose: Allocate a numops length list of opsizeXopsize matrices, name it Ops
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Ops(:)
INTEGER, INTENT(IN) :: numops,opsize
INTEGER :: i
	ALLOCATE(Ops(numops), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Ops'
			END IF 
	DO i=1,numops
			ALLOCATE(Ops(i)%m(opsize,opsize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Ops'
			END IF 
	END DO
END SUBROUTINE AllocateOps
	
SUBROUTINE DeallocateOps(Ops,numops)
!
!Purpose: Deallocate a numops length list of opsizeXopsize matrices
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Ops(:)
INTEGER, INTENT(IN) :: numops
INTEGER :: i
	DO i=1,numops
		DEALLOCATE(Ops(i)%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Ops'
			END IF 
	END DO
	DEALLOCATE(Ops, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Ops'
			END IF 
END SUBROUTINE DeallocateOps

SUBROUTINE AllocateProp(U)
!
!Purpose: Allocate the Propagator based on boundary conditions
! and order of trotter decomposition
!
IMPLICIT NONE
TYPE(matrix), POINTER :: U(:)
INTEGER :: i

!Second order trotter routines
IF(trotterOrder==2) THEN
	!Open Boundary Conditions
	IF(BoundaryCond=='O') THEN
		!Allocate 1 operator for each bond
		ALLOCATE(U(systemSize-1))
		DO i=1,systemSize-1
			ALLOCATE(U(i)%m(localSize*localSize,localSize*localSize))
		END DO
	!Periodic Boundary Conditions
	ELSE IF(BoundaryCond=='P') THEN
		!Allocate 1 operator for each bond + 1 for the L->1 bond
		ALLOCATE(U(systemSize))
		DO i=1,systemSize
			ALLOCATE(U(i)%m(localSize*localSize,localSize*localSize))
		END DO
	ELSE
		STOP "BoundaryCond specifier not recognized in AllocateProp!  Use O or P."
	END IF
!Fifth order trotter routines
ELSE IF(trotterOrder==5) THEN
	!Open Boundary Conditions
	IF(BoundaryCond=='O') THEN
		!Allocate 2 propagators for each bond
		ALLOCATE(U(2*(systemSize-1)))
		DO i=1,2*systemSize-2
			ALLOCATE(U(i)%m(localSize*localSize,localSize*localSize))
		END DO
	!Periodic Boundary Conditions
	ELSE IF(BoundaryCond=='P') THEN
		IF(MOD(systemSize,2)==0) THEN
		!Allocate 2 propagators for each bond and two for the edge
		ALLOCATE(U(2*systemSize))
		DO i=1,2*systemSize
			ALLOCATE(U(i)%m(localSize*localSize,localSize*localSize))
		END DO

		ELSE
		!Allocate 2 propagators for each odd bond, 4 for each even, and 3 for the edge
		ALLOCATE(U(3*systemSize))
		DO i=1,3*systemSize
			ALLOCATE(U(i)%m(localSize*localSize,localSize*localSize))
		END DO
		END IF

	ELSE
		STOP "BoundaryCond specifier not recognized in AllocateProp!  Use O or P."
	END IF
ELSE
	STOP "trotterOrder specifier not recognized in AllocateProp!  Use 2 or 5"
END IF

END SUBROUTINE AllocateProp

SUBROUTINE DeallocateProp(U)
!
!Purpose: Deallocate the Propagator based on boundary conditions
! and order of trotter decomposition
!
IMPLICIT NONE
TYPE(matrix), POINTER :: U(:)
INTEGER :: i

!Second order trotter routines
IF(trotterOrder==2) THEN
	!Open Boundary Conditions
	IF(BoundaryCond=='O') THEN
		!Allocate 1 operator for each bond
		DO i=1,systemSize-1
			DEALLOCATE(U(i)%m)
		END DO
		DEALLOCATE(U)
	!Periodic Boundary Conditions
	ELSE IF(BoundaryCond=='P') THEN
		!Allocate 1 operator for each bond + 1 for the L->1 bond
		DO i=1,systemSize
			DEALLOCATE(U(i)%m)
		END DO
		DEALLOCATE(U)
	ELSE
		STOP "BoundaryCond specifier not recognized in DeallocateProp!  Use O or P."
	END IF
!Fifth order trotter routines
ELSE IF(trotterOrder==5) THEN
	!Open Boundary Conditions
	IF(BoundaryCond=='O') THEN
		!Allocate 2 propagators for each bond
		DO i=1,2*systemSize-2
			DEALLOCATE(U(i)%m)
		END DO
		DEALLOCATE(U)
	!Periodic Boundary Conditions
	ELSE IF(BoundaryCond=='P') THEN
		IF(MOD(systemSize,2)==0) THEN
		DO i=1,2*systemSize
			DEALLOCATE(U(i)%m)
		END DO
		DEALLOCATE(U)
		ELSE
		DO i=1,3*systemSize
			DEALLOCATE(U(i)%m)
		END DO
		DEALLOCATE(U)
		END IF
	ELSE
		STOP "BoundaryCond specifier not recognized in DeallocateProp!  Use O or P."
	END IF
ELSE
	STOP "trotterOrder specifier not recognized in DeallocateProp!  Use 2 or 5"
END IF

END SUBROUTINE DeallocateProp

SUBROUTINE OLDConstructPropagators(H, U, dtodd, dteven)
!
!Purpose: Construct the Trotter-Suzuki propagator U from the Hamiltonian H
! using a hard-wired routine for the second order scheme
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(matrix), POINTER :: U(:)
COMPLEX(KIND=rKind), INTENT(IN) :: dtodd, dteven
INTEGER :: i

	DO i=1,(systemSize-1),2
		CALL Matrix_Exponential(H(i)%m, U(i)%m, dtodd, localSize*localSize)
	END DO
	DO i=2,(systemSize-1),2
		CALL Matrix_Exponential(H(i)%m, U(i)%m, dteven, localSize*localSize)
	END DO

END SUBROUTINE OLDConstructPropagators

SUBROUTINE NEWConstructPropagators(H, U, dt)
!
!Purpose: Construct the Trotter-Suzuki propagator U from the Hamiltonian H
!Using a routine adaptive to different boundary conditions and trotter schemes
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(matrix), POINTER :: U(:)
COMPLEX(KIND=rKind), INTENT(IN) :: dt
INTEGER :: i
REAL(KIND=rKIND) :: FRtheta

IF(trotterOrder==2) THEN
	IF(BoundaryCond=='O') THEN
		!Odd bonds
		DO i=1,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(i)%m, 0.5_rKind*dt, localSize*localSize)
		END DO
		!Even bonds
		DO i=2,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(i)%m, dt, localSize*localSize)
		END DO

	ELSE IF(BoundaryCond=='P') THEN
		! if there are an even number of lattice sites, [H_{even},H_{edge}]=0 and we can
		! take a full dt even step
		IF (MOD(systemSize,2)==0) THEN
			!Odd bonds
			DO i=1,(systemSize-1),2
				CALL Matrix_Exponential(H(i)%m, U(i)%m, 0.5_rKind*dt, localSize*localSize)
			END DO
			!Even bonds
			DO i=2,(systemSize-1),2
				CALL Matrix_Exponential(H(i)%m, U(i)%m, dt, localSize*localSize)
			END DO
			!Edge bond
			CALL Matrix_Exponential(H(systemSize)%m, U(systemSize)%m, dt, localSize*localSize)

		! if there are an even number of lattice sites, [H_{even},H_{edge}].ne.0 and we have
		! to take a half dt even steps
		ELSE
			!Odd bonds
			DO i=1,(systemSize-1),2
				CALL Matrix_Exponential(H(i)%m, U(i)%m, 0.5_rKind*dt, localSize*localSize)
			END DO
			!Even bonds
			DO i=2,(systemSize-1),2
				CALL Matrix_Exponential(H(i)%m, U(i)%m, 0.5_rKind*dt, localSize*localSize)
			END DO
			!Edge bond
			CALL Matrix_Exponential(H(systemSize)%m, U(systemSize)%m, dt, localSize*localSize)
		END IF
	ELSE
		STOP "BoundCond speicification not recognized in ConstructPropagators!"
	END IF

ELSE IF(trotterOrder==5) THEN
FRtheta=1.0_rKind/(2.0_rKind-(2.0_rKind**(1.0/3.0))) !Forest-ruth theta
	IF(BoundaryCond=='O') THEN
		!Odd bonds
		DO i=1,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+1)%m, 0.5_rKind*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+2)%m, 0.5_rKind*(1.0_rKind-FRtheta)*dt, localSize*localSize)
		END DO

		!Even bonds
		DO i=2,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+1)%m, FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+2)%m,(1.0_rKind-2.0_rKind*FRtheta)*dt, localSize*localSize)
		END DO

	ELSE IF(BoundaryCond=='P') THEN
		IF (MOD(systemSize,2)==0) THEN
		!Odd bonds
		DO i=1,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+1)%m, 0.5_rKind*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+2)%m, 0.5_rKind*(1.0_rKind-FRtheta)*dt, localSize*localSize)
		END DO

		!Even bonds
		DO i=2,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+1)%m, FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(2*(i-1)+2)%m,(1.0_rKind-2.0_rKind*FRtheta)*dt, localSize*localSize)
		END DO
			!Edge bond
			CALL Matrix_Exponential(H(systemSize)%m, U(2*(systemSize-1)+1)%m, FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(systemSize)%m, U(2*(systemSize-1)+2)%m, (1.0_rKind-2.0_rKind*FRtheta)*dt, localSize*localSize)
		ELSE

		!Odd bonds
		DO i=1,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(3*(i-1)+1)%m, 0.5_rKind*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(3*(i-1)+2)%m, 0.5_rKind*(1.0_rKind-FRtheta)*dt, localSize*localSize)
		END DO

		!Even bonds
		DO i=2,(systemSize-1),2
			CALL Matrix_Exponential(H(i)%m, U(2+3*(i-2)+1)%m, 0.5_rKind*FRtheta*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(2+3*(i-2)+2)%m,0.5_rKind*(1.0_rKind-FRtheta)*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(2+3*(i-2)+3)%m,0.5_rKind*(1.0_rKind-2.0_rKind*FRtheta)*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(i)%m, U(2+3*(i-2)+4)%m,0.5_rKind*(1.0_rKind-FRtheta)*(1.0_rKind-2.0_rKind*FRtheta)*dt, localSize*localSize)
		END DO
			!Edge bond
			CALL Matrix_Exponential(H(systemSize)%m, U(3*(systemSize-1)+1)%m, FRtheta*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(systemSize)%m, U(3*(systemSize-1)+2)%m, (1.0_rKind-2.0_rKind*FRtheta)*FRtheta*dt, localSize*localSize)
			CALL Matrix_Exponential(H(systemSize)%m, U(3*(systemSize-1)+3)%m, (1.0_rKind-2.0_rKind*FRtheta)*(1.0_rKind-2.0_rKind*FRtheta)*dt, localSize*localSize)
		END IF

	ELSE
		STOP "BoundaryCond specifier not recognized in ConstructPropagators!  Use O or P."
	END IF

ELSE
	STOP "trotterOrder specifier not recognized in ConstructPropagators!  Use 2 or 5"
END IF


END SUBROUTINE NEWConstructPropagators



END MODULE TEBDtools_module
