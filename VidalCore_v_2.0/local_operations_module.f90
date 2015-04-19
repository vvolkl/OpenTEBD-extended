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
MODULE local_operations_module
!
! Purpose: Module to perform local operations
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

IMPLICIT NONE


INTERFACE OneSiteOp
MODULE PROCEDURE OneSiteOp_r,OneSiteOp_c
END INTERFACE OneSiteOp

INTERFACE ThetaOperation
MODULE PROCEDURE ThetaOperation_r, ThetaOperation_c
END INTERFACE ThetaOperation

INTERFACE TwoSiteOp
MODULE PROCEDURE TwoSiteOp_r,TwoSiteOp_c
END INTERFACE TwoSiteOp

INTERFACE ThetaOperationNC
MODULE PROCEDURE ThetaOperationNC_r,ThetaOperationNC_c
END INTERFACE ThetaOperationNC

INTERFACE TwoSiteOpNC
MODULE PROCEDURE TwoSiteOpNC_r, TwoSiteOpNC_c
END INTERFACE TwoSiteOpNC

! These variables are set for using the subroutine "ZGESVD" in LAPACK, which performs an SVD on a general matrix.
CHARACTER(1) :: jobu_SVD, jobvt_SVD
INTEGER :: matrixSizeSM_SVD, workSizeSM_SVD, matrixSizeLG_SVD, workSizeLG_SVD, &
matrixSize_SVD, workSize_SVD, info_SVD, matrixSizeL_SVD, matrixSizeT_SVD

REAL(KIND=rKind), ALLOCATABLE :: rworkSM_SVD(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: workSM_SVD(:)
REAL(KIND=rKind), ALLOCATABLE :: rworkLG_SVD(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: workLG_SVD(:)
REAL(KIND=rKind), ALLOCATABLE :: rwork_SVD(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: work_SVD(:)	
	
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Number Non-conserving method starts !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SVDInit(chi)
!
!Purpose: Allocate variables to perform SVD using ZGESVD in LAPACK
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: chi
	jobu_SVD='A'
	jobvt_SVD='A'
	matrixSizeSM_SVD=localSize
	workSizeSM_SVD=5*matrixSizeSM_SVD
	matrixSizeLG_SVD=chi*localSize
	workSizeLG_SVD=5*matrixSizeLG_SVD
	matrixSize_SVD=chi
	workSize_SVD=5*matrixSize_SVD
	ALLOCATE(rworkSM_SVD(workSizeSM_SVD), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVD variables'
			END IF
	ALLOCATE(workSM_SVD(workSizeSM_SVD), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVD variables'
			END IF
	ALLOCATE(rworkLG_SVD(workSizeLG_SVD), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVD variables'
			END IF
	ALLOCATE(workLG_SVD(workSizeLG_SVD), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVD variables'
			END IF	
	ALLOCATE(rwork_SVD(workSize_SVD), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVD variables'
			END IF
	ALLOCATE(work_SVD(workSize_SVD), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVD variables'
			END IF	
END SUBROUTINE SVDInit
	
SUBROUTINE SVDFinish()
!
!Purpose: Deallocate SVD variables
!
IMPLICIT NONE	
	DEALLOCATE(rworkSM_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables'
			END IF
	DEALLOCATE(workSM_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables'
			END IF
	DEALLOCATE(rworkLG_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables'
			END IF
	DEALLOCATE(workLG_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables'
			END IF
	DEALLOCATE(rwork_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables'
			END IF
	DEALLOCATE(work_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables'
			END IF
END SUBROUTINE SVDFinish

SUBROUTINE OneSiteOp_r(Op1,Gamma)
!
!Purpose: Perform the real one-site operation Op1 on the Gamma specified
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:,:)
INTEGER alpha,chi
	chi = SIZE(Gamma,1);
	DO alpha = 1,chi
		Gamma(alpha,:,:) = MATMUL(Op1,Gamma(alpha,:,:));
	END DO		
END SUBROUTINE OneSiteOp_r

SUBROUTINE OneSiteOp_c(Op1,Gamma)
!
!Purpose: Perform the complex one-site operation Op1 on the Gamma specified
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma(:,:,:)
INTEGER alpha,chi
	chi = SIZE(Gamma,1);
	DO alpha = 1,chi
		Gamma(alpha,:,:) = MATMUL(Op1,Gamma(alpha,:,:));
	END DO		
END SUBROUTINE OneSiteOp_c

		
SUBROUTINE FormTheta(Theta, Lambda0, Gamma1, Lambda1, Gamma2, Lambda2)
!
!Purpose: Form Theta as defined in the Manual
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma1(:,:,:), Gamma2(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda0(:), Lambda1(:), Lambda2(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER :: chi0,chi1,chi2,alpha,beta,gamma,i,j,k,l
		chi0 = SIZE(Theta,1);
		chi1 = SIZE(Gamma1,3)
		chi2 = SIZE(Theta,4);
		DO alpha=1,chi0
			DO beta=1,chi2
				DO i=1,localSize
					DO j=1,localSize
						ThetaTemp(alpha,i,j,beta) = CMPLX(0.0,KIND=rKind);					
						DO gamma=1,chi1
							ThetaTemp(alpha,i,j,beta)=ThetaTemp(alpha,i,j,beta) & 
				         	+ Lambda0(alpha)*Gamma1(alpha,i,gamma)*Lambda1(gamma)*Gamma2(gamma,j,beta)*Lambda2(beta)
						END DO 
					END DO 
				END DO 
			END DO 
		END DO 
		Theta = ThetaTemp
END SUBROUTINE FormTheta

SUBROUTINE ThetaOperation_r(Op2,Theta)
!
!Purpose: Operate the real two-site operation Op2 on Theta
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN)  :: Op2(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER :: chi0,chi2,alpha,beta,gamma,i,j,k,l
		
	chi0 = SIZE(Theta,1)
	chi2 = SIZE(Theta,4)
		
	ThetaTemp=Theta
		
		DO alpha=1,chi0
			DO beta=1,chi2
				DO i=1,localSize
					DO j=1,localSize
						Theta(alpha,i,j,beta) = CMPLX(0.0,KIND=rKind);					
						DO k=1,localSize
						DO l=1,localSize
								Theta(alpha,i,j,beta)=Theta(alpha,i,j,beta) & 
									             + Op2((i-1)*localSize+j,(k-1)*localSize+l) &
									             * ThetaTemp(alpha,k,l,beta)
							END DO 
							END DO 
					END DO 
				END DO 
			END DO 
		END DO 
	
END SUBROUTINE ThetaOperation_r


SUBROUTINE ThetaOperation_c (Op2,Theta)
!
!Purpose: Operate the complex two-site operation Op2 on Theta
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN)  :: Op2(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER :: chi0,chi2,alpha,beta,gamma,i,j,k,l
		
	chi0 = SIZE(Theta,1)
	chi2 = SIZE(Theta,4)
		
	ThetaTemp=Theta
		
		DO alpha=1,chi0
			DO beta=1,chi2
				DO i=1,localSize
					DO j=1,localSize
						Theta(alpha,i,j,beta) = CMPLX(0.0,KIND=rKind);					
						DO k=1,localSize
						DO l=1,localSize
								Theta(alpha,i,j,beta)=Theta(alpha,i,j,beta) & 
									             + Op2((i-1)*localSize+j,(k-1)*localSize+l) &
									             * ThetaTemp(alpha,k,l,beta)
							END DO 
							END DO 
					END DO 
				END DO 
			END DO 
		END DO 
	
END SUBROUTINE ThetaOperation_c

	
SUBROUTINE ReshapeTheta(Theta,ThetaRS)
!
!Purpose: Reshape the chiXdXdXchi 4-tensor Theta into a (chi d)X(chi d) matrix ThetaRS and renormalize such that
!\sum_{aijb}|Theta_{ab}^{ij}|^2=1
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind), INTENT(OUT) :: ThetaRS(SIZE(Theta,1)*SIZE(Theta,2),SIZE(Theta,3)*SIZE(Theta,4))
COMPLEX(KIND=rKind) :: normSQcmplx
REAL(KIND=rKind) :: norm
INTEGER :: chi0,chi1,alpha,beta,gamma,i,j,k,l
		
	chi0 = SIZE(Theta,1);
	chi1 = SIZE(Theta,4);
		
	normSQcmplx=CMPLX(0.0,KIND=rKind)
		DO alpha=1,chi0
			DO beta=1,chi1
				DO i=1,localSize
					DO j=1,localSize
						normSQcmplx=normSQcmplx+CONJG(Theta(alpha,i,j,beta))*Theta(alpha,i,j,beta)
					END DO 
				END DO 
			END DO 
		END DO 
		
		norm=SQRT(REAL(normSQcmplx,KIND=rKind))
		
		DO alpha=1,chi0
			DO beta=1,chi1
				DO i=1,localSize
					DO j=1,localSize
						ThetaRS((i-1)*chi0+alpha,(j-1)*chi1+beta)=Theta(alpha,i,j,beta)/norm
					END DO 
				END DO 
			END DO 
		END DO 	
END SUBROUTINE ReshapeTheta

SUBROUTINE SVDTruncation(link, MatrixIn, S, U, V)
!
!Purpose: Perform an SVD on the reshaped Theta, keep only the largest chi singular values
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(INOUT) :: MatrixIn(:,:)       
REAL(KIND=rKind), INTENT(OUT) :: S(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: U(:,:), V(:,:)
COMPLEX(KIND=rKind) :: VT(SIZE(U,1),SIZE(U,2))
INTEGER :: i, j, dchi1, dchi2

!Call the LAPACK routine ZGESVD, which performs a SVD on a general matrix
	CALL ZGESVD(jobu_SVD, jobvt_SVD, matrixSizeLG_SVD, matrixSizeLG_SVD, MatrixIn, matrixSizeLG_SVD, S, U, & 
				matrixSizeLG_SVD, VT, matrixSizeLG_SVD, workLG_SVD, workSizeLG_SVD, rworkLG_SVD, info_SVD)

	dchi1=SIZE(U,1)
	dchi2=SIZE(U,2)
		DO i=1,dchi1
			DO j=1,dchi2
				V(i,j)=VT(i,j)
			END DO
		END DO
END SUBROUTINE SVDTruncation

SUBROUTINE FormLambda1(Lambda1, truncerr, S, chi1)
!
!Purpose: Form Lambda from the singular values
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: chi1
REAL(KIND=rKind), INTENT(OUT) :: truncerr
REAL(KIND=rKind), INTENT(INOUT) :: Lambda1(:)
REAL(KIND=rKind), INTENT(IN) :: S(:)
REAL(KIND=rKind) :: norm
INTEGER :: alpha

!Renormalize the kept singular values such that the sum of the squares is 1
			norm=0.0_rKind
		DO alpha=1, chi1
			IF(ABS(S(alpha))>1E-15) THEN
				norm=norm + S(alpha) ** 2
			ELSE
				norm=norm
			END IF
		END DO

!Schmidt error		
		truncerr =ABS( 1.0_rKind - norm)
!Redefine Lambda		
		DO alpha=1, chi1
			IF(ABS(S(alpha))>1E-15) THEN
				Lambda1(alpha)=S(alpha)/SQRT(norm)
			ELSE
				Lambda1(alpha)=0.0_rKind
			END IF
		END DO		
END SUBROUTINE FormLambda1

SUBROUTINE FormGamma1(Lambda0,Gamma1,U,chi0,chi1)
!
!Purpose: Form Gamma to the left of the link from the SVD.
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma1(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda0(:)
COMPLEX(KIND=rKind), INTENT(IN) :: U(:,:)
INTEGER, INTENT(IN) :: chi0, chi1
INTEGER :: alpha, beta, i
		
		DO beta=1,chi1
			DO alpha=1,chi0
				DO i=1,localSize
					IF(ABS(Lambda0(alpha))>1E-15) THEN
						Gamma1(alpha,i,beta)=U((i-1)*chi0+alpha,beta)/Lambda0(alpha)
					ELSE
						Gamma1(alpha,i,beta)=CMPLX(0.0,KIND=rKind)
					END IF
				END DO 
			END DO		
		END DO	
		
END SUBROUTINE FormGamma1
	
SUBROUTINE FormGamma2(Gamma2,Lambda2,V,chi1,chi2)
!
!Purpose: Form Gamma to the right of the link from the SVD.
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma2(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda2(:)
COMPLEX(KIND=rKind), INTENT(IN) :: V(:,:)
INTEGER, INTENT(IN) :: chi1, chi2
INTEGER :: alpha, beta, i
		
		DO beta=1,chi2
			DO alpha=1,chi1
				DO i=1,localSize
					IF(ABS(Lambda2(beta))>1E-15) THEN
						Gamma2(alpha,i,beta)=V(alpha,(i-1)*chi2+beta)/Lambda2(beta)
					ELSE
						Gamma2(alpha,i,beta)=CMPLX(0.0,KIND=rKind)
					END IF
				END DO 
			END DO		
		END DO	
END SUBROUTINE FormGamma2
	
SUBROUTINE TwoSiteOp_r(link,Op2,Gammas,Lambdas,truncerr)
!
!Purpose: Perform the real two-site operation Op2 on the sites neighboring link
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: ThetaRS(SIZE(Gammas(link)%t,1)*SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2)*SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: U(localSize*Size(Gammas(link)%t,1),localSize*Size(Gammas(link)%t,3))
COMPLEX(KIND=rKind) :: V(localSize*Size(Gammas(link+1)%t,1),localSize*Size(Gammas(link+1)%t,3))
REAL(KIND=rKind) :: S(localSize*Size(Gammas(link)%t,3))
INTEGER :: chi0, chi1, chi2

	chi0=Size(Gammas(link)%t,1)
	chi1=Size(Gammas(link)%t,3)
	chi2=Size(Gammas(link+1)%t,3)

	CALL FormTheta(Theta,Lambdas(link)%v,Gammas(link)%t,Lambdas(link+1)%v, &
	Gammas(link+1)%t,Lambdas(link+2)%v)
	CALL ThetaOperation(Op2,Theta)
	CALL ReshapeTheta(Theta,ThetaRS)
	CALL SVDTruncation(link,ThetaRS, S, U, V)
	CALL FormLambda1(Lambdas(link+1)%v, truncerr, S, chi1)
	CALL FormGamma1(Lambdas(link)%v,Gammas(link)%t,U,chi0,chi1)
	CALL FormGamma2(Gammas(link+1)%t,Lambdas(link+2)%v,V,chi1,chi2)
END SUBROUTINE TwoSiteOp_r

SUBROUTINE TwoSiteOp_c(link,Op2,Gammas,Lambdas,truncerr)
!
!Purpose: Perform the complex two-site operation Op2 on the sites neighboring link
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: ThetaRS(SIZE(Gammas(link)%t,1)*SIZE(Gammas(link)%t,2), &
SIZE(Gammas(link+1)%t,2)*SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: U(localSize*Size(Gammas(link)%t,1),localSize*Size(Gammas(link)%t,3))
COMPLEX(KIND=rKind) :: V(localSize*Size(Gammas(link+1)%t,1),localSize*Size(Gammas(link+1)%t,3))
REAL(KIND=rKind) :: S(localSize*Size(Gammas(link)%t,3))
INTEGER :: chi0, chi1, chi2

	chi0=Size(Gammas(link)%t,1)
	chi1=Size(Gammas(link)%t,3)
	chi2=Size(Gammas(link+1)%t,3)

	CALL FormTheta(Theta,Lambdas(link)%v,Gammas(link)%t,Lambdas(link+1)%v, &
	Gammas(link+1)%t,Lambdas(link+2)%v)
	CALL ThetaOperation(Op2,Theta)
	CALL ReshapeTheta(Theta,ThetaRS)
	CALL SVDTruncation(link,ThetaRS, S, U, V)
	CALL FormLambda1(Lambdas(link+1)%v, truncerr, S, chi1)
	CALL FormGamma1(Lambdas(link)%v,Gammas(link)%t,U,chi0,chi1)
	CALL FormGamma2(Gammas(link+1)%t,Lambdas(link+2)%v,V,chi1,chi2)
END SUBROUTINE TwoSiteOp_c
	
SUBROUTINE CanonicalForm(Lambda0,Gamma1,Lambda1,Gamma2,Lambda2)
!
!Purpose: Put the tensor network defined by the Bipartite splitting at link into canonical form.
!This amounts to reorthogonalizing the Schmidt basis.
!
!From PRA 74, 022320 "A tree [tensor network] is in the canonical form for bipartition A:B if
!(i) the wights on the connecting index \alpha correspond to the Schmidt coefficients {\lambda_{\alpha}} and
!(ii) the subtrees A and B describe a set of Schmidt bases {|\psi_{\alpha}^{[A]}>} and {|\psi_{\alpha}^{[B]}>}
!
!See manual for more detail
!

IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma1(:,:,:), Gamma2(:,:,:)
REAL(KIND=rKind), INTENT(INOUT) :: Lambda1(:), Lambda0(:), Lambda2(:)
COMPLEX(KIND=rKind) :: M(SIZE(Gamma1,3),SIZE(Gamma1,3)), X(SIZE(Gamma1,3),SIZE(Gamma1,3)), Y(SIZE(Gamma1,3),SIZE(Gamma1,3))
COMPLEX(KIND=rKind) :: U(SIZE(Gamma1,3),SIZE(Gamma1,3)), VT(SIZE(Gamma1,3),SIZE(Gamma1,3))
COMPLEX(KIND=rKind) :: Gamma1Temp(SIZE(Gamma1,1),SIZE(Gamma1,2),SIZE(Gamma1,3))
COMPLEX(KIND=rKind) :: Gamma2Temp(SIZE(Gamma2,1),SIZE(Gamma2,2),SIZE(Gamma2,3))
REAL(KIND=rKind) :: S(SIZE(Gamma1,3))
REAL(KIND=rKind) :: normsq
INTEGER :: i,j,k,l,alpha,beta,gamma,chi0,chi1,chi2
		chi0=SIZE(Gamma1,1)
		chi1=SIZE(Gamma1,3)
		chi2=SIZE(Gamma2,3)
				DO alpha=1,chi1
					DO beta=1,chi1
					M(alpha,beta)=CMPLX(0.0,KIND=rKind)
						DO i=1,localSize
							DO gamma=1,chi0
						M(alpha,beta)=M(alpha,beta)+CONJG(Gamma1(gamma,i,alpha))*Lambda0(gamma)*Lambda0(gamma)*Gamma1(gamma,i,beta)
							END DO 
						END DO	
					END DO	
				END DO	
		CALL SVD(M, U, S, VT)
				DO alpha=1,chi0
					DO beta=1,chi1
						DO i=1,localSize
							Gamma1Temp(alpha,i,beta) = CMPLX(0.0,KIND=rKind)
							IF(S(beta)>1E-8) THEN
								DO gamma=1,chi1
									Gamma1Temp(alpha,i,beta) = Gamma1Temp(alpha,i,beta) &
											     	 + Gamma1(alpha,i,gamma)*U(gamma,beta)*((1.0_rKind)/SQRT(S(beta)))
								END DO 
							END IF
						END DO 
					END DO 
				END DO 
				
				
		DO alpha=1,chi1
			DO beta=1,chi1
				X(alpha,beta)=SQRT(S(alpha))*CONJG(U(beta,alpha))
			END DO
		END DO
		
		
		DO alpha=1,chi1
			DO beta=1,chi1
				M(alpha,beta)=CMPLX(0.0,KIND=rKind)
				DO i=1,localSize
					DO gamma=1,chi2
						M(alpha,beta)=M(alpha,beta)+Gamma2(alpha,i,gamma)*Lambda2(gamma)*Lambda2(gamma)*CONJG(Gamma2(beta,i,gamma))
					END DO
				END DO
			END DO
		END DO
		
		CALL SVD(M, U, S, VT)
		
		DO alpha=1,chi1
			DO beta=1,chi2
				DO i=1,localSize
					Gamma2Temp(alpha,i,beta) = CMPLX(0.0,KIND=rKind)
					IF(S(alpha)>1E-8) THEN
						DO gamma=1,chi1
							Gamma2Temp(alpha,i,beta) = Gamma2Temp(alpha,i,beta) &
												 	+ ((1.0_rKind)/SQRT(S(alpha)))*CONJG(U(gamma,alpha))*Gamma2(gamma,i,beta)
						END DO
					END IF
				END DO
			END DO
		END DO
		
		
		DO alpha=1,chi1
			DO beta=1,chi1
				Y(alpha,beta)=Lambda1(alpha)*U(alpha,beta)*SQRT(S(beta))
			END DO
		END DO
		
		
		M=MATMUL(X,Y)
		CALL SVD(M, U, S, VT)
		normsq = 0.0_rKind
		
		DO alpha=1,chi1
			normsq = normsq + S(alpha)*S(alpha)
		END DO
		
		Lambda1 = S/SQRT(normsq)
		DO alpha=1,chi0
			DO beta=1,chi1
				DO i=1,localSize
					Gamma1(alpha,i,beta)=CMPLX(0.0,KIND=rKind)
					DO gamma=1,chi1
						Gamma1(alpha,i,beta) = Gamma1(alpha,i,beta)+Gamma1Temp(alpha,i,gamma)*U(gamma,beta)
					END DO
				END DO
			END DO
		END DO
		
		
		DO alpha=1,chi1
			DO beta=1,chi2
				DO i=1,localSize
					Gamma2(alpha,i,beta) = CMPLX(0.0,KIND=rKind)
					DO gamma=1,chi1
						Gamma2(alpha,i,beta) = Gamma2(alpha,i,beta) + VT(alpha,gamma)*Gamma2Temp(gamma,i,beta)
					END DO
				END DO
			END DO
		END DO
END SUBROUTINE CanonicalForm
	
SUBROUTINE SVD(MatrixIn, U, S, VT)
!
!Purpose: Perform an SVD using the ZGESVD LAPACK routine
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: MatrixIn(:,:)       
REAL(KIND=rKind), INTENT(OUT) :: S(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: U(:,:), VT(:,:)

		CALL ZGESVD(jobu_SVD, jobvt_SVD, matrixSize_SVD, matrixSize_SVD, MatrixIn, matrixSize_SVD, S, U, & 
					matrixSize_SVD, VT, matrixSize_SVD, work_SVD, workSize_SVD, rwork_SVD, info_SVD)
END SUBROUTINE SVD
		
SUBROUTINE SwapTheta(Theta,ThetaSW)							
!
!Purpose: Swap indices on a reshaped theta
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Theta(:,:,:,:)						
COMPLEX(KIND=rKind), INTENT(OUT) :: ThetaSW(SIZE(Theta,1)*SIZE(Theta,3),SIZE(Theta,2)*SIZE(Theta,4))						
INTEGER :: chi0,chi2,alpha,gamma,i,j						
								
	chi0 = SIZE(Theta,1);						
	chi2 = SIZE(Theta,4);						
								
DO alpha=1,chi0						
	DO gamma=1,chi2					
		DO i=1,localSize				
			DO j=1,localSize			
				ThetaSW((j-1)*chi0+alpha,(i-1)*chi2+gamma)=Theta(alpha,i,j,gamma)		
			END DO			
		END DO				
	END DO					
END DO						
END SUBROUTINE SwapTheta							
	
	
SUBROUTINE Swapping(link,truncerr,Gammas,Lambdas)
!
!Purpose: Swap the local indices of the current bipartite splitting and then put into canonical form
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
REAL(KIND=rKind), INTENT(OUT) :: truncerr
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
	SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: ThetaSW(SIZE(Gammas(link)%t,1)*SIZE(Gammas(link)%t,2), &
	SIZE(Gammas(link+1)%t,2)*SIZE(Gammas(link+1)%t,3))
COMPLEX(KIND=rKind) :: U(localSize*Size(Gammas(link)%t,1),localSize*Size(Gammas(link)%t,3))
COMPLEX(KIND=rKind) :: V(localSize*Size(Gammas(link+1)%t,1),localSize*Size(Gammas(link+1)%t,3))
REAL(KIND=rKind) :: S(localSize*Size(Gammas(link)%t,3))
INTEGER :: chi0, chi1, chi2

chi0=Size(Gammas(link)%t,1)
chi1=Size(Gammas(link)%t,3)
chi2=Size(Gammas(link+1)%t,3)
PRINT *, 'chis',chi0,chi1,chi2,localsize
CALL FormTheta(Theta,Lambdas(link)%v,Gammas(link)%t,Lambdas(link+1)%v, &
	Gammas(link+1)%t,Lambdas(link+2)%v)
PRINT *, 'theta formed'
CALL SwapTheta(Theta,ThetaSW)
PRINT *, 'theta swapped'
CALL SVDTruncation(link,ThetaSW, S, U, V)
PRINT *, 'SVD'
CALL FormLambda1(Lambdas(link+1)%v,truncerr,S,chi1)
PRINT *, 'L1'
CALL FormGamma1(Lambdas(link)%v,Gammas(link)%t,U,chi0,chi1)
PRINT *, 'G1'
CALL FormGamma2(Gammas(link+1)%t,Lambdas(link+2)%v,V,chi1,chi2)
PRINT *, 'G2'

END SUBROUTINE Swapping
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Number conserving method starts !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SVDInitNC(k, BlockSize)
!
!Purpose: Allocate variables for SVD
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(IN) :: BlockSize(:,:)
		jobu_SVD='A'
		jobvt_SVD='A'
		matrixSizeL_SVD=BlockSize(k,2)
		matrixSizeT_SVD=BlockSize(k,3)
		workSize_SVD=5*MAX(matrixSizeL_SVD,matrixSizeT_SVD)
		ALLOCATE(work_SVD(workSize_SVD), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVDNC variables'
			END IF
		ALLOCATE(rwork_SVD(5*MIN(matrixSizeL_SVD,matrixSizeT_SVD)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate SVDNC variables'
			END IF

END SUBROUTINE SVDInitNC

SUBROUTINE SVDFinishNC()
!
!Purpose: Deallocate variables for SVD
!
IMPLICIT NONE
		DEALLOCATE(work_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVDNC variables'
			END IF
		DEALLOCATE(rwork_SVD, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVDNC variables'
			END IF
END SUBROUTINE SVDFinishNC

SUBROUTINE FormThetaNC(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight, intDegFree)
!
!Purpose: Form Theta consistent with number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, i, j, alpha, beta, gamma, lftN1, lftN2, rgtN1, rgtN2, ni, nj
	
		
		chi=SIZE(Lambdas(link+1)%v,1)
! Initialize Theta
		Theta=CMPLX(0.0,KIND=rKind)
		
	
!If internal degrees of freedom are present, we must sum over them with
!a constraint that number is conserved
IF(PRESENT(intDegFree)) THEN
		
			DO alpha=1,chi
			DO gamma=1,chi
				DO beta=1,chi
!We still have to sum over internal degrees of freedom
				DO i=1,localSize,1
!number on left site
				ni=Conserv%vi(i)
				DO j=1,localSize,1
!number on right site
				nj=Conserv%vi(j)
! Number of particles on the left side of 'link'.
					lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
					lftN2=LabelLeft(link+1)%vi(beta)
! Number of particles on the right side of 'link+2'.
					rgtN1=LabelRight(link+2)%vi(gamma)
! Number of particles on the right side of 'link+1'.
					rgtN2=LabelRight(link+1)%vi(beta)
! So that lambda is not equal to zero
!Label*=1000 means that it has been left uninitialized (incoming state does not have enough entanglement)
					IF((MAX(lftN1,lftN2,rgtN1,rgtN2)<1000) &
! 1<=i<=d
					.AND.(lftN2-lftN1<=maxFilling).AND.(lftN2-lftN1>=0) &
! 1<=j<=d
					.AND.(rgtN2-rgtN1<=maxFilling).AND.(rgtN2-rgtN1>=0) &
! So that lftN1+n_i=lftN2 and rgtN1+n_j=rgtN2
					.AND.(lftN1+ni==lftN2) &
					.AND.(rgtN1+nj==rgtN2)) THEN

						Theta(alpha,i,j,gamma) = &
						Theta(alpha,i,j,gamma) + &
						Lambdas(link)%v(alpha) * Gammas(link)%t(alpha,i,beta) * &
						Lambdas(link+1)%v(beta) * Gammas(link+1)%t(beta,j,gamma) * &
						Lambdas(link+2)%v(gamma)

					ELSE
					END IF
				END DO
			END DO
		END DO
		END DO
		END DO	
		
!Else number completely specifies the state, and we avoid summing over local dimension		
ELSE
		
		DO alpha=1,chi
			DO gamma=1,chi
				DO beta=1,chi
! Number of particles on the left side of 'link'.
					lftN1=LabelLeft(link)%vi(alpha)
! Number of particles on the left side of 'link+1'.
					lftN2=LabelLeft(link+1)%vi(beta)
! Number of particles on the right side of 'link+2'.
					rgtN1=LabelRight(link+2)%vi(gamma)
! Number of particles on the right side of 'link+1'.
					rgtN2=LabelRight(link+1)%vi(beta)
! So that lambda is not equal to zero
					IF((MAX(lftN1,lftN2,rgtN1,rgtN2)<1000) &
! 1<=i<=d
					.AND.(lftN2-lftN1+1<=localSize).AND.(lftN2-lftN1+1>=1) &
! 1<=j<=d
					.AND.(rgtN2-rgtN1+1<=localSize).AND.(rgtN2-rgtN1+1>=1)) THEN
! So that lftN1+i-1=lftN2 and rgtN1+j-1=rgtN2
						Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) = &
						Theta(alpha,lftN2-lftN1+1,rgtN2-rgtN1+1,gamma) + &
						Lambdas(link)%v(alpha) * Gammas(link)%t(alpha,lftN2-lftN1+1,beta) * &
						Lambdas(link+1)%v(beta) * Gammas(link+1)%t(beta,rgtN2-rgtN1+1,gamma) * &
						Lambdas(link+2)%v(gamma)
					ELSE
					END IF
				END DO
			END DO
		END DO
		
END IF		
END SUBROUTINE FormThetaNC


SUBROUTINE ThetaOperationNC_r(link, Op2, Theta, LabelLeft, LabelRight, intDegFree)
!
!Purpose:Operate on Theta with the real two-site operation Op2 consistent with number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, i,j,l, k, alpha, gamma, lftN, rgtN,nk,nj,nl,ni
		
	chi=SIZE(Theta,1)
	ThetaTemp=CMPLX(0.0,KIND=rKind)
		
!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
		DO alpha=1,chi
			DO k=1,localSize
			nk=Conserv%vi(k)
			DO j=1, localSize
			nj=Conserv%vi(j)
				DO gamma=1,chi
! Number of particles on the left side of 'link'.
					lftN=LabelLeft(link)%vi(alpha)
! Number of particles on the right side of 'link+2'.
					rgtN=LabelRight(link+2)%vi(gamma)
! So that corresponding lambda is not equal to zero.
					IF((MAX(lftN,rgtN)<1000).AND. &
! 0<=n_k<=maxFilling
					(totNum-lftN-rgtN-nk>=0).AND.(totNum-lftN-rgtN-nk<=maxFilling) &
					.AND.(nj==totNum-lftN-rgtN-nk)) THEN
						ThetaTemp(alpha,k,j,gamma)=CMPLX(0.0,KIND=rKind)
						DO i=1,localSize
						ni=Conserv%vi(i)
						DO l=1,localSize
						nl=Conserv%vi(l)
! 0<=n_i<=maxfilling
							IF((totNum-lftN-rgtN-ni>=0).AND.(totNum-lftN-rgtN-ni<=maxFilling) &
							.AND.(nl==totNum-lftN-rgtN-ni)) THEN

								ThetaTemp(alpha,k,j,gamma) = &
								ThetaTemp(alpha,k,j,gamma) + &
								Op2((k-1)*localSize+j, (i-1)*localSize+l) * &
								Theta(alpha,i,l,gamma)
							ELSE
							END IF
						END DO
						END DO
					ELSE
					END IF
				END DO
			END DO
			END DO
		END DO

!Internal degrees of freedom absent	
ELSE		
		
		DO alpha=1,chi
			DO k=1,localSize
				DO gamma=1,chi
! Number of particles on the left side of 'link'.
					lftN=LabelLeft(link)%vi(alpha)
! Number of particles on the right side of 'link+2'.
					rgtN=LabelRight(link+2)%vi(gamma)
! So that corresponding lambda is not equal to zero.
					IF((MAX(lftN,rgtN)<1000).AND. &
! 0<=k<=maxfilling
					(totNum-lftN-rgtN-k+1>=0).AND.(totNum-lftN-rgtN-k+1<=maxFilling)) THEN
						ThetaTemp(alpha,k,totNum-lftN-rgtN-k+2,gamma)=CMPLX(0.0,KIND=rKind)
						DO i=1,localSize
! 0<=i<=maxfilling
							IF((totNum-lftN-rgtN-i+1>=0).AND.(totNum-lftN-rgtN-i+1<=maxFilling)) THEN
! Due to number conservation, l and j have to be equal to totNum-lftN-rgtN-k+1 and totNum-lftN-rgtN-i+1.
								ThetaTemp(alpha,k,totNum-lftN-rgtN-k+2,gamma) = &
								ThetaTemp(alpha,k,totNum-lftN-rgtN-k+2,gamma) + &
								Op2((k-1)*localSize+totNum-lftN-rgtN-k+2, (i-1)*localSize+totNum-lftN-rgtN-i+2) * &
								Theta(alpha,i,totNum-lftN-rgtN-i+2,gamma)
							ELSE
							END IF
						END DO
					ELSE
					END IF
				END DO
			END DO
		END DO
		
END IF
		Theta=ThetaTemp
END SUBROUTINE ThetaOperationNC_r


SUBROUTINE ThetaOperationNC_c(link, Op2, Theta, LabelLeft, LabelRight, intDegFree)
!
!Purpose:Operate on Theta with the complex two-site operation Op2 consistent with number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, i,j,l, k, alpha, gamma, lftN, rgtN,nk,nj,nl,ni
		
	chi=SIZE(Theta,1)
	ThetaTemp=CMPLX(0.0,KIND=rKind)
		
!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
		DO alpha=1,chi
			DO k=1,localSize
			nk=Conserv%vi(k)
			DO j=1, localSize
			nj=Conserv%vi(j)
				DO gamma=1,chi
! Number of particles on the left side of 'link'.
					lftN=LabelLeft(link)%vi(alpha)
! Number of particles on the right side of 'link+2'.
					rgtN=LabelRight(link+2)%vi(gamma)
! So that corresponding lambda is not equal to zero.
					IF((MAX(lftN,rgtN)<1000).AND. &
! 0<=n_k<=maxFilling
					(totNum-lftN-rgtN-nk>=0).AND.(totNum-lftN-rgtN-nk<=maxFilling) &
					.AND.(nj==totNum-lftN-rgtN-nk)) THEN
						ThetaTemp(alpha,k,j,gamma)=CMPLX(0.0,KIND=rKind)
						DO i=1,localSize
						ni=Conserv%vi(i)
						DO l=1,localSize
						nl=Conserv%vi(l)
! 0<=n_i<=maxfilling
							IF((totNum-lftN-rgtN-ni>=0).AND.(totNum-lftN-rgtN-ni<=maxFilling) &
							.AND.(nl==totNum-lftN-rgtN-ni)) THEN

								ThetaTemp(alpha,k,j,gamma) = &
								ThetaTemp(alpha,k,j,gamma) + &
								Op2((k-1)*localSize+j, (i-1)*localSize+l) * &
								Theta(alpha,i,l,gamma)
							ELSE
							END IF
						END DO
						END DO
					ELSE
					END IF
				END DO
			END DO
			END DO
		END DO

!Internal degrees of freedom absent	
ELSE		
		
		DO alpha=1,chi
			DO k=1,localSize
				DO gamma=1,chi
! Number of particles on the left side of 'link'.
					lftN=LabelLeft(link)%vi(alpha)
! Number of particles on the right side of 'link+2'.
					rgtN=LabelRight(link+2)%vi(gamma)
! So that corresponding lambda is not equal to zero.
					IF((MAX(lftN,rgtN)<1000).AND. &
! 0<=k<=maxfilling
					(totNum-lftN-rgtN-k+1>=0).AND.(totNum-lftN-rgtN-k+1<=maxFilling)) THEN
						ThetaTemp(alpha,k,totNum-lftN-rgtN-k+2,gamma)=CMPLX(0.0,KIND=rKind)
						DO i=1,localSize
! 0<=i<=maxfilling
							IF((totNum-lftN-rgtN-i+1>=0).AND.(totNum-lftN-rgtN-i+1<=maxFilling)) THEN
! Due to number conservation, l and j have to be equal to totNum-lftN-rgtN-k+1 and totNum-lftN-rgtN-i+1.
								ThetaTemp(alpha,k,totNum-lftN-rgtN-k+2,gamma) = &
								ThetaTemp(alpha,k,totNum-lftN-rgtN-k+2,gamma) + &
								Op2((k-1)*localSize+totNum-lftN-rgtN-k+2, (i-1)*localSize+totNum-lftN-rgtN-i+2) * &
								Theta(alpha,i,totNum-lftN-rgtN-i+2,gamma)
							ELSE
							END IF
						END DO
					ELSE
					END IF
				END DO
			END DO
		END DO
		
END IF
		Theta=ThetaTemp
END SUBROUTINE ThetaOperationNC_c
	
SUBROUTINE RenormThetaNC(link, Theta, LabelLeft, LabelRight, intDegFree)
!
!Purpose:Renormalize Theta so that \sum_{alpha,i,j,gamma}|Theta(alpha,i,j,gamma)|^2 is unity.
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: normSQcmplx
REAL(KIND=rKind) :: norm
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, alpha, gamma, i,j,  lftN, rgtN,ni,nj
		chi=SIZE(Theta,1)
		
		normSQcmplx=CMPLX(0.0,KIND=rKind)

!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
		DO alpha=1,chi
			DO gamma=1,chi
				lftN=LabelLeft(link)%vi(alpha)
				rgtN=LabelRight(link+2)%vi(gamma)
				DO i=1,localSize
				ni=Conserv%vi(i)
				DO j=1,localSize
				nj=Conserv%vi(j)
					IF((MAX(lftN,rgtN)<1000).AND. &
					(totNum-lftN-rgtN-ni>=0).AND.(totNum-lftN-rgtN-ni<=maxFilling) &
					.AND.(nj==totNum-lftN-rgtN-ni)) THEN
					normSQcmplx=normSQcmplx+CONJG(Theta(alpha,i,j,gamma)) * &
					Theta(alpha,i,j,gamma)
					ELSE
					END IF
				END DO
				END DO
			END DO
		END DO


		norm=SQRT(REAL(normSQcmplx,KIND=rKind))
		DO alpha=1,chi
			DO gamma=1,chi
				lftN=LabelLeft(link)%vi(alpha)
				rgtN=LabelRight(link+2)%vi(gamma)
				DO i=1,localSize
				ni=Conserv%vi(i)
				DO j=1,localSize
				nj=Conserv%vi(j)
					IF((MAX(lftN,rgtN)<1000).AND. &
					(totNum-lftN-rgtN-ni>=0).AND.(totNum-lftN-rgtN-ni<=maxFilling) &
					.AND.(nj==totNum-lftN-rgtN-ni)) THEN
					Theta(alpha,i,j,gamma) = &
					Theta(alpha,i,j,gamma)/norm
					ELSE
					END IF
				END DO
				END DO
			END DO
		END DO


		
!Internal degrees of freedom absent	
ELSE

		DO alpha=1,chi
			DO gamma=1,chi
				lftN=LabelLeft(link)%vi(alpha)
				rgtN=LabelRight(link+2)%vi(gamma)
				DO i=1,localSize
					IF((MAX(lftN,rgtN)<1000).AND. &
					(totNum-lftN-rgtN-i+1>=0).AND.(totNum-lftN-rgtN-i+1<=maxFilling)) THEN
					normSQcmplx=normSQcmplx+CONJG(Theta(alpha,i,totNum-lftN-rgtN-i+2,gamma)) * &
					Theta(alpha,i,totNum-lftN-rgtN-i+2,gamma)
					ELSE
					END IF
				END DO
			END DO
		END DO
		norm=SQRT(REAL(normSQcmplx,KIND=rKind))
		DO alpha=1,chi
			DO gamma=1,chi
				lftN=LabelLeft(link)%vi(alpha)
				rgtN=LabelRight(link+2)%vi(gamma)
				DO i=1,localSize
					IF((MAX(lftN,rgtN)<1000).AND. &
					(totNum-lftN-rgtN-i+1>=0).AND.(totNum-lftN-rgtN-i+1<=maxFilling)) THEN
					Theta(alpha,i,totNum-lftN-rgtN-i+2,gamma) = &
					Theta(alpha,i,totNum-lftN-rgtN-i+2,gamma)/norm
					ELSE
					END IF
				END DO
			END DO
		END DO
END IF
END SUBROUTINE RenormThetaNC
	
SUBROUTINE SwapThetaNC(link, Theta, LabelLeft, LabelRight, intDegFree)
!
!Purpose:Swap the two local indices on Theta consistent with number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
INTEGER, INTENT(IN) :: link
COMPLEX(KIND=rKind) :: SwappedTheta(SIZE(Theta,1),SIZE(Theta,2),SIZE(Theta,3),SIZE(Theta,4))
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, alpha, i, j, gamma, lftN, rgtN,ni,nj
		chi=SIZE(Theta,1)


!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
		DO alpha=1,chi
			DO gamma=1,chi
				lftN=LabelLeft(link)%vi(alpha)
				rgtN=LabelRight(link+2)%vi(gamma)
				DO i=1,localSize
				ni=Conserv%vi(i)
				DO j=1,localSize
				nj=Conserv%vi(j)
					IF((MAX(lftN,rgtN)<1000).AND. &
					(totNum-lftN-rgtN-ni>=0).AND.(totNum-lftN-rgtN-ni<=maxFilling) &
					.AND.(nj==totNum-lftN-rgtN-ni)) THEN
						SwappedTheta(alpha,j,i,gamma) = &
						Theta(alpha,i,j,gamma)
					ELSE
					END IF
				END DO
				END DO
			END DO
		END DO
!Internal degrees of freedom absent	
ELSE
		SwappedTheta=CMPLX(0.0,KIND=rKind)
		DO alpha=1,chi
			DO gamma=1,chi
				lftN=LabelLeft(link)%vi(alpha)
				rgtN=LabelRight(link+2)%vi(gamma)
				DO i=1,localSize
					IF((MAX(lftN,rgtN)<1000).AND. &
					(totNum-lftN-rgtN-i+1>=0).AND.(totNum-lftN-rgtN-i+1<=maxFilling)) THEN
						SwappedTheta(alpha,totNum-lftN-rgtN-i+2,i,gamma) = &
						Theta(alpha,i,totNum-lftN-rgtN-i+2,gamma)
					ELSE
					END IF
				END DO
			END DO
		END DO
END IF
		Theta=SwappedTheta
END SUBROUTINE SwapThetaNC
	
SUBROUTINE minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
!
!Purpose:Find the minimum and maximum values of the number on the left=N_L(\alpha)+N_S(i) and on the right=N_R(\alpha )+N_S(j)
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
INTEGER, INTENT(OUT) :: minNL, maxNL, minNR, maxNR
INTEGER :: chi, alpha, gamma
INTEGER :: Label0(SIZE(LabelLeft(link)%vi,1))
		
		chi=SIZE(LabelLeft(link)%vi,1)
		!as few as 0 particles on site i
		minNL=MINVAL(LabelLeft(link)%vi)
		DO alpha=1,chi
		!Label=10000 means the state is notinitialized (not enough entanglement to need it)
			IF(LabelLeft(link)%vi(alpha)==10000) THEN
				Label0(alpha)=0
			ELSE
				Label0(alpha)=LabelLeft(link)%vi(alpha)
			END IF
		END DO
		!up to maxfilling particles allowed on site i
		maxNL=MAXVAL(Label0)+maxFilling
		!as few as 0 particles on site j
		minNR=MINVAL(LabelRight(link+2)%vi)
		DO gamma=1,chi
		!Label=10000 means the state is notinitialized (not enough entanglement to need it)
			IF(LabelRight(link+2)%vi(gamma)==10000) THEN
				Label0(gamma)=0
			ELSE
				Label0(gamma)=LabelRight(link+2)%vi(gamma)
			END IF
		END DO
		!up to maxfilling particles allowed on site j
		maxNR=MAXVAL(Label0)+maxFilling
		
		!If the onsite maximum number puts us over the total number, make the total number our bound
		IF(maxNL>totNum-minNR) THEN
			maxNL=totNum-minNR
		ELSE
			minNR=totNum-maxNL
		END IF
		
		IF(maxNR>totNum-minNL) THEN
			maxNR=totNum-minNL
		ELSE
			minNL=totNum-maxNR
		END IF
END SUBROUTINE minmaxNLR
	
SUBROUTINE SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight, intDegFree)
!
!Purpose:The reshaped Theta will be block diagonal due to number conservation, 
!with each block corresponding to a fixed number on the left.  This subroutine finds how many
!blocks exist and their sizes
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(OUT) :: BlockSize(:,:)
INTEGER, INTENT(IN) :: link, minNL, maxNL, minNR, maxNR
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, alpha, gamma, i, j, k, NLeft, NRight, ni,nj
		chi=SIZE(LabelLeft(link)%vi,1)
		BlockSize=0
		!least number on the left corresponds to index 1
		!Blocksize(i,1)=(number on the left)_i
		DO k=minNL,maxNL,1
			BlockSize(k-minNL+1,1)=k
		END DO

!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
!Count the size of each block (of fixed number on the left) on the left		
		DO alpha=1,chi
			DO i=1,localSize
			ni=Conserv%vi(i)
				If(LabelLeft(link)%vi(alpha)<1000) THEN
					NLeft=LabelLeft(link)%vi(alpha)+ni
					IF((minNL<=NLeft).AND.(NLeft<=maxNL)) THEN
						BlockSize(Nleft-minNL+1,2)=BlockSize(Nleft-minNL+1,2)+1
					ELSE
					END IF
				ELSE		
				END IF					
			END DO
		END DO

!Count the size of each block (of fixed number on the left) on the right		
		DO gamma=1,chi
			DO j=1,localSize
			nj=Conserv%vi(j)
				IF(LabelRight(link+2)%vi(gamma)<1000) THEN
					NRight=LabelRight(link+2)%vi(gamma)+nj
					IF((minNR<=NRight).AND.(NRight<=maxNR)) THEN
						BlockSize(totNum-NRight-minNL+1,3)=BlockSize(totNum-NRight-minNL+1,3)+1
					ELSE
					END IF
				ELSE
				END IF
			END DO
		END DO

!Internal degrees of freedom absent	
ELSE

!Count the size of each block (of fixed number on the left) on the left				
		DO alpha=1,chi
			DO i=1,localSize
				If(LabelLeft(link)%vi(alpha)<1000) THEN
					NLeft=LabelLeft(link)%vi(alpha)+i-1
					IF((minNL<=NLeft).AND.(NLeft<=maxNL)) THEN
						BlockSize(Nleft-minNL+1,2)=BlockSize(Nleft-minNL+1,2)+1
					ELSE
					END IF
				ELSE		
				END IF					
			END DO
		END DO
!Count the size of each block (of fixed number on the left) on the right		
		DO gamma=1,chi
			DO j=1,localSize
				IF(LabelRight(link+2)%vi(gamma)<1000) THEN
					NRight=LabelRight(link+2)%vi(gamma)+j-1
					IF((minNR<=NRight).AND.(NRight<=maxNR)) THEN
						BlockSize(totNum-NRight-minNL+1,3)=BlockSize(totNum-NRight-minNL+1,3)+1
					ELSE
					END IF
				ELSE
				END IF
			END DO
		END DO
		
END IF
END SUBROUTINE SizeOfBlocks
	
	
SUBROUTINE IndexLeft(link, indL, minNL, maxNL, LabelLeft, intDegFree)
!
!Purpose: Find the on-site indices i and left schmidt indices alpha that 
!correspond to the fixed number on the left indexed by indind
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:)
TYPE(matrixInt), POINTER :: indL(:)
INTEGER, INTENT(IN) :: link, minNL, maxNL
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, i, alpha, indind,ni
INTEGER :: turnL(maxNL-minNL+1)
	chi=SIZE(LabelLeft(link)%vi,1)
	turnL=1

!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
		DO alpha=1,chi
			DO i=1,localSize
			ni=Conserv%vi(i)
!Check that the number on the left is an allowed value
				IF((LabelLeft(link)%vi(alpha)<1000).AND.(minNL<=LabelLeft(link)%vi(alpha)+ni) &
					.AND.(LabelLeft(link)%vi(alpha)+ni<=maxNL)) THEN
!This index is the number
				indind=LabelLeft(link)%vi(alpha)+ni-minNL+1
!Index of the schmidt vector
				indL(indind)%mi(turnL(indind),1)=alpha
!Index of the on-site state
				indL(indind)%mi(turnL(indind),2)=i
!Next allowed state index
				turnL(indind)=turnL(indind)+1
				ELSE
				END IF
			END DO
		END DO
!Internal degrees of freedom absent
ELSE	

	DO alpha=1,chi
		DO i=1,localSize
			IF((LabelLeft(link)%vi(alpha)<1000).AND.(minNL<=LabelLeft(link)%vi(alpha)+i-1) &
				.AND.(LabelLeft(link)%vi(alpha)+i-1<=maxNL)) THEN
			indind=LabelLeft(link)%vi(alpha)+i-1-minNL+1
			indL(indind)%mi(turnL(indind),1)=alpha
			indL(indind)%mi(turnL(indind),2)=i
			turnL(indind)=turnL(indind)+1
			ELSE
			END IF
		END DO
	END DO
END IF
END SUBROUTINE IndexLeft
	
SUBROUTINE IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,intDegFree)
!
!Purpose: Find the on-site indices j and right schmidt gamma alpha that 
!correspond to the fixed number on the left indexed by indind
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelRight(:)
TYPE(matrixInt), POINTER :: indR(:)
INTEGER, INTENT(IN) :: link, minNL, maxNL, minNR, maxNR
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: chi, j, gamma, indind,nj
INTEGER :: turnR(maxNL-minNL+1)
	chi=SIZE(LabelRight(link+2)%vi,1)
	turnR=1

!Internal degrees of freedom present	
IF(PRESENT(intDegFree)) THEN
		DO gamma=1,chi
			DO j=1,localSize
			nj=Conserv%vi(j)
!Same as above, now for the right
				IF((LabelRight(link+2)%vi(gamma)<1000).AND.(minNR<=LabelRight(link+2)%vi(gamma)+nj) &
					.AND.(LabelRight(link+2)%vi(gamma)+nj<=maxNR)) THEN
				indind=totNum-(LabelRight(link+2)%vi(gamma)+nj)-minNL+1
				indR(indind)%mi(turnR(indind),1)=gamma
				indR(indind)%mi(turnR(indind),2)=j
				turnR(indind)=turnR(indind)+1
				ELSE
				END IF
			END DO
		END DO
!Internal degrees of freedom absent
ELSE

	DO gamma=1,chi
		DO j=1,localSize
			IF((LabelRight(link+2)%vi(gamma)<1000).AND.(minNR<=LabelRight(link+2)%vi(gamma)+j-1) &
				.AND.(LabelRight(link+2)%vi(gamma)+j-1<=maxNR)) THEN
				!totNum-number on right
			indind=totNum-(LabelRight(link+2)%vi(gamma)+j-1)-minNL+1
			!indind is block number, turnR(indind) is the position in the block
			indR(indind)%mi(turnR(indind),1)=gamma
			indR(indind)%mi(turnR(indind),2)=j
			turnR(indind)=turnR(indind)+1
			ELSE
			END IF
		END DO
	END DO
END IF
END SUBROUTINE IndexRight
	
SUBROUTINE FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
!
!Purpose: Form the NumBlocks Blocks that make up the full theta
!
IMPLICIT NONE
TYPE(matrix), POINTER :: BlockTheta(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Theta(:,:,:,:)
INTEGER :: chi, NumOfBlocks, k, falp, fgamm, lgt, trv
	chi=SIZE(Theta,1)
	NumOfBlocks=SIZE(BlockTheta,1)
		
	DO k=1,NumOfBlocks
		lgt=BlockSize(k,2)
		trv=BlockSize(k,3)
		DO falp=1,lgt
			DO fgamm=1,trv
				BlockTheta(k)%m(falp,fgamm)= &
				Theta(indL(k)%mi(falp,1),indL(k)%mi(falp,2),indR(k)%mi(fgamm,2),indR(k)%mi(fgamm,1))
			END DO
		END DO
	END DO
END SUBROUTINE FormBlockTheta
	
SUBROUTINE SVDNC(US, SS, VS, BlockTheta, BlockSize)
!
!Purpose: Perform an SVD on each one of the Blocks from FormBlockTheta
!
IMPLICIT NONE
TYPE(matrix), POINTER :: US(:), VS(:), BlockTheta(:)
TYPE(vector), POINTER :: SS(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, NumOfBlocks
	NumOfBlocks=SIZE(BlockTheta,1)
		
	DO k=1,NumOfBlocks
	
	!Allocate SVD variables for the block
		CALL SVDInitNC(k, BlockSize)
	!Perform SVD on the block
		CALL ZGESVD(jobu_SVD, jobvt_SVD, matrixSizeL_SVD, matrixSizeT_SVD, BlockTheta(k)%m, &
		matrixSizeL_SVD, SS(k)%v, US(k)%m, matrixSizeL_SVD, VS(k)%m, matrixSizeT_SVD, &
		work_SVD, workSize_SVD, rwork_SVD, info_SVD)
		CALL SVDFinishNC()
	END DO
	
END SUBROUTINE SVDNC
	
SUBROUTINE FlattenSS(SS, ssfl, BlockSize)
!
!Purpose: Compile Singular values from all blocks into one vector ssfl
!
IMPLICIT NONE
TYPE(vector), POINTER :: SS(:)
REAL(KIND=rKind), INTENT(OUT) :: ssfl(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: k, fbeta, NumOfBlocks, SizeOfSS, turn
	NumOfBlocks=SIZE(BlockSize,1)
	SizeOfSS=SIZE(ssfl,1)
	turn=1
	DO k=1,NumOfBlocks
		!There are min(m,n) singular values of an mXn matrix
		DO fbeta=1,MIN(BlockSize(k,2),BlockSize(k,3))
		ssfl(turn)=SS(k)%v(fbeta)
		turn=turn+1
		END DO
	END DO  
END SUBROUTINE FlattenSS
	
SUBROUTINE Ordering(RealArray, order)
!
!Purpose: Create an array of the indices of the singular values from greatest to least
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: RealArray(:)
INTEGER, INTENT(OUT) :: order(:)
INTEGER :: SizeOfSS, beta, MaxElement
REAL(KIND=rKind) :: DeformedArray(SIZE(RealArray,1))
INTEGER :: MINLOC_array(1)
	order=0
	SizeOfSS=SIZE(RealArray,1)
	MaxElement=MAXVAL(RealArray)
	DeformedArray=RealArray
	!Work backwards
	DO beta=SizeOfSS,1,-1
		!Find position of minimum value
		MINLOC_array=MINLOC(DeformedArray)
		!Place that in the next available position
		order(beta)=MINLOC_array(1)
		!Remove that element from consideration
		DeformedArray(MINLOC_array(1))=MaxElement+10.0
	END DO
END SUBROUTINE Ordering
	
SUBROUTINE JudgePosition(Position, order, BlockSize)
!
!Purpose: Find the "inverse map" to ordering above i.e. find the block index and
! the position within the block for each value in the new singular value ordering
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(OUT) :: Position(:,:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: BlockSize(:,:)
INTEGER :: ordertemp(SIZE(order,1))
INTEGER :: k, beta, NumOfBlocks, sz
	Position=0
	sz=SIZE(order,1)
	NumOfBlocks=SIZE(BlockSize,1)
	ordertemp=order
	DO beta=1,sz
		DO k=1,NumOfBlocks
			IF(ordertemp(beta)<=0) EXIT
			ordertemp(beta)=ordertemp(beta)-MIN(BlockSize(k,2),BlockSize(k,3))
			Position(beta,1)=k
			Position(beta,2)=ordertemp(beta)+MIN(BlockSize(k,2),BlockSize(k,3))
		END DO
	END DO
END SUBROUTINE JudgePosition
	
SUBROUTINE UpdateLabelLeft(link, LabelLeft, minNL, Position, ssfl, order)
!
!Purpose: Update LabelLeft based on the new signular value ordering
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:)
INTEGER, INTENT(IN) :: link, minNL
INTEGER, INTENT(IN) :: Position(:,:)
REAL(KIND=rKind), INTENT(IN) :: ssfl(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER :: order2(SIZE(order,1))
INTEGER :: chi, beta, sz, l1
	sz=SIZE(order,1)
	chi=SIZE(LabelLeft(link)%vi,1)
	DO l1=1,sz
		order2(order(l1))=l1
	END DO
	DO beta=1,chi
		IF(beta<=sz) THEN
			!Truncate Schmidt index at chi
			IF(order2(order(beta))<=chi) THEN
				!Keep only the singular values greater than 10**(-10)
				IF(ssfl(order(beta))>10.0**(-10)) THEN
					LabelLeft(link+1)%vi(beta)=minNL+Position(beta,1)-1
				ELSE
					!If the state is simply unused, make LabelLeft 10000
					LabelLeft(link+1)%vi(beta)=10000
				END IF
			ELSE
			END IF
		ELSE
			!If the state is simply unused, make LabelLeft 10000
			LabelLeft(link+1)%vi(beta)=10000
		END IF

	END DO	
END SUBROUTINE UpdateLabelLeft
	
SUBROUTINE UpdateLabelRight(link, LabelLeft, LabelRight)
!
!Purpose: Update LabelRight based on the new signular value ordering
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
INTEGER :: chi, beta
	chi=SIZE(LabelLeft(link+1)%vi,1)
	DO beta=1,chi
		IF(LabelLeft(link+1)%vi(beta)<1000) THEN
			LabelRight(link+1)%vi(beta)=totNum-LabelLeft(link+1)%vi(beta)
		ELSE
			!If the state is simply unused, make LabelRight 10000
			LabelRight(link+1)%vi(beta)=10000
		END IF

	END DO
END SUBROUTINE UpdateLabelRight
	
SUBROUTINE FormLambdaNC(Lambda1, truncerr, ssfl, order)
!
!Purpose: Form Lambdas from the ordered singular values
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(INOUT) :: Lambda1(:)
REAL(KIND=rKind), INTENT(OUT) :: truncerr
REAL(KIND=rKind), INTENT(IN) :: ssfl(:)
INTEGER, INTENT(IN) :: order(:)
REAL(KIND=rKind) :: tempLam(SIZE(Lambda1,1))
INTEGER :: chi, beta, sz
REAL(KIND=rKind) :: renorm
	chi=SIZE(Lambda1,1)
	sz=SIZE(order,1)
	tempLam=Lambda1
	DO beta=1,chi
		IF((beta<=sz).AND.(ssfl(order(beta))>10.0**(-10))) THEN
			tempLam(beta)=ssfl(order(beta))
		ELSE
			tempLam(beta)=0.0
		END IF
	END DO
	renorm = 0.0
	DO beta=1,chi
		renorm=renorm+tempLam(beta)**(2)
	END DO
	truncerr = ABS(1.0_rKind - renorm)
	DO beta=1,chi
		Lambda1(beta)=tempLam(beta)/SQRT(renorm)
	END DO
END SUBROUTINE FormLambdaNC
	
SUBROUTINE FormGamma1NC(Lambda0, Gamma1, US, indL, order, Position, BlockSize)
!
!Purpose: Form Left Gamma from the left SVD matrix
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma1(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda0(:)
TYPE(matrix), POINTER :: US(:)
TYPE(matrixInt), POINTER :: indL(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: Position(:,:), BlockSize(:,:)
INTEGER :: sz, chi, falp, beta, bn, psn
	sz=SIZE(order,1)
	chi=SIZE(Lambda0,1)
	Gamma1=CMPLX(0.0,KIND=rKind)
	DO beta=1,chi
		IF(beta<=sz) THEN
			bn=Position(beta,1)
			psn=Position(beta,2)
			DO falp=1,BlockSize(bn,2)
				IF(Lambda0(indL(bn)%mi(falp,1))>10.0**(-10)) THEN
					Gamma1(indL(bn)%mi(falp,1),indL(bn)%mi(falp,2),beta)= &
					US(bn)%m(falp,psn)/Lambda0(indL(bn)%mi(falp,1))
				ELSE
					Gamma1(indL(bn)%mi(falp,1),indL(bn)%mi(falp,2),beta)=CMPLX(0.0,KIND=rKind)
				END IF
			END DO
		ELSE
		END IF
	END DO
END SUBROUTINE FormGamma1NC

SUBROUTINE FormGamma2NC(Gamma2, Lambda2, VS, indR, order, Position, BlockSize)
!
!Purpose: Form right Gamma from the right SVD matrix
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Gamma2(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda2(:)
TYPE(matrix), POINTER :: VS(:)
TYPE(matrixInt), POINTER :: indR(:)
INTEGER, INTENT(IN) :: order(:)
INTEGER, INTENT(IN) :: Position(:,:), BlockSize(:,:)
INTEGER :: sz, chi, beta, fgamm, bn, psn
	sz=SIZE(order,1)
	chi=SIZE(Lambda2,1)
	Gamma2=CMPLX(0.0,KIND=rKind)
	DO beta=1,chi
		IF(beta<=sz) THEN
			bn=Position(beta,1)
			psn=Position(beta,2)
			DO fgamm=1,BlockSize(bn,3)
				IF(Lambda2(indR(bn)%mi(fgamm,1))>10.0**(-10)) THEN
					Gamma2(beta,indR(bn)%mi(fgamm,2),indR(bn)%mi(fgamm,1))= &
					VS(bn)%m(psn,fgamm)/Lambda2(indR(bn)%mi(fgamm,1))
				ELSE
					Gamma2(beta,indR(bn)%mi(fgamm,2),indR(bn)%mi(fgamm,1))=CMPLX(0.0,KIND=rKind)
				END IF
			END DO
		ELSE
		END IF
	END DO
END SUBROUTINE FormGamma2NC
	
SUBROUTINE TwoSiteOpNC_r(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr,intDegFree)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
					SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: minNL, maxNL, minNR, maxNR

!Internal degrees of freedom present			
IF(PRESENT(intDegFree)) THEN

	CALL FormThetaNC(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight,intDegFree)
	CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight,intDegFree)
!		CALL TruncationForNC(link, Theta, LabelLeft, LabelRight,intDegFree)
	CALL RenormThetaNC(link, Theta, LabelLeft, LabelRight,intDegFree)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
			END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight,intDegFree)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft,intDegFree)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,intDegFree)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
			END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
			END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL UpdateLabelLeft(link, LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	CALL FormLambdaNC(Lambdas(link+1)%v, truncerr, ssflat%v, ord%vi)
	CALL FormGamma1NC(Lambdas(link)%v, Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
	CALL FormGamma2NC(Gammas(link+1)%t,Lambdas(link+2)%v, VS, indR, ord%vi, Position, BlockSize)
!Internal degrees of freedom absent				
ELSE
	CALL FormThetaNC(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
	CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight)
!		CALL TruncationForNC(link, Theta, LabelLeft, LabelRight)
	CALL RenormThetaNC(link, Theta, LabelLeft, LabelRight)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
			END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
			END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
			END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL UpdateLabelLeft(link, LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	CALL FormLambdaNC(Lambdas(link+1)%v, truncerr, ssflat%v, ord%vi)
	CALL FormGamma1NC(Lambdas(link)%v, Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
	CALL FormGamma2NC(Gammas(link+1)%t,Lambdas(link+2)%v, VS, indR, ord%vi, Position, BlockSize)
END IF
		
	DEALLOCATE(BlockSize, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
			END IF
	CALL DeallocateIndexLR(indL,indR)
	CALL DeallocateBlockTheta(BlockTheta)
	CALL DeallocateUSV(US, SS, VS)
	DEALLOCATE(ssflat%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
			END IF
	DEALLOCATE(ord%vi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
			END IF
	DEALLOCATE(Position, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
			END IF
END SUBROUTINE TwoSiteOpNC_r


SUBROUTINE TwoSiteOpNC_c(link, Op2, Gammas, Lambdas, LabelLeft, LabelRight, truncerr,intDegFree)
!
!Purpose: Perform the two-site operation Op2 on the sites divided by link
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: link
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
REAL(KIND=rKind), INTENT(INOUT) :: truncerr
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
					SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: minNL, maxNL, minNR, maxNR

!Internal degrees of freedom present			
IF(PRESENT(intDegFree)) THEN

	CALL FormThetaNC(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight,intDegFree)
	CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight,intDegFree)
!		CALL TruncationForNC(link, Theta, LabelLeft, LabelRight,intDegFree)
	CALL RenormThetaNC(link, Theta, LabelLeft, LabelRight,intDegFree)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
			END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight,intDegFree)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft,intDegFree)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,intDegFree)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
			END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
			END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL UpdateLabelLeft(link, LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	CALL FormLambdaNC(Lambdas(link+1)%v, truncerr, ssflat%v, ord%vi)
	CALL FormGamma1NC(Lambdas(link)%v, Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
	CALL FormGamma2NC(Gammas(link+1)%t,Lambdas(link+2)%v, VS, indR, ord%vi, Position, BlockSize)
!Internal degrees of freedom absent				
ELSE
	CALL FormThetaNC(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
	CALL ThetaOperationNC(link, Op2, Theta, LabelLeft, LabelRight)
!		CALL TruncationForNC(link, Theta, LabelLeft, LabelRight)
	CALL RenormThetaNC(link, Theta, LabelLeft, LabelRight)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockSize in TwoSiteOpNC'
			END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ord in TwoSiteOpNC'
			END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Position in TwoSiteOpNC'
			END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL UpdateLabelLeft(link, LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	CALL FormLambdaNC(Lambdas(link+1)%v, truncerr, ssflat%v, ord%vi)
	CALL FormGamma1NC(Lambdas(link)%v, Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
	CALL FormGamma2NC(Gammas(link+1)%t,Lambdas(link+2)%v, VS, indR, ord%vi, Position, BlockSize)
END IF
		
	DEALLOCATE(BlockSize, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate BlockSize in TwoSiteOpNC'
			END IF
	CALL DeallocateIndexLR(indL,indR)
	CALL DeallocateBlockTheta(BlockTheta)
	CALL DeallocateUSV(US, SS, VS)
	DEALLOCATE(ssflat%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate ssflat in TwoSiteOpNC'
			END IF
	DEALLOCATE(ord%vi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate ord in TwoSiteOpNC'
			END IF
	DEALLOCATE(Position, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate position in TwoSiteOpNC'
			END IF
END SUBROUTINE TwoSiteOpNC_c
	
SUBROUTINE SwappingNC(link, Gammas, Lambdas, LabelLeft, LabelRight, truncerr,intDegFree)
!
!Purpose: Swap the sites divided by link, then reindex to match
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(vector), POINTER :: Lambdas(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: link
TYPE(vector), POINTER :: SS(:)
TYPE(vector) :: ssflat
TYPE(vectorInt) :: ord
TYPE(matrix), POINTER :: BlockTheta(:), US(:), VS(:)
TYPE(matrixInt), POINTER :: indL(:), indR(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(link)%t,1),SIZE(Gammas(link)%t,2), &
						SIZE(Gammas(link+1)%t,2),SIZE(Gammas(link+1)%t,3))
INTEGER, DIMENSION(:,:), ALLOCATABLE :: BlockSize, Position
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: minNL, maxNL, minNR, maxNR
REAL(KIND=rKind) :: truncerr

!Internal degrees of freedom present		
IF(PRESENT(intDegFree)) THEN		
				
	CALL FormThetaNC(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight,intDegFree)
	CALL SwapThetaNC(link, Theta, LabelLeft, LabelRight,intDegFree)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockSize in SwappingNC'
			END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight,intDegFree)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft,intDegFree)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight,intDegFree)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ord in SwappingNC'
			END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Position in SwappingNC'
			END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL UpdateLabelLeft(link, LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	CALL FormLambdaNC(Lambdas(link+1)%v, truncerr, ssflat%v, ord%vi)
	CALL FormGamma1NC(Lambdas(link)%v, Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
	CALL FormGamma2NC(Gammas(link+1)%t,Lambdas(link+2)%v, VS, indR, ord%vi, Position, BlockSize)

!Internal degrees of freedom absent
ELSE

	CALL FormThetaNC(link, Theta, Gammas, Lambdas, LabelLeft, LabelRight)
	CALL SwapThetaNC(link, Theta, LabelLeft, LabelRight)
	CALL minmaxNLR(link, LabelLeft, LabelRight, minNL, maxNL, minNR, maxNR)
	ALLOCATE(BlockSize(maxNL-minNL+1,3), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate BlockSize in SwappingNC'
			END IF
	CALL SizeOfBlocks(link, BlockSize, minNL, maxNL, minNR, maxNR, LabelLeft, LabelRight)
	CALL AllocateIndexLR(indL, indR , BlockSize)
	CALL IndexLeft(link, indL, minNL, maxNL, LabelLeft)
	CALL IndexRight(link, indR, minNL, maxNL, minNR, maxNR, LabelRight)
	CALL AllocateBlockTheta(BlockTheta, BlockSize)
	CALL FormBlockTheta(BlockTheta, indL, indR, BlockSize, Theta)
	CALL AllocateUSV(US, SS, VS, BlockSize)
	CALL SVDNC(US, SS, VS, BlockTheta, BlockSize)
	CALL AllocateSSflat(ssflat, BlockSize)
	CALL FlattenSS(SS, ssflat%v, BlockSize)
	ALLOCATE(ord%vi(SIZE(ssflat%v,1)), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate ord in SwappingNC'
			END IF
	CALL Ordering(ssflat%v, ord%vi)
	ALLOCATE(Position(SIZE(ssflat%v,1),2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate Position in SwappingNC'
			END IF
	CALL JudgePosition(Position, ord%vi, BlockSize)
	CALL UpdateLabelLeft(link, LabelLeft, minNL, Position, ssflat%v, ord%vi)
	CALL UpdateLabelRight(link, LabelLeft, LabelRight)
	CALL FormLambdaNC(Lambdas(link+1)%v, truncerr, ssflat%v, ord%vi)
	CALL FormGamma1NC(Lambdas(link)%v, Gammas(link)%t, US, indL, ord%vi, Position, BlockSize)
	CALL FormGamma2NC(Gammas(link+1)%t,Lambdas(link+2)%v, VS, indR, ord%vi, Position, BlockSize)
END IF

	DEALLOCATE(BlockSize, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate BlockSize in SwappingNC'
			END IF
	CALL DeallocateIndexLR(indL,indR)
	CALL DeallocateBlockTheta(BlockTheta)
	CALL DeallocateUSV(US, SS, VS)
	DEALLOCATE(ssflat%v, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate ssflat in SwappingNC'
			END IF
	DEALLOCATE(ord%vi, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate ord in SwappingNC'
			END IF
	DEALLOCATE(Position, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate Position in SwappingNC'
			END IF
END SUBROUTINE SwappingNC

SUBROUTINE SpecialState(Gammas, Lambdas, stateChar)
!
!Purpose: Construct the Vidal decomposition of a GHZ, W, or cluster state
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: i, j,  chi
CHARACTER(len=*), INTENT(IN) :: stateChar
TYPE(matrix) ::  osOp, expOp, expOph
COMPLEX(KIND=rKind) :: angle
REAL(KIND=rKind) :: trunctemp

IF(stateChar=='W') THEN	
	IF(SIZE(Lambdas(2)%v).lt.2) THEN
		STOP 'Chi must be at least 2 to represent a W state!'
	END IF

	Lambdas(1)%v = 0.0_rKind
	Lambdas(1)%v(1)=1.0_rKind

	DO i = 2, systemSize
		Lambdas(i)%v = 0.0_rKind
		Lambdas(i)%v(1) = SQRT((systemSize-i+1)*1.0_rKind/(systemSize*1.0_rKind))
		Lambdas(i)%v(2) = SQRT((i-1)*1.0_rKind/(systemSize*1.0_rKind))
	END DO

	Lambdas(systemSize+1)%v = 0.0_rKind
	Lambdas(systemSize+1)%v(1)=1.0_rKind

	Gammas(1)%t = CMPLX(0.0, KIND=rKind)
	Gammas(1)%t(1,1,1) = 1.0_rKind
	Gammas(1)%t(1,2,2) = 1.0_rKind

	DO i = 2, systemSize-1
		Gammas(i)%t = CMPLX(0.0, KIND=rKind)
		Gammas(i)%t(1, 1, 1) = SQRT(systemSize*1.0_rKind/((systemSize-i+1)*1.0_rKind))
		Gammas(i)%t(1, 2, 2) = SQRT(systemSize*1.0_rKind/(i*(systemSize-i+1)*1.0_rKind))
		Gammas(i)%t(2, 1, 2) = SQRT(systemSize*1.0_rKind/(i*1.0_rKind))
	END DO

	Gammas(systemSize)%t = CMPLX(0.0, KIND=rKind)
	Gammas(systemSize)%t(1,2,1) = 1.0_rKind
	Gammas(systemSize)%t(2,1,1) = 1.0_rKind


ELSE IF(stateChar=='GHZ') THEN
	IF(SIZE(Lambdas(2)%v).lt.localSize) THEN
		STOP 'Chi must be at least localSize to represent a GHZ state!'
	END IF
	DO i = 2, systemSize-1
	Gammas(i)%t = CMPLX(0.0, KIND=rKind)
	Lambdas(i)%v = 0.0_rKind
		DO j = 1, localSize
		Lambdas(i)%v(j) = 1.0_rKind/SQRT(localSize*1.0_rKind) 
		Gammas(i)%t(j, j, j) =SQRT(localSize*1.0_rKind) 
		END DO
	END DO
	Lambdas(1)%v = 0.0_rKind
	Gammas(1)%t = CMPLX(0.0, KIND=rKind)
	Gammas(systemSize)%t = CMPLX(0.0, KIND=rKind)
	Lambdas(systemSize)%v = 0.0_rKind
	Lambdas(systemSize+1)%v = 0.0_rKind
	DO j = 1, localSize
	Gammas(1)%t(1, j, j) =1.0_rKind 
	Gammas(systemSize)%t(j, j, 1) =1.0_rKind 
	Lambdas(systemSize)%v(j)= 1.0_rKind/SQRT(localSize*1.0_rKind) 
	END DO
	Lambdas(systemSize+1)%v(1)=1.0_rKind
	Lambdas(1)%v(1)=1.0_rKind
ELSE IF(stateChar=='Cluster') THEN
	IF(localSize.ne.2) THEN
		STOP "Only hard-core boson cluster states supported!"
	END IF
	!Begin in the completely "up" polarized state
	DO i=1,systemSize
		Lambdas(i)%v=0.0_rKind
		Gammas(i)%t=CMPLX(0.0, KIND=rKind)
		Lambdas(i)%v(1)=1.0_rKind
		Gammas(i)%t(1,1,1)=1.0_rKind
	END DO

	lambdas(systemSize+1)%v=0.0_rKind
	Lambdas(systemSize+1)%v(1)=1.0_rKind
	!Rotate by an angle \pi/4 wrt the y axis

	ALLOCATE(osOp%m(localSize,localSize), expOp%m(localSize,localSize))
	!Form sigma_y in the Schwinger representation
	osOp%m=-CMPLX(0.0,1.0,KIND=rKind)*(a_op%mr-TRANSPOSE(a_op%mr))

	!Exponentiate
	angle=3.1415926535897_rKind*0.25_rKind
	CALL Matrix_Exponential(osOp%m,expOp%m,angle,localSize)

	!Rotate
	DO i=1,systemSize
		CALL OneSiteOp(expOp%m,Gammas(i)%t)
	END DO
	DEALLOCATE(osOp%m, expOp%m)
	!Evolve under the ising Hamiltonian for t=\pi/4
	!Create the Hamiltonian in the Schwinger representation
	ALLOCATE(osOp%m(localSize*localSize,localSize*localSize), expOp%m(localSize*localSize,localSize*localSize),expOph%m(localSize*localSize,localSize*localSize))
	osOp%m=-TensorProd(one_op%mr-2.0*MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr-2.0*MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
	!Exponentiate with deltat=1/100 the final time
	angle=3.1415926535897_rKind*0.25_rKind*0.01_rKind
	CALL Matrix_Exponential(osOp%m,expOp%m,angle,localSize*localSize)
	angle=3.1415926535897_rKind*0.25_rKind*0.005_rKind
	CALL Matrix_Exponential(osOp%m,expOph%m,angle,localSize*localSize)

	!Time step 100 time steps
	CALL SVDInit(SIZE(Lambdas(FLOOR(systemSize*0.5_rKind))%v))
	DO j=1,100
	!!! Operate Exp(-i Hodd dt/2)	
		DO i=1,(systemSize-1),2
			CALL TwoSiteOp(i, expOph%m, Gammas, Lambdas, trunctemp)
		END DO
	!!! Operate Exp(-i Heven dt)
		DO i=2,(systemSize-1),2
			CALL TwoSiteOp(i, expOp%m, Gammas, Lambdas, trunctemp)
		END DO
	!!! Operate Exp(-i Hodd dt/2)			
		DO i=1,(systemSize-1),2
			CALL TwoSiteOp(i, expOph%m, Gammas, Lambdas, trunctemp )
		END DO
	END DO
	CALL SVDFinish()
	DEALLOCATE(osOp%m, expOp%m, expOph%m)
	!Rotate by an angle \pi/4 wrt the z axis

	ALLOCATE(osOp%m(localSize,localSize), expOp%m(localSize,localSize))
	!Form sigma_z in the Schwinger representation
	osOp%m=one_op%mr-2.0*MATMUL(TRANSPOSE(a_op%mr),a_op%mr)

	!Exponentiate
	angle=3.1415926535897_rKind*0.25_rKind
	CALL Matrix_Exponential(osOp%m,expOp%m,angle,localSize)

	!Rotate
	DO i=1,systemSize
		CALL OneSiteOp(expOp%m,Gammas(i)%t)
	END DO
	DEALLOCATE(osOp%m, expOp%m)

ELSE
	STOP 'Character argument of SpecialState not recognized!'
END IF
END SUBROUTINE SpecialState
	
END MODULE local_operations_module
