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
MODULE observables_module
!
! Purpose: Module to compute observables
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
IMPLICIT NONE


!Type for all collected measures
TYPE measure
	REAL(KIND=rKind) :: en
	LOGICAL :: localpres
	TYPE(mlocal), POINTER :: local(:)
	LOGICAL :: avgpres
	TYPE(mavg), POINTER :: avg(:)
	INTEGER :: entpres
	TYPE(entropy) :: ent
	LOGICAL :: corrpres
	TYPE(mcorr), POINTER :: corr(:)
	LOGICAL :: fermicorrpres
	TYPE(mcorrf), POINTER :: fermicorr(:)
END TYPE measure

!*** INTERFACES
INTERFACE TotalOneSite
MODULE PROCEDURE TotalOneSite_mr, &
				TotalOneSite_m, TotalOneSite_r,&
				TotalOneSite_c
END INTERFACE TotalOneSite

INTERFACE OneSiteExpVal
MODULE PROCEDURE OneSiteExpVal_mr,&
				 OneSiteExpVal_m, OneSiteExpVal_r,&
				 OneSiteExpVal_c
END INTERFACE  OneSiteExpVal

INTERFACE OneSiteVar
MODULE PROCEDURE OneSiteVar_mr,&
				 OneSiteVar_m, OneSiteVar_r,&
				 OneSiteVar_c
END INTERFACE  OneSiteVar

INTERFACE TwoSiteExpVal
MODULE PROCEDURE TwoSiteExpVal_r,&
				 TwoSiteExpVal_c, &
				 TwoSiteExpVal_rc,&
				 TwoSiteExpVal_cr
END INTERFACE  TwoSiteExpVal

INTERFACE TwoSiteExpValG
MODULE PROCEDURE TwoSiteExpValG_r,&
				 TwoSiteExpValG_c, &
				 TwoSiteExpValG_rc,&
				 TwoSiteExpValG_cr
END INTERFACE  TwoSiteExpValG

INTERFACE TwoPointExpValNC
MODULE PROCEDURE TwoPointExpValNC_r,&
				 TwoPointExpValNC_c, &
				 TwoPointExpValNC_rc,&
				 TwoPointExpValNC_cr
END INTERFACE  TwoPointExpValNC
	
CONTAINS

SUBROUTINE FormSingleSiteRho(rho1, Lambda0, Gamma1, Lambda1)
!
!Purpose: Calculate the single site density matrix on the site where Gamma resides, store in rho1
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: rho1(:,:)
REAL(KIND=rKind) :: Lambda0(:), Lambda1(:)
COMPLEX(KIND=rKind) :: Gamma1(:,:,:)
INTEGER :: chi0, chi1, i, j, alpha, beta
	chi0 = SIZE(Gamma1,1)
	chi1 = SIZE(Gamma1,3)
	DO i=1,localSize
		DO j=1,localSize
			rho1(i,j)=CMPLX(0.0,KIND=rKind)
			DO alpha=1,chi0
				DO beta=1,chi1
					rho1(i,j) = rho1(i,j)+Lambda0(alpha)*Gamma1(alpha,i,beta)*Lambda1(beta) &
								*Lambda1(beta)*CONJG(Gamma1(alpha,j,beta))*Lambda0(alpha)
				END DO
			END DO
		END DO
	END DO
END SUBROUTINE FormSingleSiteRho
		
SUBROUTINE SingleSiteDensityMatrix(rho,Gammas,Lambdas)
!
!Purpose: Calculate the single site density matrix on every site, store in rho
!
IMPLICIT NONE
TYPE(matrix), INTENT(OUT) :: rho(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: i
	DO i=1,systemSize
		CALL FormSingleSiteRho(rho(i)%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
	END DO
END SUBROUTINE SingleSiteDensityMatrix

!!!!!!!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE OneSiteExpVal!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneSiteExpVal_mr(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the real operator Op on every site
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: expList(systemSize)
TYPE(matrixReal), INTENT(IN) :: Op
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_r'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		expList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			expList(i)=REAL(TraceMatmul(Op%mr,rho%m),KIND=rKind)
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF
END SUBROUTINE OneSiteExpVal_mr

SUBROUTINE OneSiteExpVal_m(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the complex operator Op on every site
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: expList(systemSize)
TYPE(matrix), INTENT(IN) :: Op
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_c'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		expList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			expList(i)=TraceMatmul(Op%m,rho%m)
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_c'
			END IF
END SUBROUTINE OneSiteExpVal_m

SUBROUTINE OneSiteExpVal_r(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the real operator Op on every site
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: expList(systemSize)
REAL(KIND=rKind), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_r'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		expList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			expList(i)=REAL(TraceMatmul(Op,rho%m),KIND=rKind)
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF
END SUBROUTINE OneSiteExpVal_r

SUBROUTINE OneSiteExpVal_c(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the complex operator Op on every site
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: expList(systemSize)
COMPLEX(KIND=rKIND), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_c'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		expList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			expList(i)=TraceMatmul(Op,rho%m)
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_c'
			END IF
END SUBROUTINE OneSiteExpVal_c
!!!!!!!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE OneSiteExpVal!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE OneSiteVar!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneSiteVar_mr(varList,Op, Gammas, Lambdas)
!
!Purpose: Calculate variance of the real operator Op on every site
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: varList(systemSize)
TYPE(matrixReal), INTENT(IN) :: Op
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_r'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		varList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			varList(i)=REAL(TraceMatmul(MATMUL(Op%mr,Op%mr),rho%m),KIND=rKind)-REAL(TraceMatmul(Op%mr,rho%m),KIND=rKind)**2
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF
END SUBROUTINE OneSiteVar_mr

SUBROUTINE OneSiteVar_m(varList,Op, Gammas, Lambdas)
!
!Purpose: Calculate variance of the complex operator Op on every site
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: varList(systemSize)
TYPE(matrix), INTENT(IN) :: Op
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_c'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		varList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			varList(i)=TraceMatmul(MATMUL(Op%m,Op%m),rho%m)-(TraceMatmul(Op%m,rho%m))**2
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_c'
			END IF
END SUBROUTINE OneSiteVar_m

SUBROUTINE OneSiteVar_r(varList,Op, Gammas, Lambdas)
!
!Purpose: Calculate variance of the real operator Op on every site
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: varList(systemSize)
REAL(KIND=rKind), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_r'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		varList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			varList(i)=REAL(TraceMatmul(MATMUL(Op,Op),rho%m),KIND=rKind)-REAL(TraceMatmul(Op,rho%m),KIND=rKind)**2
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF
END SUBROUTINE OneSiteVar_r

SUBROUTINE OneSiteVar_c(varList,Op, Gammas, Lambdas)
!
!Purpose: Calculate variance of the complex operator Op on every site
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: varList(systemSize)
COMPLEX(KIND=rKind), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpVal_c'
			END IF
		
	DO i=1,systemSize,1
		rho%m=0.0_rKind
		varList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			varList(i)=TraceMatmul(MATMUL(Op,Op),rho%m)-(TraceMatmul(Op,rho%m))**2
	END DO
				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_c'
			END IF
END SUBROUTINE OneSiteVar_c
!!!!!!!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE OneSiteVar!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GKernal(gee,Gamma,GammaP,Lambda)
!
!Purpose: First step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: gee(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:),GammaP(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda(:)
INTEGER :: alpha, beta, betapr,i

gee=0.0_rKind
DO alpha=1,SIZE(Gamma,3)
	DO beta=1,SIZE(Gamma,3)
		DO betapr=1,SIZE(Gamma,3)
			DO i=1,localSize,1
			gee(alpha,beta)=gee(alpha,beta)+(Lambda(betapr)**2)*CONJG(GammaP(beta,i,betapr))*Gamma(alpha,i,betapr)
			END DO
		END DO
	END DO
END DO

END SUBROUTINE GKernal

SUBROUTINE GNext(gee,Gamma,GammaP,Lambda)
!
!Purpose:  Recursive step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: gee(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:),GammaP(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: geeTemp(:,:)
INTEGER :: alpha, beta, alphapr,betapr,i

ALLOCATE(geetemp(SIZE(Gamma,3),SIZE(Gamma,3)))

geetemp=0.0_rKind
DO alpha=1,SIZE(Gamma,3)
	DO beta=1,SIZE(Gamma,3)
		DO betapr=1,SIZE(Gamma,3)
			DO alphapr=1,SIZE(Gamma,3)
				DO i=1,localSize,1
			geetemp(alpha,beta)=geetemp(alpha,beta)+Lambda(alphapr)*Lambda(betapr)&
			*CONJG(GammaP(beta,i,betapr))*Gamma(alpha,i,alphapr)*gee(alphapr,betapr)
				END DO
			END DO
		END DO
	END DO
END DO

gee=geetemp

DEALLOCATE(geetemp)

END SUBROUTINE GNext

SUBROUTINE GContraction(obsv,gee,Gamma,GammaP,Lambda1,Lambda2)
!
!Purpose: Final step in calculating the two-site observable
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: obsv
COMPLEX(KIND=rKind), INTENT(IN) :: gee(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:),GammaP(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda1(:),Lambda2(:)
INTEGER :: alpha, beta, betapr,i

obsv=0.0_rKind
DO alpha=1,SIZE(Gamma,3)
	DO beta=1,SIZE(Gamma,3)
		DO betapr=1,SIZE(Gamma,3)
			DO i=1,localSize,1
			obsv=obsv+(Lambda1(betapr)**2)*CONJG(GammaP(betapr,i,beta))*Gamma(betapr,i,alpha)*&
			Lambda2(alpha)*Lambda2(beta)*gee(alpha,beta)
			END DO
		END DO
	END DO
END DO

END SUBROUTINE GContraction

!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE TwoSiteExpValG !!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TwoSiteExpValG_r(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the real two-site operator Op1XOp2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)		
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2, i,alpha, beta, betapr
INTEGER, INTENT(IN), OPTIONAL :: phaseStat

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
            observable(l1,l1)=TraceMatmul(MATMUL(Op1,Op2),rho1)
		END DO				
		
		DO l2=systemSize,2,-1		
		GammaP = Gammas(l2)%t
!Compute the initial G 
CALL OneSiteOp(Op2,GammaP)
CALL GKernal(gee,Gammas(l2)%t,GammaP,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(l1,l2),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
				observable(l2,l1)=CONJG(observable(l1,l2))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO
		END DO
DEALLOCATE(gee)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpValG_r

SUBROUTINE TwoSiteExpValG_c(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the complex two-site operator Op1XOp2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2, i,alpha, beta, betapr
INTEGER, INTENT(IN), OPTIONAL :: phaseStat

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(MATMUL(Op1,Op2),rho1)
		END DO				
		
		DO l2=systemSize,2,-1		
		GammaP = Gammas(l2)%t
!Compute the initial G 
CALL OneSiteOp(Op2,GammaP)
CALL GKernal(gee,Gammas(l2)%t,GammaP,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(l1,l2),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
				observable(l2,l1)=CONJG(observable(l1,l2))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO
		END DO
DEALLOCATE(gee)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpValG_c

SUBROUTINE TwoSiteExpValG_rc(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the real/complex two-site operator Op1XOp2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)		
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2, i,alpha, beta, betapr
INTEGER, INTENT(IN), OPTIONAL :: phaseStat

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(MATMUL(Op1,Op2),rho1)
		END DO				
		
		DO l2=systemSize,2,-1		
		GammaP = Gammas(l2)%t
!Compute the initial G 
CALL OneSiteOp(Op2,GammaP)
CALL GKernal(gee,Gammas(l2)%t,GammaP,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(l1,l2),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
				observable(l2,l1)=CONJG(observable(l1,l2))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO
		END DO
DEALLOCATE(gee)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpValG_rc

SUBROUTINE TwoSiteExpValG_cr(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the complex/real two-site operator Op1XOp2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2, i,alpha, beta, betapr
INTEGER, INTENT(IN), OPTIONAL :: phaseStat

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(MATMUL(Op1,Op2),rho1)
		END DO				
		
		DO l2=systemSize,2,-1		
		GammaP = Gammas(l2)%t
!Compute the initial G 
CALL OneSiteOp(Op2,GammaP)
CALL GKernal(gee,Gammas(l2)%t,GammaP,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(l1,l2),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
				observable(l2,l1)=CONJG(observable(l1,l2))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO
		END DO
DEALLOCATE(gee)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpValG_cr

!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE TwoSiteExpValG !!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ThetaKernal(Theta,Lambda1,Gamma,Lambda2)
!
!Purpose: Initial step in calculating the two-site density matrix
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda1(:), Lambda2(:)
INTEGER :: alpha, beta, eta, i, j

	DO alpha=1,SIZE(Gamma,1)
		DO beta=1,SIZE(Gamma,1)
			DO i=1,localSize
				DO j=1,localSize
					Theta(alpha,i,j,beta)=CMPLX(0.0,KIND=rKind)
					DO eta=1,SIZE(Gamma,3)
						Theta(alpha,i,j,beta)=Theta(alpha,i,j,beta) &
											  +Lambda1(alpha)*Gamma(alpha,i,eta)*Lambda2(eta) &
											  *Lambda2(eta)*CONJG(Gamma(beta,j,eta))*Lambda1(beta)
					END DO
				END DO
			END DO
		END DO
	END DO
END SUBROUTINE ThetaKernal
	
SUBROUTINE ThetaNext(Theta, Gamma, GammaP, Lambda)
!
!Purpose: Recursive step in calculating the two-site density matrix
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(INOUT) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: Gamma(:,:,:), GammaP(:,:,:)
REAL(KIND=rKind) :: Lambda(:)
COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),localSize,localSize,localSize,SIZE(Gamma,1))
INTEGER :: alpha, beta, eta, i, j, k, l, d, chi
	chi=SIZE(Gamma,3)
	d=localSize
	DO alpha=1,SIZE(Theta,1)
		DO beta=1,SIZE(Gamma,1)
			DO i=1,d
				DO j=1,d
					DO k=1,d
						ThetaTemp(alpha,i,j,k,beta)=CMPLX(0.0,KIND=rKind)
						DO eta=1,chi
							ThetaTemp(alpha,i,j,k,beta) = ThetaTemp(alpha,i,j,k,beta) &
														+ Theta(alpha,i,j,eta)*CONJG(Gamma(beta,k,eta)) &
														* Lambda(beta)
						END DO
					END DO
				END DO
			END DO
		END DO
	END DO
	
	DO alpha=1,SIZE(GammaP,1)
		DO beta=1,SIZE(Gamma,1)
			DO i=1,d
				DO j=1,d
					Theta(alpha,i,j,beta)=CMPLX(0.0,KIND=rKind)
					DO k=1,d
						DO eta=1,chi
							Theta(alpha,i,j,beta) = Theta(alpha,i,j,beta) &
													+ Lambda(alpha)*GammaP(alpha,k,eta)*ThetaTemp(eta,i,j,k,beta)
						END DO
					END DO
				END DO
			END DO
		END DO
	END DO
END SUBROUTINE ThetaNext

SUBROUTINE TwoSiteRho(rho2, Theta, Gamma, GammaP, Lambda)
!
!Purpose: Calculate the two-site density matrix
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: rho2(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: Gamma(:,:,:), GammaP(:,:,:)
REAL(KIND=rKind) :: Lambda(:)
!COMPLEX(KIND=rKind) :: ThetaTemp(SIZE(Theta,1),localSize,localSize,localSize,SIZE(Gamma,1))
COMPLEX(KIND=rKind), ALLOCATABLE :: ThetaTemp(:,:,:,:,:)
INTEGER :: alpha, beta, eta, i, j, k, l, d, l1
ALLOCATE(ThetaTemp(SIZE(Theta,1),localSize,localSize,localSize,SIZE(Gamma,1)))

	d=localSize
	DO alpha=1,SIZE(Theta,1)
		DO beta=1,SIZE(Gamma,1)
			DO i=1,d
				DO j=1,d
					DO k=1,d
						ThetaTemp(alpha,i,j,k,beta)=CMPLX(0.0,KIND=rKind)
						DO eta=1,SIZE(Gamma,3)
							ThetaTemp(alpha,i,j,k,beta) = ThetaTemp(alpha,i,j,k,beta) &
														+ Theta(alpha,i,j,eta)*CONJG(Gamma(beta,k,eta)) &
														* Lambda(beta)
						END DO
					END DO
				END DO
			END DO
		END DO
	END DO
			
	DO i=1,d
		DO j=1,d
			DO k=1,d
				DO l=1,d
					rho2((i-1)*d+j,(l-1)*d+k)=CMPLX(0.0,KIND=rKind)
					DO alpha=1,SIZE(Gamma,1)
						DO beta=1,SIZE(Gamma,3)
							rho2((i-1)*d+j,(l-1)*d+k) = rho2((i-1)*d+j,(l-1)*d+k) &
									+ Lambda(alpha)*GammaP(alpha,i,beta)*ThetaTemp(beta,j,k,l,alpha)
						END DO
					END DO
				END DO
			END DO
		END DO
	END DO

DEALLOCATE(ThetaTemp)

END SUBROUTINE TwoSiteRho

!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE TwoSiteExpVal !!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TwoSiteExpVal_r(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the real two-site operator Op2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)		
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: Theta(:,:,:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
INTEGER, INTENT(IN), OPTIONAL :: phaseStat

ALLOCATE(Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(Op1,rho1)
		END DO
		DO l2=systemSize,2,(-1)
		!Compute the initial theta (Theta at site l2)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
!Fermi Phase added Here
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
				!Form the two-site density matrix from the cumulative Theta
				CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
				!Calculate the two-site density matrix
				observable(l1,l2)=TraceMatmul(Op2,rho2)
				observable(l2,l1)=CONJG(observable(l1,l2))
				!Add the site=l1 local tensors to the cumulative Theta
				CALL ThetaNext(Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
			END DO
		END DO
DEALLOCATE(Theta)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpVal_r

SUBROUTINE TwoSiteExpVal_c(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the complex two-site operator Op2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: Theta(:,:,:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
INTEGER, INTENT(IN), OPTIONAL :: phaseStat
ALLOCATE(Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(Op1,rho1)
		END DO
		DO l2=systemSize,2,(-1)
		!Compute the initial theta (Theta at site l2)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
!Fermi Phase added Here
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
				!Form the two-site density matrix from the cumulative Theta
				CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
				!Calculate the two-site density matrix
				observable(l1,l2)=TraceMatmul(Op2,rho2)
				observable(l2,l1)=CONJG(observable(l1,l2))
				!Add the site=l1 local tensors to the cumulative Theta
				CALL ThetaNext(Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
			END DO
		END DO
DEALLOCATE(Theta)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpVal_c

SUBROUTINE TwoSiteExpVal_rc(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the real two-site operator Op2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)		
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: Theta(:,:,:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
INTEGER, INTENT(IN), OPTIONAL :: phaseStat
ALLOCATE(Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(Op1,rho1)
		END DO
		DO l2=systemSize,2,(-1)
		!Compute the initial theta (Theta at site l2)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
!Fermi Phase added Here
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
				!Form the two-site density matrix from the cumulative Theta
				CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
				!Calculate the two-site density matrix
				observable(l1,l2)=TraceMatmul(Op2,rho2)
				observable(l2,l1)=CONJG(observable(l1,l2))
				!Add the site=l1 local tensors to the cumulative Theta
				CALL ThetaNext(Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
			END DO
		END DO
DEALLOCATE(Theta)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpVal_rc

SUBROUTINE TwoSiteExpVal_cr(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the complex two-site operator Op2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: Theta(:,:,:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
INTEGER, INTENT(IN), OPTIONAL :: phaseStat
ALLOCATE(Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

		DO l1=1,systemSize
		!On-diagonal elements use Op1
			CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
			observable(l1,l1)=TraceMatmul(Op1,rho1)
		END DO
		DO l2=systemSize,2,(-1)
		!Compute the initial theta (Theta at site l2)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				GammaP = Gammas(l1)%t
!Fermi Phase added Here
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
				!Form the two-site density matrix from the cumulative Theta
				CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
				!Calculate the two-site density matrix
				observable(l1,l2)=TraceMatmul(Op2,rho2)
				observable(l2,l1)=CONJG(observable(l1,l2))
				!Add the site=l1 local tensors to the cumulative Theta
				CALL ThetaNext(Theta,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
			END DO
		END DO
DEALLOCATE(Theta)
DEALLOCATE(GammaP)
END SUBROUTINE TwoSiteExpVal_cr

!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE TwoSiteExpVal !!!!!!!!!!!!!!!!!!!!!!!!

COMPLEX(KIND=rKind) FUNCTION InnerProduct(GammasL, LambdasL, GammasR, LambdasR)
!
!Purpose: Calculate the inner product of the wavefunction by the Gammas and Lambdas.
!This routine can also be used to compute the Fidelity or Loschmidt echo if GammasL/LambdasL
!is some initial state and GammasR/LambdasR some time evolved final state
!
IMPLICIT NONE
TYPE(tensor), POINTER :: GammasL(:), GammasR(:)
TYPE(vector), POINTER :: LambdasL(:), LambdasR(:)
COMPLEX(KIND=rKind) :: temp
COMPLEX(KIND=rKind) :: normKernal(SIZE(GammasL(1)%t,3),SIZE(GammasL(1)%t,3))
COMPLEX(KIND=rKind) :: GammaTemp(SIZE(GammasL(1)%t,3),localSize,SIZE(GammasL(1)%t,3))
INTEGER :: alpha,beta,gamma,eta,i,n,chi
		chi=SIZE(GammasL(1)%t,3)
		DO alpha=1,chi
			DO beta=1,chi
				normKernal(alpha,beta)=CMPLX(0.0,KIND=rKind)
				DO i=1,localSize
					normKernal(alpha,beta)=normKernal(alpha,beta)+LambdasL(2)%v(alpha)*CONJG(GammasL(1)%t(1,i,alpha)) &
											*GammasR(1)%t(1,i,beta)*LambdasR(2)%v(beta)
				END DO
			END DO
		END DO
		DO n=2,(systemSize-1)
			DO gamma=1,chi
				DO i=1,localSize
					DO beta=1,chi
						GammaTemp(gamma,i,beta)=CMPLX(0.0,KIND=rKind)
						DO eta=1,chi
							GammaTemp(gamma,i,beta) = GammaTemp(gamma,i,beta) &
												    + normKernal(gamma,eta)*GammasR(n)%t(eta,i,beta)
						END DO
					END DO
				END DO
			END DO
			DO alpha=1,chi
				DO beta=1,chi
					normKernal(alpha,beta)=CMPLX(0.0,KIND=rKind)
					DO gamma=1,chi
						DO i=1,localSize
							normKernal(alpha,beta) = normKernal(alpha,beta) &
												   + LambdasL(n+1)%v(alpha)*CONJG(GammasL(n)%t(gamma,i,alpha)) &
												   * GammaTemp(gamma,i,beta)*LambdasR(n+1)%v(beta)
						END DO
					END DO
				END DO
			END DO
		END DO
		temp=CMPLX(0.0,KIND=rKind);
		DO alpha=1,chi
			DO beta=1,chi
				DO i=1,localSize
					temp = temp+CONJG(GammasL(systemSize)%t(alpha,i,1))*normKernal(alpha,beta) &
							*GammasR(systemSize)%t(beta,i,1)
				END DO
			END DO
		END DO
		InnerProduct=temp
END FUNCTION InnerProduct

SUBROUTINE OnSiteNumber(number, Gammas, Lambdas,siteIndex, comPonent)
!
!Purpose: Calculate the total number on site
!OPTIONAL argument comPonent specifies the presence of internal degrees of freedom
!Calculates the average number on site in the component comPonent
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: number
INTEGER, INTENT(IN) :: siteIndex
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER, INTENT(IN), OPTIONAL :: comPonent 
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OnSiteNumber'
			END IF
		number=0.0_rKind
		
!Internal degrees of freedom present
IF(PRESENT(comPonent)) THEN
		
			CALL FormSingleSiteRho(rho%m, Lambdas(siteIndex)%v, Gammas(siteIndex)%t, Lambdas(siteIndex+1)%v)	
			dumNum=TraceMatmul(MATMUL(Transpose(a_opS(comPonent)%mr),a_opS(comPonent)%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
!Internal degreesof freedom absent
ELSE
			CALL FormSingleSiteRho(rho%m, Lambdas(siteIndex)%v, Gammas(siteIndex)%t, Lambdas(siteIndex+1)%v)	
			dumNum=TraceMatmul(MATMUL(Transpose(a_op%mr),a_op%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
END IF

				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OnSiteNumber'
			END IF
END SUBROUTINE OnSiteNumber

SUBROUTINE TotalNumber(number, Gammas, Lambdas, comPonent)
!
!Purpose: Calculate the total number for spinless code
!OPTIONAL argument comPonent specifies the presence of internal degrees of freedom
!Calculates the average number in the component comPonent
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: number
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER, INTENT(IN), OPTIONAL :: comPonent
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i
		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in TotalNumber'
			END IF
		number=0.0_rKind

!Internal degrees of freedom present
IF(PRESENT(comPonent)) THEN
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			dumNum=TraceMatmul(MATMUL(Transpose(a_opS(comPonent)%mr),a_opS(comPonent)%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
		END DO
		number = number/REAL(systemSize, KIND=rKind)
!Internal degrees of freedom absent
ELSE
		
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			dumNum=TraceMatmul(MATMUL(Transpose(a_op%mr),a_op%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
		END DO
		number = number/REAL(systemSize, KIND=rKind)
END IF		
		
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in TotalNumber'
			END IF
END SUBROUTINE TotalNumber

SUBROUTINE TotalOneSite_r(total,Op, Gammas, Lambdas)
!
!Purpose: Calculate the total number value of Op across all lattice sites: real version
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: total
REAL(KIND=rKind) ,INTENT(INOUT) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i
		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in TotalOneSite'
			END IF
		total=0.0_rKind
		
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			dumNum=TraceMatmul(Op,rho%m)
				total=total+REAL(dumNum, KIND=rKind)
		END DO
		
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in TotalNumber'
			END IF
END SUBROUTINE TotalOneSite_r

SUBROUTINE TotalOneSite_c(total,Op, Gammas, Lambdas)
!
!Purpose: Calculate the total number value of Op across all lattice sites: complex version
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: total
COMPLEX(KIND=rKind) ,INTENT(INOUT) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i
		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in TotalOneSite'
			END IF
		total=0.0_rKind
		
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			dumNum=TraceMatmul(Op,rho%m)
				total=total+dumNum
		END DO
		
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in TotalNumber'
			END IF
END SUBROUTINE TotalOneSite_c

SUBROUTINE TotalOneSite_mr(total,Op, Gammas, Lambdas)
!
!Purpose: Calculate the total number value of Op across all lattice sites: matrixreal version
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: total
TYPE(matrixReal) :: Op
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i
		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in TotalOneSite'
			END IF
		total=0.0_rKind
		
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			dumNum=TraceMatmul(Op%mr,rho%m)
				total=total+REAL(dumNum, KIND=rKind)
		END DO
		
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in TotalNumber'
			END IF
END SUBROUTINE TotalOneSite_mr

SUBROUTINE TotalOneSite_m(total,Op, Gammas, Lambdas)
!
!Purpose: Calculate the total number value of Op across all lattice sites: matrix version
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: total
TYPE(matrix) :: Op
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i
		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in TotalOneSite'
			END IF
		total=0.0_rKind
		
		DO i=1,systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			dumNum=TraceMatmul(Op%m,rho%m)
				total=total+dumNum
		END DO
		
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in TotalNumber'
			END IF
END SUBROUTINE TotalOneSite_m

SUBROUTINE LocalNumDev(numbers, deviations, Gammas, Lambdas, comPonent)
!
!Purpose: Calculate the local number deviations
!OPTIONAL argument comPonent specifies the presence of internal degrees of freedom
!Calculates the Local number deviations in the component comPonent
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: numbers(:), deviations(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER, OPTIONAL, INTENT(IN) :: comPonent
COMPLEX(KIND=rKind) :: dumNum
INTEGER :: i

ALLOCATE(rho%m(localSize,localSize))
!Internal degrees of freedom present
IF(PRESENT(comPonent)) THEN
	DO i=1,systemSize
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
		!Compute the expectation of number on each site
		dumNum=TraceMatmul(MATMUL(Transpose(a_opS(comPonent)%mr),a_opS(comPonent)%mr),rho%m)
		numbers(i) = REAL(dumNum, KIND=rKind)
		!Compute the expectation of number squared on each site		
		dumNum=TraceMatmul(MATMUL(MATMUL(TRANSPOSE(a_opS(comPonent)%mr),a_opS(comPonent)%mr),&
			MATMUL(TRANSPOSE(a_opS(comPonent)%mr),a_opS(comPonent)%mr)),rho%m)
!Compute the standard deviation
		deviations(i) = REAL(dumNum, KIND=rKind) - numbers(i)**2
		deviations(i) = SQRT(deviations(i))
	END DO
ELSE
	DO i=1,systemSize
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
		!Compute the expectation of number on each site
		dumNum=TraceMatmul(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),rho%m)
		numbers(i) = REAL(dumNum, KIND=rKind)
		!Compute the expectation of number squared on each site		
        dumNum=TraceMatmul(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),&
			MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),rho%m)
!Compute the standard deviation
		deviations(i) = REAL(dumNum, KIND=rKind) - numbers(i)**2
		deviations(i) = SQRT(deviations(i))
	END DO
END IF

DEALLOCATE(rho%m)
END SUBROUTINE LocalNumDev

SUBROUTINE LocalEnergy(energy, H, Gammas, Lambdas)
!
!Purpose: Calculate the energy associated with each lattice bond.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy(:)
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(1)%t, 3), localSize, localSize, SIZE(Gammas(1)%t, 3))
COMPLEX(KIND=rKind) :: rho2(localSize*localSize, localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: i, l1, l2	
	energy = 0.0_rKind
	DO l2 = 2, systemSize
		CALL ThetaKernal(Theta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v)
		l1 = l2 - 1
		CALL TwoSiteRho(rho2, Theta, Gammas(l1)%t, Gammas(l1)%t, Lambdas(l1)%v)
        dumEn=TraceMatmul(H(l1)%m, rho2)
		energy(l1) = REAL(dumEn, KIND=rKind)
	END DO

	IF(BoundaryCond=='P') THEN
		l2=systemSize
		!Compute the initial theta (Theta at site l2)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
		IF(l1==1) THEN
				!Form the two-site density matrix from the cumulative Theta
				CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t, Gammas(l1)%t,Lambdas(l1)%v)
				!Calculate the two-site density matrix
				dumEn=TraceMatmul(H(systemSize)%m,rho2)
				energy(systemSize)= REAL(dumEn, KIND=rKind)
		ELSE
				!Add the site=l1 local tensors to the cumulative Theta
				CALL ThetaNext(Theta,Gammas(l1)%t, Gammas(l1)%t,Lambdas(l1)%v)
		END IF
			END DO
	END IF


END SUBROUTINE LocalEnergy


SUBROUTINE LocalSpin(spinVec, Gammas, Lambdas)
!
!Purpose: Calculate the total spin (S_i+S_{i+1})^2 associated with each lattice bond.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: spinVec(:)
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(1)%t, 3), localSize, localSize, SIZE(Gammas(1)%t, 3))
COMPLEX(KIND=rKind) :: rho2(localSize*localSize, localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: i, l1, l2	
	spinVec = 0.0_rKind
	DO l2 = 2, systemSize
		CALL ThetaKernal(Theta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v)
		l1 = l2 - 1
		CALL TwoSiteRho(rho2, Theta, Gammas(l1)%t, Gammas(l1)%t, Lambdas(l1)%v)
	spinVec(l1)=2.0_rKind*(REAL(TraceMatmul(TensorProd(Sx_opS%mr,Sx_opS%mr), rho2), KIND=rKind) +&
			REAL(TraceMatmul(TensorProd(Sy_opS%m,Sy_opS%m), rho2), KIND=rKind) +&
			REAL(TraceMatmul(TensorProd(Sz_opS%mr,Sz_opS%mr), rho2), KIND=rKind))+&
			REAL(TraceMatmul(TensorProd(Ssq_opS%mr,one_op%mr),rho2),KIND=rKind)+&
			REAL(TraceMatmul(TensorProd(one_op%mr,Ssq_opS%mr),rho2),KIND=rKind)
	END DO

IF(BoundaryCond=='P') THEN
	l2=systemSize
	CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
	DO l1=l2-1, 2, (-1)
		CALL ThetaNext(Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
	END DO
			CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
spinVec(systemSize)=2.0_rKind*(REAL(TraceMatmul(TensorProd(Sx_opS%mr,Sx_opS%mr), rho2), KIND=rKind) +&
			REAL(TraceMatmul(TensorProd(Sy_opS%m,Sy_opS%m), rho2), KIND=rKind) +&
			REAL(TraceMatmul(TensorProd(Sz_opS%mr,Sz_opS%mr), rho2), KIND=rKind))+&
			REAL(TraceMatmul(TensorProd(Ssq_opS%mr,one_op%mr),rho2),KIND=rKind)+&
			REAL(TraceMatmul(TensorProd(one_op%mr,Ssq_opS%mr),rho2),KIND=rKind)
END IF


END SUBROUTINE LocalSpin
	
SUBROUTINE TotalEnergy(energy, H, Gammas, Lambdas)
!
!Purpose: Calculate the energy eigenvalue associated with the Hamiltonian H
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind) :: Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: l1, l2	
	energy=0.0_rKind
	DO l2=2,systemSize
		CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
		l1=l2-1
		CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
		dumEn=TraceMatmul(H(l1)%m,rho2)
		energy=energy+REAL(dumEn, KIND=rKind)
	END DO

IF(BoundaryCond=='P') THEN
	l2=systemSize
	CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
	DO l1=l2-1, 2, (-1)
		CALL ThetaNext(Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
	END DO
			CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
			dumEn=TraceMatmul(H(systemSize)%m,rho2)
			energy=energy+REAL(dumEn, KIND=rKind)
END IF


END SUBROUTINE TotalEnergy
	
SUBROUTINE Qdepletion(depletion, rho, population, CW)
!
!Purpose: Calculate the Quantum depletion
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT), OPTIONAL :: CW(:)
COMPLEX(KIND=rKind), INTENT(IN) :: rho(:,:)
REAL(KIND=rKind), INTENT(IN) :: population
REAL(KIND=rKind), INTENT(OUT) :: depletion
REAL(KIND=rKind), ALLOCATABLE :: S(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: U(:,:), VT(:,:), work(:,:), rwork(:,:), rhoCopy(:,:)
INTEGER :: workSize, info, l1, l2
CHARACTER(1) :: jobu, jobvt

IF(population==0.0) THEN
depletion=0.0_rKind
	IF(PRESENT(CW)) THEN
		DO l1=1,systemSize
			CW(l1)=0.0_rKind
		END DO
	END IF
ELSE

!Allocate the variables needed to perform an SVD on rho		
	jobu='A'
	jobvt='A'
	workSize=5*systemSize
	ALLOCATE(S(systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate S in Qdepletion'
			END IF
	ALLOCATE(U(systemSize,systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate U in Qdepletion'
			END IF
	ALLOCATE(VT(systemSize,systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate VT in Qdepletion'
			END IF
	ALLOCATE(rhoCopy(systemSize,systemSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rhoCopy in Qdepletion'
			END IF
	ALLOCATE(work(workSize,workSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate work in Qdepletion'
			END IF
	ALLOCATE(rwork(workSize,workSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rwork in Qdepletion'
			END IF
	DO l1=1,systemSize
		DO l2=1,systemSize
			rhoCopy(l1,l2)=rho(l1,l2)
		END DO
	END DO
	!Perform an SVD on rho
	CALL ZGESVD(jobu, jobvt, systemSize, systemSize, rhoCopy, systemSize, S, U, systemSize, VT, systemSize, &
				work, workSize, rwork, info)


	depletion=1.0_rKind-S(1)/population

IF(PRESENT(CW)) THEN
	DO l1=1,systemSize
		CW(l1)=SQRT(S(1))*VT(1,l1)
	END DO
END IF

	!Deallocate the unnecessary variables
	DEALLOCATE(U,VT,work,rwork,S, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate SVD variables in Qdepletion'
			END IF
	DEALLOCATE(rhoCopy, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rhocopy in Qdepletion'
			END IF
END IF
END SUBROUTINE Qdepletion

REAL(KIND=rKind) FUNCTION MeyerQmeasure(Gammas, Lambdas)
!
!Purpose: Calculate the Meyer Q-measure (average local impurity)
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumPur
INTEGER :: i
REAL(KIND=rKind) :: totalPurity

	ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in MeyerQMeasure'
			END IF
	totalPurity = 0.0_rKind
	! Calculate total purity, i.e., sum of tr(rho**2).
	DO i = 1, systemSize
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
		dumPur=TraceMatmul(rho%m, rho%m)
		totalPurity = totalPurity + dumPur
	END DO
MeyerQmeasure = localSize*1.0_rKind/(localSize*1.0_rKind-1.0_rKind) * (1.0_rKind - totalPurity/systemSize*1.0_rKind) ! Calculate average impurity.
	DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in MeyerQMeasure'
			END IF
END FUNCTION MeyerQmeasure

REAL(KIND=rKind) FUNCTION ChainEntropy(link,Lambdas)
!
!Purpose: Calculate the entropy of entanglement of the chain to the left of link
!with the chain to the right on link (in the MPS approximation)
!
INTEGER, INTENT(IN) :: link
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind) :: temp
INTEGER :: chi, i

chi=SIZE(Lambdas(link)%v)
temp = 0.0_rKind
DO i=1,chi
IF(Lambdas(link)%v(i).ne.0.0_rKind) THEN
temp=temp-Lambdas(link)%v(i)*Lambdas(link)%v(i)*LOG(ABS(Lambdas(link)%v(i)*Lambdas(link)%v(i)))/LOG(1.0_rKind*localSize)
END IF
END DO
ChainEntropy=temp

END FUNCTION ChainEntropy



SUBROUTINE LocalEntropyDist(entDist, Gammas, Lambdas,tsalliSq)
!
!Purpose: Calculate the local entropy at each lattice site.
!If the OPTIONAL argument tsalliSq is present, compute the Tsallis entropy
!S_q=(1/(1-q))(Tr( rho^q) -1), else compute the von Neumann entropy
!S=-Tr( rho Log_d rho)
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: entDist(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: tsallisQ
INTEGER :: i, j, workSize, info
REAL(KIND=rKind) :: evals(localSize), temp
COMPLEX(KIND=rKind), ALLOCATABLE :: workArr(:), rworkArr(:)


!Allocate single-site density matrix
	ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in LocalEntropyDist'
			END IF
!Allocate workspace for ZHEEV
	workSize = 2*localSize-1
	ALLOCATE(workarr(workSize), rworkarr(3*localSize-2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
			END IF

	DO i = 1, systemSize
		temp = 0.0_rKind
!Form single-site density matrix
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
!Diagonalize density matrix
		CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info) ! Diagonalize single site density matrix.

IF(PRESENT(tsalliSq)) THEN
!Compute the Tsallis entropy
temp=-(1.0_rKind/(1.0_rKind-tsallisQ))
DO j=1,localSize
	IF(evals(j).ne.0.0_rKind) THEN
	temp=temp+(1.0_rKind/(1.0_rKind-tsalliSq))*(ABS(evals(j))**tsalliSq)
	ELSE
	temp=temp+0.0_rKind
	END IF
END DO
entDist(i) = temp
ELSE
!Compute vonNeumann entropy
			DO j = 1, localSize
	IF(evals(j).ne.0.0_rKind) THEN
        temp = temp + evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
	ELSE
	temp=temp+0.0_rKind
	END IF

			END DO
		entDist(i) = -1.0_rKind*temp ! Set value of entropy on-site.
END IF

	END DO

	DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in LocalEntropyDist'
			END IF

	DEALLOCATE(workarr, rworkarr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
			END IF
END SUBROUTINE LocalEntropyDist

SUBROUTINE TwoBodyEntropyDist(entDist, Gammas, Lambdas,tsalliSq)
!
!Purpose: Calculate the 2-site entropy for each site pair.
!If the OPTIONAL argument tsalliSq is present, compute the Tsallis entropy
!S_q=(1/(1-q))(Tr( rho^q) -1), else compute the von Neumann entropy
!S=-Tr( rho Log_d rho)
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: entDist(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: tsallisQ
TYPE(matrix) :: rho2
COMPLEX(KIND=rKind),ALLOCATABLE :: Theta(:,:,:,:)
INTEGER :: i, j, workSize, info, l1, l2
REAL(KIND=rKind) :: evals(localSize*localSize), temp
COMPLEX(KIND=rKind), ALLOCATABLE :: workArr(:), rworkArr(:)
ALLOCATE(Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))

!Allocate single-site density matrix
	ALLOCATE(rho2%m(localSize*localSize, localSize*localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in LocalEntropyDist'
			END IF
			
!Allocate workspace for ZHEEV
	workSize = 2*localSize*localSize-1
	ALLOCATE(workarr(workSize), rworkarr(3*localSize*localSize-2), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
			END IF

		DO l2=systemSize,2,(-1)
		!Compute the initial theta (Theta at site l2)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
			DO l1=(l2-1),1,(-1)
				entDist(l1,l2)=0.0_rKind
				temp = 0.0_rKind
!Fermi Phase added Here if needed
				!Form the two-site density matrix from the cumulative Theta
				CALL TwoSiteRho(rho2%m,Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
				!Diagonalize the two-site density matrix
				CALL ZHEEV('N', 'U', localSize*localSize, rho2%m, localSize*localSize, evals, workArr, workSize, rworkArr, info)

IF(PRESENT(tsalliSq)) THEN
!Compute the Tsallis entropy
temp=-(1.0_rKind/(1.0_rKind-tsallisQ))
DO j=1,localSize*localSize
	IF(evals(j).ne.0.0_rKind) THEN
	temp=temp+(1.0_rKind/(1.0_rKind-tsalliSq))*(ABS(evals(j))**tsalliSq)
	ELSE
	temp=temp+0.0_rKind
	END IF
END DO
entDist(l1,l2) = temp
ELSE
!Compute vonNeumann entropy
			DO j = 1, localSize*localSize
	IF(evals(j).ne.0.0_rKind) THEN
        temp = temp + evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
	ELSE
	temp=temp+0.0_rKind
	END IF
			END DO
		entDist(l1,l2) = -1.0_rKind*temp ! Set value of entropy on-site.
END IF
entDist(l2,l1)=entDist(l1,l2)
		!Add the site=l1 local tensors to the cumulative Theta
		CALL ThetaNext(Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
			END DO
		END DO

	DEALLOCATE(rho2%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in LocalEntropyDist'
			END IF

	DEALLOCATE(workarr, rworkarr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
			END IF
DEALLOCATE(Theta)
END SUBROUTINE TwoBodyEntropyDist

SUBROUTINE PBphaseDist(phaseDist, Gammas, Lambdas)
!
!Purpose: Calculate the expectation of the Pegg-Barnett phase operator at each lattice site.
!spinless codes only
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: phaseDist(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind) :: dumPha
INTEGER :: i
	ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in PBphaseDist'
			END IF
	phaseDist = (/  (0.0_rKind, i = 1, systemSize) /)
		DO i = 1, systemSize
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			dumPha=TraceMatmul(PBphase_op%m, rho%m)
			phaseDist(i) = REAL(dumPha, KIND=rKind) ! Compute expectation value of phase operator.
		END DO
	DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in PBphaseDist'
			END IF
END SUBROUTINE PBphaseDist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Number conserving method starts !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!Spinless code only!!!!!!!!!!!!!!!!!!!
!!!!!!IDF to be supported in a later version!!!!!
SUBROUTINE ZetaKernelNC(Zeta, Lambda0, Gamma1, Lambda1, LabelL0, LabelL1)
!
!Purpose: Number conserving initial step in calculating the two-site density matrix
!Currently only supports spinless codes.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma1(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda0(:), Lambda1(:)
INTEGER, INTENT(IN) :: LabelL0(:), LabelL1(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: Zeta(:,:,:,:)
INTEGER :: NL1, NL2, NLlink, chi, i1, i2, alp1, alp2, beta
	chi=SIZE(Lambda0,1)
	Zeta=CMPLX(0.0,KIND=rKind)
	DO alp1=1,chi
		DO alp2=1,chi
			DO beta=1,chi
				NL1=LabelL0(alp1)
				NL2=LabelL0(alp2)
				NLlink=LabelL1(beta)
				IF((MAX(NL1,NL2,NLlink)<1000).AND.(NLlink-NL1>=0).AND.(NLlink-NL1<=maxFilling) &
				.AND.(NLlink-NL2>=0).AND.(NLlink-NL2<=maxFilling)) THEN
					Zeta(alp1,NLlink-NL1+1,NLlink-NL2+1,alp2) = &
					Zeta(alp1,NLlink-NL1+1,NLlink-NL2+1,alp2) + &
					Lambda0(alp1)*Gamma1(alp1,NLlink-NL1+1,beta)*Lambda1(beta)* &
					Lambda1(beta)*CONJG(Gamma1(alp2,NLlink-NL2+1,beta))*Lambda0(alp2)
				ELSE
				END IF
			END DO
		END DO
	END DO 
END SUBROUTINE ZetaKernelNC

	
SUBROUTINE ZetaNextNC(Zeta, Lambda0, Gamma1, LabelL0, LabelL1)
!
!Purpose: Number conserving recursive step in calculating the two-site density matrix
!Currently only supports spinless codes.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma1(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda0(:)
INTEGER, INTENT(IN) :: LabelL0(:), LabelL1(:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: Zeta(:,:,:,:)
COMPLEX(KIND=rKind) :: ZetaTemp(SIZE(Lambda0,1),localSize,localSize,localSize,SIZE(Lambda0,1))
INTEGER :: NL1, NL0, chi, i1, i2, alp1, alp2, beta1, beta2
	chi=SIZE(Lambda0,1)
	ZetaTemp=CMPLX(0.0,KIND=rKind)
	DO alp2=1,chi
		DO i2=1,localSize
			DO i1=1,localSize
				DO beta1=1,chi
					DO beta2=1,chi
						NL1=LabelL1(beta2)
						NL0=LabelL0(alp2)
						IF((MAX(NL1,NL0)<1000).AND.(NL1-NL0>=0).AND.(NL1-NL0<=maxFilling)) THEN
							ZetaTemp(beta1,i1,i2,NL1-NL0+1,alp2)= &
							ZetaTemp(beta1,i1,i2,NL1-NL0+1,alp2)+Zeta(beta1,i1,i2,beta2)* &
							CONJG(Gamma1(alp2,NL1-NL0+1,beta2))*Lambda0(alp2)
						ELSE
						END IF
					END DO
				END DO
			END DO
		END DO
	END DO
	Zeta=CMPLX(0.0,KIND=rKind)
	DO alp2=1,chi
		DO i2=1,localSize
			DO i1=1,localSize
				DO alp1=1,chi
					DO beta1=1,chi
						NL1=LabelL1(beta1)
						NL0=LabelL0(alp1)
						IF((MAX(NL1,NL0)<1000).AND.(NL1-NL0>=0).AND.(NL1-NL0<=maxFilling)) THEN
							Zeta(alp1,i1,i2,alp2)=Zeta(alp1,i1,i2,alp2)+Lambda0(alp1)* &
							Gamma1(alp1,NL1-NL0+1,beta1)*ZetaTemp(beta1,i1,i2,NL1-NL0+1,alp2)
						ELSE
						END IF
					END DO
				END DO
			END DO
		END DO
	END DO
END SUBROUTINE ZetaNextNC
	
SUBROUTINE TwoPointRhoNC(rho2, Zeta, Lambda0, Gamma1, LabelL0, LabelL1)
!
!Purpose: Calculate the two-site density matrix keeping number conservation
!Currently only supports spinless codes.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: rho2(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Zeta(:,:,:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Gamma1(:,:,:)
REAL(KIND=rKind), INTENT(IN) :: Lambda0(:)
INTEGER, INTENT(IN) :: LabelL0(:), LabelL1(:)
COMPLEX(KIND=rKind) :: ZetaTemp(SIZE(Lambda0,1),localSize,localSize,localSize,SIZE(Lambda0,1))
INTEGER :: NL0, NL1, chi, i1, i2, j2, alp1, alp2, beta1, beta2, k, l1
	chi=SIZE(Lambda0,1)
	ZetaTemp=CMPLX(0.0,KIND=rKind)
	DO alp2=1,chi
		DO i2=1,localSize
			DO i1=1,localSize
				DO beta1=1,chi
					DO beta2=1,chi
						NL1=LabelL1(beta2)
						NL0=LabelL0(alp2)
						IF((MAX(NL1,NL0)<1000).AND.(NL1-NL0>=0).AND.(NL1-NL0<=maxFilling)) THEN
							ZetaTemp(beta1,i1,i2,NL1-NL0+1,alp2)= &
							ZetaTemp(beta1,i1,i2,NL1-NL0+1,alp2)+Zeta(beta1,i1,i2,beta2)* &
							CONJG(Gamma1(alp2,NL1-NL0+1,beta2))*Lambda0(alp2)
						ELSE
						END IF
					END DO
				END DO
			END DO
		END DO
	END DO
			
	rho2=CMPLX(0.0,KIND=rKind)
	DO j2=1,localSize
		DO i2=1,localSize
			DO i1=1,localSize
				DO alp1=1,chi
					DO beta1=1,chi
						NL1=LabelL1(beta1)
						NL0=LabelL0(alp1)
						IF((MAX(NL1,NL0)<1000).AND.(NL1-NL0>=0).AND.(NL1-NL0<=maxFilling)) THEN
							rho2((NL1-NL0+1-1)*localSize+i1, (j2-1)*localSize+i2)= &
							rho2((NL1-NL0+1-1)*localSize+i1, (j2-1)*localSize+i2)+ &
							Lambda0(alp1)*Gamma1(alp1,NL1-NL0+1,beta1)*ZetaTemp(beta1,i1,i2,j2,alp1)
						ELSE
						END IF
					END DO
				END DO
			END DO
		END DO
	END DO
END SUBROUTINE TwoPointRhoNC

!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE TwoPointExpValNC !!!!!!!!!!!!!!!!!!!!!!!!
	
SUBROUTINE TwoPointExpValNC_r(observable, Op1, Op2, Lambdas, Gammas, LabelLeft)
!
!Purpose: Calculate the expectation value of the real two-site operator Op2 keeping number conservation
!Currently only supports spinless codes.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
COMPLEX(KIND=rKind) :: Zeta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
	observable=CMPLX(0.0,KIND=rKind)
	DO l1=1,systemSize
		CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
		observable(l1,l1)=TraceMatmul(Op1,rho1)
	END DO
			
	DO l2=systemSize,2,(-1)
		CALL ZetaKernelNC(Zeta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v, &
			LabelLeft(l2)%vi, LabelLeft(l2+1)%vi)
		DO l1=(l2-1),1,(-1)
			CALL TwoPointRhoNC(rho2, Zeta, Lambdas(l1)%v, Gammas(l1)%t, &
			LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
			observable(l1,l2)=TraceMatmul(Op2,rho2)
			observable(l2,l1)=CONJG(observable(l1,l2))
			CALL ZetaNextNC(Zeta, Lambdas(l1)%v, Gammas(l1)%t, LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
		END DO
	END DO
END SUBROUTINE TwoPointExpValNC_r
	
SUBROUTINE TwoPointExpValNC_c(observable, Op1, Op2, Lambdas, Gammas, LabelLeft)
!
!Purpose: Calculate the expectation value of the complex two-site operator Op2 keeping number conservation.
!Currently only supports spinless codes.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
COMPLEX(KIND=rKind) :: Zeta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
	observable=CMPLX(0.0,KIND=rKind)
	DO l1=1,systemSize
		CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
		observable(l1,l1)=TraceMatmul(Op1,rho1)
	END DO
			
	DO l2=systemSize,2,(-1)
		CALL ZetaKernelNC(Zeta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v, &
			LabelLeft(l2)%vi, LabelLeft(l2+1)%vi)
		DO l1=(l2-1),1,(-1)
			CALL TwoPointRhoNC(rho2, Zeta, Lambdas(l1)%v, Gammas(l1)%t, &
			LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
			observable(l1,l2)=TraceMatmul(Op2,rho2)
			observable(l2,l1)=CONJG(observable(l1,l2))
			CALL ZetaNextNC(Zeta, Lambdas(l1)%v, Gammas(l1)%t, LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
		END DO
	END DO
END SUBROUTINE TwoPointExpValNC_c


SUBROUTINE TwoPointExpValNC_rc(observable, Op1, Op2, Lambdas, Gammas, LabelLeft)
!
!Purpose: Calculate the expectation value of the real two-site operator Op2 keeping number conservation
!Currently only supports spinless codes.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op1(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
COMPLEX(KIND=rKind) :: Zeta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
	observable=CMPLX(0.0,KIND=rKind)
	DO l1=1,systemSize
		CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
		observable(l1,l1)=TraceMatmul(Op1,rho1)
	END DO
			
	DO l2=systemSize,2,(-1)
		CALL ZetaKernelNC(Zeta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v, &
			LabelLeft(l2)%vi, LabelLeft(l2+1)%vi)
		DO l1=(l2-1),1,(-1)
			CALL TwoPointRhoNC(rho2, Zeta, Lambdas(l1)%v, Gammas(l1)%t, &
			LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
			observable(l1,l2)=TraceMatmul(Op2,rho2)
			observable(l2,l1)=CONJG(observable(l1,l2))
			CALL ZetaNextNC(Zeta, Lambdas(l1)%v, Gammas(l1)%t, LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
		END DO
	END DO
END SUBROUTINE TwoPointExpValNC_rc
	
SUBROUTINE TwoPointExpValNC_cr(observable, Op1, Op2, Lambdas, Gammas, LabelLeft)
!
!Purpose: Calculate the expectation value of the complex two-site operator Op2 keeping number conservation.
!Currently only supports spinless codes.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
COMPLEX(KIND=rKind) :: Zeta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2
	observable=CMPLX(0.0,KIND=rKind)
	DO l1=1,systemSize
		CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
		observable(l1,l1)=TraceMatmul(Op1,rho1)
	END DO
			
	DO l2=systemSize,2,(-1)
		CALL ZetaKernelNC(Zeta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v, &
			LabelLeft(l2)%vi, LabelLeft(l2+1)%vi)
		DO l1=(l2-1),1,(-1)
			CALL TwoPointRhoNC(rho2, Zeta, Lambdas(l1)%v, Gammas(l1)%t, &
			LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
			observable(l1,l2)=TraceMatmul(Op2,rho2)
			observable(l2,l1)=CONJG(observable(l1,l2))
			CALL ZetaNextNC(Zeta, Lambdas(l1)%v, Gammas(l1)%t, LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
		END DO
	END DO
END SUBROUTINE TwoPointExpValNC_cr

!!!!!!!!!!!!!!!!!!! END CONTENTS OF INTERFACE TwoPointExpValNC !!!!!!!!!!!!!!!!!!!!!!!!

	
SUBROUTINE TotalEnergyNC(energy, H, Gammas, Lambdas, LabelLeft)
!
!Purpose: Calculate the energy eigenvalue associated with the Hamiltonian H keeping number conservation
!Currently only supports spinless codes.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
COMPLEX(KIND=rKind) :: Zeta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3))
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: l1, l2
	energy=0.0_rKind
	DO l2=2,systemSize
		CALL ZetaKernelNC(Zeta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v, &
			LabelLeft(l2)%vi, LabelLeft(l2+1)%vi)
		l1=l2-1
		CALL TwoPointRhoNC(rho2, Zeta, Lambdas(l1)%v, Gammas(l1)%t, &
			LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
			dumEn=TraceMatmul(H(l1)%m,rho2)
		energy=energy+REAL(dumEn, KIND=rKind)
	END DO
END SUBROUTINE TotalEnergyNC

SUBROUTINE LocalEnergyNC(energy, H, Gammas, Lambdas, LabelLeft)
!
!Purpose: Calculate the energy associated with each lattice bond keeping number conservation
!Currently only supports spinless codes.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy(systemSize-1)
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:)
COMPLEX(KIND=rKind) :: Zeta(SIZE(Gammas(1)%t, 3), localSize, localSize, SIZE(Gammas(1)%t, 3))
COMPLEX(KIND=rKind) :: rho2(localSize*localSize, localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: l1, l2	
	energy = 0.0_rKind
	DO l2 = 2, systemSize
		CALL ZetaKernelNC(Zeta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v, LabelLeft(l2)%vi, LabelLeft(l2+1)%vi)
		l1 = l2 - 1
		CALL TwoPointRhoNC(rho2, Zeta, Lambdas(l1)%v, Gammas(l1)%t, LabelLeft(l1)%vi, LabelLeft(l1+1)%vi)
		dumEn=TraceMatmul(H(l1)%m,rho2)
		energy(l1) = REAL(dumEn, KIND=rKind)
	END DO
END SUBROUTINE LocalEnergyNC


SUBROUTINE AllocateMeasures(Measures,numLocal, numAvg, numCorr, numFermiCorr, numEnt)
!
!Purpose: Allocate the measure derived type to hold num* measures of * kind
!
IMPLICIT NONE
TYPE(measure) :: Measures
INTEGER, INTENT(IN) :: numLocal, numAvg, numCorr, numFermiCorr, numEnt
INTEGER :: i

IF(numlocal.gt.0) THEN
	ALLOCATE(Measures%local(numLocal))
	DO i=1,numLocal
		ALLOCATE(Measures%local(i)%Op(localSize,localSize))
		ALLOCATE(Measures%local(i)%value(systemSize))
	END DO
	Measures%localpres=.TRUE.
ELSE
	Measures%localpres=.FALSE.
END IF

IF(numavg.gt.0) THEN
	ALLOCATE(Measures%avg(numAvg))
	DO i=1,numAvg
		ALLOCATE(Measures%avg(i)%Op(localSize,localSize))
	END DO
	Measures%avgpres=.TRUE.
ELSE
	Measures%avgpres=.FALSE.
END IF

SELECT CASE (numEnt)
CASE (1)
	Measures%entpres=1
	ALLOCATE(Measures%ent%vN(systemSize))
	ALLOCATE(Measures%ent%chain(systemSize+1))
CASE (2)
	Measures%entpres=2
	ALLOCATE(Measures%ent%vN(systemSize))
	ALLOCATE(Measures%ent%chain(systemSize+1))
	ALLOCATE(Measures%ent%tbvN(systemSize,systemSize))
CASE DEFAULT
	Measures%entpres=0
END SELECT

IF(numCorr.gt.0) THEN
	ALLOCATE(Measures%corr(numCorr))
	DO i=1,numCorr
		ALLOCATE(Measures%corr(i)%Op(localSize*localSize,localSize*localSize))	
		ALLOCATE(Measures%corr(i)%value(systemSize,systemSize))
	END DO
	Measures%corrpres=.TRUE.
ELSE
	Measures%corrpres=.FALSE.
END IF

IF(numfermicorr.gt.0) THEN
	ALLOCATE(Measures%fermicorr(numFermiCorr))
	DO i=1,numFermiCorr
		ALLOCATE(Measures%fermicorr(i)%Op(localSize*localSize,localSize*localSize))	
		ALLOCATE(Measures%fermicorr(i)%value(systemSize,systemSize))
	END DO
	Measures%fermicorrpres=.TRUE.
ELSE
	Measures%fermicorrpres=.FALSE.
END IF

END SUBROUTINE AllocateMeasures

SUBROUTINE DeallocateMeasures(Measures)
!
!Purpose: Deallocate the measure derived type
!
IMPLICIT NONE
TYPE(measure) :: Measures
INTEGER :: i

IF(Measures%localpres) THEN
	DO i=1,SIZE(Measures%local)
		DEALLOCATE(Measures%local(i)%Op)
		DEALLOCATE(Measures%local(i)%value)
	END DO
	DEALLOCATE(Measures%local)
END IF


IF(Measures%avgpres) THEN
	DO i=1,SIZE(Measures%avg)
		DEALLOCATE(Measures%avg(i)%Op)
	END DO
	DEALLOCATE(Measures%avg)
END IF

SELECT CASE (Measures%entpres)
CASE (1)
	DEALLOCATE(Measures%ent%vN)
	DEALLOCATE(Measures%ent%chain)
CASE (2)
	DEALLOCATE(Measures%ent%vN)
	DEALLOCATE(Measures%ent%chain)
	DEALLOCATE(Measures%ent%tbvN)
CASE DEFAULT
END SELECT

IF(Measures%corrpres) THEN
	DO i=1,SIZE(Measures%corr)
		DEALLOCATE(Measures%corr(i)%value)	
		DEALLOCATE(Measures%corr(i)%Op)
	END DO
	DEALLOCATE(Measures%corr)
END IF

IF(Measures%fermicorrpres) THEN
	DO i=1,SIZE(Measures%fermicorr)
		DEALLOCATE(Measures%fermicorr(i)%value)	
		DEALLOCATE(Measures%fermicorr(i)%Op)
	END DO
	DEALLOCATE(Measures%fermicorr)
END IF

END SUBROUTINE DeallocateMeasures

SUBROUTINE EvaluateMeasures(Measures, Gammas, Lambdas, H)
!
!Purpose: Evaluate the measures stored in the measure derived type
!
IMPLICIT NONE
TYPE(measure) :: Measures
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix), POINTER :: H(:)
TYPE(matrix) :: rho
COMPLEX(KIND=rKind), ALLOCATABLE :: Theta(:,:,:,:),Thetafermi(:,:,:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
REAL(KIND=rKIND) :: temp
REAL(KIND=rKIND), ALLOCATABLE :: evals(:)
INTEGER :: i,j, l1,l2, numO
INTEGER :: workSize, info
COMPLEX(KIND=rKind), ALLOCATABLE :: workArr(:), rworkArr(:)


!!!!!!!!!!!!Local Observables!!!!!!!!!!!!!!
ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in EvaluateMeasures'
END IF

!Zero out average measures
IF(Measures%avgpres) THEN
	DO i=1,SIZE(Measures%avg)
		Measures%avg(i)%value=0.0_rKind
	END DO
END IF
!Initialize single-site entropies
IF(Measures%entpres.ge.1) THEN
	Measures%ent%qme=localSize*1.0_rKind/(localSize*1.0_rKind-1.0_rKind)
!Allocate workspace for ZHEEV
	workSize = 2*localSize-1
	ALLOCATE(workarr(workSize), rworkarr(3*localSize-2),evals(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
			END IF
END IF


!Construct the single-site density matrix at every site
DO i=1,systemSize,1
	rho%m=0.0_rKind
	CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)

!local measures
IF(Measures%localpres) THEN
	DO numO=1,SIZE(Measures%local)
		Measures%local(numO)%value(i)=0.0_rKind
		Measures%local(numO)%value(i)=TraceMatmul(Measures%local(numO)%Op,rho%m)
	END DO
END IF

!average measures
IF(Measures%avgpres) THEN
	DO numO=1,SIZE(Measures%avg)
		Measures%avg(numO)%value=Measures%avg(numO)%value+TraceMatmul(Measures%avg(numO)%Op,rho%m)/(systemSize*1.0_rKind)
	END DO

END IF


!!Diagonal elements of correlation functions
!IF(Measures%corrpres) THEN
!	DO numO=1,SIZE(Measures%corr)
!		Measures%corr(numO)%value(i,i)=TraceMatmul(MATMUL(Measures%corr(numO)%Op(1)%m,Measures%corr(numO)%Op(2)%m),rho%m)
!	END DO
!END IF

!!Diagonal elements of correlation functions with fermi phases
!IF(Measures%fermicorrpres) THEN
!	DO numO=1,SIZE(Measures%fermicorr)
!		Measures%fermicorr(numO)%value(i,i)=TraceMatmul(MATMUL(Measures%fermicorr(numO)%Op(1)%m,Measures%fermicorr(numO)%Op(2)%m),rho%m)
!	END DO
!END IF

!One-body entropies
IF(Measures%entpres.ge.1) THEN
	Measures%ent%qme=Measures%ent%qme-TraceMatmul(rho%m, rho%m)*localSize/(1.0_rKind*systemSize*(localSize*1.0_rKind-1.0_rKind))
	CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info)
	Measures%ent%vN(i)=0.0_rKind
	DO j = 1, localSize
		IF(evals(j).ne.0.0_rKind) THEN
        	Measures%ent%vN(i) = Measures%ent%vN(i)-evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
		ELSE
			Measures%ent%vN(i)=Measures%ent%vN(i)+0.0_rKind
		END IF
	END DO
END IF



IF(Measures%entpres==2) THEN
	Measures%ent%tbvN(i,i)=2.0_rKind*Measures%ent%vN(i)
END IF

END DO				

IF(Measures%entpres.ge.1) THEN
	DO i=1,systemSize+1
		Measures%ent%chain(i)=ChainEntropy(i,Lambdas)
	END DO
END IF
		
DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in EvaluateMeasures'
END IF

IF(Measures%entpres.ge.1) THEN
	DEALLOCATE(workarr, rworkarr,evals, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
			END IF
END IF

!!!!!!!!!!!!Correlation Functions, energy, and two-body entropy!!!!!!!!!!!!!!

IF(Measures%corrpres.or.Measures%fermicorrpres.or.(Measures%entpres==2)) THEN

Measures%en=0.0_rKind

ALLOCATE(rho%m(localSize*localSize,localSize*localSize), STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to allocate rho in EvaluateMeasures'
END IF

!Initialize two-site entropies
IF(Measures%entpres.ge.2) THEN
!Allocate workspace for ZHEEV
	workSize = 2*localSize*localSize-1
	ALLOCATE(workarr(workSize), rworkarr(3*localSize*localSize-2),evals(localSize*localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
			END IF
END IF

	ALLOCATE(Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))
IF(Measures%fermicorrpres) THEN
	ALLOCATE(Thetafermi(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))
	ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))
END IF
	Measures%en=0.0_rKind
		DO l2=systemSize,2,(-1)
		!Compute the initial theta (Theta at site l2)
			CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
		IF(Measures%fermicorrpres) THEN
			CALL ThetaKernal(Thetafermi,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
		END IF
			DO l1=(l2-1),1,(-1)
				!Form the two-site density matrix from the cumulative Theta
		CALL TwoSiteRho(rho%m,Theta,Gammas(l1)%t, Gammas(l1)%t,Lambdas(l1)%v)

	!Total energy
		IF(l1==l2-1) THEN
			Measures%en=Measures%en+TraceMatmul(H(l1)%m,rho%m)
		END IF
	!Boundary term
		IF((l2==systemSize).and.(l1==1).and.(BoundaryCond=='P')) THEN
			Measures%en=Measures%en+TraceMatmul(H(systemSize)%m,rho%m)
		END IF
	

!Phaseless correlation functions
		IF(Measures%corrpres) THEN
			DO i=1,SIZE(Measures%corr)
				Measures%corr(i)%value(l1,l2)=TraceMatmul(Measures%corr(i)%Op,rho%m)
				Measures%corr(i)%value(l2,l1)=CONJG(Measures%corr(i)%value(l1,l2))
			END DO
		END IF

!Two-site entropy		
		IF(Measures%entpres==2) THEN
	CALL ZHEEV('N', 'U', localSize*localSize, rho%m, localSize*localSize, evals, workArr, workSize, rworkArr, info)
		Measures%ent%tbvN(l1,l2)=0.0_rKind
		DO j = 1, localSize*localSize
			IF(evals(j).ne.0.0_rKind) THEN
		        	Measures%ent%tbvN(l1,l2) = Measures%ent%tbvN(l1,l2)-evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) 
			ELSE
				Measures%ent%tbvN(l1,l2)=Measures%ent%tbvN(l1,l2)+0.0_rKind
			END IF
		END DO
		Measures%ent%tbvN(l2,l1)=Measures%ent%tbvN(l1,l2)

		END IF

!Correlation functions with a fermi phase
		IF(Measures%fermicorrpres) THEN
				GammaP = Gammas(l1)%t
				CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				CALL TwoSiteRho(rho%m,Thetafermi,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
			DO i=1,SIZE(Measures%fermicorr)
				Measures%fermicorr(i)%value(l1,l2)=TraceMatmul(Measures%fermicorr(i)%Op,rho%m)
				Measures%fermicorr(i)%value(l2,l1)=CONJG(Measures%fermicorr(i)%value(l1,l2))
			END DO
		END IF

				!Add the site=l1 local tensors to the cumulative Theta
		CALL ThetaNext(Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
	
		IF(Measures%fermicorrpres) THEN
				CALL ThetaNext(Thetafermi,Gammas(l1)%t,GammaP,Lambdas(l1)%v)
		END IF
			END DO
		END DO


	DEALLOCATE(Theta)

IF(Measures%fermicorrpres) THEN
	DEALLOCATE(Thetafermi)
	DEALLOCATE(GammaP)
END IF

!Deallocate variables needed for two-site entropies
IF(Measures%entpres.ge.2) THEN
!Allocate workspace for ZHEEV
	DEALLOCATE(workarr, rworkarr,evals, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate workarr in LocalEntropyDist'
			END IF
END IF

DEALLOCATE(rho%m, STAT=statInt)
IF(statInt.ne.0) THEN
	PRINT *, 'Failed to deallocate rho in EvaluateMeasures'
END IF

!If no correlation functions or two-body entropy are needed, 
!just compute the reduced density matrices needed for the energy
ELSE

CALL TotalEnergy(Measures%en, H, Gammas, Lambdas)

END IF


END SUBROUTINE EvaluateMeasures
	
END MODULE observables_module
