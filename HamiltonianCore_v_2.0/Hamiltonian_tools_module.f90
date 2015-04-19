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
MODULE Hamiltonian_tools_module
!
! Purpose: Module which contains operators to define various Hamiltonians, 
!procedures to generate Hamiltonians,  procedures to generate fock spaces
!and initial states with and without internal degrees of freedom, and vector 
!coupling coefficients for OpenSourceTEBD v1.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   9/28/09   M. L. Wall	v2.0 release
!
USE system_parameters	
USE TEBDtools_module

IMPLICIT NONE

!Hubbard parameters
REAL(KIND=rKind) :: jTunn !Tunneling energy
REAL(KIND=rKind) :: U0 !on-site repulsion
REAL(KIND=rKind) :: mu0 !chemical potential
REAL(KIND=rKind) :: V0 !nearest-neighbor repulsion
REAL(KIND=rKind), ALLOCATABLE :: extPot(:) !Optional site-dependent potential


!Spin system parameters
REAL(KIND=rKind) :: spin !Spin of particles
INTEGER :: spinSize !Shorthand for 2*spin+1

! *** General operators ***
TYPE(matrixReal) :: one_op !Unit operator

! *** Spinless operators ***
TYPE(matrixReal) :: a_op !Destruction operator
TYPE(matrixReal) :: t_op !Tunneling operator
TYPE(matrix) :: PBphase_op !Pegg-Barnett phase operator.

! *** Fermionic phase operator ***
TYPE(matrixReal) :: fermiPhase_op

! *** Internal degrees of freedom variables ***
TYPE(vectorInt) :: Conserv !Conserv holds an Abelian conserved quantity
REAL(KIND=rkind) :: lFac(200) !Array containing the logarithms of the first 200 factorials

! *** Internal degrees of freedom operators ***
!The list index is the internal degrees of freedom index
TYPE(matrixReal), POINTER :: a_opS(:) !Destruction operator.  Index is idof component
TYPE(matrixReal) :: Sx_opS !Operator for the x-component of Spin
TYPE(matrix) :: Sy_opS !Operator for the y-component of Spin
TYPE(matrixReal) :: Sz_opS !Operator for the z-component of Spin
TYPE(matrixReal) :: Ssq_opS !Total spin squared operator
TYPE(matrixReal) :: ntot_opS !Total nuber operator

! *** Spinor operators ***
TYPE(matrixReal) :: VB_opS !Zeeman interaction operator
TYPE(matrixReal) :: ttot_opS !Tunneling operator


! Interfaces for generic procedures
INTERFACE HamiOneSite
MODULE PROCEDURE HamiOneSite_r, HamiOneSite_c
END INTERFACE HamiOneSite

INTERFACE HamiLeft
MODULE PROCEDURE HamiLeft_r, HamiLeft_c
END INTERFACE HamiLeft

INTERFACE HamiRight
MODULE PROCEDURE HamiRight_r, HamiRight_c
END INTERFACE HamiRight

CONTAINS

SUBROUTINE ProductStateMPD(Gammas, Lambdas, carray)
!
!Purpose: Construct the Vidal decomposition of a product state whose coefficients
!		  are stored in carray
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: i, j, m, d, chi
COMPLEX(KIND=rKind), INTENT(IN) :: carray(:, :)
	m = SIZE(carray, 2)
	d = SIZE(carray, 1)

	!Check to be sure that M and d read in agree with the global variables
	IF (m /= systemSize) THEN
	PRINT *, 'systemSize parameter conflicts with input data in ProductStateMPD' 
	END IF
	
	IF (d /= localSize) THEN
	PRINT *, 'localSize parameter conflicts with input data in ProductStateMPD'
	END IF
	
	DO i = 1, systemSize
	Gammas(i)%t = CMPLX(0.0, KIND=rKind)
	Lambdas(i)%v = 0.0_rKind
	Lambdas(i)%v(1) = 1.0_rKind ! Assign the first component of each lambda the value 1, as this is a product state.
		DO j = 1, localSize
		Gammas(i)%t(1, j, 1) = carray(j, i) ! Assign the alpha=1 values of Gammas to be the coefficients of each on-site state.
		END DO
	END DO
	Lambdas(systemSize+1)%v = 0.0_rKind
	Lambdas(systemSize+1)%v(1) = 1.0_rKind
END SUBROUTINE ProductStateMPD


SUBROUTINE ProductStateLabels(LabelLeft, LabelRight, carray,intDegFree)
!
!Purpose: Construct the lists of number conserving labels of a product state whose coefficients
!		  are stored in carray
!
!WARNING: The input state should be an eigenstate of total number.  This routine will not check for this!
!
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: i, j, m, d, chi
COMPLEX(KIND=rKind), INTENT(IN) :: carray(:, :)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
	m = SIZE(carray, 2)
	d = SIZE(carray, 1)

	!Check to be sure that M and d read in agree with the global variables
	IF (m /= systemSize) THEN
	PRINT *, 'systemSize parameter conflicts with input data in ProductStateMPD' 
	END IF
	
	IF (d /= localSize) THEN
	PRINT *, 'localSize parameter conflicts with input data in ProductStateMPD'
	END IF

!If internal degrees of freedom are present, find the occupation of each component
IF(PRESENT(intDegFree)) THEN
LabelLeft(1)%vi(1)=0.0_rKind
DO i=1,systemSize
	DO j=1,localSize
		IF(ABS(carray(j,i)).ne.0.0_rKind) THEN
		LabelLeft(i+1)%vi(1)=LabelLeft(i)%vi(1)+Conserv%vi(j)
		EXIT
		END IF
	END DO
END DO

ELSE
LabelLeft(1)%vi(1)=0.0_rKind
DO i=1,systemSize
	DO j=1,localSize
		IF(ABS(carray(j,i)).ne.0.0_rKind) THEN
		LabelLeft(i+1)%vi(1)=LabelLeft(i)%vi(1)+j-1
		END IF
	END DO
END DO

END IF

!Construct LabelRight from LabelLeft
DO i=1,(systemSize+1),1
	LabelRight(i)%vi(1)=totNum-LabelLeft(i)%vi(1)
END DO
	
END SUBROUTINE ProductStateLabels


SUBROUTINE AllStates(Gammas, Lambdas)
!
!Purpose: Creates an initial state that is a product of local states 
! which contain all possible states in the same amount.  Used as initial ITP
! state for number non-conserving code.
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER :: i, j
	DO i=1,systemSize
	Gammas(i)%t=CMPLX(0.0, KIND=rKind)
	Lambdas(i)%v=0.0_rKind
	Lambdas(i)%v(1)=1.0_rKind
	!Each state is weighted equally by normalization by 1/sqrt(d)
		DO j=1,localSize
		Gammas(i)%t(1,j,1)=CMPLX((1.0_rKind)/SQRT(localSize*1.0_rKind),KIND=rKind)
		END DO
	END DO
	Lambdas(systemSize+1)%v=0.0_rKind
	Lambdas(systemSize+1)%v(1)=1.0_rKind
END SUBROUTINE AllStates


SUBROUTINE onsiteStateListIdof(list, idofSize, fermiSwitch)
!
!Purpose: Outer routine of a subroutine that returns the list of on-site states for  
!an idofSize-size internal Hilbert space truncated at maxfilling particles per site.
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: list(:, :)
INTEGER, INTENT(INOUT) :: idofSize
INTEGER i,j,n(idofSize),counter,k,l, dum3, dum4
INTEGER, INTENT(INOUT), OPTIONAL :: fermiSwitch

	list = 0.0_rKind
	n=0.0_rKind
	counter=0
IF(PRESENT(fermiSwitch)) THEN
!Loop over total number of particles nmax, beginning with nMax=0	
	DO k=1,maxFilling+1,1
	
	n=0.0_rKind
	n(1)=k-1
!FORTRAN calls by reference, and so all
!objects passed into a recursive function
!need to be variables-NOT indicies
		dum3=k-1
		dum4=1
		CALL onsiteIdofInner(list,dum3,counter,dum4,n, idofSize, fermiSwitch)
	END DO
ELSE
!Loop over total number of particles nmax, beginning with nMax=0	
	DO k=1,maxFilling+1,1
	
	n=0.0_rKind
	n(1)=k-1
!FORTRAN calls by reference, and so all
!objects passed into a recursive function
!need to be variables-NOT indicies
		dum3=k-1
		dum4=1
		CALL onsiteIdofInner(list,dum3,counter,dum4,n, idofSize)
	END DO
END IF
END SUBROUTINE onsiteStateListIdof

RECURSIVE SUBROUTINE onsiteIdofInner(list, nmax,  counter, m,n, idofSize, fermiSwitch)
!
!Purpose: Inner routine of a subroutine that returns the list of on-site states for  
!an idofSize-size internal Hilbert space truncated at maxfilling particles per site.
!
!See manual for more detail
!
IMPLICIT NONE
!All passed variables in a recursive function HAVE to be INOUTs
INTEGER, INTENT(INOUT) :: nmax
INTEGER, INTENT(INOUT) :: m
INTEGER, INTENT(INOUT) :: counter
INTEGER, INTENT(INOUT) :: idofSize
INTEGER, INTENT(INOUT) :: n(:)
COMPLEX(KIND=rKind), INTENT(INOUT) :: list(:, :)
INTEGER, INTENT(INOUT), OPTIONAL :: fermiSwitch
INTEGER i,j,ntemp(idofSize),k,l,dum1, dum2

!!!!!!!!!!!!!!!!Fermi Statistics!!!!!!!!!!!!!!
IF(PRESENT(fermiSwitch)) THEN

!If the state is stuff,0,0,...,0 just count it
IF (nmax==0) THEN
	IF(MAXVAL(n).le.1) THEN
		counter=counter+1
		list(counter, :)=n(:)
	END IF
!If the state is stuff,1,0,...,0 then
!move the 1 over sequentially, counting
!each time
ELSE IF (nmax==1) THEN
	IF(MAXVAL(n).le.1) THEN
		DO i=m,idofSize,1
!define a temporary array
		ntemp=n
!move a particle out of the ith state into
!the mth state
		ntemp(m)=ntemp(m)-1
		ntemp(i)=ntemp(i)+1
		counter=counter+1
		list(counter, :)=ntemp(:)
		END DO
!If the state is stuff,k,0,...,0 then
!perform the algorithm on the subspace beginning
!with k-see manual for more detail
	END IF
ELSE
!Loop over putting all k particles in the ith component
	DO i=m,idofSize,1
!Zero out the subspace array	
		DO j=m,idofSize,1
		n(j)=0
		END DO
		n(i)=nmax
!Loop over recursive calls
		DO l=1,nmax,1
		ntemp=n
		ntemp(i)=ntemp(i)+1-l
		ntemp(i+1)=ntemp(i+1)+l-1
!FORTRAN calls by reference, and so all
!things passed into a recursive function
!need to be variables-NOT indices
		dum1=l-1
		dum2=i+1
		CALL onsiteIdofInner(list,dum1,counter,dum2,ntemp, idofSize, fermiSwitch)
		END DO
	END DO
END IF

!!!!!!!!!!!!!!!!Bose Statistics!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE
!If the state is stuff,0,0,...,0 just count it
IF (nmax==0) THEN
	counter=counter+1
	list(counter, :)=n(:)
!If the state is stuff,1,0,...,0 then
!move the 1 over sequentially, counting
!each time
ELSE IF (nmax==1) THEN

		DO i=m,idofSize,1
!define a temporary array
		ntemp=n
!move a particle out of the ith state into
!the mth state
		ntemp(m)=ntemp(m)-1
		ntemp(i)=ntemp(i)+1
		counter=counter+1
		list(counter, :)=ntemp(:)
		END DO
!If the state is stuff,k,0,...,0 then
!perform the algorithm on the subspace beginning
!with k-see manual for more detail
ELSE

!Loop over putting all k particles in the ith component
	DO i=m,idofSize,1
!Zero out the subspace array	
		DO j=m,idofSize,1
		n(j)=0
		END DO
		n(i)=nmax
!Loop over recursive calls
		DO l=1,nmax,1
		ntemp=n
		ntemp(i)=ntemp(i)+1-l
		ntemp(i+1)=ntemp(i+1)+l-1
!FORTRAN calls by reference, and so all
!things passed into a recursive function
!need to be variables-NOT indices
		dum1=l-1
		dum2=i+1
		CALL onsiteIdofInner(list,dum1,counter,dum2,ntemp, idofSize)
		END DO
	END DO
END IF
END IF
END SUBROUTINE onsiteIdofInner

SUBROUTINE InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight, intDegFree)
!
!Purpose: Creates an initial state consistent with number conservation
!
!		  The algorithm places the particles in a "wedding cake" structure
!		  with the greatest number of particles in the center of the cake ("tops")
!		  A lesser number surrounding this peak ("center"), and no particles in the gap between
!		  the cake and the edge ("hole").  This mimics the ground state in a harmonic trap.
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!		  See the manual for details/examples concerning the algorithm
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
REAL(KIND=rKind) :: minc, maxc
INTEGER :: i, j, k, l, hole, tops, center,nj

IF(totNum==0) THEN
	DO i=1,(systemSize+1),1
		LabelLeft(i)%vi(1)=0
		LabelRight(i)%vi(1)=0
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	DO i=1,systemSize
		Gammas(i)%t(1,1,1)=1.0_rKind
	END DO


ELSE
! If the number of sites with the least number of particles (hole) is even
	If((MOD(systemSize-MOD(totNum,systemSize),2)==0).OR. &
! or there is no hole (integer filling)
	(MOD(totNum,systemSize)==0)) THEN
		!!If we have unit filling there is no hole
		If(MOD(totNum,systemSize)==0) THEN
		hole=0
		!for non-unit filling the hole is
		ELSE
		hole=(systemSize-MOD(totNum,systemSize))/2
		END IF
! Else if the number of holes is odd
	ELSE
	hole=(systemSize-MOD(totNum,systemSize)+1)/2
	END IF
				
! Number of sites that have the greatest number of particles
! number in the "top of the wedding cake"
	tops=MOD(totNum-1,systemSize)+1
		
! Number of sites that have the lesser number of particles
! number in the "center of the wedding cake"
!Floor ensures that when N=M we fill uniformly
	center=FLOOR(REAL(totNum,KIND=8)/REAL(systemSize,KIND=8)-10.0**(-8))+1

! LabelLeft in 0<=link<=hole
	DO i=1,(hole+1),1
		IF(i==1) THEN
		!There are never any particles to the left of the system
		LabelLeft(i)%vi(1)=0
		ELSE
		!Count the cumulative number to the left
		LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center-1
		END IF
	END DO
	
! LabelLeft in hole+1<=link<=hole+top
	DO i=(hole+2),(hole+tops+1),1
		LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center
	END DO
		
! LabelLeft in hole+top+1<=link<=systemSize+1
	DO i=(hole+tops+2),(systemSize+1),1
		LabelLeft(i)%vi(1)=LabelLeft(i-1)%vi(1)+center-1
	END DO
! Construct LabelRight from LabelLeft
	DO i=1,(systemSize+1),1
		LabelRight(i)%vi(1)=totNum-LabelLeft(i)%vi(1)
	END DO
						
! Construct Lambdas
	DO i=1,(systemSize+1),1
		Lambdas(i)%v(1)=1.0_rKind
	END DO


!Internal degree(s) of freedom present	
IF(PRESENT(intDegFree)) THEN

!Find the number of states with number=center-1 and center particles
!Store these in minc and maxc
		minc=0.0_rKind
		maxc=0.0_rKind
		DO j=1,localSize,1
		nj=Conserv%vi(j)
			IF(nj==center-1) THEN
			minc=minc+1.0_rKind
			ELSE IF(nj==center) THEN
			maxc=maxc+1.0_rKind
			END IF
		END DO

! Construct Gammas
! Gammas in 1<=site<=hole
		DO i=1,hole,1
!Sum over internal degree of freedom, weighting each democratically
		DO j=1,localSize,1
		nj=Conserv%vi(j)
		IF(nj==center-1) THEN
			Gammas(i)%t(1,j,1)=CMPLX(1.0/SQRT(minc),KIND=rKind)
		END IF
		END DO		
		END DO
! Gammas in hole+1<=site<=hole+top
		DO i=hole+1,(hole+tops),1
!Sum over internal degree of freedom, weighting each democratically
		DO j=1,localSize,1
		nj=Conserv%vi(j)
		IF(nj==center) THEN
			Gammas(i)%t(1,j,1)=CMPLX(1.0/SQRT(maxc),KIND=rKind)
		END IF
		END DO	
		END DO
! Gammas in hole+top+1<=site<=systemSize
		DO i=(hole+tops+1),systemSize,1
!Sum over internal degree of freedom, weighting each democratically
		DO j=1,localSize,1
		nj=Conserv%vi(j)
		IF(nj==center-1) THEN
			Gammas(i)%t(1,j,1)=CMPLX(1.0/SQRT(minc),KIND=rKind)
		END IF
		END DO		
		END DO

!Internal degree(s) of freedom absent
ELSE

! Construct Gammas
! Gammas in 1<=site<=hole
	DO i=1,hole,1
		Gammas(i)%t(1,center,1)=CMPLX(1.0,KIND=rKind)
	END DO
! Gammas in hole+1<=site<=hole+top
	DO i=hole+1,(hole+tops),1
		Gammas(i)%t(1,center+1,1)=CMPLX(1.0,KIND=rKind)
	END DO
! Gammas in hole+top+1<=site<=systemSize
	DO i=(hole+tops+1),systemSize,1
		Gammas(i)%t(1,center,1)=CMPLX(1.0,KIND=rKind)
	END DO
END IF
END IF		
END SUBROUTINE InitialSetNC

SUBROUTINE SetupLogFac()
!
!Purpose: Initialize lFac by finding the logarithms of the first 200 factorials
!         lFac(i)=Log((i-1)!).  This is used to make computation of vector coupling
!		  coefficients more efficient.
!
IMPLICIT NONE

INTEGER :: i
INTEGER, PARAMETER :: N=200
lFac(1)=0.0_rKind
	DO i=2,N,1
	lFac(i)=lFac(i-1)+LOG((i-1)*1.0_rKind)
	END DO
END SUBROUTINE SetupLogFac

REAL(KIND=rKind) FUNCTION LogTriCoef(X1,X2,X3)
!
!Purpose: This function evaluates the log of the
!triangle coefficient as defined in the manual.                                             
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: X1,X2,X3  
LogTriCoef=0.5_rKind*(lFac(FLOOR(X1+X2-X3+1))+lFac(FLOOR(X1-X2+X3+1))+lFac(FLOOR(X2-X1+X3+1))-lFac(FLOOR(X1+X2+X3+2)))

END FUNCTION LogTriCoef     
	
REAL(KIND=rKind) FUNCTION TriTest(J1,J2,J3)
!
!Purpose: This function returns 1 if triad J1, J2, J3 satisfies the triangle inequalities
! zero if it doesn't.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: J1,J2,J3
TriTest=0.0_rKind

IF((J1.ge.ABS(J2-J3)).AND.(J1.le.(J2+J3))) THEN
TriTest=1.0_rKind
END IF
END FUNCTION TriTest																																																										 

REAL(KIND=rKind) FUNCTION IntTest(J1,J2,J3)
!
!Purpose: This function returns 1 if triad J1, J2, J3 satisfies the integer rule
! zero if it doesn't
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: J1,J2,J3
IntTest=0.0_rKind

IF(FLOOR(J1+J2+J3)==CEILING(J1+J2+J3)) THEN
IntTest=1.0_rKind
END IF

END FUNCTION IntTest																																																										 

REAL(KIND=rKind) FUNCTION MTest(M1,J1)
!
!Purpose: This function returns 1 if -j<=m<=j and m+j is an integer
!zero otherwise
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: M1,J1
MTest=0.0_rKind

IF((M1.le.ABS(J1)).AND.(FLOOR(M1+J1)==CEILING(M1+J1))) THEN
MTest=1.0_rKind
END IF

END FUNCTION MTest
																																																																																																				
REAL(KIND=rKind) FUNCTION tIndTJ(X1,X2,X3,Y1,Y2,Y3)
!
!Purpose: This function computes the log of the t-independent part of the Racah formula
!for the three-J coefficient.  See the manual for more detail.                                              
!
IMPLICIT NONE

REAL(KIND=rKind), INTENT(IN) :: X1,X2,X3,Y1,Y2,Y3   

tIndTJ = LogTriCoef(X1,X2,X3) &
		+0.5_rKind*(lFac(FLOOR(X1+Y1+1))+lFac(FLOOR(X1-Y1+1))+lFac(FLOOR(X2+Y2+1))+lFac(FLOOR(X2-Y2+1))+lFac(FLOOR(X3+Y3+1))+lFac(FLOOR(X3-Y3+1)))

END FUNCTION tIndTJ

REAL(KIND=rKind) FUNCTION tIndSJ(X1,X2,X3,Y1,Y2,Y3)
!
!Purpose: This function computes the log of the t-independent part of the Racah formula
!for the six-J coefficient.  See the manual for more detail.                                          
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: X1,X2,X3,Y1,Y2,Y3   
tIndSJ = LogTriCoef(X1,X2,X3)+LogTriCoef(X1,Y2,Y3)+LogTriCoef(Y1,X2,Y3)+LogTriCoef(Y1,Y2,X3)

END FUNCTION tIndSJ

REAL(KIND=rKind) FUNCTION ThreeJ(J1D,M1D,J2D,M2D,J3D,M3D)                            
!
!Purpose: This function computes the Wigner-3J coefficient                                               
!(J1D/2  J2D/2  J3D/2)
!(M1D/2  M2D/2  M3D/2)
!using the Racah formula.  See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: J1D,M1D,J2D,M2D,J3D,M3D   
INTEGER :: tMIN,tMIN1,tMIN2,tMAX,tMAX1,tMAX2,tMAX3, t
REAL(KIND=rKind) :: X1, X2, X3, Y1, Y2, Y3 
REAL(KIND=rKind) :: term, term2, term1 
                 
X1 = J1D*0.5_rKind                                                    
X2 = J2D*0.5_rKind                                                   
X3 = J3D*0.5_rKind                                                   
Y1 = M1D*0.5_rKind                                                   
Y2 = M2D*0.5_rKind                                                   
Y3 = M3D*0.5_rKind
                                                
!Triangularity checks
if_triangularity:   IF((TriTest(X1,X2,X3)==0).OR.(M1D+M2D+M3D.NE.0).OR.(MTest(Y1,X1)==0).OR.(MTest(Y2,X2)==0).OR.(MTest(Y3,X3)==0)) THEN
					ThreeJ = 0.0_rKind 
					ELSE                       


!The t-sum includes all factorials with non-negative arguments
!This finds the smallest and largest values

!Factorials with +t
tMIN1=(J3D-J2D+M1D)/2
tMIN2=(J3D-J1D-M2D)/2
tMIN=MIN(tMIN1,tMIN2)

!If the smallest value allowed by the above terms is greater than zero, then
!the lower bound is 0 from the t! term
	IF(tMIN.ge.0) THEN
	tMIN=0
	ELSE
!If the smallest value is negative, we choose the least t such that
!t+KMIN=0
	tMIN=(-1)*tMIN
	END IF

!-t terms
tMAX1=FLOOR(X1+X2-X3)
tMAX2=FLOOR(X1-Y1)
tMAX3=FLOOR(X2+Y2)
tMAX=MIN(tMAX1,tMAX2,tMAX3)
        
!compute the t-independent part of the Racah formula                                                                    
      term1 = tIndTJ(X1,X2,X3,Y1,Y2,Y3)                                 

!Compute the summation over t
term = 0.0_rKind
DO t = tMIN,tMAX,1                                             
	term2 = lFac(t+1)+lFac(tMIN1+t+1)+lFac(tMIN2+t+1)+lFac(tMAX1-t+1)+lFac(tMAX2-t+1)+lFac(tMAX3-t+1)
	term = term + ((-1)**(IABS((J1D-J2D-M3D)/2)))*EXP(term1-term2)*((-1)**(t))                     
END DO

     ThreeJ=term    
					END IF if_triangularity  
					     
END FUNCTION ThreeJ  

REAL(KIND=rKind) FUNCTION Clebsch(J1D,M1D,J2D,M2D,JD,MD)    
!
!Purpose: This function computes the Clebsch-Gordon coefficient                                               
!<J1D/2, M1D/2; J2D/2, M2D/2| J3D/2, M3D/2>   
!using the 3-J coefficient.  See manual for more detail
!
IMPLICIT NONE  
                      
INTEGER, INTENT(IN) :: J1D,M1D,J2D,M2D,JD,MD
REAL(KIND=rKind) ::  Q, PHASE

Clebsch = 0.0_rKind                                                     

!Triangularity check
IF(M1D+M2D.NE.MD) THEN 
CLEBSCH = 0.0_rKind          
ELSE                                                                     
      Q = ThreeJ(J1D,M1D,J2D,M2D,JD,-MD) 
      PHASE = ((-1)**( (J2D-J1D-MD)/2))  
      Clebsch = Q*PHASE*SQRT((JD+1.0)*1.0_rKind)      
END IF                                                          
END FUNCTION Clebsch                                                          


REAL(KIND=rKind) FUNCTION SixJ(J11D,J21D,J12D,J22D,J13D,J23D)                            
!
!Purpose: This function computes the Wigner-6J coefficient                                               
!{J11D/2  J12D/2  J13D/2}
!{J21D/2  J22D/2  J23D/2}
!using the Racah formula.  See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: J11D,J21D,J12D,J22D,J13D,J23D  
REAL(KIND=rKind) :: X1, X2, X3, Y1, Y2, Y3, term1,term2,term
INTEGER :: t, tMIN, tMIN1, tMIN2, tMIN3, tMIN4, tMAX, tMAX1, tMAX2, tMAX3

!Convert to actual values
X1 = J11D*0.5_rKind                                                    
X2 = J12D*0.5_rKind                                                   
X3 = J13D*0.5_rKind                                                   
Y1 = J21D*0.5_rKind                                                   
Y2 = J22D*0.5_rKind                                                   
Y3 = J23D*0.5_rKind

!Selection rules
IF((TriTest(X1,X2,X3).ne.0.0_rKind).AND.(TriTest(X1,Y2,Y3).ne.0.0_rKind).AND.(TriTest(Y1,X2,Y3).ne.0.0_rKind).AND.(TriTest(Y1,Y2,X3).ne.0.0_rKind)&
	.AND.(IntTest(X1,X2,X3).ne.0.0_rKind).AND.(IntTest(X1,Y2,Y3).ne.0.0_rKind).AND.(IntTest(Y1,X2,Y3).ne.0.0_rKind).AND.(IntTest(Y1,Y2,X3).ne.0.0_rKind)) THEN


!Minimum allowed summation index
tMIN1=(J11D+J12D+J13D)/2
tMIN2=(J11D+J22D+J23D)/2
tMIN3=(J21D+J12D+J23D)/2
tMIN4=(J21D+J22D+J13D)/2
tMIN=MAX(tMIN1,tMIN2,tMIN3,tMIN4)

!Maximum allowed summation index
tMAX1=(J11D+J12D+J21D+J22D)/2
tMAX2=(J12D+J13D+J22D+J23D)/2
tMAX3=(J13D+J11D+J23D+J21D)/2
tMAX=MIN(tMAX1,tMAX2,tMAX3)

!compute the t-independent part of the Racah formula                                                                    
term1=tIndSJ(X1,X2,X3,Y1,Y2,Y3)

!Compute the summation over t
term = 0.0_rKind
DO t = tMIN,tMAX,1                                             
	term2 = lFac(tMAX1-t+1)+lFac(tMAX2-t+1)+lFac(tMAX3-t+1) &
			+ lFac(t-tMIN1+1)+lFac(t-tMIN2+1)+lFac(t-tMIN3+1)+lFac(t-tMIN4+1)
	term = term + EXP(lFAC(t+2)+term1-term2)*((-1)**(t))                     
END DO

SixJ=term    

ELSE
SixJ=0.0_rKind
END IF


END FUNCTION SixJ


REAL(KIND=rKind) FUNCTION NineJ(J11D,J21D,J31D,J12D,J22D,J32D,J13D,J23D,J33D)                            
!
!Purpose: This function computes the Wigner-9J coefficient                                               
!{J11D/2  J12D/2  J13D/2}
!{J21D/2  J22D/2  J23D/2}
!{J31D/2  J32D/2  J33D/2}
!using the 6-j contraction formula.  See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: J11D,J21D,J31D,J12D,J22D,J32D,J13D,J23D,J33D  
INTEGER :: t, tMIN, tMAX
REAL(KIND=rKind) :: X1, X2, X3, Y1, Y2, Y3, Z1,Z2,Z3,term2,term

!Convert to actual values
X1 = J11D*0.5_rKind                                                    
X2 = J12D*0.5_rKind                                                   
X3 = J13D*0.5_rKind                                                   
Y1 = J21D*0.5_rKind                                                   
Y2 = J22D*0.5_rKind                                                   
Y3 = J23D*0.5_rKind
Z1 = J31D*0.5_rKind
Z2 = J32D*0.5_rKind
Z3 = J33D*0.5_rKind


IF((TriTest(X1,Y1,Z1).NE.0).AND.(TriTest(X2,Y2,Z2).NE.0).AND.(TriTest(X3,Y3,Z3).NE.0)&
	.AND.(IntTest(X1,Y1,Z1).NE.0).AND.(IntTest(X2,Y2,Z2).NE.0).AND.(IntTest(X3,Y3,Z3).NE.0)) THEN

tMIN=MAX(ABS((J11D-J33D)),ABS((J32D-J21D)),ABS((J23D-J12D)))
tMAX=MIN((J11D+J33D),(J32D+J21D),(J23D+J12D))

term=0.0_rKind

DO t=tMIN,tMAX,2
term2=SixJ(J11D,J32D,J21D,J33D,J31D,t)*SixJ(J12D,J21D,J22D,t,J32D,J23D)*SixJ(J13D,t,J23D,J11D,J33D,J12D)
term=term+term2*(t+1)*((-1)**t)
END DO

NineJ=term

ELSE
NineJ=0.0_rKind
END IF


END FUNCTION NineJ

FUNCTION HamiOneSite_r(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrixReal), INTENT(IN) :: Op
REAL(KIND=rKind) :: HamiOneSite_r(localSize*localSize,localSize*localSize)

HamiOneSite_r=0.5_rKind*(TensorProd(Op%mr,one_op%mr)+TensorProd(one_op%mr,Op%mr))

END FUNCTION HamiOneSite_r

FUNCTION HamiOneSite_c(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrix), INTENT(IN) :: Op
COMPLEX(KIND=rKind) :: HamiOneSite_c(localSize*localSize,localSize*localSize)

HamiOneSite_c=0.5_rKind*(TensorProd(Op%m,one_op%mr)+TensorProd(one_op%mr,Op%m))

END FUNCTION HamiOneSite_c

FUNCTION HamiLeft_r(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrixReal), INTENT(IN) :: Op
REAL(KIND=rKind) :: HamiLeft_r(localSize*localSize,localSize*localSize)

HamiLeft_r=0.5_rKind*TensorProd(Op%mr,one_op%mr)

END FUNCTION HamiLeft_r

FUNCTION HamiLeft_c(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrix), INTENT(IN) :: Op
COMPLEX(KIND=rKind) :: HamiLeft_c(localSize*localSize,localSize*localSize)

HamiLeft_c=0.5_rKind*TensorProd(Op%m,one_op%mr)

END FUNCTION HamiLeft_c

FUNCTION HamiRight_r(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrixReal), INTENT(IN) :: Op
REAL(KIND=rKind) :: HamiRight_r(localSize*localSize,localSize*localSize)

HamiRight_r=0.5_rKind*TensorProd(one_op%mr,Op%mr)

END FUNCTION HamiRight_r

FUNCTION HamiRight_c(Op)
!
!Purpose: Return the contribution of the one-site operator Op to a Hamiltonian in TEBD form
!
IMPLICIT NONE
TYPE(matrix), INTENT(IN) :: Op
COMPLEX(KIND=rKind) :: HamiRight_c(localSize*localSize,localSize*localSize)

HamiRight_c=0.5_rKind*TensorProd(one_op%mr,Op%m)

END FUNCTION HamiRight_c

END MODULE Hamiltonian_tools_module
