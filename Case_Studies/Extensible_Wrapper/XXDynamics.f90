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
PROGRAM XXDynamics
!
! Purpose: Main program to compute the dynamics of the Bose
! Hubbard model during a parameter quench for OpenSourceTEBD v1.0 Case Study.
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   2/24/09   M. L. Wall	v1.0 release
!

USE system_parameters
USE TEBDtools_module
USE io_module
USE Hamiltonian_tools_module
USE bose_hubbard_module
USE fermi_hubbard_module
USE heisenberg_module
USE spinS_module
USE rotation_module
USE local_operations_module
USE observables_module
USE propagation_module


IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:),Gammas0(:) !List of gamma tensors
TYPE(vector), POINTER :: Lambdas(:),Lambdas0(:) !List of Lambda vectors
TYPE(matrix), POINTER :: H(:) !Hamiltonian
TYPE(matrix), POINTER :: Urtp(:) !Real time propagator
REAL(KIND=rKind) ::  time, totaltime, totalTruncerr, localTruncerr !Time, total time, and truncation errors
REAL(KIND=rKind) :: energy, number, dtrtpin !Observables
REAL(KIND=rKIND) :: tick, tock !Timing Variables
COMPLEX(KIND=rKIND), ALLOCATABLE :: coefArray(:,:)
REAL(KIND=RKIND), ALLOCATABLE :: mag(:)
CHARACTER(len=132) :: localName, avgName, CorrName, entName, stub !File names
INTEGER :: i,j,k !Dummy integers
!Read in input parameters
NAMELIST /SystemSettings/ systemSize,spin,  BoundaryCond, trotterOrder
NAMELIST /RTPParams/ chiMax,dtrtpin,totaltime,  stepsForStore

OPEN(138,file='XXDynamics.nml')
READ(138,SystemSettings)
READ(138,RTPParams)
CLOSE(138)

!Print output to screen
print_switch=.TRUE.
ncSwitch=.FALSE.

spin=spin*1.0_rKind
spinSize=FLOOR(2.0_rKind*spin+1.0_rKind)
maxFilling=spinSize
localSize=heisenbergLocalDim()

!Define rtp parameters
totaltime=totaltime*1.0_rKind
dtRTP=CMPLX(dtrtpin)
totalStep=FLOOR(totaltime/REAL(dtRTP))


IF(print_Switch) THEN
PRINT *, 'Beginning XX dynamics Case study.'
PRINT *, systemSize,'sites',spin,'spin'
PRINT *, 'Order of trotter expansion',trotterOrder
IF(BoundaryCond=='O') THEN
	PRINT *, 'Open boundary conditions'
ELSE IF(BoundaryCond=='P') THEN
	PRINT *, 'Periodic boundary conditions'
END IF
PRINT *, 'dt for RTP', REAL(dtRTP), 'Totaltime', totaltime
END IF
!Begin timing
CALL CPU_TIME(tick)


!Allocate Hamiltonian
IF(BoundaryCond=='O') THEN
	CALL AllocateOps(H,systemSize-1,localSize*localSize)
ELSE
	CALL AllocateOps(H,systemSize,localSize*localSize)
END IF

CALL CreateHeisenbergOps()
!Create the Hamiltonian
CALL HamiltonianHeisenberg(H , 1.0_rKind, 1.0_rKind, 0.0_rKind)

write(*,*) SHAPE(H)
write(*,*) 'dumping H matrix ...'
!DO i=1,systemSize,1
write(*,*) REAL(H(2)%m)
!end do

!Allocate the Gammas, Labdas, and labels
CALL AllocateGamLam(Gammas, Lambdas, chiMax)

!Allocate a matrix to imprint the initial state
ALLOCATE(coefArray(localSize,systemSize))
ALLOCATE(mag(systemSize))

!Define the step state consisting of spins aligned with local magnetization 0.5
!to the left of the center and spins aligned with local magnetization -0.5
!to the right of the center 
DO i=1,systemSize,1
	coefArray(:,i)=0.0_rKind
	IF(i.le.FLOOR(0.5_rKind*systemSize)) THEN
	coefArray(1,i)=1.0_rKind
	ELSE
	coefArray(2,i)=1.0_rKind
	END IF
END DO

!write(*,*) coefArray(1,20)
!write(*,*) coefArray(2,20)
!Imprint the state on Lambdas & Gammas
CALL ProductStateMPD(Gammas, Lambdas, coefArray)

!DO i=1,systemSize,1
!write(*,*) Gammas(i)%t(1,1,1)
!write(*,*) Gammas(i)%t(1,2,1)
!end do
write(*,*) Gammas(49)%t(1,1,1)
write(*,*) Gammas(49)%t(1,2,1)
write(*,*) Gammas(50)%t(1,1,1)
write(*,*) Gammas(50)%t(1,2,1)
write(*,*) Gammas(51)%t(1,1,1)
write(*,*) Gammas(51)%t(1,2,1)
!Initialize time and cumulative truncation error
time=0.0_rKind
totalTruncerr=0.0_rKind

!Allocate and construct real time propagator
CALL AllocateProp(Urtp)
CALL ConstructPropagators(H, Urtp, dtrtp)

!Set up i/o
CALL createFileName(avgname,rtpDir)
CALL appendBaseName(avgname,'XXDyn_')
CALL appendBaseName(avgname,'L',systemSize)
CALL appendBaseName(avgname,'Chi',chiMax)
avgname=TRIM(avgname)//TRIM(ADJUSTL(BoundaryCond))//'BC'
CALL appendBaseName(avgname,'dt',3,dtrtpin)
CALL appendBaseName(avgname,'.dat')
CALL openUnit(avgname,100)

CALL TotalEnergy(energy,H, Gammas, Lambdas)
CALL OneSiteExpVal(mag,Sz_opS, Gammas, Lambdas)

WRITE(100,*), time, MeyerQMeasure(Gammas,Lambdas), totalTruncerr, (mag(j),j=FLOOR(0.5_rKind*systemSize),systemSize)
CLOSE(100)

! Initialize variables for SVD routine.
	CALL SVDInit(chiMax)

!Time propagation loop
DO i=1,totalStep

CALL TrotterStep(Urtp, Gammas, Lambdas, localTruncerr)
totalTruncerr=totalTruncerr+localTruncerr

!Increase time step
time=time+dtRTP
!Write out data	
IF(MOD(i,stepsForStore)==0) THEN

CALL openUnit(avgName,100,'A')

CALL TotalEnergy(energy,H, Gammas, Lambdas)
CALL OneSiteExpVal(mag,Sz_opS, Gammas, Lambdas)

WRITE(100,*), time, MeyerQMeasure(Gammas,Lambdas), totalTruncerr, (mag(j),j=FLOOR(0.5_rKind*systemSize),systemSize)
CLOSE(100)

IF(print_switch) THEN
PRINT *, 'RTP step',i,' time =',time
PRINT *, ' truncation error this step is', localTruncerr, 'cumulative truncation error is', totalTruncerr
END IF
		END IF
		
END DO

! Deallocate SVD variables.
CALL SVDFinish()

!Clean up
CALL DeallocateGamLam(Gammas, Lambdas)
CALL DestroyHeisenbergOps()
IF(BoundaryCond=='O') THEN
	CALL DeallocateOps(H,systemSize-1)
ELSE
	CALL DeallocateOps(H,systemSize)
END IF
CALL DeallocateProp(Urtp)
DEALLOCATE(coefArray)
DEALLOCATE(mag)

!End the timing routine
CALL CPU_TIME(tock)
PRINT *, 'XX dynamics Case study exited normally!'
PRINT *, 'Time taken=',tock-tick,'seconds!'

END PROGRAM XXDynamics
