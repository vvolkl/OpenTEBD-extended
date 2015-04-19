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
PROGRAM BoseHubbard_Quench
!
! Purpose: Main program to compute the dynamics of the Bose
! Hubbard model during a parameter quench for OpenSourceTEBD v2.0 Case Study.
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
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:) !Lists of number conserving vectors
TYPE(matrix), POINTER :: H(:) !Hamiltonian
TYPE(matrix), POINTER :: Urtp(:) !Real time propagator
REAL(KIND=rKind) ::  rampTime, time, totalTruncerr, localTruncerr !Length of linear ramp, time, and truncation errors
REAL(KIND=rKind) :: energy, number, dtrtpin !Observables
REAL(KIND=rKIND) :: tick, tock !Timing Variables
CHARACTER(len=132) :: localName, avgName, CorrName, entName, stub !File names
INTEGER :: i,j,k !Dummy integers
!Read in input parameters
NAMELIST /SystemSettings/ systemSize, maxFilling, totNum,BoundaryCond, trotterOrder
NAMELIST /BHParams/ jTunn, U0, V0, mu0
NAMELIST /ITPParams/ chiMin, chiMax, convCriterion1, convCriterion2, stepsForJudge, dtITP, maxITPsteps, itpDir
NAMELIST /RTPParams/ dtrtpin, ramptime, stepsForStore

OPEN(138,file='BH_Quench.nml')
READ(138,SystemSettings)
READ(138,BHParams)
READ(138,ITPParams)
READ(138,RTPParams)
CLOSE(138)

!Print output to screen
print_switch=.TRUE.
!Number conservation
IF(totNum==0) THEN
	ncswitch=.FALSE.
ELSE
	ncswitch=.TRUE.
END IF

!Convert to rKind reals
jTunn=jTunn*1.0_rKind
U0=U0*1.0_rKind
V0=V0*1.0_rKind
mu0=mu0*1.0_rKind
convCriterion1=convCriterion1*1.0_rKind
convCriterion2=convCriterion2*1.0_rKind
dtrtpin=dtrtpin*1.0_rKind
dtRTP=CMPLX(dtrtpin)
ramptime=ramptime*1.0_rKind

!Define local dimension
localSize=boseHubbardLocalDim()

!Define rtp parameters
totalStep=FLOOR(rampTime/REAL(dtRTP))

IF(print_Switch) THEN
PRINT *, 'Beginning Bose-Hubbard Quench Case study.'
PRINT *, systemSize,'sites',totNum,'particles',maxFilling,'particles allowed per site'
PRINT *, 'Hubbard params: J',jTunn,'U',U0,'V',V0,'mu',mu0
PRINT *, 'ITP parameters: chi1,2',chiMin,chiMax
PRINT *, 'convergence criteria',convCriterion1,convCriterion2
PRINT *, 'dt for ITP',dtITP,'maxITPsteps',maxITPsteps
PRINT *, 'Order of trotter expansion',trotterOrder
IF(BoundaryCond=='O') THEN
	PRINT *, 'Open boundary conditions'
ELSE IF(BoundaryCond=='P') THEN
	PRINT *, 'Periodic boundary conditions'
END IF
PRINT *, 'dt for RTP', REAL(dtRTP), 'ramp time', rampTime
END IF
!Begin timing
CALL CPU_TIME(tick)

!Allocate Hamiltonian
IF(BoundaryCond=='O') THEN
	CALL AllocateOps(H,systemSize-1,localSize*localSize)
ELSE
	CALL AllocateOps(H,systemSize,localSize*localSize)
END IF

!Create BH operators and Hamiltonian
CALL CreateFieldOps()

IF(ALLOCATED(extPot)) THEN
CALL HamiltonianBoseHubbard(H , jTunn, U0,mu0, V0, extPot)
ELSE
CALL HamiltonianBoseHubbard(H , jTunn, U0,mu0, V0)
END IF

!Allocate Gammas and Lambdas
CALL AllocateGamLam(Gammas, Lambdas, chiMin)

!Define initial state and propagate in imaginary time to find the ground state
IF(ncswitch) THEN
	CALL AllocateLabel(LabelLeft, LabelRight, chiMin)
	!Define the initial state consistent with number conservation
	CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
	!Propagate in imaginary time
	CALL ImagTimePropNC(H, Gammas, Lambdas, LabelLeft, LabelRight, chiMin)
	!Number non-conserving routines
ELSE
!Define the initial state
CALL AllStates(Gammas,Lambdas)
!Propagate in imaginary time
CALL ImagTimeProp(H, Gammas, Lambdas, chiMin)
END IF

IF(print_switch) THEN
PRINT *, 'ITP finished!'
END IF

!Copy the initial state to a set of dummy Gammas/Labdas
CALL AllocateGamLam(Gammas0, Lambdas0, chiMax)
CALL CopyGamLam(Gammas0, Lambdas0, Gammas, Lambdas)

!Initialize time and cumulative truncation error
time=0.0_rKind
totalTruncerr=0.0_rKind

!Allocate and construct real time propagator
CALL AllocateProp(Urtp)
CALL ConstructPropagators(H, Urtp, dtrtp)

!Set up i/o
CALL createFileName(avgname,rtpDir)
CALL appendBaseName(avgname,'BHQuench_')
CALL appendBaseName(avgname,'tR',1,rampTime)
CALL appendBaseName(avgname,'L',systemSize)
CALL appendBaseName(avgname,'N',totNum)
CALL appendBaseName(avgname,'Chi',chiMax)
CALL appendBaseName(avgname,'.dat')
CALL openUnit(avgname,100)

CALL TotalEnergy(energy,H, Gammas, Lambdas)
CALL TotalNumber(number, Gammas, Lambdas)

WRITE(100,*), time, U0,energy, MeyerQMeasure(Gammas,Lambdas), ABS(InnerProduct(Gammas,Lambdas,Gammas0,Lambdas0))**2, totalTruncerr
CLOSE(100)

IF(.NOT.ncswitch) THEN
! Initialize variables for SVD routine.
	CALL SVDInit(chiMax)
END IF

!Time propagation loop
DO i=1,totalStep

IF(.NOT.ncswitch) THEN
!Trotter step one dt
	CALL TrotterStep(Urtp, Gammas, Lambdas, localTruncerr)
	totalTruncerr=totalTruncerr+localTruncerr
ELSE
	CALL TrotterStepNC(Urtp, Gammas, Lambdas, LabelLeft, LabelRight, localTruncerr)
	totalTruncerr=totalTruncerr+localTruncerr
END IF

!Increase time step
time=time+dtRTP
!Ramp, redefine Hamiltonian and Propagator
IF(time.le.(0.5_rKind*rampTime)) THEN
U0=U0-dtrtp*18.0_rKind/rampTime
CALL HamiltonianBoseHubbard(H , jTunn, U0,mu0, V0)
CALL ConstructPropagators(H, Urtp, dtrtp)
ELSE
U0=U0+dtrtp*18.0_rKind/rampTime
CALL HamiltonianBoseHubbard(H , jTunn, U0,mu0, V0)
CALL ConstructPropagators(H, Urtp, dtrtp)
END IF

!Write out data	
			IF(MOD(i,stepsForStore)==0) THEN

CALL openUnit(avgName,100,'A')

CALL TotalEnergy(energy,H, Gammas, Lambdas)
CALL TotalNumber(number, Gammas, Lambdas)

WRITE(100,*), time, U0,energy, MeyerQMeasure(Gammas,Lambdas), ABS(InnerProduct(Gammas,Lambdas,Gammas0,Lambdas0))**2, totalTruncerr
CLOSE(100)

IF(print_switch) THEN
PRINT *, 'RTP step',i,' time =',time,'number is', number, ' energy is', energy
PRINT *, ' truncation error this step is', localTruncerr, 'cumulative truncation error is', totalTruncerr
END IF
		END IF
		
END DO

IF(.NOT.ncswitch) THEN		
! Deallocate SVD variables.
          CALL SVDFinish()
END IF

!Clean up
CALL DeallocateGamLam(Gammas, Lambdas)
CALL DeallocateGamLam(Gammas0, Lambdas0)
IF(ncswitch) THEN
CALL DeallocateLabel(LabelLeft, LabelRight)
END IF
CALL DestroyFieldOps()
IF(BoundaryCond=='O') THEN
	CALL DeallocateOps(H,systemSize-1)
ELSE
	CALL DeallocateOps(H,systemSize)
END IF
CALL DeallocateProp(Urtp)


!End the timing routine
CALL CPU_TIME(tock)
PRINT *, 'Bose-Hubbard Quench Case study exited normally!'
PRINT *, 'Time taken=',tock-tick,'seconds!'

END PROGRAM BoseHubbard_Quench
