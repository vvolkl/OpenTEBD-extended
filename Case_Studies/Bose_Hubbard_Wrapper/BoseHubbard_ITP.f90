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
PROGRAM BoseHubbard_ITP
!
! Purpose: Main program to compute the ground state properties of the Bose
! Hubbard model for OpenSourceTEBD v2.0 Case Study.
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
TYPE(tensor), POINTER :: Gammas(:) !List of gamma tensors
TYPE(vector), POINTER :: Lambdas(:) !List of Lambda vectors
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:) !Lists of number conserving vectors
TYPE(matrix), POINTER :: H(:) !Hamiltonian
REAL(KIND=rKind) ::  trapOmega,swaperr !Amplitude of optional external trap
TYPE(measure) :: Measures !Measures derived type
REAL(KIND=rKIND) :: tick, tock !Timing Variables
CHARACTER(len=132) :: localName, avgName, CorrName, entName, stub !File names
INTEGER :: i,j,k !Dummy integers
INTEGER :: j0,jblock,kblock
!Read in input parameters
NAMELIST /SystemSettings/ systemSize, maxFilling, totNum,  BoundaryCond, trotterOrder
NAMELIST /BHParams/ jTunn, U0, V0, mu0
NAMELIST /ITPParams/ chiMin, chiMax, convCriterion1, convCriterion2, stepsForJudge, dtITP, maxITPsteps, itpDir

OPEN(138,file='BH_ITP.nml')
READ(138,SystemSettings)
READ(138,BHParams)
READ(138,ITPParams)
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

!Define local dimension
localSize=boseHubbardLocalDim()

!!If desired, put on  a trap
!trapOmega=0.1*jTunn
!ALLOCATE(extPot(systemSize))
!DO i=1,systemSize,1
!extPot(i)=trapOmega*(i*1.0_rKind-0.5_rKind*(systemSize+1))**2
!END DO

IF(print_Switch) THEN
PRINT *, 'Beginning Bose-Hubbard ITP Case study.'
PRINT *, systemSize,'sites',totNum,'particles',maxFilling,'particles allowed per site'
PRINT *, 'Hubbard params: J',jTunn,'U',U0,'V',V0,'mu',mu0
PRINT *, 'ITP parameters: chi1,2',chiMin,chiMax
PRINT *, 'convergence criteria',convCriterion1,convCriterion2
PRINT *, 'dt',dtITP,'maxITPsteps',maxITPsteps
PRINT *, 'Order of trotter expansion',trotterOrder
IF(BoundaryCond=='O') THEN
	PRINT *, 'Open boundary conditions'
ELSE IF(BoundaryCond=='P') THEN
	PRINT *, 'Periodic boundary conditions'
END IF
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

!Allocate necessary measures
	CALL AllocateMeasures(Measures,2, 1, 2, 0, 2)
!Define local measures
!Number
Measures%local(1)%op=MATMUL(TRANSPOSE(a_op%mr),a_op%mr)
!Number squared
Measures%local(2)%op=MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Define average measures
!Number
Measures%avg(1)%op=MATMUL(TRANSPOSE(a_op%mr),a_op%mr)

!Define correlation functions
!single particle density matrix
Measures%corr(1)%op=TensorProd(TRANSPOSE(a_op%mr),a_op%mr)
!<nn>
Measures%corr(2)%op=TensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Evaluate the given measures
	CALL EvaluateMeasures(Measures, Gammas, Lambdas, H)

DO i=1,systemSize
	!Form fluctuations from <n^2> and <n>^2
	Measures%local(2)%value(i)=SQRT(Measures%local(2)%value(i)-Measures%local(1)%value(i)**2)
	!Add diagonal elements to spdm
	Measures%corr(1)%value(i,i)=Measures%local(1)%value(i)
	Measures%corr(2)%value(i,i)=0.0_rKind
END DO

!Create file names
CALL SetupBHName(stub,itpDir)
CALL copyName(stub,localName)
CALL copyName(stub,avgName)
CALL copyName(stub,CorrName)
CALL copyName(stub,entName)
CALL appendBaseName(localName,'localmeasures.dat')
CALL appendBaseName(avgName,'averagemeasures.dat')
CALL appendBaseName(CorrName,'correlationfunctions.dat')
CALL appendBaseName(entName,'entropies.dat')

!Open files
CALL openunit(localName,170)
CALL openunit(avgName,171)
CALL openunit(CorrName,172)
CALL openunit(entName,173)

!Record measures
CALL RecordOneSiteOb(170, REAL(Measures%local(1)%value))
CALL RecordOneSiteOb(170, REAL(Measures%local(2)%value))
WRITE(171, '(2E30.15)') Real(Measures%avg(1)%value), Measures%en
CALL RecordTwoSiteOb(172, Measures%corr(1)%value)
CALL RecordTwoSiteOb(172, Measures%corr(2)%value)
WRITE(173, '(1E30.15)') Measures%ent%qme
CALL RecordOneSiteOb(173, Measures%ent%vn)
CALL RecordOneSiteOb(173, Measures%ent%chain(1:systemSize))
CALL RecordTwoSiteOb(173, Measures%ent%tbvn)

!Close files
CLOSE(170)
CLOSE(171)
CLOSE(172)
CLOSE(173)

!Clean up
IF(ncswitch) THEN
CALL DeallocateLabel(LabelLeft, LabelRight)
END IF

CALL DeallocateGamLam(Gammas, Lambdas)
CALL DeallocateMeasures(Measures)
CALL DestroyFieldOps()
IF(BoundaryCond=='O') THEN
	CALL DeallocateOps(H,systemSize-1)
ELSE
	CALL DeallocateOps(H,systemSize)
END IF

IF(ALLOCATED(extPot)) THEN
DEALLOCATE(extPot)
END IF

!End the timing routine
CALL CPU_TIME(tock)
PRINT *, 'Bose-Hubbard ITP Case study exited normally!'
PRINT *, 'Time taken=',tock-tick,'seconds!'

END PROGRAM BoseHubbard_ITP
