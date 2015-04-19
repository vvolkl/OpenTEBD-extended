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
PROGRAM XXZ_ITP_block
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
TYPE(tensor), POINTER :: Gammas(:) !List of Gamma tensors
TYPE(vector), POINTER :: Lambdas(:) !List of Lambda vectors
!TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:) !Lists of number conserving vectors
TYPE(matrix), POINTER :: H(:) !Hamiltonian
TYPE(measure) :: Measures !Measures derived type
REAL(KIND=rKIND) :: tick, tock !Timing Variables
REAL(KIND=rKind) :: energy, number, dtrtpin !Observables
CHARACTER(len=132) :: localName, avgName, CorrName, entName, stub !File names
REAL(KIND=RKIND), ALLOCATABLE :: mag(:)
INTEGER :: i,j,k !Dummy integers
REAL(KIND=rKIND) :: JL, JR, Delta
character ( len = 80 ) arg
!integer ( kind = 4 ) i
integer ( kind = 4 ) iargc
integer ( kind = 4 ) numarg
INTEGER :: jblock, kblock, j0
REAL(KIND=rKind) :: swaperr
!Read in input parameters
NAMELIST /SystemSettings/ systemSize, spin,  BoundaryCond, trotterOrder
NAMELIST /XXZParams/ JR, JL, Delta
NAMELIST /ITPParams/ chiMin, chiMax, convCriterion1, convCriterion2, stepsForJudge, dtITP, maxITPsteps, itpDir

OPEN(138,file='myparams.nml')
!REWIND(138)
READ(138, SystemSettings)
READ(138,XXZParams)
READ(138, ITPParams)
CLOSE(138)

!SystemSettings
!systemSize=25 !preferably odd  
!spin=0.5
!BoundaryCond='O'
!trotterOrder=2


!XXZParams
!Jz1=0.5
!Jz2=0.5

!i=1
!call getarg ( i, arg )
!Jz1=real(arg)
!write ( *, '(2x,i3,2x,a20)' ) i, Jz1
!i=2
!call getarg ( i, arg )
!Jz2=real(arg)
!write ( *, '(2x,i3,2x,a20)' ) i, Jz2

!ITPParams

!chiMin=15 
!chiMax=20
!convCriterion1=0.00001
!convCriterion2=0.0000001
!stepsForJudge=100
!dtITP=0.05
!maxITPsteps=4000
!itpDir='ITPDATA/'

!Print output to screen
print_switch=.TRUE.
!We're dealing with spins so no number conversation
ncSwitch=.FALSE.
!JL=0.5
!JR=0.5
!Delta=0.8
!Convert to rKind reals
JL=JL*1.0_rKind
JR=JR*1.0_rKind
Delta=Delta*1.0_rKind
spin=spin*1.0_rKind
spinSize=FLOOR(2.0_rKind*spin+1.0_rKind)
maxFilling=spinSize
!Define local dimension
localSize=heisenbergLocalDim()
ALLOCATE(mag(systemSize))
!Define rtp parameters
!totaltime=totaltime*1.0_rKind
!dtRTP=CMPLX(dtrtpin)
!totalStep=FLOOR(totaltime/REAL(dtRTP))

IF(print_Switch) THEN
PRINT *, 'Beginning search for ground state of Jprime XXZ chain.'
PRINT *, systemSize,'sites',spin,'spin'
PRINT *, 'XXZ params: JL', JL, 'JR', JR, 'DELTA', Delta
PRINT *, 'ITP parameters: chi1,2', chiMin, chiMax
PRINT *, 'convergence criteria', convCriterion1, convCriterion2
PRINT *, 'dt', dtITP, 'maxITPsteps', maxITPsteps
PRINT *, 'Order of trotter expansion', trotterOrder
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

CALL CreateHeisenbergOps()
!Create the Hamiltonian
!CALL HamiltonianHeisenberg(H , 1.0_rKind, 1.0_rKind, 0.0_rKind)
CALL perturbedChain(H , JL, JR, Delta)

!Allocate Gammas and Lambdas
CALL AllocateGamLam(Gammas, Lambdas, chiMax)

!Define initial state and propagate in imaginary time to find the ground state
!IF(ncswitch) THEN
!	CALL AllocateLabel(LabelLeft, LabelRight, chiMin)
!	!Define the initial state consistent with number conservation
!	CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
!	!Propagate in imaginary time
!	CALL ImagTimePropNC(H, Gammas, Lambdas, LabelLeft, LabelRight, chiMin)
!	!Number non-conserving routines
!ELSE
!Define the initial state
CALL AllStates(Gammas,Lambdas)

!Propagate in imaginary time
PRINT *, 'Propagating in imaginary time ...'
CALL ImagTimePropSpin(H, Gammas, Lambdas, chiMin)

PRINT *, 'Imag Time Propag completed.'
!Create file names
!CALL SetupBHName(stub,itpDir)
!CALL copyName(stub,localName)
!CALL copyName(stub,avgName)
!CALL copyName(stub,CorrName)
!CALL copyName(stub,entName)
!CALL appendBaseName(localName,'localmeasures.dat')
!CALL appendBaseName(avgName,'averagemeasures.dat')
!CALL appendBaseName(CorrName,'correlationfunctions.dat')
!CALL appendBaseName(entName,'entropies.dat')
CALL createFileName(avgname,itpDir)
CALL appendBaseName(avgname,'Jprime')
!CALL appendBaseName(avgname,'J1',Jz1)
!CALL appendBaseName(avgname,'J2',Jz2)
CALL appendBaseName(avgname,'L',systemSize)
CALL appendBaseName(avgname,'Chi',chiMax)
avgname=TRIM(avgname)//TRIM(ADJUSTL(BoundaryCond))//'BC.dat'
!CALL appendBaseName(avgname,'idt',3,dtITP)
!CALL appendBaseName(avgname,'.dat')
CALL openUnit(avgname,170)
CALL createFileName(entname,itpDir)
CALL appendBaseName(entname,'Entropy')
!CALL appendBaseName(avgname,'JL',JL)
!CALL appendBaseName(avgname,'JR',JR)
CALL appendBaseName(entname,'L',systemSize)
CALL appendBaseName(entname,'Chi',chiMax)
entname=TRIM(entname)//TRIM(ADJUSTL(BoundaryCond))//'BC.dat'
!CALL appendBaseName(avgname,'idt',3,dtITP)
!CALL appendBaseName(avgname,'.dat')
CALL openUnit(entname,171)
PRINT *, 'file created'

CALL TotalEnergy(energy,H, Gammas, Lambdas)

PRINT *, 'energy computed'
CALL OneSiteExpVal(mag,Sz_opS, Gammas, Lambdas)

PRINT *, 'magnetization computed'
!Open files
!CALL openunit(avgName,170)
!CALL openunit(entName,171)

!Record measures
PRINT*, 'E0(J1J2):', energy
DO j=FLOOR(0.5*systemSize)+1,systemSize
   PRINT*, j,mag(j), ChainEntropy(j,Lambdas)
   WRITE(170,*), j,mag(j),ChainEntropy(j,Lambdas)
END DO

!Swap block j...j+k-1 to 1...k
j0 = FLOOR(0.5*systemSize)+1
kblock = 2
CALL SVDInit(chiMax)
PRINT *, 'start the swap of the block ...'
DO jblock=j0,j0+kblock-1
   DO k=0,j0-2
      PRINT*, 'swapping ', jblock-k
      CALL Swapping(jblock-k,swaperr,Gammas,Lambdas)
      PRINT*, 'swap (Schmmidt) error = ', swaperr
   END DO
END DO
PRINT *, 'swapping completed.'
CALL SVDFinish()

PRINT *, 'Entropy of the right sub-chain:', log(2.)*ChainEntropy(FLOOR(0.5*systemSize)+1,Lambdas)
WRITE (171,*) JR, JL, Delta, systemSize, log(2.)*ChainEntropy(FLOOR(0.5*systemSize)+1,Lambdas)
!CALL RecordOneSiteOb(170, REAL(Measures%local(1)%value))
!CALL RecordOneSiteOb(170, REAL(Measures%local(2)%value))
!WRITE(171, '(2E30.15)') Real(Measures%avg(1)%value), Measures%en
!CALL RecordTwoSiteOb(172, Measures%corr(1)%value)
!CALL RecordTwoSiteOb(172, Measures%corr(2)%value)
!WRITE(173, '(1E30.15)') Measures%ent%qme
!CALL RecordOneSiteOb(173, Measures%ent%vn)
!CALL RecordOneSiteOb(173, Measures%ent%chain(1:systemSize))
!CALL RecordTwoSiteOb(173, Measures%ent%tbvn)

!Close files
CLOSE(170)
CLOSE(171)
!CLOSE(172)
!CLOSE(173)

!Clean up
CALL DeallocateGamLam(Gammas, Lambdas)

CALL DestroyHeisenbergOps()
IF(BoundaryCond=='O') THEN
	CALL DeallocateOps(H,systemSize-1)
ELSE
	CALL DeallocateOps(H,systemSize)
END IF
!CALL DeallocateProp(Urtp)
DEALLOCATE(mag)



!End the timing routine
CALL CPU_TIME(tock)
PRINT *, 'XXZ ITP Case study exited normally!'
PRINT *, 'Time taken=',tock-tick,'seconds!'

END PROGRAM XXZ_ITP_block

