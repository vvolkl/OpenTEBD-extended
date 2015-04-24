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
PROGRAM TEBDSolve
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
! helper variables
INTEGER :: ii, jj, kk
!REAL :: HtempIm(4,4)
!REAL :: HtempRe(4,4)
REAL, allocatable :: HtempIm(:,:)
REAL, allocatable :: HtempRe(:,:)
!Read in input parameters
NAMELIST /SystemSettings/ systemSize, spin, BoundaryCond, trotterOrder
NAMELIST /RTPParams/ chiMax, dtrtpin,totaltime, stepsForStore

OPEN(138,file='TEBDSolve.nml')
READ(138,SystemSettings)
READ(138,RTPParams)
CLOSE(138)

!Print output to screen
print_switch=.TRUE.

spin=spin*1.0_rKind
spinSize=FLOOR(2.0_rKind*spin+1.0_rKind)
maxFilling=spinSize
localSize=heisenbergLocalDim()

!Define rtp parameters
totaltime=totaltime*1.0_rKind
dtRTP=CMPLX(dtrtpin)
totalStep=FLOOR(totaltime/REAL(dtRTP))


IF(print_Switch) THEN
PRINT *, 'Beginning TEBD simulation.'
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

allocate(HtempIm( localSize*localSize,  localSize*localSize))
allocate(HtempRe( localSize*localSize,  localSize*localSize))

CALL CreateHeisenbergOps()
!Create the Hamiltonian
!CALL HamiltonianHeisenberg(H , 1.0_rKind, 1.0_rKind, 0.0_rKind)


open(unit=11, file='Liouvillean_LEFT_IMAG.in')
read(11,*) HtempRe
close(11)
open(unit=11, file='Liouvillean_LEFT_REAL.in')
read(11,*) HtempIm
close(11)
H(1)%m = HtempRe + COMPLEX(0,1) * HtempIm

do ii=2,systemSize-2
open(unit=11, file='Liouvillean_BULK_IMAG.in')
read(11,*) HtempRe
close(11)
open(unit=11, file='Liouvillean_BULK_REAL.in')
read(11,*) HtempIm
close(11)
H(ii)%m = HtempRe + COMPLEX(0,1) * HtempIm
end do

open(unit=11, file='Liouvillean_RIGHT_IMAG.in')
read(11,*) HtempRe
close(11)
open(unit=11, file='Liouvillean_RIGHT_REAL.in')
read(11,*) HtempIm
close(11)
H(systemSize-1)%m = HtempRe + COMPLEX(0,1) * HtempIm

!write(*,*) REAL(H(1)%m)
!write(*,*) AIMAG(H(1)%m)
write(*,*) 'dumping H matrix ...'
!DO i=1,systemSize,1
!end do
open(unit=11, file='hamtest.dat')

kk =1 
write(*,*) 'dumping H matrix ... REAL'
do ii=1,localSize*localSize
write(*,'(16(F10.5))') (REAL(H(kk)%m(jj,ii)), jj=1,localSize*localSize)
end do

write(*,*) 'dumping H matrix ... IMAG'
do ii=1,localSize*localSize
write(*,'(16(F10.5))') (AIMAG(H(kk)%m(jj,ii)), jj=1,localSize*localSize)
end do
close(unit=11)


write(*,*) 'H matrix dumped ...'
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
write(*, *) 'dumping Gammas'
write(*,*) (shape(Gammas(ii)%t), ii=1,systemSize)


write(*, *) 'dumping Lambdas'
write(*,*)  (Lambdas(floor( systemSize * 0.5))%v(kk), kk=1,chiMax)
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

OPEN(unit=12, file='schmidttest.dat', status='replace')
close(12)
OPEN(unit=13, file='gammatest.dat', status='replace')
close(13)
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
!CALL TotalEnergy(energy,H, Gammas, Lambdas)
!CALL OneSiteExpVal(mag,Sz_opS, Gammas, Lambdas)

OPEN(unit=12, file='schmidttest.dat',status='old', Access='append', form='unformatted')
WRITE(12) (Lambdas(ii)%v, ii = 1, systemSize)
CLOSE(12)
OPEN(unit=13, file='gammatest.dat',status='old', Access='append', form='unformatted')
WRITE(13) (Gammas(ii)%t , ii = 1,systemSize)
close(13)
!WRITE(100,*), time, MeyerQMeasure(Gammas,Lambdas), totalTruncerr, (mag(j),j=FLOOR(0.5_rKind*systemSize),systemSize)
!CLOSE(100)

IF(print_switch) THEN
PRINT *, 'RTP step',i,' time =',time
PRINT *, ' truncation error this step is', localTruncerr, 'cumulative truncation error is', totalTruncerr
PRINT *,  (Lambdas(floor(systemSize * 0.5))%v(kk), kk=1,10)

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

END PROGRAM TEBDSolve
