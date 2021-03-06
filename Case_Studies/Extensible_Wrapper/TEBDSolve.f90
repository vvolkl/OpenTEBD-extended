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
! Purpose: Main program to compute the dynamics of a general Model with an
! Evolution Operator contained in the Propagator*.dat files.
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   2/24/09   M. L. Wall	v1.0 release
!
! Modified 2015  by V. Volkl 

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
REAL(KIND=rKIND), ALLOCATABLE :: coefArrayRe(:,:)
REAL(KIND=RKIND), ALLOCATABLE :: mag(:)
CHARACTER(len=132) :: localName, avgName, CorrName, entName, stub !File names
INTEGER :: i,j,k !Dummy integers
! helper variables
INTEGER :: ii, jj, kk
! matrices to hold real and imaginary part of the evolution operator read from
! file.
REAL, allocatable :: UtempIm(:,:)
REAL, allocatable :: UtempRe(:,:)
real :: tr, Mz



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


!CALL CreateHeisenbergOps()
!Create the Hamiltonian
!CALL HamiltonianHeisenberg(H , 1.0_rKind, 1.0_rKind, 0.0_rKind)

!Allocate the Gammas, Labdas, and labels
CALL AllocateGamLam(Gammas, Lambdas, chiMax)

!Allocate a matrix to imprint the initial state
ALLOCATE(coefArray(localSize,systemSize))
ALLOCATE(coefArrayRe(localSize,systemSize))
ALLOCATE(mag(systemSize))

!Define the step state consisting of spins aligned with local magnetization 0.5
!to the left of the center and spins aligned with local magnetization -0.5
!to the right of the center 
DO i=1,systemSize,1
	coefArray(:,i)=0.0_rKind
	!IF(i.le.FLOOR(0.5_rKind*systemSize)) THEN
	coefArray(1,i)=1.0_rKind
	!ELSE
	!coefArray(2,i)=1.0_rKind
	!END IF
END DO
!OPEN(unit=12, file='initialStatCoeffTest.dat')
!do ii=1,localsize
!READ(12, * ) coefArrayRe(ii, :)
!end do
!CLOSE(12)

!coefArray = coefArray + coefArrayRe
!write(*,*) 'coefArray as read from file'
write(*,*) Real(coefArray(1,:))

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
!CALL ConstructPropagators(H, Urtp, dtrtp)


allocate(UtempIm( localSize*localSize,  localSize*localSize))
allocate(UtempRe( localSize*localSize,  localSize*localSize))
open(unit=11, file='Propagator_0.005_left_IMAG.nml')
read(11,*) UtempIm
close(11)
open(unit=11, file='Propagator_0.005_left_REAL.nml')
read(11,*) UtempRe
close(11)
Urtp(1)%m =  UtempRe + COMPLEX(0,1) * UtempIm
do ii=3,systemSize-2,2
open(unit=11, file='Propagator_0.005_bulk_IMAG.nml')
read(11,*) UtempIm
close(11)
open(unit=11, file='Propagator_0.005_bulk_REAL.nml')
read(11,*) UtempRe
close(11)
Urtp(ii)%m =  UtempRe + COMPLEX(0,1) * UtempIm
end do
do ii=2,systemSize-2,2
open(unit=11, file='Propagator_0.01_bulk_IMAG.nml')
read(11,*) UtempIm
close(11)
open(unit=11, file='Propagator_0.01_bulk_REAL.nml')
read(11,*) UtempRe
close(11)
Urtp(ii)%m =  UtempRe + COMPLEX(0,1) * UtempIm
end do

open(unit=11, file='Propagator_0.005_right_IMAG.nml')
read(11,*) UtempIm
close(11)
open(unit=11, file='Propagator_0.005_right_REAL.nml')
read(11,*) UtempRe
close(11)
Urtp(systemSize-1)%m = UtempRe + COMPLEX(0,1) * UtempIm

!open(unit=11, file='hamtestreal.dat')
!write(*,*) 'dumping U matrix ...'
!kk =1 
!write(*,*) 'dumping U matrix ... REAL'
!do ii=1,localSize*localSize
!write(11,'(16(F10.5))') (REAL(Urtp(kk)%m(jj,ii)), jj=1,localSize*localSize)
!end do
!close(unit=11)
!open(unit=11, file='hamtestimag.dat')
!write(*,*) 'dumping U matrix ... IMAG'
!do ii=1,localSize*localSize
!write(11,'(16(F10.5))') (AIMAG(Urtp(kk)%m(jj,ii)), jj=1,localSize*localSize)
!end do
!close(unit=11)
!write(*,*) 'U matrix dumped ...'

!Set up i/o
CALL createFileName(avgname,rtpDir)
CALL appendBaseName(avgname,'Dyn_')
CALL appendBaseName(avgname,'L',systemSize)
CALL appendBaseName(avgname,'Chi',chiMax)
avgname=TRIM(avgname)//TRIM(ADJUSTL(BoundaryCond))//'BC'
CALL appendBaseName(avgname,'dt',3,dtrtpin)
CALL appendBaseName(avgname,'.dat')
CALL openUnit(avgname,100)

!CALL TotalEnergy(energy,H, Gammas, Lambdas)
!CALL OneSiteExpVal(mag,Sz_opS, Gammas, Lambdas)

WRITE(100,*), time, MeyerQMeasure(Gammas,Lambdas), totalTruncerr, (mag(j),j=FLOOR(0.5_rKind*systemSize),systemSize)
CLOSE(100)

OPEN(unit=12, file='RTPDATA/schmidttest.dat', status='replace')
close(12)
OPEN(unit=13, file='RTPDATA/Gammatest.dat', status='replace')
close(13)
OPEN(unit=14, file='RTPDATA/Mztest.dat', status='replace')
close(14)
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

OPEN(unit=12, file='RTPDATA/schmidttest.dat',status='old', Access='append')
!WRITE(12, *) (Lambdas(ii)%v, ii = 1, systemSize)
WRITE(12, *) (Lambdas(systemSize / 2)%v)
CLOSE(12)
!OPEN(unit=13, file='gammatest.dat',status='old', Access='append', form='unformatted')
!WRITE(13) (Gammas(ii)%t , ii = 1,systemSize)
!close(13)
!WRITE(100,*), time, MeyerQMeasure(Gammas,Lambdas), totalTruncerr, (mag(j),j=FLOOR(0.5_rKind*systemSize),systemSize)
!CLOSE(100)

IF(print_switch) THEN
PRINT *, 'RTP step',i,' time =',time
PRINT *, ' truncation error this step is', localTruncerr, 'cumulative truncation error is', totalTruncerr
PRINT *,  (Lambdas(floor(systemSize * 0.5))%v(kk), kk=1,10)

PRINT *, 'test measurement routine'
PRINT *, shape(Gammas(1)%t)
call  mpstrace(Gammas, Lambdas, systemSize, localSize, chiMax,chiMax, tr)
call  measureMz(Gammas, Lambdas, systemSize, localSize, chiMax,chiMax, Mz)
PRINT *, tr, Mz,  Mz / tr
OPEN(unit=14, file='RTPDATA/Mztest.dat',status='old', Access='append')
WRITE(14, *)  tr, Mz
CLOSE(14)

END IF
		END IF
		
END DO

! Deallocate SVD variables.
CALL SVDFinish()

!Clean up
!CALL DeallocateGamLam(Gammas, Lambdas)
!CALL DestroyHeisenbergOps()
!IF(BoundaryCond=='O') THEN
!	CALL DeallocateOps(H,systemSize-1)
!ELSE
!	CALL DeallocateOps(H,systemSize)
!END IF
CALL DeallocateProp(Urtp)
DEALLOCATE(coefArray)
DEALLOCATE(mag)

!End the timing routine
CALL CPU_TIME(tock)
PRINT *, 'Case study exited normally!'
PRINT *, 'Time taken=',tock-tick,'seconds!'

contains
subroutine mpstrace(Gammas, Lambdas, length, d2, chi1, chi2, tr)

    implicit none
    integer, intent(in)                         ::  length, d2, chi1, chi2
    TYPE(tensor), POINTER, intent(in) :: Gammas(:)
    TYPE(vector), POINTER, intent(in) :: Lambdas(:)
    complex    ::  mtx3(chi1, chi2)
    integer :: i,j, k
    real, intent(out) :: tr
    tr = 0
    do i=1,chi1
    do j=1,chi2
        mtx3(i,j)=0
        if (i == j) then
        mtx3(i,j) = 1
        end if
    end do
    end do

    !do i=1,chi1
    !write(*,'(5F5.2)') (mtx3(i,j), j=1,chi2)
    !end do
    !if(size(mtx1,dim=2) /= size(mtx2,dim=1)) stop "input array sizenot match"
    do i=1,length
    do j=1,chi1
    mtx3(:,j) = mtx3(:,j) * Lambdas(i)%v(j)
    end do
    !write(*,*) 'site ', i
    mtx3 =   matmul(mtx3,Gammas(i)%t(:,1,:) + Gammas(i)%t(:,d2,:))
    !write(*,*) 'multiplication done!'
    !do j=1,m3
    !write(*,'(5F5.2)') (mtx3(j,k), k=1,m4)
    !end do
    end do
do i=1,chi1
tr = tr + mtx3(i,i)
end do

    !do i=1,m3
    !write(*,'(5F5.2)') (mtx3(i,j), j=1,m4)
    !end do
!write(*,*) tr
end subroutine
subroutine measureMz(Gammas, Lambdas, length, d2, chi1, chi2, tr)
    implicit none
    integer, intent(in)                         ::  length, d2, chi1, chi2
    TYPE(tensor), POINTER, intent(in) :: Gammas(:)
    TYPE(vector), POINTER, intent(in) :: Lambdas(:)
    complex    ::  mtx3(chi1, chi2)
    integer :: i,j, k
    real, intent(out) :: tr
    tr = 0
    do i=1,chi1
    do j=1,chi2
        mtx3(i,j)=0
        if (i == j) then
        mtx3(i,j) = 1
        end if
    end do
    end do

    !do i=1,chi1
    !write(*,'(5F5.2)') (mtx3(i,j), j=1,chi2)
    !end do
    !if(size(mtx1,dim=2) /= size(mtx2,dim=1)) stop "input array sizenot match"
    do i=1,length
    !write(*,*) 'site ', i
    do j=1,chi1
    mtx3(:,j) = mtx3(:,j) * Lambdas(i)%v(j)
    end do
    if (i == 3 ) then
    mtx3 =   matmul(mtx3,Gammas(i)%t(:,1,:) - Gammas(i)%t(:,d2,:))
    else 
    mtx3 =   matmul(mtx3,Gammas(i)%t(:,1,:) + Gammas(i)%t(:,d2,:))
    end if
    !write(*,*) 'multiplication done!'
    !do j=1,m3
    !write(*,'(5F5.2)') (mtx3(j,k), k=1,m4)
    !end do
    end do
do i=1,chi1
tr = tr + mtx3(i,i)
end do

    !do i=1,m3
    !write(*,'(5F5.2)') (mtx3(i,j), j=1,m4)
    !end do
!write(*,*) tr
end subroutine
END PROGRAM TEBDSolve
