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
MODULE SitesParallelSetup_Module
!
! Purpose: Module to include MPI and define global variables for parallel implementations
!for OpenSourceTEBD v1.0
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
USE heisenberg_module
USE spinS_module
USE rotation_module
USE local_operations_module
USE observables_module
USE propagation_module
USE timing_module

IMPLICIT NONE
! *** MPI Parameters
INCLUDE "mpif.h"
INTEGER :: ierror !Error flag for MPI
INTEGER :: errorcode !Error code for MPI_Abort
INTEGER :: my_rank !processor number
INTEGER :: num_cores !Total number of processors
INTEGER :: my_tag !Send/Recv tag
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: my_status !Status flag for MPI_Recv
INTEGER :: my_request !request tag for nonblocking sends
INTEGER :: my_local_dim !Number of Gammas a processor owns
INTEGER :: my_bounds(2) !Bounds of Gammas the processor owns
INTEGER :: my_mpi_rKind, my_mpi_cKind, my_count
INTEGER, ALLOCATABLE :: countvec(:),offsetVec(:)


INTERFACE OneSiteExpValParallel
MODULE PROCEDURE OneSiteExpValParallel_r,OneSiteExpValParallel_rc,&
				 OneSiteExpValParallel_cr,OneSiteExpValParallel_c
END INTERFACE  OneSiteExpValParallel

INTERFACE OneSiteVarParallel
MODULE PROCEDURE OneSiteVarParallel_r,&
				 OneSiteVarParallel_c
END INTERFACE  OneSiteVarParallel

INTERFACE TwoSiteExpValParallelG
MODULE PROCEDURE TwoSiteExpValParallelG_r,&
				 TwoSiteExpValParallelG_c, &
				 TwoSiteExpValParallelG_rc,&
				 TwoSiteExpValParallelG_cr
END INTERFACE  TwoSiteExpValParallelG


CONTAINS

SUBROUTINE Initialize_MPI()
!
!Purpose: Set up the MPI communicator world and derived types
!
IMPLICIT NONE

CALL MPI_Init(ierror) !initialize MPI
CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,ierror) !Find the rank of each processor
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_cores,ierror) !Find the total number of cores

!Define the proper KIND for REAL rKind MPI messages
CALL MPI_Type_create_f90_real(precis,range,my_mpi_rKind,ierror)
CALL MPI_Type_commit(my_mpi_rKind,ierror)

!Define the proper KIND for COMPLEX rKind MPI messages
CALL MPI_Type_create_f90_complex(precis,range,my_mpi_cKind,ierror)
CALL MPI_Type_commit(my_mpi_cKind,ierror)

ALLOCATE(countvec(0:num_cores-1))
ALLOCATE(offsetVec(0:num_cores-1))

END SUBROUTINE Initialize_MPI

SUBROUTINE Finalize_MPI()
!
!Purpose: Make sure all processors are finished and clean up MPI
!
IMPLICIT NONE

!Make sure that all processors are finished before exiting
	CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

CALL MPI_TYPE_FREE(my_mpi_rKind, ierror)
CALL MPI_TYPE_FREE(my_mpi_cKind, ierror)

CALL MPI_Finalize(ierror) !Clean up MPI
END SUBROUTINE Finalize_MPI


SUBROUTINE Setup_sites_MPI()
!
!Purpose: Divide the Gammas evenly among all processors for the site-parallel code
!
IMPLICIT NONE
INTEGER :: errorcode, i


IF(num_cores.gt.FLOOR(0.5*systemSize)) THEN
PRINT *, 'Number of processors cannot exceed systemSize/2!'
PRINT *, 'Program is aborting now!'
CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
END IF

!The number of solely owned Gammas is (L+1)/P-1
my_local_dim=FLOOR((systemSize+1)*1.0/(num_cores*1.0))-1

!The boundaries have only one shared Gamma, so add only 1
IF((my_rank==(num_cores-1)).or.(my_rank==0)) THEN
my_local_dim=my_local_dim+1
ELSE
!All other processors have two shared Gammas
my_local_dim=my_local_dim+2
END IF

!Find the bounds on this uniform grid
IF(my_rank==0) THEN
my_bounds(1)=1
my_bounds(2)=my_local_dim
ELSE IF(my_rank==(num_cores-1)) THEN
my_bounds(2)=systemSize
my_bounds(1)=systemSIze-my_local_dim+1
ELSE
!rank 0 portion
my_bounds(1)=my_local_dim-1
!Intermediate portions
DO i=1,my_rank-1
my_bounds(1)=my_bounds(1)+my_local_dim-1
END DO
my_bounds(2)=my_bounds(1)+my_local_dim-1
END IF


!Distribute the remainder
IF(my_rank.le.(MOD(systemSize+1,num_cores)-1)) THEN
my_local_dim=my_local_dim+1
END IF

!Re-evaluate the bounds
IF(my_rank==0) THEN
	my_bounds(1)=1
	my_bounds(2)=my_local_dim
ELSE IF(my_rank==(num_cores-1)) THEN
	my_bounds(2)=systemSize
	my_bounds(1)=systemSIze-my_local_dim+1
ELSE
!Add 1 for each processor that gained a Gamma
	DO i=0,my_rank-1
		IF(i.le.(MOD(systemSize+1,num_cores)-1)) THEN
		my_bounds(1)=my_bounds(1)+1
		my_bounds(2)=my_bounds(2)+1
		END IF
	END DO

	IF(my_rank.le.(MOD(systemSize+1,num_cores)-1)) THEN
	my_bounds(2)=my_bounds(2)+1
	END IF

END IF

END SUBROUTINE Setup_sites_MPI

SUBROUTINE AllocateGamLamParallel(Gammas, Lambdas, chi)
!
!Purpose: Allocate gammas and Lambdas based on a value of chi for the site-parallel code
!
IMPLICIT NONE  
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER, INTENT(IN) :: chi
INTEGER :: i	

ALLOCATE(Gammas(my_local_dim))
!Lambdas live on links, the extra 2 assist in computing two-site observables
ALLOCATE(Lambdas(my_local_dim+1) , STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate lambdas'
			END IF 
	DO i=1,my_local_dim
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
ALLOCATE(Lambdas(my_local_dim+1)%v(chi), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate lambdas'
			END IF 
Lambdas(my_local_dim+1)%v=0.0_rKind

END SUBROUTINE AllocateGamLamParallel

SUBROUTINE CopyGamLamParallel(GammasCopy, LambdasCopy, GammasOrig, LambdasOrig)
!
!Purpose: Copy Gammas and Lambdas from Orig to Copy for the site-parallel code
!
IMPLICIT NONE
TYPE(tensor), POINTER :: GammasCopy(:), GammasOrig(:)
TYPE(vector), POINTER :: LambdasCopy(:), LambdasOrig(:)
INTEGER :: i, alpha, j, beta, chimin
!If one of the Gammas is larger than the other, only sum up to the 
!maximum indices of the smaller

 
	chimin=MIN(SIZE(GammasCopy(1)%t,3),SIZE(GammasOrig(1)%t,3))
	DO i=1,(my_local_dim+1)
	LambdasCopy(i)%v=0.0_rKind
		DO alpha=1,chimin
		LambdasCopy(i)%v(alpha)=LambdasOrig(i)%v(alpha)
		END DO
	END DO
	
	DO i=1,my_local_dim
	GammasCopy(i)%t=0.0_rKind
		DO alpha=1,chimin
			DO beta=1,chimin
				DO j=1,localSize
					GammasCopy(i)%t(alpha,j,beta)=GammasOrig(i)%t(alpha,j,beta)
				END DO
			END DO
		END DO
	END DO
	
END SUBROUTINE CopyGamLamParallel


SUBROUTINE DeallocateGamLamParallel(Gammas, Lambdas)
!
!Purpose: Deallocate gammas and Lambdas for the site-parallel code
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER :: i
!Deallocate each site/link object
	DO	i=1,(my_local_dim)
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
	DEALLOCATE(Lambdas(my_local_dim+1)%v, STAT=statInt)
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
END SUBROUTINE DeallocateGamLamParallel

SUBROUTINE ConvertGamLamParallel(SerialGams, ParaGams, SerialLams, ParaLams)
!
!Purpose: Convert the serial Gammas and Lamdas into Parallel form
!
IMPLICIT NONE
TYPE(tensor), POINTER :: SerialGams(:), ParaGams(:)
TYPE(vector), POINTER :: SerialLams(:),ParaLams(:)	
INTEGER :: i


DO i=1,my_local_dim
ParaGams(i)%t=SerialGams(my_bounds(1)+i-1)%t
ParaLams(i)%v=SerialLams(my_bounds(1)+i-1)%v
END DO

ParaLams(my_local_dim+1)%v=SerialLams(my_bounds(1)+my_local_dim)%v

END SUBROUTINE ConvertGamLamParallel

SUBROUTINE ConvertGamLamLabParallel(SerialGams, ParaGams, SerialLams, ParaLams,SerialLabL, ParaLabL, SerialLabR, ParaLabR)
!
!Purpose: Convert the serial Gammas, Lamdas, and Labels into Parallel form
!
IMPLICIT NONE
TYPE(tensor), POINTER :: SerialGams(:), ParaGams(:)
TYPE(vector), POINTER :: SerialLams(:),ParaLams(:)
TYPE(vectorInt), POINTER :: SerialLabL(:),ParaLabL(:), SerialLabR(:),ParaLabR(:)	
INTEGER :: i


DO i=1,my_local_dim
ParaGams(i)%t=SerialGams(my_bounds(1)+i-1)%t
ParaLams(i)%v=SerialLams(my_bounds(1)+i-1)%v
ParaLabL(i)%vi=SerialLabL(my_bounds(1)+i-1)%vi
ParaLabR(i)%vi=SerialLabR(my_bounds(1)+i-1)%vi
END DO

ParaLams(my_local_dim+1)%v=SerialLams(my_bounds(1)+my_local_dim)%v
ParaLabL(my_local_dim+1)%vi=SerialLabL(my_bounds(1)+my_local_dim)%vi
ParaLabR(my_local_dim+1)%vi=SerialLabR(my_bounds(1)+my_local_dim)%vi

END SUBROUTINE ConvertGamLamLabParallel

SUBROUTINE RecordLambdasParallel(fileid, Lambdas,openKind)
!
!Purpose: Record Lambdas on a file whose ID is fileid.
!
IMPLICIT NONE
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER, INTENT(IN) :: fileid
INTEGER :: i, j, chi
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName
CHARACTER, INTENT(IN) :: openKind

!Write chi into a string
chi = SIZE(Lambdas(2)%v)
WRITE(mString,'(I4)') chi

!'S' specifies scientific notation output-human readable
IF(openKind=='S') THEN
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
!'B' specifies binary output-machine readable
ELSE IF(openKind=='B') THEN
fmtName='('//TRIM(ADJUSTL(mString))//'B80)'
ELSE
STOP "Unknown option in RecordLambdas"
END IF

	DO i=1,(my_local_dim+1)
		WRITE(UNIT=fileid, FMT=fmtname) (Lambdas(i)%v(j), j=1,chi)
	END DO
END SUBROUTINE RecordLambdasParallel

SUBROUTINE RecordGammasParallel(fileid, Gammas, openKind)
!
!Purpose: Record Gammas on a file whose ID is fileid.
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)	
INTEGER, INTENT(IN) :: fileid
INTEGER :: i,j,k,l,chi
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName
CHARACTER, INTENT(IN) :: openKind

chi=SIZE(Gammas(1)%t,3)
!Write localSize into a string
WRITE(mString,'(I4)') localSize

!'S' specifies scientific notation output-human readable
IF(openKind=='S') THEN
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
!'B' specifies binary output-machine readable
ELSE IF(openKind=='B') THEN
fmtName='('//TRIM(ADJUSTL(mString))//'B80)'
ELSE
STOP "Unknown option in RecordGammas"
END IF

	DO l=1,my_local_dim
		DO i=1,chi
			DO j=1,chi
				WRITE(UNIT=fileid, FMT=fmtName) (REAL(Gammas(l)%t(i,k,j)), k=1,localSize)
			END DO
		END DO
	END DO
	
	DO l=1,my_local_dim
		DO i=1,chi
			DO j=1,chi
				WRITE(UNIT=fileid, FMT=fmtName) (AIMAG(Gammas(l)%t(i,k,j)), k=1,localSize)
			END DO
		END DO
	END DO
		
END SUBROUTINE RecordGammasParallel

SUBROUTINE RecordLabelParallel(fileid, LabelLorR)
!
!Purpose: Record Parallel Label on a file whose ID is fileid.
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLorR(:)
INTEGER, INTENT(IN) :: fileid
INTEGER :: i, j, chi
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

chi = SIZE(LabelLorR(2)%vi)
WRITE(mString,'(I4)') chi

!Always integers	
fmtName='('//TRIM(ADJUSTL(mString))//'I16)'

	DO i=1,(my_local_dim+1)
		WRITE(UNIT=fileid, FMT=fmtname) (LabelLorR(i)%vi(j), j=1,chi)
	END DO
END SUBROUTINE RecordLabelParallel

SUBROUTINE readGammaLambdaParallel(lambdafileID, gammafileID,Gammas,Lambdas, openKind, chiNow)
!
!Purpose: Read Gammas and lambdas from file UNIT=*fileID into the allocated local tensors.
!
IMPLICIT NONE
INTEGER, INTENT(IN) ::  gammaFileID, lambdaFileID
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER, INTENT(IN) :: chiNow
CHARACTER, INTENT(IN) :: openKind
CHARACTER(4) :: mString
CHARACTER(64) :: fmtNameL, fmtNameG
COMPLEX(KIND=rKind) :: eye
REAL(KIND=rKind), ALLOCATABLE :: gammatemp(:)
INTEGER i,j,k,l


eye = CMPLX(0.0,1.0,KIND=rKind)

WRITE(mString,'(I4)') chiNow

!'S' specifies scientific notation output-human readable
IF(openKind=='S') THEN
fmtNameL='('//TRIM(ADJUSTL(mString))//'E30.15)'
!'B' specifies binary output-machine readable
ELSE IF(openKind=='B') THEN
fmtNameL='('//TRIM(ADJUSTL(mString))//'B80)'
ELSE
STOP "Unknown option in readGammaLambda"
END IF	

	WRITE(mString,'(I4)') localSize

IF(openKind=='S') THEN
fmtNameG='('//TRIM(ADJUSTL(mString))//'E30.15)'
!'B' specifies binary output-machine readable
ELSE IF(openKind=='B') THEN
fmtNameG='('//TRIM(ADJUSTL(mString))//'B80)'
ELSE
STOP "Unknown option in readGammaLambda"
END IF	


! Read Lambdas
	    DO i=1,(my_local_dim+1)
	    	READ(UNIT=lambdaFileID, FMT=fmtNamel) (Lambdas(i)%v(j), j=1,chiNow)
	    END DO

		ALLOCATE(gammatemp(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate gammatemp'
			END IF

! Read real part of Gammas
		DO l=1,my_local_dim
			DO i=1,chiNow
				DO j=1,chiNow
					READ(UNIT=gammaFileID, FMT=fmtNameG) (gammatemp(k), k=1,localSize)
					DO k=1,localSize
						Gammas(l)%t(i,k,j)=gammatemp(k)
					END DO
				END DO
			END DO
		END DO


! Read imaginary part of Gammas
		DO l=1,my_local_dim
			DO i=1,chiNow
				DO j=1,chiNow
					READ(UNIT=gammaFileID, FMT=fmtNameG) (gammatemp(k), k=1,localSize)
					DO k=1,localSize
						Gammas(l)%t(i,k,j)=Gammas(l)%t(i,k,j)+eye*gammatemp(k)
					END DO
				END DO
			END DO
		END DO


		DEALLOCATE(gammatemp, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate gammatemp'
			END IF
			
END SUBROUTINE readGammaLambdaParallel


SUBROUTINE readGammaLambdaLabelsParallel(lambdafileID, gammafileID,labelleftFileID, labelrightFileID,&
									Gammas,Lambdas,LabelLeft, LabelRight, openKind, chiNow)
!
!Purpose: Read Gammas and lambdas from file UNIT=fileID into the allocated local tensors.
!
IMPLICIT NONE
INTEGER, INTENT(IN) ::  gammaFileID, lambdaFileID,labelleftFileID, labelrightFileID
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: chiNow
CHARACTER, INTENT(IN) :: openKind
CHARACTER(4) :: mString
CHARACTER(64) :: fmtNameL, fmtNameG, fmtNameLabel
COMPLEX(KIND=rKind) :: eye
REAL(KIND=rKind), ALLOCATABLE :: gammatemp(:)
INTEGER i,j,k,l

eye = CMPLX(0.0,1.0,KIND=rKind)

WRITE(mString,'(I4)') chiNow

!'S' specifies scientific notation output-human readable
IF(openKind=='S') THEN
fmtNameL='('//TRIM(ADJUSTL(mString))//'E30.15)'
!'B' specifies binary output-machine readable
ELSE IF(openKind=='B') THEN
fmtNameL='('//TRIM(ADJUSTL(mString))//'B80)'
ELSE
STOP "Unknown option in readGammaLambda"
END IF	

fmtNameLabel='('//TRIM(ADJUSTL(mString))//'I16)'


	WRITE(mString,'(I4)') localSize

!'S' specifies scientific notation output-human readable
IF(openKind=='S') THEN
fmtNameG='('//TRIM(ADJUSTL(mString))//'E30.15)'
!'B' specifies binary output-machine readable
ELSE IF(openKind=='B') THEN
fmtNameG='('//TRIM(ADJUSTL(mString))//'B80)'
ELSE
STOP "Unknown option in readGammaLambda"
END IF	

! Read Lambdas
	    DO i=1,(my_local_dim+1)
	    	READ(UNIT=lambdaFileID, FMT=fmtNamel) (Lambdas(i)%v(j), j=1,chiNow)
	    END DO

		ALLOCATE(gammatemp(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate gammatemp'
			END IF

! Read real part of Gammas
		DO l=1,my_local_dim
			DO i=1,chiNow
				DO j=1,chiNow
					READ(UNIT=gammaFileID, FMT=fmtNameG) (gammatemp(k), k=1,localSize)
					DO k=1,localSize
						Gammas(l)%t(i,k,j)=gammatemp(k)
					END DO
				END DO
			END DO
		END DO

! Read imaginary part of Gammas
		DO l=1,my_local_dim
			DO i=1,chiNow
				DO j=1,chiNow
					READ(UNIT=gammaFileID, FMT=fmtNameg) (gammatemp(k), k=1,localSize)
					DO k=1,localSize
						Gammas(l)%t(i,k,j)=Gammas(l)%t(i,k,j)+eye*gammatemp(k)
					END DO
				END DO
			END DO
		END DO

		DEALLOCATE(gammatemp, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate gammatemp'
			END IF


! Read LabelLeft
		DO i=1,(my_local_dim+1)
	    	READ(UNIT=labelleftFileID, FMT=fmtNameLabel) (LabelLeft(i)%vi(j), j=1,chiNow)
	    END DO
! Read LabelRight
		DO i=1,(my_local_dim+1)
	    	READ(UNIT=labelrightFileID, FMT=fmtNameLabel) (LabelRight(i)%vi(j), j=1,chiNow)
	    END DO
			
END SUBROUTINE readGammaLambdaLabelsParallel

SUBROUTINE CheckpointParallel(basename,indexNow,Gammas,Lambdas,wherechar)
!
!Purpose: Save the current index, Gammas, and Lambdas so that a Parallel job may be restarted
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
CHARACTER(len=*), INTENT(IN) :: basename
INTEGER, INTENT(IN) :: indexNow
INTEGER :: my_unitg, my_unitl, my_unitc
CHARACTER(132) :: gamsname, lamsname, ckptname
CHARACTER(2), INTENT(IN) :: wherechar

my_unitg=1900+my_rank
my_unitl=2900+my_rank
my_unitc=3900+my_rank

CALL copyname(basename,gamsname)
CALL copyname(basename,lamsname)
CALL copyname(basename,ckptname)

CALL appendbasename(gamsname,'gams.dat')
CALL appendbasename(lamsname,'lams.dat')
CALL appendbasename(ckptname,'ckpt.dat')

CALL openUnit(gamsname,my_unitg)
CALL openUnit(lamsname,my_unitl)
CALL openUnit(ckptname,my_unitc)

CALL RecordLambdasParallel(my_unitl, Lambdas,ITPopenKind)
CALL RecordGammasParallel(my_unitg, Gammas,ITPopenKind)

WRITE(my_unitc,'(A2)') wherechar
WRITE(my_unitc,'(I4)') num_cores
WRITE(my_unitc,'(I4)') SIZE(Gammas(1)%t,1)
WRITE(my_unitc,'(I4)') indexNow


CLOSE(my_unitg)
CLOSE(my_unitl)
CLOSE(my_unitc)

END SUBROUTINE CheckpointParallel

CHARACTER(2) FUNCTION CheckCheckpoint(basename)
!
!Purpose: Test where the Checkpointed file basename was stopped
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: basename
CHARACTER(132) :: ckptname
INTEGER :: my_unitc
CHARACTER(2) wherechar

my_unitc=3900+my_rank
CALL copyname(basename,ckptname)
CALL appendbasename(ckptname,'ckpt.dat')

IF(CheckName(ckptname)) THEN
CALL openUnit(ckptname,my_unitc,'O')
REWIND(my_unitc)
READ(my_unitc,'(A2)') wherechar
ELSE
wherechar='NN'
END IF

CheckCheckpoint=wherechar

END FUNCTION CheckCheckpoint

SUBROUTINE RestartParallel(basename,indexNow,Gammas,Lambdas, exitstatus)
!
!Purpose: Restart a job from the checkpointed files basename*
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
CHARACTER(len=*), INTENT(IN) :: basename
INTEGER, INTENT(OUT) :: indexNow, exitstatus
INTEGER :: my_unitg, my_unitl, my_unitc,corescheck, chicheck
CHARACTER(132) :: gamsname, lamsname, ckptname
CHARACTER(2) :: wherechar

my_unitg=1900+my_rank
my_unitl=2900+my_rank
my_unitc=3900+my_rank

CALL copyname(basename,gamsname)
CALL copyname(basename,lamsname)
CALL copyname(basename,ckptname)

CALL appendbasename(gamsname,'gams.dat')
CALL appendbasename(lamsname,'lams.dat')
CALL appendbasename(ckptname,'ckpt.dat')

IF(CheckName(gamsname).and.CheckName(lamsName).and.CheckName(ckptname)) THEN
CALL openUnit(gamsname,my_unitg,'O')
CALL openUnit(lamsname,my_unitl,'O')
CALL openUnit(ckptname,my_unitc,'O')
REWIND(my_unitc)
READ(my_unitc,'(A2)') wherechar
READ(my_unitc,'(I4)') corescheck
READ(my_unitc,'(I4)') chicheck
READ(my_unitc,'(I4)') indexnow


	IF((corescheck==num_cores).and.(chicheck==SIZE(Gammas(1)%t,1))) THEN
	CALL readGammaLambdaParallel(my_unitl, my_unitg,Gammas,Lambdas, ITPopenKind, chicheck)
	ELSE
	PRINT *, 'Checkpointing parameters do not match.'
	PRINT *,  'exitstatus=1 returned!'
	PRINT *, ''
	exitstatus=1
	END IF
CLOSE(my_unitg)
CLOSE(my_unitl)
CLOSE(my_unitc)
exitstatus=0
ELSE
	PRINT *, 'Checkpointing files not found.'
	PRINT *,  'exitstatus=1 returned!'
	PRINT *, ''

exitstatus=1
END IF

END SUBROUTINE RestartParallel

SUBROUTINE CheckpointNCParallel(basename,indexNow,Gammas,Lambdas, LabelLeft,LabelRight,wherechar)
!
!Purpose: Save the current index, Gammas, and Lambdas so that a Parallel job may be restarted
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
CHARACTER(len=*), INTENT(IN) :: basename
INTEGER, INTENT(IN) :: indexNow
INTEGER :: my_unitg, my_unitl, my_unitc, my_unitll,my_unitlr
CHARACTER(132) :: gamsname, lamsname, ckptname, llname, lrname
CHARACTER(2), INTENT(IN) :: wherechar

my_unitg=1900+my_rank
my_unitl=2900+my_rank
my_unitc=3900+my_rank
my_unitll=4900+my_rank
my_unitlr=5900+my_rank

CALL copyname(basename,gamsname)
CALL copyname(basename,lamsname)
CALL copyname(basename,ckptname)
CALL copyname(basename,llname)
CALL copyname(basename, lrname)

CALL appendbasename(gamsname,'gams.dat')
CALL appendbasename(lamsname,'lams.dat')
CALL appendbasename(lamsname,'labell.dat')
CALL appendbasename(lamsname,'labelr.dat')
CALL appendbasename(ckptname,'ckpt.dat')

CALL openUnit(gamsname,my_unitg)
CALL openUnit(lamsname,my_unitl)
CALL openUnit(llname,my_unitll)
CALL openUnit(lrname,my_unitlr)
CALL openUnit(ckptname,my_unitc)

CALL RecordLambdasParallel(my_unitl, Lambdas,ITPopenKind)
CALL RecordGammasParallel(my_unitg, Gammas,ITPopenKind)
CALL RecordLabelParallel(my_unitll, LabelLeft)
CALL RecordLabelParallel(my_unitlr, LabelLeft)

WRITE(my_unitc,'(A2)') wherechar
WRITE(my_unitc,'(I4)') num_cores
WRITE(my_unitc,'(I4)') SIZE(Gammas(1)%t,1)
WRITE(my_unitc,'(I4)') indexNow

CLOSE(my_unitg)
CLOSE(my_unitl)
CLOSE(my_unitll)
CLOSE(my_unitlr)
CLOSE(my_unitc)

END SUBROUTINE CheckpointNCParallel

SUBROUTINE RestartNCParallel(basename,indexNow,Gammas,Lambdas,LabelLeft,LabelRight, exitstatus)
!
!Purpose: Restart a job from the checkpointed files basename*
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
CHARACTER(len=*), INTENT(IN) :: basename
INTEGER, INTENT(OUT) :: indexNow, exitstatus
INTEGER :: my_unitg, my_unitl, my_unitc,  my_unitll,my_unitlr,corescheck,chicheck
CHARACTER(132) :: gamsname, lamsname, ckptname, llname, lrname
CHARACTER(2) wherechar

my_unitg=1900+my_rank
my_unitl=2900+my_rank
my_unitc=3900+my_rank
my_unitll=4900+my_rank
my_unitlr=5900+my_rank

CALL copyname(basename,gamsname)
CALL copyname(basename,lamsname)
CALL copyname(basename,ckptname)
CALL copyname(basename,llname)
CALL copyname(basename,lrname)

CALL appendbasename(gamsname,'gams.dat')
CALL appendbasename(lamsname,'lams.dat')
CALL appendbasename(llname,'labell.dat')
CALL appendbasename(lrname,'labelr.dat')
CALL appendbasename(ckptname,'ckpt.dat')

IF(CheckName(gamsname).and.CheckName(lamsName).and.CheckName(ckptname).and.&
CheckName(llname).and.CheckName(lrname)) THEN

CALL openUnit(gamsname,my_unitg,'O')
CALL openUnit(lamsname,my_unitl,'O')
CALL openUnit(ckptname,my_unitc,'O')
CALL openUnit(llname,my_unitll, 'O')
CALL openUnit(lrname,my_unitlr, 'O')
REWIND(my_unitc)
READ(my_unitc,'(A2)') wherechar
READ(my_unitc,'(I4)') corescheck
READ(my_unitc,'(I4)') chicheck
READ(my_unitc,'(I4)') indexnow


	IF((corescheck==num_cores).and.(chicheck==SIZE(Gammas(1)%t,1))) THEN
	CALL readGammaLambdaLabelsParallel(my_unitl, my_unitg,my_unitll, my_unitlr,&
				Gammas,Lambdas,LabelLeft, LabelRight, ITPopenKind, chicheck)
	ELSE
	PRINT *, 'Checkpointing parameters do not match.'
	PRINT *,  'exitstatus=1 returned!'
	PRINT *, ''
	exitstatus=1
	END IF
CLOSE(my_unitg)
CLOSE(my_unitl)
CLOSE(my_unitll)
CLOSE(my_unitlr)
CLOSE(my_unitc)
exitstatus=0
ELSE
	PRINT *, 'Checkpointing files not found.'
	PRINT *,  'exitstatus=1 returned!'
	PRINT *, ''

exitstatus=1
END IF

END SUBROUTINE RestartNCParallel

SUBROUTINE AllocateLabelParallel(LabelLeft, LabelRight, chi)
!
!Purpose: Allocate labels needed for conserving a single Abelian symmetry for the site-parallel code
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN) :: chi
INTEGER :: i
!LableLeft(l)%vi(alpha) gives the conserved quantity associated with the alpha^th 
!left Schmidt vector for a bipartite splitting at link l, likewise for LabelRight
ALLOCATE(LabelLeft(my_local_dim+1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate labelLeft'
			END IF 
ALLOCATE(LabelRight(my_local_dim+1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate labelRight'
			END IF 
	DO i=1, (my_local_dim+1)
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
END SUBROUTINE AllocateLabelParallel	

SUBROUTINE CopyLabelParallel(LabLCopy, LabRCopy, LabLOrig, LabROrig)
!
!Purpose: Copy Labels from Orig to Copy for the site-parallel code
!		  Used in single Abelian symmetry conserving codes
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabLCopy(:), LabRCopy(:), LabLOrig(:), LabROrig(:)
INTEGER :: i, alpha, chimin
!If one of the Labels is larger than the other, only sum up to the 
!maximum indices of the smaller 	
chimin=MIN(SIZE(LabLCopy(1)%vi,1),SIZE(LabLOrig(1)%vi,1))
	DO i=1,(my_local_dim+1)
	LabLCopy(i)%vi=10000
	LabRCopy(i)%vi=10000
		DO alpha=1,chimin
		LabLCopy(i)%vi(alpha)=LabLOrig(i)%vi(alpha)
		LabRCopy(i)%vi(alpha)=LabROrig(i)%vi(alpha)
		END DO
	END DO
END SUBROUTINE CopyLabelParallel

SUBROUTINE DeallocateLabelParallel(LabelLeft, LabelRight)
!
!Purpose: Deallocate labels needed for conserving a single Abelian symmetry for the site-parallel code
!
IMPLICIT NONE
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER :: i
	DO i=1, (my_local_dim+1)
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
END SUBROUTINE DeallocateLabelParallel


SUBROUTINE ConstructPropagatorsParallel(H, U, dtodd, dteven)
!
!Purpose: Construct the Trotter-Suzuki propagator U from the Hamiltonian H for the site-parallel code
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(matrix), POINTER :: U(:)
COMPLEX(KIND=rKind), INTENT(IN) :: dtodd, dteven
INTEGER :: i

	DO i=1,(my_local_dim-1),2
IF(MOD(my_bounds(1),2).ne.0) THEN			
		CALL Matrix_Exponential(H(i)%m, U(i)%m, dtodd, localSize*localSize)!Same parity as lattice
ELSE
		CALL Matrix_Exponential(H(i)%m, U(i)%m, dteven, localSize*localSize)!Opposite parity
END IF
	END DO

	DO i=2,(my_local_dim-1),2
IF(MOD(my_bounds(1),2).ne.0) THEN			
		CALL Matrix_Exponential(H(i)%m, U(i)%m, dteven, localSize*localSize)!Same parity as lattice
ELSE
		CALL Matrix_Exponential(H(i)%m, U(i)%m, dtodd, localSize*localSize)!Opposite parity
END IF
	END DO
END SUBROUTINE ConstructPropagatorsParallel


SUBROUTINE ProductStateMPDParallel(Gammas, Lambdas, carray)
!
!Purpose: Construct the Vidal decomposition of a product state whose coefficients
!		  are stored in carray for the site-parallel code
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
	
	DO i=1,my_local_dim,1
	Gammas(i)%t = CMPLX(0.0, KIND=rKind)
	Lambdas(i)%v = 0.0_rKind
	Lambdas(i)%v(1) = 1.0_rKind ! Assign the first component of each lambda the value 1, as this is a product state.
		DO j = 1, localSize
		Gammas(i)%t(1, j, 1) = carray(j, i+my_bounds(1)-1) ! Assign the alpha=1 values of Gammas to be the coefficients of each on-site state.
		END DO
	END DO
	Lambdas(my_local_dim+1)%v = 0.0_rKind
	Lambdas(my_local_dim+1)%v(1) = 1.0_rKind
END SUBROUTINE ProductStateMPDParallel


SUBROUTINE ProductStateLabelsParallel(LabelLeft, LabelRight, carray,intDegFree)
!
!Purpose: Construct the lists of number conserving labels of a product state whose coefficients
!		  are stored in carray for the site-parallel code
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
DO i=1,my_local_dim+1
	DO j=1,localSize
		IF(ABS(carray(j,i+my_bounds(1)-1)).ne.0.0_rKind) THEN
		LabelLeft(i+1)%vi(1)=LabelLeft(i)%vi(1)+Conserv%vi(j)
		EXIT
		END IF
	END DO
END DO

ELSE
LabelLeft(1)%vi(1)=0.0_rKind
DO i=1,my_local_dim+1
	DO j=1,localSize
		IF(ABS(carray(j, i+my_bounds(1)-1)).ne.0.0_rKind) THEN
		LabelLeft(i+1)%vi(1)=LabelLeft(i)%vi(1)+j-1
		END IF
	END DO
END DO

END IF

!Construct LabelRight from LabelLeft
DO i=1,my_local_dim+1
	LabelRight(i)%vi(1)=totNum-LabelLeft(i)%vi(1)
END DO
	
END SUBROUTINE ProductStateLabelsParallel


SUBROUTINE AllStatesParallel(Gammas, Lambdas)
!
!Purpose: Creates an initial state that is a product state of local states 
! which contain all possible states in the same amount for site-parallel codes.  
! Used as initial ITP state for number non-conserving code.
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)	
INTEGER :: i, j
!The state is a product state, so the lambdas corresponding to only one Schmidt index are nonzero
	DO i=1,my_local_dim
	Lambdas(i)%v(1)=1.0_rKind
	!Each state is weighted equally-by normalization by 1/sqrt(d)
		DO j=1,localSize
		Gammas(i)%t(1,j,1)=CMPLX((1.0_rKind)/SQRT(localSize*1.0_rKind),KIND=rKind)
		END DO
	END DO
	Lambdas(my_local_dim+1)%v(1)=1.0_rKind
END SUBROUTINE AllStatesParallel

SUBROUTINE InitialSetNCParallel(Gammas, Lambdas, LabelLeft, LabelRight, intDegFree)
!
!Purpose: Creates an initial state consistent with number conservation for site-parallel codes
!
!		  The algorithm places the particles in a "wedding cake" structure
!		  with the greatest number of particles in the center of the cake ("tops")
!		  A lesser number surrounding this peak ("center"), and no particles in the gap between
!		  the cake and the edge ("hole")
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
REAL(KIND=rKind) :: minc, maxc
INTEGER :: i, j, k, l, hole, tops, center,nj, dummypass


IF(totNum==0) THEN
	DO i=1,my_local_dim+1,1
		LabelLeft(i)%vi(1)=0
		LabelRight(i)%vi(1)=0
		Lambdas(i)%v(1)=1.0_rKind
	END DO
	DO i=1,my_local_dim,1
		Gammas(i)%t(1,1,1)=1.0_rKind
	END DO


ELSE

! If the number of sites with the least number of particles (hole) is even
	IF((MOD(systemSize-MOD(totNum,systemSize),2)==0).OR. &
! or there is no hole (uniform filling)
	(MOD(totNum,systemSize)==0)) THEN
		!!If we have unit filling there is no hole
		IF(MOD(totNum,systemSize)==0) THEN
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
	center=FLOOR(REAL(totNum,KIND=rKind)/REAL(systemSize,KIND=rKind)-10.0**(-8))+1


!Loop over processors to initialize LabelLeft
DO k=0,num_cores-1
IF(my_rank==k) THEN

	IF((my_rank.ne.0).AND.(num_cores.ne.1)) THEN
	!Recieve cumulative value from the last core
	CALL MPI_Recv(dummypass,1,my_mpi_rKind,k-1,my_rank-1, MPI_COMM_WORLD,my_status,ierror)
	END IF
	
			
! LabelLeft in 0<=link<=hole
	DO i=1,(hole+1),1
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.(my_bounds(2)+1))) THEN
		j=i-my_bounds(1)+1
			IF(my_rank==0) THEN
				!Boundary condition
				IF(j==1) THEN
				LabelLeft(j)%vi(1)=0
				ELSE
				LabelLeft(j)%vi(1)=LabelLeft(j-1)%vi(1)+center-1
				END IF
			ELSE
			!Count the cumulative number to the left
				!If my_rank owns the previous label then sum it
				IF(j.ne.1) THEN
				LabelLeft(j)%vi(1)=LabelLeft(j-1)%vi(1)+center-1
				ELSE
				!Else use the result from the previous processor
				LabelLeft(j)%vi(1)=dummypass
				END IF
			END IF
		END IF
	END DO


	
! LabelLeft in hole+1<=link<=hole+top
	DO i=(hole+2),(hole+tops+1),1
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.(my_bounds(2)+1))) THEN
		j=i-my_bounds(1)+1
			IF(my_rank.ne.0) THEN
				IF(j.ne.1) THEN
				LabelLeft(j)%vi(1)=LabelLeft(j-1)%vi(1)+center
				ELSE
				LabelLeft(j)%vi(1)=dummypass
				END IF
			ELSE
			LabelLeft(j)%vi(1)=LabelLeft(j-1)%vi(1)+center
			END IF
		END IF
	END DO
		
! LabelLeft in hole+top+1<=link<=systemSize+1
	DO i=(hole+tops+2),(systemSize+1),1
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.(my_bounds(2)+1))) THEN
		j=i-my_bounds(1)+1
			IF(my_rank.ne.0) THEN
				IF(j.ne.1) THEN
				LabelLeft(j)%vi(1)=LabelLeft(j-1)%vi(1)+center-1
				ELSE
				LabelLeft(j)%vi(1)=dummypass
				END IF
			ELSE 
			LabelLeft(j)%vi(1)=LabelLeft(j-1)%vi(1)+center-1
			END IF
		END IF
	END DO
		

	IF((my_rank.ne.(num_cores-1)).AND.(num_cores.ne.1)) THEN
	dummypass=LabelLeft(my_local_dim)%vi(1)
	!Send cumulative value to the next core
	CALL MPI_Send(dummypass,1,my_mpi_rKind,k+1,my_rank, MPI_COMM_WORLD,ierror)
	END IF

END IF
END DO


! Construct LabelRight from LabelLeft
! Construct Lambdas
	DO i=1,(systemSize+1),1
	j=i-my_bounds(1)+1
		IF((j.ge.1).AND.(j.le.(my_local_dim+1))) THEN
		LabelRight(j)%vi(1)=totNum-LabelLeft(j)%vi(1)
		Lambdas(j)%v(1)=1.0_rKind
		END IF
	END DO


!Internal degree(s) of freedom present	
IF(PRESENT(intDegFree)) THEN

!Find the number of states with number=center-1 and center particles
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
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.my_bounds(2))) THEN
		l=i-my_bounds(1)+1
			DO j=1,localSize,1
			nj=Conserv%vi(j)
				IF(nj==center-1) THEN
				Gammas(l)%t(1,j,1)=CMPLX(1.0/SQRT(minc),KIND=rKind)
				END IF
			END DO		
		END IF
		END DO
! Gammas in hole+1<=site<=hole+top
		DO i=hole+1,(hole+tops),1
		!Sum over internal degree of freedom, weighting each democratically
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.my_bounds(2))) THEN
		l=i-my_bounds(1)+1
			DO j=1,localSize,1
			nj=Conserv%vi(j)
				IF(nj==center) THEN
				Gammas(l)%t(1,j,1)=CMPLX(1.0/SQRT(maxc),KIND=rKind)
				END IF
			END DO	
		END IF
		END DO
! Gammas in hole+top+1<=site<=systemSize
		DO i=(hole+tops+1),systemSize,1
		!Sum over internal degree of freedom, weighting each democratically
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.my_bounds(2))) THEN
		l=i-my_bounds(1)+1
			DO j=1,localSize,1
			nj=Conserv%vi(j)
				IF(nj==center-1) THEN
				Gammas(l)%t(1,j,1)=CMPLX(1.0/SQRT(minc),KIND=rKind)
				END IF
			END DO		
		END IF	
		END DO


!Internal degree(s) of freedom absent
ELSE

! Construct Gammas
! Gammas in 1<=site<=hole
	DO i=1,hole,1
	!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.my_bounds(2))) THEN
		j=i-my_bounds(1)+1
		Gammas(j)%t(1,center,1)=CMPLX(1.0,KIND=rKind)
		END IF
	END DO
! Gammas in hole+1<=site<=hole+top
	DO i=hole+1,(hole+tops),1
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.my_bounds(2))) THEN
		j=i-my_bounds(1)+1
		Gammas(j)%t(1,center+1,1)=CMPLX(1.0,KIND=rKind)
		END IF
	END DO
! Gammas in hole+top+1<=site<=systemSize
	DO i=(hole+tops+1),systemSize,1
		!If my_rank owns the current site then perform the operations
		IF((i.ge.my_bounds(1)).AND.(i.le.my_bounds(2))) THEN
		j=i-my_bounds(1)+1
		Gammas(j)%t(1,center,1)=CMPLX(1.0,KIND=rKind)
		END IF
	END DO
END IF

END IF
END SUBROUTINE InitialSetNCParallel


SUBROUTINE HamiltonianBoseHubbardParallel(H, jTunn, U0, mu0,  V0, extPot)
!
!Purpose: Construct the (extended) Bose-Hubbard Hamiltonian (see manual) in TEBD form
!Here jTunn is the tunneling energy, U0 is the on-site repulsion, mu0 is the chemical potential, 
!V0 is the nearest neighbor repulsion, and extPot is an optional vector argument specifing an
!arbitrary site-dependent trap.
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) ::  jTunn, U0, mu0, V0
REAL(KIND=rKind), INTENT(IN), OPTIONAL :: extPot(:)
INTEGER :: i

	DO i=1,(my_local_dim-1)
H(i)%m=0.0_rKind
!One site operations
!"Left" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of U/2(n(n-1))
H(i)%m=H(i)%m+0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!"Left" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Right" part of mu*n
H(i)%m=H(i)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))

!Two site operations
!Tunneling part
H(i)%m=H(i)%m-jTunn*t_op%mr

!Nearest-neighbor part
H(i)%m=H(i)%m+V0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
	END DO

!Fox box boundary conditions first and last site get additional left and right
!one site operations, respectively
IF(my_rank==0) THEN
!"Left" part of U/2(n(n-1)) for first site
H(1)%m=H(1)%m+0.25_rKind*U0*tensorProd(MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)),one_op%mr)&
-0.25_rKind*U0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
!"Left" part of mu*n for first site
H(1)%m=H(1)%m-0.5_rKind*mu0*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
END IF

IF(my_rank==num_cores-1) THEN
!"Right" part of U/2(n(n-1)) for last site
H(my_local_dim-1)%m=H(my_local_dim-1)%m+0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))&
-0.25_rKind*U0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
!"Right" part of mu*n for last site
H(my_local_dim-1)%m=H(my_local_dim-1)%m-0.5_rKind*mu0*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
END IF

!Application of an arbitrary external potential
IF(PRESENT(extPot)) THEN

DO i=1,(my_local_dim-1)
H(i)%m=H(i)%m+0.5_rKind*extPot(i+my_bounds(1)-1)*(tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr))+0.5_rKind*extPot(i+my_bounds(1))*(tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr)))
END DO
	IF(my_rank==0) THEN
		H(1)%m=H(1)%m+0.5_rKind*extPot(1)*tensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),one_op%mr)
	END IF
	IF(my_rank==num_cores-1) THEN
		H(my_local_dim-1)%m=H(my_local_dim-1)%m+0.5_rKind*extPot(systemSize)*tensorProd(one_op%mr,MATMUL(TRANSPOSE(a_op%mr),a_op%mr))
	END IF

END IF

END SUBROUTINE HamiltonianBoseHubbardParallel


SUBROUTINE HamiltonianRotationTIParallel(H)
!
!Purpose: Construct the TEBD form of Rotational Hamiltonian for site-parallel code: Time-independent version for ITP
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
INTEGER :: i

DO i = 1, (my_local_dim-1)
	H(i)%m =-ttot_opR%mr+dipdip_opR%mr+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr)
END DO
IF(my_rank==0) THEN
	H(1)%m = H(1)%m+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)
END IF
IF(my_rank==num_cores-1) THEN
	H(my_local_dim-1)%m = H(my_local_dim-1)%m+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr)
END IF
END SUBROUTINE HamiltonianRotationTIParallel


SUBROUTINE HamiltonianRotationTDParallel(H,time)
!
!Purpose: Construct the TEBD form of Rotational Hamiltonian for site-parallel code: Time-dependent version for RTP
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKind), INTENT(IN) ::   time
REAL (KIND=rKind) :: ACfac
INTEGER :: i
	ACfac=2.0_rKind*eAC*sin(omega*time)

	DO i = 1, (my_local_dim-1)
		H(i)%m =-ttot_opR%mr+dipdip_opR%mr+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr) &
		+0.5_rKind*TensorProd(ACfac*EAC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(one_op%mr,ACfac*EAC_opR%mr)
	END DO
IF(my_rank==0) THEN
	H(1)%m = H(1)%m+0.5_rKind*TensorProd(EDC_opR%mr,one_op%mr)+0.5_rKind*TensorProd(ACfac*EAC_opR%mr,one_op%mr)
END IF
IF(my_rank==num_cores-1) THEN
	H(my_local_dim-1)%m = H(my_local_dim-1)%m+0.5_rKind*TensorProd(one_op%mr,EDC_opR%mr)+0.5_rKind*TensorProd(one_op%mr,ACfac*EAC_opR%mr)
END IF

END SUBROUTINE HamiltonianRotationTDParallel


SUBROUTINE TrotterStep2ndOrderParallel(Udt, Gammas, Lambdas, totalTruncerr)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2) in the site-parallel code
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
REAL(KIND=rKind) :: trunctemp
COMPLEX(KIND=rKind) :: dumGam(SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1)), dumLam(SIZE(Lambdas(1)%v,1))
INTEGER :: i,j,k, chiHere, gamSize, lamSize
		totalTruncerr = 0.0_rKind

gamSize=SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1)
chiHere=SIZE(Gammas(1)%t,1)

!!!BEGIN Operate Exp(-i Hodd dt/2)
	IF(MOD(my_bounds(1),2).ne.0) THEN
	DO i=1,my_local_dim-1,2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!Case 1:If we are updating odd sites and begin on an odd site, send the first Gamma and second Lambda to the preceding processor
IF(my_rank.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,100+my_rank, MPI_COMM_WORLD,ierror)
END IF

!Case 2:If we are updating odd sites and begin on an even site, receive the last Gamma and next to last lambda from the preceding processor	
	ELSE
	DO i=2,my_local_dim-1,2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
	IF(my_rank.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,700+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		END DO
				
	END IF

END IF

!Receive Case 1
IF((MOD(my_bounds(2),2).ne.0).AND.(my_rank.ne.num_cores-1)) THEN
CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,100+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		END DO
		
END IF

!Send Case 2
IF((MOD(my_bounds(2),2)==0).AND.(my_rank.ne.num_cores-1)) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,700+my_rank, MPI_COMM_WORLD,ierror)
END IF


!!!END Operate Exp(-i Hodd dt/2)


!!! Operate Exp(-i Heven dt)
	IF(MOD(my_bounds(1),2)==0) THEN
	DO i=1,my_local_dim-1,2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
	
!Case 1:If we are updating even sites and begin on an even site, send the first Gamma and second Lambda to the preceding processor
IF(my_rank.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)

END IF


!Case 2:If we are updating even sites and begin on an odd site, receive the last Gamma and next to last lambda from the preceding processor	
	ELSE
	DO i=2,my_local_dim-1,2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
		IF(my_rank.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,800+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		END DO
	END IF

	
	END IF



!Receive Case 1
IF((MOD(my_bounds(2),2)==0).AND.(my_rank.ne.num_cores-1)) THEN
CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		END DO
	
END IF


!Send Case 2
IF((MOD(my_bounds(2),2).ne.0).AND.(my_rank.ne.num_cores-1)) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,800+my_rank, MPI_COMM_WORLD,ierror)
END IF

!!!END Operate Exp(-i Heven dt)


!!!BEGIN Operate Exp(-i Hodd dt/2)
	IF(MOD(my_bounds(1),2).ne.0) THEN
	DO i=1,my_local_dim-1,2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!Case 1:If we are updating odd sites and begin on an odd site, send the first Gamma and second Lambda to the preceding processor
IF(my_rank.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,300+my_rank, MPI_COMM_WORLD,ierror)
END IF

!Case 2:If we are updating odd sites and begin on an even site, receive the last Gamma and next to last lambda from the preceding processor	
	ELSE
	DO i=2,my_local_dim-1,2
		CALL TwoSiteOp(i, Udt(i)%m, Gammas, Lambdas, trunctemp)
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
	IF(my_rank.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,900+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		END DO
				
	END IF

	END IF

!Receive Case 1
IF((MOD(my_bounds(2),2).ne.0).AND.(my_rank.ne.num_cores-1)) THEN
CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,300+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		END DO
		
END IF

!Send Case 2
IF((MOD(my_bounds(2),2)==0).AND.(my_rank.ne.num_cores-1)) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,900+my_rank, MPI_COMM_WORLD,ierror)
END IF
trunctemp=totaltruncerr
	CALL MPI_REDUCE(trunctemp,totaltruncerr,1,my_mpi_rkind,MPI_SUM,0,MPI_COMM_WORLD,ierror)

END SUBROUTINE TrotterStep2ndOrderParallel


SUBROUTINE CanonicalFormAllParallel(Gammas,Lambdas)
!
!Purpose: Make all Bipartite splittings canonical. 
!Used to reorthogonalize the Schmidt basis after an ITP timestep.
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind) :: dumGam(SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1))
INTEGER :: k,i,j,siteInd, gamSize, chiHere

gamSize=SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1)
chiHere=SIZE(Gammas(1)%t,1)


!Sweep forwards
DO siteInd=0,num_cores-1

IF(my_rank==siteInd) THEN

!If you are not the first processor, recieve your updated leftmost local tensors from the previous processor
	IF(siteInd.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
		my_rank-1,100+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
			DO i=1,SIZE(Gammas(1)%t,1)
				DO j=1,SIZE(Gammas(1)%t,2)
					DO k=1,SIZE(Gammas(1)%t,3)
				Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
					END DO
				END DO
			END DO
			DO i=1,SIZE(Lambdas(1)%v,1)
			Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
			END DO
	END IF


!put all owned sites into canonical form
		DO k=1,(my_local_dim-1)
	CALL CanonicalForm(Lambdas(k)%v,Gammas(k)%t,Lambdas(k+1)%v,Gammas(k+1)%t,Lambdas(k+2)%v)
		END DO

!If you are not the last processor, send your rightmost local tensors to the next processor
	IF(siteInd.ne.num_cores-1) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		END DO
	CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,100+my_rank, MPI_COMM_WORLD,ierror)
	END IF


END IF



END DO


!Sweep backwards
DO siteInd=(num_cores-1),0,(-1)

IF(my_rank==siteInd) THEN

!If you are not the last processor, recieve your updated rightmost local tensors from the next processor
	IF(siteInd.ne.num_cores-1) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
		my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
			DO i=1,SIZE(Gammas(1)%t,1)
				DO j=1,SIZE(Gammas(1)%t,2)
					DO k=1,SIZE(Gammas(1)%t,3)
				Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
					END DO
				END DO
			END DO
			DO i=1,SIZE(Lambdas(1)%v,1)
			Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
			END DO
	END IF

!put all owned sites into canonical form
	DO k=(my_local_dim-1),1,(-1)
		CALL CanonicalForm(Lambdas(k)%v,Gammas(k)%t,Lambdas(k+1)%v,Gammas(k+1)%t,Lambdas(k+2)%v)
	END DO

!If you are not the first processor, send your leftmost local tensors to the previous processor
	IF(siteInd.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		END DO
	CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)
	END IF

!If you are not the last processor, send your updated rightmost tensors to the next processor
	IF(siteInd.ne.num_cores-1) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		END DO
	CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,300+my_rank, MPI_COMM_WORLD,ierror)
	END IF

!If you are not the first site, receive the your updated leftmost tensors from the previous processor
	IF(siteInd.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
		my_rank-1,300+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
			DO i=1,SIZE(Gammas(1)%t,1)
				DO j=1,SIZE(Gammas(1)%t,2)
					DO k=1,SIZE(Gammas(1)%t,3)
				Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
					END DO
				END DO
			END DO
			DO i=1,SIZE(Lambdas(1)%v,1)
			Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
			END DO
	END IF


END IF


END DO

END SUBROUTINE CanonicalFormAllParallel

SUBROUTINE ImagTimePropParallel(H, GammasOuter, LambdasOuter, chiIn, intDegFree,basename)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!and is the internal degree of freedom component whose number we wish to see (should be the ground state)
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: chiIn
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: GammasOuter(:), GammasInner(:)
TYPE(vector), POINTER :: LambdasOuter(:), LambdasInner(:)
TYPE(matrix), POINTER :: Uitp(:)
INTEGER, INTENT(IN),OPTIONAL :: intDegFree
COMPLEX(KIND=rKind) :: eye, dtodd, dteven
REAL(KIND=rKind) :: LrgstLamsBfr(my_local_dim+1), LrgstLamsAft(my_local_dim+1), GapBfrAft(my_local_dim+1) 
REAL(KIND=rKind) :: numberInner, cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l,siteInd, chi, chiInc
REAL(KIND=8) :: my_max(2), global_max(2)
CHARACTER(len=*),INTENT(INOUT), OPTIONAL :: basename
INTEGER :: my_loc(1), jstart, exitstatus, chicheckout
CHARACTER(2) :: wherechar

!!! Define the trotter time steps for ITP. 
	eye=CMPLX(0.0,1.0,KIND=rKind)
	dtodd=-eye*dtITP/2.0
	dteven=-eye*dtITP
!!! Construct the Imaginary time propagator
	CALL AllocateOps(Uitp,my_local_dim-1,localSize*localSize)
	CALL ConstructPropagatorsParallel(H, Uitp, dtodd, dteven)
!!! This 'if' statement is for preventing the chi increment to be 0 when chiIn=chiMax.		
	IF(chiIn==chiMax) THEN
		chiInc=1
	ELSE
		chiInc=chiMax-chiIn
	END IF

!!! We conduct the ITP for both chi=chiMin and chiMax for efficiency.
!!! The ITP for chi=chiMin usually creates a good initial state for the ITP for chi=chiMax. 
jstart=1
!Restart if possible
IF(restartSwitch) THEN
IF(.NOT.PRESENT(basename)) THEN
PRINT *,'A restart base filename must be supplied for restart option!'
STOP
END IF

IF(CheckCheckpoint(basename)=='I1') THEN
	CALL RestartParallel(basename,jstart,GammasOuter,LambdasOuter, exitstatus)
		IF(exitstatus==1) THEN	
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
		jstart=1
		END IF
ELSE IF(CheckCheckpoint(basename)=='I2') THEN
	CALL AllocateGamLamParallel(GammasInner, LambdasInner, chiMax)
	CALL CopyGamLamParallel(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
	CALL DeallocateGamLamParallel(GammasOuter, LambdasOuter)
	CALL RestartParallel(basename,jstart,GammasInner,LambdasInner, exitstatus)
		CALL AllocateGamLamParallel(GammasOuter, LambdasOuter, chiMax)
		CALL CopyGamLamParallel(GammasOuter, LambdasOuter, GammasInner, LambdasInner)
		CALL DeallocateGamLamParallel(GammasInner, LambdasInner)
		IF(exitstatus==1) THEN	
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
		jstart=1
		ELSE
		chiIn=chiMax
		END IF
ELSE
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
jstart=1
END IF

END IF

			
	DO i=chiIn,chiMax,chiInc
		chi=i

IF(print_switch) THEN
	PRINT *, 'chi', chi
END IF

!!! The convergence criterion in the ITP is tighter for chi=chiMax than for chi=chiMin.
	IF((chiIn==chiMax).OR.(chi==chiMax)) THEN
		convCr=convCriterion2
	ELSE
		convCr=convCriterion1
	END IF
	
!!! Initialize the parameters for SVD.
CALL SVDInit(chi)

!!! Initialize MPS for ITP.			
		CALL AllocateGamLamParallel(GammasInner, LambdasInner, chi)
		CALL CopyGamLamParallel(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
		CALL DeallocateGamLamParallel(GammasOuter, LambdasOuter)

!Store the largest lambda of each splitting initially
LrgstLamsBfr=0.0_rKind
		DO l=1,my_local_dim-1
			LrgstLamsBfr(l)=LambdasInner(l)%v(1)
		END DO

IF(my_rank==num_cores-1) THEN
		DO l=my_local_dim,my_local_dim+1
			LrgstLamsBfr(l)=LambdasInner(l)%v(1)
		END DO
END IF

!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
!!! or until the convergence criterion is satisfied.
		DO j=jstart,maxITPsteps
			IF(ckptSwitch) THEN
				IF(MOD(j,stepsforckpt)==0) THEN
			IF(chi==chiMax) THEN
			CALL CheckpointParallel(basename,j,GammasInner,LambdasInner, 'I2')
			ELSE
			CALL CheckpointParallel(basename,j,GammasInner,LambdasInner, 'I1')
			END IF
				END IF
			END IF
		
!!! This 'if' statement is for judging the convergence.
!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
			IF(MOD(j,stepsForJudge)==0) THEN

!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN
			CALL TotalNumberParallel(numberInner, GammasInner, LambdasInner, intDegFree)
!Internal degree(s) of freedom absent			
ELSE
			CALL TotalNumberParallel(numberInner, GammasInner, LambdasInner)
END IF
			CALL TotalEnergyParallel(energy,H, GammasInner, LambdasInner)



LrgstLamsAft=0.0_rKind
GapBfrAft=0.0_rKind
!Store the largest lambda of each splitting after stepsForJudge time steps
				DO l=1,my_local_dim-1
					LrgstLamsAft(l)=LambdasInner(l)%v(1)
!Find the percent differences of each lambda
					GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
				END DO

IF(my_rank==num_cores-1) THEN
		DO l=my_local_dim,my_local_dim+1
			LrgstLamsAft(l)=LambdasInner(l)%v(1)
			GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
		END DO
END IF

!Have all processors compute their maximum lambda difference and its location
my_loc=MAXLOC(GapBfrAft)
my_max(1)=REAL(MAXVAL(GapBfrAft), KIND=8)
my_max(2)=my_rank

!Have the master find the global max and its location
global_max=0.0_8
CALL MPI_REDUCE(my_max,global_max,2,MPI_2DOUBLE_PRECISION,&
MPI_MAXLOC,0,MPI_COMM_WORLD,ierror)


!Have the master broadcast this information to all cores
CALL MPI_BCast(global_max,2,my_mpi_rKind,0,MPI_COMM_WORLD,ierror)
testTol=global_max(1)

IF(my_rank==0) THEN
IF(PRESENT(intDegFree)) THEN
				PRINT *, 'number in the',intDegFree,'mode', numberInner, 'Energy',energy
ELSE
				PRINT *, 'number', numberInner, 'Energy',energy
END IF
END IF

IF(my_rank==FLOOR(global_max(2))) THEN
	PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(my_loc(1))%v(1), &
				'found at position', my_loc(1), 'on core #',my_rank
	PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr

END IF

				IF(testTol < convCr) EXIT
!Store the new largest lambda of each splitting
LrgstLamsBfr=0.0_rKind
		DO l=1,my_local_dim-1
			LrgstLamsBfr(l)=LrgstLamsAft(l)
		END DO

IF(my_rank==num_cores-1) THEN
		DO l=my_local_dim,my_local_dim+1
			LrgstLamsBfr(l)=LrgstLamsAft(l)
		END DO
END IF


			END IF

!Time step
			CALL TrotterStep2ndOrderParallel(Uitp, GammasInner, LambdasInner, totalTruncerr)
!Reorthogonalize
			CALL CanonicalFormAllParallel(GammasInner,LambdasInner)

		END DO

!!! Deallocate the parameters for SVD.
CALL SVDFinish()

!!! Reset GammasOuter and LambdasOuter.		
		CALL AllocateGamLamParallel(GammasOuter, LambdasOuter, chi)
		CALL CopyGamLamParallel(GammasOuter, LambdasOuter, GammasInner, LambdasInner)
!Clean up
		CALL DeallocateGamLamParallel(GammasInner, LambdasInner)
	jstart=1
	END DO
		
!	chiIn=chi	
	CALL DeallocateOps(Uitp,my_local_dim-1) 

END SUBROUTINE ImagTimePropParallel

SUBROUTINE TrotterStep2ndOrderNCParallel(Udt, Gammas, Lambdas, LabelLeft, LabelRight, totalTruncerr, intDegFree)
!
!Purpose: Propagate a TEBD form wavefunction via the 2nd order Trotter expansion:
! Exp(-i Hodd dt/2) Exp(-i Heven dt) Exp(-i Hodd dt/2)
!
!See manual for more detail
!
IMPLICIT NONE
!!! Udt is a site-indexed array of time-propagation operators. Udt(i) is Exp(-i Hodd dt/2) when i is odd,
!!! while Udt(i) is Exp(-i Heven dt) when i is even
TYPE(matrix), POINTER :: Udt(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
REAL(KIND=rKind), INTENT(OUT) :: totalTruncerr
REAL(KIND=rKind) :: trunctemp
COMPLEX(KIND=rKind) :: dumGam(SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1)), dumLam(SIZE(Lambdas(1)%v,1))
INTEGER :: dumLab(2*SIZE(LabelLeft(1)%vi,1))
INTEGER :: i,j,k, chiHere, gamSize, labSize
INTEGER, INTENT(IN), OPTIONAL :: intDegFree

		totalTruncerr = 0.0_rKind

gamSize=SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1)
chiHere=SIZE(Gammas(1)%t,1)
labSize=2*SIZE(LabelLeft(1)%vi,1)


!!!BEGIN Operate Exp(-i Hodd dt/2)
	IF(MOD(my_bounds(1),2).ne.0) THEN
	DO i=1,my_local_dim-1,2
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF	
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!Case 1:If we are updating odd sites and begin on an odd site, send the first Gamma and second Lambda to the preceding processor
IF(my_rank.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		dumLab(i)=LabelLeft(2)%vi(i)
		dumLab(chiHere+i)=LabelRight(2)%vi(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,100+my_rank, MPI_COMM_WORLD,ierror)

CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank-1,1000+my_rank, MPI_COMM_WORLD,ierror)
END IF

!Case 2:If we are updating odd sites and begin on an even site, receive the last Gamma and next to last lambda from the preceding processor	
	ELSE
	DO i=2,my_local_dim-1,2
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF	
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
	IF(my_rank.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,700+my_rank-1, MPI_COMM_WORLD,my_status,ierror)

	CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
	my_rank-1,7000+my_rank-1, MPI_COMM_WORLD,my_status,ierror)

		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		LabelLeft(1)%vi(i)=dumLab(i)
		LabelRight(1)%vi(i)=dumLab(chiHere+i)
		END DO
				
	END IF

	END IF

!Receive Case 1
IF((MOD(my_bounds(2),2).ne.0).AND.(my_rank.ne.num_cores-1)) THEN
CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,100+my_rank+1, MPI_COMM_WORLD,my_status,ierror)

CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,1000+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		LabelLeft(my_local_dim+1)%vi(i)=dumLab(i)
		LabelRight(my_local_dim+1)%vi(i)=dumLab(chiHere+i)
		END DO
		
END IF

!Send Case 2
IF((MOD(my_bounds(2),2)==0).AND.(my_rank.ne.num_cores-1)) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
!					PRINT *, 'proc', my_rank,'dumgam', i,j,k,dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		dumLab(i)=LabelLeft(my_local_dim)%vi(i)
		dumLab(chiHere+i)=LabelRight(my_local_dim)%vi(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,700+my_rank, MPI_COMM_WORLD,ierror)

CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,7000+my_rank, MPI_COMM_WORLD,ierror)
END IF


!!!END Operate Exp(-i Hodd dt/2)


!!! Operate Exp(-i Heven dt)
	IF(MOD(my_bounds(1),2)==0) THEN
	DO i=1,my_local_dim-1,2
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF	
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
	
!Case 1:If we are updating even sites and begin on an even site, send the first Gamma and second Lambda to the preceding processor
IF(my_rank.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		dumLab(i)=LabelLeft(2)%vi(i)
		dumLab(chiHere+i)=LabelRight(2)%vi(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)

CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
my_rank-1,2000+my_rank, MPI_COMM_WORLD,ierror)

END IF


!Case 2:If we are updating even sites and begin on an odd site, receive the last Gamma and next to last lambda from the preceding processor	
	ELSE
	DO i=2,my_local_dim-1,2
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF	
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
		IF(my_rank.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,800+my_rank-1, MPI_COMM_WORLD,my_status,ierror)

	CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
	my_rank-1,8000+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		LabelLeft(1)%vi(i)=dumLab(i)
		LabelRight(1)%vi(i)=dumLab(chiHere+i)
		END DO
	END IF

	
	END IF



!Receive Case 1
IF((MOD(my_bounds(2),2)==0).AND.(my_rank.ne.num_cores-1)) THEN
CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)

CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,2000+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		LabelLeft(my_local_dim+1)%vi(i)=dumLab(i)
		LabelRight(my_local_dim+1)%vi(i)=dumLab(chiHere+i)
		END DO
	
END IF


!Send Case 2
IF((MOD(my_bounds(2),2).ne.0).AND.(my_rank.ne.num_cores-1)) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		dumLab(i)=LabelLeft(my_local_dim)%vi(i)
		dumLab(chiHere+i)=LabelRight(my_local_dim)%vi(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,800+my_rank, MPI_COMM_WORLD,ierror)

CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,8000+my_rank, MPI_COMM_WORLD,ierror)
END IF

!!!END Operate Exp(-i Heven dt)


!!!BEGIN Operate Exp(-i Hodd dt/2)
	IF(MOD(my_bounds(1),2).ne.0) THEN
	DO i=1,my_local_dim-1,2
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF	
		totalTruncerr = totalTruncerr + trunctemp
	END DO
!Case 1:If we are updating odd sites and begin on an odd site, send the first Gamma and second Lambda to the preceding processor
IF(my_rank.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		dumLab(i)=LabelLeft(2)%vi(i)
		dumLab(chiHere+i)=LabelRight(2)%vi(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,300+my_rank, MPI_COMM_WORLD,ierror)

CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank-1,3000+my_rank, MPI_COMM_WORLD,ierror)
END IF

!Case 2:If we are updating odd sites and begin on an even site, receive the last Gamma and next to last lambda from the preceding processor	
	ELSE
	DO i=2,my_local_dim-1,2
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(i, Udt(i)%m, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF	
		totalTruncerr = totalTruncerr + trunctemp
	END DO
	
	IF(my_rank.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,900+my_rank-1, MPI_COMM_WORLD,my_status,ierror)

	CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
	my_rank-1,9000+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		LabelLeft(1)%vi(i)=dumLab(i)
		LabelRight(1)%vi(i)=dumLab(chiHere+i)
		END DO
				
	END IF

	END IF

!Receive Case 1
IF((MOD(my_bounds(2),2).ne.0).AND.(my_rank.ne.num_cores-1)) THEN
CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,300+my_rank+1, MPI_COMM_WORLD,my_status,ierror)

CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,3000+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
		LabelLeft(my_local_dim+1)%vi(i)=dumLab(i)
		LabelRight(my_local_dim+1)%vi(i)=dumLab(chiHere+i)
		END DO
		
END IF

!Send Case 2
IF((MOD(my_bounds(2),2)==0).AND.(my_rank.ne.num_cores-1)) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		dumLab(i)=LabelLeft(my_local_dim)%vi(i)
		dumLab(chiHere+i)=LabelRight(my_local_dim)%vi(i)
		END DO
CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,900+my_rank, MPI_COMM_WORLD,ierror)

CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,9000+my_rank, MPI_COMM_WORLD,ierror)
END IF
trunctemp=totaltruncerr
	CALL MPI_REDUCE(trunctemp,totaltruncerr,1,my_mpi_rkind,MPI_SUM,0,MPI_COMM_WORLD,ierror)

END SUBROUTINE TrotterStep2ndOrderNCParallel


SUBROUTINE CanonicalFormAllNCParallel(Gammas, Lambdas, LabelLeft, LabelRight,intDegFree)
!
!Purpose: Make all bipartite splittings canonical keeping number conservation
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
COMPLEX(KIND=rKind) :: idenMat(localSize*localSize,localSize*localSize)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
INTEGER :: i, k, j, siteInd
REAL(KIND=rKind) :: trunctemp
COMPLEX(KIND=rKind) :: dumGam(SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1))
INTEGER :: dumLab(2*SIZE(LabelLeft(1)%vi,1))
INTEGER :: gamSize, chiHere,  labSize

gamSize=SIZE(Gammas(1)%t,1)*SIZE(Gammas(1)%t,2)*SIZE(Gammas(1)%t,3)+SIZE(Gammas(1)%t,1)
chiHere=SIZE(Gammas(1)%t,1)
labSize=2*SIZE(LabelLeft(1)%vi,1)

		idenMat=CMPLX(0.0,KIND=rKind)
		DO k=1,localSize*localSize
			idenMat(k,k)=CMPLX(1.0,KIND=rKind)
		END DO


DO siteInd=0,num_cores-1

IF(my_rank==siteInd) THEN


	IF(siteInd.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
		my_rank-1,100+my_rank-1, MPI_COMM_WORLD,my_status,ierror)

	CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
		my_rank-1,1000+my_rank-1, MPI_COMM_WORLD,my_status,ierror)

			DO i=1,SIZE(Gammas(1)%t,1)
				DO j=1,SIZE(Gammas(1)%t,2)
					DO k=1,SIZE(Gammas(1)%t,3)
				Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
					END DO
				END DO
			END DO
			DO i=1,SIZE(Lambdas(1)%v,1)
			Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
			LabelLeft(1)%vi(i)=dumLab(i)
			LabelRight(1)%vi(i)=dumLab(chiHere+i)
			END DO
	END IF



		DO k=1,(my_local_dim-1)
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(k, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(k, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF
		END DO


	IF(siteInd.ne.num_cores-1) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		dumLab(i)=LabelLeft(my_local_dim)%vi(i)
		dumLab(chiHere+i)=LabelRight(my_local_dim)%vi(i)
		END DO
	CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,100+my_rank, MPI_COMM_WORLD,ierror)

	CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,1000+my_rank, MPI_COMM_WORLD,ierror)
	END IF


END IF



END DO



DO siteInd=(num_cores-1),0,(-1)

IF(my_rank==siteInd) THEN

	IF(siteInd.ne.num_cores-1) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
		my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)

	CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
		my_rank+1,2000+my_rank+1, MPI_COMM_WORLD,my_status,ierror)
			DO i=1,SIZE(Gammas(1)%t,1)
				DO j=1,SIZE(Gammas(1)%t,2)
					DO k=1,SIZE(Gammas(1)%t,3)
				Gammas(my_local_dim)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
					END DO
				END DO
			END DO
			DO i=1,SIZE(Lambdas(1)%v,1)
			Lambdas(my_local_dim+1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
			LabelLeft(my_local_dim+1)%vi(i)=dumLab(i)
			LabelRight(my_local_dim+1)%vi(i)=dumLab(chiHere+i)
			END DO
	END IF


	DO k=(my_local_dim-1),1,(-1)
IF(PRESENT(intDegFree)) THEN
			CALL TwoSiteOpNC(k, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp,intDegFree)
ELSE
			CALL TwoSiteOpNC(k, idenMat, Gammas, Lambdas, LabelLeft, LabelRight, trunctemp)
END IF

	END DO


	IF(siteInd.ne.0) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(1)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(2)%v(i)
		dumLab(i)=LabelLeft(2)%vi(i)
		dumLab(chiHere+i)=LabelRight(2)%vi(i)
		END DO
	CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)

	CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank-1,2000+my_rank, MPI_COMM_WORLD,ierror)
	END IF

	IF(siteInd.ne.num_cores-1) THEN
		DO i=1,SIZE(Gammas(1)%t,1)
			DO j=1,SIZE(Gammas(1)%t,2)
				DO k=1,SIZE(Gammas(1)%t,3)
					dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)=Gammas(my_local_dim)%t(i,j,k)
				END DO
			END DO
		END DO
		DO i=1,SIZE(Lambdas(1)%v,1)
		dumGam(chiHere*chiHere*localSize+i)=Lambdas(my_local_dim)%v(i)
		dumLab(i)=LabelLeft(my_local_dim)%vi(i)
		dumLab(chiHere+i)=LabelRight(my_local_dim)%vi(i)
		END DO
	CALL MPI_Send(dumGam,gamSize,my_mpi_cKind,&
	my_rank+1,300+my_rank, MPI_COMM_WORLD,ierror)

	CALL MPI_Send(dumLab,labSize,MPI_INTEGER,&
	my_rank+1,3000+my_rank, MPI_COMM_WORLD,ierror)
	END IF


	IF(siteInd.ne.0) THEN
	CALL MPI_Recv(dumGam,gamSize,my_mpi_cKind,&
		my_rank-1,300+my_rank-1, MPI_COMM_WORLD,my_status,ierror)

	CALL MPI_Recv(dumLab,labSize,MPI_INTEGER,&
		my_rank-1,3000+my_rank-1, MPI_COMM_WORLD,my_status,ierror)
			DO i=1,SIZE(Gammas(1)%t,1)
				DO j=1,SIZE(Gammas(1)%t,2)
					DO k=1,SIZE(Gammas(1)%t,3)
				Gammas(1)%t(i,j,k)=dumGam((i-1)*localSize*chiHere+(j-1)*chiHere+k)
					END DO
				END DO
			END DO
			DO i=1,SIZE(Lambdas(1)%v,1)
			Lambdas(1)%v(i)=REAL(dumGam(chiHere*chiHere*localSize+i),KIND=rKind)
			LabelLeft(1)%vi(i)=dumLab(i)
			LabelRight(1)%vi(i)=dumLab(chiHere+i)
			END DO
	END IF


END IF


END DO


END SUBROUTINE CanonicalFormAllNCParallel


SUBROUTINE ImagTimePropNCParallel(H, GammasOuter, LambdasOuter, LabelLeftOuter, LabelRightOuter, chiIn, intDegFree, basename)
!
!Purpose: Imaginary time propagaton algorithm-find the ground state keeping number conservaton
!OPTIONAL argument intDegFree specifies the presence of internal degrees of freedom
!and is the internal degree of freedom component whose number we wish to see (should be the ground state)
!
!See manual for more detail
!
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: chiIn
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: GammasOuter(:)
TYPE(vector), POINTER :: LambdasOuter(:)
TYPE(vectorInt), POINTER :: LabelLeftOuter(:), LabelRightOuter(:)
TYPE(tensor), POINTER :: GammasInner(:)
TYPE(vector), POINTER :: LambdasInner(:)
TYPE(vectorInt), POINTER :: LabelLeftInner(:), LabelRightInner(:)
TYPE(matrix), POINTER :: Uitp(:)
INTEGER, INTENT(IN), OPTIONAL :: intDegFree
COMPLEX(KIND=rKind) :: eye, dtodd, dteven
REAL(KIND=rKind) :: LrgstLamsBfr(systemSize+1), LrgstLamsAft(systemSize+1), GapBfrAft(systemSize+1) 
REAL(KIND=rKind) :: numberInner, cenLamBfr, cenLamAft, testTol, convCr, totalTruncerr, energy
INTEGER, DIMENSION(1) :: lrgstPoint
INTEGER :: i, j, l, chi, chiInc
REAL(KIND=8) :: my_max(2), global_max(2)
CHARACTER(len=*), OPTIONAL :: basename
INTEGER :: my_loc(1), jstart, exitstatus, chicheckout


!!! Define the time steps for ITP. 
		eye=CMPLX(0.0,1.0,KIND=rKind)
		dtodd=-eye*dtITP/2.0
		dteven=-eye*dtITP
!!! Construct the propagation operator for ITP.
		CALL AllocateOps(Uitp,my_local_dim-1,localSize*localSize)
		CALL ConstructPropagatorsParallel(H, Uitp, dtodd, dteven)

jstart=1
!Restart if possible
IF(restartSwitch) THEN
IF(.NOT.PRESENT(basename)) THEN
PRINT *,'A restart base filename must be supplied for restart option!'
STOP
END IF

IF(CheckCheckpoint(basename)=='I1') THEN
	CALL RestartNCParallel(basename,jstart,GammasOuter,LambdasOuter, LabelLeftOuter, LabelRightOuter,exitstatus)
		IF(exitstatus==1) THEN	
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
		jstart=1
		END IF
ELSE IF(CheckCheckpoint(basename)=='I2') THEN
			CALL AllocateGamLamParallel(GammasInner, LambdasInner, chiMax)
			CALL AllocateLabelParallel(LabelLeftInner, LabelRightInner, chiMax)
			CALL CopyGamLamParallel(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
			CALL CopyLabelParallel(LabelLeftInner, LabelRightInner, LabelLeftOuter, LabelRightOuter)
			CALL DeallocateGamLamParallel(GammasOuter, LambdasOuter)
			CALL DeallocateLabelParallel(LabelLeftOuter, LabelRightOuter)
	CALL RestartNCParallel(basename,jstart,GammasInner,LambdasInner, LabelLeftInner, LabelRightInner,exitstatus)
			CALL AllocateGamLamParallel(GammasOuter, LambdasOuter, chiMax)
			CALL AllocateLabelParallel(LabelLeftOuter, LabelRightOuter, chiMax)
			CALL CopyGamLamParallel(GammasOuter, LambdasOuter, GammasInner, LambdasInner)
			CALL CopyLabelParallel(LabelLeftOuter, LabelRightOuter, LabelLeftInner, LabelRightInner)
			CALL DeallocateGamLamParallel(GammasInner, LambdasInner)
			CALL DeallocateLabelParallel(LabelLeftInner, LabelRightInner)

		IF(exitstatus==1) THEN	
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
		jstart=1
		ELSE
		chiIn=chiMax
		END IF
ELSE
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
jstart=1
END IF

END IF



!!! This 'if' statement is for preventing The increment of Chi to be 0 when chiIn=chiMax.		
		IF(chiIn==chiMax) THEN
			chiInc=1
		ELSE
			chiInc=chiMax-chiIn
		END IF

!!! We conduct the ITP for both chi=chiMin and chiMax for efficiency.
!!! The ITP for chi=chiMin usually creates a good initial state for the ITP for chi=chiMax. 			
		DO i=chiIn,chiMax,chiInc

		chi=i
		PRINT *, chi
			
!!! The convergence criterion in the ITP is tighter for chi=chiMax than for chi=chiMin.
			IF((chiIn==chiMax).OR.(chi==chiMax)) THEN
				convCr=convCriterion2
			ELSE
				convCr=convCriterion1
			END IF

!!! Initialize MPS for ITP.			
			CALL AllocateGamLamParallel(GammasInner, LambdasInner, chi)
			CALL AllocateLabelParallel(LabelLeftInner, LabelRightInner, chi)
			CALL CopyGamLamParallel(GammasInner, LambdasInner, GammasOuter, LambdasOuter)
			CALL CopyLabelParallel(LabelLeftInner, LabelRightInner, LabelLeftOuter, LabelRightOuter)
			CALL DeallocateGamLamParallel(GammasOuter, LambdasOuter)
			CALL DeallocateLabelParallel(LabelLeftOuter, LabelRightOuter)

!Store the largest lambda of each splitting initially
LrgstLamsBfr=0.0_rKind
		DO l=1,my_local_dim-1
			LrgstLamsBfr(l)=LambdasInner(l)%v(1)
		END DO

IF(my_rank==num_cores-1) THEN
		DO l=my_local_dim,my_local_dim+1
			LrgstLamsBfr(l)=LambdasInner(l)%v(1)
		END DO
END IF

!!! We iterate the time step until the iteration time reaches 'maxITPsteps' 
!!! or until the convergence criterion is satisfied.
		DO j=jstart,maxITPsteps

			IF(ckptSwitch) THEN
				IF(MOD(j,stepsforckpt)==0) THEN
			IF(chi==chiMax) THEN
			CALL CheckpointNCParallel(basename,j,GammasInner,LambdasInner, LabelLeftInner, LabelRightInner, 'I2')
			ELSE
			CALL CheckpointNCParallel(basename,j,GammasInner,LambdasInner, LabelLeftInner, LabelRightInner,'I1')
			END IF
				END IF
			END IF

!!! This 'if' statement is for judging the convergence.
!!! We judge the convergence when the number of iterations j is a multiple of 'stepsForJudge'.
			IF(MOD(j,stepsForJudge)==0) THEN
			
!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN
			CALL TotalNumberParallel(numberInner, GammasInner, LambdasInner, intDegFree)
!Internal degree(s) of freedom absent			
ELSE
			CALL TotalNumberParallel(numberInner, GammasInner, LambdasInner)
END IF
			CALL TotalEnergyParallel(energy,H, GammasInner, LambdasInner)


LrgstLamsAft=0.0_rKind
!Store the largest lambda of each splitting after stepsForJudge time steps
				DO l=1,my_local_dim-1
					LrgstLamsAft(l)=LambdasInner(l)%v(1)
!Find the percent differences of each lambda
					GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
				END DO

IF(my_rank==num_cores-1) THEN
		DO l=my_local_dim,my_local_dim+1
			LrgstLamsAft(l)=LambdasInner(l)%v(1)
			GapBfrAft(l)=ABS((LrgstLamsBfr(l)-LrgstLamsAft(l))/LrgstLamsBfr(l))
		END DO
END IF

!Have all processors compute their maximum lambda difference and its location
my_loc=MAXLOC(GapBfrAft)
my_max(1)=REAL(MAXVAL(GapBfrAft), KIND=8)
my_max(2)=my_rank

!Have the master find the global max and its location
global_max=0.0_8
CALL MPI_REDUCE(my_max,global_max,2,MPI_2DOUBLE_PRECISION,&
MPI_MAXLOC,0,MPI_COMM_WORLD,ierror)

!Have the master broadcast this information to all cores
CALL MPI_BCast(global_max,2,my_mpi_rKind,0,MPI_COMM_WORLD,ierror)

testTol=global_max(1)

IF(my_rank==0) THEN
IF(PRESENT(intDegFree)) THEN
				PRINT *, 'number in the',intDegFree,'mode', numberInner, 'Energy',energy
ELSE
				PRINT *, 'number', numberInner, 'Energy',energy
END IF
END IF

IF(my_rank==FLOOR(global_max(2))) THEN
	PRINT *, 'ITP step j', j, 'lambda with largest percent difference', LambdasInner(my_loc(1))%v(1), &
				'found at position', my_loc(1), 'on core #',my_rank
	PRINT *, 'Percent difference', testTol,'convergence Criterion', convCr

END IF

				IF(testTol < convCr) EXIT
!Store the new largest lambda of each splitting
LrgstLamsBfr=0.0_rKind
		DO l=1,my_local_dim-1
			LrgstLamsBfr(l)=LrgstLamsAft(l)
		END DO

IF(my_rank==num_cores-1) THEN
		DO l=my_local_dim,my_local_dim+1
			LrgstLamsBfr(l)=LrgstLamsAft(l)
		END DO
END IF

			END IF

!Internal degree(s) of freedom present			
IF(PRESENT(intDegFree)) THEN
!Time step				
				CALL TrotterStep2ndOrderNCParallel(Uitp, GammasInner, LambdasInner, LabelLeftInner, LabelRightInner, totalTruncerr, intDegFree)
!Reorthogonalize
				CALL CanonicalFormAllNCParallel(GammasInner, LambdasInner, LabelLeftInner, LabelRightInner,intDegFree)
!Internal degree(s) of freedom absent			
ELSE
!Time step				
				CALL TrotterStep2ndOrderNCParallel(Uitp, GammasInner, LambdasInner, LabelLeftInner, LabelRightInner, totalTruncerr)
!Reorthogonalize
				CALL CanonicalFormAllNCParallel(GammasInner, LambdasInner, LabelLeftInner, LabelRightInner)
END IF


			END DO
!!! Reset GammasOuter and LambdasOuter.		
			CALL AllocateGamLamParallel(GammasOuter, LambdasOuter, chi)
			CALL AllocateLabelParallel(LabelLeftOuter, LabelRightOuter, chi)
			CALL CopyGamLamParallel(GammasOuter, LambdasOuter, GammasInner, LambdasInner)
			CALL CopyLabelParallel(LabelLeftOuter, LabelRightOuter, LabelLeftInner, LabelRightInner)
			CALL DeallocateGamLamParallel(GammasInner, LambdasInner)
			CALL DeallocateLabelParallel(LabelLeftInner, LabelRightInner)
		jstart=1
		END DO

		CALL DeallocateOps(Uitp,my_local_dim-1) 
END SUBROUTINE ImagTimePropNCParallel

SUBROUTINE SingleSiteDensityMatrixParallel(rho,Gammas,Lambdas)
!
!Purpose: Calculate the single site density matrix on every site, store in rho
!
IMPLICIT NONE
TYPE(matrix), INTENT(OUT) :: rho(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
INTEGER :: i
	DO i=1,my_local_dim-1
		CALL FormSingleSiteRho(rho(i)%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
	END DO


IF(my_rank==num_cores-1) THEN
CALL FormSingleSiteRho(rho(my_local_dim)%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)
END IF

END SUBROUTINE SingleSiteDensityMatrixParallel


SUBROUTINE TotalNumberParallel(number, Gammas, Lambdas, comPonent)
!
!Purpose: Calculate the total number for spinless site-parallel code
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
REAL(KIND=rKind) :: sendNum
INTEGER :: i
		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in TotalNumber'
			END IF
		number=0.0_rKind

!Internal degrees of freedom present
IF(PRESENT(comPonent)) THEN
		DO i=1,my_local_dim-1
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			CALL TraceMatmul(dumNum,MATMUL(Transpose(a_opS(comPonent)%mr),a_opS(comPonent)%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
		END DO

!Add in last site		
		IF(my_rank==num_cores-1) THEN
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			CALL TraceMatmul(dumNum,MATMUL(Transpose(a_opS(comPonent)%mr),a_opS(comPonent)%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
		END IF		
!Internal degrees of freedom absent
ELSE
		
		DO i=1,my_local_dim-1
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
			CALL TraceMatmul(dumNum,MATMUL(Transpose(a_op%mr),a_op%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
		END DO

!Add in last site
		IF(my_rank==num_cores-1) THEN
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			CALL TraceMatmul(dumNum,MATMUL(Transpose(a_op%mr),a_op%mr),rho%m)
				number=number+REAL(dumNum, KIND=rKind)
		END IF
END IF		

!Have all cores add the numbers together and give to master
sendNum=number	
CALL MPI_REDUCE(sendNum,number,1,my_mpi_rKind, MPI_SUM,0,MPI_COMM_WORLD,ierror)
number = number/REAL(systemSize, KIND=rKind)

		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in TotalNumber'
			END IF
END SUBROUTINE TotalNumberParallel

SUBROUTINE MeyerQmeasureParallel(Qout,Gammas, Lambdas)
!
!Purpose: Calculate the Meyer Q-measure (average local impurity)
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
REAL(KIND=rKind), INTENT(OUT) :: Qout
COMPLEX(KIND=rKind) :: dumPur
INTEGER :: i
REAL(KIND=rKind) :: totalPurity, locPurity

	ALLOCATE(rho%m(localSize, localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in MeyerQMeasure'
			END IF
	locPurity = 0.0_rKind
	! Calculate total purity, i.e., sum of tr(rho**2).
	DO i = 1, my_local_dim-1
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
		CALL TraceMatmul(dumPur,rho%m, rho%m)
		locPurity = locPurity + REAL(dumPur)
	END DO

	IF(my_rank==num_cores-1) THEN
		CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)
		CALL TraceMatmul(dumPur,rho%m, rho%m)
		locPurity = locPurity + REAL(dumPur)
	END IF

totalPurity=0.0_rKind

CALL MPI_REDUCE(locPurity,totalPurity,1,my_mpi_rKind, MPI_SUM,0,MPI_COMM_WORLD,ierror)

IF(my_rank==0) THEN
Qout = localSize*1.0_rKind/(localSize*1.0_rKind-1.0_rKind) * (1.0_rKind - totalPurity/systemSize*1.0_rKind) ! Calculate average impurity.
ELSE
Qout=0.0_rKind
END IF

	DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in MeyerQMeasure'
			END IF
END SUBROUTINE MeyerQmeasureParallel


SUBROUTINE TotalEnergyParallel(energy, H, Gammas, Lambdas)
!
!Purpose: Calculate the energy eigenvalue associated with the Hamiltonian H for site-parallel code
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: energy
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: Theta(:,:,:,:)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: l1, l2

ALLOCATE(Theta(SIZE(Gammas(1)%t,3),localSize,localSize,SIZE(Gammas(1)%t,3)))	
	energy=0.0_rKind
	
	DO l2=2,my_local_dim
		CALL ThetaKernal(Theta,Lambdas(l2)%v,Gammas(l2)%t,Lambdas(l2+1)%v)
		l1=l2-1
		CALL TwoSiteRho(rho2,Theta,Gammas(l1)%t,Gammas(l1)%t,Lambdas(l1)%v)
		CALL TraceMatmul(dumEn,H(l1)%m,rho2)
		energy=energy+REAL(dumEn, KIND=rKind)
	END DO

DEALLOCATE(Theta)
dumEn=energy	
!Compute the total energy and give it to master
CALL MPI_REDUCE(dumEn,energy,1,my_mpi_rKind, MPI_SUM,0,MPI_COMM_WORLD,ierror)
	
END SUBROUTINE TotalEnergyParallel

SUBROUTINE DefineCountVec()
!
!Purpose: Define how many local variables each core has, give it to master
!
IMPLICIT NONE
INTEGER :: i

IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
END IF

CALL MPI_Gather(my_count,1,MPI_INTEGER, countVec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

IF(my_rank==0) THEN
offsetVec(0)=0

DO i=1,num_cores-1
	offsetVec(i)=countVec(i-1)+offsetVec(i-1)
END DO
END IF

END SUBROUTINE DefineCountVec

SUBROUTINE LocalEnergyParallel(enList, H, Gammas, Lambdas)
!
!Purpose: Calculate the energy associated with each lattice bond.
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: enList(systemSize-1)
REAL(KIND=rKind), ALLOCATABLE :: myList(:)
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: Theta(:, :, :, :)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize, localSize*localSize)
COMPLEX(KIND=rKind) :: dumEn
INTEGER :: i, l1, l2

ALLOCATE(Theta(SIZE(Gammas(1)%t, 3), localSize, localSize, SIZE(Gammas(1)%t, 3)))

	enList = 0.0_rKind
my_count=my_local_dim-1

ALLOCATE(myList(my_count))

	DO l2 = 2, my_local_dim
		CALL ThetaKernal(Theta, Lambdas(l2)%v, Gammas(l2)%t, Lambdas(l2+1)%v)
		l1 = l2 - 1
		CALL TwoSiteRho(rho2, Theta, Gammas(l1)%t, Gammas(l1)%t, Lambdas(l1)%v)
		CALL TraceMatmul(dumEn,H(l1)%m, rho2)
		myList(l1)=REAL(dumEn, KIND=rKind)
		enList(l1) = myList(l1)
	END DO

DEALLOCATE(Theta)

countVec(num_cores-1)=countVec(num_cores-1)-1
offsetVec(num_cores-1)=offsetVec(num_cores-1)-1
CALL MPI_Gatherv(myList, my_count,my_mpi_rKind, enList, countVec, offsetVec,my_mpi_rKind, 0, MPI_COMM_WORLD, ierror)
countVec(num_cores-1)=countVec(num_cores-1)+1
offsetVec(num_cores-1)=offsetVec(num_cores-1)+1

DEALLOCATE(myList)

END SUBROUTINE LocalEnergyParallel


SUBROUTINE LocalEntropyDistParallel(entDist, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the real operator Op on every site, send full list to master
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: entDist(systemSize)
REAL(KIND=rKind), ALLOCATABLE :: myList(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j
INTEGER :: workSize, info
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

IF(my_rank.ne.num_cores-1) THEN
ALLOCATE(myList(my_local_dim-1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in LocalEntropyDist'
			END IF
ELSE
ALLOCATE(myList(my_local_dim), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in LocalEntropyDist'
			END IF
END IF


	DO i = 1, my_local_dim-1
		temp = 0.0_rKind
!Form single-site density matrix
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
!Diagonalize density matrix
		CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info) ! Diagonalize single site density matrix.

		DO j = 1, localSize
	IF(evals(j).ne.0.0_rKind) THEN
        temp = temp + evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
	ELSE
	temp=temp+0.0_rKind
	END IF

		END DO
	myList(i) = -1.0_rKind*temp ! Set value of entropy on-site.
	entDist(i)=myList(i)
	END DO


IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
i=my_local_dim
		temp = 0.0_rKind
!Form single-site density matrix
		CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)
!Diagonalize density matrix
		CALL ZHEEV('N', 'U', localSize, rho%m, localSize, evals, workArr, workSize, rworkArr, info) ! Diagonalize single site density matrix.

		DO j = 1, localSize
	IF(evals(j).ne.0.0_rKind) THEN
        temp = temp + evals(j)*LOG(ABS(evals(j)))/LOG(1.0_rKind*localSize) ! Compute tr(rho*log_d(rho)).
	ELSE
	temp=temp+0.0_rKind
	END IF

		END DO
	myList(i) = -1.0_rKind*temp ! Set value of entropy on-site.
	entDist(i)=myList(i)
END IF


CALL MPI_Gatherv(myList, my_count,my_mpi_rKind, entDist, countVec, offsetVec,my_mpi_rKind, 0, MPI_COMM_WORLD, ierror)

	DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in LocalEntropyDist'
			END IF

		DEALLOCATE(myList, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in LocalEntropyDist'
			END IF

	DEALLOCATE(workarr, rworkarr, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate variables in LocalEntropyDist'
			END IF

END SUBROUTINE LocalEntropyDistParallel

SUBROUTINE ChainEntropyParallel(entVec,Lambdas)
!
!Purpose: Calculate the entropy of entanglement of the chain to the left of link
!with the chain to the right on link (in the MPS approximation)
!
REAL(KIND=rKInd) :: entVec(systemSize)
TYPE(vector), POINTER :: Lambdas(:)
REAL(KIND=rKind), ALLOCATABLE :: myVec(:)
INTEGER :: chi,i,j


IF(my_rank.ne.num_cores-1) THEN
my_Count=my_local_dim-1
ALLOCATE(myVec(my_Count))
ELSE
my_Count=my_local_dim
ALLOCATE(myVec(my_Count))
END IF

chi=SIZE(Lambdas(1)%v)
myVec = 0.0_rKind
DO j=1,my_local_dim-1
	DO i=1,chi
		IF(Lambdas(j)%v(i).ne.0.0_rKind) THEN
		myVec(j)=myVec(j)-Lambdas(j)%v(i)*LOG(ABS(Lambdas(j)%v(i)))/LOG(1.0_rKind*localSize)
		END IF
	END DO
END DO


IF(my_rank==num_cores-1) THEN
j=my_local_dim
	DO i=1,chi
		IF(Lambdas(j)%v(i).ne.0.0_rKind) THEN
		myVec(j)=myVec(j)-Lambdas(j)%v(i)*LOG(ABS(Lambdas(j)%v(i)))/LOG(1.0_rKind*localSize)
		END IF
	END DO
END IF


CALL MPI_Gatherv(myVec, my_count,my_mpi_rKind, entVec, countVec, offsetVec,my_mpi_rKind, 0, MPI_COMM_WORLD, ierror)

DEALLOCATE(myVec)

END SUBROUTINE ChainEntropyParallel


!!!!!!!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE OneSiteExpValParallel!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneSiteExpValParallel_r(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the real operator Op on every site, send full list to master
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: expList(systemSize)
REAL(KIND=rKind), ALLOCATABLE :: myList(:)
REAL(KIND=rKind), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpValParallel'
			END IF

IF(my_rank.ne.num_cores-1) THEN
ALLOCATE(myList(my_local_dim-1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
ELSE
ALLOCATE(myList(my_local_dim), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
END IF

expList=0.0_rKind		
	DO i=1,my_local_dim-1,1
		rho%m=0.0_rKind
		myList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			myList(i)=REAL(TraceMatmul(Op,rho%m),KIND=rKind)
			expList(i)=myList(i)
	END DO


IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
		rho%m=0.0_rKind
		myList(my_local_dim)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			myList(my_local_dim)=REAL(TraceMatmul(Op,rho%m),KIND=rKind)
			expList(i)=myList(i)
END IF


CALL MPI_Gatherv(myList, my_count,my_mpi_rKind, expList, countVec, offsetVec,my_mpi_rKind, 0, MPI_COMM_WORLD, ierror)

	
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF

		DEALLOCATE(myList, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF


END SUBROUTINE OneSiteExpValParallel_r


SUBROUTINE OneSiteExpValParallel_rc(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the real operator Op on every site, send full list to master
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: expList(systemSize)
COMPLEX(KIND=rKind), ALLOCATABLE :: myList(:)
REAL(KIND=rKind), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpValParallel'
			END IF

IF(my_rank.ne.num_cores-1) THEN
ALLOCATE(myList(my_local_dim-1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
ELSE
ALLOCATE(myList(my_local_dim), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
END IF

expList=0.0_rKind		
	DO i=1,my_local_dim-1,1
		rho%m=0.0_rKind
		myList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			myList(i)=TraceMatmul(Op,rho%m)
			expList(i)=myList(i)
	END DO


IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
		rho%m=0.0_rKind
		myList(my_local_dim)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			myList(my_local_dim)=TraceMatmul(Op,rho%m)
			expList(i)=myList(i)
END IF


CALL MPI_Gatherv(myList, my_count,my_mpi_cKind, expList, countVec, offsetVec,my_mpi_cKind, 0, MPI_COMM_WORLD, ierror)

	
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF

		DEALLOCATE(myList, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF


END SUBROUTINE OneSiteExpValParallel_rc


SUBROUTINE OneSiteExpValParallel_cr(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the real operator Op on every site, send full list to master
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: expList(systemSize)
REAL(KIND=rKind), ALLOCATABLE :: myList(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpValParallel'
			END IF

IF(my_rank.ne.num_cores-1) THEN
ALLOCATE(myList(my_local_dim-1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
ELSE
ALLOCATE(myList(my_local_dim), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
END IF

expList=0.0_rKind		
	DO i=1,my_local_dim-1,1
		rho%m=0.0_rKind
		myList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			myList(i)=REAL(TraceMatmul(Op,rho%m))
			expList(i)=myList(i)
	END DO


IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
		rho%m=0.0_rKind
		myList(my_local_dim)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			myList(my_local_dim)=REAL(TraceMatmul(Op,rho%m))
			expList(i)=myList(i)
END IF


CALL MPI_Gatherv(myList, my_count,my_mpi_rKind, expList, countVec, offsetVec,my_mpi_rKind, 0, MPI_COMM_WORLD, ierror)

	
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF

		DEALLOCATE(myList, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF


END SUBROUTINE OneSiteExpValParallel_cr

SUBROUTINE OneSiteExpValParallel_c(expList,Op, Gammas, Lambdas)
!
!Purpose: Calculate expectation of the real operator Op on every site, send full list to master
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: expList(systemSize)
COMPLEX(KIND=rKind), ALLOCATABLE :: myList(:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpValParallel'
			END IF

IF(my_rank.ne.num_cores-1) THEN
ALLOCATE(myList(my_local_dim-1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
ELSE
ALLOCATE(myList(my_local_dim), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
END IF
expList=0.0_rKind
		
	DO i=1,my_local_dim-1,1
		rho%m=0.0_rKind
		myList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			myList(i)=TraceMatmul(Op,rho%m)
			expList(i)=myList(i)
	END DO


IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
		rho%m=0.0_rKind
		myList(my_local_dim)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			myList(my_local_dim)=TraceMatmul(Op,rho%m)
			expList(i)=myList(i)
END IF


CALL MPI_Gatherv(myList, my_count,my_mpi_cKind, expList, countVec, offsetVec,my_mpi_cKind, 0, MPI_COMM_WORLD, ierror)

				
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF

		DEALLOCATE(myList, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF


END SUBROUTINE OneSiteExpValParallel_c
!!!!!!!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE OneSiteExpValParallel!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!BEGIN CONTENTS OF INTERFACE OneSiteVarParallel!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneSiteVarParallel_r(varList,Op, Gammas, Lambdas)
!
!Purpose: Calculate variance of the real operator Op on every site
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(OUT) :: varList(systemSize)
REAL(KIND=rKind), ALLOCATABLE :: myList(:)
REAL(KIND=rKInd), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpValParallel'
			END IF

IF(my_rank.ne.num_cores-1) THEN
ALLOCATE(myList(my_local_dim-1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
ELSE
ALLOCATE(myList(my_local_dim), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
END IF

varList=0.0_rKind		
	DO i=1,my_local_dim-1,1
		rho%m=0.0_rKind
		myList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			myList(i)=REAL(TraceMatmul(MATMUL(Op,Op),rho%m),KIND=rKind)-REAL(TraceMatmul(Op,rho%m),KIND=rKind)**2
			varList(i)=myList(i)
	END DO


IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
		rho%m=0.0_rKind
		myList(my_local_dim)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			myList(my_local_dim)=REAL(TraceMatmul(MATMUL(Op,Op),rho%m),KIND=rKind)-REAL(TraceMatmul(Op,rho%m),KIND=rKind)**2
			varList(i)=myList(i)
END IF


CALL MPI_Gatherv(myList, my_count,my_mpi_rKind, varList, countVec, offsetVec,my_mpi_rKind, 0, MPI_COMM_WORLD, ierror)

	
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF

		DEALLOCATE(myList, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF

END SUBROUTINE OneSiteVarParallel_r


SUBROUTINE OneSiteVarParallel_c(varList,Op, Gammas, Lambdas)
!
!Purpose: Calculate variance of the real operator Op on every site
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: varList(systemSize)
COMPLEX(KIND=rKind), ALLOCATABLE :: myList(:)
COMPLEX(KIND=rKInd), INTENT(IN) :: Op(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix) :: rho
INTEGER :: i,j

		ALLOCATE(rho%m(localSize,localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate rho in OneSiteExpValParallel'
			END IF

IF(my_rank.ne.num_cores-1) THEN
ALLOCATE(myList(my_local_dim-1), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
ELSE
ALLOCATE(myList(my_local_dim), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF
END IF

varList=0.0_rKind		
	DO i=1,my_local_dim-1,1
		rho%m=0.0_rKind
		myList(i)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(i)%v, Gammas(i)%t, Lambdas(i+1)%v)	
			myList(i)=TraceMatmul(MATMUL(Op,Op),rho%m)-TraceMatmul(Op,rho%m)**2
			varList(i)=myList(i)
	END DO


IF(my_rank.ne.num_cores-1) THEN
my_count=my_local_dim-1
ELSE
my_count=my_local_dim
		rho%m=0.0_rKind
		myList(my_local_dim)=0.0_rKind
			CALL FormSingleSiteRho(rho%m, Lambdas(my_local_dim)%v, Gammas(my_local_dim)%t, Lambdas(my_local_dim+1)%v)	
			myList(my_local_dim)=TraceMatmul(MATMUL(Op,Op),rho%m)-TraceMatmul(Op,rho%m)**2
			varList(i)=myList(i)
END IF


CALL MPI_Gatherv(myList, my_count,my_mpi_rKind, varList, countVec, offsetVec,my_mpi_rKind, 0, MPI_COMM_WORLD, ierror)

	
		DEALLOCATE(rho%m, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to deallocate rho in OneSiteExpVal_r'
			END IF

		DEALLOCATE(myList, STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate my_List in OneSiteExpValParallel'
			END IF

END SUBROUTINE OneSiteVarParallel_c

!!!!!!!!!!!!!!!!!!!!!END CONTENTS OF INTERFACE OneSiteVarParallel!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!! BEGIN CONTENTS OF INTERFACE TwoSiteExpValParallelG !!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TwoSiteExpValParallelG_r(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
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
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:), geeSave(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: dumgee(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2,i,j,k,l, rankInt, geeSize, dumBound, outerInd
INTEGER, INTENT(IN), OPTIONAL :: phaseStat
INTEGER :: my_max(2), global_max(2), doneFlag
INTEGER ::  fullSize, counter
COMPLEX(KIND=rKind), ALLOCATABLE :: recvRho(:)

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(geeSave(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

geeSize=SIZE(Gammas(1)%t,3)*SIZE(Gammas(1)%t,3)
ALLOCATE(dumgee(geeSize))

DO l1=1,my_local_dim
!On-diagonal elements use Op1
	CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
	CALL TraceMatmul(observable(my_bounds(1)+l1-1,my_bounds(1)+l1-1),MATMUL(Op1,Op2),rho1)
END DO


!Have all processors compute number of owned sites
my_max(1)=my_local_dim
my_max(2)=my_rank

!Have the master find the global max and its location
global_max=0.0_8
CALL MPI_REDUCE(my_max,global_max,2,MPI_2INTEGER,&
MPI_MAXLOC,0,MPI_COMM_WORLD,ierror)

!Have the master broadcast this information to all cores
CALL MPI_BCast(global_max,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

DO l2=global_max(1),2,(-1)

IF(l2.le.my_local_dim) THEN
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
CALL GContraction(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(my_bounds(1)+l2-1,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

	END DO

END IF


geeSave=gee

DO outerInd=num_cores-1,1,(-1)
gee=geeSave

IF(my_rank==outerInd) THEN
dumBound=my_bounds(1)+l2-1

	IF(dumbound.gt.my_bounds(2)) THEN
	dumbound=-1
	END IF
END IF

!Have the owner broadcast this information to all cores
CALL MPI_BCast(dumBound,1,MPI_INTEGER,outerInd,MPI_COMM_WORLD,ierror)

IF(dumBound.gt.0) THEN


	DO rankInt=outerInd,1,(-1)

	IF(my_rank==rankInt) THEN


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)=gee(i,l)
		END DO
		END DO


			CALL MPI_Send(dumgee,geeSize,my_mpi_cKind,&
				my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)


	ELSE IF(my_rank==rankInt-1) THEN
			CALL MPI_Recv(dumgee,geeSize,my_mpi_cKind,&
				my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		gee(i,l)=dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)
		END DO
		END DO

			DO l1=(my_local_dim-1),1,(-1)

				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(my_bounds(1)+l1-1,dumBound),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(dumBound,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,dumBound))

				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO

	END IF



	END DO
END IF

END DO

END DO


DEALLOCATE(dumgee)


IF(my_rank.ne.num_Cores-1) THEN
geeSize=FLOOR(0.5*(my_bounds(2)-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)+1))
ELSE
geeSize=FLOOR(0.5*(my_bounds(2)+1-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)))
END IF


CALL MPI_Gather(geeSize,1,MPI_INTEGER, countvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)


IF(my_rank==0) THEN
offsetVec(0)=0

DO i=1,num_cores-1
	offsetVec(i)=countVec(i-1)+offsetVec(i-1)
END DO
END IF

ALLOCATE(dumgee(geeSize))

IF(my_rank.ne.num_Cores-1) THEN
counter=1
DO i=my_bounds(1),my_bounds(2)-1
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO
ELSE
counter=1
DO i=my_bounds(1),my_bounds(2)
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO

END IF

IF(my_rank==0) THEN
fullSize=SUM(countvec)
ALLOCATE(recvRho(fullSize))
END IF


CALL MPI_Gatherv(dumgee, geeSize,my_mpi_cKind, recvRho, countVec, offsetVec,my_mpi_cKind, 0, MPI_COMM_WORLD, ierror)

DEALLOCATE(dumgee)

IF(my_rank==0) THEN

counter=1
DO i=1,systemSize
	DO j=i,systemSize
	observable(i,j)=recvRho(counter)
	observable(j,i)=CONJG(observable(i,j))
	counter=counter+1
	END DO
END DO

DEALLOCATE(recvRho)
END IF

DEALLOCATE(gee) 
DEALLOCATE(geeSave)
DEALLOCATE(GammaP)

END SUBROUTINE TwoSiteExpValParallelG_r

SUBROUTINE TwoSiteExpValParallelG_c(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the real two-site operator Op2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
COMPLEX(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:), geeSave(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: dumgee(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2,i,j,k,l, rankInt, geeSize, dumBound, outerInd
INTEGER, INTENT(IN), OPTIONAL :: phaseStat
INTEGER :: my_max(2), global_max(2), doneFlag
INTEGER ::  fullSize, counter
COMPLEX(KIND=rKind), ALLOCATABLE :: recvRho(:)

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(geeSave(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

geeSize=SIZE(Gammas(1)%t,3)*SIZE(Gammas(1)%t,3)
ALLOCATE(dumgee(geeSize))

DO l1=1,my_local_dim
!On-diagonal elements use Op1
	CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
	CALL TraceMatmul(observable(my_bounds(1)+l1-1,my_bounds(1)+l1-1),MATMUL(Op1,Op2),rho1)
END DO


!Have all processors compute number of owned sites
my_max(1)=my_local_dim
my_max(2)=my_rank

!Have the master find the global max and its location
global_max=0.0_8
CALL MPI_REDUCE(my_max,global_max,2,MPI_2INTEGER,&
MPI_MAXLOC,0,MPI_COMM_WORLD,ierror)

!Have the master broadcast this information to all cores
CALL MPI_BCast(global_max,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

DO l2=global_max(1),2,(-1)

IF(l2.le.my_local_dim) THEN
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
CALL GContraction(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(my_bounds(1)+l2-1,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

	END DO

END IF


geeSave=gee

DO outerInd=num_cores-1,1,(-1)
gee=geeSave

IF(my_rank==outerInd) THEN
dumBound=my_bounds(1)+l2-1

	IF(dumbound.gt.my_bounds(2)) THEN
	dumbound=-1
	END IF
END IF

!Have the owner broadcast this information to all cores
CALL MPI_BCast(dumBound,1,MPI_INTEGER,outerInd,MPI_COMM_WORLD,ierror)

IF(dumBound.gt.0) THEN


	DO rankInt=outerInd,1,(-1)

	IF(my_rank==rankInt) THEN


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)=gee(i,l)
		END DO
		END DO


			CALL MPI_Send(dumgee,geeSize,my_mpi_cKind,&
				my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)


	ELSE IF(my_rank==rankInt-1) THEN
			CALL MPI_Recv(dumgee,geeSize,my_mpi_cKind,&
				my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		gee(i,l)=dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)
		END DO
		END DO

			DO l1=(my_local_dim-1),1,(-1)

				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(my_bounds(1)+l1-1,dumBound),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(dumBound,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,dumBound))

				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO

	END IF



	END DO
END IF

END DO

END DO


DEALLOCATE(dumgee)


IF(my_rank.ne.num_Cores-1) THEN
geeSize=FLOOR(0.5*(my_bounds(2)-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)+1))
ELSE
geeSize=FLOOR(0.5*(my_bounds(2)+1-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)))
END IF


CALL MPI_Gather(geeSize,1,MPI_INTEGER, countvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)


IF(my_rank==0) THEN
offsetVec(0)=0

DO i=1,num_cores-1
	offsetVec(i)=countVec(i-1)+offsetVec(i-1)
END DO
END IF

ALLOCATE(dumgee(geeSize))

IF(my_rank.ne.num_Cores-1) THEN
counter=1
DO i=my_bounds(1),my_bounds(2)-1
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO
ELSE
counter=1
DO i=my_bounds(1),my_bounds(2)
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO

END IF

IF(my_rank==0) THEN
fullSize=SUM(countvec)
ALLOCATE(recvRho(fullSize))
END IF


CALL MPI_Gatherv(dumgee, geeSize,my_mpi_cKind, recvRho, countVec, offsetVec,my_mpi_cKind, 0, MPI_COMM_WORLD, ierror)

DEALLOCATE(dumgee)

IF(my_rank==0) THEN

counter=1
DO i=1,systemSize
	DO j=i,systemSize
	observable(i,j)=recvRho(counter)
	observable(j,i)=CONJG(observable(i,j))
	counter=counter+1
	END DO
END DO

DEALLOCATE(recvRho)
END IF

DEALLOCATE(gee) 
DEALLOCATE(geeSave)
DEALLOCATE(GammaP)

END SUBROUTINE TwoSiteExpValParallelG_c

SUBROUTINE TwoSiteExpValParallelG_rc(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
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
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:), geeSave(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: dumgee(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2,i,j,k,l, rankInt, geeSize, dumBound, outerInd
INTEGER, INTENT(IN), OPTIONAL :: phaseStat
INTEGER :: my_max(2), global_max(2), doneFlag
INTEGER ::  fullSize, counter
COMPLEX(KIND=rKind), ALLOCATABLE :: recvRho(:)

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(geeSave(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

geeSize=SIZE(Gammas(1)%t,3)*SIZE(Gammas(1)%t,3)
ALLOCATE(dumgee(geeSize))

DO l1=1,my_local_dim
!On-diagonal elements use Op1
	CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
	CALL TraceMatmul(observable(my_bounds(1)+l1-1,my_bounds(1)+l1-1),MATMUL(Op1,Op2),rho1)
END DO


!Have all processors compute number of owned sites
my_max(1)=my_local_dim
my_max(2)=my_rank

!Have the master find the global max and its location
global_max=0.0_8
CALL MPI_REDUCE(my_max,global_max,2,MPI_2INTEGER,&
MPI_MAXLOC,0,MPI_COMM_WORLD,ierror)

!Have the master broadcast this information to all cores
CALL MPI_BCast(global_max,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

DO l2=global_max(1),2,(-1)

IF(l2.le.my_local_dim) THEN
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
CALL GContraction(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(my_bounds(1)+l2-1,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

	END DO

END IF


geeSave=gee

DO outerInd=num_cores-1,1,(-1)
gee=geeSave

IF(my_rank==outerInd) THEN
dumBound=my_bounds(1)+l2-1

	IF(dumbound.gt.my_bounds(2)) THEN
	dumbound=-1
	END IF
END IF

!Have the owner broadcast this information to all cores
CALL MPI_BCast(dumBound,1,MPI_INTEGER,outerInd,MPI_COMM_WORLD,ierror)

IF(dumBound.gt.0) THEN


	DO rankInt=outerInd,1,(-1)

	IF(my_rank==rankInt) THEN


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)=gee(i,l)
		END DO
		END DO


			CALL MPI_Send(dumgee,geeSize,my_mpi_cKind,&
				my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)


	ELSE IF(my_rank==rankInt-1) THEN
			CALL MPI_Recv(dumgee,geeSize,my_mpi_cKind,&
				my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		gee(i,l)=dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)
		END DO
		END DO

			DO l1=(my_local_dim-1),1,(-1)

				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(my_bounds(1)+l1-1,dumBound),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(dumBound,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,dumBound))

				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO

	END IF



	END DO
END IF

END DO

END DO


DEALLOCATE(dumgee)


IF(my_rank.ne.num_Cores-1) THEN
geeSize=FLOOR(0.5*(my_bounds(2)-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)+1))
ELSE
geeSize=FLOOR(0.5*(my_bounds(2)+1-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)))
END IF


CALL MPI_Gather(geeSize,1,MPI_INTEGER, countvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)


IF(my_rank==0) THEN
offsetVec(0)=0

DO i=1,num_cores-1
	offsetVec(i)=countVec(i-1)+offsetVec(i-1)
END DO
END IF

ALLOCATE(dumgee(geeSize))

IF(my_rank.ne.num_Cores-1) THEN
counter=1
DO i=my_bounds(1),my_bounds(2)-1
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO
ELSE
counter=1
DO i=my_bounds(1),my_bounds(2)
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO

END IF

IF(my_rank==0) THEN
fullSize=SUM(countvec)
ALLOCATE(recvRho(fullSize))
END IF


CALL MPI_Gatherv(dumgee, geeSize,my_mpi_cKind, recvRho, countVec, offsetVec,my_mpi_cKind, 0, MPI_COMM_WORLD, ierror)

DEALLOCATE(dumgee)

IF(my_rank==0) THEN

counter=1
DO i=1,systemSize
	DO j=i,systemSize
	observable(i,j)=recvRho(counter)
	observable(j,i)=CONJG(observable(i,j))
	counter=counter+1
	END DO
END DO

DEALLOCATE(recvRho)
END IF

DEALLOCATE(gee) 
DEALLOCATE(geeSave)
DEALLOCATE(GammaP)

END SUBROUTINE TwoSiteExpValParallelG_rc

SUBROUTINE TwoSiteExpValParallelG_cr(observable, Op1, Op2, Gammas, Lambdas, phaseStat)
!
!Purpose: Calculate the expectation value of the real two-site operator Op2 at every pair of sites
!
!See manual for more detail
!
IMPLICIT NONE
COMPLEX(KIND=rKind), INTENT(OUT) :: observable(:,:)
COMPLEX(KIND=rKind), INTENT(IN) :: Op1(:,:)		
REAL(KIND=rKind), INTENT(IN) :: Op2(:,:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:), geeSave(:,:)
COMPLEX(KIND=rKind), ALLOCATABLE :: dumgee(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: GammaP(:,:,:)
COMPLEX(KIND=rKind) :: rho1(localSize,localSize)
COMPLEX(KIND=rKind) :: rho2(localSize*localSize,localSize*localSize)
INTEGER :: l1, l2,i,j,k,l, rankInt, geeSize, dumBound, outerInd
INTEGER, INTENT(IN), OPTIONAL :: phaseStat
INTEGER :: my_max(2), global_max(2), doneFlag
INTEGER ::  fullSize, counter
COMPLEX(KIND=rKind), ALLOCATABLE :: recvRho(:)

ALLOCATE(gee(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(geeSave(SIZE(Gammas(1)%t,3),SIZE(Gammas(1)%t,3)))
ALLOCATE(GammaP(SIZE(Gammas(1)%t,3),localSize,SIZE(Gammas(1)%t,3)))

geeSize=SIZE(Gammas(1)%t,3)*SIZE(Gammas(1)%t,3)
ALLOCATE(dumgee(geeSize))

DO l1=1,my_local_dim
!On-diagonal elements use Op1
	CALL FormSingleSiteRho(rho1, Lambdas(l1)%v, Gammas(l1)%t, Lambdas(l1+1)%v)
	CALL TraceMatmul(observable(my_bounds(1)+l1-1,my_bounds(1)+l1-1),MATMUL(Op1,Op2),rho1)
END DO


!Have all processors compute number of owned sites
my_max(1)=my_local_dim
my_max(2)=my_rank

!Have the master find the global max and its location
global_max=0.0_8
CALL MPI_REDUCE(my_max,global_max,2,MPI_2INTEGER,&
MPI_MAXLOC,0,MPI_COMM_WORLD,ierror)

!Have the master broadcast this information to all cores
CALL MPI_BCast(global_max,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

DO l2=global_max(1),2,(-1)

IF(l2.le.my_local_dim) THEN
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
CALL GContraction(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(my_bounds(1)+l2-1,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,my_bounds(1)+l2-1))
				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

	END DO

END IF


geeSave=gee

DO outerInd=num_cores-1,1,(-1)
gee=geeSave

IF(my_rank==outerInd) THEN
dumBound=my_bounds(1)+l2-1

	IF(dumbound.gt.my_bounds(2)) THEN
	dumbound=-1
	END IF
END IF

!Have the owner broadcast this information to all cores
CALL MPI_BCast(dumBound,1,MPI_INTEGER,outerInd,MPI_COMM_WORLD,ierror)

IF(dumBound.gt.0) THEN


	DO rankInt=outerInd,1,(-1)

	IF(my_rank==rankInt) THEN


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)=gee(i,l)
		END DO
		END DO


			CALL MPI_Send(dumgee,geeSize,my_mpi_cKind,&
				my_rank-1,200+my_rank, MPI_COMM_WORLD,ierror)


	ELSE IF(my_rank==rankInt-1) THEN
			CALL MPI_Recv(dumgee,geeSize,my_mpi_cKind,&
				my_rank+1,200+my_rank+1, MPI_COMM_WORLD,my_status,ierror)


		DO i=1,SIZE(Gammas(1)%t,3)
		DO l=1,SIZE(Gammas(1)%t,3)
		gee(i,l)=dumgee((i-1)*SIZE(Gammas(1)%t,3)+l)
		END DO
		END DO

			DO l1=(my_local_dim-1),1,(-1)

				GammaP = Gammas(l1)%t
				
!Fermi Phase for final
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF

!Compute final G
CALL OneSiteOp(Op1,GammaP)
CALL GContraction(observable(my_bounds(1)+l1-1,dumBound),gee,Gammas(l1)%t,GammaP,Lambdas(l1)%v,Lambdas(l1+1)%v)
observable(dumBound,my_bounds(1)+l1-1)=CONJG(observable(my_bounds(1)+l1-1,dumBound))

				GammaP = Gammas(l1)%t
!Fermi Phase for next
				IF(PRESENT(phaseStat)) THEN
					CALL OneSiteOp(fermiPhase_op%mr,GammaP)
				END IF
!Next G
CALL GNext(gee,Gammas(l1)%t,GammaP,Lambdas(l1+1)%v)

			END DO

	END IF



	END DO
END IF

END DO

END DO


DEALLOCATE(dumgee)


IF(my_rank.ne.num_Cores-1) THEN
geeSize=FLOOR(0.5*(my_bounds(2)-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)+1))
ELSE
geeSize=FLOOR(0.5*(my_bounds(2)+1-my_bounds(1))*(2*(1+systemSize)-my_bounds(1)-my_bounds(2)))
END IF


CALL MPI_Gather(geeSize,1,MPI_INTEGER, countvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)


IF(my_rank==0) THEN
offsetVec(0)=0

DO i=1,num_cores-1
	offsetVec(i)=countVec(i-1)+offsetVec(i-1)
END DO
END IF

ALLOCATE(dumgee(geeSize))

IF(my_rank.ne.num_Cores-1) THEN
counter=1
DO i=my_bounds(1),my_bounds(2)-1
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO
ELSE
counter=1
DO i=my_bounds(1),my_bounds(2)
	DO j=i,systemSize
	dumgee(counter)=observable(i,j)
	counter=counter+1
	END DO
END DO

END IF

IF(my_rank==0) THEN
fullSize=SUM(countvec)
ALLOCATE(recvRho(fullSize))
END IF


CALL MPI_Gatherv(dumgee, geeSize,my_mpi_cKind, recvRho, countVec, offsetVec,my_mpi_cKind, 0, MPI_COMM_WORLD, ierror)

DEALLOCATE(dumgee)

IF(my_rank==0) THEN

counter=1
DO i=1,systemSize
	DO j=i,systemSize
	observable(i,j)=recvRho(counter)
	observable(j,i)=CONJG(observable(i,j))
	counter=counter+1
	END DO
END DO

DEALLOCATE(recvRho)
END IF

DEALLOCATE(gee) 
DEALLOCATE(geeSave)
DEALLOCATE(GammaP)

END SUBROUTINE TwoSiteExpValParallelG_cr



SUBROUTINE InnerProductParallel(inProd,GammasL, LambdasL, GammasR, LambdasR)
!
!Purpose: Calculate the inner product of the wavefunction by the Gammas and Lambdas.
!This routine can also be used to compute the Fidelity or Loschmidt echo if GammasL/LambdasL
!is some initial state and GammasR/LambdasR some time evolved final state
!
IMPLICIT NONE
TYPE(tensor), POINTER :: GammasL(:), GammasR(:)
TYPE(vector), POINTER :: LambdasL(:), LambdasR(:)
COMPLEX(KIND=rKind), INTENT(OUT) :: inProd
COMPLEX(KIND=rKind) :: temp
COMPLEX(KIND=rKind), ALLOCATABLE :: gee(:,:), geetemp(:,:), dumGee(:)
INTEGER :: alpha,beta,betapr,alphapr,i,l,n, geeSize, chi

chi=SIZE(GammasL(1)%t,3)
ALLOCATE(gee(SIZE(GammasL(1)%t,3),SIZE(GammasL(1)%t,3)))
ALLOCATE(geetemp(SIZE(GammasL(1)%t,3),SIZE(GammasL(1)%t,3)))

geeSize=SIZE(GammasL(1)%t,3)*SIZE(GammasL(1)%t,3)
ALLOCATE(dumgee(geeSize))

IF(my_rank==0) THEN

gee=0.0_rKind
DO alpha=1,chi
	DO beta=1,chi
		DO i=1,localSize,1
			gee(alpha,beta)=gee(alpha,beta)+LambdasL(2)%v(beta)*CONJG(GammasL(1)%t(1,i,beta))*GammasR(1)%t(1,i,alpha)*LambdasR(2)%v(alpha)
		END DO
	END DO
END DO

DO n=2,my_local_dim-1
geetemp=0.0_rKind
DO alpha=1,chi
	DO beta=1,chi
		DO betapr=1,chi
			DO alphapr=1,chi
				DO i=1,localSize,1
			geetemp(alpha,beta)=geetemp(alpha,beta)+LambdasL(n+1)%v(beta)*LambdasR(n+1)%v(alpha)&
			*CONJG(GammasL(n)%t(betapr,i,beta))*GammasR(n)%t(alphapr,i,alpha)*gee(alphapr,betapr)
				END DO
			END DO
		END DO
	END DO
END DO
gee=geetemp
END DO

		DO i=1,chi
		DO l=1,chi
		dumgee((i-1)*chi+l)=gee(i,l)
		END DO
		END DO


			CALL MPI_Send(dumgee,geeSize,my_mpi_cKind,&
				my_rank+1,200+my_rank, MPI_COMM_WORLD,ierror)

			CALL MPI_Recv(temp,1,my_mpi_cKind,&
				num_cores-1,1000, MPI_COMM_WORLD,my_status,ierror)
inProd=temp


END IF


IF((my_rank.ne.0).AND.(my_rank.ne.num_cores-1)) THEN


			CALL MPI_Recv(dumgee,geeSize,my_mpi_cKind,&
				my_rank-1,200+my_rank-1, MPI_COMM_WORLD,my_status,ierror)


		DO i=1,chi
		DO l=1,chi
		gee(i,l)=dumgee((i-1)*chi+l)
		END DO
		END DO


DO n=1,my_local_dim-1
geetemp=0.0_rKind
DO alpha=1,chi
	DO beta=1,chi
		DO betapr=1,chi
			DO alphapr=1,chi
				DO i=1,localSize,1
			geetemp(alpha,beta)=geetemp(alpha,beta)+LambdasL(n+1)%v(beta)*LambdasR(n+1)%v(alpha)&
			*CONJG(GammasL(n)%t(betapr,i,beta))*GammasR(n)%t(alphapr,i,alpha)*gee(alphapr,betapr)
				END DO
			END DO
		END DO
	END DO
END DO
gee=geetemp
END DO

		DO i=1,chi
		DO l=1,chi
		dumgee((i-1)*chi+l)=gee(i,l)
		END DO
		END DO


			CALL MPI_Send(dumgee,geeSize,my_mpi_cKind,&
				my_rank+1,200+my_rank, MPI_COMM_WORLD,ierror)
inProd=0.0_rKind
END IF


IF(my_rank==num_cores-1) THEN

			CALL MPI_Recv(dumgee,geeSize,my_mpi_cKind,&
				my_rank-1,200+my_rank-1, MPI_COMM_WORLD,my_status,ierror)


		DO i=1,chi
		DO l=1,chi
		gee(i,l)=dumgee((i-1)*chi+l)
		END DO
		END DO


DO n=1,my_local_dim-1
geetemp=0.0_rKind
DO alpha=1,chi
	DO beta=1,chi
		DO betapr=1,chi
			DO alphapr=1,chi
				DO i=1,localSize,1
			geetemp(alpha,beta)=geetemp(alpha,beta)+LambdasL(n+1)%v(beta)*LambdasR(n+1)%v(alpha)&
			*CONJG(GammasL(n)%t(betapr,i,beta))*GammasR(n)%t(alphapr,i,alpha)*gee(alphapr,betapr)
				END DO
			END DO
		END DO
	END DO
END DO
gee=geetemp
END DO


		temp=CMPLX(0.0,KIND=rKind);
		DO alpha=1,chi
			DO beta=1,chi
				DO i=1,localSize
		temp = temp+CONJG(GammasL(my_local_dim)%t(alpha,i,1))*gee(alpha,beta) &
							*GammasR(my_local_dim)%t(beta,i,1)
				END DO
			END DO
		END DO


			CALL MPI_Send(temp,1,my_mpi_cKind,&
				0,1000, MPI_COMM_WORLD,ierror)

inProd=temp
END IF


DEALLOCATE(gee) 
DEALLOCATE(dumgee) 
DEALLOCATE(geetemp) 

END SUBROUTINE InnerProductParallel


SUBROUTINE RealTimePropParallel(H, Gammas, Lambdas, chi,basename)
!
!Purpose: This subroutine performs real time propagation on the inputs Gammas and Lambdas.
! H is a list of systemSize-1 matrices of dimension localSize**2 by localSize**2.
! chi is the number of Schmidt coefficients kept at each splitting during real time propagation.
IMPLICIT NONE
INTEGER, INTENT(IN) :: chi
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(matrix), POINTER :: Urtp(:)
COMPLEX(KIND=rKind) :: dtodd, dteven
REAL(KIND=rKind) :: dt, bigDt, timeToRun
REAL(KIND=rKind) :: time, energy, number, totalTruncerr, numList(systemSize), numvarList(systemSize), bhatRealList(systemSize), bhatImagList(systemSize),  localEnergyList(systemSize-1), Q, entropyList(systemSize), chainEnt(systemSize), phaseList(systemSize)
TYPE(tensor), POINTER :: GammasInitial(:)
TYPE(vector), POINTER :: LambdasInitial(:)
COMPLEX(KIND=rKind) :: bhatList(systemSize), innerProd
TYPE(matrix) :: oneBodyDensMat, densDensCorrMat
INTEGER :: ic_ind, phase_ind, i, j, totalSteps, smallSteps, istart, exitstatus
CHARACTER(32) :: paramsName, timesName, energyName, numberName, numvarName, bhatName, phaseName, localEnergyName, QName, entropyName, OBDMName, densDensName, lambdasName, truncerrName, innerProdName, ChainEntName, progressName
CHARACTER(32) :: GammasTempName, LambdasTempName, timeTempName
CHARACTER(len=*), INTENT(IN), OPTIONAL :: basename

istart=1
!Restart if possible
IF(restartSwitch) THEN
IF(.NOT.PRESENT(basename)) THEN
PRINT *,'A restart base filename must be supplied for restart option!'
STOP
END IF

IF(CheckCheckpoint(basename)=='R') THEN
	CALL RestartParallel(basename,istart,Gammas,Lambdas, exitstatus)
		IF(exitstatus==1) THEN	
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
		istart=1
		END IF
ELSE
istart=1
END IF

END IF


!Have Master define File Names
IF(my_rank==0) THEN
CALL createFileName(paramsName,rtpDir)
CALL appendBaseName(paramsName,'',simlabel)
CALL copyName(paramsName,timesName)
CALL copyName(paramsName,energyName)
CALL copyName(paramsName,numberName)
CALL copyName(paramsName,numvarName)
CALL copyName(paramsName,bhatName)
CALL copyName(paramsName,phaseName)
CALL copyName(paramsName,localEnergyName)
CALL copyName(paramsName,QName)
CALL copyName(paramsName,entropyName)
IF(tsSwitch) THEN
CALL copyName(paramsName,densDensName)
CALL copyName(paramsName,OBDMName)
CALL appendBaseName(densDensName,'densDensCorr.dat')
CALL appendBaseName(truncerrName,'truncerr.dat')
END IF
CALL copyName(paramsName,truncerrName)
CALL copyName(paramsName,ChainEntName)
CALL copyName(paramsName,progressName)
CALL copyName(paramsName,timeTempName)

CALL appendBaseName(paramsName,'systemParams.dat')
CALL appendBaseName(timesName,'times.dat')
CALL appendBaseName(energyName,'averageEnergy.dat')
CALL appendBaseName(numberName,'averageNumber.dat')
CALL appendBaseName(numvarName,'numberVariance.dat')
CALL appendBaseName(bhatName,'orderParameter.dat')
CALL appendBaseName(phaseName,'PBphase.dat')
CALL appendBaseName(localEnergyName,'localEnergy.dat')
CALL appendBaseName(QName,'MeyerQmeasure.dat')
CALL appendBaseName(entropyName,'localEntropy.dat')
CALL appendBaseName(ChainEntName,'BlockEntropy.dat')
CALL appendBaseName(progressName,'RTPprogress.dat')
CALL appendBaseName(timeTempName,'timeTemp.dat')
CALL appendBaseName(OBDMName,'oneBodyDensMats.dat')

END IF

! Allocate single-particle density matrix, density-density correlations, and single-site density matrices.
IF(tsSwitch) THEN
	ALLOCATE(oneBodyDensMat%m(systemSize, systemSize))
	ALLOCATE(densDensCorrMat%m(systemSize, systemSize))
END IF

!Have all cores compute observables at t=0
	totalTruncerr = 0.0_8
	time=0.0_rKind
IF(istart==1) THEN
	CALL DefineCountVec()
        CALL TotalNumberParallel(number, Gammas, Lambdas)
        CALL TotalEnergyParallel(energy, H, Gammas, Lambdas)
	CALL OneSiteExpValParallel(numList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteVarParallel(numvarList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteExpValParallel(bhatList,a_op%mr, Gammas, Lambdas)
	CALL OneSiteExpValParallel(phaseList,PBphase_op%m, Gammas, Lambdas)
	CALL LocalEnergyParallel(localEnergyList, H, Gammas, Lambdas)
	CALL MeyerQMeasureParallel(Q, Gammas,Lambdas)
        CALL LocalEntropyDistParallel(entropyList, Gammas, Lambdas)
	CALL ChainEntropyParallel(chainEnt,Lambdas)
IF(tsSwitch) THEN
        CALL TwoSiteExpValParallelG(oneBodyDensMat%m,  MATMUL(TRANSPOSE(a_op%mr),a_op%mr), TensorProd(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
        CALL TwoSiteExpValParallelG(densDensCorrMat%m,  MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), TensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), Gammas, Lambdas)
END IF

!Have master write output
IF(my_rank==0) THEN
CALL openUnit(paramsName,3)
	WRITE(UNIT=3, FMT=*) num_cores,systemSize, localSize, chiRTPdnls, totNum, U0, jTunn
CLOSE(3)
CALL openUnit(timesName,5)
	WRITE(UNIT=5, FMT=*) time
CLOSE(5)
CALL openUnit(energyName,7)
	WRITE(UNIT=7, FMT=*) energy
CLOSE(7)
CALL openUnit(numberName,9)
	WRITE(UNIT=9, FMT=*) numList
CLOSE(9)
CALL openUnit(numvarName,11)
	WRITE(UNIT=11, FMT=*) numvarList
CLOSE(11)
CALL openUnit(bhatName,13)
	WRITE(UNIT=13, FMT=*) REAL(bhatList)
	WRITE(UNIT=13, FMT=*) AIMAG(bhatList)
CLOSE(13)
CALL openUnit(phaseName,15)
	WRITE(UNIT=15, FMT=*) phaseList
CLOSE(15)
CALL openUnit(localEnergyName,17)
	WRITE(UNIT=17, FMT=*) localEnergyList
CLOSE(17)
CALL openUnit(QName,19)
	WRITE(UNIT=19, FMT=*) Q
CLOSE(19)
CALL openUnit(entropyName,21)
          WRITE(UNIT=21, FMT=*) entropyList
CLOSE(21)
IF(tsSwitch) THEN
CALL openUnit(OBDMName,23)
          CALL RecordTwoSiteOb(23, oneBodyDensMat%m)
CLOSE(23)
CALL openUnit(densDensName,25)
	  CALL RecordTwoSiteOb(25, densDensCorrMat%m)
CLOSE(25)
END IF

CALL openUnit(truncerrName,27)
          WRITE(UNIT=27, FMT=*) totalTruncerr
CLOSE(27)
CALL openUnit(ChainEntName,29)
	  WRITE(UNIT=29, FMT=*) chainEnt
CLOSE(29)
CALL openUnit(progressName,31)
CLOSE(31)
CALL openUnit(timeTempName,33)
	  WRITE(UNIT=33, FMT=*) time
CLOSE(33)
END IF

END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)


! Calculate appropriate time steps based on number of data stores.
	  timeToRun = tfinalRTP - time
          bigDt = timeToRun/(1.0*stepsForStore-1) ! Compute time between data stores.
          smallSteps = CEILING(bigDt/dtRTP) ! This is the number of steps between data stores.
          dt = bigDt/(1.0*smallSteps-1.0) ! Adjust time step so that we have an integer number of time steps between stores.

IF(istart.ne.1) THEN
time=istart*(smallSteps-1)*dt
END IF

          dtodd = CMPLX(dt/2.0) ! Set odd dt.
          dteven = CMPLX(dt) ! Set even dt.

! Allocate unitaries.
          CALL AllocateOps(Urtp, my_local_dim-1, localSize*localSize)

! Construct unitary evolution operators from input Hamiltonian.  This function simply takes the matrix exponential of the Hamiltonian.
          CALL ConstructPropagatorsParallel(H, Urtp, dtodd, dteven)

! Initialize variables for SVD routine.
          CALL SVDInit(chi)

! Perform real time propagation.
IF(my_rank==0) THEN
	CALL openUnit(progressName,55,'A')
        WRITE(55, *) 'RTP is now starting!!!!'
        CLOSE(55)
	PRINT *, 'RTP is now starting!!!!'

	CALL openUnit(progressName,55,'A')
        WRITE(55, *) 'Running RTP: t =', time, ', total number =', number, ', energy =', energy, ' truncation error =', totalTruncerr
	CLOSE(55)
END IF

	IF (verboseSwitch == 1) THEN
	PRINT *, 'Running RTP: t = ', time, ', total number =', number, ', energy =', energy, ' truncation error =', totalTruncerr
	END IF

          DO i = istart, stepsForStore-1 ! Outer loop where data is stored.

			IF(ckptSwitch) THEN
				IF(MOD(i,stepsforckpt)==0) THEN
			CALL CheckpointParallel(basename,i,Gammas,Lambdas, 'R')
				END IF
			END IF


          DO j = 1, smallSteps-1 ! Inner loop where Gammas and Lambdas are being updated.
                CALL TrotterStep2ndOrderParallel(Urtp, Gammas, Lambdas, totalTruncerr) ! Update Gammas and Lambdas and compute truncation error.
                time = time + dt ! Update time variable.
             END DO

	CALL DefineCountVec()
        CALL TotalNumberParallel(number, Gammas, Lambdas)
        CALL TotalEnergyParallel(energy, H, Gammas, Lambdas)
	CALL OneSiteExpValParallel(numList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteVarParallel(numvarList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteExpValParallel(bhatList,a_op%mr, Gammas, Lambdas)
	CALL OneSiteExpValParallel(phaseList,PBphase_op%m, Gammas, Lambdas)
	CALL LocalEnergyParallel(localEnergyList, H, Gammas, Lambdas)
	CALL MeyerQMeasureParallel(Q,Gammas,Lambdas)
        CALL LocalEntropyDistParallel(entropyList, Gammas, Lambdas)
	CALL ChainEntropyParallel(chainEnt,Lambdas)
IF(tsSwitch) THEN
        CALL TwoSiteExpValParallelG(oneBodyDensMat%m,  MATMUL(TRANSPOSE(a_op%mr),a_op%mr), TensorProd(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
        CALL TwoSiteExpValParallelG(densDensCorrMat%m,  MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), TensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), Gammas, Lambdas)
END IF

!Have master write output
IF(my_rank==0) THEN
CALL openUnit(paramsName,3,'A')
	WRITE(UNIT=3, FMT=*) num_cores,systemSize, localSize, chiRTPdnls, totNum, U0, jTunn
CLOSE(3)

CALL openUnit(timesName,5,'A')
	WRITE(UNIT=5, FMT=*) time
CLOSE(5)

CALL openUnit(energyName,7,'A')
	WRITE(UNIT=7, FMT=*) energy
CLOSE(7)

CALL openUnit(numberName,9,'A')
	WRITE(UNIT=9, FMT=*) numList
CLOSE(9)
CALL openUnit(numvarName,11,'A')
	WRITE(UNIT=11, FMT=*) numvarList
CLOSE(11)
CALL openUnit(bhatName,13,'A')
	WRITE(UNIT=13, FMT=*) REAL(bhatList)
	WRITE(UNIT=13, FMT=*) AIMAG(bhatList)
CLOSE(13)
CALL openUnit(phaseName,15,'A')
	WRITE(UNIT=15, FMT=*) phaseList
CLOSE(15)
CALL openUnit(localEnergyName,17,'A')
	WRITE(UNIT=17, FMT=*) localEnergyList
CLOSE(17)
CALL openUnit(QName,19,'A')
	WRITE(UNIT=19, FMT=*) Q
CLOSE(19)
CALL openUnit(entropyName,21,'A')
          WRITE(UNIT=21, FMT=*) entropyList
CLOSE(21)
IF(tsSwitch) THEN
CALL openUnit(OBDMName,23,'A')
          CALL RecordTwoSiteOb(23, oneBodyDensMat%m)
CLOSE(23)
CALL openUnit(densDensName,25,'A')
	  CALL RecordTwoSiteOb(25, densDensCorrMat%m)
CLOSE(25)
END IF
CALL openUnit(truncerrName,27,'A')
          WRITE(UNIT=27, FMT=*) totalTruncerr
CLOSE(27)

CALL openUnit(ChainEntName,29,'A')
	  WRITE(UNIT=29, FMT=*) chainEnt
CLOSE(29)

CALL openUnit(timeTempName,33,'A')
	  WRITE(UNIT=33, FMT=*) time
CLOSE(33)
END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
	     ! Write simulation progress to screen and to file.
                IF (MOD(i+1, printPeriod) == 0) THEN
IF(my_rank==0) THEN
		   CALL Openunit(progressName,55,'A')
                   WRITE(55, *) 'Running RTP: t =', time, ', total number =', number, ' energy =', energy,  &
                        ' truncation error =', totalTruncerr
		   CLOSE(55)
		   IF (verboseSwitch == 1) THEN
                   	PRINT *, 'Running RTP: t =', time, ', total number =', number, ' energy =', energy,  &
                        ' truncation error =', totalTruncerr
		   END IF
END IF
                END IF
          END DO
IF(my_rank==0) THEN
	CALL Openunit(progressName,55,'A')

          WRITE(55, *) 'RTP is now finished!!!!'
          CLOSE(55)
END IF
	  PRINT *, 'RTP is now finished!!!!'

! Deallocate SVD variables.
          CALL SVDFinish()

! Deallocate unitaries, single-particle density matrices, single-site density matrices, and initial MPD, etc.
          CALL DeallocateOps(Urtp, my_local_dim-1)
          DEALLOCATE(oneBodyDensMat%m, densDensCorrMat%m)

END SUBROUTINE RealTimePropParallel




SUBROUTINE RealTimePropNCParallel(H, Gammas, Lambdas, LabelLeft, LabelRight, basename)
!
!Purpose: This subroutine performs real time propagation on the inputs Gammas and Lambdas.
! H is a list of systemSize-1 matrices of dimension localSize**2 by localSize**2.
! chi is the number of Schmidt coefficients kept at each splitting during real time propagation.
! initialTime is the 
IMPLICIT NONE
INTEGER :: chi
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
TYPE(matrix), POINTER :: Urtp(:)
COMPLEX(KIND=rKind) :: dtodd, dteven
REAL(KIND=rKind) :: dt, bigDt, timeToRun
REAL(KIND=rKind) :: time, energy, number, totalTruncerr, numList(systemSize), numvarList(systemSize), bhatRealList(systemSize), bhatImagList(systemSize), localEnergyList(systemSize-1), Q, entropyList(systemSize), chainEnt(systemSize), phaseList(systemSize)
TYPE(tensor), POINTER :: GammasInitial(:)
TYPE(vector), POINTER :: LambdasInitial(:)
COMPLEX(KIND=rKind) :: bhatList(systemSize), innerProd
TYPE(matrix) :: oneBodyDensMat, densDensCorrMat
INTEGER :: ic_ind, phase_ind, i, j, totalSteps, smallSteps, istart,exitstatus
CHARACTER(32) :: paramsName, timesName, energyName, numberName, numvarName, bhatName, phaseName, localEnergyName, QName, entropyName, OBDMName, densDensName, lambdasName, truncerrName, innerProdName, ChainEntName, progressName
CHARACTER(32) :: GammasTempName, LambdasTempName, timeTempName
CHARACTER(len=*), INTENT(IN), OPTIONAL :: basename

istart=1
!Restart if possible
IF(restartSwitch) THEN
IF(.NOT.PRESENT(basename)) THEN
PRINT *,'A restart base filename must be supplied for restart option!'
STOP
END IF

IF(CheckCheckpoint(basename)=='R') THEN
	CALL RestartNCParallel(basename,istart,Gammas,Lambdas,LabelLeft, LabelRight, exitstatus)
		IF(exitstatus==1) THEN	
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
		istart=1
		END IF
ELSE
istart=1
END IF

END IF

!Have Master define File Names
IF(my_rank==0) THEN
CALL createFileName(paramsName,rtpDir)
CALL appendBaseName(paramsName,'',simlabel)
CALL copyName(paramsName,timesName)
CALL copyName(paramsName,energyName)
CALL copyName(paramsName,numberName)
CALL copyName(paramsName,numvarName)
CALL copyName(paramsName,bhatName)
CALL copyName(paramsName,phaseName)
CALL copyName(paramsName,localEnergyName)
CALL copyName(paramsName,QName)
CALL copyName(paramsName,entropyName)
IF(tsSwitch) THEN
CALL copyName(paramsName,OBDMName)
CALL copyName(paramsName,densDensName)
CALL appendBaseName(densDensName,'densDensCorr.dat')
CALL appendBaseName(OBDMName,'oneBodyDensMats.dat')
END IF

CALL copyName(paramsName,truncerrName)
CALL copyName(paramsName,ChainEntName)
CALL copyName(paramsName,progressName)
CALL copyName(paramsName,timeTempName)


CALL appendBaseName(paramsName,'systemParams.dat')
CALL appendBaseName(timesName,'times.dat')
CALL appendBaseName(energyName,'averageEnergy.dat')
CALL appendBaseName(numberName,'averageNumber.dat')
CALL appendBaseName(numvarName,'numberVariance.dat')
CALL appendBaseName(bhatName,'orderParameter.dat')
CALL appendBaseName(phaseName,'PBphase.dat')
CALL appendBaseName(localEnergyName,'localEnergy.dat')
CALL appendBaseName(QName,'MeyerQmeasure.dat')
CALL appendBaseName(entropyName,'localEntropy.dat')
CALL appendBaseName(truncerrName,'truncerr.dat')
CALL appendBaseName(ChainEntName,'BlockEntropy.dat')
CALL appendBaseName(progressName,'RTPprogress.dat')
CALL appendBaseName(timeTempName,'timeTemp.dat')


END IF




! Allocate single-particle density matrix, density-density correlations, and single-site density matrices.
IF(tsSwitch) THEN
	ALLOCATE(oneBodyDensMat%m(systemSize, systemSize))
	ALLOCATE(densDensCorrMat%m(systemSize, systemSize))
END IF

!Have all cores compute observables at t=0
	totalTruncerr = 0.0_8
	time=0.0_rKind
IF(istart==1) THEN
	CALL DefineCountVec()
        CALL TotalNumberParallel(number, Gammas, Lambdas)
        CALL TotalEnergyParallel(energy, H, Gammas, Lambdas)
	CALL OneSiteExpValParallel(numList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteVarParallel(numvarList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteExpValParallel(bhatList,a_op%mr, Gammas, Lambdas)
	CALL OneSiteExpValParallel(phaseList,PBphase_op%m, Gammas, Lambdas)
        CALL LocalEntropyDistParallel(entropyList, Gammas, Lambdas)
	CALL ChainEntropyParallel(chainEnt,Lambdas)
	CALL MeyerQMeasureParallel(Q,Gammas,Lambdas)
	CALL LocalEnergyParallel(localEnergyList, H, Gammas, Lambdas)
IF(tsSwitch) THEN
        CALL TwoSiteExpValParallelG(oneBodyDensMat%m,  MATMUL(TRANSPOSE(a_op%mr),a_op%mr), TensorProd(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
        CALL TwoSiteExpValParallelG(densDensCorrMat%m,  MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), TensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), Gammas, Lambdas)
END IF


!Have master write output
IF(my_rank==0) THEN
CALL openUnit(paramsName,3)
	WRITE(UNIT=3, FMT=*) num_cores,systemSize, localSize, chiRTPdnls, totNum, U0, jTunn
CLOSE(3)

CALL openUnit(timesName,5)
	WRITE(UNIT=5, FMT=*) time
CLOSE(5)

CALL openUnit(energyName,7)
	WRITE(UNIT=7, FMT=*) energy
CLOSE(7)

CALL openUnit(numberName,9)
	WRITE(UNIT=9, FMT=*) numList
CLOSE(9)
CALL openUnit(numvarName,11)
	WRITE(UNIT=11, FMT=*) numvarList
CLOSE(11)
CALL openUnit(bhatName,13)
	WRITE(UNIT=13, FMT=*) REAL(bhatList)
	WRITE(UNIT=13, FMT=*) AIMAG(bhatList)
CLOSE(13)
CALL openUnit(phaseName,15)
	WRITE(UNIT=15, FMT=*) phaseList
CLOSE(15)
CALL openUnit(localEnergyName,17)
	WRITE(UNIT=17, FMT=*) localEnergyList
CLOSE(17)
CALL openUnit(QName,19)
	WRITE(UNIT=19, FMT=*) Q
CLOSE(19)
CALL openUnit(entropyName,21)
          WRITE(UNIT=21, FMT=*) entropyList
CLOSE(21)
IF(tsSwitch) THEN
CALL openUnit(OBDMName,23)
          CALL RecordTwoSiteOb(23, oneBodyDensMat%m)
CLOSE(23)
CALL openUnit(densDensName,25)
	  CALL RecordTwoSiteOb(25, densDensCorrMat%m)
CLOSE(25)
END IF
CALL openUnit(truncerrName,27)
          WRITE(UNIT=27, FMT=*) totalTruncerr
CLOSE(27)

CALL openUnit(ChainEntName,29)
	  WRITE(UNIT=29, FMT=*) chainEnt
CLOSE(29)

CALL openUnit(progressName,31)
CLOSE(31)

CALL openUnit(timeTempName,33)
	  WRITE(UNIT=33, FMT=*) time
CLOSE(33)
END IF
END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
! Calculate appropriate time steps based on number of data stores.
	  timeToRun = tfinalRTP - time
          bigDt = timeToRun/(1.0*stepsForStore-1) ! Compute time between data stores.
          smallSteps = CEILING(bigDt/dtRTP) ! This is the number of steps between data stores.
          dt = bigDt/(1.0*smallSteps-1) ! Adjust time step so that we have an integer number of time steps between stores.

IF(istart.ne.1) THEN
time=istart*(smallSteps-1)*dt
END IF
          dtodd = CMPLX(dt/2.0) ! Set odd dt.
          dteven = CMPLX(dt) ! Set even dt.

! Allocate unitaries.
          CALL AllocateOps(Urtp, my_local_dim-1, localSize*localSize)

! Construct unitary evolution operators from input Hamiltonian.  This function simply takes the matrix exponential of the Hamiltonian.
          CALL ConstructPropagatorsParallel(H, Urtp, dtodd, dteven)

! Compute chi of input state.
          chi = SIZE(Gammas(1)%t, 3)


! Perform real time propagation.
IF(my_rank==0) THEN
	CALL openUnit(progressName,55,'A')
        WRITE(55, *) 'RTP is now starting!!!!'
        CLOSE(55)
	PRINT *, 'RTP is now starting!!!!'

	CALL openUnit(progressName,55,'A')
        WRITE(55, *) 'Running RTP: t =', time, ', total number =', number, ', energy =', energy, ' truncation error =', totalTruncerr
	CLOSE(55)

	IF (verboseSwitch == 1) THEN
	PRINT *, 'Running RTP: t = ', time, ', total number =', number, ', energy =', energy, ' truncation error =', totalTruncerr
	END IF
END IF

          DO i = istart, stepsForStore-1 ! Outer loop where data is stored.

			IF(ckptSwitch) THEN
				IF(MOD(i,stepsforckpt)==0) THEN
		CALL CheckpointNCParallel(basename,i,Gammas,Lambdas,LabelLeft,LabelRight, 'R')
				END IF
			END IF
 
            DO j = 1, smallSteps-1 ! Inner loop where Gammas and Lambdas are being updated.
                CALL TrotterStep2ndOrderNCParallel(Urtp, Gammas, Lambdas,LabelLeft, LabelRight, totalTruncerr) ! Update Gammas and Lambdas and compute truncation error.
                time = time + dt ! Update time variable.
             END DO




	CALL DefineCountVec()
        CALL TotalNumberParallel(number, Gammas, Lambdas)
        CALL TotalEnergyParallel(energy, H, Gammas, Lambdas)
	CALL OneSiteExpValParallel(numList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteVarParallel(numvarList,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteExpValParallel(bhatList,a_op%mr, Gammas, Lambdas)
	CALL OneSiteExpValParallel(phaseList,PBphase_op%m, Gammas, Lambdas)
	CALL MeyerQMeasureParallel(Q,Gammas,Lambdas)
        CALL LocalEntropyDistParallel(entropyList, Gammas, Lambdas)
	CALL ChainEntropyParallel(chainEnt,Lambdas)
	CALL LocalEnergyParallel(localEnergyList, H, Gammas, Lambdas)
IF(tsSwitch) THEN
        CALL TwoSiteExpValParallelG(oneBodyDensMat%m,  MATMUL(TRANSPOSE(a_op%mr),a_op%mr), TensorProd(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
        CALL TwoSiteExpValParallelG(densDensCorrMat%m,  MATMUL(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), TensorProd(MATMUL(TRANSPOSE(a_op%mr),a_op%mr),MATMUL(TRANSPOSE(a_op%mr),a_op%mr)), Gammas, Lambdas)
END IF
!Have master write output
IF(my_rank==0) THEN
CALL openUnit(paramsName,3,'A')
	WRITE(UNIT=3, FMT=*) num_cores,systemSize, localSize, chiRTPdnls, totNum, U0, jTunn
CLOSE(3)

CALL openUnit(timesName,5,'A')
	WRITE(UNIT=5, FMT=*) time
CLOSE(5)

CALL openUnit(energyName,7,'A')
	WRITE(UNIT=7, FMT=*) energy
CLOSE(7)

CALL openUnit(numberName,9,'A')
	WRITE(UNIT=9, FMT=*) numList
CLOSE(9)
CALL openUnit(numvarName,11,'A')
	WRITE(UNIT=11, FMT=*) numvarList
CLOSE(11)
CALL openUnit(bhatName,13,'A')
	WRITE(UNIT=13, FMT=*) REAL(bhatList)
	WRITE(UNIT=13, FMT=*) AIMAG(bhatList)
CLOSE(13)
CALL openUnit(phaseName,15,'A')
	WRITE(UNIT=15, FMT=*) phaseList
CLOSE(15)
CALL openUnit(localEnergyName,17,'A')
	WRITE(UNIT=17, FMT=*) localEnergyList
CLOSE(17)
CALL openUnit(QName,19,'A')
	WRITE(UNIT=19, FMT=*) Q
CLOSE(19)
CALL openUnit(entropyName,21,'A')
          WRITE(UNIT=21, FMT=*) entropyList
CLOSE(21)
IF(tsSwitch) THEN
CALL openUnit(OBDMName,23,'A')
          CALL RecordTwoSiteOb(23, oneBodyDensMat%m)
CLOSE(23)
CALL openUnit(densDensName,25,'A')
	  CALL RecordTwoSiteOb(25, densDensCorrMat%m)
CLOSE(25)
END IF
CALL openUnit(truncerrName,27,'A')
          WRITE(UNIT=27, FMT=*) totalTruncerr
CLOSE(27)

CALL openUnit(ChainEntName,29,'A')
	  WRITE(UNIT=29, FMT=*) chainEnt
CLOSE(29)

CALL openUnit(progressName,31,'A')
CLOSE(31)

CALL openUnit(timeTempName,33,'A')
	  WRITE(UNIT=33, FMT=*) time
CLOSE(33)
END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

	     ! Write simulation progress to screen and to file.
                IF (MOD(i+1, printPeriod) == 0) THEN
IF(my_rank==0) THEN
		   CALL Openunit(progressName,55,'A')
                   WRITE(55, *) 'Running RTP: t =', time, ', total number =', number, ' energy =', energy,  &
                        ' truncation error =', totalTruncerr
		   CLOSE(55)
		   IF (verboseSwitch == 1) THEN
                   	PRINT *, 'Running RTP: t =', time, ', total number =', number, ' energy =', energy,  &
                        ' truncation error =', totalTruncerr
		   END IF
END IF
                END IF
          END DO
IF(my_rank==0) THEN
	CALL Openunit(progressName,55,'A')

          WRITE(55, *) 'RTP is now finished!!!!'
          CLOSE(55)
END IF
	  PRINT *, 'RTP is now finished!!!!'


! Deallocate unitaries, single-particle density matrices, single-site density matrices, and initial MPD, etc.
          CALL DeallocateOps(Urtp, my_local_dim-1)
          DEALLOCATE(oneBodyDensMat%m, densDensCorrMat%m)


END SUBROUTINE RealTimePropNCParallel


SUBROUTINE BoseHubbardQuenchParallel(rampTime)
!
!Purpose: Quench dynamics in parallel
!
!See manual for more detail
!
IMPLICIT NONE
REAL(KIND=rKind), INTENT(IN) :: rampTime
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:), Gammas0(:)
TYPE(vector), POINTER :: Lambdas(:), Lambdas0(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
TYPE(matrix), POINTER :: Urtp(:)
REAL(KIND=rKind), ALLOCATABLE :: ssOpDevs(:), ssOp(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: ssOpDevsC(:), ssOpC(:)
TYPE(matrix) :: rho
REAL(KIND=rKind) :: time, energy, number, localTruncerr, totalTruncerr,numbRamps, Qmeas, chainEnt(systemSIze)
COMPLEX(KIND=rKind) :: eye, dtodd, dteven, inProd
CHARACTER(132) :: testName,numbersName,spdmName,  Energiesname,EntropyName,ChainEntropyName,MutEntropyName
INTEGER :: i,j,l
REAL(KIND=rKind) :: start_time,end_time

IF(my_rank==0) THEN
start_time=MPI_WTIME()
END IF

!Allocate Hamiltonian
CALL AllocateOps(H,my_local_dim-1,localSize*localSize)
!Create BH operators and Hamiltonian
CALL CreateFieldOps()

CALL HamiltonianBoseHubbardParallel(H , jTunn, U0,mu0, V0)

!Allocate Gammas and Lambdas
CALL AllocateGamLamParallel(Gammas, Lambdas, chiMin)
!Numer conserving routines
IF(ncswitch) THEN
	CALL AllocateLabelParallel(LabelLeft, LabelRight, chiMin)
	!Define the initial state consistent with number conservation
	CALL InitialSetNCParallel(Gammas, Lambdas, LabelLeft, LabelRight)
	!Propagate in imaginary time
	CALL ImagTimePropNCParallel(H, Gammas, Lambdas, LabelLeft, LabelRight, chiMin)
	!Number non-conserving routines
ELSE
!Define the initial state
CALL AllStatesParallel(Gammas,Lambdas)
!Propagate in imaginary time
CALL ImagTimePropParallel(H, Gammas, Lambdas, chiMin)
END IF


IF((print_switch).and.(my_rank==0)) THEN
end_time=MPI_WTIME()
PRINT *, 'ITP finished in BoseHubbardQuenchParallel'
PRINT *, 'It took',end_time-start_time,'seconds'
start_time=MPI_WTIME()
END IF

!Copy the initial state to a set of dummy Gammas/Labdas
CALL AllocateGamLamParallel(Gammas0, Lambdas0, chiMax)
CALL CopyGamLamParallel(Gammas0, Lambdas0, Gammas, Lambdas)

!!! Define the trotter time steps for RTP. 
	dtodd = CMPLX(dtRTP/2.0) ! Set odd dt.
	dteven = CMPLX(dtRTP) ! Set even dt.

!Initialize time and cumulative truncation error
time=0.0_rKind
totalTruncerr=0.0_rKind

numbRamps=rampTime/dtRTP

!!! Construct the Real time propagator
	CALL AllocateOps(Urtp,my_local_dim-1,localSize*localSize)
	CALL ConstructPropagatorsParallel(H, Urtp, dtodd, dteven)

IF(my_rank==0) THEN
!Set up i/o
CALL createFileName(testName,rtpDir)
CALL appendBaseName(testName,'BHQuench_')
CALL appendBaseName(testName,'tR',1,rampTime)
CALL copyName(testName,numbersName)
CALL copyName(testName,Energiesname)
CALL copyName(testName,spdmName)
CALL copyName(testName,EntropyName)
CALL copyName(testName,ChainEntropyName)
CALL appendBaseName(testName,'LE.dat')
CALL openUnit(testName,100,'N')
END IF

CALL TotalEnergyParallel(energy,H, Gammas, Lambdas)
CALL TotalNumberParallel(number, Gammas, Lambdas)
CALL MeyerQMeasureParallel(Qmeas,Gammas,Lambdas)
CALL InnerProductParallel(inProd,Gammas,Lambdas,Gammas0,Lambdas0)

IF(my_rank==0) THEN
WRITE(100,*), time,rampTime, U0,number,energy, Qmeas, ABS(inProd)**2, totalTruncerr
CLOSE(100)
END IF

IF(.NOT.ncswitch) THEN
! Initialize variables for SVD routine.
	CALL SVDInit(chiMax)
END IF

!Time propagation loop
DO i=1,totalStep

IF(.NOT.ncswitch) THEN
!Trotter step one dt
	CALL TrotterStep2ndOrderParallel(Urtp, Gammas, Lambdas, localTruncerr)
	totalTruncerr=totalTruncerr+localTruncerr
ELSE
			CALL TrotterStep2ndOrderNCParallel(Urtp, Gammas, Lambdas, LabelLeft, LabelRight, localTruncerr)
	totalTruncerr=totalTruncerr+localTruncerr
END IF

!Increase time step
	time=time+dtRTP
!Ramp, redefine Hamiltonian and Propagator
IF(time.le.(0.5_rKind*rampTime)) THEN
U0=U0-18.0_rKind/numbRamps
CALL HamiltonianBoseHubbardParallel(H , jTunn, U0,mu0, V0)
CALL ConstructPropagatorsParallel(H, Urtp, dtodd, dteven)
ELSE
U0=U0+18.0_rKind/numbRamps
CALL HamiltonianBoseHubbardParallel(H , jTunn, U0,mu0, V0)
CALL ConstructPropagatorsParallel(H, Urtp, dtodd, dteven)
END IF

!Write out data	
			IF(MOD(i,stepsForStore)==0) THEN



CALL TotalEnergyParallel(energy,H, Gammas, Lambdas)
CALL TotalNumberParallel(number, Gammas, Lambdas)
CALL MeyerQMeasureParallel(Qmeas,Gammas,Lambdas)
CALL InnerProductParallel(inProd,Gammas,Lambdas,Gammas0,Lambdas0)

IF(my_rank==0) THEN
CALL openUnit(testName,100,'A')

WRITE(100,*), time,rampTime, U0,number,energy, Qmeas, ABS(inProd)**2, totalTruncerr
CLOSE(100)
END IF

IF((print_switch).and.(my_rank==0)) THEN
end_time=MPI_WTIME()
PRINT *, 'RTP step',i,' time =',time,'number is', number, ' energy is', energy
PRINT *, ' truncation error this step is', localTruncerr, 'cumulative truncation error is', totalTruncerr
PRINT *, 'Time taken so far is',end_time-start_time
END IF

			END IF
END DO

IF(.NOT.ncswitch) THEN		
! Deallocate SVD variables.
          CALL SVDFinish()
END IF

!Allocate Observables, spdm, and mutual entropy
ALLOCATE(ssOp(systemSize), ssOpDevs(systemSize))
ALLOCATE(ssOpC(systemSize), ssOpDevsC(systemSize))
ALLOCATE(rho%m(systemSize,systemSize))

IF(my_rank==0) THEN
!Allocate Observables
CALL appendBaseName(numbersName,'_Nums')
CALL appendBaseName(numbersName,'.dat')
CALL openUnit(numbersName,100,'N')
END IF
	CALL DefineCountVec()
!Compute number and deviations on site
	CALL OneSiteExpValParallel(ssOp,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
	CALL OneSiteVarParallel(ssOpDevs,MATMUL(TRANSPOSE(a_op%mr),a_op%mr), Gammas, Lambdas)
IF(my_rank==0) THEN
!Write number to file
CALL RecordOneSiteOb(100, ssOp)
!Write number deviations to file
CALL RecordOneSiteOb(100, ssOpDevs)
CLOSE(100)

!Energies file
CALL appendBaseName(EnergiesName,'_Energy')
CALL appendBaseName(EnergiesName,'.dat')
CALL openUnit(EnergiesName,101,'N')
END IF

!Compute energy if each link
CALL LocalEnergyParallel(ssOp, H, Gammas, Lambdas)
IF(my_rank==0) THEN
!Write spin projection to file
CALL RecordOneSiteOb(101, ssOp)
CLOSE(101)

!spdm file
CALL appendBaseName(spdmName,'_spdm')
CALL appendBaseName(spdmName,'.dat')
CALL openUnit(spdmName,102,'N')
END IF

        CALL TwoSiteExpValParallelG(rho%m,  TRANSPOSE(a_op%mr), a_op%mr, Gammas, Lambdas)

IF(my_rank==0) THEN
CALL RecordOp(102, rho)

CLOSE(102)

!entropy file
CALL appendBaseName(EntropyName,'_Entropy')
CALL appendBaseName(EntropyName,'.dat')
CALL openUnit(EntropyName,103,'N')

END IF
	CALL DefineCountVec()
!Compute Local vN entropy
CALL LocalEntropyDistParallel(ssOp, Gammas, Lambdas)

IF(my_rank==0) THEN
CALL RecordOneSiteOb(103, ssOp)
CLOSE(103)


!Chain entropy file
CALL appendBaseName(ChainEntropyName,'_ChainEntropy')
CALL appendBaseName(ChainEntropyName,'.dat')
CALL openUnit(ChainEntropyName,203,'N')
END IF
	CALL ChainEntropyParallel(chainEnt,Lambdas)

IF(my_rank==0) THEN
!Compute the entropy of one side of each bipartite splitting with the other
DO i=1,systemSize
WRITE(203,*) i,chainEnt(i)
END DO
CLOSE(203)

END IF



DEALLOCATE(ssOp, ssOpDevs)
DEALLOCATE(rho%m)


!Clean up
CALL DeallocateOps(Urtp,my_local_dim-1)
CALL DeallocateGamLamParallel(Gammas, Lambdas)
CALL DeallocateGamLamParallel(Gammas0, Lambdas0)
IF(ncswitch) THEN
CALL DeallocateLabelParallel(LabelLeft, LabelRight)
END IF
CALL DestroyFieldOps()
CALL DeallocateOps(H,my_local_dim-1)
		
END SUBROUTINE BoseHubbardQuenchParallel

SUBROUTINE RealTimePropRotationNCParallel(H, Gammas, Lambdas, LabelLeft, LabelRight, basename)
!
!Purpose: Propagate MHH system in real time
!
!See manual for more detail
!
IMPLICIT NONE
TYPE(matrix), POINTER :: H(:)
TYPE(tensor), POINTER :: Gammas(:), Gammas0(:)
TYPE(vector), POINTER :: Lambdas(:), Lambdas0(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
TYPE(matrix), POINTER :: Urtp(:)
TYPE(matrixReal), POINTER :: numberOp(:)
REAL(KIND=rKind) :: time, energy, numberz,numbero,numbert, localTruncerr
REAL(KIND=rKind) :: totalTruncerr, Qmeas, chainEnt(systemSize)
COMPLEX(KIND=rKind) :: eye, dtodd, dteven, inProd
TYPE(vector), POINTER ::  ssOp(:)
COMPLEX(KIND=rKind), ALLOCATABLE :: ssOpDevsC(:,:), ssOpC(:,:)
CHARACTER(132) :: globalQuantitiesName, numbersName, entropyName
INTEGER :: i,j,l, istart, exitstatus
CHARACTER(len=*), INTENT(IN), OPTIONAL :: basename

istart=1
!Restart if possible
IF(restartSwitch) THEN
IF(.NOT.PRESENT(basename)) THEN
PRINT *,'A restart base filename must be supplied for restart option!'
STOP
END IF

IF(CheckCheckpoint(basename)=='R') THEN
	CALL RestartNCParallel(basename,istart,Gammas,Lambdas,LabelLeft, LabelRight, exitstatus)
		IF(exitstatus==1) THEN	
		PRINT *,'RestartParallel Failed!'
		PRINT *,'ignoring restart request!'
		istart=1
		END IF
ELSE
istart=1
END IF

END IF

IF(my_rank==0) THEN
!Set up i/o
CALL SetupRotName(globalQuantitiesName,rtpDir)
!Copy base name to other files
CALL copyName(globalQuantitiesName,numbersName)
CALL copyName(globalQuantitiesName,entropyName)

CALL appendBaseName(numbersName,'_NumsRTP.dat')
CALL appendBaseName(globalQuantitiesName,'_GpRTP.dat')
CALL appendBaseName(entropyName,'_EntropyRTP.dat')
END IF


!Copy the gammas and lambdas for Loschmidt echo
CALL AllocateGamLamParallel(Gammas0, Lambdas0, chiMax)
CALL CopyGamLamParallel(Gammas0, Lambdas0, Gammas, Lambdas)


!!! Define the trotter time steps for RTP. 
	dtodd = CMPLX(dtRTP/2.0) ! Set odd dt.
	dteven = CMPLX(dtRTP) ! Set even dt.

!Initialize time and cumulative truncation error
time=0.0_rKind
totalTruncerr=0.0_rKind

!!! Allocate the Real time propagator
	CALL AllocateOps(Urtp,my_local_dim-1,localSize*localSize)
!Allocate Observables
ALLOCATE(ssOp(3))
ALLOCATE(numberOp(3))
DO i=1,3
ALLOCATE(numberOp(i)%mr(localSize,localSize))
numberOp(i)%mr=MATMUL(TRANSPOSE(a_opS(i)%mr),a_opS(i)%mr)
ALLOCATE(ssOp(i)%v(systemSize))
END DO

IF(istart==1) THEN
CALL DefineCountVec()
DO j=1,3,1
!Compute j^th component number and deviations on site
	CALL OneSiteExpValParallel(ssOp(j)%v,numberOp(j)%mr, Gammas, Lambdas)
END DO
CALL TotalNumberParallel(numberz, Gammas, Lambdas, 1)
CALL TotalNumberParallel(numbero, Gammas, Lambdas, 2)
CALL TotalNumberParallel(numbert, Gammas, Lambdas, 3)
CALL TotalEnergyParallel(energy, H, Gammas, Lambdas)
CALL MeyerQMeasureParallel(Qmeas,Gammas,Lambdas)
CALL InnerProductParallel(inProd,Gammas,Lambdas,Gammas0,Lambdas0)
	CALL ChainEntropyParallel(chainEnt,Lambdas)

IF(my_rank==0) THEN
CALL openUnit(numbersName,100,'N')

DO j=1,3,1
	!Write initial number to file
	CALL RecordOneSiteOb(100, ssOp(j)%v, time)
END DO
CLOSE(100)
CALL openUnit(globalQuantitiesName,101,'N')

WRITE(101,'(8E30.15)') time, totalTruncerr, energy,ABS(inProd)**2, numberz, numbero, numbert,Qmeas
CLOSE(101)

CALL openUnit(entropyName,102,'N')

WRITE(102,*) time, (chainEnt(j),j=2,systemSize)
CLOSE(102)
END IF

END IF

IF(istart.ne.1) THEN
time=istart*dtRtp
END IF

!Time propagation loop
DO i=istart,totalStep

			IF(ckptSwitch) THEN
				IF(MOD(i,stepsforckpt)==0) THEN
		CALL CheckpointNCParallel(basename,i,Gammas,Lambdas,LabelLeft,LabelRight, 'R')
				END IF
			END IF
 
	CALL HamiltonianRotationTDParallel(H,time)
	CALL ConstructPropagatorsParallel(H, Urtp, dtodd, dteven)

	CALL TrotterStep2ndOrderNCParallel(Urtp, Gammas, Lambdas, LabelLeft, LabelRight, localTruncerr, 1)
	totalTruncerr=totalTruncerr+localTruncerr
	time=time+dtRTP
!Write out data	
		IF(MOD(i,stepsForStore)==0) THEN

CALL DefineCountVec()
DO j=1,3,1
!Compute j^th component number and deviations on site
	CALL OneSiteExpValParallel(ssOp(j)%v,numberOp(j)%mr, Gammas, Lambdas)
END DO
CALL TotalNumberParallel(numberz, Gammas, Lambdas, 1)
CALL TotalNumberParallel(numbero, Gammas, Lambdas, 2)
CALL TotalNumberParallel(numbert, Gammas, Lambdas, 3)
CALL TotalEnergyParallel(energy, H, Gammas, Lambdas)
CALL MeyerQMeasureParallel(Qmeas,Gammas,Lambdas)
CALL InnerProductParallel(inProd,Gammas,Lambdas,Gammas0,Lambdas0)
	CALL ChainEntropyParallel(chainEnt,Lambdas)

IF(my_rank==0) THEN
CALL openUnit(numbersName,100,'A')
DO j=1,1+rotLevel,1
	!Write initial number to file
	CALL RecordOneSiteOb(100, ssOp(j)%v, time)
END DO
CLOSE(100)
!Global properties at t=0
CALL openUnit(globalQuantitiesName,101,'A')

WRITE(101,'(8E30.15)') time, totalTruncerr, energy,ABS(inProd)**2, numberz, numbero, numbert,Qmeas
CLOSE(101)

!Chain entropies at t=0
CALL openUnit(entropyName,102,'A')

WRITE(102,*) time, (chainEnt(j),j=2,systemSize)
CLOSE(102)

IF(print_Switch) THEN
PRINT *, 'RTP step',i,' time =',time,'number in the ',1,' component is', numberz, ' energy is', energy
PRINT *, ' truncation error this step is', localTruncerr, 'cumulative truncation error is', totalTruncerr
END IF
END IF

			END IF
		END DO

!Clean up		
CALL DeallocateOps(Urtp,my_local_dim-1)
DO i=1,3
DEALLOCATE(numberOp(i)%mr)
DEALLOCATE(ssOp(i)%v)
END DO
DEALLOCATE(ssOp)
DEALLOCATE(numberOp)

CALL DeallocateGamLamParallel(Gammas0, Lambdas0)
		
END SUBROUTINE RealTimePropRotationNCParallel


END MODULE SitesParallelSetup_Module
