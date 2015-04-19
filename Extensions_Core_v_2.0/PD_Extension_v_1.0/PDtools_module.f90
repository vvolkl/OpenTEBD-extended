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
MODULE PDtools_module

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
USE timing_module
!
! Purpose: Module to compute Bose-Hubbard and Fermi Hubbard phase diagrams
!	   in parallel for OpenSourceTEBD v2.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   9/28/09   M. L. Wall	v2.0 release
!

IMPLICIT NONE
! *** MPI Parameters
INCLUDE "mpif.h"
INTEGER :: ierror !Error flag for MPI
INTEGER :: errorcode !Error code for MPI_Abort
INTEGER :: my_rank !processor number
INTEGER :: num_cores !Total number of processors
INTEGER :: my_tag !Send/Recv tag
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: my_status !Status flag for MPI_Recv
INTEGER :: my_local_dim !Number of Gammas a processor owns
INTEGER :: my_bounds(2) !Bounds of Gammas the processor owns
INTEGER :: my_mpi_rKind, my_mpi_cKind !selected kinds for mpi

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

END SUBROUTINE Initialize_MPI

SUBROUTINE Finalize_MPI()
!
!Purpose: Free derived types and finalize MPI
!
IMPLICIT NONE

!Make sure that all processors are finished before exiting
	CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

CALL MPI_TYPE_FREE(my_mpi_rKind, ierror)
CALL MPI_TYPE_FREE(my_mpi_cKind, ierror)

CALL MPI_Finalize(ierror) !Finalize MPI
END SUBROUTINE Finalize_MPI

SUBROUTINE PhaseDiagramMaster(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
!
!Purpose: Master's Routine for computing a phase diagram in parallel
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, j_Min, j_Max
INTEGER, INTENT(IN) :: mu_Res, j_Res
REAL(KIND=rKind) :: varPass(2), dummies(2)
REAL(KIND=rKind) :: obsvec(7),obsvecSave(7)
REAL(KIND=8) :: start_time, end_time
INTEGER :: i,j,k,i1, counter, nPass, whoRank, restartCounter
INTEGER :: ioSw1, ioSw2
CHARACTER(len=132) :: paramsName, outputName,progName,jpName, tempName

!Begin timer
start_time=MPI_WTIME()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Job setup begins!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Set up i/o headings
CALL createFileName(paramsName,pdDir)

!NC switch decision 
!0 mu resolution implies number conserving code
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
	IF((FLOOR(mu_Min).ne.CEILING(mu_Min)).OR.(FLOOR(mu_Max).ne.CEILING(mu_Max))) THEN
		PRINT *, 'Warning: The number conserving method is being used but either mu_Min or mu_Max is not integer.'
		PRINT *,  'The floor is being taken!'
		PRINT *, ''
	END IF
	CALL appendBaseName(paramsName,'PDNC_')
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
	CALL appendBaseName(paramsName,'PD_')
END IF


CALL appendBaseName(paramsName,'L',systemSize)
CALL appendBaseName(paramsName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(paramsName,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(paramsName,'N',4,mu_Min)
END IF

CALL appendBaseName(paramsName,'t',4,mu_Max)
CALL appendBaseName(paramsName,'jR',4,j_Min)
CALL appendBaseName(paramsName,'t',4,j_Max)

!Copy file heading to all files
CALL copyName(paramsName,outputName)
CALL copyName(paramsName,progName)
CALL copyName(paramsName,jpName)
!For number conservation, define a temp file for postprocessing
IF(ncSwitch) THEN
CALL copyName(paramsName,tempName)
CALL appendBaseName(tempName,'output.dat')
CALL appendBaseName(outputName,'temp.dat')
ELSE
CALL appendBaseName(outputName,'output.dat')
END IF
CALL appendBaseName(paramsName,'params.dat')
CALL appendBaseName(progName,'Complete.dat')
CALL appendBaseName(jpName,'jobPool.dat')

!restartSwitch=1 restarts an existing job if possible
IF(restartSwitch==1) THEN
	!Check for an existing job
	IF(CheckName(progName)) THEN
		IF(print_Switch) THEN
			PRINT *,'Restarting phase diagram from old values!'
		END IF
	!Create a job pool from all incomplete jobs
	CALL openUnit(progName,539,'O')
	ioSw2=0
	DO WHILE(ioSw2==0)
		!Read in a completed job's pair
		READ(539,'(2E30.15)',IOSTAT=ioSw2) varPass(1),varPass(2)
		IF(ioSw2.ne.0) EXIT
		!Open job pool and scratch
		CALL openUnit(jpName,541,'O')
		OPEN(UNIT=542,STATUS='Scratch')	
		ioSw1=0
		DO WHILE(ioSw1==0)
			!Read in pairs from the job pool
			READ(541,'(2E30.15)',IOSTAT=ioSw1) dummies(1),dummies(2)
			IF(ioSw1.ne.0) EXIT
		!If the pair does not match the completed pair, write it to scratch
			IF((ABS(dummies(1)-varPass(1)).le.10.0_rKind**(-precis)).AND.&
			(ABS(dummies(2)-varPass(2)).le.10.0_rKind**(-precis))) THEN
			ELSE
			WRITE(542,'(2E30.15)') dummies(1),dummies(2)
			END IF
		END DO
		!Delete the job pool and rewind the scratch
		REWIND(542)
		CLOSE(541,STATUS='DELETE')
		!Write the scratch (possibly minus the completed job) to the job pool
		CALL openUnit(jpName,541)
		ioSw1=0
		DO WHILE(ioSw1==0)
			!Read in pairs from the scratch
			READ(542,'(2E30.15)',IOSTAT=ioSw1) dummies(1),dummies(2)
			IF(ioSw1.ne.0) EXIT
			!Write to job pool
			WRITE(541,'(2E30.15)') dummies(1),dummies(2)
		END DO
		!Delete scratch and rewind job pool
		REWIND(541)
		CLOSE(541)
		CLOSE(542)
	END DO

	CLOSE(539)

!If a checkpoint file does not exist
	ELSE
		PRINT *, 'File named '//progname//' not found!'
		PRINT *, 'Ignoring restart option!'
!If checkpoint file does not exist but a parameters file does
!(such as will happen if an identical job has been completed)
!then abort execution of the program
	IF(CheckName(paramsName)) THEN
	PRINT *, 'File named '//paramsname//' already exists!'
	PRINT *, 'Job has already been performed!  Program is exiting NOW!'
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF

	CALL openUnit(paramsName,537)

	!Write parameters to file
	WRITE(537,*) 'Name of File=',paramsName
	WRITE(537,*) 'systemSize=',systemSize
	WRITE(537,*) 'maxFilling=',maxFilling
	WRITE(537,*) 'Tunneling ranges from',j_Min, ' to ',j_Max,' with ', j_Res, ' points'
	IF(.NOT.ncSwitch) THEN
		WRITE(537,*) 'Chemical potential ranges from',mu_Min, ' to ',mu_Max,' with ', mu_Res, ' points'
	ELSE
		WRITE(537,*) 'Number ranges from',MAX(FLOOR(mu_Min)-1,0), ' to ',FLOOR(mu_Max+1),' with ',FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1, ' points'
	END IF
	WRITE(537,*) 'dtITP=',dtITP
	WRITE(537,*) 'maxITPsteps=',maxITPsteps
	WRITE(537,*) 'stepsForJudge=',stepsForJudge
	WRITE(537,*) 'minimum Chi=',chiMin
	WRITE(537,*) 'maximum Chi=',chiMax
	WRITE(537,*) 'convergence criteria=', convCriterion1,convCriterion2
	WRITE(537,*) 'number of cores', num_cores

	CLOSE(537)

	CALL openUnit(outputName,538)
	CLOSE(538)


!Open file to hold completed jobs
		CALL openUnit(progName,539,'N')
		CLOSE(539)
!Open file to hold job pool-those that are both incomplete 
!and not currently being worked on
		CALL openUnit(jpName,541,'N')
		IF(.NOT.ncSwitch) THEN
		!loop over j
			DO j=1,j_Res,1
			varPass(1)=j_Min+(j-1)*(J_max-J_min)/(j_res-1)
	
			!Loop over mu
				DO i=1,mu_res,1
				varPass(2)=mu_min+(i-1)*(mu_max-mu_min)/(mu_res-1)
				!Record (j,mu) pairs
				WRITE(541,'(2E30.15)') varPass(1), varPass(2)
				END DO
			END DO
		!Rewind and close the file
		REWIND(541)
		CLOSE(541)

		ELSE
		!loop over j
		DO j=1,j_Res,1
		varPass(1)=j_Min+(j-1)*(J_max-J_min)/(j_res-1)
			!Loop over number
			DO i=MAX(FLOOR(mu_Min)-1,0),FLOOR(mu_Max+1),1
			varPass(2)=i*1.0_rKind
			!Record (j,N) pairs
			WRITE(541,'(2E30.15)') varPass(1), varPass(2)
			END DO
		END DO
		!Rewind and close the file
		REWIND(541)
		CLOSE(541)
		!END NC IF
		END IF
	END IF

ELSE

!Even if restart not claimed, be sure that another PD with the same parameters doesn't exist
	IF(CheckName(paramsName)) THEN
	PRINT *, 'File named '//paramsname//' already exists!'
	PRINT *, 'Job has possibly already been performed!  Program is exiting NOW!'
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
	CALL openUnit(paramsName,537,'N')

	!Write parameters to file
	WRITE(537,*) 'Name of File=',paramsName
	WRITE(537,*) 'systemSize=',systemSize
	WRITE(537,*) 'maxFilling=',maxFilling
	WRITE(537,*) 'Tunneling ranges from',j_Min, ' to ',j_Max,' with ', j_Res, ' points'
	IF(.NOT.ncSwitch) THEN
		WRITE(537,*) 'Chemical potential ranges from',mu_Min, ' to ',mu_Max,' with ', mu_Res, ' points'
	ELSE
		WRITE(537,*) 'Number ranges from',MAX(FLOOR(mu_Min)-1,0), ' to ',FLOOR(mu_Max+1),' with ',FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1, ' points'
	END IF
	WRITE(537,*) 'dtITP=',dtITP
	WRITE(537,*) 'maxITPsteps=',maxITPsteps
	WRITE(537,*) 'stepsForJudge=',stepsForJudge
	WRITE(537,*) 'minimum Chi=',chiMin
	WRITE(537,*) 'maximum Chi=',chiMax
	WRITE(537,*) 'convergence criteria=', convCriterion1,convCriterion2
	WRITE(537,*) 'number of cores', num_cores

	CLOSE(537)

	CALL openUnit(outputName,538,'N')
	CLOSE(538)

!Open file to hold incomplete jobs
		CALL openUnit(progName,539,'N')
		CLOSE(539)
!Open file to hold job pool-those that are both incomplete 
!and not currently being worked on
		CALL openUnit(jpName,541,'N')

IF(.NOT.ncSwitch) THEN
!loop over j
	DO j=1,j_Res,1
	varPass(1)=j_Min+(j-1)*(J_max-J_min)/(j_res-1)

	!Loop over mu
		DO i=1,mu_res,1
		varPass(2)=mu_min+(i-1)*(mu_max-mu_min)/(mu_res-1)
	!Record (j,mu) pairs
		WRITE(541,'(2E30.15)') varPass(1), varPass(2)
		END DO
	END DO
	!Rewind and close the files
	REWIND(541)
	CLOSE(541)

ELSE
	!loop over j
	DO j=1,j_Res,1
	varPass(1)=j_Min+(j-1)*(J_max-J_min)/(j_res-1)
		!Loop over number
		DO i=MAX(FLOOR(mu_Min)-1,0),FLOOR(mu_Max+1),1
		varPass(2)=i*1.0_rKind
		!Record (j,N) pairs
		WRITE(541,'(2E30.15)') varPass(1), varPass(2)
		END DO
	END DO
	!Rewind and close the file
	REWIND(541)
	CLOSE(541)

END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Job setup ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Phase diagram computation begins!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ioSw1=0
counter=1
!Open Job pool
	CALL openUnit(jpName,541, 'O')

!Job delegation loop
DO WHILE(ioSw1==0)

!Get a (j,*) pair
	READ(541,'(2E30.15)', IOSTAT=ioSw1) varPass(1), varPass(2)
!If it was the end of the file, exit the loop
IF(ioSw1.ne.0) EXIT

!If counter<=P, send out an initial batch of jobs
	IF(counter.le.num_cores-1) THEN
		!Send a (j,*) pair to processor # counter
		CALL MPI_Send(varPass,2,my_mpi_rKind,counter,&
		1000+counter, MPI_COMM_WORLD,ierror)
	ELSE
		!Subsequently, receive the rank of any core that has finished its task
		CALL MPI_Recv(whoRank,1,MPI_INTEGER,MPI_ANY_SOURCE,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
		!Ask that core for its observables
		CALL MPI_Recv(obsvec,7,my_mpi_rKind,whoRank,whoRank,&
		MPI_COMM_WORLD,my_status,ierror)
		!Write those observables to a file
		CALL openUnit(outputName,538,'A')
		WRITE(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		CLOSE(538)

		!Record the finished pair
		CALL openUnit(progName,539,'A')
		WRITE(539,'(2E30.15)') obsvec(1),obsvec(2)
		CLOSE(539)

		!Send the core a new job
		CALL MPI_Send(varPass,2,my_mpi_rKind,whoRank,&
		1000+whoRank, MPI_COMM_WORLD,ierror)
	END IF

	counter=counter+1
END DO	

!Close the job pool
CLOSE(UNIT=541)

!Receive the last running jobs and send a stop command
DO i=1,num_cores-1
		!Receive the rank of a core who has finished their task
		CALL MPI_Recv(whoRank,1,MPI_INTEGER,MPI_ANY_SOURCE,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
		!Ask that core for its observables
		CALL MPI_Recv(obsvec,7,my_mpi_rKind,whoRank,whoRank,&
		MPI_COMM_WORLD,my_status,ierror)
		!Write those observables to a file
		CALL openUnit(outputName,538,'A')
		WRITE(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		CLOSE(538)

		!Record the finished pair
		CALL openUnit(progName,539,'A')
		WRITE(539,'(2E30.15)') obsvec(1),obsvec(2)
		CLOSE(539)

		!Define the stop command
		varPass(1)=-100.0_rKind
		varPass(2)=-200.0_rKind
		!Send the core the stop command
		CALL MPI_Send(varPass,2,my_mpi_rKind,whoRank,&
		1000+whoRank, MPI_COMM_WORLD,ierror)
counter=counter+1
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Phase diagram computation ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Postprocessing begins!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Postprocess the number conserving file to extract the chemical potential
IF(ncSwitch) THEN

!Open a new output file
CALL openUnit(tempName,540)			
CLOSE(540)

!Loop over all parameter values
!Loop over j
jloop: DO j=1,j_Res,1
varPass(1)=j_Min+(j-1)*(J_max-J_min)/(j_res-1)

!Loop over N
iloop:	DO i=MAX(FLOOR(mu_Min)-1,0),FLOOR(mu_Max+1)-1,1
	varPass(2)=i*1.0_rKind

!Open the temp file
	CALL openUnit(outputName,538,'O')
i1loop:		DO i1=1,j_Res*(FLOOR(mu_Max+1)+1-MAX(FLOOR(mu_Min)-1,0))
			!Read in data from the temp file

			READ(538, '(7E30.15)') obsvecsave(1),obsvecsave(2),obsvecsave(3),&
			obsvecsave(4),obsvecsave(5),obsvecsave(6),obsvecsave(7) 

		!If the read in data matches the parameter set of the loop, store it
		IF((ABS(obsvecSave(1)-varPass(1)).le.10.0_rKind**(-precis+1)).AND.(ABS(obsvecSave(2)-varPass(2)).le.10.0_rKind**(-precis+1))) THEN
		!Rewind the file to the beginning
				REWIND(538)

			!Loop over reading in data
kloop:			DO k=1,j_Res*(FLOOR(mu_Max+1)+1-MAX(FLOOR(mu_Min)-1,0))
		READ(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
			!If the read in matches j and has N'=N+1, keep it
			IF((ABS(obsvecSave(1)-obsvec(1)).le.10.0_rKind**(-precis+1)).AND.(ABS(obsvecSave(2)+1.0_rKind-obsvec(2)).le.10.0_rKind**(-precis+1))) THEN

			CALL openUnit(tempName,540,'A')			
!Compute the chemical potential from \mu(N+1)=E(N+1)-E(N), write data to file
			WRITE(540,'(7E30.15)') obsvec(1), obsvec(3)-obsvecSave(3), obsvec(3),obsvec(4),obsvec(5),obsvec(6),obsvec(7)
			CLOSE(540)
!Stop reading in files for this parameter set
EXIT kloop
			END IF

			
			END DO kloop
!Stop reading in  files for this parameter set
EXIT i1loop
		END IF


		END DO i1loop
!Rewind file for the next parameter set
			REWIND(538)
			CLOSE(538)

	END DO iloop

END DO jloop

!Get rid of the temporary file
CALL openUnit(outputName,538,'O')
CLOSE(UNIT=538, STATUS='DELETE')

!Postprocess the number non-conserving file to get parameters in the right order
ELSE

!Rewind the output file
CALL openUnit(outputName,538,'O')		
REWIND(538)

!Open a scratch file
OPEN(UNIT=542,STATUS='Scratch')
DO j=1,j_Res,1
varPass(1)=j_Min+(j-1)*(J_max-J_min)/(j_res-1)
	!Loop over mu
	DO i=1,mu_res,1
	varPass(2)=mu_min+(i-1)*(mu_max-mu_min)/(mu_res-1)

	ioSw1=0
		DO WHILE(ioSw1==0)
!Read in observables from output
		READ(538, '(7E30.15)',IOSTAT=ioSw1) obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		IF(ioSw1.ne.0) EXIT
!If the obseravbles match the parameters of the loop, write them to scratch
			IF((ABS(obsvec(1)-varPass(1)).le.10.0_rKind**(-precis+2)).AND.(ABS(obsvec(2)-varPass(2)).le.10.0_rKind**(-precis+2))) THEN
			WRITE(542,'(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
			obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
			EXIT
			END IF
		END DO
!Rewind observables file
		REWIND(538)
	END DO
END DO

REWIND(542)
CLOSE(538,STATUS='DELETE')
CALL openUnit(outputName,538)

!Write scratch to observables file
ioSw1=0
DO WHILE(ioSw1==0)
	READ(542, '(7E30.15)',IOSTAT=ioSw1) obsvec(1),obsvec(2),obsvec(3),&
	obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		IF(ioSw1.ne.0) EXIT
	WRITE(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
	obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
END DO

CLOSE(542)
CLOSE(538)


END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Postprocessing ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!end timer
end_time=MPI_WTIME()

!Output success message and timing data
	CALL openUnit(paramsName,537,'A')
	WRITE(537,*) 'Program has exited normally'
	WRITE(537,*) 'Total time taken=',end_time-start_time,' seconds'
	CLOSE(537)

!Get rid of the job pool and progress file
		CALL openUnit(progName,539)
		CLOSE(UNIT=539, STATUS='DELETE')

		CALL openUnit(jpName,541)
		CLOSE(UNIT=541, STATUS='DELETE')

END SUBROUTINE PhaseDiagramMaster


SUBROUTINE PhaseDiagramMasterFromU(mu_Min, mu_Max, mu_Res, U_Min, U_Max, U_Res)
!
!Purpose: Master's Routine for computing a phase diagram in parallel changing U and fixing j.
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, U_Min, U_Max
INTEGER, INTENT(IN) :: mu_Res, U_Res
REAL(KIND=rKind) :: varPass(2), dummies(2)
REAL(KIND=rKind) :: obsvec(7),obsvecSave(7)
REAL(KIND=8) :: start_time, end_time
INTEGER :: i,j,k,i1, counter, nPass, whoRank, restartCounter
INTEGER :: ioSw1, ioSw2
CHARACTER(len=132) :: paramsName, outputName,progName,jpName, tempName

!Begin timer
start_time=MPI_WTIME()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Job setup begins!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Set up i/o headings
CALL createFileName(paramsName,pdDir)

!NC switch decision 
!0 mu resolution implies number conserving code
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
	IF((FLOOR(mu_Min).ne.CEILING(mu_Min)).OR.(FLOOR(mu_Max).ne.CEILING(mu_Max))) THEN
		PRINT *, 'Warning: The number conserving method is being used but either mu_Min or mu_Max is not integer.'
		PRINT *,  'The floor is being taken!'
		PRINT *, ''
	END IF
	CALL appendBaseName(paramsName,'PDNC_')
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
	CALL appendBaseName(paramsName,'PD_')
END IF


CALL appendBaseName(paramsName,'L',systemSize)
CALL appendBaseName(paramsName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(paramsName,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(paramsName,'N',4,mu_Min)
END IF

CALL appendBaseName(paramsName,'t',4,mu_Max)
CALL appendBaseName(paramsName,'uR',4,U_Min)
CALL appendBaseName(paramsName,'t',4,U_Max)

!Copy file heading to all files
CALL copyName(paramsName,outputName)
CALL copyName(paramsName,progName)
CALL copyName(paramsName,jpName)
!For number conservation, define a temp file for postprocessing
IF(ncSwitch) THEN
CALL copyName(paramsName,tempName)
CALL appendBaseName(tempName,'output.dat')
CALL appendBaseName(outputName,'temp.dat')
ELSE
CALL appendBaseName(outputName,'output.dat')
END IF
CALL appendBaseName(paramsName,'params.dat')
CALL appendBaseName(progName,'Complete.dat')
CALL appendBaseName(jpName,'jobPool.dat')

!restartSwitch=1 restarts an existing job if possible
IF(restartSwitch==1) THEN
	!Check for an existing job
	IF(CheckName(progName)) THEN
		IF(print_Switch) THEN
			PRINT *,'Restarting phase diagram from old values!'
		END IF
	!Create a job pool from all incomplete jobs
	CALL openUnit(progName,539,'O')
	ioSw2=0
	DO WHILE(ioSw2==0)
		!Read in a completed job's pair
		READ(539,'(2E30.15)',IOSTAT=ioSw2) varPass(1),varPass(2)
		IF(ioSw2.ne.0) EXIT
		!Open job pool and scratch
		CALL openUnit(jpName,541,'O')
		OPEN(UNIT=542,STATUS='Scratch')	
		ioSw1=0
		DO WHILE(ioSw1==0)
			!Read in pairs from the job pool
			READ(541,'(2E30.15)',IOSTAT=ioSw1) dummies(1),dummies(2)
			IF(ioSw1.ne.0) EXIT
		!If the pair does not match the completed pair, write it to scratch
			IF((ABS(dummies(1)-varPass(1)).le.10.0_rKind**(-precis)).AND.&
			(ABS(dummies(2)-varPass(2)).le.10.0_rKind**(-precis))) THEN
			ELSE
			WRITE(542,'(2E30.15)') dummies(1),dummies(2)
			END IF
		END DO
		!Delete the job pool and rewind the scratch
		REWIND(542)
		CLOSE(541,STATUS='DELETE')
		!Write the scratch (possibly minus the completed job) to the job pool
		CALL openUnit(jpName,541)
		ioSw1=0
		DO WHILE(ioSw1==0)
			!Read in pairs from the scratch
			READ(542,'(2E30.15)',IOSTAT=ioSw1) dummies(1),dummies(2)
			IF(ioSw1.ne.0) EXIT
			!Write to job pool
			WRITE(541,'(2E30.15)') dummies(1),dummies(2)
		END DO
		!Delete scratch and rewind job pool
		REWIND(541)
		CLOSE(541)
		CLOSE(542)
	END DO

	CLOSE(539)

!If a checkpoint file does not exist
	ELSE
		PRINT *, 'File named '//progname//' not found!'
		PRINT *, 'Ignoring restart option!'
!If checkpoint file does not exist but a parameters file does
!(such as will happen if an identical job has been completed)
!then abort execution of the program
	IF(CheckName(paramsName)) THEN
	PRINT *, 'File named '//paramsname//' already exists!'
	PRINT *, 'Job has already been performed!  Program is exiting NOW!'
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF

	CALL openUnit(paramsName,537)

	!Write parameters to file
	WRITE(537,*) 'Name of File=',paramsName
	WRITE(537,*) 'systemSize=',systemSize
	WRITE(537,*) 'maxFilling=',maxFilling
	WRITE(537,*) 'U ranges from',U_Min, ' to ',U_Max,' with ', U_Res, ' points'
	IF(.NOT.ncSwitch) THEN
		WRITE(537,*) 'Chemical potential ranges from',mu_Min, ' to ',mu_Max,' with ', mu_Res, ' points'
	ELSE
		WRITE(537,*) 'Number ranges from',MAX(FLOOR(mu_Min)-1,0), ' to ',FLOOR(mu_Max+1),' with ',FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1, ' points'
	END IF
	WRITE(537,*) 'dtITP=',dtITP
	WRITE(537,*) 'maxITPsteps=',maxITPsteps
	WRITE(537,*) 'stepsForJudge=',stepsForJudge
	WRITE(537,*) 'minimum Chi=',chiMin
	WRITE(537,*) 'maximum Chi=',chiMax
	WRITE(537,*) 'convergence criteria=', convCriterion1,convCriterion2
	WRITE(537,*) 'number of cores', num_cores

	CLOSE(537)

	CALL openUnit(outputName,538)
	CLOSE(538)


!Open file to hold completed jobs
		CALL openUnit(progName,539,'N')
		CLOSE(539)
!Open file to hold job pool-those that are both incomplete 
!and not currently being worked on
		CALL openUnit(jpName,541,'N')
		IF(.NOT.ncSwitch) THEN
		!loop over j
			DO j=1,U_Res,1
			varPass(1)=U_Min+(j-1)*(U_max-U_min)/(U_res-1)
	
			!Loop over mu
				DO i=1,mu_res,1
				varPass(2)=mu_min+(i-1)*(mu_max-mu_min)/(mu_res-1)
				!Record (j,mu) pairs
				WRITE(541,'(2E30.15)') varPass(1), varPass(2)
				END DO
			END DO
		!Rewind and close the file
		REWIND(541)
		CLOSE(541)

		ELSE
		!loop over j
		DO j=1,U_Res,1
		varPass(1)=U_Min+(j-1)*(U_max-U_min)/(U_res-1)
			!Loop over number
			DO i=MAX(FLOOR(mu_Min)-1,0),FLOOR(mu_Max+1),1
			varPass(2)=i*1.0_rKind
			!Record (j,N) pairs
			WRITE(541,'(2E30.15)') varPass(1), varPass(2)
			END DO
		END DO
		!Rewind and close the file
		REWIND(541)
		CLOSE(541)
		!END NC IF
		END IF
	END IF

ELSE

!Even if restart not claimed, be sure that another PD with the same parameters doesn't exist
	IF(CheckName(paramsName)) THEN
	PRINT *, 'File named '//paramsname//' already exists!'
	PRINT *, 'Job has possibly already been performed!  Program is exiting NOW!'
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
	CALL openUnit(paramsName,537,'N')

	!Write parameters to file
	WRITE(537,*) 'Name of File=',paramsName
	WRITE(537,*) 'systemSize=',systemSize
	WRITE(537,*) 'maxFilling=',maxFilling
	WRITE(537,*) 'U ranges from',U_Min, ' to ',U_Max,' with ', U_Res, ' points'
	IF(.NOT.ncSwitch) THEN
		WRITE(537,*) 'Chemical potential ranges from',mu_Min, ' to ',mu_Max,' with ', mu_Res, ' points'
	ELSE
		WRITE(537,*) 'Number ranges from',MAX(FLOOR(mu_Min)-1,0), ' to ',FLOOR(mu_Max+1),' with ',FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1, ' points'
	END IF
	WRITE(537,*) 'dtITP=',dtITP
	WRITE(537,*) 'maxITPsteps=',maxITPsteps
	WRITE(537,*) 'stepsForJudge=',stepsForJudge
	WRITE(537,*) 'minimum Chi=',chiMin
	WRITE(537,*) 'maximum Chi=',chiMax
	WRITE(537,*) 'convergence criteria=', convCriterion1,convCriterion2
	WRITE(537,*) 'number of cores', num_cores

	CLOSE(537)

	CALL openUnit(outputName,538,'N')
	CLOSE(538)

!Open file to hold incomplete jobs
		CALL openUnit(progName,539,'N')
		CLOSE(539)
!Open file to hold job pool-those that are both incomplete 
!and not currently being worked on
		CALL openUnit(jpName,541,'N')

IF(.NOT.ncSwitch) THEN
!loop over U
	DO j=1,U_Res,1
	varPass(1)=U_Min+(j-1)*(U_max-U_min)/(U_res-1)

	!Loop over mu
		DO i=1,mu_res,1
		varPass(2)=mu_min+(i-1)*(mu_max-mu_min)/(mu_res-1)
	!Record (j,mu) pairs
		WRITE(541,'(2E30.15)') varPass(1), varPass(2)
		END DO
	END DO
	!Rewind and close the files
	REWIND(541)
	CLOSE(541)

ELSE
	!loop over j
	DO j=1,U_Res,1
	varPass(1)=U_Min+(j-1)*(U_max-U_min)/(U_res-1)
		!Loop over number
		DO i=MAX(FLOOR(mu_Min)-1,0),FLOOR(mu_Max+1),1
		varPass(2)=i*1.0_rKind
		!Record (j,N) pairs
		WRITE(541,'(2E30.15)') varPass(1), varPass(2)
		END DO
	END DO
	!Rewind and close the file
	REWIND(541)
	CLOSE(541)

END IF

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Job setup ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Phase diagram computation begins!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ioSw1=0
counter=1
!Open Job pool
	CALL openUnit(jpName,541, 'O')

!Job delegation loop
DO WHILE(ioSw1==0)

!Get a (U,*) pair
	READ(541,'(2E30.15)', IOSTAT=ioSw1) varPass(1), varPass(2)
!If it was the end of the file, exit the loop
IF(ioSw1.ne.0) EXIT

!If counter<=P, send out an initial batch of jobs
	IF(counter.le.num_cores-1) THEN
		!Send a (U,*) pair to processor # counter
		CALL MPI_Send(varPass,2,my_mpi_rKind,counter,&
		1000+counter, MPI_COMM_WORLD,ierror)
	ELSE
		!Subsequently, receive the rank of any core that has finished its task
		CALL MPI_Recv(whoRank,1,MPI_INTEGER,MPI_ANY_SOURCE,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
		!Ask that core for its observables
		CALL MPI_Recv(obsvec,7,my_mpi_rKind,whoRank,whoRank,&
		MPI_COMM_WORLD,my_status,ierror)
		!Write those observables to a file
		CALL openUnit(outputName,538,'A')
		WRITE(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		CLOSE(538)

		!Record the finished pair
		CALL openUnit(progName,539,'A')
		WRITE(539,'(2E30.15)') obsvec(1),obsvec(2)
		CLOSE(539)

		!Send the core a new job
		CALL MPI_Send(varPass,2,my_mpi_rKind,whoRank,&
		1000+whoRank, MPI_COMM_WORLD,ierror)
	END IF

	counter=counter+1
END DO	

!Close the job pool
CLOSE(UNIT=541)

!Receive the last running jobs and send a stop command
DO i=1,num_cores-1
		!Receive the rank of a core who has finished their task
		CALL MPI_Recv(whoRank,1,MPI_INTEGER,MPI_ANY_SOURCE,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
		!Ask that core for its observables
		CALL MPI_Recv(obsvec,7,my_mpi_rKind,whoRank,whoRank,&
		MPI_COMM_WORLD,my_status,ierror)
		!Write those observables to a file
		CALL openUnit(outputName,538,'A')
		WRITE(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		CLOSE(538)

		!Record the finished pair
		CALL openUnit(progName,539,'A')
		WRITE(539,'(2E30.15)') obsvec(1),obsvec(2)
		CLOSE(539)

		!Define the stop command
		varPass(1)=-100.0_rKind
		varPass(2)=-200.0_rKind
		!Send the core the stop command
		CALL MPI_Send(varPass,2,my_mpi_rKind,whoRank,&
		1000+whoRank, MPI_COMM_WORLD,ierror)
counter=counter+1
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Phase diagram computation ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Postprocessing begins!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Postprocess the number conserving file to extract the chemical potential
IF(ncSwitch) THEN

!Open a new output file
CALL openUnit(tempName,540)			
CLOSE(540)

!Loop over all parameter values
!Loop over U
jloop: DO j=1,U_Res,1
varPass(1)=U_Min+(j-1)*(U_max-U_min)/(U_res-1)

!Loop over N
iloop:	DO i=MAX(FLOOR(mu_Min)-1,0),FLOOR(mu_Max+1)-1,1
	varPass(2)=i*1.0_rKind

!Open the temp file
	CALL openUnit(outputName,538,'O')
i1loop:		DO i1=1,U_Res*(FLOOR(mu_Max+1)+1-MAX(FLOOR(mu_Min)-1,0))
			!Read in data from the temp file

			READ(538, '(7E30.15)') obsvecsave(1),obsvecsave(2),obsvecsave(3),&
			obsvecsave(4),obsvecsave(5),obsvecsave(6),obsvecsave(7) 

		!If the read in data matches the parameter set of the loop, store it
		IF((ABS(obsvecSave(1)-varPass(1)).le.10.0_rKind**(-precis+1)).AND.(ABS(obsvecSave(2)-varPass(2)).le.10.0_rKind**(-precis+1))) THEN
		!Rewind the file to the beginning
				REWIND(538)

			!Loop over reading in data
kloop:			DO k=1,U_Res*(FLOOR(mu_Max+1)+1-MAX(FLOOR(mu_Min)-1,0))
		READ(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
			!If the read in matches j and has N'=N+1, keep it
			IF((ABS(obsvecSave(1)-obsvec(1)).le.10.0_rKind**(-precis+1)).AND.(ABS(obsvecSave(2)+1.0_rKind-obsvec(2)).le.10.0_rKind**(-precis+1))) THEN

			CALL openUnit(tempName,540,'A')			
!Compute the chemical potential from \mu(N+1)=E(N+1)-E(N), write data to file
			WRITE(540,'(7E30.15)') obsvec(1), obsvec(3)-obsvecSave(3), obsvec(3),obsvec(4),obsvec(5),obsvec(6),obsvec(7)
			CLOSE(540)
!Stop reading in files for this parameter set
EXIT kloop
			END IF

			
			END DO kloop
!Stop reading in  files for this parameter set
EXIT i1loop
		END IF


		END DO i1loop
!Rewind file for the next parameter set
			REWIND(538)
			CLOSE(538)

	END DO iloop

END DO jloop

!Get rid of the temporary file
CALL openUnit(outputName,538,'O')
CLOSE(UNIT=538, STATUS='DELETE')

!Postprocess the number non-conserving file to get parameters in the right order
ELSE

!Rewind the output file
CALL openUnit(outputName,538,'O')		
REWIND(538)

!Open a scratch file
OPEN(UNIT=542,STATUS='Scratch')
DO j=1,U_Res,1
varPass(1)=U_Min+(j-1)*(U_max-U_min)/(U_res-1)
	!Loop over mu
	DO i=1,mu_res,1
	varPass(2)=mu_min+(i-1)*(mu_max-mu_min)/(mu_res-1)

	ioSw1=0
		DO WHILE(ioSw1==0)
!Read in observables from output
		READ(538, '(7E30.15)',IOSTAT=ioSw1) obsvec(1),obsvec(2),obsvec(3),&
		obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		IF(ioSw1.ne.0) EXIT
!If the obseravbles match the parameters of the loop, write them to scratch
			IF((ABS(obsvec(1)-varPass(1)).le.10.0_rKind**(-precis+2)).AND.(ABS(obsvec(2)-varPass(2)).le.10.0_rKind**(-precis+2))) THEN
			WRITE(542,'(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
			obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
			EXIT
			END IF
		END DO
!Rewind observables file
		REWIND(538)
	END DO
END DO

REWIND(542)
CLOSE(538,STATUS='DELETE')
CALL openUnit(outputName,538)

!Write scratch to observables file
ioSw1=0
DO WHILE(ioSw1==0)
	READ(542, '(7E30.15)',IOSTAT=ioSw1) obsvec(1),obsvec(2),obsvec(3),&
	obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
		IF(ioSw1.ne.0) EXIT
	WRITE(538, '(7E30.15)') obsvec(1),obsvec(2),obsvec(3),&
	obsvec(4),obsvec(5),obsvec(6),obsvec(7) 
END DO

CLOSE(542)
CLOSE(538)


END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Postprocessing ends!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!end timer
end_time=MPI_WTIME()

!Output success message and timing data
	CALL openUnit(paramsName,537,'A')
	WRITE(537,*) 'Program has exited normally'
	WRITE(537,*) 'Total time taken=',end_time-start_time,' seconds'
	CLOSE(537)

!Get rid of the job pool and progress file
		CALL openUnit(progName,539)
		CLOSE(UNIT=539, STATUS='DELETE')

		CALL openUnit(jpName,541)
		CLOSE(UNIT=541, STATUS='DELETE')

END SUBROUTINE PhaseDiagramMasterFromU


SUBROUTINE PhaseDiagramBoseHubbardWorker(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
!
!Purpose: Worker's Routine for computing the phase diagram of the BH model in parallel
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
TYPE(matrix) :: oneBodyRho
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, j_Min, j_Max
INTEGER, INTENT(IN) :: mu_Res, j_Res
REAL(KIND=rKind) :: varPass(2)
REAL(KIND=rKind) :: obsvec(7), entList(systemSize)
CHARACTER(len=132) :: itpGName,itpLName,itpLLName,itpLRName
INTEGER :: myUnitG, myUnitL,myUnitLL, myUnitLR
INTEGER :: i,j ,k


myUnitG=790
myUnitL=791
myUnitLL=792
myUnitLR=793

!NC switch decision 
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
END IF

!Allocate Hamiltonian
CALL AllocateOps(H,systemSize-1,localSize*localSize)
!Allocate one-body rho
ALLOCATE(oneBodyRho%m(systemSize,systemSize))
!Create BH operators
CALL CreateFieldOps()

!Number non-conserving version
IF(.NOT.ncSwitch) THEN

DO i=1,j_Res*mu_Res
!Get a (j,mu) pair from master
		CALL MPI_Recv(varPass,2,my_mpi_rKind,0,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
jTunn=varPass(1)
mu0=varPass(2)

!Exit if the stop command is given
IF(varPass(1)==-100.0_rKind) EXIT

		!Define the Hamiltonian
		CALL HamiltonianBoseHubbard(H , jTunn, U0,mu0, V0)

		!Allocate Gammas and Lambdas
		CALL AllocateGamLam(Gammas, Lambdas, chiMin)

!Read in a previously generated ground state
IF(ITPreadMPDswitch==1) THEN
!Create the file name corresponding to the given parameters
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'PD_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(itpGName,'mu',4,mu0)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpGName).AND.CheckName(itpLName)) THEN
	CALL openUnit(itpGName,myUnitG,'O')
	CALL openUnit(itpLName,myUnitL,'O')
	!Read in the MPS
	CALL readGammaLambda(myUnitL, myUnitG,Gammas,Lambdas, ITPopenKind, chiMin)
CLOSE(myUnitG)
CLOSE(myUnitL)
	ELSE
	PRINT *, 'Pre-existing MPS not found!  Computing ground state!'
	CALL AllStates(Gammas,Lambdas)
	END IF
ELSE
		!Define the initial state
		CALL AllStates(Gammas,Lambdas)
END IF
		!Propagate in imaginary time
		CALL ImagTimeProp(H, Gammas, Lambdas, chiMin)

	obsvec=0.0_rKind
	obsvec(1)=jTunn
	obsvec(2)=mu0
		!Compute the total energy
		CALL TotalEnergy(obsvec(3), H, Gammas, Lambdas)
		!Compute the average filling
		CALL TotalNumber(obsvec(4), Gammas, Lambdas)
		!Compute the entropy at each site
		CALL LocalEntropyDist(entList, Gammas, Lambdas)
		!Find the site-averaged entropy
			DO k=1,systemSize
			obsvec(5)=obsvec(5)+entList(k)/(systemSize*1.0_rKind)
			END DO
		!Compute the Q-measure
		obsvec(6)=MeyerQmeasure(Gammas, Lambdas)
!!! Calculate the single-particle density matrix.
CALL TwoSiteExpVal(oneBodyRho%m,MATMUL(TRANSPOSE(a_op%mr),a_op%mr),TensorProd(TRANSPOSE(a_op%mr),a_op%mr),Gammas,Lambdas)
		!Calculate the depletion from the spdm
		CALL Qdepletion(obsvec(7), oneBodyRho%m, REAL(obsvec(4)*systemSize,KIND=rKind))

!Output the ground state for later use
IF(ITPwriteMPDswitch==1) THEN
!Redefine name with chiMax
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'PD_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMax) !Note this is chiMax for writing!
CALL appendBaseName(itpGName,'mu',4,mu0)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
!Open the files (prepared to write over old data)
	CALL openUnit(itpGName,myUnitG)
	CALL openUnit(itpLName,myUnitL)
!Write out the local tensors
CALL RecordLambdas(myUnitL, Lambdas,ITPopenKind)
CALL RecordGammas(myUnitG, Gammas, ITPopenKind)
CLOSE(myUnitG)
CLOSE(myUnitL)
END IF

	!Clean up
	CALL DeallocateGamLam(Gammas, Lambdas)

	!Send rank to the master to let it know you are finished
	CALL MPI_Send(my_rank,1,MPI_INTEGER,0,&
		2000+my_rank, MPI_COMM_WORLD,ierror)

	!When master responds, send observables
	CALL MPI_Send(obsvec,7,my_mpi_rKind,0,&
		my_rank, MPI_COMM_WORLD,ierror)

END DO


!Number conserving version
ELSE

mu0=0.0_rKind

DO i=1,j_Res*(FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1)
!Get a (j,N) pair
		CALL MPI_Recv(varPass,2,my_mpi_rKind,0,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
jTunn=varPass(1)
totNum=FLOOR(varPass(2))

!Exit if the stop command is given
IF(varPass(1)==-100.0_rKind) EXIT

		!Define the Hamiltonian
		CALL HamiltonianBoseHubbard(H , jTunn, U0,mu0, V0)

		!Allocate Gammas, Lambdas, and Labels
		CALL AllocateGamLam(Gammas, Lambdas, chiMin)
		CALL AllocateLabel(LabelLeft, LabelRight, chiMin)

!Read in a previously generated ground state
IF(ITPreadMPDswitch==1) THEN
!Create the file name corresponding to the given parameters
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'PDNC_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(itpGName,'N',totNum)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL copyName(itpGName,itpLLName)
CALL copyName(itpGName,itpLRName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
CALL appendBaseName(itpLLName,'itpLL.dat')
CALL appendBaseName(itpLRName,'itpLR.dat')

!Check to see if a state exists
	IF(CheckName(itpGName).AND.CheckName(itpLName).AND.CheckName(itpLLName).AND.CheckName(itpLRName)) THEN
	CALL openUnit(itpGName,myUnitG,'O')
	CALL openUnit(itpLName,myUnitL,'O')
	CALL openUnit(itpLLName,myUnitLL,'O')
	CALL openUnit(itpLRName,myUnitLR,'O')
	!Read in the MPS
CALL readGammaLambdaLabels(myUnitL, myUnitG,myUnitLL, myUnitLR,Gammas,Lambdas,LabelLeft, LabelRight, ITPopenKind, chiMin)
CLOSE(myUnitG)
CLOSE(myUnitL)
CLOSE(myUnitLL)
CLOSE(myUnitLR)
	ELSE
PRINT *, 'Pre-existing MPS not found!  Computing ground state!'
		CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
	END IF
ELSE
	!Define the initial state consistent with number conservation
	CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
END IF

		!Propagate in imaginary time
		CALL ImagTimePropNC(H, Gammas, Lambdas, LabelLeft, LabelRight, chiMin)

	obsvec=0.0_rKind
	obsvec(1)=jTunn
	obsvec(2)=totNum*1.0_rKind
		!Compute the total energy
		CALL TotalEnergy(obsvec(3), H, Gammas, Lambdas)
		!Compute the average filling
		obsvec(4)=totNum*1.0_rKind/(systemSize*1.0_rKind)
		!Compute the entropy at each site
		CALL LocalEntropyDist(entList, Gammas, Lambdas)
		!Find the site-averaged entropy
			DO k=1,systemSize
			obsvec(5)=obsvec(5)+entList(k)/(systemSize*1.0_rKind)
			END DO
		!Compute the Q-measure
		obsvec(6)=MeyerQmeasure(Gammas, Lambdas)
!!! Calculate the single-particle density matrix.
CALL TwoSiteExpVal(oneBodyRho%m,MATMUL(TRANSPOSE(a_op%mr),a_op%mr),TensorProd(TRANSPOSE(a_op%mr),a_op%mr),Gammas,Lambdas)
		!Calculate the depletion from the spdm
		CALL Qdepletion(obsvec(7), oneBodyRho%m, REAL(obsvec(4)*systemSize,KIND=rKind))

!Output the ground state for later use
IF(ITPwriteMPDswitch==1) THEN
!Redefine name with chiMax
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'PDNC_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMax) !Note this is chiMax for writing!
CALL appendBaseName(itpGName,'N',totNum)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL copyName(itpGName,itpLLName)
CALL copyName(itpGName,itpLRName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
CALL appendBaseName(itpLLName,'itpLL.dat')
CALL appendBaseName(itpLRName,'itpLR.dat')

!Open the files (prepared to write over old data)
	CALL openUnit(itpGName,myUnitG)
	CALL openUnit(itpLName,myUnitL)
	CALL openUnit(itpLLName,myUnitLL)
	CALL openUnit(itpLRName,myUnitLR)
!Write out the local tensors
CALL RecordLambdas(myUnitL, Lambdas,ITPopenKind)
CALL RecordGammas(myUnitG, Gammas, ITPopenKind)
CALL RecordLabel(myUnitLL, LabelLeft)
CALL RecordLabel(myUnitLR, LabelRight)
CLOSE(myUnitG)
CLOSE(myUnitL)
CLOSE(myUnitLL)
CLOSE(myUnitLR)
END IF

		!Clean up
		CALL DeallocateLabel(LabelLeft, LabelRight)
		CALL DeallocateGamLam(Gammas, Lambdas)

	!Send rank to the master to let it know you are finished
	CALL MPI_Send(my_rank,1,MPI_INTEGER,0,&
		2000+my_rank, MPI_COMM_WORLD,ierror)

	!When master responds, send observables
	CALL MPI_Send(obsvec,7,my_mpi_rKind,0,&
		my_rank, MPI_COMM_WORLD,ierror)

END DO
END IF

!Clean up
CALL DestroyFieldOps()
CALL DeallocateOps(H,systemSize-1)
DEALLOCATE(oneBodyRho%m)

END SUBROUTINE PhaseDiagramBoseHubbardWorker

SUBROUTINE PhaseDiagramFermiHubbardWorker(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
!
!Purpose: Worker's Routine for computing the phase diagram of the FH model in parallel
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
TYPE(matrix) :: oneBodyRho
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, j_Min, j_Max
INTEGER, INTENT(IN) :: mu_Res, j_Res
REAL(KIND=rKind) :: varPass(2), numVec(systemSize) 
COMPLEX(KIND=rKind) :: COvec(systemSize+1), Mvec(systemSize+1), singSupOP
COMPLEX(KIND=rKind) :: COparam, COparamvec(systemSize,systemSize)
COMPLEX(KIND=rKind) :: pairingParam, pairingParamvec(systemSize,systemSize)
REAL(KIND=rKind) :: obsvec(7), entList(systemSize), dumNum
CHARACTER(len=132) :: itpGName,itpLName,itpLLName,itpLRName, OPName
INTEGER :: myUnitG, myUnitL,myUnitLL, myUnitLR
INTEGER :: i,j ,k, ind1, ind2, ind3
TYPE(matrixReal) :: num_op,num_opP
COMPLEX(KIND=rKind) :: eye

eye=CMPLX(0.0,1.0,KIND=rKind)

myUnitG=790
myUnitL=791
myUnitLL=792
myUnitLR=793

!NC switch decision 
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
END IF

!Allocate Hamiltonian
CALL AllocateOps(H,systemSize-1,localSize*localSize)
!Allocate one-body rho
ALLOCATE(oneBodyRho%m(systemSize,systemSize))
!Create FH operators
CALL CreateFermiSops()

ALLOCATE(num_op%mr(localSize, localSize))
ALLOCATE(num_opP%mr(localSize, localSize))

!Number non-conserving version
IF(.NOT.ncSwitch) THEN

!Create the file name corresponding to the given parameters
CALL createFileName(OPname,pdDir)
CALL appendBaseName(OPname,'FHPD_')
CALL appendBaseName(OPname,'L',systemSize)
CALL appendBaseName(OPname,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(OPname,'mu',4,mu0)
CALL appendBaseName(OPname,'j',4,jtunn)
CALL appendBaseName(OPname,'rank',my_rank)
CALL appendBaseName(OPname,'OP.dat')

!CALL openUnit(OPname,831)

DO i=1,j_Res*mu_Res
!Get a (j,mu) pair
		CALL MPI_Recv(varPass,2,my_mpi_rKind,0,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
jTunn=varPass(1)
mu0=varPass(2)

!Exit if the stop command is given
IF(varPass(1)==-100.0_rKind) EXIT

		!Redefine the Hamiltonian
		CALL HamiltonianHubbard(H,jTunn,U0,mu0)

		!Allocate Gammas and Lambdas
		CALL AllocateGamLam(Gammas, Lambdas, chiMin)

!Read in a previously generated ground state
IF(ITPreadMPDswitch==1) THEN
!Create the file name corresponding to the given parameters
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPD_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(itpGName,'mu',4,mu0)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpGName).AND.CheckName(itpLName)) THEN
	CALL openUnit(itpGName,myUnitG,'O')
	CALL openUnit(itpLName,myUnitL,'O')
	!Read in the MPS
	CALL readGammaLambda(myUnitL, myUnitG,Gammas,Lambdas, ITPopenKind, chiMin)
CLOSE(myUnitG)
CLOSE(myUnitL)
	ELSE
	PRINT *, 'Pre-existing MPS not found!  Computing ground state!'
	CALL AllStates(Gammas,Lambdas)
	END IF
ELSE
		!Define the initial state
		CALL AllStates(Gammas,Lambdas)
END IF
		!Propagate in imaginary time
		CALL ImagTimeProp(H, Gammas, Lambdas, chiMin)

	obsvec=0.0_rKind
	obsvec(1)=jTunn
	obsvec(2)=mu0
		!Compute the average filling
		CALL TotalNumber(dumNum, Gammas, Lambdas,1)
		obsVec(3)=dumNum
		obsVec(4)=dumNum
		CALL TotalNumber(dumNum, Gammas, Lambdas,2)
		obsVec(3)=obsVec(3)+dumNum
!Magnetization
		obsvec(4)=obsVec(4)-dumNum
!Q-measure
		obsvec(5)=MeyerQmeasure(Gammas, Lambdas)
!Charge-density correlation evaluated at q=Pi
obsvec(6)=0.0_rKind
CALL TwoSiteExpValG(COparamvec,&
MATMUL(TRANSPOSE(a_opS(1)%mr),a_opS(1)%mr)+MATMUL(TRANSPOSE(a_opS(2)%mr),a_opS(2)%mr)&
,MATMUL(TRANSPOSE(a_opS(1)%mr),a_opS(1)%mr)+MATMUL(TRANSPOSE(a_opS(2)%mr),a_opS(2)%mr)&
,Gammas,Lambdas,1)
DO ind1=1,systemSize
	DO ind2=1,systemSize
	obsvec(6)=obsvec(6)+((-1.0_rKind)**(ind1-ind2))*COparamvec(ind1,ind2)
	END DO
END DO

obsvec(7)=0.0_rKind
!Pairing order parameter evaluated at q=Pi
CALL TwoSiteExpValG(COparamvec,&
MATMUL(TRANSPOSE(a_opS(1)%mr),TRANSPOSE(a_opS(2)%mr))&
,MATMUL(a_opS(2)%mr,a_opS(1)%mr)&
,Gammas,Lambdas,1)
DO ind1=1,systemSize
	DO ind2=1,systemSize
	obsvec(7)=obsvec(7)+COparamvec(ind1,ind2)
	END DO
END DO

CALL TwoSiteExpValG(COparamvec,&
MATMUL(a_opS(2)%mr,a_opS(1)%mr)&
,MATMUL(TRANSPOSE(a_opS(1)%mr),TRANSPOSE(a_opS(2)%mr))&
,Gammas,Lambdas,1)
DO ind1=1,systemSize
	DO ind2=1,systemSize
	obsvec(7)=obsvec(7)+COparamvec(ind1,ind2)
	END DO
END DO




!Output the ground state for later use
IF(ITPwriteMPDswitch==1) THEN
!Redefine name with chiMax
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPD_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMax) !Note this is chiMax for writing!
CALL appendBaseName(itpGName,'mu',4,mu0)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
!Open the files (prepared to write over old data)
	CALL openUnit(itpGName,myUnitG)
	CALL openUnit(itpLName,myUnitL)
!Write out the local tensors
CALL RecordLambdas(myUnitL, Lambdas,ITPopenKind)
CALL RecordGammas(myUnitG, Gammas, ITPopenKind)
CLOSE(myUnitG)
CLOSE(myUnitL)
END IF

	!Clean up
	CALL DeallocateGamLam(Gammas, Lambdas)

	!Send rank to the master to let it know you are finished
	CALL MPI_Send(my_rank,1,MPI_INTEGER,0,&
		2000+my_rank, MPI_COMM_WORLD,ierror)

	!When master responds, send observables
	CALL MPI_Send(obsvec,7,my_mpi_rKind,0,&
		my_rank, MPI_COMM_WORLD,ierror)

END DO

!CLOSE(831)

!Number conserving version
ELSE

mu0=0.0_rKind

DO i=1,j_Res*(FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1)
!Get a (j,N) pair
		CALL MPI_Recv(varPass,2,my_mpi_rKind,0,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
jTunn=varPass(1)
totNum=FLOOR(varPass(2))

!Exit if the stop command is given
IF(varPass(1)==-100.0_rKind) EXIT

		!Redefine the Hamiltonian
		CALL HamiltonianHubbard(H ,jTunn,U0, mu0)

		!Allocate Gammas, Lambdas, and Labels
		CALL AllocateGamLam(Gammas, Lambdas, chiMin)
		CALL AllocateLabel(LabelLeft, LabelRight, chiMin)

!Read in a previously generated ground state
IF(ITPreadMPDswitch==1) THEN
!Create the file name corresponding to the given parameters
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPDNC_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(itpGName,'N',totNum)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL copyName(itpGName,itpLLName)
CALL copyName(itpGName,itpLRName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
CALL appendBaseName(itpLLName,'itpLL.dat')
CALL appendBaseName(itpLRName,'itpLR.dat')

!Check to see if a state exists
	IF(CheckName(itpGName).AND.CheckName(itpLName).AND.CheckName(itpLLName).AND.CheckName(itpLRName)) THEN
	CALL openUnit(itpGName,myUnitG,'O')
	CALL openUnit(itpLName,myUnitL,'O')
	CALL openUnit(itpLLName,myUnitLL,'O')
	CALL openUnit(itpLRName,myUnitLR,'O')
	!Read in the MPS
CALL readGammaLambdaLabels(myUnitL, myUnitG,myUnitLL, myUnitLR,Gammas,Lambdas,LabelLeft, LabelRight, ITPopenKind, chiMin)
CLOSE(myUnitG)
CLOSE(myUnitL)
CLOSE(myUnitLL)
CLOSE(myUnitLR)
	ELSE
PRINT *, 'Pre-existing MPS not found!  Computing ground state!'
		CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
	END IF
ELSE
	!Define the initial state consistent with number conservation
	CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
END IF

		!Propagate in imaginary time
		CALL ImagTimePropNC(H, Gammas, Lambdas, LabelLeft, LabelRight, chiMin)

	obsvec=0.0_rKind
	obsvec(1)=jTunn
	obsvec(2)=totNum*1.0_rKind
		!Compute the total energy
		CALL TotalEnergy(obsvec(3), H, Gammas, Lambdas)
		!Compute the average filling
		obsvec(4)=totNum*1.0_rKind/(systemSize*1.0_rKind)
		!Compute the entropy at each site
		CALL LocalEntropyDist(entList, Gammas, Lambdas)
		!Find the site-average entropy
			DO k=1,systemSize
			obsvec(5)=obsvec(5)+entList(k)/(systemSize*1.0_rKind)
			END DO
		obsvec(6)=MeyerQmeasure(Gammas, Lambdas)
!!! Calculate the single-particle density matrix with fermi phase.
CALL TwoSiteExpVal(oneBodyRho%m,MATMUL(TRANSPOSE(a_op%mr),a_op%mr),TensorProd(TRANSPOSE(a_op%mr),a_op%mr),Gammas,Lambdas,1)
!Calculate the depletion from the spdm
		CALL Qdepletion(obsvec(7), oneBodyRho%m, REAL(obsvec(4)*systemSize,KIND=rKind))

!Output the ground state for later use
IF(ITPwriteMPDswitch==1) THEN
!Redefine name with chiMax
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPDNC_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMax) !Note this is chiMax for writing!
CALL appendBaseName(itpGName,'N',totNum)
CALL appendBaseName(itpGName,'j',4,jtunn)
CALL copyName(itpGName,itpLName)
CALL copyName(itpGName,itpLLName)
CALL copyName(itpGName,itpLRName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
CALL appendBaseName(itpLLName,'itpLL.dat')
CALL appendBaseName(itpLRName,'itpLR.dat')

!Open the files (prepared to write over old data)
	CALL openUnit(itpGName,myUnitG)
	CALL openUnit(itpLName,myUnitL)
	CALL openUnit(itpLLName,myUnitLL)
	CALL openUnit(itpLRName,myUnitLR)
!Write out the local tensors
CALL RecordLambdas(myUnitL, Lambdas,ITPopenKind)
CALL RecordGammas(myUnitG, Gammas, ITPopenKind)
CALL RecordLabel(myUnitLL, LabelLeft)
CALL RecordLabel(myUnitLR, LabelRight)
CLOSE(myUnitG)
CLOSE(myUnitL)
CLOSE(myUnitLL)
CLOSE(myUnitLR)
END IF

		!Clean up
		CALL DeallocateLabel(LabelLeft, LabelRight)
		CALL DeallocateGamLam(Gammas, Lambdas)

	!Send rank to the master to let it know you are finished
	CALL MPI_Send(my_rank,1,MPI_INTEGER,0,&
		2000+my_rank, MPI_COMM_WORLD,ierror)

	!When master responds, send observables
	CALL MPI_Send(obsvec,7,my_mpi_rKind,0,&
		my_rank, MPI_COMM_WORLD,ierror)

END DO
END IF

!Clean up
CALL DestroyFermiSops()
CALL DeallocateOps(H,systemSize-1)
DEALLOCATE(oneBodyRho%m)

END SUBROUTINE PhaseDiagramFermiHubbardWorker


SUBROUTINE PhaseDiagramFermiHubbardWorkerFromU(mu_Min, mu_Max, mu_Res, U_Min, U_Max, U_Res)
!
!Purpose: Worker's Routine for computing the phase diagram of the FH model in parallel
!
IMPLICIT NONE
TYPE(tensor), POINTER :: Gammas(:)
TYPE(vector), POINTER :: Lambdas(:)
TYPE(vectorInt), POINTER :: LabelLeft(:), LabelRight(:)
TYPE(matrix) :: oneBodyRho
TYPE(matrix), POINTER :: H(:)
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, U_Min, U_Max
INTEGER, INTENT(IN) :: mu_Res, U_Res
REAL(KIND=rKind) :: varPass(2), numVec(systemSize) 
COMPLEX(KIND=rKind) :: COvec(systemSize+1), Mvec(systemSize+1), singSupOP
COMPLEX(KIND=rKind) :: COparam, COparamvec(systemSize,systemSize)
COMPLEX(KIND=rKind) :: pairingParam, pairingParamvec(systemSize,systemSize)
REAL(KIND=rKind) :: obsvec(7), entList(systemSize), dumNum
CHARACTER(len=132) :: itpGName,itpLName,itpLLName,itpLRName, OPName
INTEGER :: myUnitG, myUnitL,myUnitLL, myUnitLR
INTEGER :: i,j ,k, ind1, ind2, ind3
TYPE(matrixReal) :: num_op,num_opP
COMPLEX(KIND=rKind) :: eye

eye=CMPLX(0.0,1.0,KIND=rKind)

myUnitG=790
myUnitL=791
myUnitLL=792
myUnitLR=793

!NC switch decision 
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
END IF

!Allocate Hamiltonian
CALL AllocateOps(H,systemSize-1,localSize*localSize)
!Allocate one-body rho
ALLOCATE(oneBodyRho%m(systemSize,systemSize))
!Create FH operators
CALL CreateFermiSops()

ALLOCATE(num_op%mr(localSize, localSize))
ALLOCATE(num_opP%mr(localSize, localSize))

!Number non-conserving version
IF(.NOT.ncSwitch) THEN

!Create the file name corresponding to the given parameters
CALL createFileName(OPname,pdDir)
CALL appendBaseName(OPname,'FHPD_')
CALL appendBaseName(OPname,'L',systemSize)
CALL appendBaseName(OPname,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(OPname,'mu',4,mu0)
CALL appendBaseName(OPname,'U',4,U0)
CALL appendBaseName(OPname,'rank',my_rank)
CALL appendBaseName(OPname,'OP.dat')

!CALL openUnit(OPname,831)

DO i=1,U_Res*mu_Res
!Get a (j,mu) pair
		CALL MPI_Recv(varPass,2,my_mpi_rKind,0,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
U0=varPass(1)
mu0=varPass(2)

!Exit if the stop command is given
IF(varPass(1)==-100.0_rKind) EXIT

		!Redefine the Hamiltonian
		CALL HamiltonianHubbard(H,jTunn,U0,mu0)

		!Allocate Gammas and Lambdas
		CALL AllocateGamLam(Gammas, Lambdas, chiMin)

!Read in a previously generated ground state
IF(ITPreadMPDswitch==1) THEN
!Create the file name corresponding to the given parameters
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPD_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(itpGName,'mu',4,mu0)
CALL appendBaseName(itpGName,'U',4,U0)
CALL copyName(itpGName,itpLName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpGName).AND.CheckName(itpLName)) THEN
	CALL openUnit(itpGName,myUnitG,'O')
	CALL openUnit(itpLName,myUnitL,'O')
	!Read in the MPS
	CALL readGammaLambda(myUnitL, myUnitG,Gammas,Lambdas, ITPopenKind, chiMin)
CLOSE(myUnitG)
CLOSE(myUnitL)
	ELSE
	PRINT *, 'Pre-existing MPS not found!  Computing ground state!'
	CALL AllStates(Gammas,Lambdas)
	END IF
ELSE
		!Define the initial state
		CALL AllStates(Gammas,Lambdas)
END IF
		!Propagate in imaginary time
		CALL ImagTimeProp(H, Gammas, Lambdas, chiMin)

	obsvec=0.0_rKind
	obsvec(1)=U0
	obsvec(2)=mu0
		!Compute the average filling
		CALL TotalNumber(dumNum, Gammas, Lambdas,1)
		obsVec(3)=dumNum
		obsVec(4)=dumNum
		CALL TotalNumber(dumNum, Gammas, Lambdas,2)
		obsVec(3)=obsVec(3)+dumNum
!Magnetization
		obsvec(4)=obsVec(4)-dumNum
!Q-measure
		obsvec(5)=MeyerQmeasure(Gammas, Lambdas)
!Charge-density correlation evaluated at q=Pi
obsvec(6)=0.0_rKind
CALL TwoSiteExpValG(COparamvec,&
MATMUL(TRANSPOSE(a_opS(1)%mr),a_opS(1)%mr)+MATMUL(TRANSPOSE(a_opS(2)%mr),a_opS(2)%mr)&
,MATMUL(TRANSPOSE(a_opS(1)%mr),a_opS(1)%mr)+MATMUL(TRANSPOSE(a_opS(2)%mr),a_opS(2)%mr)&
,Gammas,Lambdas,1)
DO ind1=1,systemSize
	DO ind2=1,systemSize
	obsvec(6)=obsvec(6)+((-1.0_rKind)**(ind1-ind2))*COparamvec(ind1,ind2)
	END DO
END DO

obsvec(7)=0.0_rKind
!Pairing order parameter evaluated at q=Pi
CALL TwoSiteExpValG(COparamvec,&
MATMUL(TRANSPOSE(a_opS(1)%mr),TRANSPOSE(a_opS(2)%mr))&
,MATMUL(a_opS(2)%mr,a_opS(1)%mr)&
,Gammas,Lambdas,1)
DO ind1=1,systemSize
	DO ind2=1,systemSize
	obsvec(7)=obsvec(7)+COparamvec(ind1,ind2)
	END DO
END DO

CALL TwoSiteExpValG(COparamvec,&
MATMUL(a_opS(2)%mr,a_opS(1)%mr)&
,MATMUL(TRANSPOSE(a_opS(1)%mr),TRANSPOSE(a_opS(2)%mr))&
,Gammas,Lambdas,1)
DO ind1=1,systemSize
	DO ind2=1,systemSize
	obsvec(7)=obsvec(7)+COparamvec(ind1,ind2)
	END DO
END DO


!Output the ground state for later use
IF(ITPwriteMPDswitch==1) THEN
!Redefine name with chiMax
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPD_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMax) !Note this is chiMax for writing!
CALL appendBaseName(itpGName,'mu',4,mu0)
CALL appendBaseName(itpGName,'U',4,U0)
CALL copyName(itpGName,itpLName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
!Open the files (prepared to write over old data)
	CALL openUnit(itpGName,myUnitG)
	CALL openUnit(itpLName,myUnitL)
!Write out the local tensors
CALL RecordLambdas(myUnitL, Lambdas,ITPopenKind)
CALL RecordGammas(myUnitG, Gammas, ITPopenKind)
CLOSE(myUnitG)
CLOSE(myUnitL)
END IF

	!Clean up
	CALL DeallocateGamLam(Gammas, Lambdas)

	!Send rank to the master to let it know you are finished
	CALL MPI_Send(my_rank,1,MPI_INTEGER,0,&
		2000+my_rank, MPI_COMM_WORLD,ierror)

	!When master responds, send observables
	CALL MPI_Send(obsvec,7,my_mpi_rKind,0,&
		my_rank, MPI_COMM_WORLD,ierror)

END DO

!CLOSE(831)

!Number conserving version
ELSE

mu0=0.0_rKind

DO i=1,U_Res*(FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1)
!Get a (j,N) pair
		CALL MPI_Recv(varPass,2,my_mpi_rKind,0,&
		 MPI_ANY_TAG, MPI_COMM_WORLD, my_status,ierror)
U0=varPass(1)
totNum=FLOOR(varPass(2))

!Exit if the stop command is given
IF(varPass(1)==-100.0_rKind) EXIT

		!Redefine the Hamiltonian
		CALL HamiltonianHubbard(H ,jTunn,U0, mu0)

		!Allocate Gammas, Lambdas, and Labels
		CALL AllocateGamLam(Gammas, Lambdas, chiMin)
		CALL AllocateLabel(LabelLeft, LabelRight, chiMin)

!Read in a previously generated ground state
IF(ITPreadMPDswitch==1) THEN
!Create the file name corresponding to the given parameters
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPDNC_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMin) !Note this is chiMin for reading!
CALL appendBaseName(itpGName,'N',totNum)
CALL appendBaseName(itpGName,'U',4,U0)
CALL copyName(itpGName,itpLName)
CALL copyName(itpGName,itpLLName)
CALL copyName(itpGName,itpLRName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
CALL appendBaseName(itpLLName,'itpLL.dat')
CALL appendBaseName(itpLRName,'itpLR.dat')

!Check to see if a state exists
	IF(CheckName(itpGName).AND.CheckName(itpLName).AND.CheckName(itpLLName).AND.CheckName(itpLRName)) THEN
	CALL openUnit(itpGName,myUnitG,'O')
	CALL openUnit(itpLName,myUnitL,'O')
	CALL openUnit(itpLLName,myUnitLL,'O')
	CALL openUnit(itpLRName,myUnitLR,'O')
	!Read in the MPS
CALL readGammaLambdaLabels(myUnitL, myUnitG,myUnitLL, myUnitLR,Gammas,Lambdas,LabelLeft, LabelRight, ITPopenKind, chiMin)
CLOSE(myUnitG)
CLOSE(myUnitL)
CLOSE(myUnitLL)
CLOSE(myUnitLR)
	ELSE
PRINT *, 'Pre-existing MPS not found!  Computing ground state!'
		CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
	END IF
ELSE
	!Define the initial state consistent with number conservation
	CALL InitialSetNC(Gammas, Lambdas, LabelLeft, LabelRight)
END IF

		!Propagate in imaginary time
		CALL ImagTimePropNC(H, Gammas, Lambdas, LabelLeft, LabelRight, chiMin)

	obsvec=0.0_rKind
	obsvec(1)=U0
	obsvec(2)=totNum*1.0_rKind
		!Compute the total energy
		CALL TotalEnergy(obsvec(3), H, Gammas, Lambdas)
		!Compute the average filling
		obsvec(4)=totNum*1.0_rKind/(systemSize*1.0_rKind)
		!Compute the entropy at each site
		CALL LocalEntropyDist(entList, Gammas, Lambdas)
		!Find the site-average entropy
			DO k=1,systemSize
			obsvec(5)=obsvec(5)+entList(k)/(systemSize*1.0_rKind)
			END DO
		obsvec(6)=MeyerQmeasure(Gammas, Lambdas)
!!! Calculate the single-particle density matrix with fermi phase.
CALL TwoSiteExpVal(oneBodyRho%m,MATMUL(TRANSPOSE(a_op%mr),a_op%mr),TensorProd(TRANSPOSE(a_op%mr),a_op%mr),Gammas,Lambdas,1)
!Calculate the depletion from the spdm
		CALL Qdepletion(obsvec(7), oneBodyRho%m, REAL(obsvec(4)*systemSize,KIND=rKind))

!Output the ground state for later use
IF(ITPwriteMPDswitch==1) THEN
!Redefine name with chiMax
CALL createFileName(itpGName,itpDir)
CALL appendBaseName(itpGName,'FHPDNC_')
CALL appendBaseName(itpGName,'L',systemSize)
CALL appendBaseName(itpGName,'Chi',chiMax) !Note this is chiMax for writing!
CALL appendBaseName(itpGName,'N',totNum)
CALL appendBaseName(itpGName,'U',4,U0)
CALL copyName(itpGName,itpLName)
CALL copyName(itpGName,itpLLName)
CALL copyName(itpGName,itpLRName)
CALL appendBaseName(itpGName,'itpG.dat')
CALL appendBaseName(itpLName,'itpL.dat')
CALL appendBaseName(itpLLName,'itpLL.dat')
CALL appendBaseName(itpLRName,'itpLR.dat')

!Open the files (prepared to write over old data)
	CALL openUnit(itpGName,myUnitG)
	CALL openUnit(itpLName,myUnitL)
	CALL openUnit(itpLLName,myUnitLL)
	CALL openUnit(itpLRName,myUnitLR)
!Write out the local tensors
CALL RecordLambdas(myUnitL, Lambdas,ITPopenKind)
CALL RecordGammas(myUnitG, Gammas, ITPopenKind)
CALL RecordLabel(myUnitLL, LabelLeft)
CALL RecordLabel(myUnitLR, LabelRight)
CLOSE(myUnitG)
CLOSE(myUnitL)
CLOSE(myUnitLL)
CLOSE(myUnitLR)
END IF

		!Clean up
		CALL DeallocateLabel(LabelLeft, LabelRight)
		CALL DeallocateGamLam(Gammas, Lambdas)

	!Send rank to the master to let it know you are finished
	CALL MPI_Send(my_rank,1,MPI_INTEGER,0,&
		2000+my_rank, MPI_COMM_WORLD,ierror)

	!When master responds, send observables
	CALL MPI_Send(obsvec,7,my_mpi_rKind,0,&
		my_rank, MPI_COMM_WORLD,ierror)

END DO
END IF

!Clean up
CALL DestroyFermiSops()
CALL DeallocateOps(H,systemSize-1)
DEALLOCATE(oneBodyRho%m)

END SUBROUTINE PhaseDiagramFermiHubbardWorkerFromU

SUBROUTINE FidelityAnalysis(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
!
!Purpose: Read in Gammas and Lambdas from output files to compute
!Fidelities for Bose code
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, j_Min, j_Max
INTEGER, INTENT(IN) :: mu_Res, j_Res
TYPE(tensor), POINTER :: Gammas1(:),Gammas2(:)
TYPE(vector), POINTER :: Lambdas1(:), Lambdas2(:)
INTEGER i,j,k, counter, ioSw1, ioSw2, counterVec(0:num_cores-1), myUnit
REAL(KIND=rKind) :: jSave(2), muSave(2), varSave(2),dummies(2)
INTEGER :: numSave(2), numactive, my_b(2)
REAL(KIND=rKind) :: mu_low, mu_high, mu_diff
REAL(KIND=rKind) :: obseVec(j_Res), jDummy
COMPLEX(KIND=rKind) :: fidel
CHARACTER(len=132) :: outputmuName,outputJName, FinaloutputName,FinaloutputNameJ
CHARACTER(len=132) :: itpG1name, itpG2name, itpL1name, itpL2name
CHARACTER(len=132) :: outputName(0:num_cores-1)

!We cannot use more cores than we have lines
IF(num_cores.gt.mu_Res) THEN
numactive=mu_Res
ELSE
numactive=num_cores
END IF

IF(my_rank.le.numactive-1) THEN

CALL createFileName(outputmuName,pdDir)

!NC switch decision 
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
	IF((FLOOR(mu_Min).ne.CEILING(mu_Min)).OR.(FLOOR(mu_Max).ne.CEILING(mu_Max))) THEN
		PRINT *, 'Warning: The number conserving method is being used but either mu_Min or mu_Max is not integer.'
		PRINT *,  'The floor is being taken!'
		PRINT *, ''
	END IF
	CALL appendBaseName(outputmuName,'FANC_')
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
	CALL appendBaseName(outputmuName,'FA_')
END IF


CALL appendBaseName(outputmuName,'L',systemSize)
CALL appendBaseName(outputmuName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(outputmuName,'muR',4,mu_Min)
	ELSE
	CALL appendBaseName(outputmuName,'N',4,mu_Min)
END IF
CALL appendBaseName(outputmuName,'t',4,mu_Max)
CALL appendBaseName(outputmuName,'jR',4,j_Min)
CALL appendBaseName(outputmuName,'t',4,j_Max)
CALL appendBaseName(outputmuName,'rank',my_rank)

!Copy i/o to all files
CALL copyName(outputmuName,outputJName)
CALL appendBaseName(outputmuName,'FixedMuOutput.dat')
CALL appendBaseName(outputJName,'FixedJOutput.dat')

CALL openUnit(outputmuName,538)
CALL openUnit(outputJName,539)
CLOSE(538)
CLOSE(539)

!!!!!!!!!!!!!!FIXED MU FIDELITY ANALYSIS!!!!!!!!!!!!!!!!!1
CALL AllocateGamLam(Gammas1, Lambdas1, chiMax)
CALL AllocateGamLam(Gammas2, Lambdas2, chiMax)

IF(.NOT.ncSwitch) THEN

	my_b(1)=1+my_rank*FLOOR(mu_Res*1.0_rKind/((numactive)*1.0_rKind))
	my_b(2)=(my_rank+1)*FLOOR(mu_Res*1.0_rKind/((numactive)*1.0_rKind))

	IF(my_rank==numactive-1) THEN
		my_b(2)=my_b(2)+MOD(mu_Res,numactive)
	END IF

	DO i=my_b(1),my_b(2),1
		mu0=mu_Min+(i-1)*(mu_max-mu_min)/((mu_Res-1)*1.0_rKind)
		DO j=1,j_Res-1

		jSave(1)=j_Min+(j-1)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)
		jSave(2)=j_Min+(j)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'PD_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'mu',4,mu0)
CALL appendBaseName(itpG1Name,'j',4,jSave(1))
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'PD_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'mu',4,mu0)
CALL appendBaseName(itpG2Name,'j',4,jSave(2))
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
		CALL openUnit(itpG1Name,540,'O')
		CALL openUnit(itpL1Name,541,'O')
		!Read in the MPS
		CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
		CLOSE(540)
		CLOSE(541)
	ELSE
		PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
		CALL openUnit(itpG2Name,542,'O')
		CALL openUnit(itpL2Name,543,'O')
		!Read in the MPS
		CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
		CLOSE(542)
		CLOSE(543)
	ELSE
		PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap
fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)

CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jSave(1),mu0,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((J_Max-j_Min)/((j_res-1)*1.0_rKind))**2))

CLOSE(538)

		END DO
	END DO

ELSE

mu_low=(Min(FLOOR(mu_Min)-1,0))*1.0_rKind
mu_high=(FLOOR(mu_Max)+1)*1.0_rKind
mu_diff=mu_low-mu_High+1

my_b(1)=1+my_rank*FLOOR((mu_diff*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((mu_diff*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(FLOOR(mu_diff),numactive)
END IF

DO i=my_b(1),my_b(2),1

	totNum=i	
	DO j=1,j_Res-1

		jSave(1)=j_Min+(j-1)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)
		jSave(2)=j_Min+(j)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'PDNC_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'N',totNum)
CALL appendBaseName(itpG1Name,'j',4,jSave(1))
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'PDNC_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'N',totNum)
CALL appendBaseName(itpG2Name,'j',4,jSave(2))
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
		CALL openUnit(itpG1Name,540,'O')
		CALL openUnit(itpL1Name,541,'O')
		!Read in the MPS
		CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
		CLOSE(540)
		CLOSE(541)
	ELSE
		PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
		CALL openUnit(itpG2Name,542,'O')
		CALL openUnit(itpL2Name,543,'O')
		!Read in the MPS
		CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
		CLOSE(542)
		CLOSE(543)
	ELSE
		PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap
fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)

CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jSave(1),totNum*1.0_rKind,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((J_Max-j_Min)/((j_res-1)*1.0_rKind))**2))

CLOSE(538)

	END DO
	END DO
END IF

END IF


!Have the master consolidate all files
IF(my_rank==0) THEN
!Set up consolidated filename

CALL createFileName(FinaloutputName,pdDir)
!NC switch decision 
IF(mu_Res==0) THEN
	CALL appendBaseName(FinaloutputName,'FANC_')
ELSE
	CALL appendBaseName(FinaloutputName,'FA_')
END IF

CALL appendBaseName(FinaloutputName,'L',systemSize)
CALL appendBaseName(FinaloutputName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(FinaloutputName,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(FinaloutputName,'N',4,mu_Min)
END IF
CALL appendBaseName(FinaloutputName,'t',4,mu_Max)
CALL appendBaseName(FinaloutputName,'jR',4,j_Min)
CALL appendBaseName(FinaloutputName,'t',4,j_Max)

!Copy i/o to all files
DO i=0,numactive-1
CALL copyName(FinaloutputName,outputName(i))
CALL appendBaseName(outputName(i),'rank',i)
CALL appendBaseName(outputName(i),'FixedMuOutput.dat')
END DO

CALL appendBaseName(FinaloutputName,'MuOutput.dat')

CALL openUnit(FinaloutputName,400)

DO i=0,numactive-1
myUnit=100+i
CALL openUnit(outputName(i),myUnit,'O')
!Start from the beginning of the file
REWIND(myUnit)
	ioSw1=0
	DO WHILE(ioSw1==0)
	READ(myUnit,'(4E30.15)',IOSTAT=ioSw1) jSave(1), jSave(2), muSave(1), muSave(2)
	IF(ioSw1.ne.0) EXIT
	WRITE(400,'(4E30.15)') jSave(1), jSave(2), muSave(1), muSave(2)
	END DO
	CLOSE(UNIT=myUnit, STATUS='DELETE')
END DO

CLOSE(400)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

!!!!!!!!!!!!!!FIXED J FIDELITY ANALYSIS!!!!!!!!!!!!!!!!!!

IF(num_cores.gt.j_Res) THEN
numactive=j_Res
ELSE
numactive=num_cores
END IF

IF(my_rank.le.numactive-1) THEN


IF(.NOT.ncSwitch) THEN


my_b(1)=1+my_rank*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(j_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

jTunn=j_Min+(i-1)*(j_max-j_min)/((j_Res-1)*1.0_rKind)


DO j=1,mu_Res-1

muSave(1)=mu_Min+(j-1)*(mu_Max-mu_Min)/((mu_res-1)*1.0_rKind)
muSave(2)=mu_Min+(j)*(mu_Max-mu_Min)/((mu_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'PD_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'mu',4,muSave(1))
CALL appendBaseName(itpG1Name,'j',4,jTunn)
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'PD_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'mu',4,muSave(2))
CALL appendBaseName(itpG2Name,'j',4,jTunn)
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap


fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)

CALL openUnit(outputjName,538,'A')
WRITE(538,'(4E30.15)') jTunn,muSave(1),ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((mu_Max-mu_Min)/((mu_res-1)*1.0_rKind))**2))

CLOSE(538)

END DO
END DO

ELSE

my_b(1)=1+my_rank*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(j_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

jTunn=j_Min+(i-1)*(j_max-j_min)/((j_Res-1)*1.0_rKind)


DO j=MIN(FLOOR(mu_min)-1,0),FLOOR(mu_Max)+1

numSave(1)=j
numSave(2)=j+1


!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'PDNC_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'N',numSave(1))
CALL appendBaseName(itpG1Name,'j',4,jTunn)
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'PDNC_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'N',numSave(1))
CALL appendBaseName(itpG2Name,'j',4,jTunn)
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap
fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)


CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jTunn,numSave(1)*1.0_rKind,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize)

CLOSE(538)

END DO
END DO


END IF
CALL DeallocateGamLam(Gammas1, Lambdas1)
CALL DeallocateGamLam(Gammas2, Lambdas2)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

!Have the master consolidate all files
IF(my_rank==0) THEN
!Set up consolidated filename

CALL createFileName(FinaloutputNameJ,pdDir)
!NC switch decision 
IF(mu_Res==0) THEN
	CALL appendBaseName(FinaloutputNameJ,'FANC_')
ELSE
	CALL appendBaseName(FinaloutputNameJ,'FA_')
END IF

CALL appendBaseName(FinaloutputNameJ,'L',systemSize)
CALL appendBaseName(FinaloutputNameJ,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(FinaloutputNameJ,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(FinaloutputNameJ,'N',4,mu_Min)
END IF
CALL appendBaseName(FinaloutputNameJ,'t',4,mu_Max)
CALL appendBaseName(FinaloutputNameJ,'jR',4,j_Min)
CALL appendBaseName(FinaloutputNameJ,'t',4,j_Max)

!Copy i/o to all files
DO i=0,numactive-1
CALL copyName(FinaloutputNameJ,outputName(i))
CALL appendBaseName(outputName(i),'rank',i)
CALL appendBaseName(outputName(i),'FixedJOutput.dat')
END DO

CALL appendBaseName(FinaloutputNameJ,'JOutput.dat')

CALL openUnit(FinaloutputNameJ,400)

DO i=0,numactive-1
myUnit=100+i
CALL openUnit(outputName(i),myUnit,'O')
!Start from the beginning of the file
REWIND(myUnit)
	ioSw1=0
	DO WHILE(ioSw1==0)
	READ(myUnit,'(4E30.15)',IOSTAT=ioSw1) jSave(1), jSave(2), muSave(1), muSave(2)
	IF(ioSw1.ne.0) EXIT
	WRITE(400,'(4E30.15)') jSave(1), jSave(2), muSave(1), muSave(2)
	END DO
	CLOSE(UNIT=myUnit, STATUS='DELETE')
END DO

CLOSE(400)
END IF

END SUBROUTINE FidelityAnalysis


SUBROUTINE FidelityAnalysisFermi(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
!
!Purpose: Read in Gammas and Lambdas from output files to compute
!Fidelities for Fermi code
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, j_Min, j_Max
INTEGER, INTENT(IN) :: mu_Res, j_Res
TYPE(tensor), POINTER :: Gammas1(:),Gammas2(:)
TYPE(vector), POINTER :: Lambdas1(:), Lambdas2(:)
INTEGER i,j,k, counter, ioSw1, ioSw2, counterVec(0:num_cores-1), myUnit
REAL(KIND=rKind) :: jSave(2), muSave(2), varSave(2),dummies(2)
INTEGER :: numSave(2), numactive, my_b(2)
REAL(KIND=rKind) :: mu_low, mu_high, mu_diff
REAL(KIND=rKind) :: obseVec(j_Res), jDummy
COMPLEX(KIND=rKind) :: fidel
CHARACTER(len=132) :: outputmuName,outputJName, FinaloutputName,FinaloutputNameJ
CHARACTER(len=132) :: itpG1name, itpG2name, itpL1name, itpL2name
CHARACTER(len=132) :: outputName(0:num_cores-1)

IF(num_cores.gt.mu_Res) THEN
numactive=mu_Res
ELSE
numactive=num_cores
END IF

IF(my_rank.le.numactive-1) THEN

CALL createFileName(outputmuName,pdDir)

!NC switch decision 
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
	IF((FLOOR(mu_Min).ne.CEILING(mu_Min)).OR.(FLOOR(mu_Max).ne.CEILING(mu_Max))) THEN
		PRINT *, 'Warning: The number conserving method is being used but either mu_Min or mu_Max is not integer.'
		PRINT *,  'The floor is being taken!'
		PRINT *, ''
	END IF
	CALL appendBaseName(outputmuName,'FANC_')
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
	CALL appendBaseName(outputmuName,'FA_')
END IF


CALL appendBaseName(outputmuName,'L',systemSize)
CALL appendBaseName(outputmuName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(outputmuName,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(outputmuName,'N',4,mu_Min)
END IF
CALL appendBaseName(outputmuName,'t',4,mu_Max)
CALL appendBaseName(outputmuName,'jR',4,j_Min)
CALL appendBaseName(outputmuName,'t',4,j_Max)
CALL appendBaseName(outputmuName,'rank',my_rank)

!Copy i/o to all files
CALL copyName(outputmuName,outputJName)
CALL appendBaseName(outputmuName,'FixedMuOutput.dat')
CALL appendBaseName(outputJName,'FixedJOutput.dat')

CALL openUnit(outputmuName,538)
CALL openUnit(outputJName,539)
CLOSE(538)
CLOSE(539)

!!!!!!!!!!!!!!FIXED MU FIDELITY ANALYSIS!!!!!!!!!!!!!!!!!1
CALL AllocateGamLam(Gammas1, Lambdas1, chiMax)
CALL AllocateGamLam(Gammas2, Lambdas2, chiMax)

IF(.NOT.ncSwitch) THEN

my_b(1)=1+my_rank*FLOOR(mu_Res*1.0_rKind/((numactive)*1.0_rKind))
my_b(2)=(my_rank+1)*FLOOR(mu_Res*1.0_rKind/((numactive)*1.0_rKind))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(mu_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

mu0=mu_Min+(i-1)*(mu_max-mu_min)/((mu_Res-1)*1.0_rKind)


DO j=1,j_Res-1

jSave(1)=j_Min+(j-1)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)
jSave(2)=j_Min+(j)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPD_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'mu',4,mu0)
CALL appendBaseName(itpG1Name,'j',4,jSave(1))
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')


CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPD_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'mu',4,mu0)
CALL appendBaseName(itpG2Name,'j',4,jSave(2))
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap


fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)

CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jSave(1),mu0,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((J_Max-j_Min)/((j_res-1)*1.0_rKind))**2))

CLOSE(538)

END DO
END DO

ELSE

mu_low=(Min(FLOOR(mu_Min)-1,0))*1.0_rKind
mu_high=(FLOOR(mu_Max)+1)*1.0_rKind
mu_diff=mu_low-mu_High+1

my_b(1)=1+my_rank*FLOOR((mu_diff*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((mu_diff*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(FLOOR(mu_diff),numactive)
END IF

DO i=my_b(1),my_b(2),1

totNum=i


DO j=1,j_Res-1

jSave(1)=j_Min+(j-1)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)
jSave(2)=j_Min+(j)*(J_Max-j_Min)/((j_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPDNC_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'N',totNum)
CALL appendBaseName(itpG1Name,'j',4,jSave(1))
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPDNC_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'N',totNum)
CALL appendBaseName(itpG2Name,'j',4,jSave(2))
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap
fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)


CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jSave(1),totNum*1.0_rKind,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((J_Max-j_Min)/((j_res-1)*1.0_rKind))**2))

CLOSE(538)

END DO
END DO


END IF

END IF


!Have the master consolidate all files
IF(my_rank==0) THEN
!Set up consolidated filename

CALL createFileName(FinaloutputName,pdDir)
!NC switch decision 
IF(mu_Res==0) THEN
	CALL appendBaseName(FinaloutputName,'FANC_')
ELSE
	CALL appendBaseName(FinaloutputName,'FA_')
END IF

CALL appendBaseName(FinaloutputName,'L',systemSize)
CALL appendBaseName(FinaloutputName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(FinaloutputName,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(FinaloutputName,'N',4,mu_Min)
END IF
CALL appendBaseName(FinaloutputName,'t',4,mu_Max)
CALL appendBaseName(FinaloutputName,'jR',4,j_Min)
CALL appendBaseName(FinaloutputName,'t',4,j_Max)

!Copy i/o to all files
DO i=0,numactive-1
CALL copyName(FinaloutputName,outputName(i))
CALL appendBaseName(outputName(i),'rank',i)
CALL appendBaseName(outputName(i),'FixedMuOutput.dat')
END DO

CALL appendBaseName(FinaloutputName,'MuOutput.dat')

CALL openUnit(FinaloutputName,400)

DO i=0,numactive-1
myUnit=100+i
CALL openUnit(outputName(i),myUnit,'O')
!Start from the beginning of the file
REWIND(myUnit)
	ioSw1=0
	DO WHILE(ioSw1==0)
	READ(myUnit,'(4E30.15)',IOSTAT=ioSw1) jSave(1), jSave(2), muSave(1), muSave(2)
	IF(ioSw1.ne.0) EXIT
	WRITE(400,'(4E30.15)') jSave(1), jSave(2), muSave(1), muSave(2)
	END DO
	CLOSE(UNIT=myUnit, STATUS='DELETE')
END DO

CLOSE(400)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
!!!!!!!!!!!!!!FIXED J FIDELITY ANALYSIS!!!!!!!!!!!!!!!!!1

IF(num_cores.gt.j_Res) THEN
numactive=j_Res
ELSE
numactive=num_cores
END IF

IF(my_rank.le.numactive-1) THEN


IF(.NOT.ncSwitch) THEN


my_b(1)=1+my_rank*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(j_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

jTunn=j_Min+(i-1)*(j_max-j_min)/((j_Res-1)*1.0_rKind)


DO j=1,mu_Res-1

muSave(1)=mu_Min+(j-1)*(mu_Max-mu_Min)/((mu_res-1)*1.0_rKind)
muSave(2)=mu_Min+(j)*(mu_Max-mu_Min)/((mu_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPD_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'mu',4,muSave(1))
CALL appendBaseName(itpG1Name,'j',4,jTunn)
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPD_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'mu',4,muSave(2))
CALL appendBaseName(itpG2Name,'j',4,jTunn)
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap


fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)

CALL openUnit(outputjName,538,'A')
WRITE(538,'(4E30.15)') jTunn,muSave(1),ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((mu_Max-mu_Min)/((mu_res-1)*1.0_rKind))**2))

CLOSE(538)

END DO
END DO

ELSE

my_b(1)=1+my_rank*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((j_Res*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(j_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

jTunn=j_Min+(i-1)*(j_max-j_min)/((j_Res-1)*1.0_rKind)


DO j=MIN(FLOOR(mu_min)-1,0),FLOOR(mu_Max)+1

numSave(1)=j
numSave(2)=j+1


!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPDNC_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'N',numSave(1))
CALL appendBaseName(itpG1Name,'j',4,jTunn)
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPDNC_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'N',numSave(1))
CALL appendBaseName(itpG2Name,'j',4,jTunn)
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap
fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)


CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jTunn,numSave(1)*1.0_rKind,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize)

CLOSE(538)

END DO
END DO


END IF
CALL DeallocateGamLam(Gammas1, Lambdas1)
CALL DeallocateGamLam(Gammas2, Lambdas2)

END IF



!!!!!!!!!!!!!!!POSTPROCESSING!!!!!!!!!!!!!!!!!!!!!
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

!Have the master consolidate all files
IF(my_rank==0) THEN
!Set up consolidated filename

CALL createFileName(FinaloutputNameJ,pdDir)
!NC switch decision 
IF(mu_Res==0) THEN
	CALL appendBaseName(FinaloutputNameJ,'FANC_')
ELSE
	CALL appendBaseName(FinaloutputNameJ,'FA_')
END IF

CALL appendBaseName(FinaloutputNameJ,'L',systemSize)
CALL appendBaseName(FinaloutputNameJ,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(FinaloutputNameJ,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(FinaloutputNameJ,'N',4,mu_Min)
END IF
CALL appendBaseName(FinaloutputNameJ,'t',4,mu_Max)
CALL appendBaseName(FinaloutputNameJ,'jR',4,j_Min)
CALL appendBaseName(FinaloutputNameJ,'t',4,j_Max)

!Copy i/o to all files
DO i=0,numactive-1
CALL copyName(FinaloutputNameJ,outputName(i))
CALL appendBaseName(outputName(i),'rank',i)
CALL appendBaseName(outputName(i),'FixedJOutput.dat')
END DO

CALL appendBaseName(FinaloutputNameJ,'JOutput.dat')

CALL openUnit(FinaloutputNameJ,400)

DO i=0,numactive-1
myUnit=100+i
CALL openUnit(outputName(i),myUnit,'O')
!Start from the beginning of the file
REWIND(myUnit)
	ioSw1=0
	DO WHILE(ioSw1==0)
	READ(myUnit,'(4E30.15)',IOSTAT=ioSw1) jSave(1), jSave(2), muSave(1), muSave(2)
	IF(ioSw1.ne.0) EXIT
	WRITE(400,'(4E30.15)') jSave(1), jSave(2), muSave(1), muSave(2)
	END DO
	CLOSE(UNIT=myUnit, STATUS='DELETE')
END DO

CLOSE(400)
END IF

END SUBROUTINE FidelityAnalysisFermi



SUBROUTINE FidelityAnalysisFermiFromU(mu_Min, mu_Max, mu_Res, U_Min, U_Max, U_Res)
!
!Purpose: Read in Gammas and Lambdas from output files to compute
!Fidelities for Fermi code
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, U_Min, U_Max
INTEGER, INTENT(IN) :: mu_Res, U_Res
TYPE(tensor), POINTER :: Gammas1(:),Gammas2(:)
TYPE(vector), POINTER :: Lambdas1(:), Lambdas2(:)
INTEGER i,j,k, counter, ioSw1, ioSw2, counterVec(0:num_cores-1), myUnit
REAL(KIND=rKind) :: jSave(2), muSave(2), varSave(2),dummies(2)
INTEGER :: numSave(2), numactive, my_b(2)
REAL(KIND=rKind) :: mu_low, mu_high, mu_diff
REAL(KIND=rKind) :: obseVec(U_Res), jDummy
COMPLEX(KIND=rKind) :: fidel
CHARACTER(len=132) :: outputmuName,outputJName, FinaloutputName,FinaloutputNameJ
CHARACTER(len=132) :: itpG1name, itpG2name, itpL1name, itpL2name
CHARACTER(len=132) :: outputName(0:num_cores-1)

IF(num_cores.gt.mu_Res) THEN
numactive=mu_Res
ELSE
numactive=num_cores
END IF

IF(my_rank.le.numactive-1) THEN

CALL createFileName(outputmuName,pdDir)

!NC switch decision 
IF(mu_Res==0) THEN
	ncSwitch=.TRUE.
	IF(print_switch) THEN
		PRINT *, 'Number conserving method chosen by 0 mu resolution!'
	END IF
	IF((FLOOR(mu_Min).ne.CEILING(mu_Min)).OR.(FLOOR(mu_Max).ne.CEILING(mu_Max))) THEN
		PRINT *, 'Warning: The number conserving method is being used but either mu_Min or mu_Max is not integer.'
		PRINT *,  'The floor is being taken!'
		PRINT *, ''
	END IF
	CALL appendBaseName(outputmuName,'FANC_')
ELSE
	ncSwitch=.FALSE.
	IF(print_switch) THEN
		PRINT *, 'Number non-conserving method chosen by nonzero mu resolution!'
	END IF
	CALL appendBaseName(outputmuName,'FA_')
END IF


CALL appendBaseName(outputmuName,'L',systemSize)
CALL appendBaseName(outputmuName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(outputmuName,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(outputmuName,'N',4,mu_Min)
END IF
CALL appendBaseName(outputmuName,'t',4,mu_Max)
CALL appendBaseName(outputmuName,'uR',4,U_Min)
CALL appendBaseName(outputmuName,'t',4,U_Max)
CALL appendBaseName(outputmuName,'rank',my_rank)

!Copy i/o to all files
CALL copyName(outputmuName,outputJName)
CALL appendBaseName(outputmuName,'FixedMuOutput.dat')
CALL appendBaseName(outputJName,'FixedUOutput.dat')

CALL openUnit(outputmuName,538)
CALL openUnit(outputJName,539)
CLOSE(538)
CLOSE(539)

!!!!!!!!!!!!!!FIXED MU FIDELITY ANALYSIS!!!!!!!!!!!!!!!!!1
CALL AllocateGamLam(Gammas1, Lambdas1, chiMax)
CALL AllocateGamLam(Gammas2, Lambdas2, chiMax)

IF(.NOT.ncSwitch) THEN

my_b(1)=1+my_rank*FLOOR(mu_Res*1.0_rKind/((numactive)*1.0_rKind))
my_b(2)=(my_rank+1)*FLOOR(mu_Res*1.0_rKind/((numactive)*1.0_rKind))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(mu_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

mu0=mu_Min+(i-1)*(mu_max-mu_min)/((mu_Res-1)*1.0_rKind)


DO j=1,U_Res-1

jSave(1)=U_Min+(j-1)*(U_Max-U_Min)/((U_res-1)*1.0_rKind)
jSave(2)=U_Min+(j)*(U_Max-U_Min)/((U_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPD_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'mu',4,mu0)
CALL appendBaseName(itpG1Name,'U',4,jSave(1))
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')


CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPD_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'mu',4,mu0)
CALL appendBaseName(itpG2Name,'U',4,jSave(2))
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap


fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)

CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jSave(1),mu0,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((U_Max-U_Min)/((U_res-1)*1.0_rKind))**2))

CLOSE(538)

END DO
END DO

ELSE

mu_low=(Min(FLOOR(mu_Min)-1,0))*1.0_rKind
mu_high=(FLOOR(mu_Max)+1)*1.0_rKind
mu_diff=mu_low-mu_High+1

my_b(1)=1+my_rank*FLOOR((mu_diff*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((mu_diff*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(FLOOR(mu_diff),numactive)
END IF

DO i=my_b(1),my_b(2),1

totNum=i


DO j=1,U_Res-1

jSave(1)=U_Min+(j-1)*(U_Max-U_Min)/((U_res-1)*1.0_rKind)
jSave(2)=U_Min+(j)*(U_Max-U_Min)/((U_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPDNC_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'N',totNum)
CALL appendBaseName(itpG1Name,'U',4,jSave(1))
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPDNC_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'N',totNum)
CALL appendBaseName(itpG2Name,'U',4,jSave(2))
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap
fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)


CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jSave(1),totNum*1.0_rKind,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((U_Max-U_Min)/((U_res-1)*1.0_rKind))**2))

CLOSE(538)

END DO
END DO


END IF

END IF


!Have the master consolidate all files
IF(my_rank==0) THEN
!Set up consolidated filename

CALL createFileName(FinaloutputName,pdDir)
!NC switch decision 
IF(mu_Res==0) THEN
	CALL appendBaseName(FinaloutputName,'FANC_')
ELSE
	CALL appendBaseName(FinaloutputName,'FA_')
END IF

CALL appendBaseName(FinaloutputName,'L',systemSize)
CALL appendBaseName(FinaloutputName,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(FinaloutputName,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(FinaloutputName,'N',4,mu_Min)
END IF
CALL appendBaseName(FinaloutputName,'t',4,mu_Max)
CALL appendBaseName(FinaloutputName,'uR',4,U_Min)
CALL appendBaseName(FinaloutputName,'t',4,U_Max)

!Copy i/o to all files
DO i=0,numactive-1
CALL copyName(FinaloutputName,outputName(i))
CALL appendBaseName(outputName(i),'rank',i)
CALL appendBaseName(outputName(i),'FixedMuOutput.dat')
END DO

CALL appendBaseName(FinaloutputName,'MuOutput.dat')

CALL openUnit(FinaloutputName,400)

DO i=0,numactive-1
myUnit=100+i
CALL openUnit(outputName(i),myUnit,'O')
!Start from the beginning of the file
REWIND(myUnit)
	ioSw1=0
	DO WHILE(ioSw1==0)
	READ(myUnit,'(4E30.15)',IOSTAT=ioSw1) jSave(1), jSave(2), muSave(1), muSave(2)
	IF(ioSw1.ne.0) EXIT
	WRITE(400,'(4E30.15)') jSave(1), jSave(2), muSave(1), muSave(2)
	END DO
	CLOSE(UNIT=myUnit, STATUS='DELETE')
END DO

CLOSE(400)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
!!!!!!!!!!!!!!FIXED U FIDELITY ANALYSIS!!!!!!!!!!!!!!!!!1

IF(num_cores.gt.U_Res) THEN
numactive=U_Res
ELSE
numactive=num_cores
END IF

IF(my_rank.le.numactive-1) THEN


IF(.NOT.ncSwitch) THEN


my_b(1)=1+my_rank*FLOOR((U_Res*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((U_Res*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(U_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

jTunn=U_Min+(i-1)*(U_max-U_min)/((U_Res-1)*1.0_rKind)


DO j=1,mu_Res-1

muSave(1)=mu_Min+(j-1)*(mu_Max-mu_Min)/((mu_res-1)*1.0_rKind)
muSave(2)=mu_Min+(j)*(mu_Max-mu_Min)/((mu_res-1)*1.0_rKind)

!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPD_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'mu',4,muSave(1))
CALL appendBaseName(itpG1Name,'U',4,jTunn)
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPD_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'mu',4,muSave(2))
CALL appendBaseName(itpG2Name,'U',4,jTunn)
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap


fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)

CALL openUnit(outputjName,538,'A')
WRITE(538,'(4E30.15)') jTunn,muSave(1),ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize*(((mu_Max-mu_Min)/((mu_res-1)*1.0_rKind))**2))

CLOSE(538)

END DO
END DO

ELSE

my_b(1)=1+my_rank*FLOOR((U_Res*1.0_rKind/((numactive)*1.0_rKind)))
my_b(2)=(my_rank+1)*FLOOR((U_Res*1.0_rKind/((numactive)*1.0_rKind)))

IF(my_rank==numactive-1) THEN
my_b(2)=my_b(2)+MOD(U_Res,numactive)
END IF

DO i=my_b(1),my_b(2),1

jTunn=U_Min+(i-1)*(U_max-U_min)/((U_Res-1)*1.0_rKind)


DO j=MIN(FLOOR(mu_min)-1,0),FLOOR(mu_Max)+1

numSave(1)=j
numSave(2)=j+1


!Create the file names corresponding to the given parameters
CALL createFileName(itpG1Name,itpDir)
CALL appendBaseName(itpG1Name,'FHPDNC_')
CALL appendBaseName(itpG1Name,'L',systemSize)
CALL appendBaseName(itpG1Name,'Chi',chiMax)
CALL appendBaseName(itpG1Name,'N',numSave(1))
CALL appendBaseName(itpG1Name,'U',4,U0)
CALL copyName(itpG1Name,itpL1Name)
CALL appendBaseName(itpG1Name,'itpG.dat')
CALL appendBaseName(itpL1Name,'itpL.dat')

CALL createFileName(itpG2Name,itpDir)
CALL appendBaseName(itpG2Name,'FHPDNC_')
CALL appendBaseName(itpG2Name,'L',systemSize)
CALL appendBaseName(itpG2Name,'Chi',chiMax)
CALL appendBaseName(itpG2Name,'N',numSave(1))
CALL appendBaseName(itpG2Name,'U',4,U0)
CALL copyName(itpG2Name,itpL2Name)
CALL appendBaseName(itpG2Name,'itpG.dat')
CALL appendBaseName(itpL2Name,'itpL.dat')

!Check to see if a state exists
	IF(CheckName(itpG1Name).AND.CheckName(itpL1Name)) THEN
	CALL openUnit(itpG1Name,540,'O')
	CALL openUnit(itpL1Name,541,'O')
	!Read in the MPS
	CALL readGammaLambda(541, 540,Gammas1,Lambdas1, ITPopenKind, chiMax)
CLOSE(540)
CLOSE(541)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Check to see if a state exists
	IF(CheckName(itpG2Name).AND.CheckName(itpL2Name)) THEN
	CALL openUnit(itpG2Name,542,'O')
	CALL openUnit(itpL2Name,543,'O')
	!Read in the MPS
	CALL readGammaLambda(543, 542,Gammas2,Lambdas2, ITPopenKind, chiMax)
CLOSE(542)
CLOSE(543)
	ELSE
	PRINT *, 'Warning: MPS not found!'
	END IF

!Compute the overlap
fidel=InnerProduct(Gammas1, Lambdas1, Gammas2, Lambdas2)


CALL openUnit(outputmuName,538,'A')
WRITE(538,'(4E30.15)') jTunn,numSave(1)*1.0_rKind,ABS(fidel)**2,2.0_rKind*(1.0_rKind-ABS(fidel)**2)/(systemSize)

CLOSE(538)

END DO
END DO


END IF
CALL DeallocateGamLam(Gammas1, Lambdas1)
CALL DeallocateGamLam(Gammas2, Lambdas2)

END IF



!!!!!!!!!!!!!!!POSTPROCESSING!!!!!!!!!!!!!!!!!!!!!
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

!Have the master consolidate all files
IF(my_rank==0) THEN
!Set up consolidated filename

CALL createFileName(FinaloutputNameJ,pdDir)
!NC switch decision 
IF(mu_Res==0) THEN
	CALL appendBaseName(FinaloutputNameJ,'FANC_')
ELSE
	CALL appendBaseName(FinaloutputNameJ,'FA_')
END IF

CALL appendBaseName(FinaloutputNameJ,'L',systemSize)
CALL appendBaseName(FinaloutputNameJ,'Chi',chiMax)
IF(.NOT.ncswitch) THEN
	CALL appendBaseName(FinaloutputNameJ,'muR',4,mu_Min)
ELSE
	CALL appendBaseName(FinaloutputNameJ,'N',4,mu_Min)
END IF
CALL appendBaseName(FinaloutputNameJ,'t',4,mu_Max)
CALL appendBaseName(FinaloutputNameJ,'uR',4,U_Min)
CALL appendBaseName(FinaloutputNameJ,'t',4,U_Max)

!Copy i/o to all files
DO i=0,numactive-1
CALL copyName(FinaloutputNameJ,outputName(i))
CALL appendBaseName(outputName(i),'rank',i)
CALL appendBaseName(outputName(i),'FixedUOutput.dat')
END DO

CALL appendBaseName(FinaloutputNameJ,'UOutput.dat')

CALL openUnit(FinaloutputNameJ,400)

DO i=0,numactive-1
myUnit=100+i
CALL openUnit(outputName(i),myUnit,'O')
!Start from the beginning of the file
REWIND(myUnit)
	ioSw1=0
	DO WHILE(ioSw1==0)
	READ(myUnit,'(4E30.15)',IOSTAT=ioSw1) jSave(1), jSave(2), muSave(1), muSave(2)
	IF(ioSw1.ne.0) EXIT
	WRITE(400,'(4E30.15)') jSave(1), jSave(2), muSave(1), muSave(2)
	END DO
	CLOSE(UNIT=myUnit, STATUS='DELETE')
END DO

CLOSE(400)
END IF

END SUBROUTINE FidelityAnalysisFermiFromU


SUBROUTINE PhaseDiagramBoseHubbard(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
!
!Purpose: Compute the phase diagram of the Bose hubbard model varying mu and J while keeping U fixed.
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, j_Min, j_Max
INTEGER, INTENT(IN) :: mu_Res, j_Res


U0=1.0_rKind
V0=0.0_rKind

!NC switch decision 
IF(mu_Res==0) THEN
	IF(num_cores.gt.j_Res*(FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1)) THEN
	PRINT *, "Number of cores is greater than number of tasks!"
	PRINT *, "Program is exiting!"
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
ELSE
	IF(num_cores.gt.j_Res*mu_Res) THEN
	PRINT *, "Number of cores is greater than number of tasks!"
	PRINT *, "Program is exiting!"
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
END IF

IF(my_rank==0) THEN
CALL PhaseDiagramMaster(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
ELSE
CALL PhaseDiagramBoseHubbardWorker(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
END IF


END SUBROUTINE PhaseDiagramBoseHubbard

SUBROUTINE PhaseDiagramFermiHubbard(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
!
!Purpose: Compute the phase diagram of the Fermi hubbard model varying mu and J while keeping U fixed.
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, j_Min, j_Max
INTEGER, INTENT(IN) :: mu_Res, j_Res


U0=1.0_rKind

!NC switch decision 
IF(mu_Res==0) THEN
	IF(num_cores.gt.j_Res*(FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1)) THEN
	PRINT *, "Number of cores is greater than number of tasks!"
	PRINT *, "Program is exiting!"
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
ELSE
	IF(num_cores.gt.j_Res*mu_Res) THEN
	PRINT *, "Number of cores is greater than number of tasks!"
	PRINT *, "Program is exiting!"
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
END IF

IF(my_rank==0) THEN
CALL PhaseDiagramMaster(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
ELSE
CALL PhaseDiagramFermiHubbardWorker(mu_Min, mu_Max, mu_Res, j_Min, j_Max, j_Res)
END IF


END SUBROUTINE PhaseDiagramFermiHubbard

SUBROUTINE PhaseDiagramFermiHubbardFromU(mu_Min, mu_Max, mu_Res, U_Min, U_Max, U_Res)
!
!Purpose: Compute the phase diagram of the Fermi hubbard model varying mu and U while keeping j fixed.
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: mu_Min, mu_Max, U_Min, U_Max
INTEGER, INTENT(IN) :: mu_Res, U_Res


jTunn=1.0_rKind

!NC switch decision 
IF(mu_Res==0) THEN
	IF(num_cores.gt.U_Res*(FLOOR(mu_Max+1)-MAX(FLOOR(mu_Min)-1,0)+1)) THEN
	PRINT *, "Number of cores is greater than number of tasks!"
	PRINT *, "Program is exiting!"
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
ELSE
	IF(num_cores.gt.U_Res*mu_Res) THEN
	PRINT *, "Number of cores is greater than number of tasks!"
	PRINT *, "Program is exiting!"
	CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
	END IF
END IF

IF(my_rank==0) THEN
CALL PhaseDiagramMasterFromU(mu_Min, mu_Max, mu_Res, U_Min, U_Max, U_Res)
ELSE
CALL PhaseDiagramFermiHubbardWorkerFromU(mu_Min, mu_Max, mu_Res, U_Min, U_Max, U_Res)
END IF


END SUBROUTINE PhaseDiagramFermiHubbardFromU

END MODULE PDtools_module
