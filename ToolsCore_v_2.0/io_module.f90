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
MODULE io_module
!
! Purpose: Module to perform basic i/o operations
! for OpenSourceTEBD v1.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   9/28/09   M. L. Wall	v2.0 release
!
USE system_parameters
USE TEBDtools_module
IMPLICIT NONE

INTERFACE appendBaseName
MODULE PROCEDURE appendBaseName_r, appendBaseName_i, appendBaseName_c
END INTERFACE appendBaseName

INTERFACE RecordOp
MODULE PROCEDURE RecordOp_m,RecordOp_mr,RecordOp_c,RecordOp_r
END INTERFACE RecordOp

INTERFACE RecordOpList
MODULE PROCEDURE RecordOpList_m,RecordOpList_mr
END INTERFACE RecordOpList

INTERFACE RecordOneSiteOb
MODULE PROCEDURE RecordOneSiteOb_r,RecordOneSiteOb_c
END INTERFACE RecordOneSiteOb

INTERFACE RecordTwoSiteOb
MODULE PROCEDURE RecordTwoSiteOb_r,RecordTwoSiteOb_c
END INTERFACE RecordTwoSiteOb

CONTAINS

SUBROUTINE createFileName(basename,diRectory)
!
!Purpose: Begin a file name in the directory diRectory
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: diRectory

baseName=diRectory

END SUBROUTINE createFileName

SUBROUTINE appendBaseName_r(basename,partName,partDigs,partValue)
!
!Purpose: Append to a file name the character string partName followed by 
!the value partValue to partDigs digits for real numbers
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: partname
INTEGER, INTENT(IN) :: partDigs
REAL(KIND=rKind), INTENT(IN) :: partValue
CHARACTER(16) :: iWstring, specString

!Write the number of decimal places wanted into a string
	WRITE(iWstring,'(I4)') partDigs
!Generate the format type for a float string with partDigs decimal places	
	specString="(1f16."//TRIM(ADJUSTL(iWstring))//")"
!Write the value partValue to a string using partDigs decimal places
	WRITE(iWstring,specString) partValue
!append partname and the value to the given basename
basename=TRIM(basename)//partName//TRIM(ADJUSTL(iWstring))

END SUBROUTINE appendBaseName_r

SUBROUTINE appendBaseName_i(basename,partName,partValue)
!
!Purpose: Append to a file name the character string partName followed by 
!the value partValue to partDigs digits for integers
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: partname
INTEGER, INTENT(IN) :: partValue
CHARACTER(16) :: iWstring, specString

!Write the number of digits wanted into a string
	WRITE(iWstring,'(I4)') 16
!Generate the format type for an integer string with partDigs digits	
	specString="(I"//TRIM(ADJUSTL(iWstring))//")"
!Write the value partValue to a string using partDigs digits
	WRITE(iWstring,specString) partValue
!append partname and the value to the given basename
basename=TRIM(basename)//partName//TRIM(ADJUSTL(iWstring))

END SUBROUTINE appendBaseName_i

SUBROUTINE appendBaseName_c(basename,partName)
!
!Purpose: Append to a file name the character string partName
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(INOUT) :: baseName
CHARACTER(len=*), INTENT(IN) :: partname
CHARACTER(16) :: iWstring, specString

basename=TRIM(basename)//TRIM(partName)

END SUBROUTINE appendBaseName_c

SUBROUTINE copyName(name1,name2)
!
!Purpose:Copy name1 to name2
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: name1
CHARACTER(len=*), INTENT(OUT) :: name2

name2=TRIM(name1)

END SUBROUTINE copyName

LOGICAL FUNCTION CheckName(baseName)
!
!Purpose: Returns TRUE if file name baseName exists and FALSE otherwise
!
IMPLICIT NONE
CHARACTER(len=*), INTENT(IN) :: baseName

INQUIRE(FILE=baseName, EXIST=CheckName)

END FUNCTION CheckName

SUBROUTINE openUnit(fileName,myUnit,openKind)
!
!Purpose: Open a file 'filename' as UNIT=myUnit
!An error message is returned if the specified file cannot be opened
!
IMPLICIT NONE
INTEGER, INTENT(IN) ::  myUnit
CHARACTER(len=*), INTENT(IN) :: fileName
CHARACTER(132) ::  stopname
CHARACTER, INTENT(IN), OPTIONAL :: openKind

	stopname='*** Cannot open file named '//filename//'***'

IF(PRESENT(openKind)) THEN
IF(openKind=='N') THEN
	OPEN(UNIT=myUnit, FILE=filename, STATUS='NEW', ACTION='WRITE',IOSTAT=fileStatus)
		IF(fileStatus>0) THEN
			PRINT *,stopname
			STOP
		END IF
ELSE IF(openKind=='O') THEN
	OPEN(UNIT=myUnit, FILE=filename, STATUS='OLD', ACTION='READWRITE',IOSTAT=fileStatus)
		IF(fileStatus>0) THEN
			PRINT *,stopname
			STOP
		END IF
ELSE IF(openKind=='A') THEN
	OPEN(UNIT=myUnit, FILE=filename, ACTION='WRITE',POSITION='APPEND',IOSTAT=fileStatus)
		IF(fileStatus>0) THEN
			PRINT *,stopname
			STOP
		END IF
ELSE
STOP "Unknown option in openUnit!"
END IF

ELSE
	OPEN(UNIT=myUnit, FILE=filename, STATUS='UNKNOWN', ACTION='READWRITE',IOSTAT=fileStatus)
		IF(fileStatus>0) THEN
			PRINT *,stopname
			STOP
		END IF
END IF

END SUBROUTINE openUnit

SUBROUTINE RecordLambdas(fileid, Lambdas,openKind)
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

	DO i=1,(systemSize+1)
		WRITE(UNIT=fileid, FMT=fmtname) (Lambdas(i)%v(j), j=1,chi)
	END DO
END SUBROUTINE RecordLambdas

SUBROUTINE RecordGammas(fileid, Gammas, openKind)
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

	DO l=1,systemSize
		DO i=1,chi
			DO j=1,chi
				WRITE(UNIT=fileid, FMT=fmtName) (REAL(Gammas(l)%t(i,k,j)), k=1,localSize)
			END DO
		END DO
	END DO
	
	DO l=1,systemSize
		DO i=1,chi
			DO j=1,chi
				WRITE(UNIT=fileid, FMT=fmtName) (AIMAG(Gammas(l)%t(i,k,j)), k=1,localSize)
			END DO
		END DO
	END DO
		
END SUBROUTINE RecordGammas


SUBROUTINE RecordLabel(fileid, LabelLorR)
!
!Purpose: Record Label on a file whose ID is fileid.
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

	DO i=1,(systemSize+1)
		WRITE(UNIT=fileid, FMT=fmtname) (LabelLorR(i)%vi(j), j=1,chi)
	END DO
END SUBROUTINE RecordLabel

SUBROUTINE readGammaLambda(lambdafileID, gammafileID,Gammas,Lambdas, openKind, chiNow)
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
	    DO i=1,(systemSize+1)
	    	READ(UNIT=lambdaFileID, FMT=fmtNamel) (Lambdas(i)%v(j), j=1,chiNow)
	    END DO

		ALLOCATE(gammatemp(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate gammatemp'
			END IF

! Read real part of Gammas
		DO l=1,systemSize
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
		DO l=1,systemSize
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
			
END SUBROUTINE readGammaLambda


SUBROUTINE readGammaLambdaLabels(lambdafileID, gammafileID,labelleftFileID, labelrightFileID,&
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
	    DO i=1,(systemSize+1)
	    	READ(UNIT=lambdaFileID, FMT=fmtNamel) (Lambdas(i)%v(j), j=1,chiNow)
	    END DO

		ALLOCATE(gammatemp(localSize), STAT=statInt)
			IF(statInt.ne.0) THEN
			PRINT *, 'Failed to allocate gammatemp'
			END IF

! Read real part of Gammas
		DO l=1,systemSize
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
		DO l=1,systemSize
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
		DO i=1,(systemSize+1)
	    	READ(UNIT=labelleftFileID, FMT=fmtNameLabel) (LabelLeft(i)%vi(j), j=1,chiNow)
	    END DO
! Read LabelRight
		DO i=1,(systemSize+1)
	    	READ(UNIT=labelrightFileID, FMT=fmtNameLabel) (LabelRight(i)%vi(j), j=1,chiNow)
	    END DO
			
END SUBROUTINE readGammaLambdaLabels

SUBROUTINE RecordOp_m(fileid, Op)
!
!Purpose: Record a matrix operator on a file whose ID is fileid.	
!
IMPLICIT NONE
TYPE(matrix) :: Op
INTEGER, INTENT(IN) :: fileid
INTEGER :: i,j,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=SIZE(Op%m,2)
WRITE(mString,'(I4)') 2*dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
		
	DO i=1,dd
		WRITE(UNIT=fileid, FMT=fmtName) (REAL(Op%m(i,j)), j=1,dd), (AIMAG(Op%m(i,j)), j=1,dd)
	END DO
	
END SUBROUTINE RecordOp_m

SUBROUTINE RecordOp_mr(fileid, Op)
!
!Purpose: Record a matrixReal operator on a file whose ID is fileid.	
!
IMPLICIT NONE
TYPE(matrixReal) :: Op
INTEGER, INTENT(IN) :: fileid
INTEGER :: i,j,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=SIZE(Op%mr,2)
WRITE(mString,'(I4)') dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'

		
	DO i=1,dd
		WRITE(UNIT=fileid, FMT=fmtName) (Op%mr(i,j), j=1,dd)
	END DO
	
END SUBROUTINE RecordOp_mr

SUBROUTINE RecordOp_c(fileid, Op)
!
!Purpose: Record a matrix operator on a file whose ID is fileid.	
!
IMPLICIT NONE
COMPLEX(KIND=rKIND), INTENT(IN) :: Op(:,:)
INTEGER, INTENT(IN) :: fileid
INTEGER :: i,j,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=SIZE(Op,2)
WRITE(mString,'(I4)') 2*dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
		
	DO i=1,dd
		WRITE(UNIT=fileid, FMT=fmtName) (REAL(Op(i,j)), j=1,dd), (AIMAG(Op(i,j)), j=1,dd)
	END DO
	
END SUBROUTINE RecordOp_c

SUBROUTINE RecordOp_r(fileid, Op)
!
!Purpose: Record a matrixReal operator on a file whose ID is fileid.	
!
IMPLICIT NONE
REAL(KIND=rKIND), INTENT(IN) :: Op(:,:)
INTEGER, INTENT(IN) :: fileid
INTEGER :: i,j,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=SIZE(Op,2)
WRITE(mString,'(I4)') dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'

		
	DO i=1,dd
		WRITE(UNIT=fileid, FMT=fmtName) (Op(i,j), j=1,dd)
	END DO
	
END SUBROUTINE RecordOp_r

SUBROUTINE RecordOpList_m(fileid, Op)
!
!Purpose: Record a matrix operator List on a file whose ID is fileid.	
!
IMPLICIT NONE
TYPE(matrix), POINTER :: Op(:)
INTEGER, INTENT(IN) :: fileid
INTEGER :: i,j,k,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=SIZE(Op(1)%m,2)
WRITE(mString,'(I4)') 2*dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	
	DO k=1,SIZE(Op,1),1		
		DO i=1,dd
			WRITE(UNIT=fileid, FMT=fmtName) (REAL(Op(k)%m(i,j)), j=1,dd), (AIMAG(Op(k)%m(i,j)), j=1,dd)
		END DO
	END DO
	
END SUBROUTINE RecordOpList_m

SUBROUTINE RecordOpList_mr(fileid, Op)
!
!Purpose: Record a matrixReal operator list on a file whose ID is fileid.	
!
IMPLICIT NONE
TYPE(matrixReal), POINTER :: Op(:)
INTEGER, INTENT(IN) :: fileid
INTEGER :: i,j,k,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=SIZE(Op(1)%mr,2)
WRITE(mString,'(I4)') dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	
	DO k=1,SIZE(Op,1),1		
		DO i=1,dd
			WRITE(UNIT=fileid, FMT=fmtName) (REAL(Op(k)%mr(i,j)), j=1,dd)
		END DO
	END DO
	
END SUBROUTINE RecordOpList_mr

!!!!!!!!!!!!!BEGIN INTERFACE RecordOneSiteOb!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RecordOneSiteOb_r(fileid, Obs, timeN)
!
!Purpose: Record a matrix operator List on a file whose ID is fileid.	
!
IMPLICIT NONE
REAL(KIND=rKind) :: Obs(systemSize)
REAL(KIND=rKind), OPTIONAL :: timeN
INTEGER, INTENT(IN) :: fileid
INTEGER :: j,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=systemSize
IF(PRESENT(timeN)) THEN
WRITE(mString,'(I4)') dd+1
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	
			WRITE(UNIT=fileid, FMT=fmtName) timeN, (Obs(j), j=1,dd)
ELSE
WRITE(mString,'(I4)') dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	
			WRITE(UNIT=fileid, FMT=fmtName) (Obs(j), j=1,dd)
END IF
	
END SUBROUTINE RecordOneSiteOb_r

SUBROUTINE RecordOneSiteOb_c(fileid, Obs, timeN)
!
!Purpose: Record a matrix operator List on a file whose ID is fileid.	
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: Obs(systemSize)
REAL(KIND=rKind), OPTIONAL :: timeN
INTEGER, INTENT(IN) :: fileid
INTEGER :: j,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=systemSize
IF(PRESENT(timeN)) THEN
WRITE(mString,'(I4)') 2*dd+1
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	
					WRITE(UNIT=fileid, FMT=fmtName) timeN, (REAL(Obs(j)), j=1,dd), (AIMAG(Obs(j)), j=1,dd)
ELSE
WRITE(mString,'(I4)') 2*dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'
	
					WRITE(UNIT=fileid, FMT=fmtName) (REAL(Obs(j)), j=1,dd), (AIMAG(Obs(j)), j=1,dd)
END IF
	
END SUBROUTINE RecordOneSiteOb_c

!!!!!!!!!!!!!END INTERFACE RecordOneSiteOb!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!BEGIN INTERFACE RecordTwoSiteOb!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RecordTwoSiteOb_r(fileid, Obs)
!
!Purpose: Record a matrix operator List on a file whose ID is fileid.	
!
IMPLICIT NONE
REAL(KIND=rKind) :: Obs(systemSize,systemSize)
INTEGER, INTENT(IN) :: fileid
INTEGER :: j,k,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=systemSize
WRITE(mString,'(I4)') dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'

		DO k=1,systemSize
			WRITE(UNIT=fileid, FMT=fmtName) (Obs(k,j), j=1,dd)
		END DO
		
END SUBROUTINE RecordTwoSiteOb_r

SUBROUTINE RecordTwoSiteOb_c(fileid, Obs)
!
!Purpose: Record a matrix operator List on a file whose ID is fileid.	
!
IMPLICIT NONE
COMPLEX(KIND=rKind) :: Obs(systemSize,systemSize)
INTEGER, INTENT(IN) :: fileid
INTEGER :: j,k,dd
CHARACTER(4) :: mString
CHARACTER(64) :: fmtName

dd=systemSize
WRITE(mString,'(I4)') 2*dd
fmtName='('//TRIM(ADJUSTL(mString))//'E30.15)'

		DO k=1,systemSize	
					WRITE(UNIT=fileid, FMT=fmtName) (REAL(Obs(k,j)), j=1,dd), (AIMAG(Obs(k,j)), j=1,dd)
		END DO
		
END SUBROUTINE RecordTwoSiteOb_c

!!!!!!!!!!!!!END INTERFACE RecordTwoSiteOb!!!!!!!!!!!!!!!!!!!!!

	
END MODULE io_module
