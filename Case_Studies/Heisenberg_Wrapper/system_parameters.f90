MODULE system_parameters
!
! Purpose: Module to declare variables for OpenSourceTEBD v2.0
!
! Record of Revisions
!	Date	Programmer	Description of change
!	====	==========	=====================
!   2/22/09   M. L. Wall	v1.0 release
!   8/10/09   M. L. Wall	Boundary conditions toggle added
!				Switches changed to LOGICALS
!   9/24/09   M. L. Wall	v2.0 release
!
IMPLICIT NONE

!Select the precision and range of declared real variables
INTEGER, PARAMETER :: precis=15
INTEGER, PARAMETER :: range=30
INTEGER, PARAMETER :: rKind=SELECTED_REAL_KIND(p=precis,r=range)

! *** GLOBAL INPUT PARAMETERS
INTEGER :: systemSize=4 !Number of lattice sites
INTEGER :: maxFilling=1 !Maximum number of particles allowed per lattice site
INTEGER :: trotterOrder=2 !Order of the trotter decomposition, can be 2 or 5 for OBC and 2 for PBC
INTEGER :: chiMin=2 !Entanglement cutoff for initial sweep of ITP
INTEGER :: chiMax=5 !Entanglement cutoff for final sweep of ITP/RTP
REAL(KIND=rKind) :: dtITP=0.0001_rKind !Time step for ITP
INTEGER :: stepsForJudge=100 !Number of time steps before ITP convergence is checked
INTEGER :: maxITPsteps=2000 !greatest allowed number of ITP steps
REAL(KIND=rKind) ::  convCriterion1=0.00001_rKind ! convCriterion1 (convCriterion2) is the convergence criterion for the first (second) iteration of ITP
REAL(KIND=rKind) :: convCriterion2=0.000001_rKind 
COMPLEX(KIND=rKind) :: dtRTP=0.1_rKind !Time step for RTP
INTEGER ::  totalStep=10000 !Total number of RTP steps
INTEGER :: stepsForStore=10 !Number of RTP steps between each output write

! *** GLOBAL DERIVED PARAMETERS
INTEGER :: localSize !On-site hilbert space dimension

! *** I/O Data
CHARACTER(32) :: itpDir='ITPDATA/' !Directory where ITP data is stored
CHARACTER(32) :: itpExt='.dat' !Extension of ITP data
CHARACTER(32) :: rtpDir='RTPDATA/' !Directory where ITP data is stored
CHARACTER(32) :: rtpExt='.dat' !Extension of ITP data
CHARACTER(32) :: pdDir='PD_DATA/' !Directory where phase diagram data is stored

! *** Switches
LOGICAL :: print_switch=.TRUE. !Toggle printing to screen
LOGICAL :: ncSwitch=.FALSE. !Toggle number conservation
LOGICAL :: ITPreadMPDSwitch=.FALSE. !Toggle read in of Initial state (Phase Diagram code)
LOGICAL :: ITPwriteMPDswitch=.FALSE. !Toggle output of final state (Phase diagram Code)
CHARACTER :: BoundaryCond='O' !Boundary conditions: 'O' is open, 'P' is periodic
CHARACTER :: ITPopenKind='B' !'B' means output the MPS in binary, 'S' means output the MPS in scientific notation

! ** Status checks
INTEGER :: statInt !Status integer to ensure allocation/deallocation happens properly
INTEGER :: fileStatus !Status integer to ensure file open happens properly

! *** Symmetry conservation parameters
INTEGER :: totNum=4 !Total number of particles

END MODULE
