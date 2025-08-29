      MODULE MODEL_ATOM_COUNTS

         IMPLICIT NONE

         INTEGER, PARAMETER :: NUMB_MODEL_SPECIES =    328
         INTEGER, PARAMETER :: NUMB_MODEL_ATOMS   =   14
         TYPE MODEL_SPECIES_INFO
               CHARACTER( 16 ) :: NAME = "" 
               CHARACTER( 100) :: REPRESENTATIVE = "" 
               CHARACTER(  10) :: REPRESENTATION = "" 
               CHARACTER( 100) :: DSSTOX_ID = "" 
               CHARACTER( 100) :: SMILES = "" 
               CHARACTER(   2) :: SPECIES_TYPE = ""
               REAL            :: MOLWT = 0.0
         END TYPE 
         TYPE( MODEL_SPECIES_INFO ) ::  MODEL_SPECIES( NUMB_MODEL_SPECIES )
         REAL :: SPECIES_ATOMS_COUNTS( NUMB_MODEL_ATOMS,NUMB_MODEL_SPECIES )
         CHARACTER(2) :: MODEL_ATOMS(NUMB_MODEL_ATOMS) = (/ &
                          'CA', 'MN', 'CL', 'HG', 'BR', 'NA', 'SI', 'S ', 'TI', 'FE', 'K ', 'I ', 'N ', 'C ', '


         CONTAINS
            LOGICAL FUNCTION SET_ATOMS_COUNTS( MECHANISM_NAME )

            IMPLICIT NONE


            CHARACTER(LEN=*), INTENT( IN ) :: MECHANISM_NAME
! local:
            LOGICAL, SAVE :: FIRSTCALL = .TRUE.
            INTEGER       :: ISPECIES 

            IF( FIRSTCALL )THEN
               FIRSTCALL =  .FALSE.
               IF( TRIM( MECHANISM_NAME ) .NE. "CRACMM3M" ) SET_ATOMS_COUNTS = .FALSE.
               RETURN
            END IF
!               CA      MN      CL      HG      BR      NA      SI      S       TI      FE      K       I       N       C     
            SPECIES_ATOMS_COUNTS( :,    1) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! O3 GC namelist
            SPECIES_ATOMS_COUNTS( :,    2) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! O3P GC namelist
            SPECIES_ATOMS_COUNTS( :,    3) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! O1D GC namelist
            SPECIES_ATOMS_COUNTS( :,    4) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! H2O2 GC namelist
            SPECIES_ATOMS_COUNTS( :,    5) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! HO GC namelist
            SPECIES_ATOMS_COUNTS( :,    6) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! NO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,    7) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! NO GC namelist
            SPECIES_ATOMS_COUNTS( :,    8) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! NO3 GC namelist
            SPECIES_ATOMS_COUNTS( :,    9) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! HONO GC namelist
            SPECIES_ATOMS_COUNTS( :,   10) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! HNO3 GC namelist
            SPECIES_ATOMS_COUNTS( :,   11) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! HNO4 GC namelist
            SPECIES_ATOMS_COUNTS( :,   12) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! HO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   13) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! HCHO GC namelist
            SPECIES_ATOMS_COUNTS( :,   14) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! CO GC namelist
            SPECIES_ATOMS_COUNTS( :,   15) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ACD GC namelist
            SPECIES_ATOMS_COUNTS( :,   16) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! MO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   17) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! ALD GC namelist
            SPECIES_ATOMS_COUNTS( :,   18) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ETHP GC namelist
            SPECIES_ATOMS_COUNTS( :,   19) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! ACT GC namelist
            SPECIES_ATOMS_COUNTS( :,   20) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ACO3 GC namelist
            SPECIES_ATOMS_COUNTS( :,   21) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! UALD GC namelist
            SPECIES_ATOMS_COUNTS( :,   22) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! KET GC namelist
            SPECIES_ATOMS_COUNTS( :,   23) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MEK GC namelist
            SPECIES_ATOMS_COUNTS( :,   24) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! HKET GC namelist
            SPECIES_ATOMS_COUNTS( :,   25) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MACR GC namelist
            SPECIES_ATOMS_COUNTS( :,   26) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MACP GC namelist
            SPECIES_ATOMS_COUNTS( :,   27) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! XO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   28) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MVK GC namelist
            SPECIES_ATOMS_COUNTS( :,   29) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! GLY GC namelist
            SPECIES_ATOMS_COUNTS( :,   30) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! MGLY GC namelist
            SPECIES_ATOMS_COUNTS( :,   31) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! DCB1 GC namelist
            SPECIES_ATOMS_COUNTS( :,   32) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! DCB2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   33) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! BALD GC namelist
            SPECIES_ATOMS_COUNTS( :,   34) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! CHO GC namelist
            SPECIES_ATOMS_COUNTS( :,   35) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! OP1 GC namelist
            SPECIES_ATOMS_COUNTS( :,   36) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! OP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   37) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! OPB GC namelist
            SPECIES_ATOMS_COUNTS( :,   38) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! VOP3 GC namelist
            SPECIES_ATOMS_COUNTS( :,   39) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! PAA GC namelist
            SPECIES_ATOMS_COUNTS( :,   40) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    4.0/) ! ONIT GC namelist
            SPECIES_ATOMS_COUNTS( :,   41) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    2.0/) ! PAN GC namelist
            SPECIES_ATOMS_COUNTS( :,   42) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0,    0.0/) ! N2O5 GC namelist
            SPECIES_ATOMS_COUNTS( :,   43) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! SO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   44) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! SULF GC namelist
            SPECIES_ATOMS_COUNTS( :,   45) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! SULRXN GC namelist
            SPECIES_ATOMS_COUNTS( :,   46) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ETH GC namelist
            SPECIES_ATOMS_COUNTS( :,   47) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! HC3 GC namelist
            SPECIES_ATOMS_COUNTS( :,   48) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! HC3P GC namelist
            SPECIES_ATOMS_COUNTS( :,   49) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! HC5 GC namelist
            SPECIES_ATOMS_COUNTS( :,   50) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! HC5P GC namelist
            SPECIES_ATOMS_COUNTS( :,   51) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! HC10 GC namelist
            SPECIES_ATOMS_COUNTS( :,   52) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! HC10P GC namelist
            SPECIES_ATOMS_COUNTS( :,   53) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! HC10P2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   54) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ETE GC namelist
            SPECIES_ATOMS_COUNTS( :,   55) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ETEP GC namelist
            SPECIES_ATOMS_COUNTS( :,   56) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! OLT GC namelist
            SPECIES_ATOMS_COUNTS( :,   57) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! OLTP GC namelist
            SPECIES_ATOMS_COUNTS( :,   58) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! OLI GC namelist
            SPECIES_ATOMS_COUNTS( :,   59) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! OLIP GC namelist
            SPECIES_ATOMS_COUNTS( :,   60) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ACE GC namelist
            SPECIES_ATOMS_COUNTS( :,   61) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! ORA1 GC namelist
            SPECIES_ATOMS_COUNTS( :,   62) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! BEN GC namelist
            SPECIES_ATOMS_COUNTS( :,   63) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! BENP GC namelist
            SPECIES_ATOMS_COUNTS( :,   64) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! PHEN GC namelist
            SPECIES_ATOMS_COUNTS( :,   65) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! TOL GC namelist
            SPECIES_ATOMS_COUNTS( :,   66) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! CSL GC namelist
            SPECIES_ATOMS_COUNTS( :,   67) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! XYL GC namelist
            SPECIES_ATOMS_COUNTS( :,   68) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! EBZ GC namelist
            SPECIES_ATOMS_COUNTS( :,   69) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! STY GC namelist
            SPECIES_ATOMS_COUNTS( :,   70) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! TOLP GC namelist
            SPECIES_ATOMS_COUNTS( :,   71) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! XYLP GC namelist
            SPECIES_ATOMS_COUNTS( :,   72) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! EBZP GC namelist
            SPECIES_ATOMS_COUNTS( :,   73) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! STYP GC namelist
            SPECIES_ATOMS_COUNTS( :,   74) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! ISO GC namelist
            SPECIES_ATOMS_COUNTS( :,   75) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! ISOP GC namelist
            SPECIES_ATOMS_COUNTS( :,   76) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! API GC namelist
            SPECIES_ATOMS_COUNTS( :,   77) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! APIP1 GC namelist
            SPECIES_ATOMS_COUNTS( :,   78) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! APIP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   79) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,   10.0/) ! APINP1 GC namelist
            SPECIES_ATOMS_COUNTS( :,   80) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,   10.0/) ! APINP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   81) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! PINAL GC namelist
            SPECIES_ATOMS_COUNTS( :,   82) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! PINALP GC namelist
            SPECIES_ATOMS_COUNTS( :,   83) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! LIM GC namelist
            SPECIES_ATOMS_COUNTS( :,   84) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! LIMP1 GC namelist
            SPECIES_ATOMS_COUNTS( :,   85) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! LIMP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   86) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,   10.0/) ! LIMNP1 GC namelist
            SPECIES_ATOMS_COUNTS( :,   87) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,   10.0/) ! LIMNP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,   88) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! LIMAL GC namelist
            SPECIES_ATOMS_COUNTS( :,   89) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! LIMALP GC namelist
            SPECIES_ATOMS_COUNTS( :,   90) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! VHOM GC namelist
            SPECIES_ATOMS_COUNTS( :,   91) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   20.0/) ! VELHOM GC namelist
            SPECIES_ATOMS_COUNTS( :,   92) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! RCO3 GC namelist
            SPECIES_ATOMS_COUNTS( :,   93) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! ACTP GC namelist
            SPECIES_ATOMS_COUNTS( :,   94) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MEKP GC namelist
            SPECIES_ATOMS_COUNTS( :,   95) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! KETP GC namelist
            SPECIES_ATOMS_COUNTS( :,   96) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MCP GC namelist
            SPECIES_ATOMS_COUNTS( :,   97) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MVKP GC namelist
            SPECIES_ATOMS_COUNTS( :,   98) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! UALP GC namelist
            SPECIES_ATOMS_COUNTS( :,   99) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! DCB3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  100) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! BALP GC namelist
            SPECIES_ATOMS_COUNTS( :,  101) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! ADDC GC namelist
            SPECIES_ATOMS_COUNTS( :,  102) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! MCT GC namelist
            SPECIES_ATOMS_COUNTS( :,  103) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! MCTO GC namelist
            SPECIES_ATOMS_COUNTS( :,  104) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! MOH GC namelist
            SPECIES_ATOMS_COUNTS( :,  105) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! EOH GC namelist
            SPECIES_ATOMS_COUNTS( :,  106) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! ROH GC namelist
            SPECIES_ATOMS_COUNTS( :,  107) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ETEG GC namelist
            SPECIES_ATOMS_COUNTS( :,  108) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! ISHP GC namelist
            SPECIES_ATOMS_COUNTS( :,  109) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! IEPOX GC namelist
            SPECIES_ATOMS_COUNTS( :,  110) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! IEPOXP GC namelist
            SPECIES_ATOMS_COUNTS( :,  111) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! MAHP GC namelist
            SPECIES_ATOMS_COUNTS( :,  112) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ORA2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  113) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0/) ! ORAP GC namelist
            SPECIES_ATOMS_COUNTS( :,  114) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    3.0/) ! PPN GC namelist
            SPECIES_ATOMS_COUNTS( :,  115) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    4.0/) ! MPAN GC namelist
            SPECIES_ATOMS_COUNTS( :,  116) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    5.0/) ! INALD GC namelist
            SPECIES_ATOMS_COUNTS( :,  117) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    5.0/) ! ISONP GC namelist
            SPECIES_ATOMS_COUNTS( :,  118) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    5.0/) ! ISON GC namelist
            SPECIES_ATOMS_COUNTS( :,  119) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! IPX GC namelist
            SPECIES_ATOMS_COUNTS( :,  120) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,   10.0/) ! VTRPN GC namelist
            SPECIES_ATOMS_COUNTS( :,  121) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,   10.0/) ! VHONIT GC namelist
            SPECIES_ATOMS_COUNTS( :,  122) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! MCTP GC namelist
            SPECIES_ATOMS_COUNTS( :,  123) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    3.0/) ! OLNN GC namelist
            SPECIES_ATOMS_COUNTS( :,  124) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    3.0/) ! OLND GC namelist
            SPECIES_ATOMS_COUNTS( :,  125) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    6.0/) ! ADCN GC namelist
            SPECIES_ATOMS_COUNTS( :,  126) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! BAL1 GC namelist
            SPECIES_ATOMS_COUNTS( :,  127) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! BAL2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  128) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! ACRO GC namelist
            SPECIES_ATOMS_COUNTS( :,  129) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! BDE13 GC namelist
            SPECIES_ATOMS_COUNTS( :,  130) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! BDE13P GC namelist
            SPECIES_ATOMS_COUNTS( :,  131) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    3.0/) ! PROG GC namelist
            SPECIES_ATOMS_COUNTS( :,  132) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! FURAN GC namelist
            SPECIES_ATOMS_COUNTS( :,  133) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! FURANO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  134) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    4.0/) ! FURANONE GC namelist
            SPECIES_ATOMS_COUNTS( :,  135) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! VROCIOXY GC namelist
            SPECIES_ATOMS_COUNTS( :,  136) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    1.0/) ! SLOWROC GC namelist
            SPECIES_ATOMS_COUNTS( :,  137) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   15.0/) ! SESQ GC namelist
            SPECIES_ATOMS_COUNTS( :,  138) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   15.0/) ! SESQRO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  139) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,   15.0/) ! SESQNRO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  140) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! NAPH GC namelist
            SPECIES_ATOMS_COUNTS( :,  141) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! NAPHP GC namelist
            SPECIES_ATOMS_COUNTS( :,  142) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   14.0/) ! VROCP5ARO GC namelist
            SPECIES_ATOMS_COUNTS( :,  143) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   13.0/) ! VROCP6ARO GC namelist
            SPECIES_ATOMS_COUNTS( :,  144) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   14.0/) ! VROCP5AROP GC namelist
            SPECIES_ATOMS_COUNTS( :,  145) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   13.0/) ! VROCP6AROP GC namelist
            SPECIES_ATOMS_COUNTS( :,  146) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   30.0/) ! VROCN2ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  147) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   29.0/) ! VROCN1ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  148) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   28.0/) ! VROCP0ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  149) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   27.0/) ! VROCP1ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  150) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   24.0/) ! VROCP2ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  151) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   21.0/) ! VROCP3ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  152) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   18.0/) ! VROCP4ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  153) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   14.0/) ! VROCP5ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  154) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   13.0/) ! VROCP6ALK GC namelist
            SPECIES_ATOMS_COUNTS( :,  155) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   27.0/) ! VROCP1ALKP GC namelist
            SPECIES_ATOMS_COUNTS( :,  156) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   24.0/) ! VROCP2ALKP GC namelist
            SPECIES_ATOMS_COUNTS( :,  157) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   21.0/) ! VROCP3ALKP GC namelist
            SPECIES_ATOMS_COUNTS( :,  158) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   18.0/) ! VROCP4ALKP GC namelist
            SPECIES_ATOMS_COUNTS( :,  159) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   14.0/) ! VROCP5ALKP GC namelist
            SPECIES_ATOMS_COUNTS( :,  160) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   13.0/) ! VROCP6ALKP GC namelist
            SPECIES_ATOMS_COUNTS( :,  161) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   27.0/) ! VROCP1ALKP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  162) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   24.0/) ! VROCP2ALKP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  163) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   21.0/) ! VROCP3ALKP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  164) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   18.0/) ! VROCP4ALKP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  165) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   14.0/) ! VROCP5ALKP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  166) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   13.0/) ! VROCP6ALKP2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  167) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   17.0/) ! VROCN2OXY2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  168) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   11.0/) ! VROCN2OXY4 GC namelist
            SPECIES_ATOMS_COUNTS( :,  169) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! VROCN2OXY8 GC namelist
            SPECIES_ATOMS_COUNTS( :,  170) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   20.0/) ! VROCN1OXY1 GC namelist
            SPECIES_ATOMS_COUNTS( :,  171) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   12.0/) ! VROCN1OXY3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  172) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    8.0/) ! VROCN1OXY6 GC namelist
            SPECIES_ATOMS_COUNTS( :,  173) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   14.0/) ! VROCP0OXY2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  174) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   10.0/) ! VROCP0OXY4 GC namelist
            SPECIES_ATOMS_COUNTS( :,  175) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   17.0/) ! VROCP1OXY1 GC namelist
            SPECIES_ATOMS_COUNTS( :,  176) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   11.0/) ! VROCP1OXY3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  177) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   12.0/) ! VROCP2OXY2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  178) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   11.0/) ! VROCP3OXY2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  179) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    9.0/) ! VROCP4OXY2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  180) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,   11.0/) ! VROCP5OXY1 GC namelist
            SPECIES_ATOMS_COUNTS( :,  181) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    9.0/) ! VROCP6OXY1 GC namelist
            SPECIES_ATOMS_COUNTS( :,  182) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! ECH4 GC namelist
            SPECIES_ATOMS_COUNTS( :,  183) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! CO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  184) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    9.0/) ! VMTN1 GC namelist
            SPECIES_ATOMS_COUNTS( :,  185) = & 
            (/    0.0,    0.0,    2.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! CL2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  186) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! CL GC namelist
            SPECIES_ATOMS_COUNTS( :,  187) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! CLO GC namelist
            SPECIES_ATOMS_COUNTS( :,  188) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! OCLO GC namelist
            SPECIES_ATOMS_COUNTS( :,  189) = & 
            (/    0.0,    0.0,    2.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! CL2O2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  190) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! HOCL GC namelist
            SPECIES_ATOMS_COUNTS( :,  191) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! CLNO GC namelist
            SPECIES_ATOMS_COUNTS( :,  192) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! CLNO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  193) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! CLNO3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  194) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! HCOCL GC namelist
            SPECIES_ATOMS_COUNTS( :,  195) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! HCL GC namelist
            SPECIES_ATOMS_COUNTS( :,  196) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! CLOO GC namelist
            SPECIES_ATOMS_COUNTS( :,  197) = & 
            (/    0.0,    0.0,    0.0,    0.0,    2.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! BR2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  198) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! BR GC namelist
            SPECIES_ATOMS_COUNTS( :,  199) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! BRO GC namelist
            SPECIES_ATOMS_COUNTS( :,  200) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! OBRO GC namelist
            SPECIES_ATOMS_COUNTS( :,  201) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! HOBR GC namelist
            SPECIES_ATOMS_COUNTS( :,  202) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! BRNO GC namelist
            SPECIES_ATOMS_COUNTS( :,  203) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! BRNO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  204) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! BRNO3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  205) = & 
            (/    0.0,    0.0,    0.0,    0.0,    2.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! CH2BR2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  206) = & 
            (/    0.0,    0.0,    0.0,    0.0,    3.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! CHBR3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  207) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! HBR GC namelist
            SPECIES_ATOMS_COUNTS( :,  208) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! HCOBR GC namelist
            SPECIES_ATOMS_COUNTS( :,  209) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0,    0.0,    0.0/) ! I2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  210) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0/) ! I GC namelist
            SPECIES_ATOMS_COUNTS( :,  211) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0/) ! IO GC namelist
            SPECIES_ATOMS_COUNTS( :,  212) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0/) ! OIO GC namelist
            SPECIES_ATOMS_COUNTS( :,  213) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0,    0.0,    0.0/) ! I2O2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  214) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0/) ! HOI GC namelist
            SPECIES_ATOMS_COUNTS( :,  215) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0/) ! HI GC namelist
            SPECIES_ATOMS_COUNTS( :,  216) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    1.0,    0.0/) ! INO GC namelist
            SPECIES_ATOMS_COUNTS( :,  217) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    1.0,    0.0/) ! INO2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  218) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    1.0,    0.0/) ! INO3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  219) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    1.0/) ! CH3I GC namelist
            SPECIES_ATOMS_COUNTS( :,  220) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0,    0.0,    1.0/) ! CH2I2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  221) = & 
            (/    0.0,    0.0,    1.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! BRCL GC namelist
            SPECIES_ATOMS_COUNTS( :,  222) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0/) ! ICL GC namelist
            SPECIES_ATOMS_COUNTS( :,  223) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0/) ! IBR GC namelist
            SPECIES_ATOMS_COUNTS( :,  224) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    1.0/) ! CH2IBR GC namelist
            SPECIES_ATOMS_COUNTS( :,  225) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    1.0/) ! CH2ICL GC namelist
            SPECIES_ATOMS_COUNTS( :,  226) = & 
            (/    0.0,    0.0,    1.0,    0.0,    2.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! CHBR2CL GC namelist
            SPECIES_ATOMS_COUNTS( :,  227) = & 
            (/    0.0,    0.0,    2.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! CHBRCL2 GC namelist
            SPECIES_ATOMS_COUNTS( :,  228) = & 
            (/    0.0,    0.0,    1.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! CH2BRCL GC namelist
            SPECIES_ATOMS_COUNTS( :,  229) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! NO2PIJ GC namelist
            SPECIES_ATOMS_COUNTS( :,  230) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! NO2PK GC namelist
            SPECIES_ATOMS_COUNTS( :,  231) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0,    0.0,    0.0/) ! I2O3 GC namelist
            SPECIES_ATOMS_COUNTS( :,  232) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    2.0,    0.0,    0.0/) ! I2O4 GC namelist
            SPECIES_ATOMS_COUNTS( :,  233) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ASO4I NA namelist
            SPECIES_ATOMS_COUNTS( :,  234) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ASO4J NA namelist
            SPECIES_ATOMS_COUNTS( :,  235) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ASO4K NA namelist
            SPECIES_ATOMS_COUNTS( :,  236) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! ANH4I NA namelist
            SPECIES_ATOMS_COUNTS( :,  237) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! ANH4J NA namelist
            SPECIES_ATOMS_COUNTS( :,  238) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! ANH4K NA namelist
            SPECIES_ATOMS_COUNTS( :,  239) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! ANO3I NA namelist
            SPECIES_ATOMS_COUNTS( :,  240) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! ANO3J NA namelist
            SPECIES_ATOMS_COUNTS( :,  241) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! ANO3K NA namelist
            SPECIES_ATOMS_COUNTS( :,  242) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    7.0/) ! ASOATJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  243) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! AGLYOLIGJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  244) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AHOMJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  245) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AELHOMJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  246) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    6.0/) ! AORGCJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  247) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! AECI NA namelist
            SPECIES_ATOMS_COUNTS( :,  248) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! AECJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  249) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AOTHRI NA namelist
            SPECIES_ATOMS_COUNTS( :,  250) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AOTHRJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  251) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0/) ! AFEJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  252) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AALJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  253) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ASIJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  254) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ATIJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  255) = & 
            (/    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ACAJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  256) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AMGJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  257) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0/) ! AKJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  258) = & 
            (/    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AMNJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  259) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ACORS NA namelist
            SPECIES_ATOMS_COUNTS( :,  260) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ASOIL NA namelist
            SPECIES_ATOMS_COUNTS( :,  261) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! NUMATKN NA namelist
            SPECIES_ATOMS_COUNTS( :,  262) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! NUMACC NA namelist
            SPECIES_ATOMS_COUNTS( :,  263) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! NUMCOR NA namelist
            SPECIES_ATOMS_COUNTS( :,  264) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! SRFATKN NA namelist
            SPECIES_ATOMS_COUNTS( :,  265) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! SRFACC NA namelist
            SPECIES_ATOMS_COUNTS( :,  266) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! SRFCOR NA namelist
            SPECIES_ATOMS_COUNTS( :,  267) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AORGH2OJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  268) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AH2OI NA namelist
            SPECIES_ATOMS_COUNTS( :,  269) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AH2OJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  270) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AH2OK NA namelist
            SPECIES_ATOMS_COUNTS( :,  271) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AH3OPI NA namelist
            SPECIES_ATOMS_COUNTS( :,  272) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AH3OPJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  273) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AH3OPK NA namelist
            SPECIES_ATOMS_COUNTS( :,  274) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ANAI NA namelist
            SPECIES_ATOMS_COUNTS( :,  275) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ANAJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  276) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ACLI NA namelist
            SPECIES_ATOMS_COUNTS( :,  277) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ACLJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  278) = & 
            (/    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ACLK NA namelist
            SPECIES_ATOMS_COUNTS( :,  279) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ASEACAT NA namelist
            SPECIES_ATOMS_COUNTS( :,  280) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! APOCI NA namelist
            SPECIES_ATOMS_COUNTS( :,  281) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! APOCJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  282) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! APNCOMI NA namelist
            SPECIES_ATOMS_COUNTS( :,  283) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! APNCOMJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  284) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AOP3J NA namelist
            SPECIES_ATOMS_COUNTS( :,  285) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2ALKI NA namelist
            SPECIES_ATOMS_COUNTS( :,  286) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2ALKJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  287) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1ALKI NA namelist
            SPECIES_ATOMS_COUNTS( :,  288) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1ALKJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  289) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP0ALKI NA namelist
            SPECIES_ATOMS_COUNTS( :,  290) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP0ALKJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  291) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP1ALKI NA namelist
            SPECIES_ATOMS_COUNTS( :,  292) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP1ALKJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  293) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP2ALKJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  294) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP3ALKJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  295) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2OXY2I NA namelist
            SPECIES_ATOMS_COUNTS( :,  296) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2OXY2J NA namelist
            SPECIES_ATOMS_COUNTS( :,  297) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2OXY4I NA namelist
            SPECIES_ATOMS_COUNTS( :,  298) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2OXY4J NA namelist
            SPECIES_ATOMS_COUNTS( :,  299) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2OXY8I NA namelist
            SPECIES_ATOMS_COUNTS( :,  300) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN2OXY8J NA namelist
            SPECIES_ATOMS_COUNTS( :,  301) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1OXY1I NA namelist
            SPECIES_ATOMS_COUNTS( :,  302) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1OXY1J NA namelist
            SPECIES_ATOMS_COUNTS( :,  303) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1OXY3I NA namelist
            SPECIES_ATOMS_COUNTS( :,  304) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1OXY3J NA namelist
            SPECIES_ATOMS_COUNTS( :,  305) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1OXY6I NA namelist
            SPECIES_ATOMS_COUNTS( :,  306) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCN1OXY6J NA namelist
            SPECIES_ATOMS_COUNTS( :,  307) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP0OXY2I NA namelist
            SPECIES_ATOMS_COUNTS( :,  308) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP0OXY2J NA namelist
            SPECIES_ATOMS_COUNTS( :,  309) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP0OXY4I NA namelist
            SPECIES_ATOMS_COUNTS( :,  310) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP0OXY4J NA namelist
            SPECIES_ATOMS_COUNTS( :,  311) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP1OXY1I NA namelist
            SPECIES_ATOMS_COUNTS( :,  312) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP1OXY1J NA namelist
            SPECIES_ATOMS_COUNTS( :,  313) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP1OXY3I NA namelist
            SPECIES_ATOMS_COUNTS( :,  314) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP1OXY3J NA namelist
            SPECIES_ATOMS_COUNTS( :,  315) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP2OXY2J NA namelist
            SPECIES_ATOMS_COUNTS( :,  316) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AROCP3OXY2J NA namelist
            SPECIES_ATOMS_COUNTS( :,  317) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! AISO3NOSJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  318) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! AISO3OSJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  319) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! AISO4J NA namelist
            SPECIES_ATOMS_COUNTS( :,  320) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    5.0/) ! AISO5J NA namelist
            SPECIES_ATOMS_COUNTS( :,  321) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ATRPNJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  322) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AHONITJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  323) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ASEASTJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  324) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! AMTN1J NA namelist
            SPECIES_ATOMS_COUNTS( :,  325) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ABRJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  326) = & 
            (/    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0/) ! ABRK NA namelist
            SPECIES_ATOMS_COUNTS( :,  327) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0/) ! AHMSJ NA namelist
            SPECIES_ATOMS_COUNTS( :,  328) = & 
            (/    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    0.0/) ! NH3 NR namelist

            MODEL_SPECIES(     1 ) = MODEL_SPECIES_INFO('O3', & 
                                  'Ozone', & 
                                  'EXPLICIT', & 
                                  'DTXSID0021098', & 
                                  '[O-][O+]=O', & 
                                  'GC',  &
                                    48.000)

            MODEL_SPECIES(     2 ) = MODEL_SPECIES_INFO('O3P', & 
                                  'Ground state oxygen', & 
                                  'EXPLICIT', & 
                                  'DTXSID00170378', & 
                                  '[O]', & 
                                  'GC',  &
                                    16.000)

            MODEL_SPECIES(     3 ) = MODEL_SPECIES_INFO('O1D', & 
                                  'Excited oxygen', & 
                                  'EXPLICIT', & 
                                  'DTXSID00170378', & 
                                  '[O]', & 
                                  'GC',  &
                                    16.000)

            MODEL_SPECIES(     4 ) = MODEL_SPECIES_INFO('H2O2', & 
                                  'Hydrogen peroxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID2020715', & 
                                  'OO', & 
                                  'GC',  &
                                    34.000)

            MODEL_SPECIES(     5 ) = MODEL_SPECIES_INFO('HO', & 
                                  'Hydroxyl radical', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  '[OH]', & 
                                  'GC',  &
                                    17.000)

            MODEL_SPECIES(     6 ) = MODEL_SPECIES_INFO('NO2', & 
                                  'Nitrogen dioxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID7020974', & 
                                  'N(=O)[O]', & 
                                  'GC',  &
                                    46.000)

            MODEL_SPECIES(     7 ) = MODEL_SPECIES_INFO('NO', & 
                                  'Nitric oxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID1020938', & 
                                  '[N]=O', & 
                                  'GC',  &
                                    30.000)

            MODEL_SPECIES(     8 ) = MODEL_SPECIES_INFO('NO3', & 
                                  'Nitrate radical', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  '[O]N(=O)=O', & 
                                  'GC',  &
                                    62.000)

            MODEL_SPECIES(     9 ) = MODEL_SPECIES_INFO('HONO', & 
                                  'Nitrous acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID7064813', & 
                                  'N(=O)O', & 
                                  'GC',  &
                                    47.000)

            MODEL_SPECIES(    10 ) = MODEL_SPECIES_INFO('HNO3', & 
                                  'Nitric acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID5029685', & 
                                  '[N+](=O)(O)[O-]', & 
                                  'GC',  &
                                    63.000)

            MODEL_SPECIES(    11 ) = MODEL_SPECIES_INFO('HNO4', & 
                                  'Hydroxy nitrate', & 
                                  'EXPLICIT', & 
                                  'DTXSID201030501', & 
                                  '[N+](=O)([O-])OO', & 
                                  'GC',  &
                                    79.000)

            MODEL_SPECIES(    12 ) = MODEL_SPECIES_INFO('HO2', & 
                                  'Hydroperoxy', & 
                                  'EXPLICIT', & 
                                  'DTXSID30894777', & 
                                  'O[O-]', & 
                                  'GC',  &
                                    33.000)

            MODEL_SPECIES(    13 ) = MODEL_SPECIES_INFO('HCHO', & 
                                  'Formaldehyde', & 
                                  'EXPLICIT', & 
                                  'DTXSID7020637', & 
                                  'C=O', & 
                                  'GC',  &
                                    30.000)

            MODEL_SPECIES(    14 ) = MODEL_SPECIES_INFO('CO', & 
                                  'Carbon monoxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID5027273', & 
                                  '[C-]#[O+]', & 
                                  'GC',  &
                                    28.000)

            MODEL_SPECIES(    15 ) = MODEL_SPECIES_INFO('ACD', & 
                                  'Acetaldehyde', & 
                                  'EXPLICIT', & 
                                  'DTXSID5039224', & 
                                  'CC=O', & 
                                  'GC',  &
                                    44.000)

            MODEL_SPECIES(    16 ) = MODEL_SPECIES_INFO('MO2', & 
                                  'Methylperoxy', & 
                                  'EXPLICIT', & 
                                  'DTXSID10944007', & 
                                  'CO[O]', & 
                                  'GC',  &
                                    47.000)

            MODEL_SPECIES(    17 ) = MODEL_SPECIES_INFO('ALD', & 
                                  'Propanal', & 
                                  'LUMPED', & 
                                  'DTXSID2021658', & 
                                  'CCC=O', & 
                                  'GC',  &
                                    58.000)

            MODEL_SPECIES(    18 ) = MODEL_SPECIES_INFO('ETHP', & 
                                  'Ethylperoxy', & 
                                  'LUMPED', & 
                                  'DTXSID90953652', & 
                                  'CCO[O]', & 
                                  'GC',  &
                                    61.000)

            MODEL_SPECIES(    19 ) = MODEL_SPECIES_INFO('ACT', & 
                                  'Acetone', & 
                                  'EXPLICIT', & 
                                  'DTXSID8021482', & 
                                  'CC(C)=O', & 
                                  'GC',  &
                                    58.000)

            MODEL_SPECIES(    20 ) = MODEL_SPECIES_INFO('ACO3', & 
                                  'Acetylperoxy', & 
                                  'EXPLICIT', & 
                                  'DTXSID40957943', & 
                                  'CC(=O)O[O]', & 
                                  'GC',  &
                                    75.000)

            MODEL_SPECIES(    21 ) = MODEL_SPECIES_INFO('UALD', & 
                                  '2-Methylbut-2-enal', & 
                                  'LUMPED', & 
                                  'DTXSID00859414', & 
                                  'CC=C(C)C=O', & 
                                  'GC',  &
                                    84.100)

            MODEL_SPECIES(    22 ) = MODEL_SPECIES_INFO('KET', & 
                                  'Diethylketone', & 
                                  'LUMPED', & 
                                  'DTXSID6021820', & 
                                  'CCC(=O)CC', & 
                                  'GC',  &
                                    86.000)

            MODEL_SPECIES(    23 ) = MODEL_SPECIES_INFO('MEK', & 
                                  'Methyl ethyl ketone', & 
                                  'EXPLICIT', & 
                                  'DTXSID3021516', & 
                                  'CCC(C)=O', & 
                                  'GC',  &
                                    72.100)

            MODEL_SPECIES(    24 ) = MODEL_SPECIES_INFO('HKET', & 
                                  'Hydroxy acetone', & 
                                  'LUMPED', & 
                                  'DTXSID8051590', & 
                                  'CC(=O)CO', & 
                                  'GC',  &
                                    74.000)

            MODEL_SPECIES(    25 ) = MODEL_SPECIES_INFO('MACR', & 
                                  'Methacrolein', & 
                                  'LUMPED', & 
                                  'DTXSID0052540', & 
                                  'CC(=C)C=O', & 
                                  'GC',  &
                                    70.000)

            MODEL_SPECIES(    26 ) = MODEL_SPECIES_INFO('MACP', & 
                                  '(2-methylprop-2-enoyl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(=C)C(=O)O[O]', & 
                                  'GC',  &
                                   101.000)

            MODEL_SPECIES(    27 ) = MODEL_SPECIES_INFO('XO2', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'GC',  &
                                     1.000)

            MODEL_SPECIES(    28 ) = MODEL_SPECIES_INFO('MVK', & 
                                  'Methyl vinyl ketone', & 
                                  'EXPLICIT', & 
                                  'DTXSID3025671', & 
                                  'CC(=O)C=C', & 
                                  'GC',  &
                                    70.100)

            MODEL_SPECIES(    29 ) = MODEL_SPECIES_INFO('GLY', & 
                                  'Glyoxal', & 
                                  'LUMPED', & 
                                  'DTXSID5025364', & 
                                  'O=CC=O', & 
                                  'GC',  &
                                    58.000)

            MODEL_SPECIES(    30 ) = MODEL_SPECIES_INFO('MGLY', & 
                                  'Methyl glyoxal', & 
                                  'LUMPED', & 
                                  'DTXSID0021628', & 
                                  'CC(=O)C=O', & 
                                  'GC',  &
                                    72.000)

            MODEL_SPECIES(    31 ) = MODEL_SPECIES_INFO('DCB1', & 
                                  '2-methylbut-2-enedial', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=CC=C(C)C=O', & 
                                  'GC',  &
                                    98.000)

            MODEL_SPECIES(    32 ) = MODEL_SPECIES_INFO('DCB2', & 
                                  '2-methyl-4-oxopent-2-enal', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=CC(=CC(=O)C)C', & 
                                  'GC',  &
                                   112.100)

            MODEL_SPECIES(    33 ) = MODEL_SPECIES_INFO('BALD', & 
                                  'Benzaldehyde', & 
                                  'LUMPED', & 
                                  'DTXSID8039241', & 
                                  'O=CC1=CC=CC=C1', & 
                                  'GC',  &
                                   106.000)

            MODEL_SPECIES(    34 ) = MODEL_SPECIES_INFO('CHO', & 
                                  '(2-methylphenyl)oxidanyl', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC1=CC=CC=C1[O]', & 
                                  'GC',  &
                                   107.000)

            MODEL_SPECIES(    35 ) = MODEL_SPECIES_INFO('OP1', & 
                                  'Methyl hydroperoxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID10184401', & 
                                  'COO', & 
                                  'GC',  &
                                    48.000)

            MODEL_SPECIES(    36 ) = MODEL_SPECIES_INFO('OP2', & 
                                  'Ethyl hydroperoxide', & 
                                  'LUMPED', & 
                                  'DTXSID70184402', & 
                                  'CCOO', & 
                                  'GC',  &
                                    62.000)

            MODEL_SPECIES(    37 ) = MODEL_SPECIES_INFO('OPB', & 
                                  '2-hydroperoxy-2;6;6-trimethylbicyclo[3.1.1]heptan-3-ol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OOC1(C)C(O)CC2CC1C2(C)C', & 
                                  'GC',  &
                                   186.200)

            MODEL_SPECIES(    38 ) = MODEL_SPECIES_INFO('VOP3', & 
                                  '5-hydroperoxy-6-hydroxyoctan-3-one', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCC(=O)CC(OO)C(O)CC', & 
                                  'GC',  &
                                   176.200)

            MODEL_SPECIES(    39 ) = MODEL_SPECIES_INFO('PAA', & 
                                  'Peroxyacetic acid', & 
                                  'LUMPED', & 
                                  'DTXSID1025853', & 
                                  'CC(=O)OO', & 
                                  'GC',  &
                                    76.000)

            MODEL_SPECIES(    40 ) = MODEL_SPECIES_INFO('ONIT', & 
                                  'Butan-2-yl nitrate', & 
                                  'LUMPED', & 
                                  'DTXSID00871813', & 
                                  'CCC(C)O[N+](=O)[O-]', & 
                                  'GC',  &
                                   119.000)

            MODEL_SPECIES(    41 ) = MODEL_SPECIES_INFO('PAN', & 
                                  'Peroxyacetyl nitrate', & 
                                  'LUMPED', & 
                                  'DTXSID4062301', & 
                                  'CC(=O)OON(=O)=O', & 
                                  'GC',  &
                                   121.000)

            MODEL_SPECIES(    42 ) = MODEL_SPECIES_INFO('N2O5', & 
                                  'Dinitrogen pentoxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID90143672', & 
                                  '[N+](=O)([O-])O[N+](=O)[O-]', & 
                                  'GC',  &
                                   108.000)

            MODEL_SPECIES(    43 ) = MODEL_SPECIES_INFO('SO2', & 
                                  'Sulfur dioxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID6029672', & 
                                  'O=S=O', & 
                                  'GC',  &
                                    64.000)

            MODEL_SPECIES(    44 ) = MODEL_SPECIES_INFO('SULF', & 
                                  'Sulfuric acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID5029683', & 
                                  'OS(=O)(=O)O', & 
                                  'GC',  &
                                    98.000)

            MODEL_SPECIES(    45 ) = MODEL_SPECIES_INFO('SULRXN', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'GC',  &
                                    98.000)

            MODEL_SPECIES(    46 ) = MODEL_SPECIES_INFO('ETH', & 
                                  'Ethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID6026377', & 
                                  'CC', & 
                                  'GC',  &
                                    30.100)

            MODEL_SPECIES(    47 ) = MODEL_SPECIES_INFO('HC3', & 
                                  'Propane', & 
                                  'LUMPED', & 
                                  'DTXSID5026386', & 
                                  'CCC', & 
                                  'GC',  &
                                    44.100)

            MODEL_SPECIES(    48 ) = MODEL_SPECIES_INFO('HC3P', & 
                                  'Isopropyl peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(C)O[O]', & 
                                  'GC',  &
                                    75.000)

            MODEL_SPECIES(    49 ) = MODEL_SPECIES_INFO('HC5', & 
                                  'Pentane', & 
                                  'LUMPED', & 
                                  'DTXSID2025846', & 
                                  'CCCCC', & 
                                  'GC',  &
                                    72.100)

            MODEL_SPECIES(    50 ) = MODEL_SPECIES_INFO('HC5P', & 
                                  'Pentan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCC(O[O])CC', & 
                                  'GC',  &
                                   103.000)

            MODEL_SPECIES(    51 ) = MODEL_SPECIES_INFO('HC10', & 
                                  'Decane', & 
                                  'LUMPED', & 
                                  'DTXSID6024913', & 
                                  'CCCCCCCCCC', & 
                                  'GC',  &
                                   142.280)

            MODEL_SPECIES(    52 ) = MODEL_SPECIES_INFO('HC10P', & 
                                  'Decan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCC(CC)O[O]', & 
                                  'GC',  &
                                   173.270)

            MODEL_SPECIES(    53 ) = MODEL_SPECIES_INFO('HC10P2', & 
                                  '(8-hydroxydecan-5-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCC(O[O])CCC(O)CC', & 
                                  'GC',  &
                                   189.270)

            MODEL_SPECIES(    54 ) = MODEL_SPECIES_INFO('ETE', & 
                                  'Ethylene', & 
                                  'EXPLICIT', & 
                                  'DTXSID1026378', & 
                                  'C=C', & 
                                  'GC',  &
                                    28.100)

            MODEL_SPECIES(    55 ) = MODEL_SPECIES_INFO('ETEP', & 
                                  '(2-hydroxyethyl)peroxy', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'OCCO[O]', & 
                                  'GC',  &
                                    77.000)

            MODEL_SPECIES(    56 ) = MODEL_SPECIES_INFO('OLT', & 
                                  '1-Propene', & 
                                  'LUMPED', & 
                                  'DTXSID5021205', & 
                                  'CC=C', & 
                                  'GC',  &
                                    42.000)

            MODEL_SPECIES(    57 ) = MODEL_SPECIES_INFO('OLTP', & 
                                  '(1-hydroxypropan-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(CO)O[O]', & 
                                  'GC',  &
                                    91.000)

            MODEL_SPECIES(    58 ) = MODEL_SPECIES_INFO('OLI', & 
                                  '2-Methyl-2-butene', & 
                                  'LUMPED', & 
                                  'DTXSID8027165', & 
                                  'CC=C(C)C', & 
                                  'GC',  &
                                    70.100)

            MODEL_SPECIES(    59 ) = MODEL_SPECIES_INFO('OLIP', & 
                                  '(3-hydroxy-2-methylbutan-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC(C)(C)C(C)O', & 
                                  'GC',  &
                                   119.000)

            MODEL_SPECIES(    60 ) = MODEL_SPECIES_INFO('ACE', & 
                                  'Acetylene', & 
                                  'EXPLICIT', & 
                                  'DTXSID6026379', & 
                                  'C#C', & 
                                  'GC',  &
                                    26.000)

            MODEL_SPECIES(    61 ) = MODEL_SPECIES_INFO('ORA1', & 
                                  'Formic acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID2024115', & 
                                  'OC=O', & 
                                  'GC',  &
                                    46.000)

            MODEL_SPECIES(    62 ) = MODEL_SPECIES_INFO('BEN', & 
                                  'Benzene', & 
                                  'EXPLICIT', & 
                                  'DTXSID3039242', & 
                                  'C1=CC=CC=C1', & 
                                  'GC',  &
                                    78.110)

            MODEL_SPECIES(    63 ) = MODEL_SPECIES_INFO('BENP', & 
                                  '{8-hydroxy-6;7-dioxabicyclo[3.2.1]oct-3-en-2-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1C=CC2OOC1C2O', & 
                                  'GC',  &
                                   159.120)

            MODEL_SPECIES(    64 ) = MODEL_SPECIES_INFO('PHEN', & 
                                  'Phenol', & 
                                  'EXPLICIT', & 
                                  'DTXSID5021124', & 
                                  'OC1=CC=CC=C1', & 
                                  'GC',  &
                                    94.110)

            MODEL_SPECIES(    65 ) = MODEL_SPECIES_INFO('TOL', & 
                                  'Toluene', & 
                                  'EXPLICIT', & 
                                  'DTXSID7021360', & 
                                  'CC1=CC=CC=C1', & 
                                  'GC',  &
                                    92.140)

            MODEL_SPECIES(    66 ) = MODEL_SPECIES_INFO('CSL', & 
                                  'o-cresol', & 
                                  'LUMPED', & 
                                  'DTXSID8021808', & 
                                  'CC1=C(O)C=CC=C1', & 
                                  'GC',  &
                                   108.140)

            MODEL_SPECIES(    67 ) = MODEL_SPECIES_INFO('XYL', & 
                                  'm-Xylene', & 
                                  'LUMPED', & 
                                  'DTXSID6026298', & 
                                  'CC1=CC(C)=CC=C1', & 
                                  'GC',  &
                                   106.200)

            MODEL_SPECIES(    68 ) = MODEL_SPECIES_INFO('EBZ', & 
                                  'Ethylbenzene', & 
                                  'LUMPED', & 
                                  'DTXSID3020596', & 
                                  'CCC1=CC=CC=C1', & 
                                  'GC',  &
                                   106.200)

            MODEL_SPECIES(    69 ) = MODEL_SPECIES_INFO('STY', & 
                                  'Styrene', & 
                                  'EXPLICIT', & 
                                  'DTXSID2021284', & 
                                  'C=CC1=CC=CC=C1', & 
                                  'GC',  &
                                   104.000)

            MODEL_SPECIES(    70 ) = MODEL_SPECIES_INFO('TOLP', & 
                                  '{8-hydroxy-5-methyl-6;7-dioxabicyclo[3.2.1]oct-3-en-2-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1C=CC2(C)OOC1C2O', & 
                                  'GC',  &
                                   173.140)

            MODEL_SPECIES(    71 ) = MODEL_SPECIES_INFO('XYLP', & 
                                  '{8-hydroxy-1;5-dimethyl-6;7-dioxabicyclo[3.2.1]oct-3-en-2-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1C=CC2(C)OOC1(C)C2O', & 
                                  'GC',  &
                                   187.170)

            MODEL_SPECIES(    72 ) = MODEL_SPECIES_INFO('EBZP', & 
                                  '{5-ethyl-8-hydroxy-6;7-dioxabicyclo[3.2.1]oct-3-en-2-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1C=CC2(CC)OOC1C2O', & 
                                  'GC',  &
                                   187.170)

            MODEL_SPECIES(    73 ) = MODEL_SPECIES_INFO('STYP', & 
                                  '(2-hydroxy-2-phenylethyl)peroxy', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  '[O]OCC(O)C1=CC=CC=C1', & 
                                  'GC',  &
                                   153.000)

            MODEL_SPECIES(    74 ) = MODEL_SPECIES_INFO('ISO', & 
                                  'Isoprene', & 
                                  'EXPLICIT', & 
                                  'DTXSID2020761', & 
                                  'CC(=C)C=C', & 
                                  'GC',  &
                                    68.100)

            MODEL_SPECIES(    75 ) = MODEL_SPECIES_INFO('ISOP', & 
                                  '(1-hydroxy-3-methylbut-3-en-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OCC(O[O])C(C)=C', & 
                                  'GC',  &
                                   117.000)

            MODEL_SPECIES(    76 ) = MODEL_SPECIES_INFO('API', & 
                                  'alpha-pinene', & 
                                  'LUMPED', & 
                                  'DTXSID4026501', & 
                                  'CC1=CCC2CC1C2(C)C', & 
                                  'GC',  &
                                   136.400)

            MODEL_SPECIES(    77 ) = MODEL_SPECIES_INFO('APIP1', & 
                                  '{3-hydroxy-2;6;6-trimethylbicyclo[3.1.1]heptan-2-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1(C)C(O)CC2CC1C2(C)C', & 
                                  'GC',  &
                                   185.000)

            MODEL_SPECIES(    78 ) = MODEL_SPECIES_INFO('APIP2', & 
                                  '[4-hydroperoxy-6-(2-hydroperoxypropan-2-yl)-2-hydroxy-3-methylcyclohex-3-en-1-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC1=C(CC(C(O[O])C1O)C(C)(C)OO)OO', & 
                                  'GC',  &
                                   249.000)

            MODEL_SPECIES(    79 ) = MODEL_SPECIES_INFO('APINP1', & 
                                  '[2;6;6-trimethyl-3-(nitrooxy)bicyclo[3.1.1]heptan-2-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1(C)C(ON(=O)=O)CC2CC1C2(C)C', & 
                                  'GC',  &
                                   230.000)

            MODEL_SPECIES(    80 ) = MODEL_SPECIES_INFO('APINP2', & 
                                  '[4-hydroperoxy-6-(2-hydroperoxypropan-2-yl)-3-methyl-2-(nitrooxy)cyclohex-3-en-1-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC1=C(CC(C(O[O])C1(ON(=O)=O))C(C)(C)OO)OO', & 
                                  'GC',  &
                                   294.000)

            MODEL_SPECIES(    81 ) = MODEL_SPECIES_INFO('PINAL', & 
                                  'Pinonaldehyde', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=CCC1CC(C(=O)C)C1(C)C', & 
                                  'GC',  &
                                   168.000)

            MODEL_SPECIES(    82 ) = MODEL_SPECIES_INFO('PINALP', & 
                                  '[3-acetyl-2;2-dimethyl-1-(2-oxoethyl)cyclobutyl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=CCC1(O[O])CC(C(=O)C)C1(C)C', & 
                                  'GC',  &
                                   199.000)

            MODEL_SPECIES(    83 ) = MODEL_SPECIES_INFO('LIM', & 
                                  'D-Limonene', & 
                                  'LUMPED', & 
                                  'DTXSID1020778', & 
                                  'CC(=C)[C@@H]1CCC(C)=CC1', & 
                                  'GC',  &
                                   136.300)

            MODEL_SPECIES(    84 ) = MODEL_SPECIES_INFO('LIMP1', & 
                                  '[2-hydroxy-1-methyl-4-(prop-1-en-2-yl)cyclohexyl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1(C)CCC(CC1O)C(=C)C', & 
                                  'GC',  &
                                   185.000)

            MODEL_SPECIES(    85 ) = MODEL_SPECIES_INFO('LIMP2', & 
                                  '[6-hydroperoxy-4-(2-hydroperoxy-1-hydroxypropan-2-yl)-1-methylcyclohex-2-en-1-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C(OO)1C(O[O])(C)C=CC(C(OO)(C)CO)C1', & 
                                  'GC',  &
                                   249.000)

            MODEL_SPECIES(    86 ) = MODEL_SPECIES_INFO('LIMNP1', & 
                                  '[1-methyl-2-(nitrooxy)-4-(prop-1-en-2-yl)cyclohexyl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O-][N+](=O)OC1CC(CCC1(C)O[O])C(=C)C', & 
                                  'GC',  &
                                   230.000)

            MODEL_SPECIES(    87 ) = MODEL_SPECIES_INFO('LIMNP2', & 
                                  '{6-hydroperoxy-4-[2-hydroperoxy-1-(nitrooxy)propan-2-yl]-1-methylcyclohex-2-en-1-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C(OO)1C(O[O])(C)C=CC(C(OO)(C)C(ON(=O)=O))C1', & 
                                  'GC',  &
                                   294.000)

            MODEL_SPECIES(    88 ) = MODEL_SPECIES_INFO('LIMAL', & 
                                  'Limonaldehyde', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=CCC(CCC(=O)C)C(=C)C', & 
                                  'GC',  &
                                   168.000)

            MODEL_SPECIES(    89 ) = MODEL_SPECIES_INFO('LIMALP', & 
                                  '[1-hydroxy-2-methyl-5-oxo-3-(3-oxobutyl)pentan-2-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=CCC(CCC(=O)C)C(C)(CO)O[O]', & 
                                  'GC',  &
                                   217.000)

            MODEL_SPECIES(    90 ) = MODEL_SPECIES_INFO('VHOM', & 
                                  '6;7-dihydroperoxy-3;3;6-trimethyl-tetrahydro-3aH-1;2-benzodioxol-5-ol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC1CC2C(OOC2(C)C)C(OO)C1(C)OO', & 
                                  'GC',  &
                                   250.000)

            MODEL_SPECIES(    91 ) = MODEL_SPECIES_INFO('VELHOM', & 
                                  '6-hydroperoxy-7-({3-hydroxy-2;6;6-trimethylbicyclo[3.1.1]heptan-2-yl}peroxy)-3;3;6-trimethyl-tetrahy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC1CC2C(OOC2(C)C)C(OOC3(C)C4C(C)(C)C(C4)CC3O)C1(C)OO', & 
                                  'GC',  &
                                   402.000)

            MODEL_SPECIES(    92 ) = MODEL_SPECIES_INFO('RCO3', & 
                                  'Propanoylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCC(=O)O[O]', & 
                                  'GC',  &
                                    89.000)

            MODEL_SPECIES(    93 ) = MODEL_SPECIES_INFO('ACTP', & 
                                  '(2-oxopropyl)peroxy', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'CC(=O)CO[O]', & 
                                  'GC',  &
                                    89.000)

            MODEL_SPECIES(    94 ) = MODEL_SPECIES_INFO('MEKP', & 
                                  '(3-oxobutyl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OCCC(=O)C', & 
                                  'GC',  &
                                   103.000)

            MODEL_SPECIES(    95 ) = MODEL_SPECIES_INFO('KETP', & 
                                  '(3-oxopentan-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCC(C(C)O[O])=O', & 
                                  'GC',  &
                                   117.000)

            MODEL_SPECIES(    96 ) = MODEL_SPECIES_INFO('MCP', & 
                                  '(1-hydroxy-2-methyl-3-oxopropan-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OCC(C)(O[O])C=O', & 
                                  'GC',  &
                                   119.000)

            MODEL_SPECIES(    97 ) = MODEL_SPECIES_INFO('MVKP', & 
                                  '(2-hydroxy-3-oxobutyl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(=O)C(O)CO[O]', & 
                                  'GC',  &
                                   119.000)

            MODEL_SPECIES(    98 ) = MODEL_SPECIES_INFO('UALP', & 
                                  '(3-hydroxy-3-methyl-4-oxobutan-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(O[O])C(C)(O)C=O', & 
                                  'GC',  &
                                   133.000)

            MODEL_SPECIES(    99 ) = MODEL_SPECIES_INFO('DCB3', & 
                                  'but-2-enedial', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=CC=CC=O', & 
                                  'GC',  &
                                    84.000)

            MODEL_SPECIES(   100 ) = MODEL_SPECIES_INFO('BALP', & 
                                  'benzoylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=C(O[O])C1=CC=CC=C1', & 
                                  'GC',  &
                                   137.000)

            MODEL_SPECIES(   101 ) = MODEL_SPECIES_INFO('ADDC', & 
                                  '3-methyl-5-oxocyclohex-3-en-1-yloxidanyl', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC1=CC(O)=CC([O])C1', & 
                                  'GC',  &
                                   125.000)

            MODEL_SPECIES(   102 ) = MODEL_SPECIES_INFO('MCT', & 
                                  'Catechol', & 
                                  'LUMPED', & 
                                  'DTXSID3020257', & 
                                  'OC1=C(O)C=CC=C1', & 
                                  'GC',  &
                                   110.110)

            MODEL_SPECIES(   103 ) = MODEL_SPECIES_INFO('MCTO', & 
                                  '(2-hydroxyphenyl)oxidanyl', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]C1=CC=CC=C1O', & 
                                  'GC',  &
                                   109.000)

            MODEL_SPECIES(   104 ) = MODEL_SPECIES_INFO('MOH', & 
                                  'Methanol', & 
                                  'EXPLICIT', & 
                                  'DTXSID2021731', & 
                                  'CO', & 
                                  'GC',  &
                                    32.000)

            MODEL_SPECIES(   105 ) = MODEL_SPECIES_INFO('EOH', & 
                                  'Ethanol', & 
                                  'EXPLICIT', & 
                                  'DTXSID9020584', & 
                                  'CCO', & 
                                  'GC',  &
                                    46.100)

            MODEL_SPECIES(   106 ) = MODEL_SPECIES_INFO('ROH', & 
                                  'Propanol', & 
                                  'LUMPED', & 
                                  'DTXSID2021739', & 
                                  'CCCO', & 
                                  'GC',  &
                                    60.000)

            MODEL_SPECIES(   107 ) = MODEL_SPECIES_INFO('ETEG', & 
                                  'Ethylene glycol', & 
                                  'EXPLICIT', & 
                                  'DTXSID8020597', & 
                                  'OCCO', & 
                                  'GC',  &
                                    62.100)

            MODEL_SPECIES(   108 ) = MODEL_SPECIES_INFO('ISHP', & 
                                  '2-hydroperoxy-2-methylbut-3-en-1-ol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C=CC(OO)(CO)C', & 
                                  'GC',  &
                                   118.000)

            MODEL_SPECIES(   109 ) = MODEL_SPECIES_INFO('IEPOX', & 
                                  'Isoprene epoxydiol', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'OCC1OC1(C)CO', & 
                                  'GC',  &
                                   118.100)

            MODEL_SPECIES(   110 ) = MODEL_SPECIES_INFO('IEPOXP', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'GC',  &
                                   118.100)

            MODEL_SPECIES(   111 ) = MODEL_SPECIES_INFO('MAHP', & 
                                  '2-methylprop-2-eneperoxoic acid', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C=C(C)C(OO)=O', & 
                                  'GC',  &
                                   102.000)

            MODEL_SPECIES(   112 ) = MODEL_SPECIES_INFO('ORA2', & 
                                  'Acetic acid', & 
                                  'LUMPED', & 
                                  'DTXSID5024394', & 
                                  'CC(O)=O', & 
                                  'GC',  &
                                    60.200)

            MODEL_SPECIES(   113 ) = MODEL_SPECIES_INFO('ORAP', & 
                                  '(carboxymethyl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OCC(=O)O', & 
                                  'GC',  &
                                    91.000)

            MODEL_SPECIES(   114 ) = MODEL_SPECIES_INFO('PPN', & 
                                  'Peroxypopionyl nitrate', & 
                                  'EXPLICIT', & 
                                  'DTXSID90206675', & 
                                  'CCC(=O)OO[N+](=O)[O-]', & 
                                  'GC',  &
                                   135.000)

            MODEL_SPECIES(   115 ) = MODEL_SPECIES_INFO('MPAN', & 
                                  'Peroxymethacryloyl nitrate', & 
                                  'LUMPED', & 
                                  'DTXSID10236878', & 
                                  'O=N(=O)OOC(=O)C(=C)C', & 
                                  'GC',  &
                                   147.100)

            MODEL_SPECIES(   116 ) = MODEL_SPECIES_INFO('INALD', & 
                                  '2-hydroperoxy-4-hydroxy-3-methyl-3-(nitrooxy)butanal', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(CO)(O[N+](=O)O)C(C=O)OO', & 
                                  'GC',  &
                                   195.100)

            MODEL_SPECIES(   117 ) = MODEL_SPECIES_INFO('ISONP', & 
                                  '[(2E)-3-methyl-4-(nitrooxy)but-2-en-1-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C/C(=C\CO[O])CON(=O)=O', & 
                                  'GC',  &
                                   162.100)

            MODEL_SPECIES(   118 ) = MODEL_SPECIES_INFO('ISON', & 
                                  '2-methyl-2-(nitrooxy)but-3-en-1-ol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OCC(C)(C=C)ON(=O)=O', & 
                                  'GC',  &
                                   147.000)

            MODEL_SPECIES(   119 ) = MODEL_SPECIES_INFO('IPX', & 
                                  '3-(1-hydroperoxy-2-hydroxypropan-2-yl)oxiran-2-ol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(O)(COO)C1OC1O', & 
                                  'GC',  &
                                   150.100)

            MODEL_SPECIES(   120 ) = MODEL_SPECIES_INFO('VTRPN', & 
                                  '3-hydroxy-2;6;6-trimethylbicyclo[3.1.1]heptan-2-yl nitrate', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'O=N(=O)OC1(C)C(O)CC2CC1C2(C)C', & 
                                  'GC',  &
                                   215.000)

            MODEL_SPECIES(   121 ) = MODEL_SPECIES_INFO('VHONIT', & 
                                  '2-(3-hydroperoxy-4;5-dihydroxy-4-methylcyclohexyl)propan-2-yl nitrate', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(C)(O[N+]([O-])=O)C1CC(O)C(C)(O)C(C1)OO', & 
                                  'GC',  &
                                   265.000)

            MODEL_SPECIES(   122 ) = MODEL_SPECIES_INFO('MCTP', & 
                                  '(2-hydroxyphenyl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1=CC=CC=C1O', & 
                                  'GC',  &
                                   125.000)

            MODEL_SPECIES(   123 ) = MODEL_SPECIES_INFO('OLNN', & 
                                  '[1-(nitrooxy)propan-2-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(O[O])CO[N+]([O-])=O', & 
                                  'GC',  &
                                   136.000)

            MODEL_SPECIES(   124 ) = MODEL_SPECIES_INFO('OLND', & 
                                  '[1-(nitrooxy)propan-2-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(O[O])CO[N+]([O-])=O', & 
                                  'GC',  &
                                   136.000)

            MODEL_SPECIES(   125 ) = MODEL_SPECIES_INFO('ADCN', & 
                                  '4-oxocyclohex-2-en-1-yl nitrate', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC1=C[C]C(O[N+]([O-])=O)C=C1', & 
                                  'GC',  &
                                   155.000)

            MODEL_SPECIES(   126 ) = MODEL_SPECIES_INFO('BAL1', & 
                                  'Phenylperoxy radical', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1=CC=CC=C1', & 
                                  'GC',  &
                                   109.100)

            MODEL_SPECIES(   127 ) = MODEL_SPECIES_INFO('BAL2', & 
                                  'phenoxy radical', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]C1=CC=CC=C1', & 
                                  'GC',  &
                                    93.100)

            MODEL_SPECIES(   128 ) = MODEL_SPECIES_INFO('ACRO', & 
                                  'Acrolein', & 
                                  'EXPLICIT', & 
                                  'DTXSID5020023', & 
                                  'C=CC=O', & 
                                  'GC',  &
                                    56.100)

            MODEL_SPECIES(   129 ) = MODEL_SPECIES_INFO('BDE13', & 
                                  '1;3-Butadiene', & 
                                  'EXPLICIT', & 
                                  'DTXSID3020203', & 
                                  'C=CC=C', & 
                                  'GC',  &
                                    54.100)

            MODEL_SPECIES(   130 ) = MODEL_SPECIES_INFO('BDE13P', & 
                                  '(1-hydroxybut-3-en-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C=CC(O[O])CO', & 
                                  'GC',  &
                                   103.000)

            MODEL_SPECIES(   131 ) = MODEL_SPECIES_INFO('PROG', & 
                                  '1;2-Propylene glycol', & 
                                  'EXPLICIT', & 
                                  'DTXSID0021206', & 
                                  'CC(O)CO', & 
                                  'GC',  &
                                    76.100)

            MODEL_SPECIES(   132 ) = MODEL_SPECIES_INFO('FURAN', & 
                                  'Furfural', & 
                                  'LUMPED', & 
                                  'DTXSID1020647', & 
                                  'O=CC1=CC=CO1', & 
                                  'GC',  &
                                    96.100)

            MODEL_SPECIES(   133 ) = MODEL_SPECIES_INFO('FURANO2', & 
                                  '(2-formyl-5-hydroxy-5H-furan-2-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC1C=CC(O1)(O[O])(C=O)', & 
                                  'GC',  &
                                   145.100)

            MODEL_SPECIES(   134 ) = MODEL_SPECIES_INFO('FURANONE', & 
                                  '2-hydroxy-2H-furan-5-one', & 
                                  'LUMPED', & 
                                  'DTXSID10930763', & 
                                  'C1=CC(=O)OC1O', & 
                                  'GC',  &
                                   100.100)

            MODEL_SPECIES(   135 ) = MODEL_SPECIES_INFO('VROCIOXY', & 
                                  'Decamethylcyclopentasiloxane', & 
                                  'LUMPED', & 
                                  'DTXSID1027184', & 
                                  'C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1', & 
                                  'GC',  &
                                   247.000)

            MODEL_SPECIES(   136 ) = MODEL_SPECIES_INFO('SLOWROC', & 
                                  'Hydrogen cyanide', & 
                                  'LUMPED', & 
                                  'DTXSID9024148', & 
                                  'C#N', & 
                                  'GC',  &
                                    75.400)

            MODEL_SPECIES(   137 ) = MODEL_SPECIES_INFO('SESQ', & 
                                  'b-caryophyllene', & 
                                  'LUMPED', & 
                                  'DTXSID8024739', & 
                                  'C/C1=C/CCC(=C)C2CC(C)(C)C2CC\1', & 
                                  'GC',  &
                                   204.400)

            MODEL_SPECIES(   138 ) = MODEL_SPECIES_INFO('SESQRO2', & 
                                  '{5-hydroxy-4;11;11-trimethyl-8-methylidenebicyclo[7.2.0]undecan-4-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1(C)CCC2C(CC2(C)C)C(=C)CCC1O', & 
                                  'GC',  &
                                   253.400)

            MODEL_SPECIES(   139 ) = MODEL_SPECIES_INFO('SESQNRO2', & 
                                  '[4;11;11-trimethyl-8-methylidene-5-(nitrooxy)bicyclo[7.2.0]undecan-4-yl]peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[O]OC1(C)CCC2C(CC2(C)C)C(=C)CCC1O[N+](=O)[O-]', & 
                                  'GC',  &
                                   298.400)

            MODEL_SPECIES(   140 ) = MODEL_SPECIES_INFO('NAPH', & 
                                  'Naphthalene', & 
                                  'LUMPED', & 
                                  'DTXSID8020913', & 
                                  'C1=CC2=CC=CC=C2C=C1', & 
                                  'GC',  &
                                   128.170)

            MODEL_SPECIES(   141 ) = MODEL_SPECIES_INFO('NAPHP', & 
                                  '{8-hydroxy-10;11-dioxatricyclo[7.2.1.0^{2;7}]dodeca-2;4;6-trien-12-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C12=CC=CC=C1C3OOC(C3O[O])C2(O)', & 
                                  'GC',  &
                                   209.170)

            MODEL_SPECIES(   142 ) = MODEL_SPECIES_INFO('VROCP5ARO', & 
                                  'Benzene; octyl-', & 
                                  'LUMPED', & 
                                  'DTXSID2062240', & 
                                  'CCCCCCCCC1=CC=CC=C1', & 
                                  'GC',  &
                                   190.330)

            MODEL_SPECIES(   143 ) = MODEL_SPECIES_INFO('VROCP6ARO', & 
                                  '1-Hexyl-4-methylbenzene', & 
                                  'LUMPED', & 
                                  'DTXSID30333914', & 
                                  'CCCCCCC1=CC=C(C)C=C1', & 
                                  'GC',  &
                                   176.300)

            MODEL_SPECIES(   144 ) = MODEL_SPECIES_INFO('VROCP5AROP', & 
                                  '{8-hydroxy-5-octyl-6;7-dioxabicyclo[3.2.1]oct-3-en-2-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCC1(OO2)C=CC(O[O])C2C1O', & 
                                  'GC',  &
                                   271.330)

            MODEL_SPECIES(   145 ) = MODEL_SPECIES_INFO('VROCP6AROP', & 
                                  '{2-hexyl-8-hydroxy-5-methyl-6;7-dioxabicyclo[3.2.1]oct-3-en-2-yl}peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC1C2C(CCCCCC)(O[O])C=CC1(C)OO2', & 
                                  'GC',  &
                                   257.300)

            MODEL_SPECIES(   146 ) = MODEL_SPECIES_INFO('VROCN2ALK', & 
                                  'Triacontane', & 
                                  'LUMPED', & 
                                  'DTXSID0060935', & 
                                  'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', & 
                                  'GC',  &
                                   422.830)

            MODEL_SPECIES(   147 ) = MODEL_SPECIES_INFO('VROCN1ALK', & 
                                  '5;9-Dimethylheptacosane', & 
                                  'LUMPED', & 
                                  'DTXSID40823452', & 
                                  'CCCCCCCCCCCCCCCCCCC(C)CCCC(C)CCCC', & 
                                  'GC',  &
                                   408.800)

            MODEL_SPECIES(   148 ) = MODEL_SPECIES_INFO('VROCP0ALK', & 
                                  '11-Methylheptacosane', & 
                                  'LUMPED', & 
                                  'DTXSID40333900', & 
                                  'CCCCCCCCCCCCCCCCC(C)CCCCCCCCCC', & 
                                  'GC',  &
                                   394.770)

            MODEL_SPECIES(   149 ) = MODEL_SPECIES_INFO('VROCP1ALK', & 
                                  'Heptacosane', & 
                                  'LUMPED', & 
                                  'DTXSID6058637', & 
                                  'CCCCCCCCCCCCCCCCCCCCCCCCCCC', & 
                                  'GC',  &
                                   380.750)

            MODEL_SPECIES(   150 ) = MODEL_SPECIES_INFO('VROCP2ALK', & 
                                  'Tetracosane', & 
                                  'LUMPED', & 
                                  'DTXSID8060955', & 
                                  'CCCCCCCCCCCCCCCCCCCCCCCC', & 
                                  'GC',  &
                                   338.660)

            MODEL_SPECIES(   151 ) = MODEL_SPECIES_INFO('VROCP3ALK', & 
                                  'Heneicosane', & 
                                  'LUMPED', & 
                                  'DTXSID9047097', & 
                                  'CCCCCCCCCCCCCCCCCCCCC', & 
                                  'GC',  &
                                   296.580)

            MODEL_SPECIES(   152 ) = MODEL_SPECIES_INFO('VROCP4ALK', & 
                                  'Octadecane', & 
                                  'LUMPED', & 
                                  'DTXSID9047172', & 
                                  'CCCCCCCCCCCCCCCCCC', & 
                                  'GC',  &
                                   254.500)

            MODEL_SPECIES(   153 ) = MODEL_SPECIES_INFO('VROCP5ALK', & 
                                  'Tetradecane', & 
                                  'LUMPED', & 
                                  'DTXSID1027267', & 
                                  'CCCCCCCCCCCCCC', & 
                                  'GC',  &
                                   198.390)

            MODEL_SPECIES(   154 ) = MODEL_SPECIES_INFO('VROCP6ALK', & 
                                  'Tridecane', & 
                                  'LUMPED', & 
                                  'DTXSID6027266', & 
                                  'CCCCCCCCCCCCC', & 
                                  'GC',  &
                                   184.370)

            MODEL_SPECIES(   155 ) = MODEL_SPECIES_INFO('VROCP1ALKP', & 
                                  'heptacosan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCCCCCCCCCCCCCC(CC)O[O]', & 
                                  'GC',  &
                                   411.740)

            MODEL_SPECIES(   156 ) = MODEL_SPECIES_INFO('VROCP2ALKP', & 
                                  'tetracosan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCCCCCCCCCCC(CC)O[O]', & 
                                  'GC',  &
                                   369.650)

            MODEL_SPECIES(   157 ) = MODEL_SPECIES_INFO('VROCP3ALKP', & 
                                  'henicosan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCCCCCCCC(CC)O[O]', & 
                                  'GC',  &
                                   327.570)

            MODEL_SPECIES(   158 ) = MODEL_SPECIES_INFO('VROCP4ALKP', & 
                                  'octadecan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCCCCC(CC)O[O]', & 
                                  'GC',  &
                                   285.490)

            MODEL_SPECIES(   159 ) = MODEL_SPECIES_INFO('VROCP5ALKP', & 
                                  'tetradecan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCC(CC)O[O]', & 
                                  'GC',  &
                                   229.380)

            MODEL_SPECIES(   160 ) = MODEL_SPECIES_INFO('VROCP6ALKP', & 
                                  'tridecan-3-ylperoxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCC(CC)O[O]', & 
                                  'GC',  &
                                   215.360)

            MODEL_SPECIES(   161 ) = MODEL_SPECIES_INFO('VROCP1ALKP2', & 
                                  '(3-hydroxyheptacosan-6-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCCCCCCCCCCC(O[O])CCC(O)CC', & 
                                  'GC',  &
                                   427.730)

            MODEL_SPECIES(   162 ) = MODEL_SPECIES_INFO('VROCP2ALKP2', & 
                                  '(3-hydroxytetracosan-6-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCCCCCCCC(O[O])CCC(O)CC', & 
                                  'GC',  &
                                   385.650)

            MODEL_SPECIES(   163 ) = MODEL_SPECIES_INFO('VROCP3ALKP2', & 
                                  '(3-hydroxyhenicosan-6-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCCCCC(O[O])CCC(O)CC', & 
                                  'GC',  &
                                   343.570)

            MODEL_SPECIES(   164 ) = MODEL_SPECIES_INFO('VROCP4ALKP2', & 
                                  '(3-hydroxyoctadecan-6-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCCCCCC(O[O])CCC(O)CC', & 
                                  'GC',  &
                                   301.490)

            MODEL_SPECIES(   165 ) = MODEL_SPECIES_INFO('VROCP5ALKP2', & 
                                  '(3-hydroxytetradecan-6-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCCC(O[O])CCC(O)CC', & 
                                  'GC',  &
                                   245.380)

            MODEL_SPECIES(   166 ) = MODEL_SPECIES_INFO('VROCP6ALKP2', & 
                                  '(3-hydroxytridecan-6-yl)peroxy', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CCCCCCCC(O[O])CCC(O)CC', & 
                                  'GC',  &
                                   231.360)

            MODEL_SPECIES(   167 ) = MODEL_SPECIES_INFO('VROCN2OXY2', & 
                                  '12(S)-hydroxy-16-Heptadecynoic Acid', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C#CCCC[C@H](CCCCCCCCCCC(=O)O)O', & 
                                  'GC',  &
                                   282.400)

            MODEL_SPECIES(   168 ) = MODEL_SPECIES_INFO('VROCN2OXY4', & 
                                  '2-HYDROXYUNDECANEDIOIC ACID', & 
                                  'LUMPED', & 
                                  'DTXSID90726525', & 
                                  'C(CCCCC(=O)O)CCCC(C(=O)O)O', & 
                                  'GC',  &
                                   232.300)

            MODEL_SPECIES(   169 ) = MODEL_SPECIES_INFO('VROCN2OXY8', & 
                                  '3;4;5;6;7-Pentahydroxyheptan-2-one', & 
                                  'LUMPED', & 
                                  'DTXSID80956455', & 
                                  'CC(=O)C(C(C(C(CO)O)O)O)O', & 
                                  'GC',  &
                                   194.200)

            MODEL_SPECIES(   170 ) = MODEL_SPECIES_INFO('VROCN1OXY1', & 
                                  'Arachidic Acid', & 
                                  'LUMPED', & 
                                  'DTXSID1060134', & 
                                  'CCCCCCCCCCCCCCCCCCCC(=O)O', & 
                                  'GC',  &
                                   312.500)

            MODEL_SPECIES(   171 ) = MODEL_SPECIES_INFO('VROCN1OXY3', & 
                                  'DODECANEDIOIC ACID', & 
                                  'LUMPED', & 
                                  'DTXSID3027297', & 
                                  'C(CCCCCC(=O)O)CCCCC(=O)O', & 
                                  'GC',  &
                                   230.300)

            MODEL_SPECIES(   172 ) = MODEL_SPECIES_INFO('VROCN1OXY6', & 
                                  '2-hydroxy-octanedioic acid', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C(CCC(C(=O)O)O)CCC(=O)O', & 
                                  'GC',  &
                                   190.200)

            MODEL_SPECIES(   173 ) = MODEL_SPECIES_INFO('VROCP0OXY2', & 
                                  '3-oxotetradecanoic acid', & 
                                  'LUMPED', & 
                                  'DTXSID10332384', & 
                                  'CCCCCCCCCCCC(=O)CC(=O)O', & 
                                  'GC',  &
                                   242.400)

            MODEL_SPECIES(   174 ) = MODEL_SPECIES_INFO('VROCP0OXY4', & 
                                  'DECANEDIOIC ACID', & 
                                  'LUMPED', & 
                                  'DTXSID7026867', & 
                                  'C(CCCCC(=O)O)CCCC(=O)O', & 
                                  'GC',  &
                                   202.300)

            MODEL_SPECIES(   175 ) = MODEL_SPECIES_INFO('VROCP1OXY1', & 
                                  'HEPTADECANOIC ACID', & 
                                  'LUMPED', & 
                                  'DTXSID5021596', & 
                                  'CCCCCCCCCCCCCCCCC(=O)O', & 
                                  'GC',  &
                                   270.500)

            MODEL_SPECIES(   176 ) = MODEL_SPECIES_INFO('VROCP1OXY3', & 
                                  '11-hydroxyundecanoic acid', & 
                                  'LUMPED', & 
                                  'DTXSID40190136', & 
                                  'C(CCCCCO)CCCCC(=O)O', & 
                                  'GC',  &
                                   202.300)

            MODEL_SPECIES(   177 ) = MODEL_SPECIES_INFO('VROCP2OXY2', & 
                                  'Dodecanoic Acid', & 
                                  'LUMPED', & 
                                  'DTXSID5021590', & 
                                  'CCCCCCCCCCCC(=O)O', & 
                                  'GC',  &
                                   200.300)

            MODEL_SPECIES(   178 ) = MODEL_SPECIES_INFO('VROCP3OXY2', & 
                                  '11-hydroxyundecanal', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C(CCCCCO)CCCCC=O', & 
                                  'GC',  &
                                   186.300)

            MODEL_SPECIES(   179 ) = MODEL_SPECIES_INFO('VROCP4OXY2', & 
                                  '2-Ethylheptanoic Acid', & 
                                  'LUMPED', & 
                                  'DTXSID40880929', & 
                                  'CCCCCC(CC)C(=O)O', & 
                                  'GC',  &
                                   158.200)

            MODEL_SPECIES(   180 ) = MODEL_SPECIES_INFO('VROCP5OXY1', & 
                                  'Undecanal', & 
                                  'LUMPED', & 
                                  'DTXSID4021688', & 
                                  'CCCCCCCCCCC=O', & 
                                  'GC',  &
                                   170.300)

            MODEL_SPECIES(   181 ) = MODEL_SPECIES_INFO('VROCP6OXY1', & 
                                  'Nonanal', & 
                                  'LUMPED', & 
                                  'DTXSID9021639', & 
                                  'CCCCCCCCC=O', & 
                                  'GC',  &
                                   142.200)

            MODEL_SPECIES(   182 ) = MODEL_SPECIES_INFO('ECH4', & 
                                  'Methane', & 
                                  'EXPLICIT', & 
                                  'DTXSID8025545', & 
                                  'C', & 
                                  'GC',  &
                                    16.000)

            MODEL_SPECIES(   183 ) = MODEL_SPECIES_INFO('CO2', & 
                                  'Carbon dioxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID4027028', & 
                                  'O=C=O', & 
                                  'GC',  &
                                    44.000)

            MODEL_SPECIES(   184 ) = MODEL_SPECIES_INFO('VMTN1', & 
                                  '3-[carboxy(hydroxy)methyl]cyclobutane-1;2;2-tricarboxylic acid', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC(C1CC(C(O)=O)C1(C(O)=O)C(O)=O)C(O)=O', & 
                                  'GC',  &
                                   262.000)

            MODEL_SPECIES(   185 ) = MODEL_SPECIES_INFO('CL2', & 
                                  'Chlorine', & 
                                  'EXPLICIT', & 
                                  'DTXSID1020273', & 
                                  'ClCl', & 
                                  'GC',  &
                                    71.000)

            MODEL_SPECIES(   186 ) = MODEL_SPECIES_INFO('CL', & 
                                  'Chlorine atom', & 
                                  'EXPLICIT', & 
                                  'DTXSID801014230', & 
                                  '[Cl]', & 
                                  'GC',  &
                                    35.500)

            MODEL_SPECIES(   187 ) = MODEL_SPECIES_INFO('CLO', & 
                                  'Chlorine monoxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID8073136', & 
                                  '[O-]Cl', & 
                                  'GC',  &
                                    51.500)

            MODEL_SPECIES(   188 ) = MODEL_SPECIES_INFO('OCLO', & 
                                  'Chlorine dioxide isomar', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'O=Cl[O]', & 
                                  'GC',  &
                                    67.500)

            MODEL_SPECIES(   189 ) = MODEL_SPECIES_INFO('CL2O2', & 
                                  'Dichlorine dioxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID001310193', & 
                                  'ClOOCl', & 
                                  'GC',  &
                                   103.000)

            MODEL_SPECIES(   190 ) = MODEL_SPECIES_INFO('HOCL', & 
                                  'Hypochlorous acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID3036737', & 
                                  'OCl', & 
                                  'GC',  &
                                    52.500)

            MODEL_SPECIES(   191 ) = MODEL_SPECIES_INFO('CLNO', & 
                                  'Nitrosyl chloride', & 
                                  'EXPLICIT', & 
                                  'DTXSID5051945', & 
                                  'ClN=O', & 
                                  'GC',  &
                                    65.500)

            MODEL_SPECIES(   192 ) = MODEL_SPECIES_INFO('CLNO2', & 
                                  'Nitryl chloride', & 
                                  'EXPLICIT', & 
                                  'DTXSID601317066', & 
                                  'ClN(=O)=O', & 
                                  'GC',  &
                                    81.500)

            MODEL_SPECIES(   193 ) = MODEL_SPECIES_INFO('CLNO3', & 
                                  'Chlorine nitrate', & 
                                  'EXPLICIT', & 
                                  'DTXSID00163043', & 
                                  '[O-][N+](=O)OCl', & 
                                  'GC',  &
                                    97.500)

            MODEL_SPECIES(   194 ) = MODEL_SPECIES_INFO('HCOCL', & 
                                  'Formyl chloride', & 
                                  'EXPLICIT', & 
                                  'DTXSID30180344', & 
                                  'ClC=O', & 
                                  'GC',  &
                                    64.500)

            MODEL_SPECIES(   195 ) = MODEL_SPECIES_INFO('HCL', & 
                                  'Hydrochloric acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID2020711', & 
                                  'Cl', & 
                                  'GC',  &
                                    36.000)

            MODEL_SPECIES(   196 ) = MODEL_SPECIES_INFO('CLOO', & 
                                  'Chlorine dioxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID5023958', & 
                                  'O=Cl[O]', & 
                                  'GC',  &
                                    67.500)

            MODEL_SPECIES(   197 ) = MODEL_SPECIES_INFO('BR2', & 
                                  'Bromine', & 
                                  'EXPLICIT', & 
                                  'DTXSID1035238', & 
                                  'BrBr', & 
                                  'GC',  &
                                   159.800)

            MODEL_SPECIES(   198 ) = MODEL_SPECIES_INFO('BR', & 
                                  'Bromine atom', & 
                                  'EXPLICIT', & 
                                  'DTXSID201014232', & 
                                  '[Br]', & 
                                  'GC',  &
                                    79.900)

            MODEL_SPECIES(   199 ) = MODEL_SPECIES_INFO('BRO', & 
                                  'Bromine monoxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID201315575', & 
                                  '[O]Br', & 
                                  'GC',  &
                                    95.900)

            MODEL_SPECIES(   200 ) = MODEL_SPECIES_INFO('OBRO', & 
                                  'Bromine dioxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID401028281', & 
                                  'O=[Br]=O', & 
                                  'GC',  &
                                   111.900)

            MODEL_SPECIES(   201 ) = MODEL_SPECIES_INFO('HOBR', & 
                                  'Hypobromous acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID701024864', & 
                                  'OBr', & 
                                  'GC',  &
                                    96.900)

            MODEL_SPECIES(   202 ) = MODEL_SPECIES_INFO('BRNO', & 
                                  'Nitrosyl bromide', & 
                                  'EXPLICIT', & 
                                  'DTXSID90158701', & 
                                  'BrN=O', & 
                                  'GC',  &
                                   109.900)

            MODEL_SPECIES(   203 ) = MODEL_SPECIES_INFO('BRNO2', & 
                                  'Nitryl bromide', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'BrN(=O)=O', & 
                                  'GC',  &
                                   125.900)

            MODEL_SPECIES(   204 ) = MODEL_SPECIES_INFO('BRNO3', & 
                                  'Bromine nitrate', & 
                                  'EXPLICIT', & 
                                  'DTXSID90960821', & 
                                  'BrON(=O)=O', & 
                                  'GC',  &
                                   141.900)

            MODEL_SPECIES(   205 ) = MODEL_SPECIES_INFO('CH2BR2', & 
                                  'Dibromomethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID4021557', & 
                                  'BrCBr', & 
                                  'GC',  &
                                   173.800)

            MODEL_SPECIES(   206 ) = MODEL_SPECIES_INFO('CHBR3', & 
                                  'Bromoform', & 
                                  'EXPLICIT', & 
                                  'DTXSID1021374', & 
                                  'BrC(Br)Br', & 
                                  'GC',  &
                                   252.700)

            MODEL_SPECIES(   207 ) = MODEL_SPECIES_INFO('HBR', & 
                                  'Hydobromic acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID0029713', & 
                                  'Br', & 
                                  'GC',  &
                                    80.900)

            MODEL_SPECIES(   208 ) = MODEL_SPECIES_INFO('HCOBR', & 
                                  'Formyl bromide', & 
                                  'EXPLICIT', & 
                                  'DTXSID20227994', & 
                                  'BrC=O', & 
                                  'GC',  &
                                   108.900)

            MODEL_SPECIES(   209 ) = MODEL_SPECIES_INFO('I2', & 
                                  'Iodine', & 
                                  'EXPLICIT', & 
                                  'DTXSID7034672', & 
                                  'II', & 
                                  'GC',  &
                                   253.800)

            MODEL_SPECIES(   210 ) = MODEL_SPECIES_INFO('I', & 
                                  'Atmomic iodine', & 
                                  'EXPLICIT', & 
                                  'DTXSID501014231', & 
                                  '[I]', & 
                                  'GC',  &
                                   126.900)

            MODEL_SPECIES(   211 ) = MODEL_SPECIES_INFO('IO', & 
                                  'Iodine monoxide', & 
                                  'EXPLICIT', & 
                                  'DTXSID101316449', & 
                                  '[O]I', & 
                                  'GC',  &
                                   142.900)

            MODEL_SPECIES(   212 ) = MODEL_SPECIES_INFO('OIO', & 
                                  'Iodine dioxide', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'O=I(=O)', & 
                                  'GC',  &
                                   158.900)

            MODEL_SPECIES(   213 ) = MODEL_SPECIES_INFO('I2O2', & 
                                  'Diiodine dioxide', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'O=I(=O)I', & 
                                  'GC',  &
                                   285.800)

            MODEL_SPECIES(   214 ) = MODEL_SPECIES_INFO('HOI', & 
                                  'Hypoiodous acid', & 
                                  'EXPLICIT', & 
                                  'DTXSID8042050', & 
                                  'OI', & 
                                  'GC',  &
                                   143.900)

            MODEL_SPECIES(   215 ) = MODEL_SPECIES_INFO('HI', & 
                                  'Hydrogen iodide', & 
                                  'EXPLICIT', & 
                                  'DTXSID2044349', & 
                                  'I', & 
                                  'GC',  &
                                   127.900)

            MODEL_SPECIES(   216 ) = MODEL_SPECIES_INFO('INO', & 
                                  'Nitrosyl iodide', & 
                                  'EXPLICIT', & 
                                  'DTXSID80207297', & 
                                  'IN=O', & 
                                  'GC',  &
                                   156.900)

            MODEL_SPECIES(   217 ) = MODEL_SPECIES_INFO('INO2', & 
                                  'Iodine nitrite', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'IN(=O)=O', & 
                                  'GC',  &
                                   172.900)

            MODEL_SPECIES(   218 ) = MODEL_SPECIES_INFO('INO3', & 
                                  'Iodine nitrate', & 
                                  'EXPLICIT', & 
                                  'DTXSID001336624', & 
                                  'ION(=O)=O', & 
                                  'GC',  &
                                   188.900)

            MODEL_SPECIES(   219 ) = MODEL_SPECIES_INFO('CH3I', & 
                                  'Iodomethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID0024187', & 
                                  'CI', & 
                                  'GC',  &
                                   141.900)

            MODEL_SPECIES(   220 ) = MODEL_SPECIES_INFO('CH2I2', & 
                                  'Methylene iodide', & 
                                  'EXPLICIT', & 
                                  'DTXSID4058784', & 
                                  'ICI', & 
                                  'GC',  &
                                   267.800)

            MODEL_SPECIES(   221 ) = MODEL_SPECIES_INFO('BRCL', & 
                                  'Bromine chloride', & 
                                  'EXPLICIT', & 
                                  'DTXSID4035259', & 
                                  'BrCl', & 
                                  'GC',  &
                                   115.400)

            MODEL_SPECIES(   222 ) = MODEL_SPECIES_INFO('ICL', & 
                                  'Iodine chloride', & 
                                  'EXPLICIT', & 
                                  'DTXSID1064879', & 
                                  'ICl', & 
                                  'GC',  &
                                   162.400)

            MODEL_SPECIES(   223 ) = MODEL_SPECIES_INFO('IBR', & 
                                  'Iodine bromide', & 
                                  'EXPLICIT', & 
                                  'DTXSID2064862', & 
                                  'IBr', & 
                                  'GC',  &
                                   206.800)

            MODEL_SPECIES(   224 ) = MODEL_SPECIES_INFO('CH2IBR', & 
                                  'Bromoiodomethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID50204233', & 
                                  'ICBr', & 
                                  'GC',  &
                                   220.800)

            MODEL_SPECIES(   225 ) = MODEL_SPECIES_INFO('CH2ICL', & 
                                  'Chloroiodomethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID50208034', & 
                                  'ClCI', & 
                                  'GC',  &
                                   176.400)

            MODEL_SPECIES(   226 ) = MODEL_SPECIES_INFO('CHBR2CL', & 
                                  'Dibromochloromethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID1020300', & 
                                  'ClC(Br)Br', & 
                                  'GC',  &
                                   208.300)

            MODEL_SPECIES(   227 ) = MODEL_SPECIES_INFO('CHBRCL2', & 
                                  'Bromodichloromethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID1020198', & 
                                  'ClC(Cl)Br', & 
                                  'GC',  &
                                   163.900)

            MODEL_SPECIES(   228 ) = MODEL_SPECIES_INFO('CH2BRCL', & 
                                  'Bromochloromethane', & 
                                  'EXPLICIT', & 
                                  'DTXSID4021503', & 
                                  'ClCBr', & 
                                  'GC',  &
                                   129.400)

            MODEL_SPECIES(   229 ) = MODEL_SPECIES_INFO('NO2PIJ', & 
                                  'Dioxidonitrogen', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  '[N+](=O)=O', & 
                                  'GC',  &
                                    46.000)

            MODEL_SPECIES(   230 ) = MODEL_SPECIES_INFO('NO2PK', & 
                                  'Dioxidonitrogen', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  '[N+](=O)=O', & 
                                  'GC',  &
                                    46.000)

            MODEL_SPECIES(   231 ) = MODEL_SPECIES_INFO('I2O3', & 
                                  'Diiodine trioxide', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'O=IOI(=O)', & 
                                  'GC',  &
                                   301.800)

            MODEL_SPECIES(   232 ) = MODEL_SPECIES_INFO('I2O4', & 
                                  'Diiodine Tetroxide', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  'O=IOI(=O)=O', & 
                                  'GC',  &
                                   317.800)

            MODEL_SPECIES(   233 ) = MODEL_SPECIES_INFO('ASO4I', & 
                                  'Sulfate ion', & 
                                  'LUMPED', & 
                                  'DTXSID3042425', & 
                                  '[O-]S(=O)(=O)[O-]', & 
                                  'NA',  &
                                    96.000)

            MODEL_SPECIES(   234 ) = MODEL_SPECIES_INFO('ASO4J', & 
                                  'Sulfate ion', & 
                                  'LUMPED', & 
                                  'DTXSID3042425', & 
                                  '[O-]S(=O)(=O)[O-]', & 
                                  'NA',  &
                                    96.000)

            MODEL_SPECIES(   235 ) = MODEL_SPECIES_INFO('ASO4K', & 
                                  'Sulfate ion', & 
                                  'LUMPED', & 
                                  'DTXSID3042425', & 
                                  '[O-]S(=O)(=O)[O-]', & 
                                  'NA',  &
                                    96.000)

            MODEL_SPECIES(   236 ) = MODEL_SPECIES_INFO('ANH4I', & 
                                  'Ammonium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID5043974', & 
                                  '[NH4+]', & 
                                  'NA',  &
                                    18.000)

            MODEL_SPECIES(   237 ) = MODEL_SPECIES_INFO('ANH4J', & 
                                  'Ammonium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID5043974', & 
                                  '[NH4+]', & 
                                  'NA',  &
                                    18.000)

            MODEL_SPECIES(   238 ) = MODEL_SPECIES_INFO('ANH4K', & 
                                  'Ammonium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID5043974', & 
                                  '[NH4+]', & 
                                  'NA',  &
                                    18.000)

            MODEL_SPECIES(   239 ) = MODEL_SPECIES_INFO('ANO3I', & 
                                  'Nitrate ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID5024217', & 
                                  '[N+](=O)([O-])[O-]', & 
                                  'NA',  &
                                    62.000)

            MODEL_SPECIES(   240 ) = MODEL_SPECIES_INFO('ANO3J', & 
                                  'Nitrate ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID5024217', & 
                                  '[N+](=O)([O-])[O-]', & 
                                  'NA',  &
                                    62.000)

            MODEL_SPECIES(   241 ) = MODEL_SPECIES_INFO('ANO3K', & 
                                  'Nitrate ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID5024217', & 
                                  '[N+](=O)([O-])[O-]', & 
                                  'NA',  &
                                    62.000)

            MODEL_SPECIES(   242 ) = MODEL_SPECIES_INFO('ASOATJ', & 
                                  '3;4;5;6;7-Pentahydroxyheptan-2-one', & 
                                  'LUMPED', & 
                                  'DTXSID80956455', & 
                                  'CC(=O)C(C(C(C(CO)O)O)O)O', & 
                                  'NA',  &
                                   200.000)

            MODEL_SPECIES(   243 ) = MODEL_SPECIES_INFO('AGLYOLIGJ', & 
                                  '2-(4;5-dihydroxy-1;3-dioxolan-2-yl)-1;3-dioxolane-4;5-diol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC2OC(C1OC(O)C(O)O1)OC2O', & 
                                  'NA',  &
                                    66.400)

            MODEL_SPECIES(   244 ) = MODEL_SPECIES_INFO('AHOMJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   250.000)

            MODEL_SPECIES(   245 ) = MODEL_SPECIES_INFO('AELHOMJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   402.000)

            MODEL_SPECIES(   246 ) = MODEL_SPECIES_INFO('AORGCJ', & 
                                  '2-(4;5-dihydroxy-1;3-dioxolan-2-yl)-1;3-dioxolane-4;5-diol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'OC2OC(C1OC(O)C(O)O1)OC2O', & 
                                  'NA',  &
                                   177.000)

            MODEL_SPECIES(   247 ) = MODEL_SPECIES_INFO('AECI', & 
                                  'Carbon', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[C]', & 
                                  'NA',  &
                                    12.000)

            MODEL_SPECIES(   248 ) = MODEL_SPECIES_INFO('AECJ', & 
                                  'Carbon', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  '[C]', & 
                                  'NA',  &
                                    12.000)

            MODEL_SPECIES(   249 ) = MODEL_SPECIES_INFO('AOTHRI', & 
                                  'Unspeciated PM', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   200.000)

            MODEL_SPECIES(   250 ) = MODEL_SPECIES_INFO('AOTHRJ', & 
                                  'Unspeciated PM', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   200.000)

            MODEL_SPECIES(   251 ) = MODEL_SPECIES_INFO('AFEJ', & 
                                  'Iron ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID4042672', & 
                                  '[Fe++]', & 
                                  'NA',  &
                                    55.800)

            MODEL_SPECIES(   252 ) = MODEL_SPECIES_INFO('AALJ', & 
                                  'Aluminum ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID70912343', & 
                                  '[Al+3]', & 
                                  'NA',  &
                                    27.000)

            MODEL_SPECIES(   253 ) = MODEL_SPECIES_INFO('ASIJ', & 
                                  'Total Silicon', & 
                                  'EXPLICIT', & 
                                  'DTXSID0051441', & 
                                  '[Si]', & 
                                  'NA',  &
                                    28.100)

            MODEL_SPECIES(   254 ) = MODEL_SPECIES_INFO('ATIJ', & 
                                  'Total Titanium', & 
                                  'EXPLICIT', & 
                                  'DTXSID3047764', & 
                                  '[Ti]', & 
                                  'NA',  &
                                    47.900)

            MODEL_SPECIES(   255 ) = MODEL_SPECIES_INFO('ACAJ', & 
                                  'Calcium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID7037638', & 
                                  '[Ca++]', & 
                                  'NA',  &
                                    40.100)

            MODEL_SPECIES(   256 ) = MODEL_SPECIES_INFO('AMGJ', & 
                                  'Magnesium ion', & 
                                  'EXPLICIT', & 
                                  'NA', & 
                                  '[Mg++]', & 
                                  'NA',  &
                                    24.300)

            MODEL_SPECIES(   257 ) = MODEL_SPECIES_INFO('AKJ', & 
                                  'Potassium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID9042671', & 
                                  '[K+]', & 
                                  'NA',  &
                                    39.100)

            MODEL_SPECIES(   258 ) = MODEL_SPECIES_INFO('AMNJ', & 
                                  'Manganese ions', & 
                                  'LUMPED', & 
                                  'DTXSID00167687', & 
                                  '[Mn++]', & 
                                  'NA',  &
                                    54.900)

            MODEL_SPECIES(   259 ) = MODEL_SPECIES_INFO('ACORS', & 
                                  'Coarse PM', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   100.000)

            MODEL_SPECIES(   260 ) = MODEL_SPECIES_INFO('ASOIL', & 
                                  'Crustal species', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   100.000)

            MODEL_SPECIES(   261 ) = MODEL_SPECIES_INFO('NUMATKN', & 
                                  'Number of particles', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                     1.000)

            MODEL_SPECIES(   262 ) = MODEL_SPECIES_INFO('NUMACC', & 
                                  'Number of particles', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                     1.000)

            MODEL_SPECIES(   263 ) = MODEL_SPECIES_INFO('NUMCOR', & 
                                  'Number of particles', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                     1.000)

            MODEL_SPECIES(   264 ) = MODEL_SPECIES_INFO('SRFATKN', & 
                                  'Surface area of particles', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                     1.000)

            MODEL_SPECIES(   265 ) = MODEL_SPECIES_INFO('SRFACC', & 
                                  'Surface area of particles', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                     1.000)

            MODEL_SPECIES(   266 ) = MODEL_SPECIES_INFO('SRFCOR', & 
                                  'Surface area of particles', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                     1.000)

            MODEL_SPECIES(   267 ) = MODEL_SPECIES_INFO('AORGH2OJ', & 
                                  'Water', & 
                                  'EXPLICIT', & 
                                  'DTXSID6026296', & 
                                  'O', & 
                                  'NA',  &
                                    18.000)

            MODEL_SPECIES(   268 ) = MODEL_SPECIES_INFO('AH2OI', & 
                                  'Water', & 
                                  'EXPLICIT', & 
                                  'DTXSID6026296', & 
                                  'O', & 
                                  'NA',  &
                                    18.000)

            MODEL_SPECIES(   269 ) = MODEL_SPECIES_INFO('AH2OJ', & 
                                  'Water', & 
                                  'EXPLICIT', & 
                                  'DTXSID6026296', & 
                                  'O', & 
                                  'NA',  &
                                    18.000)

            MODEL_SPECIES(   270 ) = MODEL_SPECIES_INFO('AH2OK', & 
                                  'Water', & 
                                  'EXPLICIT', & 
                                  'DTXSID6026296', & 
                                  'O', & 
                                  'NA',  &
                                    18.000)

            MODEL_SPECIES(   271 ) = MODEL_SPECIES_INFO('AH3OPI', & 
                                  'Hydronium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID20893597', & 
                                  '[OH3+]', & 
                                  'NA',  &
                                    19.000)

            MODEL_SPECIES(   272 ) = MODEL_SPECIES_INFO('AH3OPJ', & 
                                  'Hydronium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID20893597', & 
                                  '[OH3+]', & 
                                  'NA',  &
                                    19.000)

            MODEL_SPECIES(   273 ) = MODEL_SPECIES_INFO('AH3OPK', & 
                                  'Hydronium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID20893597', & 
                                  '[OH3+]', & 
                                  'NA',  &
                                    19.000)

            MODEL_SPECIES(   274 ) = MODEL_SPECIES_INFO('ANAI', & 
                                  'Sodium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID8037671', & 
                                  '[Na+]', & 
                                  'NA',  &
                                    23.000)

            MODEL_SPECIES(   275 ) = MODEL_SPECIES_INFO('ANAJ', & 
                                  'Sodium ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID8037671', & 
                                  '[Na+]', & 
                                  'NA',  &
                                    23.000)

            MODEL_SPECIES(   276 ) = MODEL_SPECIES_INFO('ACLI', & 
                                  'Chloride ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID6043969', & 
                                  '[Cl-]', & 
                                  'NA',  &
                                    35.500)

            MODEL_SPECIES(   277 ) = MODEL_SPECIES_INFO('ACLJ', & 
                                  'Chloride ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID6043969', & 
                                  '[Cl-]', & 
                                  'NA',  &
                                    35.500)

            MODEL_SPECIES(   278 ) = MODEL_SPECIES_INFO('ACLK', & 
                                  'Chloride ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID6043969', & 
                                  '[Cl-]', & 
                                  'NA',  &
                                    35.500)

            MODEL_SPECIES(   279 ) = MODEL_SPECIES_INFO('ASEACAT', & 
                                  'Coarse sea spray cations', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                    23.750)

            MODEL_SPECIES(   280 ) = MODEL_SPECIES_INFO('APOCI', & 
                                  'Organic carbon', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   220.000)

            MODEL_SPECIES(   281 ) = MODEL_SPECIES_INFO('APOCJ', & 
                                  'Organic carbon', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   220.000)

            MODEL_SPECIES(   282 ) = MODEL_SPECIES_INFO('APNCOMI', & 
                                  'Non-carbon organic matter', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   220.000)

            MODEL_SPECIES(   283 ) = MODEL_SPECIES_INFO('APNCOMJ', & 
                                  'Non-carbon organic matter', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   220.000)

            MODEL_SPECIES(   284 ) = MODEL_SPECIES_INFO('AOP3J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   176.200)

            MODEL_SPECIES(   285 ) = MODEL_SPECIES_INFO('AROCN2ALKI', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   422.830)

            MODEL_SPECIES(   286 ) = MODEL_SPECIES_INFO('AROCN2ALKJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   422.830)

            MODEL_SPECIES(   287 ) = MODEL_SPECIES_INFO('AROCN1ALKI', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   408.800)

            MODEL_SPECIES(   288 ) = MODEL_SPECIES_INFO('AROCN1ALKJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   408.800)

            MODEL_SPECIES(   289 ) = MODEL_SPECIES_INFO('AROCP0ALKI', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   394.770)

            MODEL_SPECIES(   290 ) = MODEL_SPECIES_INFO('AROCP0ALKJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   394.770)

            MODEL_SPECIES(   291 ) = MODEL_SPECIES_INFO('AROCP1ALKI', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   380.750)

            MODEL_SPECIES(   292 ) = MODEL_SPECIES_INFO('AROCP1ALKJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   380.750)

            MODEL_SPECIES(   293 ) = MODEL_SPECIES_INFO('AROCP2ALKJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   338.660)

            MODEL_SPECIES(   294 ) = MODEL_SPECIES_INFO('AROCP3ALKJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   296.580)

            MODEL_SPECIES(   295 ) = MODEL_SPECIES_INFO('AROCN2OXY2I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   282.400)

            MODEL_SPECIES(   296 ) = MODEL_SPECIES_INFO('AROCN2OXY2J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   282.400)

            MODEL_SPECIES(   297 ) = MODEL_SPECIES_INFO('AROCN2OXY4I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   232.300)

            MODEL_SPECIES(   298 ) = MODEL_SPECIES_INFO('AROCN2OXY4J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   232.300)

            MODEL_SPECIES(   299 ) = MODEL_SPECIES_INFO('AROCN2OXY8I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   194.200)

            MODEL_SPECIES(   300 ) = MODEL_SPECIES_INFO('AROCN2OXY8J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   194.200)

            MODEL_SPECIES(   301 ) = MODEL_SPECIES_INFO('AROCN1OXY1I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   312.500)

            MODEL_SPECIES(   302 ) = MODEL_SPECIES_INFO('AROCN1OXY1J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   312.500)

            MODEL_SPECIES(   303 ) = MODEL_SPECIES_INFO('AROCN1OXY3I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   230.300)

            MODEL_SPECIES(   304 ) = MODEL_SPECIES_INFO('AROCN1OXY3J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   230.300)

            MODEL_SPECIES(   305 ) = MODEL_SPECIES_INFO('AROCN1OXY6I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   190.200)

            MODEL_SPECIES(   306 ) = MODEL_SPECIES_INFO('AROCN1OXY6J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   190.200)

            MODEL_SPECIES(   307 ) = MODEL_SPECIES_INFO('AROCP0OXY2I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   242.400)

            MODEL_SPECIES(   308 ) = MODEL_SPECIES_INFO('AROCP0OXY2J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   242.400)

            MODEL_SPECIES(   309 ) = MODEL_SPECIES_INFO('AROCP0OXY4I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   202.300)

            MODEL_SPECIES(   310 ) = MODEL_SPECIES_INFO('AROCP0OXY4J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   202.300)

            MODEL_SPECIES(   311 ) = MODEL_SPECIES_INFO('AROCP1OXY1I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   270.500)

            MODEL_SPECIES(   312 ) = MODEL_SPECIES_INFO('AROCP1OXY1J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   270.500)

            MODEL_SPECIES(   313 ) = MODEL_SPECIES_INFO('AROCP1OXY3I', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   202.300)

            MODEL_SPECIES(   314 ) = MODEL_SPECIES_INFO('AROCP1OXY3J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   202.300)

            MODEL_SPECIES(   315 ) = MODEL_SPECIES_INFO('AROCP2OXY2J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   200.300)

            MODEL_SPECIES(   316 ) = MODEL_SPECIES_INFO('AROCP3OXY2J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   186.300)

            MODEL_SPECIES(   317 ) = MODEL_SPECIES_INFO('AISO3NOSJ', & 
                                  '2-methylbutane-1;2;3;4-tetrol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C(O)C(O)(C)C(O)CO', & 
                                  'NA',  &
                                   136.200)

            MODEL_SPECIES(   318 ) = MODEL_SPECIES_INFO('AISO3OSJ', & 
                                  '(1;3;4-trihydroxy-2-methylbutan-2-yl)oxysulfonic acid', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'C(O)C(OS(O)(=O)(=O))(C)C(O)CO', & 
                                  'NA',  &
                                   216.200)

            MODEL_SPECIES(   319 ) = MODEL_SPECIES_INFO('AISO4J', & 
                                  '4-hydroperoxy-3-methylbutane-1;1;2;3-tetrol', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(O)(COO)C(O)C(O)O', & 
                                  'NA',  &
                                   168.100)

            MODEL_SPECIES(   320 ) = MODEL_SPECIES_INFO('AISO5J', & 
                                  '2-hydroperoxy-3;4-dihydroxy-3-methylbutanal', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'CC(O)(CO)C(OO)C=O', & 
                                  'NA',  &
                                   150.100)

            MODEL_SPECIES(   321 ) = MODEL_SPECIES_INFO('ATRPNJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   215.000)

            MODEL_SPECIES(   322 ) = MODEL_SPECIES_INFO('AHONITJ', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   265.000)

            MODEL_SPECIES(   323 ) = MODEL_SPECIES_INFO('ASEASTJ', & 
                                  'Sea salt tracer', & 
                                  'LUMPED', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                    31.300)

            MODEL_SPECIES(   324 ) = MODEL_SPECIES_INFO('AMTN1J', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA', & 
                                  'NA',  &
                                   262.000)

            MODEL_SPECIES(   325 ) = MODEL_SPECIES_INFO('ABRJ', & 
                                  'Bromide ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID6043967', & 
                                  '[Br-]', & 
                                  'NA',  &
                                    79.900)

            MODEL_SPECIES(   326 ) = MODEL_SPECIES_INFO('ABRK', & 
                                  'Bromide ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID6043967', & 
                                  '[Br-]', & 
                                  'NA',  &
                                    79.900)

            MODEL_SPECIES(   327 ) = MODEL_SPECIES_INFO('AHMSJ', & 
                                  'Hydroxymethanesulfonate ion', & 
                                  'EXPLICIT', & 
                                  'DTXSID20225910', & 
                                  'C(O)S(=O)(=O)[O-]', & 
                                  'NA',  &
                                   111.100)

            MODEL_SPECIES(   328 ) = MODEL_SPECIES_INFO('NH3', & 
                                  'Ammonia', & 
                                  'EXPLICIT', & 
                                  'DTXSID0023872', & 
                                  'N', & 
                                  'NR',  &
                                    17.000)
             SET_ATOMS_COUNTS = .TRUE.

         END FUNCTION SET_ATOMS_COUNTS
      END MODULE MODEL_ATOM_COUNTS
