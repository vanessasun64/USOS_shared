       MODULE RXNS_FUNCTION


       IMPLICIT NONE



! Name of Mechanism CRACMM3M

       PUBLIC             :: CALC_RCONST, SPECIAL_RATES, MAP_CHEMISTRY_SPECIES

       CONTAINS


       REAL( 8 ) FUNCTION POWER_T02( TEMPOT300,A0,B0 )
         IMPLICIT NONE
! rate constant for CMAQ Arrhenuis reaction type 2
! Arguements:
         REAL( 8 ), INTENT( IN ) :: TEMPOT300
         REAL( 8 ), INTENT( IN ) :: A0
         REAL( 8 ), INTENT( IN ) :: B0
         ! Local: None
         POWER_T02 =  A0 * TEMPOT300**B0
         RETURN
       END FUNCTION POWER_T02
       REAL( 8 ) FUNCTION ARRHENUIS_T04( INV_TEMP,TEMPOT300,A0,B0,C0 )
         IMPLICIT NONE
! rate constant for CMAQ Arrhenuis reaction type 4
! Arguements:
         REAL( 8 ), INTENT( IN ) :: INV_TEMP
         REAL( 8 ), INTENT( IN ) :: TEMPOT300
         REAL( 8 ), INTENT( IN ) :: A0
         REAL( 8 ), INTENT( IN ) :: B0
         REAL( 8 ), INTENT( IN ) :: C0
         ! Local:
         INTRINSIC DEXP
         ARRHENUIS_T04 =  A0 * DEXP( B0 * INV_TEMP ) * TEMPOT300**C0
         RETURN
       END FUNCTION ARRHENUIS_T04
       REAL( 8 ) FUNCTION ARRHENUIS_T03( INV_TEMP,A0,B0 )
! rate constant for CMAQ Arrhenuis reaction type 3
         IMPLICIT NONE
! Arguements:
         REAL( 8 ),   INTENT( IN ) ::  INV_TEMP
         REAL( 8 ),     INTENT(IN) ::  A0
         REAL( 8 ),     INTENT(IN) ::  B0
         ! Local:
         INTRINSIC DEXP
         ARRHENUIS_T03 =  A0 * DEXP( B0 * INV_TEMP )
         RETURN
       END FUNCTION ARRHENUIS_T03 
       REAL( 8 ) FUNCTION FALLOFF_T08(INV_TEMP,CAIR,A0,C0,A2,C2,A3,C3)
! rate constant for CMAQ fall off reaction type 8
         IMPLICIT NONE
! Arguements:
         REAL( 8 ), INTENT( IN ) :: INV_TEMP
         REAL( 8 ), INTENT( IN ) :: CAIR
         REAL( 8 ), INTENT( IN ) :: A0
         REAL( 8 ), INTENT( IN ) :: C0
         REAL( 8 ), INTENT( IN ) :: A2
         REAL( 8 ), INTENT( IN ) :: C2
         REAL( 8 ), INTENT( IN ) :: A3
         REAL( 8 ), INTENT( IN ) :: C3
         ! Local:
         REAL( 8 ) K0
         REAL( 8 ) K2
         REAL( 8 ) K3
         INTRINSIC DEXP
         K0 = A0 * DEXP( C0 * INV_TEMP )
         K2 = A2 * DEXP( C2 * INV_TEMP )
         K3 = A3 * DEXP( C3 * INV_TEMP )
         K3 = K3 * CAIR
         FALLOFF_T08 = K0 + K3/( 1.0D0 + K3/K2 )
         RETURN
       END FUNCTION FALLOFF_T08
       REAL( 8 ) FUNCTION FALLOFF_T09(INV_TEMP,CAIR,A1,C1,A2,C2)
! rate constant for CMAQ fall off reaction type 9
         IMPLICIT NONE
! Arguements:
         REAL( 8 ), INTENT( IN ) :: INV_TEMP
         REAL( 8 ), INTENT( IN ) :: CAIR
         REAL( 8 ), INTENT( IN ) :: A1
         REAL( 8 ), INTENT( IN ) :: C1
         REAL( 8 ), INTENT( IN ) :: A2
         REAL( 8 ), INTENT( IN ) :: C2
         !  Local:
         REAL( 8 ) K1
         REAL( 8 ) K2
         INTRINSIC DEXP
         K1 = A1 * DEXP( C1 * INV_TEMP )
         K2 = A2 * DEXP( C2 * INV_TEMP )
         FALLOFF_T09 = K1 + K2 * CAIR
         RETURN
       END FUNCTION FALLOFF_T09
       REAL( 8 ) FUNCTION FALLOFF_T10(INV_TEMP,TEMPOT300,CAIR,A0,B0,C0,A1,B1,C1,CE,CF)
         IMPLICIT NONE
! rate constant for CMAQ fall off reaction type 10
! Arguements:
         REAL( 8 ), INTENT( IN ) :: INV_TEMP
         REAL( 8 ), INTENT( IN ) :: TEMPOT300
         REAL( 8 ), INTENT( IN ) :: CAIR
         REAL( 8 ), INTENT( IN ) :: A0
         REAL( 8 ), INTENT( IN ) :: B0
         REAL( 8 ), INTENT( IN ) :: C0
         REAL( 8 ), INTENT( IN ) :: A1
         REAL( 8 ), INTENT( IN ) :: B1
         REAL( 8 ), INTENT( IN ) :: C1
         REAL( 8 ), INTENT( IN ) :: CE
         REAL( 8 ), INTENT( IN ) :: CF
         ! Local:
         REAL( 8 ) K0
         REAL( 8 ) K1
         REAL( 8 ) KEND
         K0 = A0 * CAIR * DEXP(B0*INV_TEMP)* TEMPOT300**C0
         K1 = A1 * DEXP(B1*INV_TEMP) * TEMPOT300**C1
         KEND = ( ( 1.0D0 + ( ( 1.0D0 / CE ) * DLOG10( K0 / K1 ) ) ** 2.0D0 ) )
         KEND = 1.0D0 / KEND
         FALLOFF_T10 = ( K0 / ( 1.0D0 + K0/K1 ) ) * CF ** KEND
         RETURN
       END FUNCTION FALLOFF_T10
       REAL( 8 ) FUNCTION FALLOFF_T11(INV_TEMP,TEMPOT300,CAIR,A1,B1,C1,A2, B2, C2, D1, D2)
! rate constant for CMAQ fall off reaction type 11
! actually expanded form of type 9
         IMPLICIT NONE
! Arguements:
         REAL( 8 ), INTENT( IN ) :: INV_TEMP
         REAL( 8 ), INTENT( IN ) :: TEMPOT300
         REAL( 8 ), INTENT( IN ) :: CAIR
         REAL( 8 ), INTENT( IN ) :: A1
         REAL( 8 ), INTENT( IN ) :: B1
         REAL( 8 ), INTENT( IN ) :: C1
         REAL( 8 ), INTENT( IN ) :: A2
         REAL( 8 ), INTENT( IN ) :: B2
         REAL( 8 ), INTENT( IN ) :: C2
         REAL( 8 ), INTENT( IN ) :: D1
         REAL( 8 ), INTENT( IN ) :: D2
         !  Local:
         REAL( 8 ) K1
         REAL( 8 ) K2
         REAL( 8 ) K3
         INTRINSIC DEXP
         K1 = A1 * DEXP( C1 * INV_TEMP ) * TEMPOT300**B1
         K2 = A2 * DEXP( C2 * INV_TEMP ) * TEMPOT300**B2
         K3 = D1 * DEXP( D2 * INV_TEMP )
         FALLOFF_T11 = K1 + K2 * CAIR + K3
         RETURN
       END FUNCTION FALLOFF_T11
       REAL( 8 ) FUNCTION HALOGEN_FALLOFF(PRESS,A1,B1,A2,B2,A3)
         IMPLICIT NONE
         REAL( 8 ), INTENT( IN ) :: PRESS
         REAL( 8 ), INTENT( IN ) :: A1
         REAL( 8 ), INTENT( IN ) :: B1
         REAL( 8 ), INTENT( IN ) :: A2
         REAL( 8 ), INTENT( IN ) :: B2
         REAL( 8 ), INTENT( IN ) :: A3 ! Maximum loss rate (1/sec)
         INTRINSIC DEXP
         HALOGEN_FALLOFF = A1 * DEXP( B1 * PRESS ) + A2 * DEXP( B2 * PRESS )
         HALOGEN_FALLOFF = DMIN1 (A3, HALOGEN_FALLOFF )
         RETURN
       END FUNCTION HALOGEN_FALLOFF

       SUBROUTINE SPECIAL_RATES( NUMCELLS, Y, TAIR, DENS, RKI )
! Purpose: calculate special rate operators and update
!         appropriate rate constants

       USE RXNS_DATA
       IMPLICIT NONE

! Arguments:
       INTEGER,      INTENT( IN  )   :: NUMCELLS        ! Number of cells in block 
       REAL( 8 ),    INTENT( IN )    :: Y( :, : )       ! species concs
       REAL( 8 ),    INTENT( IN )    :: TAIR( : )       ! air temperature, K 
       REAL( 8 ),    INTENT( IN )    :: DENS( : )       ! air density, Kg/m3
       REAL( 8 ),    INTENT( INOUT ) :: RKI( :, : )     ! reaction rate constant, ppm/min 
! Local:
       REAL( 8 ), PARAMETER :: DENSITY_TO_NUMBER = 2.07930D+19 ! Kg/m3 to molecules/cm3

       INTEGER   :: NCELL
       REAL( 8 ) :: TEMP
       REAL( 8 ) :: INV_TEMP
       REAL( 8 ) :: CAIR
       REAL( 8 ) :: CFACT         ! scales operator if not multiplied by RKI, cm^3/(molecule*sec) to 1/(ppm*min)
       REAL( 8 ) :: CFACT_SQU     ! scales operator if not multiplied by RKI, cm^6/(molecule^2*sec) to 1/(ppm^2*min)
! special rate operators listed below

       DO NCELL = 1, NUMCELLS
          TEMP      = TAIR( NCELL )
          INV_TEMP  = 1.0D0 / TEMP 
          CAIR      = DENSITY_TO_NUMBER * DENS( NCELL )
          CFACT     = 6.0D-05 * CAIR
          CFACT_SQU = 6.0D-11 * CAIR * CAIR


! define special rate operators


! define rate constants in terms of special rate operators 

       END DO

       RETURN
       END SUBROUTINE SPECIAL_RATES
 
       SUBROUTINE CALC_RCONST( BLKTEMP, BLKPRES, BLKH2O, RJBLK, BLKHET, LSUNLIGHT, SEAWATER, RKI, NUMCELLS )

!**********************************************************************

!  Function: To compute thermal and photolytic reaction rate
!            coefficients for each reaction.

!  Preconditions: Photolysis rates for individual species must have
!                 been calculated and stored in RJPHOT. Expects
!                 temperature in deg K, pressure in atm., water
!                 vapor in ppmV, and J-values in /min.
!  Key Subroutines/Functions Called: POWER_02, ARRHRENUIS_T0*, FALLOFF_T*, HALOGEN_FALLOFF 
!***********************************************************************




       USE RXNS_DATA

        IMPLICIT NONE  

!  Arguements: None 

        REAL( 8 ),           INTENT( IN  ) :: BLKTEMP( : )      ! temperature, deg K 
        REAL( 8 ),           INTENT( IN  ) :: BLKPRES( : )      ! pressure, Atm
        REAL( 8 ),           INTENT( IN  ) :: BLKH2O ( : )      ! water mixing ratio, ppm 
        REAL( 8 ),           INTENT( IN  ) :: RJBLK  ( :, : )   ! photolysis rates, 1/min 
        REAL( 8 ),           INTENT( IN  ) :: BLKHET ( :, : )   ! heterogeneous rate constants, ???/min
        INTEGER,             INTENT( IN  ) :: NUMCELLS          ! Number of cells in block 
        LOGICAL,             INTENT( IN  ) :: LSUNLIGHT         ! Is there sunlight? 
        REAL( 8 ),           INTENT( IN  ) :: SEAWATER( : )     ! fractional area of OPEN+SURF 
        REAL( 8 ),           INTENT( OUT ) :: RKI ( :, : )      ! reaction rate constant, ppm/min 
!..Parameters: 

        REAL( 8 ), PARAMETER :: COEF1  = 7.33981D+15     ! Molec/cc to ppm conv factor 
        REAL( 8 ), PARAMETER :: CONSTC = 0.6D+0          ! Constant for reaction type 7
        REAL( 8 ), PARAMETER :: TI300  = 1.0D+0/300.0D+0 ! reciprocal of 300 deg K
        REAL( 8 ), PARAMETER :: SFACT  = 60.D+0          ! seconds per minute 
!..External Functions: None

!..Local Variables:

        INTEGER   :: NRT           ! Loop index for reaction types 
        INTEGER   :: IRXN          ! Reaction number
        INTEGER   :: JNUM          ! J-value species # from PHOT)
        INTEGER   :: KNUM          ! Reaction # for a relative rate coeff.
        INTEGER   :: N             ! Loop index for reactions
        INTEGER   :: NCELL         ! Loop index for # of cells in the block
        REAL( 8 ) :: CAIR          ! air number density (wet) [molec/cm^3]
        REAL( 8 ) :: CFACT         ! Convertor cm^3/(molec*sec) to 1/(ppm*min)
        REAL( 8 ) :: CFACT_SQU     ! Convertor cm^6/(molec^2*sec) to 1/(ppm^2*min)
        REAL( 8 ) :: INV_CFACT     ! ppm/min to molec/(cm^3*sec)
        REAL( 8 ) :: TEMPOT300     ! temperature divided by 300 K, dimensionaless 
        REAL( 8 ) :: INV_TEMP      ! reciprocal of air temperature, K-1
        REAL( 8 ) :: INV_CAIR      ! reciprocal of air number density (wet), [cm^3/molec]
        REAL( 8 ) :: TEMP          ! air temperature, K
        REAL( 8 ) :: PRESS         ! pressure [Atm] 
        REAL( 8 ) :: INV_RFACT     ! ppm/min to molec/(cm^3*min)
        REAL( 8 ) :: RFACT_SQU     ! cm^6/(molec^2*min) to 1/(ppm^2*min)
        REAL( 8 ) :: RFACT         ! cm^3/(molec*min) to 1/(ppm*min)
        REAL( 8 ) :: H2O           ! concentration, [molec/cm^3] 

        RKI = 0.0D0 

! All rate constants converted from  molec/cm3 to ppm
! and 1/sec to 1/min

        IF( LSUNLIGHT )THEN 
            DO NCELL = 1, NUMCELLS 

!  Reaction Label R001            
                RKI( NCELL,    1) =  RJBLK( NCELL, IJ_O3O3P_NASA06 )
!  Reaction Label R002            
                RKI( NCELL,    2) =  RJBLK( NCELL, IJ_O3O1D_NASA06 )
!  Reaction Label R003            
                RKI( NCELL,    3) =  RJBLK( NCELL, IJ_H2O2_RACM2 )
!  Reaction Label R004            
                RKI( NCELL,    4) =  RJBLK( NCELL, IJ_NO2_RACM2 )
!  Reaction Label R005            
                RKI( NCELL,    5) =  RJBLK( NCELL, IJ_NO3NO_RACM2 )
!  Reaction Label R006            
                RKI( NCELL,    6) =  RJBLK( NCELL, IJ_NO3NO2_RACM2 )
!  Reaction Label R007            
                RKI( NCELL,    7) =  RJBLK( NCELL, IJ_HONO_RACM2 )
!  Reaction Label R008            
                RKI( NCELL,    8) =  RJBLK( NCELL, IJ_HNO3_RACM2 )
!  Reaction Label R009            
                RKI( NCELL,    9) =  RJBLK( NCELL, IJ_HNO4_RACM2 )
!  Reaction Label R010            
                RKI( NCELL,   10) =  RJBLK( NCELL, IJ_HCHO_MOL_JPL19 )
!  Reaction Label R011            
                RKI( NCELL,   11) =  RJBLK( NCELL, IJ_HCHO_RAD_JPL19 )
!  Reaction Label R012            
                RKI( NCELL,   12) =  RJBLK( NCELL, IJ_CH3CHO_RACM2 )
!  Reaction Label R013            
                RKI( NCELL,   13) =  RJBLK( NCELL, IJ_ALD_JPL19 )
!  Reaction Label R014            
                RKI( NCELL,   14) =  RJBLK( NCELL, IJ_CH3COCH3A_JPL19 )
!  Reaction Label R014a           
                RKI( NCELL,   15) =  RJBLK( NCELL, IJ_CH3COCH3B_JPL19 )
!  Reaction Label R015            
                RKI( NCELL,   16) =  RJBLK( NCELL, IJ_UALD_RACM2 )
!  Reaction Label TRP01           
                RKI( NCELL,   17) =  RJBLK( NCELL, IJ_ALD_JPL19 )
!  Reaction Label TRP02           
                RKI( NCELL,   18) =  RJBLK( NCELL, IJ_ALD_JPL19 )
!  Reaction Label R016            
                RKI( NCELL,   19) =  RJBLK( NCELL, IJ_MEK_JGR19 )
!  Reaction Label R017            
                RKI( NCELL,   20) =  RJBLK( NCELL, IJ_KET_JGR19 )
!  Reaction Label R018            
                RKI( NCELL,   21) =  RJBLK( NCELL, IJ_HKET_RACM2 )
!  Reaction Label R019            
                RKI( NCELL,   22) =  RJBLK( NCELL, IJ_MACR_RACM2 )
!  Reaction Label R020            
                RKI( NCELL,   23) =  RJBLK( NCELL, IJ_MVK_JPL19 )
!  Reaction Label R021            
                RKI( NCELL,   24) =  RJBLK( NCELL, IJ_GLYH2_JPL19 )
!  Reaction Label R022            
                RKI( NCELL,   25) =  RJBLK( NCELL, IJ_GLYF_JPL19 )
!  Reaction Label R023            
                RKI( NCELL,   26) =  RJBLK( NCELL, IJ_GLYHX_JPL19 )
!  Reaction Label R024            
                RKI( NCELL,   27) =  RJBLK( NCELL, IJ_MGLY_RACM2 )
!  Reaction Label R025            
                RKI( NCELL,   28) =  RJBLK( NCELL, IJ_MGLY_RACM2 )
!  Reaction Label R026            
                RKI( NCELL,   29) =  RJBLK( NCELL, IJ_MGLY_RACM2 )
!  Reaction Label R027a           
                RKI( NCELL,   30) =  RJBLK( NCELL, IJ_BALD1_CALVERT11 )
!  Reaction Label R027b           
                RKI( NCELL,   31) =  RJBLK( NCELL, IJ_BALD2_CALVERT11 )
!  Reaction Label R028            
                RKI( NCELL,   32) =  RJBLK( NCELL, IJ_OP1_RACM2 )
!  Reaction Label R029            
                RKI( NCELL,   33) =  RJBLK( NCELL, IJ_OP1_RACM2 )
!  Reaction Label TRP03           
                RKI( NCELL,   34) =  RJBLK( NCELL, IJ_OP1_RACM2 )
!  Reaction Label R029a           
                RKI( NCELL,   35) =  RJBLK( NCELL, IJ_OP1_RACM2 )
!  Reaction Label R030            
                RKI( NCELL,   36) =  RJBLK( NCELL, IJ_PAA_RACM2 )
!  Reaction Label R031            
                RKI( NCELL,   37) =  RJBLK( NCELL, IJ_ONIT_CALVERT08 )
!  Reaction Label R032            
                RKI( NCELL,   38) =  RJBLK( NCELL, IJ_PAN1_JPL19 )
!  Reaction Label R033            
                RKI( NCELL,   39) =  RJBLK( NCELL, IJ_PAN2_JPL19 )
!  Reaction Label R033a           
                RKI( NCELL,   40) =  RJBLK( NCELL, IJ_PPN1_JPL19 )
!  Reaction Label R033b           
                RKI( NCELL,   41) =  RJBLK( NCELL, IJ_PPN2_JPL19 )
!  Reaction Label TRP55           
                RKI( NCELL,   42) =  RJBLK( NCELL, IJ_TRPN_WANG2023 )
!  Reaction Label TRP56           
                RKI( NCELL,   43) =  RJBLK( NCELL, IJ_TRPN_WANG2023 )
!  Reaction Label T20             
                RKI( NCELL,  404) =  RJBLK( NCELL, IJ_ACRO_09 )
!  Reaction Label SOAphot01       
                RKI( NCELL,  542) =   1.0000D-02 * RJBLK( NCELL, IJ_NO2_RACM2 )
!  Reaction Label CLP01           
                RKI( NCELL,  543) =  RJBLK( NCELL, IJ_CL2_JPL19 )
!  Reaction Label CLP02           
                RKI( NCELL,  544) =  RJBLK( NCELL, IJ_CLO_JPL19 )
!  Reaction Label CLP03           
                RKI( NCELL,  545) =  RJBLK( NCELL, IJ_OCLO_JPL19 )
!  Reaction Label CLP07           
                RKI( NCELL,  546) =  RJBLK( NCELL, IJ_CL2O2_JPL19 )
!  Reaction Label CLP10           
                RKI( NCELL,  547) =  RJBLK( NCELL, IJ_HOCL_JPL19 )
!  Reaction Label CLP12           
                RKI( NCELL,  548) =  RJBLK( NCELL, IJ_CLNO_JPL19 )
!  Reaction Label CLP14           
                RKI( NCELL,  549) =  RJBLK( NCELL, IJ_CLNO2_JPL19 )
!  Reaction Label CLP17           
                RKI( NCELL,  550) =  RJBLK( NCELL, IJ_CLNO3_R_JPL19 )
!  Reaction Label CLP18           
                RKI( NCELL,  551) =  RJBLK( NCELL, IJ_CLNO3_M_JPL19 )
!  Reaction Label CLP19           
                RKI( NCELL,  552) =  RJBLK( NCELL, IJ_HCOCL_JPL19 )
!  Reaction Label BRP01           
                RKI( NCELL,  740) =  RJBLK( NCELL, IJ_BR2_JPL19 )
!  Reaction Label BRP02           
                RKI( NCELL,  741) =  RJBLK( NCELL, IJ_BRO_JPL19 )
!  Reaction Label BRP03           
                RKI( NCELL,  742) =  RJBLK( NCELL, IJ_OBRO_JPL19 )
!  Reaction Label BRP04           
                RKI( NCELL,  743) =  RJBLK( NCELL, IJ_HOBR_JPL19 )
!  Reaction Label BRP06           
                RKI( NCELL,  744) =  RJBLK( NCELL, IJ_BRNO_JPL19 )
!  Reaction Label BRP08           
                RKI( NCELL,  745) =  RJBLK( NCELL, IJ_BRNO2_JPL19 )
!  Reaction Label BRP12           
                RKI( NCELL,  746) =  RJBLK( NCELL, IJ_BRNO3_R_JPL19 )
!  Reaction Label BRP13           
                RKI( NCELL,  747) =  RJBLK( NCELL, IJ_BRNO3_M_JPL19 )
!  Reaction Label BRP15           
                RKI( NCELL,  748) =  RJBLK( NCELL, IJ_CH2BR2_JPL19 )
!  Reaction Label BRP16           
                RKI( NCELL,  749) =  RJBLK( NCELL, IJ_CHBR3_JPL19 )
!  Reaction Label BRP17           
                RKI( NCELL,  750) =  RJBLK( NCELL, IJ_HCOBR_JPL19 )
!  Reaction Label IP01            
                RKI( NCELL,  811) =  RJBLK( NCELL, IJ_I2_JPL19 )
!  Reaction Label IP02            
                RKI( NCELL,  812) =  RJBLK( NCELL, IJ_IO_JPL19 )
!  Reaction Label IP03            
                RKI( NCELL,  813) =  RJBLK( NCELL, IJ_OIO_JPL19 )
!  Reaction Label IP05            
                RKI( NCELL,  814) =  RJBLK( NCELL, IJ_INO3_06 )
!  Reaction Label IP09            
                RKI( NCELL,  815) =  RJBLK( NCELL, IJ_HOI_JPL19 )
!  Reaction Label IP10            
                RKI( NCELL,  816) =  RJBLK( NCELL, IJ_HI_JPL19 )
!  Reaction Label IP11            
                RKI( NCELL,  817) =  RJBLK( NCELL, IJ_INO_JPL19 )
!  Reaction Label IP12            
                RKI( NCELL,  818) =  RJBLK( NCELL, IJ_INO2_JPL19 )
!  Reaction Label IP16            
                RKI( NCELL,  819) =  RJBLK( NCELL, IJ_INO3_06 )
!  Reaction Label IP19            
                RKI( NCELL,  820) =  RJBLK( NCELL, IJ_CH3I_JPL19 )
!  Reaction Label IP20            
                RKI( NCELL,  821) =  RJBLK( NCELL, IJ_CH2I2_JPL19 )
!  Reaction Label IO21            
                RKI( NCELL,  822) =  RJBLK( NCELL, IJ_INO3_06 )
!  Reaction Label IO22            
                RKI( NCELL,  823) =  RJBLK( NCELL, IJ_INO3_06 )
!  Reaction Label MP01            
                RKI( NCELL,  851) =  RJBLK( NCELL, IJ_BRCL_JPL19 )
!  Reaction Label MP03            
                RKI( NCELL,  852) =  RJBLK( NCELL, IJ_ICL_JPL19 )
!  Reaction Label MP04            
                RKI( NCELL,  853) =  RJBLK( NCELL, IJ_IBR_JPL19 )
!  Reaction Label MP06            
                RKI( NCELL,  854) =  RJBLK( NCELL, IJ_CH2IBR_JPL19 )
!  Reaction Label MP07            
                RKI( NCELL,  855) =  RJBLK( NCELL, IJ_CH2ICL_JPL19 )
!  Reaction Label MP08            
                RKI( NCELL,  856) =  RJBLK( NCELL, IJ_CHBR2CL_JPL19 )
!  Reaction Label MP09            
                RKI( NCELL,  857) =  RJBLK( NCELL, IJ_CHBRCL2_JPL19 )
            END DO 
       END IF 

        DO NCELL = 1, NUMCELLS 
!  Set-up conversion factors 
             INV_TEMP  = 1.0D+00 / BLKTEMP( NCELL ) 
             CAIR      = 1.0D+06 * COEF1 * BLKPRES( NCELL ) * INV_TEMP 
             CFACT     = 6.0D-05 * CAIR
             CFACT_SQU = 6.0D-11 * CAIR * CAIR 
             INV_CAIR  = 1.0D0 / CAIR 
             INV_CFACT = 6.0D+07 * INV_CAIR 
             TEMP      = BLKTEMP( NCELL ) 
             TEMPOT300 = BLKTEMP( NCELL ) * TI300 
             RFACT     = 1.0D+06 * INV_CAIR 
             RFACT_SQU = 1.0D+12 * INV_CAIR * INV_CAIR 

!  Reaction Label R034            
             RKI( NCELL,   44) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7000D-12,  -9.4000D+02 )
!  Reaction Label R035            
             RKI( NCELL,   45) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.0000D-14,  -4.9000D+02 )
!  Reaction Label R036            
             RKI( NCELL,   46) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.0000D-12,  -1.5000D+03 )
!  Reaction Label R037            
             RKI( NCELL,   47) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.2000D-13,  -2.4500D+03 )
!  Reaction Label R038            
             RKI( NCELL,   48) =  CFACT_SQU * POWER_T02( TEMPOT300,   6.1000D-34,  -2.4000D+00 )
!  Reaction Label R039            
             RKI( NCELL,   49) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.0000D-12,  -2.0600D+03 )
!  Reaction Label R040            
             RKI( NCELL,   50) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.3000D-11,   5.5000D+01 )
!  Reaction Label R041            
             RKI( NCELL,   51) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.1500D-11,   1.1000D+02 )
!  Reaction Label R042            
             RKI( NCELL,   52) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6300D-10,   6.0000D+01 )
!  Reaction Label R043            
             RKI( NCELL,   53) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-12,  -1.8000D+03 )
!  Reaction Label R044            
             RKI( NCELL,   54) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.8000D-11,   2.5000D+02 )
!  Reaction Label R045            
             RKI( NCELL,   55) =  CFACT * FALLOFF_T09( INV_TEMP,  CAIR, & 
     &                                                 3.0000D-13,   4.6000D+02,   2.1000D-33,  & 
     &                                                 9.2000D+02 )
!  Reaction Label R046            
             RKI( NCELL,   56) =  CFACT_SQU * FALLOFF_T09( INV_TEMP,  CAIR, & 
     &                                                 4.2000D-34,   2.6600D+03,   2.9400D-54,  & 
     &                                                 3.1200D+03 )
!  Reaction Label R047            
             RKI( NCELL,   57) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.8000D-12,   0.0000D+00 )
!  Reaction Label R048            
             RKI( NCELL,   58) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 9.1000D-32,   0.0000D+00,  -1.5000D+00,  & 
     &                                                 3.0000D-11,   0.0000D+00,   0.0000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R049            
             RKI( NCELL,   59) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 7.1000D-31,   0.0000D+00,  -2.6000D+00,  & 
     &                                                 3.6000D-11,   0.0000D+00,  -1.0000D-01,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R050            
             RKI( NCELL,   60) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4400D-12,   2.6000D+02 )
!  Reaction Label R051            
             RKI( NCELL,   61) =  CFACT * FALLOFF_T11( INV_TEMP,TEMPOT300,CAIR, & 
     &                                                 6.0950D-14,  -1.0000D+00,   2.7000D+02, &
     &                                                 6.8570D-34,   1.0000D+00,   2.7000D+02,  & 
     &                                                -5.9680D-14,   2.7000D+02 )
!  Reaction Label R052            
             RKI( NCELL,   62) =  CFACT_SQU * ARRHENUIS_T03( INV_TEMP,  4.2500D-39,   6.6350D+02 )
!  Reaction Label R053            
             RKI( NCELL,   63) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.0000D-12,   2.5000D+02 )
!  Reaction Label R054            
             RKI( NCELL,   64) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.3000D-12,   2.0000D+02 )
!  Reaction Label R055            
             RKI( NCELL,   65) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 3.4000D-31,   0.0000D+00,  -1.6000D+00,  & 
     &                                                 2.3000D-11,   0.0000D+00,  -2.0000D-01,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R056            
             RKI( NCELL,   66) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 1.8000D-30,   0.0000D+00,  -3.0000D+00,  & 
     &                                                 2.8000D-11,   0.0000D+00,   0.0000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R057            
             RKI( NCELL,   67) =  CFACT * FALLOFF_T08( INV_TEMP,  CAIR, & 
     &                                                 2.4000D-14,   4.6000D+02,   2.7000D-17,  & 
     &                                                 2.1990D+03,   6.5000D-34,   1.3350D+03 )
!  Reaction Label R058            
             RKI( NCELL,   68) =   2.0000D-11 * CFACT 
!  Reaction Label R059            
             RKI( NCELL,   69) =   3.5000D-12 * CFACT 
!  Reaction Label R060            
             RKI( NCELL,   70) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7000D-11,   1.2500D+02 )
!  Reaction Label R061            
             RKI( NCELL,   71) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.3500D-14,  -1.3350D+03 )
!  Reaction Label R062            
             RKI( NCELL,   72) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.5000D-13,  -2.4500D+03 )
!  Reaction Label R063            
             RKI( NCELL,   73) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 2.4000D-30,   0.0000D+00,  -3.0000D+00,  & 
     &                                                 1.6000D-12,   0.0000D+00,   1.0000D-01,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R064            
             RKI( NCELL,   74) =  RFACT * RKI( NCELL,   73 ) & 
     &                         * (  1.7241D+26 * DEXP( -1.0840D+04 * INV_TEMP) ) 
!  Reaction Label R065            
             RKI( NCELL,   75) =   1.0000D-22 * CFACT 
!  Reaction Label R066            
             RKI( NCELL,   76) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 1.9000D-31,   0.0000D+00,  -3.4000D+00,  & 
     &                                                 4.0000D-12,   0.0000D+00,  -3.0000D-01,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R067            
             RKI( NCELL,   77) =  RFACT * RKI( NCELL,   76 ) & 
     &                         * (  4.7619D+26 * DEXP( -1.0900D+04 * INV_TEMP) ) 
!  Reaction Label R068            
             RKI( NCELL,   78) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.5000D-13,   6.1000D+02 )
!  Reaction Label R069            
             RKI( NCELL,   79) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 2.9000D-31,   0.0000D+00,  -4.1000D+00,  & 
     &                                                 1.7000D-12,   0.0000D+00,   2.0000D-01,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R070            
             RKI( NCELL,   80) =  CFACT * FALLOFF_T09( INV_TEMP,  CAIR, & 
     &                                                 1.4400D-13,   0.0000D+00,   2.7400D-33,  & 
     &                                                 0.0000D+00 )
!  Reaction Label R071            
             RKI( NCELL,   81) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.4500D-12,  -1.7750D+03 )
!  Reaction Label R072            
             RKI( NCELL,   82) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.6600D-12,  -1.0200D+03 )
!  Reaction Label R073            
             RKI( NCELL,   83) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.6800D-12,  -3.7000D+02 )
!  Reaction Label R074            
             RKI( NCELL,   84) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.0100D-11,  -2.4500D+02 )
!  Reaction Label R076            
             RKI( NCELL,   85) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 1.0000D-28,   0.0000D+00,  -4.5000D+00,  & 
     &                                                 8.8000D-12,   0.0000D+00,  -8.5000D-01,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R077            
             RKI( NCELL,   86) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.7200D-12,   5.0000D+02 )
!  Reaction Label R078            
             RKI( NCELL,   87) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.3300D-11,   5.0000D+02 )
!  Reaction Label R080            
             RKI( NCELL,   88) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 5.5000D-30,   0.0000D+00,   0.0000D+00,  & 
     &                                                 8.3000D-13,   0.0000D+00,   2.0000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label ROCARO31        
             RKI( NCELL,   89) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.3300D-12,  -1.9300D+02 )
!  Reaction Label ROCARO41        
             RKI( NCELL,   90) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.8100D-12,   3.5400D+02 )
!  Reaction Label ROCARO51        
             RKI( NCELL,   91) =   2.3300D-11 * CFACT 
!  Reaction Label ROCARO61        
             RKI( NCELL,   92) =   7.0000D-12 * CFACT 
!  Reaction Label ROCARO61a       
             RKI( NCELL,   93) =   5.7000D-16 * CFACT 
!  Reaction Label RAM01           
             RKI( NCELL,   94) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.5800D-14,  -2.0000D+03 )
!  Reaction Label RAM02           
             RKI( NCELL,   95) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9500D-12,  -4.5000D+02 )
!  Reaction Label RAM03           
             RKI( NCELL,   96) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6900D-11,   3.9000D+02 )
!  Reaction Label RAM04           
             RKI( NCELL,   97) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.5000D-13,   1.3000D+03 )
!  Reaction Label RAM05           
             RKI( NCELL,   98) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.0000D-12,   3.5000D+02 )
!  Reaction Label RAM06           
             RKI( NCELL,   99) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.6000D-12,   2.0000D+02 )
!  Reaction Label RAM07           
             RKI( NCELL,  100) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9700D-11,   3.9000D+02 )
!  Reaction Label RAM08           
             RKI( NCELL,  101) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.1400D-11,   5.8000D+02 )
!  Reaction Label RAM09           
             RKI( NCELL,  102) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.4200D-12,   5.8000D+02 )
!  Reaction Label RAM10           
             RKI( NCELL,  103) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.4000D-11,   3.9000D+02 )
!  Reaction Label RAM11           
             RKI( NCELL,  104) =   1.5000D-11 * CFACT 
!  Reaction Label RAM12           
             RKI( NCELL,  105) =   4.0000D-05 * SFACT 
!  Reaction Label RAM13           
             RKI( NCELL,  106) =   3.0000D-12 * CFACT 
!  Reaction Label R087            
             RKI( NCELL,  107) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.2100D-11,   4.4000D+02 )
!  Reaction Label R088            
             RKI( NCELL,  108) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.2000D-11,   4.0100D+02 )
!  Reaction Label TRP04           
             RKI( NCELL,  109) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.2000D-12,   6.0000D+02 )
!  Reaction Label TRP05           
             RKI( NCELL,  110) =   1.1000D-10 * CFACT 
!  Reaction Label R089            
             RKI( NCELL,  111) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.5000D-12,   1.2500D+02 )
!  Reaction Label R090            
             RKI( NCELL,  112) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.7000D-12,   3.4500D+02 )
!  Reaction Label R091            
             RKI( NCELL,  113) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.9000D-12,   4.0500D+02 )
!  Reaction Label R092            
             RKI( NCELL,  114) =  CFACT * ARRHENUIS_T04( INV_TEMP,  TEMPOT300, & 
     &                                                   4.5600D-14,  -4.2700D+02,   3.6500D+00 )
!  Reaction Label R093            
             RKI( NCELL,  115) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.5000D-12,  -9.0000D+01 )
!  Reaction Label R094            
             RKI( NCELL,  116) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-12,   1.0000D+01 )
!  Reaction Label R095            
             RKI( NCELL,  117) =   3.0000D-12 * CFACT 
!  Reaction Label R096            
             RKI( NCELL,  118) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.0000D-12,   3.8000D+02 )
!  Reaction Label R097            
             RKI( NCELL,  119) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-12,   6.1000D+02 )
!  Reaction Label R098            
             RKI( NCELL,  120) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.7700D-12,   5.3300D+02 )
!  Reaction Label R099            
             RKI( NCELL,  121) =   1.1000D-11 * CFACT 
!  Reaction Label R100            
             RKI( NCELL,  122) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.2600D-13,   8.3000D+02 )
!  Reaction Label R101            
             RKI( NCELL,  123) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-11,   1.7500D+02 )
!  Reaction Label R102            
             RKI( NCELL,  124) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-11,   1.7500D+02 )
!  Reaction Label R103            
             RKI( NCELL,  125) =   1.0000D-11 * CFACT 
!  Reaction Label R104            
             RKI( NCELL,  126) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.3200D-12,   2.4300D+02 )
!  Reaction Label R105            
             RKI( NCELL,  127) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.7500D-12,   4.0500D+02 )
!  Reaction Label R106            
             RKI( NCELL,  128) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.6500D-11,   0.0000D+00 )
!  Reaction Label R108            
             RKI( NCELL,  129) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0500D-10,   0.0000D+00 )
!  Reaction Label R109            
             RKI( NCELL,  130) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8500D-12,  -3.4500D+02 )
!  Reaction Label R110            
             RKI( NCELL,  131) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.0000D-12,   2.0000D+01 )
!  Reaction Label R111            
             RKI( NCELL,  132) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-12,   2.0000D+02 )
!  Reaction Label R112            
             RKI( NCELL,  133) =   1.4700D-11 * CFACT 
!  Reaction Label R113            
             RKI( NCELL,  134) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9000D-12,   1.9000D+02 )
!  Reaction Label R114            
             RKI( NCELL,  135) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-12,   1.9000D+02 )
!  Reaction Label TRP06           
             RKI( NCELL,  136) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-12,   1.9000D+02 )
!  Reaction Label R116            
             RKI( NCELL,  137) =   3.0000D-11 * CFACT 
!  Reaction Label R117            
             RKI( NCELL,  138) =   4.5000D-13 * CFACT 
!  Reaction Label R118            
             RKI( NCELL,  139) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.0000D-14,   8.5000D+02 )
!  Reaction Label R119            
             RKI( NCELL,  140) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9300D-12,   1.9000D+02 )
!  Reaction Label R120            
             RKI( NCELL,  141) =   4.0000D-14 * CFACT 
!  Reaction Label R121            
             RKI( NCELL,  142) =   4.0000D-14 * CFACT 
!  Reaction Label R122            
             RKI( NCELL,  143) =   3.2000D-11 * CFACT 
!  Reaction Label R123            
             RKI( NCELL,  144) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.3100D-12,  -2.6000D+02 )
!  Reaction Label TRP07           
             RKI( NCELL,  145) =   4.8000D-12 * CFACT 
!  Reaction Label TRP57           
             RKI( NCELL,  146) =  CFACT * FALLOFF_T08( INV_TEMP,  CAIR, & 
     &                                                 2.4000D-14,   4.6000D+02,   2.7000D-17,  & 
     &                                                 2.1990D+03,   6.5000D-34,   1.3350D+03 )
!  Reaction Label R126            
             RKI( NCELL,  147) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.1400D-15,  -2.5800D+03 )
!  Reaction Label R127            
             RKI( NCELL,  148) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.3300D-15,  -1.8000D+03 )
!  Reaction Label R128            
             RKI( NCELL,  149) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.4000D-15,  -8.4500D+02 )
!  Reaction Label R131            
             RKI( NCELL,  150) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.0500D-16,  -6.4000D+02 )
!  Reaction Label R132            
             RKI( NCELL,  151) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-15,  -7.7000D+02 )
!  Reaction Label TRP08           
             RKI( NCELL,  152) =   8.3000D-18 * CFACT 
!  Reaction Label TRP09           
             RKI( NCELL,  153) =   1.6700D-16 * CFACT 
!  Reaction Label R133            
             RKI( NCELL,  154) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.3600D-15,  -2.1120D+03 )
!  Reaction Label R134            
             RKI( NCELL,  155) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.5000D-16,  -1.5200D+03 )
!  Reaction Label R135            
             RKI( NCELL,  156) =   1.6600D-18 * CFACT 
!  Reaction Label R136            
             RKI( NCELL,  157) =   2.0000D-16 * CFACT 
!  Reaction Label R137            
             RKI( NCELL,  158) =   2.0000D-16 * CFACT 
!  Reaction Label R138            
             RKI( NCELL,  159) =   9.0000D-17 * CFACT 
!  Reaction Label R140            
             RKI( NCELL,  160) =   2.8600D-13 * CFACT 
!  Reaction Label R141            
             RKI( NCELL,  161) =  CFACT * ARRHENUIS_T04( INV_TEMP,  TEMPOT300, & 
     &                                                   4.3920D-13,  -2.2820D+03,   2.0000D+00 )
!  Reaction Label R142            
             RKI( NCELL,  162) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7900D-13,  -4.5000D+02 )
!  Reaction Label R143            
             RKI( NCELL,  163) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.6400D-13,   4.5000D+02 )
!  Reaction Label R146            
             RKI( NCELL,  164) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1900D-12,   4.9000D+02 )
!  Reaction Label R147            
             RKI( NCELL,  165) =   1.2200D-11 * CFACT 
!  Reaction Label TRP10           
             RKI( NCELL,  166) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.1500D-13,  -4.4800D+02 )
!  Reaction Label R148            
             RKI( NCELL,  167) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-12,  -2.4400D+03 )
!  Reaction Label R149            
             RKI( NCELL,  168) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.4000D-12,  -1.9000D+03 )
!  Reaction Label R150            
             RKI( NCELL,  169) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.7600D-12,  -1.9000D+03 )
!  Reaction Label R151            
             RKI( NCELL,  170) =   3.4000D-15 * CFACT 
!  Reaction Label R152            
             RKI( NCELL,  171) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.0200D-13,  -1.0760D+03 )
!  Reaction Label R153            
             RKI( NCELL,  172) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9000D-12,  -1.9000D+03 )
!  Reaction Label R154            
             RKI( NCELL,  173) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.7600D-12,  -1.9000D+03 )
!  Reaction Label RBAL1           
             RKI( NCELL,  174) =   2.4000D-15 * CFACT 
!  Reaction Label R155            
             RKI( NCELL,  175) =   3.7800D-12 * CFACT 
!  Reaction Label R156            
             RKI( NCELL,  176) =   1.0600D-12 * CFACT 
!  Reaction Label R158            
             RKI( NCELL,  177) =   2.0100D-10 * CFACT 
!  Reaction Label R159            
             RKI( NCELL,  178) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.2000D-14,  -5.0000D+02 )
!  Reaction Label RBAL2           
             RKI( NCELL,  179) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 9.7000D-29,   0.0000D+00,  -5.6000D+00,  & 
     &                                                 9.3000D-12,   0.0000D+00,  -1.5000D+00,  & 
     &                                                 3.0000D-01,   1.4200D+00 )
!  Reaction Label TRP11           
             RKI( NCELL,  180) =   2.9000D-02 * SFACT 
!  Reaction Label TRP12           
             RKI( NCELL,  181) =   2.4000D-02 * SFACT 
!  Reaction Label R166            
             RKI( NCELL,  182) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 9.7000D-29,   0.0000D+00,  -5.6000D+00,  & 
     &                                                 9.3000D-12,   0.0000D+00,  -1.5000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R167            
             RKI( NCELL,  183) =  RFACT * RKI( NCELL,  182 ) & 
     &                         * (  1.1111D+28 * DEXP( -1.4000D+04 * INV_TEMP) ) 
!  Reaction Label R168            
             RKI( NCELL,  184) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 9.7000D-29,   0.0000D+00,  -5.6000D+00,  & 
     &                                                 9.3000D-12,   0.0000D+00,  -1.5000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label R169            
             RKI( NCELL,  185) =  RFACT * RKI( NCELL,  184 ) & 
     &                         * (  1.1111D+28 * DEXP( -1.4000D+04 * INV_TEMP) ) 
!  Reaction Label R170            
             RKI( NCELL,  186) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-12,   1.8100D+02 )
!  Reaction Label R171            
             RKI( NCELL,  187) =  SFACT * ARRHENUIS_T03( INV_TEMP,  1.6000D+16,  -1.3486D+04 )
!  Reaction Label R172            
             RKI( NCELL,  188) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-12,   3.0000D+02 )
!  Reaction Label R173            
             RKI( NCELL,  189) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-12,   3.6500D+02 )
!  Reaction Label R174            
             RKI( NCELL,  190) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R175            
             RKI( NCELL,  191) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R177            
             RKI( NCELL,  192) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R178            
             RKI( NCELL,  193) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R179            
             RKI( NCELL,  194) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO33        
             RKI( NCELL,  195) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO43        
             RKI( NCELL,  196) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO53        
             RKI( NCELL,  197) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO63        
             RKI( NCELL,  198) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R189            
             RKI( NCELL,  199) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP13           
             RKI( NCELL,  200) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP14           
             RKI( NCELL,  201) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP15           
             RKI( NCELL,  202) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R190            
             RKI( NCELL,  203) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP16           
             RKI( NCELL,  204) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP17           
             RKI( NCELL,  205) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP18           
             RKI( NCELL,  206) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP19           
             RKI( NCELL,  207) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label TRP20           
             RKI( NCELL,  208) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R191            
             RKI( NCELL,  209) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.1000D-12,   2.7000D+02 )
!  Reaction Label R192            
             RKI( NCELL,  210) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.1000D-12,   2.7000D+02 )
!  Reaction Label R193            
             RKI( NCELL,  211) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9000D-12,   3.0000D+02 )
!  Reaction Label R194            
             RKI( NCELL,  212) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R195            
             RKI( NCELL,  213) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R196            
             RKI( NCELL,  214) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5400D-12,   3.6000D+02 )
!  Reaction Label R197            
             RKI( NCELL,  215) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5400D-12,   3.6000D+02 )
!  Reaction Label R198            
             RKI( NCELL,  216) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5400D-12,   3.6000D+02 )
!  Reaction Label R199            
             RKI( NCELL,  217) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5400D-12,   3.6000D+02 )
!  Reaction Label R200            
             RKI( NCELL,  218) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R201            
             RKI( NCELL,  219) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R202            
             RKI( NCELL,  220) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R203            
             RKI( NCELL,  221) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R204            
             RKI( NCELL,  222) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R205            
             RKI( NCELL,  223) =   4.0000D-12 * CFACT 
!  Reaction Label R206            
             RKI( NCELL,  224) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R207            
             RKI( NCELL,  225) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R208            
             RKI( NCELL,  226) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R209            
             RKI( NCELL,  227) =   2.0000D-11 * CFACT 
!  Reaction Label R210            
             RKI( NCELL,  228) =   2.0000D-11 * CFACT 
!  Reaction Label R211            
             RKI( NCELL,  229) =   2.0800D-12 * CFACT 
!  Reaction Label R212            
             RKI( NCELL,  230) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.1000D-13,   7.5000D+02 )
!  Reaction Label R213            
             RKI( NCELL,  231) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.5000D-13,   7.0000D+02 )
!  Reaction Label R214            
             RKI( NCELL,  232) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6600D-13,   1.3000D+03 )
!  Reaction Label R215            
             RKI( NCELL,  233) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6600D-13,   1.3000D+03 )
!  Reaction Label R217            
             RKI( NCELL,  234) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.9000D-13,   1.3000D+03 )
!  Reaction Label R218            
             RKI( NCELL,  235) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6600D-13,   1.3000D+03 )
!  Reaction Label R219            
             RKI( NCELL,  236) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6600D-13,   1.3000D+03 )
!  Reaction Label ROCARO32        
             RKI( NCELL,  237) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label ROCARO42        
             RKI( NCELL,  238) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label ROCARO52        
             RKI( NCELL,  239) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label ROCARO62        
             RKI( NCELL,  240) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label R229            
             RKI( NCELL,  241) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-13,   1.3000D+03 )
!  Reaction Label TRP21           
             RKI( NCELL,  242) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7300D-13,   1.3000D+03 )
!  Reaction Label TRP22           
             RKI( NCELL,  243) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7100D-13,   1.3000D+03 )
!  Reaction Label TRP23           
             RKI( NCELL,  244) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7500D-13,   1.3000D+03 )
!  Reaction Label R230            
             RKI( NCELL,  245) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-13,   1.3000D+03 )
!  Reaction Label TRP24           
             RKI( NCELL,  246) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7300D-13,   1.3000D+03 )
!  Reaction Label TRP25           
             RKI( NCELL,  247) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7100D-13,   1.3000D+03 )
!  Reaction Label TRP26           
             RKI( NCELL,  248) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7500D-13,   1.3000D+03 )
!  Reaction Label TRP27           
             RKI( NCELL,  249) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7100D-13,   1.3000D+03 )
!  Reaction Label TRP28           
             RKI( NCELL,  250) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7300D-13,   1.3000D+03 )
!  Reaction Label R231            
             RKI( NCELL,  251) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.3000D-13,   1.0400D+03 )
!  Reaction Label R232            
             RKI( NCELL,  252) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.3000D-13,   1.0400D+03 )
!  Reaction Label R233            
             RKI( NCELL,  253) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1500D-13,   1.3000D+03 )
!  Reaction Label R234            
             RKI( NCELL,  254) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1500D-13,   1.3000D+03 )
!  Reaction Label R235            
             RKI( NCELL,  255) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1500D-13,   1.3000D+03 )
!  Reaction Label R236            
             RKI( NCELL,  256) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.8200D-13,   1.3000D+03 )
!  Reaction Label R237            
             RKI( NCELL,  257) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.8200D-13,   1.3000D+03 )
!  Reaction Label R238            
             RKI( NCELL,  258) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label R239            
             RKI( NCELL,  259) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label RBAL3           
             RKI( NCELL,  260) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.3800D-13,   1.3000D+03 )
!  Reaction Label RBAL4           
             RKI( NCELL,  261) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.1100D-13,   1.3000D+03 )
!  Reaction Label RBAL5           
             RKI( NCELL,  262) =   2.8600D-13 * CFACT 
!  Reaction Label R240            
             RKI( NCELL,  263) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.7500D-13,   9.8000D+02 )
!  Reaction Label R241            
             RKI( NCELL,  264) =   1.0000D-11 * CFACT 
!  Reaction Label R242            
             RKI( NCELL,  265) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.7500D-13,   9.8000D+02 )
!  Reaction Label R243            
             RKI( NCELL,  266) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1500D-13,   1.3000D+03 )
!  Reaction Label R244            
             RKI( NCELL,  267) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6600D-13,   1.3000D+03 )
!  Reaction Label R245            
             RKI( NCELL,  268) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6600D-13,   1.3000D+03 )
!  Reaction Label R246            
             RKI( NCELL,  269) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.7500D-13,   9.8000D+02 )
!  Reaction Label R247            
             RKI( NCELL,  270) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6600D-13,   1.3000D+03 )
!  Reaction Label R248            
             RKI( NCELL,  271) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.5000D-14,   3.9000D+02 )
!  Reaction Label R249            
             RKI( NCELL,  272) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1800D-13,   1.5800D+02 )
!  Reaction Label R250            
             RKI( NCELL,  273) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.4600D-14,   4.3100D+02 )
!  Reaction Label R251            
             RKI( NCELL,  274) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.0000D-13,   4.6700D+02 )
!  Reaction Label R253            
             RKI( NCELL,  275) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7100D-13,   7.0800D+02 )
!  Reaction Label R254            
             RKI( NCELL,  276) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.4600D-13,   7.0800D+02 )
!  Reaction Label R255            
             RKI( NCELL,  277) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.1800D-14,   7.0800D+02 )
!  Reaction Label ROCARO35        
             RKI( NCELL,  278) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label ROCARO45        
             RKI( NCELL,  279) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label ROCARO55        
             RKI( NCELL,  280) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label ROCARO65        
             RKI( NCELL,  281) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label R264            
             RKI( NCELL,  282) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-14,   2.2100D+02 )
!  Reaction Label R265            
             RKI( NCELL,  283) =   2.0000D-12 * CFACT 
!  Reaction Label TRP29           
             RKI( NCELL,  284) =   1.0000D-10 * CFACT 
!  Reaction Label TRP30           
             RKI( NCELL,  285) =   2.0000D-12 * CFACT 
!  Reaction Label TRP31           
             RKI( NCELL,  286) =   1.0000D-10 * CFACT 
!  Reaction Label R266            
             RKI( NCELL,  287) =   2.0000D-12 * CFACT 
!  Reaction Label TRP32           
             RKI( NCELL,  288) =   1.0000D-10 * CFACT 
!  Reaction Label TRP33           
             RKI( NCELL,  289) =   2.0000D-12 * CFACT 
!  Reaction Label TRP34           
             RKI( NCELL,  290) =   1.0000D-10 * CFACT 
!  Reaction Label R267            
             RKI( NCELL,  291) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-11,   5.0000D+02 )
!  Reaction Label R268            
             RKI( NCELL,  292) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-11,   5.0000D+02 )
!  Reaction Label R269            
             RKI( NCELL,  293) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.5000D-13,   5.0000D+02 )
!  Reaction Label R270            
             RKI( NCELL,  294) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.9100D-13,   5.0800D+02 )
!  Reaction Label R271            
             RKI( NCELL,  295) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.9100D-13,   5.0800D+02 )
!  Reaction Label R272            
             RKI( NCELL,  296) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-14,   2.2100D+02 )
!  Reaction Label R273            
             RKI( NCELL,  297) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-14,   2.2100D+02 )
!  Reaction Label R274            
             RKI( NCELL,  298) =   8.3700D-14 * CFACT 
!  Reaction Label R275            
             RKI( NCELL,  299) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-14,   2.2100D+02 )
!  Reaction Label R276            
             RKI( NCELL,  300) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label R277            
             RKI( NCELL,  301) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label R278            
             RKI( NCELL,  302) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label R279            
             RKI( NCELL,  303) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label R280            
             RKI( NCELL,  304) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.5000D-13,   5.0000D+02 )
!  Reaction Label R281            
             RKI( NCELL,  305) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6000D-13,   7.0800D+02 )
!  Reaction Label R282            
             RKI( NCELL,  306) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.6800D-14,   7.0800D+02 )
!  Reaction Label R283            
             RKI( NCELL,  307) =   3.5600D-14 * CFACT 
!  Reaction Label R284            
             RKI( NCELL,  308) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.9900D-15,   1.5100D+03 )
!  Reaction Label R285            
             RKI( NCELL,  309) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.0300D-12,   2.1100D+02 )
!  Reaction Label R286            
             RKI( NCELL,  310) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.9000D-13,   4.6000D+02 )
!  Reaction Label R287            
             RKI( NCELL,  311) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.5900D-13,   5.2200D+02 )
!  Reaction Label R289            
             RKI( NCELL,  312) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.4800D-13,   7.6500D+02 )
!  Reaction Label R290            
             RKI( NCELL,  313) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.1100D-13,   7.6500D+02 )
!  Reaction Label R291            
             RKI( NCELL,  314) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.0900D-13,   7.6500D+02 )
!  Reaction Label ROCARO36        
             RKI( NCELL,  315) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label ROCARO46        
             RKI( NCELL,  316) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label ROCARO56        
             RKI( NCELL,  317) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label ROCARO66        
             RKI( NCELL,  318) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label R300            
             RKI( NCELL,  319) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.4000D-14,   2.2100D+02 )
!  Reaction Label R301            
             RKI( NCELL,  320) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-12,   5.0000D+02 )
!  Reaction Label TRP35           
             RKI( NCELL,  321) =   1.0000D-10 * CFACT 
!  Reaction Label TRP36           
             RKI( NCELL,  322) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-12,   5.0000D+02 )
!  Reaction Label TRP37           
             RKI( NCELL,  323) =   1.0000D-10 * CFACT 
!  Reaction Label R302            
             RKI( NCELL,  324) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-12,   5.0000D+02 )
!  Reaction Label TRP38           
             RKI( NCELL,  325) =   1.0000D-10 * CFACT 
!  Reaction Label TRP39           
             RKI( NCELL,  326) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-12,   5.0000D+02 )
!  Reaction Label TRP40           
             RKI( NCELL,  327) =   1.0000D-10 * CFACT 
!  Reaction Label R303            
             RKI( NCELL,  328) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5000D-12,   5.0000D+02 )
!  Reaction Label R304            
             RKI( NCELL,  329) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5000D-12,   5.0000D+02 )
!  Reaction Label R305            
             RKI( NCELL,  330) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.5100D-13,   5.6500D+02 )
!  Reaction Label R306            
             RKI( NCELL,  331) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.5100D-13,   5.6500D+02 )
!  Reaction Label R307            
             RKI( NCELL,  332) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.5100D-13,   5.6500D+02 )
!  Reaction Label R308            
             RKI( NCELL,  333) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.4000D-14,   2.2100D+02 )
!  Reaction Label R309            
             RKI( NCELL,  334) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.4000D-14,   2.2100D+02 )
!  Reaction Label R310            
             RKI( NCELL,  335) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6800D-12,   5.0000D+02 )
!  Reaction Label R311            
             RKI( NCELL,  336) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6800D-12,   5.0000D+02 )
!  Reaction Label R312            
             RKI( NCELL,  337) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label R313            
             RKI( NCELL,  338) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label R314            
             RKI( NCELL,  339) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.0800D+02 )
!  Reaction Label R315            
             RKI( NCELL,  340) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.0800D+02 )
!  Reaction Label R316            
             RKI( NCELL,  341) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.5100D-13,   5.6500D+02 )
!  Reaction Label R317            
             RKI( NCELL,  342) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.8500D-13,   7.6500D+02 )
!  Reaction Label R318            
             RKI( NCELL,  343) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.3700D-13,   7.6500D+02 )
!  Reaction Label R319            
             RKI( NCELL,  344) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.0800D+02 )
!  Reaction Label R320            
             RKI( NCELL,  345) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-14,   1.5600D+03 )
!  Reaction Label R321            
             RKI( NCELL,  346) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5000D-12,   5.0000D+02 )
!  Reaction Label R322            
             RKI( NCELL,  347) =   1.2000D-12 * CFACT 
!  Reaction Label R323            
             RKI( NCELL,  348) =   1.2000D-12 * CFACT 
!  Reaction Label R324            
             RKI( NCELL,  349) =   1.2000D-12 * CFACT 
!  Reaction Label R325            
             RKI( NCELL,  350) =   1.2000D-12 * CFACT 
!  Reaction Label R327            
             RKI( NCELL,  351) =   1.2000D-12 * CFACT 
!  Reaction Label R328            
             RKI( NCELL,  352) =   1.2000D-12 * CFACT 
!  Reaction Label R329            
             RKI( NCELL,  353) =   1.2000D-12 * CFACT 
!  Reaction Label ROCARO34        
             RKI( NCELL,  354) =   2.3000D-12 * CFACT 
!  Reaction Label ROCARO44        
             RKI( NCELL,  355) =   2.3000D-12 * CFACT 
!  Reaction Label ROCARO54        
             RKI( NCELL,  356) =   2.3000D-12 * CFACT 
!  Reaction Label ROCARO64        
             RKI( NCELL,  357) =   2.3000D-12 * CFACT 
!  Reaction Label R338            
             RKI( NCELL,  358) =   1.2000D-12 * CFACT 
!  Reaction Label R339            
             RKI( NCELL,  359) =   2.3000D-12 * CFACT 
!  Reaction Label R340            
             RKI( NCELL,  360) =   2.3000D-12 * CFACT 
!  Reaction Label TRP53           
             RKI( NCELL,  361) =   2.3000D-12 * CFACT 
!  Reaction Label TRP54           
             RKI( NCELL,  362) =   2.3000D-12 * CFACT 
!  Reaction Label R341            
             RKI( NCELL,  363) =   4.0000D-12 * CFACT 
!  Reaction Label R342            
             RKI( NCELL,  364) =   4.0000D-12 * CFACT 
!  Reaction Label R343            
             RKI( NCELL,  365) =   1.2000D-12 * CFACT 
!  Reaction Label R344            
             RKI( NCELL,  366) =   1.2000D-12 * CFACT 
!  Reaction Label R345            
             RKI( NCELL,  367) =   1.2000D-12 * CFACT 
!  Reaction Label R346            
             RKI( NCELL,  368) =   1.2000D-12 * CFACT 
!  Reaction Label R347            
             RKI( NCELL,  369) =   1.2000D-12 * CFACT 
!  Reaction Label R348            
             RKI( NCELL,  370) =   2.5000D-12 * CFACT 
!  Reaction Label R349            
             RKI( NCELL,  371) =   2.5000D-12 * CFACT 
!  Reaction Label R350            
             RKI( NCELL,  372) =   2.5000D-12 * CFACT 
!  Reaction Label R351            
             RKI( NCELL,  373) =   2.5000D-12 * CFACT 
!  Reaction Label R352            
             RKI( NCELL,  374) =   1.2000D-12 * CFACT 
!  Reaction Label R353            
             RKI( NCELL,  375) =   1.2000D-12 * CFACT 
!  Reaction Label R354            
             RKI( NCELL,  376) =   1.2000D-12 * CFACT 
!  Reaction Label R355            
             RKI( NCELL,  377) =   1.2000D-12 * CFACT 
!  Reaction Label R356            
             RKI( NCELL,  378) =   1.2000D-12 * CFACT 
!  Reaction Label R357            
             RKI( NCELL,  379) =   1.2000D-12 * CFACT 
!  Reaction Label R358            
             RKI( NCELL,  380) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.0000D-14,   1.0000D+03 )
!  Reaction Label R359            
             RKI( NCELL,  381) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.2500D-14,   1.0000D+03 )
!  Reaction Label R360            
             RKI( NCELL,  382) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9600D-14,   1.0000D+03 )
!  Reaction Label R361            
             RKI( NCELL,  383) =   1.2000D-12 * CFACT 
!  Reaction Label R362            
             RKI( NCELL,  384) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.5000D-12,   5.0000D+02 )
!  Reaction Label R363            
             RKI( NCELL,  385) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.1300D-17,   2.9500D+03 )
!  Reaction Label TRP41           
             RKI( NCELL,  386) =   1.0000D-10 * CFACT 
!  Reaction Label TRP42           
             RKI( NCELL,  387) =   1.0000D-10 * CFACT 
!  Reaction Label TRP43           
             RKI( NCELL,  388) =   1.0000D-10 * CFACT 
!  Reaction Label TRP44           
             RKI( NCELL,  389) =   1.0000D-10 * CFACT 
!  Reaction Label TRP45           
             RKI( NCELL,  390) =   1.0000D-10 * CFACT 
!  Reaction Label TRP46           
             RKI( NCELL,  391) =   1.0000D-10 * CFACT 
!  Reaction Label TRP47           
             RKI( NCELL,  392) =   1.0000D-10 * CFACT 
!  Reaction Label TRP48           
             RKI( NCELL,  393) =   1.0000D-10 * CFACT 
!  Reaction Label TRP49           
             RKI( NCELL,  394) =   1.0000D-10 * CFACT 
!  Reaction Label TRP50           
             RKI( NCELL,  395) =   1.0000D-10 * CFACT 
!  Reaction Label TRP51           
             RKI( NCELL,  396) =   1.0000D-10 * CFACT 
!  Reaction Label TRP52           
             RKI( NCELL,  397) =   1.0000D-10 * CFACT 
!  Reaction Label RAM17           
             RKI( NCELL,  398) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.7800D-11,  -4.0000D+02 )
!  Reaction Label R001c           
             RKI( NCELL,  399) =   6.8900D-12 * CFACT 
!  Reaction Label R002c           
             RKI( NCELL,  400) =   6.5500D-14 * CFACT 
!  Reaction Label T17             
             RKI( NCELL,  401) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.0000D-12,   3.8000D+02 )
!  Reaction Label T18             
             RKI( NCELL,  402) =   2.9000D-19 * CFACT 
!  Reaction Label T19             
             RKI( NCELL,  403) =   3.4000D-15 * CFACT 
!  Reaction Label T10             
             RKI( NCELL,  405) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.4800D-11,   4.4800D+02 )
!  Reaction Label T10a            
             RKI( NCELL,  406) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label T10b            
             RKI( NCELL,  407) =   2.3000D-12 * CFACT 
!  Reaction Label T10c            
             RKI( NCELL,  408) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.9300D-13,   1.3000D+03 )
!  Reaction Label T10d            
             RKI( NCELL,  409) =   2.3900D-12 * CFACT 
!  Reaction Label T10e            
             RKI( NCELL,  410) =   1.3700D-11 * CFACT 
!  Reaction Label T11             
             RKI( NCELL,  411) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.3400D-14,  -2.2830D+03 )
!  Reaction Label T12             
             RKI( NCELL,  412) =   1.0000D-13 * CFACT 
!  Reaction Label R003c           
             RKI( NCELL,  413) =   5.0100D-11 * CFACT 
!  Reaction Label R004c           
             RKI( NCELL,  414) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R005c           
             RKI( NCELL,  415) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.7500D-13,   9.8000D+02 )
!  Reaction Label R006c           
             RKI( NCELL,  416) =   4.4000D-11 * CFACT 
!  Reaction Label R007c           
             RKI( NCELL,  417) =   3.4300D-17 * CFACT 
!  Reaction Label R008c           
             RKI( NCELL,  418) =   8.9900D-12 * CFACT 
!  Reaction Label R010c           
             RKI( NCELL,  419) =   1.2000D-11 * CFACT 
!  Reaction Label R011c           
             RKI( NCELL,  420) =   1.9000D-11 * CFACT 
!  Reaction Label R012c           
             RKI( NCELL,  421) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8400D-13,   1.3000D+03 )
!  Reaction Label R013c           
             RKI( NCELL,  422) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label R014c           
             RKI( NCELL,  423) =   2.3000D-12 * CFACT 
!  Reaction Label R015c           
             RKI( NCELL,  424) =   1.2000D-14 * CFACT 
!  Reaction Label R016c           
             RKI( NCELL,  425) =   1.9700D-10 * CFACT 
!  Reaction Label R017c           
             RKI( NCELL,  426) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8400D-13,   1.3000D+03 )
!  Reaction Label R019c           
             RKI( NCELL,  427) =   2.3000D-12 * CFACT 
!  Reaction Label R020c           
             RKI( NCELL,  428) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label HET_GLY         
             RKI( NCELL,  429) =  BLKHET(  NCELL, IK_HETERO_GLY )
!  Reaction Label HET_MGLY        
             RKI( NCELL,  430) =  BLKHET(  NCELL, IK_HETERO_MGLY )
!  Reaction Label HET_NO2         
             RKI( NCELL,  431) =  BLKHET(  NCELL, IK_HETERO_NO2 )
!  Reaction Label HET_HO2         
             RKI( NCELL,  432) =  BLKHET(  NCELL, IK_HETERO_HO2 )
!  Reaction Label HET_NO3         
             RKI( NCELL,  433) =  BLKHET(  NCELL, IK_HETERO_NO3 )
!  Reaction Label HET_IEPOX       
             RKI( NCELL,  434) =  BLKHET(  NCELL, IK_HETERO_IEPOX )
!  Reaction Label HET_ISO3TET     
             RKI( NCELL,  435) =  BLKHET(  NCELL, IK_HETERO_ISO3NOSJ )
!  Reaction Label HET_IEPOXOS     
             RKI( NCELL,  436) =  BLKHET(  NCELL, IK_HETERO_ISO3OSJ )
!  Reaction Label HET_IPX         
             RKI( NCELL,  437) =   2.0000D+00 * BLKHET( NCELL, IK_HETERO_IEPOX )
!  Reaction Label HET_INALD       
             RKI( NCELL,  438) =   5.0000D-01 * BLKHET( NCELL, IK_HETERO_IEPOX )
!  Reaction Label ROCALK1c        
             RKI( NCELL,  439) =   1.5300D-11 * CFACT 
!  Reaction Label ROCALK2c        
             RKI( NCELL,  440) =   1.6800D-11 * CFACT 
!  Reaction Label ROCALK3c        
             RKI( NCELL,  441) =   2.2400D-11 * CFACT 
!  Reaction Label ROCALK4c        
             RKI( NCELL,  442) =   2.6700D-11 * CFACT 
!  Reaction Label ROCALK5c        
             RKI( NCELL,  443) =   3.0900D-11 * CFACT 
!  Reaction Label ROCALK6c        
             RKI( NCELL,  444) =   3.3800D-11 * CFACT 
!  Reaction Label HC1001          
             RKI( NCELL,  445) =   1.1000D-11 * CFACT 
!  Reaction Label ROCALK7c        
             RKI( NCELL,  446) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK8c        
             RKI( NCELL,  447) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK9c        
             RKI( NCELL,  448) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK10c       
             RKI( NCELL,  449) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK11c       
             RKI( NCELL,  450) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK12c       
             RKI( NCELL,  451) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label HC1002          
             RKI( NCELL,  452) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK13c       
             RKI( NCELL,  453) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK14c       
             RKI( NCELL,  454) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK15c       
             RKI( NCELL,  455) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK16c       
             RKI( NCELL,  456) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK17c       
             RKI( NCELL,  457) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK18c       
             RKI( NCELL,  458) =   2.3000D-12 * CFACT 
!  Reaction Label HC1003          
             RKI( NCELL,  459) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK19c       
             RKI( NCELL,  460) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6800D-13,   1.3000D+03 )
!  Reaction Label ROCALK20c       
             RKI( NCELL,  461) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7100D-13,   1.3000D+03 )
!  Reaction Label ROCALK21c       
             RKI( NCELL,  462) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7800D-13,   1.3000D+03 )
!  Reaction Label ROCALK22c       
             RKI( NCELL,  463) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-13,   1.3000D+03 )
!  Reaction Label ROCALK23c       
             RKI( NCELL,  464) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8100D-13,   1.3000D+03 )
!  Reaction Label ROCALK24c       
             RKI( NCELL,  465) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8100D-13,   1.3000D+03 )
!  Reaction Label HC1004          
             RKI( NCELL,  466) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6600D-13,   1.3000D+03 )
!  Reaction Label ROCALK25c       
             RKI( NCELL,  467) =   1.8800D-01 * SFACT 
!  Reaction Label ROCALK26c       
             RKI( NCELL,  468) =   1.8800D-01 * SFACT 
!  Reaction Label ROCALK27c       
             RKI( NCELL,  469) =   1.8800D-01 * SFACT 
!  Reaction Label ROCALK28c       
             RKI( NCELL,  470) =   1.8800D-01 * SFACT 
!  Reaction Label ROCALK29c       
             RKI( NCELL,  471) =   1.8800D-01 * SFACT 
!  Reaction Label ROCALK30c       
             RKI( NCELL,  472) =   1.8800D-01 * SFACT 
!  Reaction Label HC1005          
             RKI( NCELL,  473) =   1.8800D-01 * SFACT 
!  Reaction Label ROCALK31c       
             RKI( NCELL,  474) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK32c       
             RKI( NCELL,  475) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK33c       
             RKI( NCELL,  476) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK34c       
             RKI( NCELL,  477) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK35c       
             RKI( NCELL,  478) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK36c       
             RKI( NCELL,  479) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label HC1006          
             RKI( NCELL,  480) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCALK37c       
             RKI( NCELL,  481) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK38c       
             RKI( NCELL,  482) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK39c       
             RKI( NCELL,  483) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK40c       
             RKI( NCELL,  484) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK41c       
             RKI( NCELL,  485) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK42c       
             RKI( NCELL,  486) =   2.3000D-12 * CFACT 
!  Reaction Label HC1007          
             RKI( NCELL,  487) =   2.3000D-12 * CFACT 
!  Reaction Label ROCALK43c       
             RKI( NCELL,  488) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7100D-13,   1.3000D+03 )
!  Reaction Label ROCALK44c       
             RKI( NCELL,  489) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7300D-13,   1.3000D+03 )
!  Reaction Label ROCALK45c       
             RKI( NCELL,  490) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7800D-13,   1.3000D+03 )
!  Reaction Label ROCALK46c       
             RKI( NCELL,  491) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-13,   1.3000D+03 )
!  Reaction Label ROCALK47c       
             RKI( NCELL,  492) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8100D-13,   1.3000D+03 )
!  Reaction Label ROCALK48c       
             RKI( NCELL,  493) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8200D-13,   1.3000D+03 )
!  Reaction Label HC1008          
             RKI( NCELL,  494) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6600D-13,   1.3000D+03 )
!  Reaction Label ROCARO01        
             RKI( NCELL,  495) =   1.8100D-11 * CFACT 
!  Reaction Label ROCARO02        
             RKI( NCELL,  496) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label ROCARO03        
             RKI( NCELL,  497) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO04        
             RKI( NCELL,  498) =   2.3000D-12 * CFACT 
!  Reaction Label ROCARO05        
             RKI( NCELL,  499) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label ROCARO06        
             RKI( NCELL,  500) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label ROCARO11        
             RKI( NCELL,  501) =   1.8100D-11 * CFACT 
!  Reaction Label ROCARO12        
             RKI( NCELL,  502) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label ROCARO13        
             RKI( NCELL,  503) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO14        
             RKI( NCELL,  504) =   2.3000D-12 * CFACT 
!  Reaction Label ROCARO15        
             RKI( NCELL,  505) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label ROCARO16        
             RKI( NCELL,  506) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label ROCARO21        
             RKI( NCELL,  507) =   2.3100D-11 * CFACT 
!  Reaction Label ROCARO22        
             RKI( NCELL,  508) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-13,   1.3000D+03 )
!  Reaction Label ROCARO23        
             RKI( NCELL,  509) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO24        
             RKI( NCELL,  510) =   2.3000D-12 * CFACT 
!  Reaction Label ROCARO25        
             RKI( NCELL,  511) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-14,   7.0800D+02 )
!  Reaction Label ROCARO26        
             RKI( NCELL,  512) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-13,   7.6500D+02 )
!  Reaction Label ROCOXY1c        
             RKI( NCELL,  513) =   5.9000D-11 * CFACT 
!  Reaction Label ROCOXY2c        
             RKI( NCELL,  514) =   6.0700D-11 * CFACT 
!  Reaction Label ROCOXY3c        
             RKI( NCELL,  515) =   5.5400D-11 * CFACT 
!  Reaction Label ROCOXY4c        
             RKI( NCELL,  516) =   5.6300D-11 * CFACT 
!  Reaction Label ROCOXY5c        
             RKI( NCELL,  517) =   5.4600D-11 * CFACT 
!  Reaction Label ROCOXY6c        
             RKI( NCELL,  518) =   4.5000D-11 * CFACT 
!  Reaction Label ROCOXY7c        
             RKI( NCELL,  519) =   5.1700D-11 * CFACT 
!  Reaction Label ROCOXY8c        
             RKI( NCELL,  520) =   4.7300D-11 * CFACT 
!  Reaction Label ROCOXY9c        
             RKI( NCELL,  521) =   4.6000D-11 * CFACT 
!  Reaction Label ROCOXY10c       
             RKI( NCELL,  522) =   3.8000D-11 * CFACT 
!  Reaction Label ROCOXY11c       
             RKI( NCELL,  523) =   3.9300D-11 * CFACT 
!  Reaction Label ROCOXY12c       
             RKI( NCELL,  524) =   3.5200D-11 * CFACT 
!  Reaction Label ROCOXY13c       
             RKI( NCELL,  525) =   3.1200D-11 * CFACT 
!  Reaction Label ROCOXY14c       
             RKI( NCELL,  526) =   2.4000D-11 * CFACT 
!  Reaction Label ROCOXY15c       
             RKI( NCELL,  527) =   2.0500D-11 * CFACT 
!  Reaction Label ROCOXY16c       
             RKI( NCELL,  528) =   4.7400D-11 * CFACT 
!  Reaction Label R364            
             RKI( NCELL,  529) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.4500D-12,  -1.7750D+03 )
!  Reaction Label TRP58           
             RKI( NCELL,  530) =   9.2600D-05 * SFACT 
!  Reaction Label TRP59           
             RKI( NCELL,  531) =   9.2600D-05 * SFACT 
!  Reaction Label ROCARO71        
             RKI( NCELL,  532) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.0200D-11,   5.3200D+02 )
!  Reaction Label ROCARO71a       
             RKI( NCELL,  533) =   1.4000D-17 * CFACT 
!  Reaction Label ROCARO71b       
             RKI( NCELL,  534) =   1.5100D-13 * CFACT 
!  Reaction Label ROCARO72        
             RKI( NCELL,  535) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9100D-12,   1.3000D+03 )
!  Reaction Label ROCARO73        
             RKI( NCELL,  536) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   3.6000D+02 )
!  Reaction Label ROCARO74        
             RKI( NCELL,  537) =   2.3000D-12 * CFACT 
!  Reaction Label ROCARO75        
             RKI( NCELL,  538) =   2.5000D-13 * CFACT 
!  Reaction Label ROCARO76        
             RKI( NCELL,  539) =   2.5000D-13 * CFACT 
!  Reaction Label HET_ANO3I       
             RKI( NCELL,  540) =  BLKHET(  NCELL, IK_HETERO_ANO3 )
!  Reaction Label HET_ANO3J       
             RKI( NCELL,  541) =  BLKHET(  NCELL, IK_HETERO_ANO3 )
!  Reaction Label CLT01           
             RKI( NCELL,  553) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.3000D-11,  -2.0000D+02 )
!  Reaction Label CLT02           
             RKI( NCELL,  554) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.0500D-11,  -2.2700D+03 )
!  Reaction Label CLT03           
             RKI( NCELL,  555) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.6000D-11,  -3.7500D+02 )
!  Reaction Label CLT04           
             RKI( NCELL,  556) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.4000D-11,   2.7000D+02 )
!  Reaction Label CLT05           
             RKI( NCELL,  557) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1000D-11,  -9.8000D+02 )
!  Reaction Label CLT06           
             RKI( NCELL,  558) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.2000D-11,   0.0000D+00 )
!  Reaction Label CLT07           
             RKI( NCELL,  559) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.3000D-10,   0.0000D+00 )
!  Reaction Label CLT08           
             RKI( NCELL,  560) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-11,   1.6000D+02 )
!  Reaction Label CLT09           
             RKI( NCELL,  561) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.6000D-11,   6.5000D+01 )
!  Reaction Label CLT10           
             RKI( NCELL,  562) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.4000D-12,  -1.3000D+02 )
!  Reaction Label CLT11           
             RKI( NCELL,  563) =  CFACT_SQU * POWER_T02( TEMPOT300,   7.7000D-32,  -1.8000D+00 )
!  Reaction Label CLT12           
             RKI( NCELL,  564) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 1.8000D-31,   0.0000D+00,  -2.0000D+00,  & 
     &                                                 1.0000D-10,   0.0000D+00,  -1.0000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label CLT14           
             RKI( NCELL,  565) =   2.4000D-11 * CFACT 
!  Reaction Label CLT15           
             RKI( NCELL,  566) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.8000D-11,   1.0000D+02 )
!  Reaction Label CLT17           
             RKI( NCELL,  567) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.2000D-12,   1.4500D+02 )
!  Reaction Label CLT22           
             RKI( NCELL,  568) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-12,   2.9000D+02 )
!  Reaction Label CLT23           
             RKI( NCELL,  569) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.4000D-12,   2.9000D+02 )
!  Reaction Label CLT25           
             RKI( NCELL,  570) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 1.8000D-31,   0.0000D+00,  -3.4000D+00,  & 
     &                                                 1.5000D-11,   0.0000D+00,  -1.9000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label CLT27           
             RKI( NCELL,  571) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.0000D-12,  -1.5900D+03 )
!  Reaction Label CLT28           
             RKI( NCELL,  572) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5000D-13,  -1.3700D+03 )
!  Reaction Label CLT31           
             RKI( NCELL,  573) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 1.9000D-32,   0.0000D+00,  -3.6000D+00,  & 
     &                                                 3.7000D-12,   0.0000D+00,  -1.6000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label CLT32           
             RKI( NCELL,  574) =  RFACT * RKI( NCELL,  573 ) & 
     &                         * (  4.6296D+26 * DEXP( -8.5370D+03 * INV_TEMP) ) 
!  Reaction Label CLT35           
             RKI( NCELL,  575) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1000D-13,   3.5000D+02 )
!  Reaction Label CLT36           
             RKI( NCELL,  576) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-12,  -1.1000D+03 )
!  Reaction Label CLT37           
             RKI( NCELL,  577) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.0000D-13,   2.3000D+02 )
!  Reaction Label CLT38           
             RKI( NCELL,  578) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.4000D-12,   2.7000D+02 )
!  Reaction Label CLT39           
             RKI( NCELL,  579) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.8000D-12,  -2.5000D+02 )
!  Reaction Label CLT40           
             RKI( NCELL,  580) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.0000D-12,  -5.0000D+02 )
!  Reaction Label CLT41           
             RKI( NCELL,  581) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.4000D-12,  -1.2500D+03 )
!  Reaction Label CLT43           
             RKI( NCELL,  582) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.2000D-12,  -3.3000D+02 )
!  Reaction Label CLT44           
             RKI( NCELL,  583) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.4000D-12,   6.0000D+02 )
!  Reaction Label CLT45           
             RKI( NCELL,  584) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.0000D-13,   6.7000D+02 )
!  Reaction Label CLTO01          
             RKI( NCELL,  585) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.1000D-12,  -1.2700D+03 )
!  Reaction Label CLTO02          
             RKI( NCELL,  586) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.2000D-11,  -7.0000D+01 )
!  Reaction Label CLTO03          
             RKI( NCELL,  587) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.4000D-10,   0.0000D+00 )
!  Reaction Label CLTO04          
             RKI( NCELL,  588) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-10,   0.0000D+00 )
!  Reaction Label CLTO05          
             RKI( NCELL,  589) =   5.2700D-10 * CFACT 
!  Reaction Label CLTO06          
             RKI( NCELL,  590) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.1000D-11,  -3.0000D+01 )
!  Reaction Label CLTO07          
             RKI( NCELL,  591) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.0000D-11,   0.0000D+00 )
!  Reaction Label CLTO08          
             RKI( NCELL,  592) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-11,   4.5300D+02 )
!  Reaction Label CLTO09          
             RKI( NCELL,  593) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-11,  -7.4500D+02 )
!  Reaction Label CLTO10          
             RKI( NCELL,  594) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.5000D-11,   0.0000D+00 )
!  Reaction Label CLTO11          
             RKI( NCELL,  595) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.0000D-11,   1.5500D+02 )
!  Reaction Label CLTO12          
             RKI( NCELL,  596) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7000D-10,   0.0000D+00 )
!  Reaction Label CLTO13          
             RKI( NCELL,  597) =   2.0000D-13 * CFACT 
!  Reaction Label CLTO14          
             RKI( NCELL,  598) =   2.7000D-14 * CFACT 
!  Reaction Label CLTO15          
             RKI( NCELL,  599) =   2.3000D-10 * CFACT 
!  Reaction Label CTLO16          
             RKI( NCELL,  600) =   2.3000D-10 * CFACT 
!  Reaction Label CLTO17          
             RKI( NCELL,  601) =   4.1000D-11 * CFACT 
!  Reaction Label CLTO18          
             RKI( NCELL,  602) =   1.1600D-10 * CFACT 
!  Reaction Label CLTO19          
             RKI( NCELL,  603) =   5.4000D-11 * CFACT 
!  Reaction Label CLTO20          
             RKI( NCELL,  604) =   2.6000D-10 * CFACT 
!  Reaction Label CLTO21          
             RKI( NCELL,  605) =   2.0000D-10 * CFACT 
!  Reaction Label CLTO22          
             RKI( NCELL,  606) =   3.6000D-11 * CFACT 
!  Reaction Label CLTO23          
             RKI( NCELL,  607) =   5.0000D-11 * CFACT 
!  Reaction Label CLTO24          
             RKI( NCELL,  608) =   4.6000D-10 * CFACT 
!  Reaction Label CLTO25          
             RKI( NCELL,  609) =   5.9000D-11 * CFACT 
!  Reaction Label CLTO26          
             RKI( NCELL,  610) =   1.1000D-10 * CFACT 
!  Reaction Label CLTO27          
             RKI( NCELL,  611) =   1.1000D-10 * CFACT 
!  Reaction Label CLTO28          
             RKI( NCELL,  612) =   1.1000D-10 * CFACT 
!  Reaction Label CLTO29          
             RKI( NCELL,  613) =   1.3600D-09 * CFACT 
!  Reaction Label CLTO29a         
             RKI( NCELL,  614) =   1.5600D-09 * CFACT 
!  Reaction Label CLTO30          
             RKI( NCELL,  615) =   5.1000D-10 * CFACT 
!  Reaction Label CLTO31          
             RKI( NCELL,  616) =   2.2000D-13 * CFACT 
!  Reaction Label CLTO32          
             RKI( NCELL,  617) =   1.0000D-14 * CFACT 
!  Reaction Label CLTO33          
             RKI( NCELL,  618) =   1.1300D-12 * CFACT 
!  Reaction Label CLTO34          
             RKI( NCELL,  619) =   9.1000D-10 * CFACT 
!  Reaction Label CLTO35          
             RKI( NCELL,  620) =   2.2000D-11 * CFACT 
!  Reaction Label CLTO36          
             RKI( NCELL,  621) =   4.8000D-11 * CFACT 
!  Reaction Label CLTO37          
             RKI( NCELL,  622) =   9.1000D-11 * CFACT 
!  Reaction Label CLTO38          
             RKI( NCELL,  623) =   1.3000D-10 * CFACT 
!  Reaction Label CLTO39          
             RKI( NCELL,  624) =   1.8000D-09 * CFACT 
!  Reaction Label CLTO40          
             RKI( NCELL,  625) =   2.1000D-10 * CFACT 
!  Reaction Label CLTO41          
             RKI( NCELL,  626) =   9.5000D-11 * CFACT 
!  Reaction Label CLTO42          
             RKI( NCELL,  627) =   3.0000D-10 * CFACT 
!  Reaction Label CLTO43          
             RKI( NCELL,  628) =   3.5000D-10 * CFACT 
!  Reaction Label CLTO44          
             RKI( NCELL,  629) =   4.6000D-10 * CFACT 
!  Reaction Label CLTO45          
             RKI( NCELL,  630) =   2.4000D-10 * CFACT 
!  Reaction Label CLTO46          
             RKI( NCELL,  631) =   7.4000D-10 * CFACT 
!  Reaction Label CLTO47          
             RKI( NCELL,  632) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.0600D-11,   4.8000D+02 )
!  Reaction Label CLTO48          
             RKI( NCELL,  633) =   2.6000D-10 * CFACT 
!  Reaction Label CLTO49          
             RKI( NCELL,  634) =   3.7000D-10 * CFACT 
!  Reaction Label CLTO50          
             RKI( NCELL,  635) =   3.7000D-10 * CFACT 
!  Reaction Label CLTO51          
             RKI( NCELL,  636) =   7.4000D-11 * CFACT 
!  Reaction Label CLTO52          
             RKI( NCELL,  637) =   5.3000D-10 * CFACT 
!  Reaction Label CLTO53          
             RKI( NCELL,  638) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.2300D-10,   0.0000D+00 )
!  Reaction Label CLTO54          
             RKI( NCELL,  639) =   3.3000D-10 * CFACT 
!  Reaction Label CLTO55          
             RKI( NCELL,  640) =   1.5000D-16 * CFACT 
!  Reaction Label CLTO56          
             RKI( NCELL,  641) =   6.0000D-11 * CFACT 
!  Reaction Label CLTO58          
             RKI( NCELL,  642) =   1.4500D-10 * CFACT 
!  Reaction Label CLTO59          
             RKI( NCELL,  643) =   9.9000D-11 * CFACT 
!  Reaction Label CLTO60          
             RKI( NCELL,  644) =   2.0000D-10 * CFACT 
!  Reaction Label CLTO61          
             RKI( NCELL,  645) =   3.5000D-11 * CFACT 
!  Reaction Label CLTO62          
             RKI( NCELL,  646) =   1.6000D-09 * CFACT 
!  Reaction Label CLTO63          
             RKI( NCELL,  647) =   2.8000D-15 * CFACT 
!  Reaction Label CLTO64          
             RKI( NCELL,  648) =   2.0000D-10 * CFACT 
!  Reaction Label CLTO65          
             RKI( NCELL,  649) =   1.8200D-09 * CFACT 
!  Reaction Label CLTO66          
             RKI( NCELL,  650) =   1.6700D-09 * CFACT 
!  Reaction Label CLTO67          
             RKI( NCELL,  651) =   1.4500D-09 * CFACT 
!  Reaction Label CLTO68          
             RKI( NCELL,  652) =   1.2240D-09 * CFACT 
!  Reaction Label CLTO69          
             RKI( NCELL,  653) =   9.3000D-10 * CFACT 
!  Reaction Label CLTO70          
             RKI( NCELL,  654) =   8.5100D-10 * CFACT 
!  Reaction Label CLTO71          
             RKI( NCELL,  655) =   3.4300D-11 * CFACT 
!  Reaction Label CLTO72          
             RKI( NCELL,  656) =   1.3700D-10 * CFACT 
!  Reaction Label CLTO73          
             RKI( NCELL,  657) =   1.3700D-10 * CFACT 
!  Reaction Label CLTO74          
             RKI( NCELL,  658) =   3.0900D-10 * CFACT 
!  Reaction Label CLTO75          
             RKI( NCELL,  659) =   3.7700D-10 * CFACT 
!  Reaction Label CLTO76          
             RKI( NCELL,  660) =   3.8300D-10 * CFACT 
!  Reaction Label CLTO77          
             RKI( NCELL,  661) =   4.2000D-10 * CFACT 
!  Reaction Label CLTO78          
             RKI( NCELL,  662) =   3.8900D-10 * CFACT 
!  Reaction Label CLTO79          
             RKI( NCELL,  663) =   4.0800D-10 * CFACT 
!  Reaction Label CLTO80          
             RKI( NCELL,  664) =   3.2500D-10 * CFACT 
!  Reaction Label CLTO81          
             RKI( NCELL,  665) =   3.5600D-10 * CFACT 
!  Reaction Label CLTO82          
             RKI( NCELL,  666) =   2.5900D-10 * CFACT 
!  Reaction Label CLTO83          
             RKI( NCELL,  667) =   3.1600D-10 * CFACT 
!  Reaction Label CLTO84          
             RKI( NCELL,  668) =   2.6800D-10 * CFACT 
!  Reaction Label CLTO85          
             RKI( NCELL,  669) =   2.3900D-10 * CFACT 
!  Reaction Label CLTO86          
             RKI( NCELL,  670) =   2.1100D-10 * CFACT 
!  Reaction Label CLTO87          
             RKI( NCELL,  671) =   1.6000D-10 * CFACT 
!  Reaction Label CLTO88          
             RKI( NCELL,  672) =   1.3500D-10 * CFACT 
!  Reaction Label CLTO89          
             RKI( NCELL,  673) =   3.8800D-11 * CFACT 
!  Reaction Label CLTO90          
             RKI( NCELL,  674) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.5600D-11,  -1.4190D+03 )
!  Reaction Label CLTO91          
             RKI( NCELL,  675) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.1000D-11,  -3.0000D+01 )
!  Reaction Label CLT092          
             RKI( NCELL,  676) =   1.2900D-10 * CFACT 
!  Reaction Label CLT093          
             RKI( NCELL,  677) =   3.6000D-10 * CFACT 
!  Reaction Label CLOR01          
             RKI( NCELL,  678) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR02          
             RKI( NCELL,  679) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR03          
             RKI( NCELL,  680) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR04          
             RKI( NCELL,  681) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR05          
             RKI( NCELL,  682) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR06          
             RKI( NCELL,  683) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR07          
             RKI( NCELL,  684) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR08          
             RKI( NCELL,  685) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR09          
             RKI( NCELL,  686) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR10          
             RKI( NCELL,  687) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR12          
             RKI( NCELL,  688) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR13          
             RKI( NCELL,  689) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR14          
             RKI( NCELL,  690) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR15          
             RKI( NCELL,  691) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR16          
             RKI( NCELL,  692) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR17          
             RKI( NCELL,  693) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR18          
             RKI( NCELL,  694) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR19          
             RKI( NCELL,  695) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR20          
             RKI( NCELL,  696) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR21          
             RKI( NCELL,  697) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR22          
             RKI( NCELL,  698) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR23          
             RKI( NCELL,  699) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR24          
             RKI( NCELL,  700) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR25          
             RKI( NCELL,  701) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR26          
             RKI( NCELL,  702) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR27          
             RKI( NCELL,  703) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR28          
             RKI( NCELL,  704) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR29          
             RKI( NCELL,  705) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR30          
             RKI( NCELL,  706) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR31          
             RKI( NCELL,  707) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR32          
             RKI( NCELL,  708) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR33          
             RKI( NCELL,  709) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR34          
             RKI( NCELL,  710) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR35          
             RKI( NCELL,  711) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR36          
             RKI( NCELL,  712) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR37          
             RKI( NCELL,  713) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR38          
             RKI( NCELL,  714) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR39          
             RKI( NCELL,  715) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR40          
             RKI( NCELL,  716) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR41          
             RKI( NCELL,  717) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR42          
             RKI( NCELL,  718) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR43          
             RKI( NCELL,  719) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR44          
             RKI( NCELL,  720) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR45          
             RKI( NCELL,  721) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR46          
             RKI( NCELL,  722) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR47          
             RKI( NCELL,  723) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR48          
             RKI( NCELL,  724) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR49          
             RKI( NCELL,  725) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR50          
             RKI( NCELL,  726) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR51          
             RKI( NCELL,  727) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR52          
             RKI( NCELL,  728) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR53          
             RKI( NCELL,  729) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR54          
             RKI( NCELL,  730) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR55          
             RKI( NCELL,  731) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR56          
             RKI( NCELL,  732) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR57          
             RKI( NCELL,  733) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR58          
             RKI( NCELL,  734) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR59          
             RKI( NCELL,  735) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR60          
             RKI( NCELL,  736) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR61          
             RKI( NCELL,  737) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR62          
             RKI( NCELL,  738) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label CLOR63          
             RKI( NCELL,  739) =  CFACT * ARRHENUIS_T03( INV_TEMP,  3.2500D-12,  -1.1400D+02 )
!  Reaction Label BRT01           
             RKI( NCELL,  751) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.1000D-11,   2.4000D+02 )
!  Reaction Label BRT02           
             RKI( NCELL,  752) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6000D-11,  -7.8000D+02 )
!  Reaction Label BRT03           
             RKI( NCELL,  753) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.8000D-12,  -3.1000D+02 )
!  Reaction Label BRT06           
             RKI( NCELL,  754) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 4.3000D-31,   0.0000D+00,  -2.4000D+00,  & 
     &                                                 2.7000D-11,   0.0000D+00,   0.0000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label BRT07           
             RKI( NCELL,  755) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.6000D-11,   0.0000D+00 )
!  Reaction Label BRT08           
             RKI( NCELL,  756) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.9000D-11,   0.0000D+00 )
!  Reaction Label BRT10           
             RKI( NCELL,  757) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.5000D-12,   2.0000D+02 )
!  Reaction Label BRT11           
             RKI( NCELL,  758) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7000D-11,   2.5000D+02 )
!  Reaction Label BRT12           
             RKI( NCELL,  759) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.5000D-12,   4.6000D+02 )
!  Reaction Label BRT13           
             RKI( NCELL,  760) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.8000D-12,   2.6000D+02 )
!  Reaction Label BRT15           
             RKI( NCELL,  761) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 5.5000D-31,   0.0000D+00,  -3.1000D+00,  & 
     &                                                 6.6000D-12,   0.0000D+00,  -2.9000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label BRT16           
             RKI( NCELL,  762) =   1.0000D-12 * CFACT 
!  Reaction Label BRT17           
             RKI( NCELL,  763) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.4000D-12,   4.0000D+01 )
!  Reaction Label BRT18           
             RKI( NCELL,  764) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.8000D-14,   8.6000D+02 )
!  Reaction Label BRT19           
             RKI( NCELL,  765) =   1.5000D-15 * CFACT 
!  Reaction Label BRTO01          
             RKI( NCELL,  766) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7000D-11,  -8.0000D+02 )
!  Reaction Label BRTO02          
             RKI( NCELL,  767) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.8000D-11,  -4.6000D+02 )
!  Reaction Label BRTO03          
             RKI( NCELL,  768) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.3000D-11,  -5.2600D+02 )
!  Reaction Label BRTO04          
             RKI( NCELL,  769) =   8.0000D-15 * CFACT 
!  Reaction Label BRTO05          
             RKI( NCELL,  770) =   2.0000D-14 * CFACT 
!  Reaction Label BRTO06          
             RKI( NCELL,  771) =   3.3000D-14 * CFACT 
!  Reaction Label BRTO07          
             RKI( NCELL,  772) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.7000D-11,  -5.3600D+02 )
!  Reaction Label BRTO08          
             RKI( NCELL,  773) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.7000D-11,  -5.3600D+02 )
!  Reaction Label BRTO09          
             RKI( NCELL,  774) =   8.0000D-15 * CFACT 
!  Reaction Label BRTO10          
             RKI( NCELL,  775) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7000D-11,  -8.0000D+02 )
!  Reaction Label BRTO11          
             RKI( NCELL,  776) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.7000D-11,  -8.0000D+02 )
!  Reaction Label BRTO12          
             RKI( NCELL,  777) =   2.0000D-14 * CFACT 
!  Reaction Label BRTO13          
             RKI( NCELL,  778) =   6.4000D-16 * CFACT 
!  Reaction Label BRTO13a         
             RKI( NCELL,  779) =   7.3600D-15 * CFACT 
!  Reaction Label BRTO14          
             RKI( NCELL,  780) =   8.0000D-15 * CFACT 
!  Reaction Label BRTO15          
             RKI( NCELL,  781) =   1.3000D-13 * CFACT 
!  Reaction Label BRTO16          
             RKI( NCELL,  782) =   3.7000D-12 * CFACT 
!  Reaction Label BRTO17          
             RKI( NCELL,  783) =   1.0000D-11 * CFACT 
!  Reaction Label BRTO18          
             RKI( NCELL,  784) =   8.1000D-11 * CFACT 
!  Reaction Label BRTO19          
             RKI( NCELL,  785) =   2.4000D-11 * CFACT 
!  Reaction Label BRTO20          
             RKI( NCELL,  786) =   7.4000D-11 * CFACT 
!  Reaction Label BRTO21          
             RKI( NCELL,  787) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.4000D-15,   4.4000D+02 )
!  Reaction Label BRTO22          
             RKI( NCELL,  788) =   2.0000D-11 * CFACT 
!  Reaction Label BRTO23          
             RKI( NCELL,  789) =   2.6000D-11 * CFACT 
!  Reaction Label BRTO24          
             RKI( NCELL,  790) =   3.7000D-11 * CFACT 
!  Reaction Label BRTO25          
             RKI( NCELL,  791) =   3.7000D-11 * CFACT 
!  Reaction Label BRTO26          
             RKI( NCELL,  792) =   7.4000D-12 * CFACT 
!  Reaction Label BRTO27          
             RKI( NCELL,  793) =   5.3000D-11 * CFACT 
!  Reaction Label BRTO28          
             RKI( NCELL,  794) =   4.3000D-12 * CFACT 
!  Reaction Label BRTO29          
             RKI( NCELL,  795) =   2.5000D-11 * CFACT 
!  Reaction Label BRTO30          
             RKI( NCELL,  796) =   6.3000D-11 * CFACT 
!  Reaction Label BRTO31          
             RKI( NCELL,  797) =   5.0000D-16 * CFACT 
!  Reaction Label BRTO32          
             RKI( NCELL,  798) =   1.4500D-14 * CFACT 
!  Reaction Label BRTO34          
             RKI( NCELL,  799) =   7.2500D-14 * CFACT 
!  Reaction Label BRTO35          
             RKI( NCELL,  800) =   1.0000D-13 * CFACT 
!  Reaction Label BRT036          
             RKI( NCELL,  801) =   7.2500D-14 * CFACT 
!  Reaction Label BRTO37          
             RKI( NCELL,  802) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.2500D-14,   0.0000D+00 )
!  Reaction Label BRTO38          
             RKI( NCELL,  803) =  CFACT * ARRHENUIS_T03( INV_TEMP,  7.2500D-14,   0.0000D+00 )
!  Reaction Label BRTO39          
             RKI( NCELL,  804) =   3.0000D-16 * CFACT 
!  Reaction Label BRTO40          
             RKI( NCELL,  805) =   5.0100D-12 * CFACT 
!  Reaction Label BRTO41          
             RKI( NCELL,  806) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0000D-12,  -8.4000D+02 )
!  Reaction Label BRTO42          
             RKI( NCELL,  807) =   8.5000D-12 * CFACT 
!  Reaction Label BRTO43          
             RKI( NCELL,  808) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.0000D-12,  -3.6000D+02 )
!  Reaction Label BRT044          
             RKI( NCELL,  809) =   3.1400D-14 * CFACT 
!  Reaction Label BRT045          
             RKI( NCELL,  810) =   5.6000D-12 * CFACT 
!  Reaction Label IT02            
             RKI( NCELL,  824) =   1.8000D-10 * CFACT 
!  Reaction Label IT03            
             RKI( NCELL,  825) =   1.5000D-12 * CFACT 
!  Reaction Label IT04            
             RKI( NCELL,  826) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.1000D-11,  -8.3000D+02 )
!  Reaction Label IT05            
             RKI( NCELL,  827) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.5000D-11,  -1.0900D+03 )
!  Reaction Label IT06            
             RKI( NCELL,  828) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 1.8000D-32,   0.0000D+00,  -1.0000D+00,  & 
     &                                                 1.7000D-11,   0.0000D+00,   0.0000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label IT07            
             RKI( NCELL,  829) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 3.0000D-31,   0.0000D+00,  -1.0000D+00,  & 
     &                                                 6.6000D-11,   0.0000D+00,   0.0000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label IT09            
             RKI( NCELL,  830) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.1000D-11,  -1.4600D+02 )
!  Reaction Label IT10            
             RKI( NCELL,  831) =  SFACT * ARRHENUIS_T03( INV_TEMP,  9.9400D+17,  -1.1859D+04 )
!  Reaction Label IT11            
             RKI( NCELL,  832) =  SFACT * ARRHENUIS_T03( INV_TEMP,  2.1000D+15,  -1.3670D+04 )
!  Reaction Label IT12            
             RKI( NCELL,  833) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.4000D-11,  -2.6200D+03 )
!  Reaction Label IT13            
             RKI( NCELL,  834) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9000D-11,  -2.6000D+03 )
!  Reaction Label IT15            
             RKI( NCELL,  835) =   5.0000D-12 * CFACT 
!  Reaction Label IT16            
             RKI( NCELL,  836) =   3.0000D-11 * CFACT 
!  Reaction Label IT20            
             RKI( NCELL,  837) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.4000D-11,   5.4000D+02 )
!  Reaction Label IT21            
             RKI( NCELL,  838) =  CFACT * ARRHENUIS_T03( INV_TEMP,  8.6000D-12,   2.3000D+02 )
!  Reaction Label IT22            
             RKI( NCELL,  839) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.1000D-12,   5.4200D+02 )
!  Reaction Label IT23            
             RKI( NCELL,  840) =  CFACT * FALLOFF_T10( INV_TEMP,  TEMPOT300,  CAIR, & 
     &                                                 7.7000D-31,   0.0000D+00,  -3.5000D+00,  & 
     &                                                 7.7000D-12,   0.0000D+00,  -1.5000D+00,  & 
     &                                                 1.0000D+00,   6.0000D-01 )
!  Reaction Label IT25            
             RKI( NCELL,  841) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.5000D-11,   5.0000D+02 )
!  Reaction Label IT30            
             RKI( NCELL,  842) =  SFACT * ARRHENUIS_T03( INV_TEMP,  2.5000D+14,  -9.7700D+03 )
!  Reaction Label IT33            
             RKI( NCELL,  843) =  CFACT * POWER_T02( TEMPOT300,   2.7000D-12,  -2.6600D+00 )
!  Reaction Label IT34            
             RKI( NCELL,  844) =  SFACT * ARRHENUIS_T03( INV_TEMP,  1.0000D+12,  -9.7700D+03 )
!  Reaction Label IT35            
             RKI( NCELL,  845) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.5000D-10,   0.0000D+00 )
!  Reaction Label IT36            
             RKI( NCELL,  846) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.5000D-10,   0.0000D+00 )
!  Reaction Label IT37            
             RKI( NCELL,  847) =  SFACT * ARRHENUIS_T03( INV_TEMP,  3.8000D-02,   0.0000D+00 )
!  Reaction Label IT38            
             RKI( NCELL,  848) =  CFACT * ARRHENUIS_T03( INV_TEMP,  1.3000D-12,  -1.8300D+03 )
!  Reaction Label IT39            
             RKI( NCELL,  849) =   1.2500D-10 * CFACT 
!  Reaction Label ITO01           
             RKI( NCELL,  850) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.3000D-12,  -1.1200D+03 )
!  Reaction Label MT01            
             RKI( NCELL,  858) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.6000D-11,  -1.3000D+03 )
!  Reaction Label MT02            
             RKI( NCELL,  859) =  CFACT * ARRHENUIS_T03( INV_TEMP,  5.9000D-12,  -1.7000D+02 )
!  Reaction Label MT03            
             RKI( NCELL,  860) =   1.2000D-11 * CFACT 
!  Reaction Label MT04            
             RKI( NCELL,  861) =   2.7000D-11 * CFACT 
!  Reaction Label MT05            
             RKI( NCELL,  862) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.5000D-13,   5.5000D+02 )
!  Reaction Label MT06            
             RKI( NCELL,  863) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.9000D-12,   2.2000D+02 )
!  Reaction Label MT07            
             RKI( NCELL,  864) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.7000D-12,   2.8000D+02 )
!  Reaction Label MT08            
             RKI( NCELL,  865) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.4000D-12,   7.6000D+02 )
!  Reaction Label MTO01           
             RKI( NCELL,  866) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.8000D-12,  -8.3000D+02 )
!  Reaction Label MTO02           
             RKI( NCELL,  867) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.3500D-12,  -8.5000D+02 )
!  Reaction Label MTO03           
             RKI( NCELL,  868) =  CFACT * ARRHENUIS_T03( INV_TEMP,  6.8000D-12,  -8.7000D+02 )
!  Reaction Label MTO06           
             RKI( NCELL,  869) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.1000D-13,  -8.8000D+02 )
!  Reaction Label MTO07           
             RKI( NCELL,  870) =  CFACT * ARRHENUIS_T03( INV_TEMP,  4.8500D-12,  -8.5000D+02 )
!  Reaction Label MTO08           
             RKI( NCELL,  871) =  CFACT * ARRHENUIS_T03( INV_TEMP,  2.0300D-11,  -1.1100D+03 )
!  Reaction Label MTO09           
             RKI( NCELL,  872) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.0000D-13,  -4.2000D+02 )
!  Reaction Label MTO10           
             RKI( NCELL,  873) =  CFACT * ARRHENUIS_T03( INV_TEMP,  9.4000D-13,  -5.1000D+02 )
!  Reaction Label HET_N2O5IJ      
             RKI( NCELL,  874) =  BLKHET(  NCELL, IK_HETERO_N2O5IJ )
!  Reaction Label HET_N2O5K       
             RKI( NCELL,  875) =  BLKHET(  NCELL, IK_HETERO_N2O5K )
!  Reaction Label HET_H2NO3PIJA   
             RKI( NCELL,  876) =  BLKHET(  NCELL, IK_HETERO_H2NO3PAIJ )
!  Reaction Label HET_H2NO3PKA    
             RKI( NCELL,  877) =  BLKHET(  NCELL, IK_HETERO_H2NO3PAK )
!  Reaction Label HET_H2NO3PIB    
             RKI( NCELL,  878) =  BLKHET(  NCELL, IK_HETERO_H2NO3PBIJ )
!  Reaction Label HET_H2NO3PJB    
             RKI( NCELL,  879) =  BLKHET(  NCELL, IK_HETERO_H2NO3PBIJ )
!  Reaction Label HET_H2NO3PKB    
             RKI( NCELL,  880) =  BLKHET(  NCELL, IK_HETERO_H2NO3PBK )
!  Reaction Label HET_HOCL_ACLJ   
             RKI( NCELL,  881) =  BLKHET(  NCELL, IK_HETERO_HOCL_ACLJ )
!  Reaction Label HET_HOCL_ABRJ   
             RKI( NCELL,  882) =  BLKHET(  NCELL, IK_HETERO_HOCL_ABRJ )
!  Reaction Label HET_CLN2_WAI    
             RKI( NCELL,  883) =  BLKHET(  NCELL, IK_HETERO_CLN2_WAI )
!  Reaction Label HET_CLN2_WAJ    
             RKI( NCELL,  884) =  BLKHET(  NCELL, IK_HETERO_CLN2_WAJ )
!  Reaction Label HET_CLN2_ACLJ   
             RKI( NCELL,  885) =  BLKHET(  NCELL, IK_HETERO_CLN2_ACLJ )
!  Reaction Label HET_CLN2_ABRJ   
             RKI( NCELL,  886) =  BLKHET(  NCELL, IK_HETERO_CLN2_ABRJ )
!  Reaction Label HET_CLN3_WAI    
             RKI( NCELL,  887) =  BLKHET(  NCELL, IK_HETERO_CLN3_WAI )
!  Reaction Label HET_CLN3_WAJ    
             RKI( NCELL,  888) =  BLKHET(  NCELL, IK_HETERO_CLN3_WAJ )
!  Reaction Label HET_CLN3_ACLJ   
             RKI( NCELL,  889) =  BLKHET(  NCELL, IK_HETERO_CLN3_ACLJ )
!  Reaction Label HET_CLN3_ABRJ   
             RKI( NCELL,  890) =  BLKHET(  NCELL, IK_HETERO_CLN3_ABRJ )
!  Reaction Label HET_HOBR_ACLJ   
             RKI( NCELL,  891) =  BLKHET(  NCELL, IK_HETERO_HOBR_ACLJ )
!  Reaction Label HET_HOBR_ABRJ   
             RKI( NCELL,  892) =  BLKHET(  NCELL, IK_HETERO_HOBR_ABRJ )
!  Reaction Label HET_BRN2_ACLJ   
             RKI( NCELL,  893) =  BLKHET(  NCELL, IK_HETERO_BRN2_ACLJ )
!  Reaction Label HET_BRN2_ABRJ   
             RKI( NCELL,  894) =  BLKHET(  NCELL, IK_HETERO_BRN2_ABRJ )
!  Reaction Label HET_BRN3_WAI    
             RKI( NCELL,  895) =  BLKHET(  NCELL, IK_HETERO_BRN3_WAI )
!  Reaction Label HET_BRN3_WAJ    
             RKI( NCELL,  896) =  BLKHET(  NCELL, IK_HETERO_BRN3_WAJ )
!  Reaction Label HET_BRN3_ACLJ   
             RKI( NCELL,  897) =  BLKHET(  NCELL, IK_HETERO_BRN3_ACLJ )
!  Reaction Label HET_BRN3_ABRJ   
             RKI( NCELL,  898) =  BLKHET(  NCELL, IK_HETERO_BRN3_ABRJ )
!  Reaction Label HET_HBR_ABRJ    
             RKI( NCELL,  899) =  BLKHET(  NCELL, IK_HETERO_HBR_ABRJ )
!  Reaction Label HET_HOI_ACLJ    
             RKI( NCELL,  900) =  BLKHET(  NCELL, IK_HETERO_HOI_ACLJ )
!  Reaction Label HET_HOI_ABRJ    
             RKI( NCELL,  901) =  BLKHET(  NCELL, IK_HETERO_HOI_ABRJ )
!  Reaction Label HET_INO2_ACLJ   
             RKI( NCELL,  902) =  BLKHET(  NCELL, IK_HETERO_INO2_ACLJ )
!  Reaction Label HET_INO2_ABRJ   
             RKI( NCELL,  903) =  BLKHET(  NCELL, IK_HETERO_INO2_ABRJ )
!  Reaction Label HET_INO3_ACLJ   
             RKI( NCELL,  904) =  BLKHET(  NCELL, IK_HETERO_INO3_ACLJ )
!  Reaction Label HET_INO3_ABRJ   
             RKI( NCELL,  905) =  BLKHET(  NCELL, IK_HETERO_INO3_ABRJ )
!  Reaction Label HET_I2O2_AI     
             RKI( NCELL,  906) =  BLKHET(  NCELL, IK_HETERO_I2O2_AI )
!  Reaction Label HET_I2O2_AJ     
             RKI( NCELL,  907) =  BLKHET(  NCELL, IK_HETERO_I2O2_AJ )
!  Reaction Label HET_I2O3_AI     
             RKI( NCELL,  908) =  BLKHET(  NCELL, IK_HETERO_I2O3_AI )
!  Reaction Label HET_I2O3_AJ     
             RKI( NCELL,  909) =  BLKHET(  NCELL, IK_HETERO_I2O3_AJ )
!  Reaction Label HET_I2O4_AI     
             RKI( NCELL,  910) =  BLKHET(  NCELL, IK_HETERO_I2O4_AI )
!  Reaction Label HET_I2O4_AJ     
             RKI( NCELL,  911) =  BLKHET(  NCELL, IK_HETERO_I2O4_AJ )
!  Reaction Label HET_SO2_H2O2    
             RKI( NCELL,  912) =  BLKHET(  NCELL, IK_HETERO_SO2H )
!  Reaction Label HET_SO2_O3      
             RKI( NCELL,  913) =  BLKHET(  NCELL, IK_HETERO_SO2O )
!  Reaction Label HET_SO2_MEPX    
             RKI( NCELL,  914) =  BLKHET(  NCELL, IK_HETERO_SO2M )
!  Reaction Label HET_SO2_PACD    
             RKI( NCELL,  915) =  BLKHET(  NCELL, IK_HETERO_SO2P )
!  Reaction Label HET_SO2_NO2     
             RKI( NCELL,  916) =  BLKHET(  NCELL, IK_HETERO_SO2N )
!  Reaction Label HET_SO2_TMI     
             RKI( NCELL,  917) =  BLKHET(  NCELL, IK_HETERO_SO2T )
!  Reaction Label HET_SO2_HCHO    
             RKI( NCELL,  918) =  BLKHET(  NCELL, IK_HETERO_HMSP )
!  Reaction Label HET_SO2_HNO4    
             RKI( NCELL,  919) =  BLKHET(  NCELL, IK_HETERO_HNO4 )
!  Reaction Label HET_HMS_L1      
             RKI( NCELL,  920) =  BLKHET(  NCELL, IK_HETERO_HMSL1 )
!  Reaction Label HET_HMS_L2      
             RKI( NCELL,  921) =  BLKHET(  NCELL, IK_HETERO_HMSL2 )

        END DO  
!  Multiply rate constants by [M], [O2], [N2], [H2O], [H2], or [CH4]
!  where needed and return
       IF ( NWM .GT. 0 ) THEN
          DO NRT = 1, NWM
             IRXN = NRXWM( NRT )
             DO NCELL = 1, NUMCELLS
                RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_AIR
             END DO
          END DO
       END IF
       IF ( NWO2 .GT. 0 ) THEN
          DO NRT = 1, NWO2
             IRXN = NRXWO2( NRT )
             DO NCELL = 1, NUMCELLS
                RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_O2
             END DO
          END DO
       END IF
       IF ( NWN2 .GT. 0 ) THEN
          DO NRT = 1, NWN2
             IRXN = NRXWN2( NRT )
             DO NCELL = 1, NUMCELLS
                RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_N2
             END DO
          END DO
       END IF
       IF ( NWW .GT. 0 ) THEN
          DO NRT = 1, NWW
             IRXN = NRXWW( NRT )
             DO NCELL = 1, NUMCELLS
                RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * BLKH2O( NCELL )
             END DO
          END DO
       END IF
       IF ( NWH2 .GT. 0 ) THEN
          DO NRT = 1, NWH2
             IRXN = NRXWH2( NRT )
             DO NCELL = 1, NUMCELLS
                RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_H2
             END DO
          END DO
       END IF
       IF ( NWCH4 .GT. 0 ) THEN
          DO NRT = 1, NWCH4
             IRXN = NRXWCH4( NRT )
             DO NCELL = 1, NUMCELLS
                RKI( NCELL,IRXN ) = RKI( NCELL,IRXN ) * ATM_CH4
             END DO
          END DO
       END IF
       RETURN
       END SUBROUTINE CALC_RCONST
         FUNCTION MAP_CHEMISTRY_SPECIES() RESULT ( SUCCESS )

! Purpose find or test the CGRID Index, Species Type, and Conversion Factor
! for the Mechanism against the CMAQ namelists

            USE UTILIO_DEFN
            USE CGRID_SPCS
            USE RXNS_DATA

            IMPLICIT NONE

!Parameters:
            CHARACTER(  1 ), PARAMETER :: BL = ' '
            INTEGER,         PARAMETER :: SPC_DIM = 200
!Local:

            LOGICAL SUCCESS
            INTEGER I, IOS, J
            INTEGER I1, I2, I3, I4      ! SURROGATE TYPE 1 COUNTERS
            INTEGER J1, J2              ! SURROGATE TYPE 2 COUNTERS
            INTEGER K1, K2, K3, K4, K5  ! CONTROL TYPE COUNTERS
            INTEGER ICALL

            LOGICAL :: ORDER = .TRUE.
            LOGICAL :: FOUND = .TRUE.

            CHARACTER( 120 ) :: XMSG

            CHARACTER( 16 ), ALLOCATABLE     :: CGRID_SPC  ( : )
            CHARACTER( 16 ), ALLOCATABLE     :: NML_SPC    ( : )
            CHARACTER(  2 ), ALLOCATABLE     :: NML_TYPE   ( : )
            INTEGER,         ALLOCATABLE     :: NML_INDEX  ( : )
            LOGICAL,         ALLOCATABLE     :: NML_CONVERT( : )
            REAL,            ALLOCATABLE     :: NML_MOLWT  ( : )
            REAL                             :: DELTA            ! fractional difference

            LOGICAL, SAVE :: INITIALIZED = .FALSE.

            IF( INITIALIZED )RETURN

            INITIALIZED = .TRUE.
            SUCCESS     = .TRUE.


            ALLOCATE ( CGRID_SPC( NSPCSD - 1 ),    &
     &                 NML_SPC  ( NSPCSD - 1 ),    &
     &                 NML_INDEX( NSPCSD - 1 ),    &
     &                 NML_TYPE( NSPCSD - 1 ),     &
     &                 NML_CONVERT( NSPCSD - 1 ),  &
     &                 NML_MOLWT( NSPCSD - 1 ),    &
     &                 STAT = IOS )


            J = 0


            NML_INDEX     = -1
            TYPE_INDEX    = -1
            NML_TYPE      = '??'
            NML_CONVERT   = .FALSE.

            DO I = 1, N_GC_SPC ! load gc names and indices
               J = J + 1
               CGRID_SPC( I )     = GC_SPC( I )
               NML_INDEX( J )     = I + GC_STRT -1
               NML_TYPE( J )      = 'GC'
               NML_MOLWT( J )     = GC_MOLWT( I )
            END DO

            DO I = 1, N_AE_SPC ! load ae names and indices
               J = J + 1
               CGRID_SPC( J )     = AE_SPC( I )
               NML_INDEX( J )     = I + AE_STRT - 1
               NML_TYPE( J )      = 'AE'
               NML_CONVERT( J )   = .TRUE.
               NML_MOLWT( J )     = AE_MOLWT( I )
            END DO

            DO I = 1, N_NR_SPC ! load nr names and indices
               J = J + 1
               CGRID_SPC( J )     = NR_SPC( I )
               NML_INDEX( J )     = I + NR_STRT - 1
               NML_TYPE( J )      = 'NR'
               NML_MOLWT( J )     = NR_MOLWT( I )
            END DO

            DO I = 1, N_TR_SPC ! load tr names and indices
               J = J + 1
               CGRID_SPC( J )     = TR_SPC( I )
               NML_INDEX( J )     = I + TR_STRT - 1
               NML_TYPE( J )      = 'TR'
               NML_MOLWT( J )     = TR_MOLWT( I )
            END DO

            NML_SPC( 1:(NSPCSD-1) ) = CGRID_SPC( 1:(NSPCSD-1) )


! determine if mechanism species are in cgrid species

            DO I = 1, NUMB_MECH_SPC
! set species informations arrays using SPECIES_LIST array before mapping
               CHEMISTRY_SPC( I ) = SPECIES_LIST( I )%CHEMISTRY_SPC
               CGRID_INDEX  ( I ) = SPECIES_LIST( I )%CGRID_INDEX
               SPECIES_TYPE ( I ) = SPECIES_LIST( I )%SPECIES_TYPE
               CONVERT_CONC ( I ) = SPECIES_LIST( I )%CONVERT_CONC
               SPECIES_MOLWT( I ) = SPECIES_LIST( I )%SPECIES_MOLWT

               I1 = INDEX1R( CHEMISTRY_SPC( I ), (NSPCSD-1), CGRID_SPC )
               IF ( I1 .LT. 1 ) THEN
                  FOUND = .FALSE.
               ELSE
                  FOUND = .TRUE.
                  IF( .NOT. MAPPED_TO_CGRID )THEN
                      CGRID_INDEX( I )   = NML_INDEX( I1 )
                      SPECIES_TYPE( I )  = NML_TYPE ( I1 )
                      SPECIES_MOLWT( I ) = NML_MOLWT( I1 )
                      CONVERT_CONC( I )  = NML_CONVERT( I1 )
                  ELSE
                      IF(CGRID_INDEX( I ) .NE. NML_INDEX( I1 ))THEN
                         SUCCESS = .FALSE.
                         XMSG = '*** For Species ' // TRIM( CHEMISTRY_SPC( I ) ) &
    &                        // ' cgrid index does not match mechanism value.'
                         WRITE( LOGDEV,'( /5X, A )' ) TRIM( XMSG )
                         WRITE( XMSG,'(A,I3,1X,I3)')'CGRID Indices: Mechanism and NML Values are ',    &
    &                    CGRID_INDEX( I ),NML_INDEX( I1 )
                         WRITE( LOGDEV,'( 5X, A )' )XMSG
                      END IF
                      IF(CONVERT_CONC( I ) .NEQV. NML_CONVERT( I1 ))THEN
                         SUCCESS = .FALSE.
                         XMSG = '*** For Species ' // TRIM( CHEMISTRY_SPC( I ) ) &
    &                        // ' species unit conversion flag does not match mechanism value.'
                         WRITE( LOGDEV,'( /5X, A )' ) TRIM( XMSG )
                         WRITE( XMSG,'(A,1X,L21X,L2)')'CONVERSION FLAGS: Mechanism and NML Values are ', &
    &                    CONVERT_CONC( I ),NML_CONVERT( I1 )
                         WRITE( LOGDEV,'( 5X, A )' )XMSG
                         WRITE( XMSG,'(A,1X,A3,1X,A3)')'SPECIES TYPE: Mechanism and NML Values are ',    &
    &                    SPECIES_TYPE( I ),NML_TYPE( I1 )
                         WRITE( LOGDEV,'( 5X, A )' )XMSG
                      END IF
                      DELTA = ( SPECIES_MOLWT( I ) - NML_MOLWT( I1 ) )/MAX(NML_MOLWT( I1 ),1.0E-20)
                      IF( ABS( DELTA ) .GE. 0.05 )THEN
                         IF( CONVERT_CONC( I ) )SUCCESS = .FALSE.
                         XMSG = '*** For Species ' // TRIM( CHEMISTRY_SPC( I ) ) &
    &                        // ' species molecular weight does not match mechanism value.'
                         WRITE( LOGDEV,'( /5X, A )' ) TRIM( XMSG )
                         WRITE( XMSG,'(A,2(ES12.4,1X))')'Molecular Weight: Mechanism and NML Values are ', &
    &                    SPECIES_MOLWT( I ), NML_MOLWT( I1 )
                         WRITE( LOGDEV,'( 5X, A )' )XMSG
                      END IF
                 END IF
              END IF
              IF( INDEX( CHEMISTRY_SPC( I ), 'SRF') .GT. 0 )THEN
                  SUCCESS = .FALSE.
                  XMSG = '*** reactions cannot use modal aerosol surface area as species'
                  WRITE( LOGDEV,'( /5X, A )' ) TRIM( XMSG )
                  XMSG = TRIM( CHEMISTRY_SPC( I ) )
                  WRITE( LOGDEV,'( 2X, A )' ) TRIM( XMSG )
              END IF
              IF( INDEX( CHEMISTRY_SPC( I ), 'NUM') .GT. 0 )THEN
                  SUCCESS = .FALSE.
                  XMSG = '*** reactions cannot use modal aerosol number density as species'
                  WRITE( LOGDEV,'( /5X, A )' ) TRIM( XMSG )
                  XMSG = TRIM( CHEMISTRY_SPC( I ) )
                  WRITE( LOGDEV,'( 2X, A )' ) TRIM( XMSG )
              END IF
              IF ( .NOT. FOUND ) THEN
                 XMSG = 'Fatal error: Mechanism Species found not in species namelist:'
                 WRITE( LOGDEV,'( /5X, A )', ADVANCE = 'NO' ) TRIM( XMSG )
                 XMSG = TRIM( CHEMISTRY_SPC( I ) )
                 WRITE( LOGDEV,'( 2X, A )' ) TRIM( XMSG )
                 SUCCESS = .FALSE.
              END IF
            END DO

            IF( SUCCESS )RETURN

            WRITE(LOGDEV,99901)TRIM( MECHNAME )
            XMSG = 'The FATAL errors found in namelist used. Check ' &
      &          //  'the log of exiting processor if more details are needed.'
            CALL M3WARN('MAP_CHEMISTRY_SPECIES',0,0,XMSG)


99901       FORMAT( / 'FATAL error(s) found in the namelists used. Check that ' &
     &     /  'these namelists contain the above data as the respective files ' &
     &     /  'in the respository version of the mechanism: ' , A )

         RETURN

         END FUNCTION MAP_CHEMISTRY_SPECIES
!----------------------------------------------------------------------------------------
         INTEGER FUNCTION INDEX1R ( NAME, N, NLIST )
            IMPLICIT NONE
            CHARACTER( * ) NAME        ! character string being searched for
            INTEGER N                  ! length of array to be searched
            CHARACTER( * ) NLIST( : )  ! array to be searched

            INTEGER I

            DO I = 1, N
               IF ( NAME .EQ. NLIST( I ) ) THEN
                  INDEX1R = I
                  RETURN
               END IF
           END DO
           INDEX1R = 0
           RETURN

          END FUNCTION INDEX1R
          SUBROUTINE RESET_SPECIES_POINTERS( IOLD2NEW )

             USE RXNS_DATA
             IMPLICIT NONE
             INTEGER, INTENT( IN ) :: IOLD2NEW( :,: ) 


             INDEX_O3          = IOLD2NEW( INDEX_O3         , 1 )
             INDEX_O3P         = IOLD2NEW( INDEX_O3P        , 1 )
             INDEX_O1D         = IOLD2NEW( INDEX_O1D        , 1 )
             INDEX_H2O2        = IOLD2NEW( INDEX_H2O2       , 1 )
             INDEX_HO          = IOLD2NEW( INDEX_HO         , 1 )
             INDEX_NO2         = IOLD2NEW( INDEX_NO2        , 1 )
             INDEX_NO          = IOLD2NEW( INDEX_NO         , 1 )
             INDEX_NO3         = IOLD2NEW( INDEX_NO3        , 1 )
             INDEX_HONO        = IOLD2NEW( INDEX_HONO       , 1 )
             INDEX_HNO3        = IOLD2NEW( INDEX_HNO3       , 1 )
             INDEX_HNO4        = IOLD2NEW( INDEX_HNO4       , 1 )
             INDEX_HO2         = IOLD2NEW( INDEX_HO2        , 1 )
             INDEX_HCHO        = IOLD2NEW( INDEX_HCHO       , 1 )
             INDEX_CO          = IOLD2NEW( INDEX_CO         , 1 )
             INDEX_ACD         = IOLD2NEW( INDEX_ACD        , 1 )
             INDEX_MO2         = IOLD2NEW( INDEX_MO2        , 1 )
             INDEX_ALD         = IOLD2NEW( INDEX_ALD        , 1 )
             INDEX_ETHP        = IOLD2NEW( INDEX_ETHP       , 1 )
             INDEX_ACT         = IOLD2NEW( INDEX_ACT        , 1 )
             INDEX_ACO3        = IOLD2NEW( INDEX_ACO3       , 1 )
             INDEX_UALD        = IOLD2NEW( INDEX_UALD       , 1 )
             INDEX_KET         = IOLD2NEW( INDEX_KET        , 1 )
             INDEX_PINAL       = IOLD2NEW( INDEX_PINAL      , 1 )
             INDEX_PINALP      = IOLD2NEW( INDEX_PINALP     , 1 )
             INDEX_LIMAL       = IOLD2NEW( INDEX_LIMAL      , 1 )
             INDEX_LIMALP      = IOLD2NEW( INDEX_LIMALP     , 1 )
             INDEX_MEK         = IOLD2NEW( INDEX_MEK        , 1 )
             INDEX_HKET        = IOLD2NEW( INDEX_HKET       , 1 )
             INDEX_MACR        = IOLD2NEW( INDEX_MACR       , 1 )
             INDEX_MACP        = IOLD2NEW( INDEX_MACP       , 1 )
             INDEX_XO2         = IOLD2NEW( INDEX_XO2        , 1 )
             INDEX_MVK         = IOLD2NEW( INDEX_MVK        , 1 )
             INDEX_GLY         = IOLD2NEW( INDEX_GLY        , 1 )
             INDEX_MGLY        = IOLD2NEW( INDEX_MGLY       , 1 )
             INDEX_DCB1        = IOLD2NEW( INDEX_DCB1       , 1 )
             INDEX_DCB2        = IOLD2NEW( INDEX_DCB2       , 1 )
             INDEX_BALD        = IOLD2NEW( INDEX_BALD       , 1 )
             INDEX_BEN         = IOLD2NEW( INDEX_BEN        , 1 )
             INDEX_BAL1        = IOLD2NEW( INDEX_BAL1       , 1 )
             INDEX_OP1         = IOLD2NEW( INDEX_OP1        , 1 )
             INDEX_OP2         = IOLD2NEW( INDEX_OP2        , 1 )
             INDEX_OPB         = IOLD2NEW( INDEX_OPB        , 1 )
             INDEX_VOP3        = IOLD2NEW( INDEX_VOP3       , 1 )
             INDEX_PAA         = IOLD2NEW( INDEX_PAA        , 1 )
             INDEX_ONIT        = IOLD2NEW( INDEX_ONIT       , 1 )
             INDEX_PAN         = IOLD2NEW( INDEX_PAN        , 1 )
             INDEX_CO2         = IOLD2NEW( INDEX_CO2        , 1 )
             INDEX_PPN         = IOLD2NEW( INDEX_PPN        , 1 )
             INDEX_RCO3        = IOLD2NEW( INDEX_RCO3       , 1 )
             INDEX_HC3P        = IOLD2NEW( INDEX_HC3P       , 1 )
             INDEX_VTRPN       = IOLD2NEW( INDEX_VTRPN      , 1 )
             INDEX_VHONIT      = IOLD2NEW( INDEX_VHONIT     , 1 )
             INDEX_N2O5        = IOLD2NEW( INDEX_N2O5       , 1 )
             INDEX_SO2         = IOLD2NEW( INDEX_SO2        , 1 )
             INDEX_SULF        = IOLD2NEW( INDEX_SULF       , 1 )
             INDEX_SULRXN      = IOLD2NEW( INDEX_SULRXN     , 1 )
             INDEX_ETH         = IOLD2NEW( INDEX_ETH        , 1 )
             INDEX_HC3         = IOLD2NEW( INDEX_HC3        , 1 )
             INDEX_ASOATJ      = IOLD2NEW( INDEX_ASOATJ     , 1 )
             INDEX_HC5         = IOLD2NEW( INDEX_HC5        , 1 )
             INDEX_HC5P        = IOLD2NEW( INDEX_HC5P       , 1 )
             INDEX_ETE         = IOLD2NEW( INDEX_ETE        , 1 )
             INDEX_ETEP        = IOLD2NEW( INDEX_ETEP       , 1 )
             INDEX_OLT         = IOLD2NEW( INDEX_OLT        , 1 )
             INDEX_OLTP        = IOLD2NEW( INDEX_OLTP       , 1 )
             INDEX_OLI         = IOLD2NEW( INDEX_OLI        , 1 )
             INDEX_OLIP        = IOLD2NEW( INDEX_OLIP       , 1 )
             INDEX_ACE         = IOLD2NEW( INDEX_ACE        , 1 )
             INDEX_ORA1        = IOLD2NEW( INDEX_ORA1       , 1 )
             INDEX_BENP        = IOLD2NEW( INDEX_BENP       , 1 )
             INDEX_PHEN        = IOLD2NEW( INDEX_PHEN       , 1 )
             INDEX_TOL         = IOLD2NEW( INDEX_TOL        , 1 )
             INDEX_TOLP        = IOLD2NEW( INDEX_TOLP       , 1 )
             INDEX_CSL         = IOLD2NEW( INDEX_CSL        , 1 )
             INDEX_XYL         = IOLD2NEW( INDEX_XYL        , 1 )
             INDEX_XYLP        = IOLD2NEW( INDEX_XYLP       , 1 )
             INDEX_EBZ         = IOLD2NEW( INDEX_EBZ        , 1 )
             INDEX_EBZP        = IOLD2NEW( INDEX_EBZP       , 1 )
             INDEX_ISO         = IOLD2NEW( INDEX_ISO        , 1 )
             INDEX_ISON        = IOLD2NEW( INDEX_ISON       , 1 )
             INDEX_ISONP       = IOLD2NEW( INDEX_ISONP      , 1 )
             INDEX_ISOP        = IOLD2NEW( INDEX_ISOP       , 1 )
             INDEX_ISHP        = IOLD2NEW( INDEX_ISHP       , 1 )
             INDEX_IEPOX       = IOLD2NEW( INDEX_IEPOX      , 1 )
             INDEX_IPX         = IOLD2NEW( INDEX_IPX        , 1 )
             INDEX_INALD       = IOLD2NEW( INDEX_INALD      , 1 )
             INDEX_ROH         = IOLD2NEW( INDEX_ROH        , 1 )
             INDEX_API         = IOLD2NEW( INDEX_API        , 1 )
             INDEX_APIP1       = IOLD2NEW( INDEX_APIP1      , 1 )
             INDEX_APIP2       = IOLD2NEW( INDEX_APIP2      , 1 )
             INDEX_LIM         = IOLD2NEW( INDEX_LIM        , 1 )
             INDEX_LIMP1       = IOLD2NEW( INDEX_LIMP1      , 1 )
             INDEX_LIMP2       = IOLD2NEW( INDEX_LIMP2      , 1 )
             INDEX_ACTP        = IOLD2NEW( INDEX_ACTP       , 1 )
             INDEX_MEKP        = IOLD2NEW( INDEX_MEKP       , 1 )
             INDEX_KETP        = IOLD2NEW( INDEX_KETP       , 1 )
             INDEX_MCP         = IOLD2NEW( INDEX_MCP        , 1 )
             INDEX_MVKP        = IOLD2NEW( INDEX_MVKP       , 1 )
             INDEX_UALP        = IOLD2NEW( INDEX_UALP       , 1 )
             INDEX_DCB3        = IOLD2NEW( INDEX_DCB3       , 1 )
             INDEX_BALP        = IOLD2NEW( INDEX_BALP       , 1 )
             INDEX_ADDC        = IOLD2NEW( INDEX_ADDC       , 1 )
             INDEX_CHO         = IOLD2NEW( INDEX_CHO        , 1 )
             INDEX_MCT         = IOLD2NEW( INDEX_MCT        , 1 )
             INDEX_MCTO        = IOLD2NEW( INDEX_MCTO       , 1 )
             INDEX_MOH         = IOLD2NEW( INDEX_MOH        , 1 )
             INDEX_EOH         = IOLD2NEW( INDEX_EOH        , 1 )
             INDEX_ETEG        = IOLD2NEW( INDEX_ETEG       , 1 )
             INDEX_HC10P       = IOLD2NEW( INDEX_HC10P      , 1 )
             INDEX_MAHP        = IOLD2NEW( INDEX_MAHP       , 1 )
             INDEX_ORA2        = IOLD2NEW( INDEX_ORA2       , 1 )
             INDEX_ORAP        = IOLD2NEW( INDEX_ORAP       , 1 )
             INDEX_MPAN        = IOLD2NEW( INDEX_MPAN       , 1 )
             INDEX_MCTP        = IOLD2NEW( INDEX_MCTP       , 1 )
             INDEX_OLNN        = IOLD2NEW( INDEX_OLNN       , 1 )
             INDEX_OLND        = IOLD2NEW( INDEX_OLND       , 1 )
             INDEX_APINP1      = IOLD2NEW( INDEX_APINP1     , 1 )
             INDEX_APINP2      = IOLD2NEW( INDEX_APINP2     , 1 )
             INDEX_LIMNP1      = IOLD2NEW( INDEX_LIMNP1     , 1 )
             INDEX_LIMNP2      = IOLD2NEW( INDEX_LIMNP2     , 1 )
             INDEX_ADCN        = IOLD2NEW( INDEX_ADCN       , 1 )
             INDEX_VHOM        = IOLD2NEW( INDEX_VHOM       , 1 )
             INDEX_VROCP4OXY2  = IOLD2NEW( INDEX_VROCP4OXY2 , 1 )
             INDEX_VROCN1OXY6  = IOLD2NEW( INDEX_VROCN1OXY6 , 1 )
             INDEX_FURANONE    = IOLD2NEW( INDEX_FURANONE   , 1 )
             INDEX_VROCP3OXY2  = IOLD2NEW( INDEX_VROCP3OXY2 , 1 )
             INDEX_VROCP0OXY4  = IOLD2NEW( INDEX_VROCP0OXY4 , 1 )
             INDEX_BAL2        = IOLD2NEW( INDEX_BAL2       , 1 )
             INDEX_VELHOM      = IOLD2NEW( INDEX_VELHOM     , 1 )
             INDEX_VROCIOXY    = IOLD2NEW( INDEX_VROCIOXY   , 1 )
             INDEX_SLOWROC     = IOLD2NEW( INDEX_SLOWROC    , 1 )
             INDEX_ACRO        = IOLD2NEW( INDEX_ACRO       , 1 )
             INDEX_BDE13       = IOLD2NEW( INDEX_BDE13      , 1 )
             INDEX_BDE13P      = IOLD2NEW( INDEX_BDE13P     , 1 )
             INDEX_FURAN       = IOLD2NEW( INDEX_FURAN      , 1 )
             INDEX_FURANO2     = IOLD2NEW( INDEX_FURANO2    , 1 )
             INDEX_PROG        = IOLD2NEW( INDEX_PROG       , 1 )
             INDEX_SESQ        = IOLD2NEW( INDEX_SESQ       , 1 )
             INDEX_SESQNRO2    = IOLD2NEW( INDEX_SESQNRO2   , 1 )
             INDEX_VROCN2OXY2  = IOLD2NEW( INDEX_VROCN2OXY2 , 1 )
             INDEX_SESQRO2     = IOLD2NEW( INDEX_SESQRO2    , 1 )
             INDEX_VROCP0OXY2  = IOLD2NEW( INDEX_VROCP0OXY2 , 1 )
             INDEX_VROCP1OXY3  = IOLD2NEW( INDEX_VROCP1OXY3 , 1 )
             INDEX_AGLYOLIGJ   = IOLD2NEW( INDEX_AGLYOLIGJ  , 1 )
             INDEX_IEPOXP      = IOLD2NEW( INDEX_IEPOXP     , 1 )
             INDEX_AISO3NOSJ   = IOLD2NEW( INDEX_AISO3NOSJ  , 1 )
             INDEX_ASO4J       = IOLD2NEW( INDEX_ASO4J      , 1 )
             INDEX_AISO3OSJ    = IOLD2NEW( INDEX_AISO3OSJ   , 1 )
             INDEX_AISO4J      = IOLD2NEW( INDEX_AISO4J     , 1 )
             INDEX_AISO5J      = IOLD2NEW( INDEX_AISO5J     , 1 )
             INDEX_VROCP6ALK   = IOLD2NEW( INDEX_VROCP6ALK  , 1 )
             INDEX_VROCP6ALKP  = IOLD2NEW( INDEX_VROCP6ALKP , 1 )
             INDEX_VROCP5ALK   = IOLD2NEW( INDEX_VROCP5ALK  , 1 )
             INDEX_VROCP5ALKP  = IOLD2NEW( INDEX_VROCP5ALKP , 1 )
             INDEX_VROCP4ALK   = IOLD2NEW( INDEX_VROCP4ALK  , 1 )
             INDEX_VROCP4ALKP  = IOLD2NEW( INDEX_VROCP4ALKP , 1 )
             INDEX_VROCP3ALK   = IOLD2NEW( INDEX_VROCP3ALK  , 1 )
             INDEX_VROCP3ALKP  = IOLD2NEW( INDEX_VROCP3ALKP , 1 )
             INDEX_VROCP2ALK   = IOLD2NEW( INDEX_VROCP2ALK  , 1 )
             INDEX_VROCP2ALKP  = IOLD2NEW( INDEX_VROCP2ALKP , 1 )
             INDEX_VROCP1ALK   = IOLD2NEW( INDEX_VROCP1ALK  , 1 )
             INDEX_VROCP1ALKP  = IOLD2NEW( INDEX_VROCP1ALKP , 1 )
             INDEX_HC10        = IOLD2NEW( INDEX_HC10       , 1 )
             INDEX_VROCP6ALKP2 = IOLD2NEW( INDEX_VROCP6ALKP2, 1 )
             INDEX_VROCP5ALKP2 = IOLD2NEW( INDEX_VROCP5ALKP2, 1 )
             INDEX_VROCP4ALKP2 = IOLD2NEW( INDEX_VROCP4ALKP2, 1 )
             INDEX_VROCP2OXY2  = IOLD2NEW( INDEX_VROCP2OXY2 , 1 )
             INDEX_VROCP3ALKP2 = IOLD2NEW( INDEX_VROCP3ALKP2, 1 )
             INDEX_VROCP1OXY1  = IOLD2NEW( INDEX_VROCP1OXY1 , 1 )
             INDEX_VROCP2ALKP2 = IOLD2NEW( INDEX_VROCP2ALKP2, 1 )
             INDEX_VROCP1ALKP2 = IOLD2NEW( INDEX_VROCP1ALKP2, 1 )
             INDEX_VROCN1OXY1  = IOLD2NEW( INDEX_VROCN1OXY1 , 1 )
             INDEX_HC10P2      = IOLD2NEW( INDEX_HC10P2     , 1 )
             INDEX_VROCP6ARO   = IOLD2NEW( INDEX_VROCP6ARO  , 1 )
             INDEX_VROCP6AROP  = IOLD2NEW( INDEX_VROCP6AROP , 1 )
             INDEX_VROCN2OXY4  = IOLD2NEW( INDEX_VROCN2OXY4 , 1 )
             INDEX_VROCN1OXY3  = IOLD2NEW( INDEX_VROCN1OXY3 , 1 )
             INDEX_VROCP5ARO   = IOLD2NEW( INDEX_VROCP5ARO  , 1 )
             INDEX_VROCP5AROP  = IOLD2NEW( INDEX_VROCP5AROP , 1 )
             INDEX_NAPH        = IOLD2NEW( INDEX_NAPH       , 1 )
             INDEX_NAPHP       = IOLD2NEW( INDEX_NAPHP      , 1 )
             INDEX_VROCN2OXY8  = IOLD2NEW( INDEX_VROCN2OXY8 , 1 )
             INDEX_VROCP5OXY1  = IOLD2NEW( INDEX_VROCP5OXY1 , 1 )
             INDEX_VROCP6OXY1  = IOLD2NEW( INDEX_VROCP6OXY1 , 1 )
             INDEX_ECH4        = IOLD2NEW( INDEX_ECH4       , 1 )
             INDEX_ATRPNJ      = IOLD2NEW( INDEX_ATRPNJ     , 1 )
             INDEX_AHOMJ       = IOLD2NEW( INDEX_AHOMJ      , 1 )
             INDEX_AHONITJ     = IOLD2NEW( INDEX_AHONITJ    , 1 )
             INDEX_STY         = IOLD2NEW( INDEX_STY        , 1 )
             INDEX_STYP        = IOLD2NEW( INDEX_STYP       , 1 )
             INDEX_ANO3I       = IOLD2NEW( INDEX_ANO3I      , 1 )
             INDEX_ANO3J       = IOLD2NEW( INDEX_ANO3J      , 1 )
             INDEX_AMTN1J      = IOLD2NEW( INDEX_AMTN1J     , 1 )
             INDEX_CL2         = IOLD2NEW( INDEX_CL2        , 1 )
             INDEX_CL          = IOLD2NEW( INDEX_CL         , 1 )
             INDEX_CLO         = IOLD2NEW( INDEX_CLO        , 1 )
             INDEX_OCLO        = IOLD2NEW( INDEX_OCLO       , 1 )
             INDEX_CL2O2       = IOLD2NEW( INDEX_CL2O2      , 1 )
             INDEX_CLOO        = IOLD2NEW( INDEX_CLOO       , 1 )
             INDEX_HOCL        = IOLD2NEW( INDEX_HOCL       , 1 )
             INDEX_CLNO        = IOLD2NEW( INDEX_CLNO       , 1 )
             INDEX_CLNO2       = IOLD2NEW( INDEX_CLNO2      , 1 )
             INDEX_CLNO3       = IOLD2NEW( INDEX_CLNO3      , 1 )
             INDEX_HCOCL       = IOLD2NEW( INDEX_HCOCL      , 1 )
             INDEX_HCL         = IOLD2NEW( INDEX_HCL        , 1 )
             INDEX_BR2         = IOLD2NEW( INDEX_BR2        , 1 )
             INDEX_BR          = IOLD2NEW( INDEX_BR         , 1 )
             INDEX_BRO         = IOLD2NEW( INDEX_BRO        , 1 )
             INDEX_OBRO        = IOLD2NEW( INDEX_OBRO       , 1 )
             INDEX_HOBR        = IOLD2NEW( INDEX_HOBR       , 1 )
             INDEX_BRNO        = IOLD2NEW( INDEX_BRNO       , 1 )
             INDEX_BRNO2       = IOLD2NEW( INDEX_BRNO2      , 1 )
             INDEX_BRNO3       = IOLD2NEW( INDEX_BRNO3      , 1 )
             INDEX_CH2BR2      = IOLD2NEW( INDEX_CH2BR2     , 1 )
             INDEX_CHBR3       = IOLD2NEW( INDEX_CHBR3      , 1 )
             INDEX_HCOBR       = IOLD2NEW( INDEX_HCOBR      , 1 )
             INDEX_HBR         = IOLD2NEW( INDEX_HBR        , 1 )
             INDEX_I2          = IOLD2NEW( INDEX_I2         , 1 )
             INDEX_I           = IOLD2NEW( INDEX_I          , 1 )
             INDEX_IO          = IOLD2NEW( INDEX_IO         , 1 )
             INDEX_OIO         = IOLD2NEW( INDEX_OIO        , 1 )
             INDEX_I2O2        = IOLD2NEW( INDEX_I2O2       , 1 )
             INDEX_HOI         = IOLD2NEW( INDEX_HOI        , 1 )
             INDEX_HI          = IOLD2NEW( INDEX_HI         , 1 )
             INDEX_INO         = IOLD2NEW( INDEX_INO        , 1 )
             INDEX_INO2        = IOLD2NEW( INDEX_INO2       , 1 )
             INDEX_INO3        = IOLD2NEW( INDEX_INO3       , 1 )
             INDEX_CH3I        = IOLD2NEW( INDEX_CH3I       , 1 )
             INDEX_CH2I2       = IOLD2NEW( INDEX_CH2I2      , 1 )
             INDEX_I2O3        = IOLD2NEW( INDEX_I2O3       , 1 )
             INDEX_I2O4        = IOLD2NEW( INDEX_I2O4       , 1 )
             INDEX_BRCL        = IOLD2NEW( INDEX_BRCL       , 1 )
             INDEX_ICL         = IOLD2NEW( INDEX_ICL        , 1 )
             INDEX_IBR         = IOLD2NEW( INDEX_IBR        , 1 )
             INDEX_CH2IBR      = IOLD2NEW( INDEX_CH2IBR     , 1 )
             INDEX_CH2ICL      = IOLD2NEW( INDEX_CH2ICL     , 1 )
             INDEX_CHBR2CL     = IOLD2NEW( INDEX_CHBR2CL    , 1 )
             INDEX_CHBRCL2     = IOLD2NEW( INDEX_CHBRCL2    , 1 )
             INDEX_CH2BRCL     = IOLD2NEW( INDEX_CH2BRCL    , 1 )
             INDEX_NO2PIJ      = IOLD2NEW( INDEX_NO2PIJ     , 1 )
             INDEX_NO2PK       = IOLD2NEW( INDEX_NO2PK      , 1 )
             INDEX_ACLI        = IOLD2NEW( INDEX_ACLI       , 1 )
             INDEX_ACLJ        = IOLD2NEW( INDEX_ACLJ       , 1 )
             INDEX_ACLK        = IOLD2NEW( INDEX_ACLK       , 1 )
             INDEX_ABRJ        = IOLD2NEW( INDEX_ABRJ       , 1 )
             INDEX_AFEJ        = IOLD2NEW( INDEX_AFEJ       , 1 )
             INDEX_AMNJ        = IOLD2NEW( INDEX_AMNJ       , 1 )
             INDEX_AHMSJ       = IOLD2NEW( INDEX_AHMSJ      , 1 )
          END SUBROUTINE RESET_SPECIES_POINTERS
       END MODULE RXNS_FUNCTION
