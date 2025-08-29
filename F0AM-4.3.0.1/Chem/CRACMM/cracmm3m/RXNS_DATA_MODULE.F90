       MODULE RXNS_DATA


       IMPLICIT NONE



! --------- Photochemical Mechanism Reactions, Rates, etc. DAT ---------
! Source file: /work/MOD3DEV/has/2023cracmm_ages/runs/20250808f0amgen/UTIL/chemmech/input/cracmm3m/mech_cracmm3m.def
! for Mechanism Name: CRACMM3M                        

! This file is used to create mechanism data and functions

! The following are reserved symbols declared in this file:
!    MECHNAME        = Mechanism name
!    N_GAS_CHEM_SPC  = Total number of gas species in chemical mechanism
!    NUMB_MECH_SPC   = Total number of species in chemical mechanism
!    N_ACT_SP        = Number of active (determined by ODE solver) species in mechanism
!    GAS_CHEM_SPC    = Names of gas species in chemical mechanism
!    CHEMISTRY_SPC   = Names of species in chemical mechanism
!    CGRID_INDEX     = CGRID Index of species in chemical mechanism
!    SPECIES_TYPE    = Group or type of species 
!    SPECIES_MOLWT   = Molecular Weight of species (gm/mole)
!    NRXNS           = Number of mechanism reactions
!    ZERO_REACT_REACTIONS  = number zero reactant reactions
!    ONE_REACT_REACTIONS   = number one reactant reactions
!    TWO_REACT_REACTIONS   = number second order reactions
!    THREE_REACT_REACTIONS = number three reactant reactions
!    NSUNLIGHT_RXNS  = Number of mechanism reactions requiring sunlight
!    NTHERMAL_RXNS   = Number of mechanism reactions not requiring sunlight
!    KUNITS          = Units of mechanism reactions
!    KTYPE           = Reaction type
!    IRXBITS         = Bit test mask vector for selected reactions
!    IORDER          = Order of the reaction
!    NTERMS_JACOB    = Maximum number of nonzero terms in day/night Jacobian
!    MSTEPS_JACOB    = Maximum number of LU Decomposition steps to solve each Jacobian
!    KTN1            = Number of type 1 reactions
!    KRX1            = Reactions list pointer to type 1 reactions
!    KTN2            = Number of type 2 reactions
!    KRX2            = Reactions list pointer to type 2 reactions
!    KTN3            = Number of type 3 reactions
!    KRX3            = Reactions list pointer to type 3 reactions
!    KTN4            = Number of type 4 reactions
!    KRX4            = Reactions list pointer to type 4 reactions
!    KTN5            = Number of type 5 reactions
!    KRX5            = Reactions list pointer to type 5 reactions
!    KTN6            = Number of type 6 reactions
!    KRX6            = Reactions list pointer to type 6 reactions
!    KTN7            = Number of type 7 reactions
!    KRX7            = Reactions list pointer to type 7 reactions

!    NWM       = Number of air 3-body reactions
!    NRXWM     = Reactions list pointer to air 3-body reactions
!    ATM_AIR   = air 3-body reactions concentration
!    NWW       = Number of H2O 3-body reactions
!    NRXWW     = Reactions list pointer to H2O 3-body reactions
!    NWO2      = Number of reactions with O2
!    NRXWO2    = Reactions list pointer to O2 reactions
!    ATM_O2    = Oxygen reactions concentration
!    NWN2      = Number of N2 3-body reactions
!    NRXWN2    = Reactions list pointer to N2 3-body reactions
!    ATM_N2    = Nitrogen 3-body reactions concentration
!    NWCH4     = Number of reactions with CH4
!    NRXWCH4   = Reactions list pointer to CH4 reactions
!    ATM_CH4   = Methane reactions concentration
!    NWH2      = Number of reactions with H2
!    NRXWH2    = Reactions list pointer to H2 reactions
!    ATM_H2    = Hydrogen reactions concentration

!    MXPRD     = Maximum number of mechanism reaction products
!    IRR       = Reactions list pointer to reactants and products
!    RTDAT     = Kinetic reaction rates expressions components
!    NFALLOFFF = Number of falloff reactions
!    IRRFALL   = Reactions list pointer to falloff reactions
!    RFDAT     = Falloff reaction rates expressions components
!    SC        = Stoichiometric coefficients
!    NREACT    = Number of reactants in each mechanism reaction
!    NPRDCT    = Number of products in each mechanism reaction
!    RXLABEL   = Character label list for mechanism reactions
!    NMPHOT    = Number of mechanism photolytic reactions
!    NPHOTAB   = Number of photolytic reactions tables
!    IPH       = Reactions list pointer to photolytic reactions and tables
!    MHETERO   = Number of mechanism heteorogenous reactions
!    NHETERO   = Number of unique heteorogenous rate constants
!    IHETERO   = Reactions list pointer to heteorogenous reactions and tables

      CHARACTER( 32 ), PARAMETER :: MECHNAME = 'CRACMM3M'

      INTEGER, PARAMETER :: N_GAS_CHEM_SPC = 228
      INTEGER, PARAMETER :: NUMB_MECH_SPC  = 248

      CHARACTER( 16 ) :: GAS_CHEM_SPC( N_GAS_CHEM_SPC )
      CHARACTER( 16 ) :: CHEMISTRY_SPC( NUMB_MECH_SPC )
      CHARACTER( 16 ) :: SPECIES_TYPE(  NUMB_MECH_SPC )
      INTEGER         :: CGRID_INDEX (  NUMB_MECH_SPC )
      INTEGER         :: TYPE_INDEX  (  NUMB_MECH_SPC )
      REAL( 8 )       :: SPECIES_MOLWT( NUMB_MECH_SPC )
      LOGICAL         :: CONVERT_CONC(  NUMB_MECH_SPC )

! The below character and integer arrays list the model species names used in the 
! chemical mechanism. The gas species and their order should agree with 
! the GC_SPC array for the gas phase chemistry to work correctly. 
! If present, the CHEMISTRY_SPC names and species type should agree with the CGRID_SPCS module

      DATA GAS_CHEM_SPC(   1 ) / 'O3              ' /
      DATA GAS_CHEM_SPC(   2 ) / 'O3P             ' /
      DATA GAS_CHEM_SPC(   3 ) / 'O1D             ' /
      DATA GAS_CHEM_SPC(   4 ) / 'H2O2            ' /
      DATA GAS_CHEM_SPC(   5 ) / 'HO              ' /
      DATA GAS_CHEM_SPC(   6 ) / 'NO2             ' /
      DATA GAS_CHEM_SPC(   7 ) / 'NO              ' /
      DATA GAS_CHEM_SPC(   8 ) / 'NO3             ' /
      DATA GAS_CHEM_SPC(   9 ) / 'HONO            ' /
      DATA GAS_CHEM_SPC(  10 ) / 'HNO3            ' /
      DATA GAS_CHEM_SPC(  11 ) / 'HNO4            ' /
      DATA GAS_CHEM_SPC(  12 ) / 'HO2             ' /
      DATA GAS_CHEM_SPC(  13 ) / 'HCHO            ' /
      DATA GAS_CHEM_SPC(  14 ) / 'CO              ' /
      DATA GAS_CHEM_SPC(  15 ) / 'ACD             ' /
      DATA GAS_CHEM_SPC(  16 ) / 'MO2             ' /
      DATA GAS_CHEM_SPC(  17 ) / 'ALD             ' /
      DATA GAS_CHEM_SPC(  18 ) / 'ETHP            ' /
      DATA GAS_CHEM_SPC(  19 ) / 'ACT             ' /
      DATA GAS_CHEM_SPC(  20 ) / 'ACO3            ' /
      DATA GAS_CHEM_SPC(  21 ) / 'UALD            ' /
      DATA GAS_CHEM_SPC(  22 ) / 'KET             ' /
      DATA GAS_CHEM_SPC(  23 ) / 'PINAL           ' /
      DATA GAS_CHEM_SPC(  24 ) / 'PINALP          ' /
      DATA GAS_CHEM_SPC(  25 ) / 'LIMAL           ' /
      DATA GAS_CHEM_SPC(  26 ) / 'LIMALP          ' /
      DATA GAS_CHEM_SPC(  27 ) / 'MEK             ' /
      DATA GAS_CHEM_SPC(  28 ) / 'HKET            ' /
      DATA GAS_CHEM_SPC(  29 ) / 'MACR            ' /
      DATA GAS_CHEM_SPC(  30 ) / 'MACP            ' /
      DATA GAS_CHEM_SPC(  31 ) / 'XO2             ' /
      DATA GAS_CHEM_SPC(  32 ) / 'MVK             ' /
      DATA GAS_CHEM_SPC(  33 ) / 'GLY             ' /
      DATA GAS_CHEM_SPC(  34 ) / 'MGLY            ' /
      DATA GAS_CHEM_SPC(  35 ) / 'DCB1            ' /
      DATA GAS_CHEM_SPC(  36 ) / 'DCB2            ' /
      DATA GAS_CHEM_SPC(  37 ) / 'BALD            ' /
      DATA GAS_CHEM_SPC(  38 ) / 'BEN             ' /
      DATA GAS_CHEM_SPC(  39 ) / 'BAL1            ' /
      DATA GAS_CHEM_SPC(  40 ) / 'OP1             ' /
      DATA GAS_CHEM_SPC(  41 ) / 'OP2             ' /
      DATA GAS_CHEM_SPC(  42 ) / 'OPB             ' /
      DATA GAS_CHEM_SPC(  43 ) / 'VOP3            ' /
      DATA GAS_CHEM_SPC(  44 ) / 'PAA             ' /
      DATA GAS_CHEM_SPC(  45 ) / 'ONIT            ' /
      DATA GAS_CHEM_SPC(  46 ) / 'PAN             ' /
      DATA GAS_CHEM_SPC(  47 ) / 'CO2             ' /
      DATA GAS_CHEM_SPC(  48 ) / 'PPN             ' /
      DATA GAS_CHEM_SPC(  49 ) / 'RCO3            ' /
      DATA GAS_CHEM_SPC(  50 ) / 'HC3P            ' /
      DATA GAS_CHEM_SPC(  51 ) / 'VTRPN           ' /
      DATA GAS_CHEM_SPC(  52 ) / 'VHONIT          ' /
      DATA GAS_CHEM_SPC(  53 ) / 'N2O5            ' /
      DATA GAS_CHEM_SPC(  54 ) / 'SO2             ' /
      DATA GAS_CHEM_SPC(  55 ) / 'SULF            ' /
      DATA GAS_CHEM_SPC(  56 ) / 'SULRXN          ' /
      DATA GAS_CHEM_SPC(  57 ) / 'ETH             ' /
      DATA GAS_CHEM_SPC(  58 ) / 'HC3             ' /
      DATA GAS_CHEM_SPC(  59 ) / 'HC5             ' /
      DATA GAS_CHEM_SPC(  60 ) / 'HC5P            ' /
      DATA GAS_CHEM_SPC(  61 ) / 'ETE             ' /
      DATA GAS_CHEM_SPC(  62 ) / 'ETEP            ' /
      DATA GAS_CHEM_SPC(  63 ) / 'OLT             ' /
      DATA GAS_CHEM_SPC(  64 ) / 'OLTP            ' /
      DATA GAS_CHEM_SPC(  65 ) / 'OLI             ' /
      DATA GAS_CHEM_SPC(  66 ) / 'OLIP            ' /
      DATA GAS_CHEM_SPC(  67 ) / 'ACE             ' /
      DATA GAS_CHEM_SPC(  68 ) / 'ORA1            ' /
      DATA GAS_CHEM_SPC(  69 ) / 'BENP            ' /
      DATA GAS_CHEM_SPC(  70 ) / 'PHEN            ' /
      DATA GAS_CHEM_SPC(  71 ) / 'TOL             ' /
      DATA GAS_CHEM_SPC(  72 ) / 'TOLP            ' /
      DATA GAS_CHEM_SPC(  73 ) / 'CSL             ' /
      DATA GAS_CHEM_SPC(  74 ) / 'XYL             ' /
      DATA GAS_CHEM_SPC(  75 ) / 'XYLP            ' /
      DATA GAS_CHEM_SPC(  76 ) / 'EBZ             ' /
      DATA GAS_CHEM_SPC(  77 ) / 'EBZP            ' /
      DATA GAS_CHEM_SPC(  78 ) / 'ISO             ' /
      DATA GAS_CHEM_SPC(  79 ) / 'ISON            ' /
      DATA GAS_CHEM_SPC(  80 ) / 'ISONP           ' /
      DATA GAS_CHEM_SPC(  81 ) / 'ISOP            ' /
      DATA GAS_CHEM_SPC(  82 ) / 'ISHP            ' /
      DATA GAS_CHEM_SPC(  83 ) / 'IEPOX           ' /
      DATA GAS_CHEM_SPC(  84 ) / 'IPX             ' /
      DATA GAS_CHEM_SPC(  85 ) / 'INALD           ' /
      DATA GAS_CHEM_SPC(  86 ) / 'ROH             ' /
      DATA GAS_CHEM_SPC(  87 ) / 'API             ' /
      DATA GAS_CHEM_SPC(  88 ) / 'APIP1           ' /
      DATA GAS_CHEM_SPC(  89 ) / 'APIP2           ' /
      DATA GAS_CHEM_SPC(  90 ) / 'LIM             ' /
      DATA GAS_CHEM_SPC(  91 ) / 'LIMP1           ' /
      DATA GAS_CHEM_SPC(  92 ) / 'LIMP2           ' /
      DATA GAS_CHEM_SPC(  93 ) / 'ACTP            ' /
      DATA GAS_CHEM_SPC(  94 ) / 'MEKP            ' /
      DATA GAS_CHEM_SPC(  95 ) / 'KETP            ' /
      DATA GAS_CHEM_SPC(  96 ) / 'MCP             ' /
      DATA GAS_CHEM_SPC(  97 ) / 'MVKP            ' /
      DATA GAS_CHEM_SPC(  98 ) / 'UALP            ' /
      DATA GAS_CHEM_SPC(  99 ) / 'DCB3            ' /
      DATA GAS_CHEM_SPC( 100 ) / 'BALP            ' /
      DATA GAS_CHEM_SPC( 101 ) / 'ADDC            ' /
      DATA GAS_CHEM_SPC( 102 ) / 'CHO             ' /
      DATA GAS_CHEM_SPC( 103 ) / 'MCT             ' /
      DATA GAS_CHEM_SPC( 104 ) / 'MCTO            ' /
      DATA GAS_CHEM_SPC( 105 ) / 'MOH             ' /
      DATA GAS_CHEM_SPC( 106 ) / 'EOH             ' /
      DATA GAS_CHEM_SPC( 107 ) / 'ETEG            ' /
      DATA GAS_CHEM_SPC( 108 ) / 'HC10P           ' /
      DATA GAS_CHEM_SPC( 109 ) / 'MAHP            ' /
      DATA GAS_CHEM_SPC( 110 ) / 'ORA2            ' /
      DATA GAS_CHEM_SPC( 111 ) / 'ORAP            ' /
      DATA GAS_CHEM_SPC( 112 ) / 'MPAN            ' /
      DATA GAS_CHEM_SPC( 113 ) / 'MCTP            ' /
      DATA GAS_CHEM_SPC( 114 ) / 'OLNN            ' /
      DATA GAS_CHEM_SPC( 115 ) / 'OLND            ' /
      DATA GAS_CHEM_SPC( 116 ) / 'APINP1          ' /
      DATA GAS_CHEM_SPC( 117 ) / 'APINP2          ' /
      DATA GAS_CHEM_SPC( 118 ) / 'LIMNP1          ' /
      DATA GAS_CHEM_SPC( 119 ) / 'LIMNP2          ' /
      DATA GAS_CHEM_SPC( 120 ) / 'ADCN            ' /
      DATA GAS_CHEM_SPC( 121 ) / 'VHOM            ' /
      DATA GAS_CHEM_SPC( 122 ) / 'VROCP4OXY2      ' /
      DATA GAS_CHEM_SPC( 123 ) / 'VROCN1OXY6      ' /
      DATA GAS_CHEM_SPC( 124 ) / 'FURANONE        ' /
      DATA GAS_CHEM_SPC( 125 ) / 'VROCP3OXY2      ' /
      DATA GAS_CHEM_SPC( 126 ) / 'VROCP0OXY4      ' /
      DATA GAS_CHEM_SPC( 127 ) / 'BAL2            ' /
      DATA GAS_CHEM_SPC( 128 ) / 'VELHOM          ' /
      DATA GAS_CHEM_SPC( 129 ) / 'VROCIOXY        ' /
      DATA GAS_CHEM_SPC( 130 ) / 'SLOWROC         ' /
      DATA GAS_CHEM_SPC( 131 ) / 'ACRO            ' /
      DATA GAS_CHEM_SPC( 132 ) / 'BDE13           ' /
      DATA GAS_CHEM_SPC( 133 ) / 'BDE13P          ' /
      DATA GAS_CHEM_SPC( 134 ) / 'FURAN           ' /
      DATA GAS_CHEM_SPC( 135 ) / 'FURANO2         ' /
      DATA GAS_CHEM_SPC( 136 ) / 'PROG            ' /
      DATA GAS_CHEM_SPC( 137 ) / 'SESQ            ' /
      DATA GAS_CHEM_SPC( 138 ) / 'SESQNRO2        ' /
      DATA GAS_CHEM_SPC( 139 ) / 'VROCN2OXY2      ' /
      DATA GAS_CHEM_SPC( 140 ) / 'SESQRO2         ' /
      DATA GAS_CHEM_SPC( 141 ) / 'VROCP0OXY2      ' /
      DATA GAS_CHEM_SPC( 142 ) / 'VROCP1OXY3      ' /
      DATA GAS_CHEM_SPC( 143 ) / 'IEPOXP          ' /
      DATA GAS_CHEM_SPC( 144 ) / 'VROCP6ALK       ' /
      DATA GAS_CHEM_SPC( 145 ) / 'VROCP6ALKP      ' /
      DATA GAS_CHEM_SPC( 146 ) / 'VROCP5ALK       ' /
      DATA GAS_CHEM_SPC( 147 ) / 'VROCP5ALKP      ' /
      DATA GAS_CHEM_SPC( 148 ) / 'VROCP4ALK       ' /
      DATA GAS_CHEM_SPC( 149 ) / 'VROCP4ALKP      ' /
      DATA GAS_CHEM_SPC( 150 ) / 'VROCP3ALK       ' /
      DATA GAS_CHEM_SPC( 151 ) / 'VROCP3ALKP      ' /
      DATA GAS_CHEM_SPC( 152 ) / 'VROCP2ALK       ' /
      DATA GAS_CHEM_SPC( 153 ) / 'VROCP2ALKP      ' /
      DATA GAS_CHEM_SPC( 154 ) / 'VROCP1ALK       ' /
      DATA GAS_CHEM_SPC( 155 ) / 'VROCP1ALKP      ' /
      DATA GAS_CHEM_SPC( 156 ) / 'HC10            ' /
      DATA GAS_CHEM_SPC( 157 ) / 'VROCP6ALKP2     ' /
      DATA GAS_CHEM_SPC( 158 ) / 'VROCP5ALKP2     ' /
      DATA GAS_CHEM_SPC( 159 ) / 'VROCP4ALKP2     ' /
      DATA GAS_CHEM_SPC( 160 ) / 'VROCP2OXY2      ' /
      DATA GAS_CHEM_SPC( 161 ) / 'VROCP3ALKP2     ' /
      DATA GAS_CHEM_SPC( 162 ) / 'VROCP1OXY1      ' /
      DATA GAS_CHEM_SPC( 163 ) / 'VROCP2ALKP2     ' /
      DATA GAS_CHEM_SPC( 164 ) / 'VROCP1ALKP2     ' /
      DATA GAS_CHEM_SPC( 165 ) / 'VROCN1OXY1      ' /
      DATA GAS_CHEM_SPC( 166 ) / 'HC10P2          ' /
      DATA GAS_CHEM_SPC( 167 ) / 'VROCP6ARO       ' /
      DATA GAS_CHEM_SPC( 168 ) / 'VROCP6AROP      ' /
      DATA GAS_CHEM_SPC( 169 ) / 'VROCN2OXY4      ' /
      DATA GAS_CHEM_SPC( 170 ) / 'VROCN1OXY3      ' /
      DATA GAS_CHEM_SPC( 171 ) / 'VROCP5ARO       ' /
      DATA GAS_CHEM_SPC( 172 ) / 'VROCP5AROP      ' /
      DATA GAS_CHEM_SPC( 173 ) / 'NAPH            ' /
      DATA GAS_CHEM_SPC( 174 ) / 'NAPHP           ' /
      DATA GAS_CHEM_SPC( 175 ) / 'VROCN2OXY8      ' /
      DATA GAS_CHEM_SPC( 176 ) / 'VROCP5OXY1      ' /
      DATA GAS_CHEM_SPC( 177 ) / 'VROCP6OXY1      ' /
      DATA GAS_CHEM_SPC( 178 ) / 'ECH4            ' /
      DATA GAS_CHEM_SPC( 179 ) / 'STY             ' /
      DATA GAS_CHEM_SPC( 180 ) / 'STYP            ' /
      DATA GAS_CHEM_SPC( 181 ) / 'CL2             ' /
      DATA GAS_CHEM_SPC( 182 ) / 'CL              ' /
      DATA GAS_CHEM_SPC( 183 ) / 'CLO             ' /
      DATA GAS_CHEM_SPC( 184 ) / 'OCLO            ' /
      DATA GAS_CHEM_SPC( 185 ) / 'CL2O2           ' /
      DATA GAS_CHEM_SPC( 186 ) / 'CLOO            ' /
      DATA GAS_CHEM_SPC( 187 ) / 'HOCL            ' /
      DATA GAS_CHEM_SPC( 188 ) / 'CLNO            ' /
      DATA GAS_CHEM_SPC( 189 ) / 'CLNO2           ' /
      DATA GAS_CHEM_SPC( 190 ) / 'CLNO3           ' /
      DATA GAS_CHEM_SPC( 191 ) / 'HCOCL           ' /
      DATA GAS_CHEM_SPC( 192 ) / 'HCL             ' /
      DATA GAS_CHEM_SPC( 193 ) / 'BR2             ' /
      DATA GAS_CHEM_SPC( 194 ) / 'BR              ' /
      DATA GAS_CHEM_SPC( 195 ) / 'BRO             ' /
      DATA GAS_CHEM_SPC( 196 ) / 'OBRO            ' /
      DATA GAS_CHEM_SPC( 197 ) / 'HOBR            ' /
      DATA GAS_CHEM_SPC( 198 ) / 'BRNO            ' /
      DATA GAS_CHEM_SPC( 199 ) / 'BRNO2           ' /
      DATA GAS_CHEM_SPC( 200 ) / 'BRNO3           ' /
      DATA GAS_CHEM_SPC( 201 ) / 'CH2BR2          ' /
      DATA GAS_CHEM_SPC( 202 ) / 'CHBR3           ' /
      DATA GAS_CHEM_SPC( 203 ) / 'HCOBR           ' /
      DATA GAS_CHEM_SPC( 204 ) / 'HBR             ' /
      DATA GAS_CHEM_SPC( 205 ) / 'I2              ' /
      DATA GAS_CHEM_SPC( 206 ) / 'I               ' /
      DATA GAS_CHEM_SPC( 207 ) / 'IO              ' /
      DATA GAS_CHEM_SPC( 208 ) / 'OIO             ' /
      DATA GAS_CHEM_SPC( 209 ) / 'I2O2            ' /
      DATA GAS_CHEM_SPC( 210 ) / 'HOI             ' /
      DATA GAS_CHEM_SPC( 211 ) / 'HI              ' /
      DATA GAS_CHEM_SPC( 212 ) / 'INO             ' /
      DATA GAS_CHEM_SPC( 213 ) / 'INO2            ' /
      DATA GAS_CHEM_SPC( 214 ) / 'INO3            ' /
      DATA GAS_CHEM_SPC( 215 ) / 'CH3I            ' /
      DATA GAS_CHEM_SPC( 216 ) / 'CH2I2           ' /
      DATA GAS_CHEM_SPC( 217 ) / 'I2O3            ' /
      DATA GAS_CHEM_SPC( 218 ) / 'I2O4            ' /
      DATA GAS_CHEM_SPC( 219 ) / 'BRCL            ' /
      DATA GAS_CHEM_SPC( 220 ) / 'ICL             ' /
      DATA GAS_CHEM_SPC( 221 ) / 'IBR             ' /
      DATA GAS_CHEM_SPC( 222 ) / 'CH2IBR          ' /
      DATA GAS_CHEM_SPC( 223 ) / 'CH2ICL          ' /
      DATA GAS_CHEM_SPC( 224 ) / 'CHBR2CL         ' /
      DATA GAS_CHEM_SPC( 225 ) / 'CHBRCL2         ' /
      DATA GAS_CHEM_SPC( 226 ) / 'CH2BRCL         ' /
      DATA GAS_CHEM_SPC( 227 ) / 'NO2PIJ          ' /
      DATA GAS_CHEM_SPC( 228 ) / 'NO2PK           ' /




      LOGICAL   :: HALOGEN_PARMAETER = .FALSE. 


! The below character and integer arrays list the model species names used in the 
! chemical mechanism. The gas species and their order should agree with 
! the GC_SPC array for the gas phase chemistry to work correctly. 
! If present, the CHEMISTRY_SPC names and species type should agree with the CGRID_SPCS module


! MAPPED_TO_CGRID declares whether CMAQ namelists were used to determine 
! the below values of CGRID_INDEX, SPECIES_TYPE, SPECIES_MOLWT, and CONVERT_CONC
      LOGICAL, PARAMETER, PRIVATE :: F = .FALSE.
      LOGICAL, PARAMETER, PRIVATE :: T = .TRUE.


      LOGICAL   :: MAPPED_TO_CGRID   = .FALSE. 

      TYPE MEMBER
         CHARACTER( 16 ) :: CHEMISTRY_SPC
         INTEGER         :: CGRID_INDEX
         CHARACTER(  2 ) :: SPECIES_TYPE
         REAL( 8 )       :: SPECIES_MOLWT
         LOGICAL         :: CONVERT_CONC
      END TYPE MEMBER
      TYPE( MEMBER ) ::  SPECIES_LIST( NUMB_MECH_SPC ) = (/ &
      & MEMBER("O3              ",    1, "GC",   48.00D0, F), &
      & MEMBER("O3P             ",    2, "GC",   16.00D0, F), &
      & MEMBER("O1D             ",    3, "GC",   16.00D0, F), &
      & MEMBER("H2O2            ",    4, "GC",   34.00D0, F), &
      & MEMBER("HO              ",    5, "GC",   17.00D0, F), &
      & MEMBER("NO2             ",    6, "GC",   46.00D0, F), &
      & MEMBER("NO              ",    7, "GC",   30.00D0, F), &
      & MEMBER("NO3             ",    8, "GC",   62.00D0, F), &
      & MEMBER("HONO            ",    9, "GC",   47.00D0, F), &
      & MEMBER("HNO3            ",   10, "GC",   63.00D0, F), &
      & MEMBER("HNO4            ",   11, "GC",   79.00D0, F), &
      & MEMBER("HO2             ",   12, "GC",   33.00D0, F), &
      & MEMBER("HCHO            ",   13, "GC",   30.00D0, F), &
      & MEMBER("CO              ",   14, "GC",   28.00D0, F), &
      & MEMBER("ACD             ",   15, "GC",   44.00D0, F), &
      & MEMBER("MO2             ",   16, "GC",   47.00D0, F), &
      & MEMBER("ALD             ",   17, "GC",   58.00D0, F), &
      & MEMBER("ETHP            ",   18, "GC",   61.00D0, F), &
      & MEMBER("ACT             ",   19, "GC",   58.00D0, F), &
      & MEMBER("ACO3            ",   20, "GC",   75.00D0, F), &
      & MEMBER("UALD            ",   21, "GC",   84.10D0, F), &
      & MEMBER("KET             ",   22, "GC",   86.00D0, F), &
      & MEMBER("PINAL           ",   81, "GC",  168.00D0, F), &
      & MEMBER("PINALP          ",   82, "GC",  199.00D0, F), &
      & MEMBER("LIMAL           ",   88, "GC",  168.00D0, F), &
      & MEMBER("LIMALP          ",   89, "GC",  217.00D0, F), &
      & MEMBER("MEK             ",   23, "GC",   72.10D0, F), &
      & MEMBER("HKET            ",   24, "GC",   74.00D0, F), &
      & MEMBER("MACR            ",   25, "GC",   70.00D0, F), &
      & MEMBER("MACP            ",   26, "GC",  101.00D0, F), &
      & MEMBER("XO2             ",   27, "GC",    1.00D0, F), &
      & MEMBER("MVK             ",   28, "GC",   70.10D0, F), &
      & MEMBER("GLY             ",   29, "GC",   58.00D0, F), &
      & MEMBER("MGLY            ",   30, "GC",   72.00D0, F), &
      & MEMBER("DCB1            ",   31, "GC",   98.00D0, F), &
      & MEMBER("DCB2            ",   32, "GC",  112.10D0, F), &
      & MEMBER("BALD            ",   33, "GC",  106.00D0, F), &
      & MEMBER("BEN             ",   62, "GC",   78.11D0, F), &
      & MEMBER("BAL1            ",  126, "GC",  109.10D0, F), &
      & MEMBER("OP1             ",   35, "GC",   48.00D0, F), &
      & MEMBER("OP2             ",   36, "GC",   62.00D0, F), &
      & MEMBER("OPB             ",   37, "GC",  186.20D0, F), &
      & MEMBER("VOP3            ",   38, "GC",  176.20D0, F), &
      & MEMBER("PAA             ",   39, "GC",   76.00D0, F), &
      & MEMBER("ONIT            ",   40, "GC",  119.00D0, F), &
      & MEMBER("PAN             ",   41, "GC",  121.00D0, F), &
      & MEMBER("CO2             ",  183, "GC",   44.00D0, F), &
      & MEMBER("PPN             ",  114, "GC",  135.00D0, F), &
      & MEMBER("RCO3            ",   92, "GC",   89.00D0, F), &
      & MEMBER("HC3P            ",   48, "GC",   75.00D0, F), &
      & MEMBER("VTRPN           ",  120, "GC",  215.00D0, F), &
      & MEMBER("VHONIT          ",  121, "GC",  265.00D0, F), &
      & MEMBER("N2O5            ",   42, "GC",  108.00D0, F), &
      & MEMBER("SO2             ",   43, "GC",   64.00D0, F), &
      & MEMBER("SULF            ",   44, "GC",   98.00D0, F), &
      & MEMBER("SULRXN          ",   45, "GC",   98.00D0, F), &
      & MEMBER("ETH             ",   46, "GC",   30.10D0, F), &
      & MEMBER("HC3             ",   47, "GC",   44.10D0, F), &
      & MEMBER("ASOATJ          ",  243, "AE",  200.00D0, T), &
      & MEMBER("HC5             ",   49, "GC",   72.10D0, F), &
      & MEMBER("HC5P            ",   50, "GC",  103.00D0, F), &
      & MEMBER("ETE             ",   54, "GC",   28.10D0, F), &
      & MEMBER("ETEP            ",   55, "GC",   77.00D0, F), &
      & MEMBER("OLT             ",   56, "GC",   42.00D0, F), &
      & MEMBER("OLTP            ",   57, "GC",   91.00D0, F), &
      & MEMBER("OLI             ",   58, "GC",   70.10D0, F), &
      & MEMBER("OLIP            ",   59, "GC",  119.00D0, F), &
      & MEMBER("ACE             ",   60, "GC",   26.00D0, F), &
      & MEMBER("ORA1            ",   61, "GC",   46.00D0, F), &
      & MEMBER("BENP            ",   63, "GC",  159.12D0, F), &
      & MEMBER("PHEN            ",   64, "GC",   94.11D0, F), &
      & MEMBER("TOL             ",   65, "GC",   92.14D0, F), &
      & MEMBER("TOLP            ",   70, "GC",  173.14D0, F), &
      & MEMBER("CSL             ",   66, "GC",  108.14D0, F), &
      & MEMBER("XYL             ",   67, "GC",  106.20D0, F), &
      & MEMBER("XYLP            ",   71, "GC",  187.17D0, F), &
      & MEMBER("EBZ             ",   68, "GC",  106.20D0, F), &
      & MEMBER("EBZP            ",   72, "GC",  187.17D0, F), &
      & MEMBER("ISO             ",   74, "GC",   68.10D0, F), &
      & MEMBER("ISON            ",  118, "GC",  147.00D0, F), &
      & MEMBER("ISONP           ",  117, "GC",  162.10D0, F), &
      & MEMBER("ISOP            ",   75, "GC",  117.00D0, F), &
      & MEMBER("ISHP            ",  108, "GC",  118.00D0, F), &
      & MEMBER("IEPOX           ",  109, "GC",  118.10D0, F), &
      & MEMBER("IPX             ",  119, "GC",  150.10D0, F), &
      & MEMBER("INALD           ",  116, "GC",  195.10D0, F), &
      & MEMBER("ROH             ",  106, "GC",   60.00D0, F), &
      & MEMBER("API             ",   76, "GC",  136.40D0, F), &
      & MEMBER("APIP1           ",   77, "GC",  185.00D0, F), &
      & MEMBER("APIP2           ",   78, "GC",  249.00D0, F), &
      & MEMBER("LIM             ",   83, "GC",  136.30D0, F), &
      & MEMBER("LIMP1           ",   84, "GC",  185.00D0, F), &
      & MEMBER("LIMP2           ",   85, "GC",  249.00D0, F), &
      & MEMBER("ACTP            ",   93, "GC",   89.00D0, F), &
      & MEMBER("MEKP            ",   94, "GC",  103.00D0, F), &
      & MEMBER("KETP            ",   95, "GC",  117.00D0, F), &
      & MEMBER("MCP             ",   96, "GC",  119.00D0, F), &
      & MEMBER("MVKP            ",   97, "GC",  119.00D0, F), &
      & MEMBER("UALP            ",   98, "GC",  133.00D0, F), &
      & MEMBER("DCB3            ",   99, "GC",   84.00D0, F), &
      & MEMBER("BALP            ",  100, "GC",  137.00D0, F), &
      & MEMBER("ADDC            ",  101, "GC",  125.00D0, F), &
      & MEMBER("CHO             ",   34, "GC",  107.00D0, F), &
      & MEMBER("MCT             ",  102, "GC",  110.11D0, F), &
      & MEMBER("MCTO            ",  103, "GC",  109.00D0, F), &
      & MEMBER("MOH             ",  104, "GC",   32.00D0, F), &
      & MEMBER("EOH             ",  105, "GC",   46.10D0, F), &
      & MEMBER("ETEG            ",  107, "GC",   62.10D0, F), &
      & MEMBER("HC10P           ",   52, "GC",  173.27D0, F), &
      & MEMBER("MAHP            ",  111, "GC",  102.00D0, F), &
      & MEMBER("ORA2            ",  112, "GC",   60.20D0, F), &
      & MEMBER("ORAP            ",  113, "GC",   91.00D0, F), &
      & MEMBER("MPAN            ",  115, "GC",  147.10D0, F), &
      & MEMBER("MCTP            ",  122, "GC",  125.00D0, F), &
      & MEMBER("OLNN            ",  123, "GC",  136.00D0, F), &
      & MEMBER("OLND            ",  124, "GC",  136.00D0, F), &
      & MEMBER("APINP1          ",   79, "GC",  230.00D0, F), &
      & MEMBER("APINP2          ",   80, "GC",  294.00D0, F), &
      & MEMBER("LIMNP1          ",   86, "GC",  230.00D0, F), &
      & MEMBER("LIMNP2          ",   87, "GC",  294.00D0, F), &
      & MEMBER("ADCN            ",  125, "GC",  155.00D0, F), &
      & MEMBER("VHOM            ",   90, "GC",  250.00D0, F), &
      & MEMBER("VROCP4OXY2      ",  179, "GC",  158.20D0, F), &
      & MEMBER("VROCN1OXY6      ",  172, "GC",  190.20D0, F), &
      & MEMBER("FURANONE        ",  134, "GC",  100.10D0, F), &
      & MEMBER("VROCP3OXY2      ",  178, "GC",  186.30D0, F), &
      & MEMBER("VROCP0OXY4      ",  174, "GC",  202.30D0, F), &
      & MEMBER("BAL2            ",  127, "GC",   93.10D0, F), &
      & MEMBER("VELHOM          ",   91, "GC",  402.00D0, F), &
      & MEMBER("VROCIOXY        ",  135, "GC",  247.00D0, F), &
      & MEMBER("SLOWROC         ",  136, "GC",   75.40D0, F), &
      & MEMBER("ACRO            ",  128, "GC",   56.10D0, F), &
      & MEMBER("BDE13           ",  129, "GC",   54.10D0, F), &
      & MEMBER("BDE13P          ",  130, "GC",  103.00D0, F), &
      & MEMBER("FURAN           ",  132, "GC",   96.10D0, F), &
      & MEMBER("FURANO2         ",  133, "GC",  145.10D0, F), &
      & MEMBER("PROG            ",  131, "GC",   76.10D0, F), &
      & MEMBER("SESQ            ",  137, "GC",  204.40D0, F), &
      & MEMBER("SESQNRO2        ",  139, "GC",  298.40D0, F), &
      & MEMBER("VROCN2OXY2      ",  167, "GC",  282.40D0, F), &
      & MEMBER("SESQRO2         ",  138, "GC",  253.40D0, F), &
      & MEMBER("VROCP0OXY2      ",  173, "GC",  242.40D0, F), &
      & MEMBER("VROCP1OXY3      ",  176, "GC",  202.30D0, F), &
      & MEMBER("AGLYOLIGJ       ",  244, "AE",   66.40D0, T), &
      & MEMBER("IEPOXP          ",  110, "GC",  118.10D0, F), &
      & MEMBER("AISO3NOSJ       ",  318, "AE",  136.20D0, T), &
      & MEMBER("ASO4J           ",  235, "AE",   96.00D0, T), &
      & MEMBER("AISO3OSJ        ",  319, "AE",  216.20D0, T), &
      & MEMBER("AISO4J          ",  320, "AE",  168.10D0, T), &
      & MEMBER("AISO5J          ",  321, "AE",  150.10D0, T), &
      & MEMBER("VROCP6ALK       ",  154, "GC",  184.37D0, F), &
      & MEMBER("VROCP6ALKP      ",  160, "GC",  215.36D0, F), &
      & MEMBER("VROCP5ALK       ",  153, "GC",  198.39D0, F), &
      & MEMBER("VROCP5ALKP      ",  159, "GC",  229.38D0, F), &
      & MEMBER("VROCP4ALK       ",  152, "GC",  254.50D0, F), &
      & MEMBER("VROCP4ALKP      ",  158, "GC",  285.49D0, F), &
      & MEMBER("VROCP3ALK       ",  151, "GC",  296.58D0, F), &
      & MEMBER("VROCP3ALKP      ",  157, "GC",  327.57D0, F), &
      & MEMBER("VROCP2ALK       ",  150, "GC",  338.66D0, F), &
      & MEMBER("VROCP2ALKP      ",  156, "GC",  369.65D0, F), &
      & MEMBER("VROCP1ALK       ",  149, "GC",  380.75D0, F), &
      & MEMBER("VROCP1ALKP      ",  155, "GC",  411.74D0, F), &
      & MEMBER("HC10            ",   51, "GC",  142.28D0, F), &
      & MEMBER("VROCP6ALKP2     ",  166, "GC",  231.36D0, F), &
      & MEMBER("VROCP5ALKP2     ",  165, "GC",  245.38D0, F), &
      & MEMBER("VROCP4ALKP2     ",  164, "GC",  301.49D0, F), &
      & MEMBER("VROCP2OXY2      ",  177, "GC",  200.30D0, F), &
      & MEMBER("VROCP3ALKP2     ",  163, "GC",  343.57D0, F), &
      & MEMBER("VROCP1OXY1      ",  175, "GC",  270.50D0, F), &
      & MEMBER("VROCP2ALKP2     ",  162, "GC",  385.65D0, F), &
      & MEMBER("VROCP1ALKP2     ",  161, "GC",  427.73D0, F), &
      & MEMBER("VROCN1OXY1      ",  170, "GC",  312.50D0, F), &
      & MEMBER("HC10P2          ",   53, "GC",  189.27D0, F), &
      & MEMBER("VROCP6ARO       ",  143, "GC",  176.30D0, F), &
      & MEMBER("VROCP6AROP      ",  145, "GC",  257.30D0, F), &
      & MEMBER("VROCN2OXY4      ",  168, "GC",  232.30D0, F), &
      & MEMBER("VROCN1OXY3      ",  171, "GC",  230.30D0, F), &
      & MEMBER("VROCP5ARO       ",  142, "GC",  190.33D0, F), &
      & MEMBER("VROCP5AROP      ",  144, "GC",  271.33D0, F), &
      & MEMBER("NAPH            ",  140, "GC",  128.17D0, F), &
      & MEMBER("NAPHP           ",  141, "GC",  209.17D0, F), &
      & MEMBER("VROCN2OXY8      ",  169, "GC",  194.20D0, F), &
      & MEMBER("VROCP5OXY1      ",  180, "GC",  170.30D0, F), &
      & MEMBER("VROCP6OXY1      ",  181, "GC",  142.20D0, F), &
      & MEMBER("ECH4            ",  182, "GC",   16.00D0, F), &
      & MEMBER("ATRPNJ          ",  322, "AE",  215.00D0, T), &
      & MEMBER("AHOMJ           ",  245, "AE",  250.00D0, T), &
      & MEMBER("AHONITJ         ",  323, "AE",  265.00D0, T), &
      & MEMBER("STY             ",   69, "GC",  104.00D0, F), &
      & MEMBER("STYP            ",   73, "GC",  153.00D0, F), &
      & MEMBER("ANO3I           ",  240, "AE",   62.00D0, T), &
      & MEMBER("ANO3J           ",  241, "AE",   62.00D0, T), &
      & MEMBER("AMTN1J          ",  325, "AE",  262.00D0, T), &
      & MEMBER("CL2             ",  185, "GC",   71.00D0, F), &
      & MEMBER("CL              ",  186, "GC",   35.50D0, F), &
      & MEMBER("CLO             ",  187, "GC",   51.50D0, F), &
      & MEMBER("OCLO            ",  188, "GC",   67.50D0, F), &
      & MEMBER("CL2O2           ",  189, "GC",  103.00D0, F), &
      & MEMBER("CLOO            ",  196, "GC",   67.50D0, F), &
      & MEMBER("HOCL            ",  190, "GC",   52.50D0, F), &
      & MEMBER("CLNO            ",  191, "GC",   65.50D0, F), &
      & MEMBER("CLNO2           ",  192, "GC",   81.50D0, F), &
      & MEMBER("CLNO3           ",  193, "GC",   97.50D0, F), &
      & MEMBER("HCOCL           ",  194, "GC",   64.50D0, F), &
      & MEMBER("HCL             ",  195, "GC",   36.00D0, F), &
      & MEMBER("BR2             ",  197, "GC",  159.80D0, F), &
      & MEMBER("BR              ",  198, "GC",   79.90D0, F), &
      & MEMBER("BRO             ",  199, "GC",   95.90D0, F), &
      & MEMBER("OBRO            ",  200, "GC",  111.90D0, F), &
      & MEMBER("HOBR            ",  201, "GC",   96.90D0, F), &
      & MEMBER("BRNO            ",  202, "GC",  109.90D0, F), &
      & MEMBER("BRNO2           ",  203, "GC",  125.90D0, F), &
      & MEMBER("BRNO3           ",  204, "GC",  141.90D0, F), &
      & MEMBER("CH2BR2          ",  205, "GC",  173.80D0, F), &
      & MEMBER("CHBR3           ",  206, "GC",  252.70D0, F), &
      & MEMBER("HCOBR           ",  208, "GC",  108.90D0, F), &
      & MEMBER("HBR             ",  207, "GC",   80.90D0, F), &
      & MEMBER("I2              ",  209, "GC",  253.80D0, F), &
      & MEMBER("I               ",  210, "GC",  126.90D0, F), &
      & MEMBER("IO              ",  211, "GC",  142.90D0, F), &
      & MEMBER("OIO             ",  212, "GC",  158.90D0, F), &
      & MEMBER("I2O2            ",  213, "GC",  285.80D0, F), &
      & MEMBER("HOI             ",  214, "GC",  143.90D0, F), &
      & MEMBER("HI              ",  215, "GC",  127.90D0, F), &
      & MEMBER("INO             ",  216, "GC",  156.90D0, F), &
      & MEMBER("INO2            ",  217, "GC",  172.90D0, F), &
      & MEMBER("INO3            ",  218, "GC",  188.90D0, F), &
      & MEMBER("CH3I            ",  219, "GC",  141.90D0, F), &
      & MEMBER("CH2I2           ",  220, "GC",  267.80D0, F), &
      & MEMBER("I2O3            ",  231, "GC",  301.80D0, F), &
      & MEMBER("I2O4            ",  232, "GC",  317.80D0, F), &
      & MEMBER("BRCL            ",  221, "GC",  115.40D0, F), &
      & MEMBER("ICL             ",  222, "GC",  162.40D0, F), &
      & MEMBER("IBR             ",  223, "GC",  206.80D0, F), &
      & MEMBER("CH2IBR          ",  224, "GC",  220.80D0, F), &
      & MEMBER("CH2ICL          ",  225, "GC",  176.40D0, F), &
      & MEMBER("CHBR2CL         ",  226, "GC",  208.30D0, F), &
      & MEMBER("CHBRCL2         ",  227, "GC",  163.90D0, F), &
      & MEMBER("CH2BRCL         ",  228, "GC",  129.40D0, F), &
      & MEMBER("NO2PIJ          ",  229, "GC",   46.00D0, F), &
      & MEMBER("NO2PK           ",  230, "GC",   46.00D0, F), &
      & MEMBER("ACLI            ",  277, "AE",   35.50D0, T), &
      & MEMBER("ACLJ            ",  278, "AE",   35.50D0, T), &
      & MEMBER("ACLK            ",  279, "AE",   35.50D0, T), &
      & MEMBER("ABRJ            ",  326, "AE",   79.90D0, T), &
      & MEMBER("AFEJ            ",  252, "AE",   55.80D0, T), &
      & MEMBER("AMNJ            ",  259, "AE",   54.90D0, T), &
      & MEMBER("AHMSJ           ",  328, "AE",  111.10D0, T) /)

      DATA CHEMISTRY_SPC(   1 ), SPECIES_MOLWT(   1 ) / 'O3              ',   48.00D0 /
      DATA CHEMISTRY_SPC(   2 ), SPECIES_MOLWT(   2 ) / 'O3P             ',   16.00D0 /
      DATA CHEMISTRY_SPC(   3 ), SPECIES_MOLWT(   3 ) / 'O1D             ',   16.00D0 /
      DATA CHEMISTRY_SPC(   4 ), SPECIES_MOLWT(   4 ) / 'H2O2            ',   34.00D0 /
      DATA CHEMISTRY_SPC(   5 ), SPECIES_MOLWT(   5 ) / 'HO              ',   17.00D0 /
      DATA CHEMISTRY_SPC(   6 ), SPECIES_MOLWT(   6 ) / 'NO2             ',   46.00D0 /
      DATA CHEMISTRY_SPC(   7 ), SPECIES_MOLWT(   7 ) / 'NO              ',   30.00D0 /
      DATA CHEMISTRY_SPC(   8 ), SPECIES_MOLWT(   8 ) / 'NO3             ',   62.00D0 /
      DATA CHEMISTRY_SPC(   9 ), SPECIES_MOLWT(   9 ) / 'HONO            ',   47.00D0 /
      DATA CHEMISTRY_SPC(  10 ), SPECIES_MOLWT(  10 ) / 'HNO3            ',   63.00D0 /
      DATA CHEMISTRY_SPC(  11 ), SPECIES_MOLWT(  11 ) / 'HNO4            ',   79.00D0 /
      DATA CHEMISTRY_SPC(  12 ), SPECIES_MOLWT(  12 ) / 'HO2             ',   33.00D0 /
      DATA CHEMISTRY_SPC(  13 ), SPECIES_MOLWT(  13 ) / 'HCHO            ',   30.00D0 /
      DATA CHEMISTRY_SPC(  14 ), SPECIES_MOLWT(  14 ) / 'CO              ',   28.00D0 /
      DATA CHEMISTRY_SPC(  15 ), SPECIES_MOLWT(  15 ) / 'ACD             ',   44.00D0 /
      DATA CHEMISTRY_SPC(  16 ), SPECIES_MOLWT(  16 ) / 'MO2             ',   47.00D0 /
      DATA CHEMISTRY_SPC(  17 ), SPECIES_MOLWT(  17 ) / 'ALD             ',   58.00D0 /
      DATA CHEMISTRY_SPC(  18 ), SPECIES_MOLWT(  18 ) / 'ETHP            ',   61.00D0 /
      DATA CHEMISTRY_SPC(  19 ), SPECIES_MOLWT(  19 ) / 'ACT             ',   58.00D0 /
      DATA CHEMISTRY_SPC(  20 ), SPECIES_MOLWT(  20 ) / 'ACO3            ',   75.00D0 /
      DATA CHEMISTRY_SPC(  21 ), SPECIES_MOLWT(  21 ) / 'UALD            ',   84.10D0 /
      DATA CHEMISTRY_SPC(  22 ), SPECIES_MOLWT(  22 ) / 'KET             ',   86.00D0 /
      DATA CHEMISTRY_SPC(  23 ), SPECIES_MOLWT(  23 ) / 'PINAL           ',  168.00D0 /
      DATA CHEMISTRY_SPC(  24 ), SPECIES_MOLWT(  24 ) / 'PINALP          ',  199.00D0 /
      DATA CHEMISTRY_SPC(  25 ), SPECIES_MOLWT(  25 ) / 'LIMAL           ',  168.00D0 /
      DATA CHEMISTRY_SPC(  26 ), SPECIES_MOLWT(  26 ) / 'LIMALP          ',  217.00D0 /
      DATA CHEMISTRY_SPC(  27 ), SPECIES_MOLWT(  27 ) / 'MEK             ',   72.10D0 /
      DATA CHEMISTRY_SPC(  28 ), SPECIES_MOLWT(  28 ) / 'HKET            ',   74.00D0 /
      DATA CHEMISTRY_SPC(  29 ), SPECIES_MOLWT(  29 ) / 'MACR            ',   70.00D0 /
      DATA CHEMISTRY_SPC(  30 ), SPECIES_MOLWT(  30 ) / 'MACP            ',  101.00D0 /
      DATA CHEMISTRY_SPC(  31 ), SPECIES_MOLWT(  31 ) / 'XO2             ',    1.00D0 /
      DATA CHEMISTRY_SPC(  32 ), SPECIES_MOLWT(  32 ) / 'MVK             ',   70.10D0 /
      DATA CHEMISTRY_SPC(  33 ), SPECIES_MOLWT(  33 ) / 'GLY             ',   58.00D0 /
      DATA CHEMISTRY_SPC(  34 ), SPECIES_MOLWT(  34 ) / 'MGLY            ',   72.00D0 /
      DATA CHEMISTRY_SPC(  35 ), SPECIES_MOLWT(  35 ) / 'DCB1            ',   98.00D0 /
      DATA CHEMISTRY_SPC(  36 ), SPECIES_MOLWT(  36 ) / 'DCB2            ',  112.10D0 /
      DATA CHEMISTRY_SPC(  37 ), SPECIES_MOLWT(  37 ) / 'BALD            ',  106.00D0 /
      DATA CHEMISTRY_SPC(  38 ), SPECIES_MOLWT(  38 ) / 'BEN             ',   78.11D0 /
      DATA CHEMISTRY_SPC(  39 ), SPECIES_MOLWT(  39 ) / 'BAL1            ',  109.10D0 /
      DATA CHEMISTRY_SPC(  40 ), SPECIES_MOLWT(  40 ) / 'OP1             ',   48.00D0 /
      DATA CHEMISTRY_SPC(  41 ), SPECIES_MOLWT(  41 ) / 'OP2             ',   62.00D0 /
      DATA CHEMISTRY_SPC(  42 ), SPECIES_MOLWT(  42 ) / 'OPB             ',  186.20D0 /
      DATA CHEMISTRY_SPC(  43 ), SPECIES_MOLWT(  43 ) / 'VOP3            ',  176.20D0 /
      DATA CHEMISTRY_SPC(  44 ), SPECIES_MOLWT(  44 ) / 'PAA             ',   76.00D0 /
      DATA CHEMISTRY_SPC(  45 ), SPECIES_MOLWT(  45 ) / 'ONIT            ',  119.00D0 /
      DATA CHEMISTRY_SPC(  46 ), SPECIES_MOLWT(  46 ) / 'PAN             ',  121.00D0 /
      DATA CHEMISTRY_SPC(  47 ), SPECIES_MOLWT(  47 ) / 'CO2             ',   44.00D0 /
      DATA CHEMISTRY_SPC(  48 ), SPECIES_MOLWT(  48 ) / 'PPN             ',  135.00D0 /
      DATA CHEMISTRY_SPC(  49 ), SPECIES_MOLWT(  49 ) / 'RCO3            ',   89.00D0 /
      DATA CHEMISTRY_SPC(  50 ), SPECIES_MOLWT(  50 ) / 'HC3P            ',   75.00D0 /
      DATA CHEMISTRY_SPC(  51 ), SPECIES_MOLWT(  51 ) / 'VTRPN           ',  215.00D0 /
      DATA CHEMISTRY_SPC(  52 ), SPECIES_MOLWT(  52 ) / 'VHONIT          ',  265.00D0 /
      DATA CHEMISTRY_SPC(  53 ), SPECIES_MOLWT(  53 ) / 'N2O5            ',  108.00D0 /
      DATA CHEMISTRY_SPC(  54 ), SPECIES_MOLWT(  54 ) / 'SO2             ',   64.00D0 /
      DATA CHEMISTRY_SPC(  55 ), SPECIES_MOLWT(  55 ) / 'SULF            ',   98.00D0 /
      DATA CHEMISTRY_SPC(  56 ), SPECIES_MOLWT(  56 ) / 'SULRXN          ',   98.00D0 /
      DATA CHEMISTRY_SPC(  57 ), SPECIES_MOLWT(  57 ) / 'ETH             ',   30.10D0 /
      DATA CHEMISTRY_SPC(  58 ), SPECIES_MOLWT(  58 ) / 'HC3             ',   44.10D0 /
      DATA CHEMISTRY_SPC(  59 ), SPECIES_MOLWT(  59 ) / 'ASOATJ          ',  200.00D0 /
      DATA CHEMISTRY_SPC(  60 ), SPECIES_MOLWT(  60 ) / 'HC5             ',   72.10D0 /
      DATA CHEMISTRY_SPC(  61 ), SPECIES_MOLWT(  61 ) / 'HC5P            ',  103.00D0 /
      DATA CHEMISTRY_SPC(  62 ), SPECIES_MOLWT(  62 ) / 'ETE             ',   28.10D0 /
      DATA CHEMISTRY_SPC(  63 ), SPECIES_MOLWT(  63 ) / 'ETEP            ',   77.00D0 /
      DATA CHEMISTRY_SPC(  64 ), SPECIES_MOLWT(  64 ) / 'OLT             ',   42.00D0 /
      DATA CHEMISTRY_SPC(  65 ), SPECIES_MOLWT(  65 ) / 'OLTP            ',   91.00D0 /
      DATA CHEMISTRY_SPC(  66 ), SPECIES_MOLWT(  66 ) / 'OLI             ',   70.10D0 /
      DATA CHEMISTRY_SPC(  67 ), SPECIES_MOLWT(  67 ) / 'OLIP            ',  119.00D0 /
      DATA CHEMISTRY_SPC(  68 ), SPECIES_MOLWT(  68 ) / 'ACE             ',   26.00D0 /
      DATA CHEMISTRY_SPC(  69 ), SPECIES_MOLWT(  69 ) / 'ORA1            ',   46.00D0 /
      DATA CHEMISTRY_SPC(  70 ), SPECIES_MOLWT(  70 ) / 'BENP            ',  159.12D0 /
      DATA CHEMISTRY_SPC(  71 ), SPECIES_MOLWT(  71 ) / 'PHEN            ',   94.11D0 /
      DATA CHEMISTRY_SPC(  72 ), SPECIES_MOLWT(  72 ) / 'TOL             ',   92.14D0 /
      DATA CHEMISTRY_SPC(  73 ), SPECIES_MOLWT(  73 ) / 'TOLP            ',  173.14D0 /
      DATA CHEMISTRY_SPC(  74 ), SPECIES_MOLWT(  74 ) / 'CSL             ',  108.14D0 /
      DATA CHEMISTRY_SPC(  75 ), SPECIES_MOLWT(  75 ) / 'XYL             ',  106.20D0 /
      DATA CHEMISTRY_SPC(  76 ), SPECIES_MOLWT(  76 ) / 'XYLP            ',  187.17D0 /
      DATA CHEMISTRY_SPC(  77 ), SPECIES_MOLWT(  77 ) / 'EBZ             ',  106.20D0 /
      DATA CHEMISTRY_SPC(  78 ), SPECIES_MOLWT(  78 ) / 'EBZP            ',  187.17D0 /
      DATA CHEMISTRY_SPC(  79 ), SPECIES_MOLWT(  79 ) / 'ISO             ',   68.10D0 /
      DATA CHEMISTRY_SPC(  80 ), SPECIES_MOLWT(  80 ) / 'ISON            ',  147.00D0 /
      DATA CHEMISTRY_SPC(  81 ), SPECIES_MOLWT(  81 ) / 'ISONP           ',  162.10D0 /
      DATA CHEMISTRY_SPC(  82 ), SPECIES_MOLWT(  82 ) / 'ISOP            ',  117.00D0 /
      DATA CHEMISTRY_SPC(  83 ), SPECIES_MOLWT(  83 ) / 'ISHP            ',  118.00D0 /
      DATA CHEMISTRY_SPC(  84 ), SPECIES_MOLWT(  84 ) / 'IEPOX           ',  118.10D0 /
      DATA CHEMISTRY_SPC(  85 ), SPECIES_MOLWT(  85 ) / 'IPX             ',  150.10D0 /
      DATA CHEMISTRY_SPC(  86 ), SPECIES_MOLWT(  86 ) / 'INALD           ',  195.10D0 /
      DATA CHEMISTRY_SPC(  87 ), SPECIES_MOLWT(  87 ) / 'ROH             ',   60.00D0 /
      DATA CHEMISTRY_SPC(  88 ), SPECIES_MOLWT(  88 ) / 'API             ',  136.40D0 /
      DATA CHEMISTRY_SPC(  89 ), SPECIES_MOLWT(  89 ) / 'APIP1           ',  185.00D0 /
      DATA CHEMISTRY_SPC(  90 ), SPECIES_MOLWT(  90 ) / 'APIP2           ',  249.00D0 /
      DATA CHEMISTRY_SPC(  91 ), SPECIES_MOLWT(  91 ) / 'LIM             ',  136.30D0 /
      DATA CHEMISTRY_SPC(  92 ), SPECIES_MOLWT(  92 ) / 'LIMP1           ',  185.00D0 /
      DATA CHEMISTRY_SPC(  93 ), SPECIES_MOLWT(  93 ) / 'LIMP2           ',  249.00D0 /
      DATA CHEMISTRY_SPC(  94 ), SPECIES_MOLWT(  94 ) / 'ACTP            ',   89.00D0 /
      DATA CHEMISTRY_SPC(  95 ), SPECIES_MOLWT(  95 ) / 'MEKP            ',  103.00D0 /
      DATA CHEMISTRY_SPC(  96 ), SPECIES_MOLWT(  96 ) / 'KETP            ',  117.00D0 /
      DATA CHEMISTRY_SPC(  97 ), SPECIES_MOLWT(  97 ) / 'MCP             ',  119.00D0 /
      DATA CHEMISTRY_SPC(  98 ), SPECIES_MOLWT(  98 ) / 'MVKP            ',  119.00D0 /
      DATA CHEMISTRY_SPC(  99 ), SPECIES_MOLWT(  99 ) / 'UALP            ',  133.00D0 /
      DATA CHEMISTRY_SPC( 100 ), SPECIES_MOLWT( 100 ) / 'DCB3            ',   84.00D0 /
      DATA CHEMISTRY_SPC( 101 ), SPECIES_MOLWT( 101 ) / 'BALP            ',  137.00D0 /
      DATA CHEMISTRY_SPC( 102 ), SPECIES_MOLWT( 102 ) / 'ADDC            ',  125.00D0 /
      DATA CHEMISTRY_SPC( 103 ), SPECIES_MOLWT( 103 ) / 'CHO             ',  107.00D0 /
      DATA CHEMISTRY_SPC( 104 ), SPECIES_MOLWT( 104 ) / 'MCT             ',  110.11D0 /
      DATA CHEMISTRY_SPC( 105 ), SPECIES_MOLWT( 105 ) / 'MCTO            ',  109.00D0 /
      DATA CHEMISTRY_SPC( 106 ), SPECIES_MOLWT( 106 ) / 'MOH             ',   32.00D0 /
      DATA CHEMISTRY_SPC( 107 ), SPECIES_MOLWT( 107 ) / 'EOH             ',   46.10D0 /
      DATA CHEMISTRY_SPC( 108 ), SPECIES_MOLWT( 108 ) / 'ETEG            ',   62.10D0 /
      DATA CHEMISTRY_SPC( 109 ), SPECIES_MOLWT( 109 ) / 'HC10P           ',  173.27D0 /
      DATA CHEMISTRY_SPC( 110 ), SPECIES_MOLWT( 110 ) / 'MAHP            ',  102.00D0 /
      DATA CHEMISTRY_SPC( 111 ), SPECIES_MOLWT( 111 ) / 'ORA2            ',   60.20D0 /
      DATA CHEMISTRY_SPC( 112 ), SPECIES_MOLWT( 112 ) / 'ORAP            ',   91.00D0 /
      DATA CHEMISTRY_SPC( 113 ), SPECIES_MOLWT( 113 ) / 'MPAN            ',  147.10D0 /
      DATA CHEMISTRY_SPC( 114 ), SPECIES_MOLWT( 114 ) / 'MCTP            ',  125.00D0 /
      DATA CHEMISTRY_SPC( 115 ), SPECIES_MOLWT( 115 ) / 'OLNN            ',  136.00D0 /
      DATA CHEMISTRY_SPC( 116 ), SPECIES_MOLWT( 116 ) / 'OLND            ',  136.00D0 /
      DATA CHEMISTRY_SPC( 117 ), SPECIES_MOLWT( 117 ) / 'APINP1          ',  230.00D0 /
      DATA CHEMISTRY_SPC( 118 ), SPECIES_MOLWT( 118 ) / 'APINP2          ',  294.00D0 /
      DATA CHEMISTRY_SPC( 119 ), SPECIES_MOLWT( 119 ) / 'LIMNP1          ',  230.00D0 /
      DATA CHEMISTRY_SPC( 120 ), SPECIES_MOLWT( 120 ) / 'LIMNP2          ',  294.00D0 /
      DATA CHEMISTRY_SPC( 121 ), SPECIES_MOLWT( 121 ) / 'ADCN            ',  155.00D0 /
      DATA CHEMISTRY_SPC( 122 ), SPECIES_MOLWT( 122 ) / 'VHOM            ',  250.00D0 /
      DATA CHEMISTRY_SPC( 123 ), SPECIES_MOLWT( 123 ) / 'VROCP4OXY2      ',  158.20D0 /
      DATA CHEMISTRY_SPC( 124 ), SPECIES_MOLWT( 124 ) / 'VROCN1OXY6      ',  190.20D0 /
      DATA CHEMISTRY_SPC( 125 ), SPECIES_MOLWT( 125 ) / 'FURANONE        ',  100.10D0 /
      DATA CHEMISTRY_SPC( 126 ), SPECIES_MOLWT( 126 ) / 'VROCP3OXY2      ',  186.30D0 /
      DATA CHEMISTRY_SPC( 127 ), SPECIES_MOLWT( 127 ) / 'VROCP0OXY4      ',  202.30D0 /
      DATA CHEMISTRY_SPC( 128 ), SPECIES_MOLWT( 128 ) / 'BAL2            ',   93.10D0 /
      DATA CHEMISTRY_SPC( 129 ), SPECIES_MOLWT( 129 ) / 'VELHOM          ',  402.00D0 /
      DATA CHEMISTRY_SPC( 130 ), SPECIES_MOLWT( 130 ) / 'VROCIOXY        ',  247.00D0 /
      DATA CHEMISTRY_SPC( 131 ), SPECIES_MOLWT( 131 ) / 'SLOWROC         ',   75.40D0 /
      DATA CHEMISTRY_SPC( 132 ), SPECIES_MOLWT( 132 ) / 'ACRO            ',   56.10D0 /
      DATA CHEMISTRY_SPC( 133 ), SPECIES_MOLWT( 133 ) / 'BDE13           ',   54.10D0 /
      DATA CHEMISTRY_SPC( 134 ), SPECIES_MOLWT( 134 ) / 'BDE13P          ',  103.00D0 /
      DATA CHEMISTRY_SPC( 135 ), SPECIES_MOLWT( 135 ) / 'FURAN           ',   96.10D0 /
      DATA CHEMISTRY_SPC( 136 ), SPECIES_MOLWT( 136 ) / 'FURANO2         ',  145.10D0 /
      DATA CHEMISTRY_SPC( 137 ), SPECIES_MOLWT( 137 ) / 'PROG            ',   76.10D0 /
      DATA CHEMISTRY_SPC( 138 ), SPECIES_MOLWT( 138 ) / 'SESQ            ',  204.40D0 /
      DATA CHEMISTRY_SPC( 139 ), SPECIES_MOLWT( 139 ) / 'SESQNRO2        ',  298.40D0 /
      DATA CHEMISTRY_SPC( 140 ), SPECIES_MOLWT( 140 ) / 'VROCN2OXY2      ',  282.40D0 /
      DATA CHEMISTRY_SPC( 141 ), SPECIES_MOLWT( 141 ) / 'SESQRO2         ',  253.40D0 /
      DATA CHEMISTRY_SPC( 142 ), SPECIES_MOLWT( 142 ) / 'VROCP0OXY2      ',  242.40D0 /
      DATA CHEMISTRY_SPC( 143 ), SPECIES_MOLWT( 143 ) / 'VROCP1OXY3      ',  202.30D0 /
      DATA CHEMISTRY_SPC( 144 ), SPECIES_MOLWT( 144 ) / 'AGLYOLIGJ       ',   66.40D0 /
      DATA CHEMISTRY_SPC( 145 ), SPECIES_MOLWT( 145 ) / 'IEPOXP          ',  118.10D0 /
      DATA CHEMISTRY_SPC( 146 ), SPECIES_MOLWT( 146 ) / 'AISO3NOSJ       ',  136.20D0 /
      DATA CHEMISTRY_SPC( 147 ), SPECIES_MOLWT( 147 ) / 'ASO4J           ',   96.00D0 /
      DATA CHEMISTRY_SPC( 148 ), SPECIES_MOLWT( 148 ) / 'AISO3OSJ        ',  216.20D0 /
      DATA CHEMISTRY_SPC( 149 ), SPECIES_MOLWT( 149 ) / 'AISO4J          ',  168.10D0 /
      DATA CHEMISTRY_SPC( 150 ), SPECIES_MOLWT( 150 ) / 'AISO5J          ',  150.10D0 /
      DATA CHEMISTRY_SPC( 151 ), SPECIES_MOLWT( 151 ) / 'VROCP6ALK       ',  184.37D0 /
      DATA CHEMISTRY_SPC( 152 ), SPECIES_MOLWT( 152 ) / 'VROCP6ALKP      ',  215.36D0 /
      DATA CHEMISTRY_SPC( 153 ), SPECIES_MOLWT( 153 ) / 'VROCP5ALK       ',  198.39D0 /
      DATA CHEMISTRY_SPC( 154 ), SPECIES_MOLWT( 154 ) / 'VROCP5ALKP      ',  229.38D0 /
      DATA CHEMISTRY_SPC( 155 ), SPECIES_MOLWT( 155 ) / 'VROCP4ALK       ',  254.50D0 /
      DATA CHEMISTRY_SPC( 156 ), SPECIES_MOLWT( 156 ) / 'VROCP4ALKP      ',  285.49D0 /
      DATA CHEMISTRY_SPC( 157 ), SPECIES_MOLWT( 157 ) / 'VROCP3ALK       ',  296.58D0 /
      DATA CHEMISTRY_SPC( 158 ), SPECIES_MOLWT( 158 ) / 'VROCP3ALKP      ',  327.57D0 /
      DATA CHEMISTRY_SPC( 159 ), SPECIES_MOLWT( 159 ) / 'VROCP2ALK       ',  338.66D0 /
      DATA CHEMISTRY_SPC( 160 ), SPECIES_MOLWT( 160 ) / 'VROCP2ALKP      ',  369.65D0 /
      DATA CHEMISTRY_SPC( 161 ), SPECIES_MOLWT( 161 ) / 'VROCP1ALK       ',  380.75D0 /
      DATA CHEMISTRY_SPC( 162 ), SPECIES_MOLWT( 162 ) / 'VROCP1ALKP      ',  411.74D0 /
      DATA CHEMISTRY_SPC( 163 ), SPECIES_MOLWT( 163 ) / 'HC10            ',  142.28D0 /
      DATA CHEMISTRY_SPC( 164 ), SPECIES_MOLWT( 164 ) / 'VROCP6ALKP2     ',  231.36D0 /
      DATA CHEMISTRY_SPC( 165 ), SPECIES_MOLWT( 165 ) / 'VROCP5ALKP2     ',  245.38D0 /
      DATA CHEMISTRY_SPC( 166 ), SPECIES_MOLWT( 166 ) / 'VROCP4ALKP2     ',  301.49D0 /
      DATA CHEMISTRY_SPC( 167 ), SPECIES_MOLWT( 167 ) / 'VROCP2OXY2      ',  200.30D0 /
      DATA CHEMISTRY_SPC( 168 ), SPECIES_MOLWT( 168 ) / 'VROCP3ALKP2     ',  343.57D0 /
      DATA CHEMISTRY_SPC( 169 ), SPECIES_MOLWT( 169 ) / 'VROCP1OXY1      ',  270.50D0 /
      DATA CHEMISTRY_SPC( 170 ), SPECIES_MOLWT( 170 ) / 'VROCP2ALKP2     ',  385.65D0 /
      DATA CHEMISTRY_SPC( 171 ), SPECIES_MOLWT( 171 ) / 'VROCP1ALKP2     ',  427.73D0 /
      DATA CHEMISTRY_SPC( 172 ), SPECIES_MOLWT( 172 ) / 'VROCN1OXY1      ',  312.50D0 /
      DATA CHEMISTRY_SPC( 173 ), SPECIES_MOLWT( 173 ) / 'HC10P2          ',  189.27D0 /
      DATA CHEMISTRY_SPC( 174 ), SPECIES_MOLWT( 174 ) / 'VROCP6ARO       ',  176.30D0 /
      DATA CHEMISTRY_SPC( 175 ), SPECIES_MOLWT( 175 ) / 'VROCP6AROP      ',  257.30D0 /
      DATA CHEMISTRY_SPC( 176 ), SPECIES_MOLWT( 176 ) / 'VROCN2OXY4      ',  232.30D0 /
      DATA CHEMISTRY_SPC( 177 ), SPECIES_MOLWT( 177 ) / 'VROCN1OXY3      ',  230.30D0 /
      DATA CHEMISTRY_SPC( 178 ), SPECIES_MOLWT( 178 ) / 'VROCP5ARO       ',  190.33D0 /
      DATA CHEMISTRY_SPC( 179 ), SPECIES_MOLWT( 179 ) / 'VROCP5AROP      ',  271.33D0 /
      DATA CHEMISTRY_SPC( 180 ), SPECIES_MOLWT( 180 ) / 'NAPH            ',  128.17D0 /
      DATA CHEMISTRY_SPC( 181 ), SPECIES_MOLWT( 181 ) / 'NAPHP           ',  209.17D0 /
      DATA CHEMISTRY_SPC( 182 ), SPECIES_MOLWT( 182 ) / 'VROCN2OXY8      ',  194.20D0 /
      DATA CHEMISTRY_SPC( 183 ), SPECIES_MOLWT( 183 ) / 'VROCP5OXY1      ',  170.30D0 /
      DATA CHEMISTRY_SPC( 184 ), SPECIES_MOLWT( 184 ) / 'VROCP6OXY1      ',  142.20D0 /
      DATA CHEMISTRY_SPC( 185 ), SPECIES_MOLWT( 185 ) / 'ECH4            ',   16.00D0 /
      DATA CHEMISTRY_SPC( 186 ), SPECIES_MOLWT( 186 ) / 'ATRPNJ          ',  215.00D0 /
      DATA CHEMISTRY_SPC( 187 ), SPECIES_MOLWT( 187 ) / 'AHOMJ           ',  250.00D0 /
      DATA CHEMISTRY_SPC( 188 ), SPECIES_MOLWT( 188 ) / 'AHONITJ         ',  265.00D0 /
      DATA CHEMISTRY_SPC( 189 ), SPECIES_MOLWT( 189 ) / 'STY             ',  104.00D0 /
      DATA CHEMISTRY_SPC( 190 ), SPECIES_MOLWT( 190 ) / 'STYP            ',  153.00D0 /
      DATA CHEMISTRY_SPC( 191 ), SPECIES_MOLWT( 191 ) / 'ANO3I           ',   62.00D0 /
      DATA CHEMISTRY_SPC( 192 ), SPECIES_MOLWT( 192 ) / 'ANO3J           ',   62.00D0 /
      DATA CHEMISTRY_SPC( 193 ), SPECIES_MOLWT( 193 ) / 'AMTN1J          ',  262.00D0 /
      DATA CHEMISTRY_SPC( 194 ), SPECIES_MOLWT( 194 ) / 'CL2             ',   71.00D0 /
      DATA CHEMISTRY_SPC( 195 ), SPECIES_MOLWT( 195 ) / 'CL              ',   35.50D0 /
      DATA CHEMISTRY_SPC( 196 ), SPECIES_MOLWT( 196 ) / 'CLO             ',   51.50D0 /
      DATA CHEMISTRY_SPC( 197 ), SPECIES_MOLWT( 197 ) / 'OCLO            ',   67.50D0 /
      DATA CHEMISTRY_SPC( 198 ), SPECIES_MOLWT( 198 ) / 'CL2O2           ',  103.00D0 /
      DATA CHEMISTRY_SPC( 199 ), SPECIES_MOLWT( 199 ) / 'CLOO            ',   67.50D0 /
      DATA CHEMISTRY_SPC( 200 ), SPECIES_MOLWT( 200 ) / 'HOCL            ',   52.50D0 /
      DATA CHEMISTRY_SPC( 201 ), SPECIES_MOLWT( 201 ) / 'CLNO            ',   65.50D0 /
      DATA CHEMISTRY_SPC( 202 ), SPECIES_MOLWT( 202 ) / 'CLNO2           ',   81.50D0 /
      DATA CHEMISTRY_SPC( 203 ), SPECIES_MOLWT( 203 ) / 'CLNO3           ',   97.50D0 /
      DATA CHEMISTRY_SPC( 204 ), SPECIES_MOLWT( 204 ) / 'HCOCL           ',   64.50D0 /
      DATA CHEMISTRY_SPC( 205 ), SPECIES_MOLWT( 205 ) / 'HCL             ',   36.00D0 /
      DATA CHEMISTRY_SPC( 206 ), SPECIES_MOLWT( 206 ) / 'BR2             ',  159.80D0 /
      DATA CHEMISTRY_SPC( 207 ), SPECIES_MOLWT( 207 ) / 'BR              ',   79.90D0 /
      DATA CHEMISTRY_SPC( 208 ), SPECIES_MOLWT( 208 ) / 'BRO             ',   95.90D0 /
      DATA CHEMISTRY_SPC( 209 ), SPECIES_MOLWT( 209 ) / 'OBRO            ',  111.90D0 /
      DATA CHEMISTRY_SPC( 210 ), SPECIES_MOLWT( 210 ) / 'HOBR            ',   96.90D0 /
      DATA CHEMISTRY_SPC( 211 ), SPECIES_MOLWT( 211 ) / 'BRNO            ',  109.90D0 /
      DATA CHEMISTRY_SPC( 212 ), SPECIES_MOLWT( 212 ) / 'BRNO2           ',  125.90D0 /
      DATA CHEMISTRY_SPC( 213 ), SPECIES_MOLWT( 213 ) / 'BRNO3           ',  141.90D0 /
      DATA CHEMISTRY_SPC( 214 ), SPECIES_MOLWT( 214 ) / 'CH2BR2          ',  173.80D0 /
      DATA CHEMISTRY_SPC( 215 ), SPECIES_MOLWT( 215 ) / 'CHBR3           ',  252.70D0 /
      DATA CHEMISTRY_SPC( 216 ), SPECIES_MOLWT( 216 ) / 'HCOBR           ',  108.90D0 /
      DATA CHEMISTRY_SPC( 217 ), SPECIES_MOLWT( 217 ) / 'HBR             ',   80.90D0 /
      DATA CHEMISTRY_SPC( 218 ), SPECIES_MOLWT( 218 ) / 'I2              ',  253.80D0 /
      DATA CHEMISTRY_SPC( 219 ), SPECIES_MOLWT( 219 ) / 'I               ',  126.90D0 /
      DATA CHEMISTRY_SPC( 220 ), SPECIES_MOLWT( 220 ) / 'IO              ',  142.90D0 /
      DATA CHEMISTRY_SPC( 221 ), SPECIES_MOLWT( 221 ) / 'OIO             ',  158.90D0 /
      DATA CHEMISTRY_SPC( 222 ), SPECIES_MOLWT( 222 ) / 'I2O2            ',  285.80D0 /
      DATA CHEMISTRY_SPC( 223 ), SPECIES_MOLWT( 223 ) / 'HOI             ',  143.90D0 /
      DATA CHEMISTRY_SPC( 224 ), SPECIES_MOLWT( 224 ) / 'HI              ',  127.90D0 /
      DATA CHEMISTRY_SPC( 225 ), SPECIES_MOLWT( 225 ) / 'INO             ',  156.90D0 /
      DATA CHEMISTRY_SPC( 226 ), SPECIES_MOLWT( 226 ) / 'INO2            ',  172.90D0 /
      DATA CHEMISTRY_SPC( 227 ), SPECIES_MOLWT( 227 ) / 'INO3            ',  188.90D0 /
      DATA CHEMISTRY_SPC( 228 ), SPECIES_MOLWT( 228 ) / 'CH3I            ',  141.90D0 /
      DATA CHEMISTRY_SPC( 229 ), SPECIES_MOLWT( 229 ) / 'CH2I2           ',  267.80D0 /
      DATA CHEMISTRY_SPC( 230 ), SPECIES_MOLWT( 230 ) / 'I2O3            ',  301.80D0 /
      DATA CHEMISTRY_SPC( 231 ), SPECIES_MOLWT( 231 ) / 'I2O4            ',  317.80D0 /
      DATA CHEMISTRY_SPC( 232 ), SPECIES_MOLWT( 232 ) / 'BRCL            ',  115.40D0 /
      DATA CHEMISTRY_SPC( 233 ), SPECIES_MOLWT( 233 ) / 'ICL             ',  162.40D0 /
      DATA CHEMISTRY_SPC( 234 ), SPECIES_MOLWT( 234 ) / 'IBR             ',  206.80D0 /
      DATA CHEMISTRY_SPC( 235 ), SPECIES_MOLWT( 235 ) / 'CH2IBR          ',  220.80D0 /
      DATA CHEMISTRY_SPC( 236 ), SPECIES_MOLWT( 236 ) / 'CH2ICL          ',  176.40D0 /
      DATA CHEMISTRY_SPC( 237 ), SPECIES_MOLWT( 237 ) / 'CHBR2CL         ',  208.30D0 /
      DATA CHEMISTRY_SPC( 238 ), SPECIES_MOLWT( 238 ) / 'CHBRCL2         ',  163.90D0 /
      DATA CHEMISTRY_SPC( 239 ), SPECIES_MOLWT( 239 ) / 'CH2BRCL         ',  129.40D0 /
      DATA CHEMISTRY_SPC( 240 ), SPECIES_MOLWT( 240 ) / 'NO2PIJ          ',   46.00D0 /
      DATA CHEMISTRY_SPC( 241 ), SPECIES_MOLWT( 241 ) / 'NO2PK           ',   46.00D0 /
      DATA CHEMISTRY_SPC( 242 ), SPECIES_MOLWT( 242 ) / 'ACLI            ',   35.50D0 /
      DATA CHEMISTRY_SPC( 243 ), SPECIES_MOLWT( 243 ) / 'ACLJ            ',   35.50D0 /
      DATA CHEMISTRY_SPC( 244 ), SPECIES_MOLWT( 244 ) / 'ACLK            ',   35.50D0 /
      DATA CHEMISTRY_SPC( 245 ), SPECIES_MOLWT( 245 ) / 'ABRJ            ',   79.90D0 /
      DATA CHEMISTRY_SPC( 246 ), SPECIES_MOLWT( 246 ) / 'AFEJ            ',   55.80D0 /
      DATA CHEMISTRY_SPC( 247 ), SPECIES_MOLWT( 247 ) / 'AMNJ            ',   54.90D0 /
      DATA CHEMISTRY_SPC( 248 ), SPECIES_MOLWT( 248 ) / 'AHMSJ           ',  111.10D0 /


      DATA CGRID_INDEX(   1 ), SPECIES_TYPE(   1 ), CONVERT_CONC(   1 ) /    1, 'GC', F /  ! O3
      DATA CGRID_INDEX(   2 ), SPECIES_TYPE(   2 ), CONVERT_CONC(   2 ) /    2, 'GC', F /  ! O3P
      DATA CGRID_INDEX(   3 ), SPECIES_TYPE(   3 ), CONVERT_CONC(   3 ) /    3, 'GC', F /  ! O1D
      DATA CGRID_INDEX(   4 ), SPECIES_TYPE(   4 ), CONVERT_CONC(   4 ) /    4, 'GC', F /  ! H2O2
      DATA CGRID_INDEX(   5 ), SPECIES_TYPE(   5 ), CONVERT_CONC(   5 ) /    5, 'GC', F /  ! HO
      DATA CGRID_INDEX(   6 ), SPECIES_TYPE(   6 ), CONVERT_CONC(   6 ) /    6, 'GC', F /  ! NO2
      DATA CGRID_INDEX(   7 ), SPECIES_TYPE(   7 ), CONVERT_CONC(   7 ) /    7, 'GC', F /  ! NO
      DATA CGRID_INDEX(   8 ), SPECIES_TYPE(   8 ), CONVERT_CONC(   8 ) /    8, 'GC', F /  ! NO3
      DATA CGRID_INDEX(   9 ), SPECIES_TYPE(   9 ), CONVERT_CONC(   9 ) /    9, 'GC', F /  ! HONO
      DATA CGRID_INDEX(  10 ), SPECIES_TYPE(  10 ), CONVERT_CONC(  10 ) /   10, 'GC', F /  ! HNO3
      DATA CGRID_INDEX(  11 ), SPECIES_TYPE(  11 ), CONVERT_CONC(  11 ) /   11, 'GC', F /  ! HNO4
      DATA CGRID_INDEX(  12 ), SPECIES_TYPE(  12 ), CONVERT_CONC(  12 ) /   12, 'GC', F /  ! HO2
      DATA CGRID_INDEX(  13 ), SPECIES_TYPE(  13 ), CONVERT_CONC(  13 ) /   13, 'GC', F /  ! HCHO
      DATA CGRID_INDEX(  14 ), SPECIES_TYPE(  14 ), CONVERT_CONC(  14 ) /   14, 'GC', F /  ! CO
      DATA CGRID_INDEX(  15 ), SPECIES_TYPE(  15 ), CONVERT_CONC(  15 ) /   15, 'GC', F /  ! ACD
      DATA CGRID_INDEX(  16 ), SPECIES_TYPE(  16 ), CONVERT_CONC(  16 ) /   16, 'GC', F /  ! MO2
      DATA CGRID_INDEX(  17 ), SPECIES_TYPE(  17 ), CONVERT_CONC(  17 ) /   17, 'GC', F /  ! ALD
      DATA CGRID_INDEX(  18 ), SPECIES_TYPE(  18 ), CONVERT_CONC(  18 ) /   18, 'GC', F /  ! ETHP
      DATA CGRID_INDEX(  19 ), SPECIES_TYPE(  19 ), CONVERT_CONC(  19 ) /   19, 'GC', F /  ! ACT
      DATA CGRID_INDEX(  20 ), SPECIES_TYPE(  20 ), CONVERT_CONC(  20 ) /   20, 'GC', F /  ! ACO3
      DATA CGRID_INDEX(  21 ), SPECIES_TYPE(  21 ), CONVERT_CONC(  21 ) /   21, 'GC', F /  ! UALD
      DATA CGRID_INDEX(  22 ), SPECIES_TYPE(  22 ), CONVERT_CONC(  22 ) /   22, 'GC', F /  ! KET
      DATA CGRID_INDEX(  23 ), SPECIES_TYPE(  23 ), CONVERT_CONC(  23 ) /   81, 'GC', F /  ! PINAL
      DATA CGRID_INDEX(  24 ), SPECIES_TYPE(  24 ), CONVERT_CONC(  24 ) /   82, 'GC', F /  ! PINALP
      DATA CGRID_INDEX(  25 ), SPECIES_TYPE(  25 ), CONVERT_CONC(  25 ) /   88, 'GC', F /  ! LIMAL
      DATA CGRID_INDEX(  26 ), SPECIES_TYPE(  26 ), CONVERT_CONC(  26 ) /   89, 'GC', F /  ! LIMALP
      DATA CGRID_INDEX(  27 ), SPECIES_TYPE(  27 ), CONVERT_CONC(  27 ) /   23, 'GC', F /  ! MEK
      DATA CGRID_INDEX(  28 ), SPECIES_TYPE(  28 ), CONVERT_CONC(  28 ) /   24, 'GC', F /  ! HKET
      DATA CGRID_INDEX(  29 ), SPECIES_TYPE(  29 ), CONVERT_CONC(  29 ) /   25, 'GC', F /  ! MACR
      DATA CGRID_INDEX(  30 ), SPECIES_TYPE(  30 ), CONVERT_CONC(  30 ) /   26, 'GC', F /  ! MACP
      DATA CGRID_INDEX(  31 ), SPECIES_TYPE(  31 ), CONVERT_CONC(  31 ) /   27, 'GC', F /  ! XO2
      DATA CGRID_INDEX(  32 ), SPECIES_TYPE(  32 ), CONVERT_CONC(  32 ) /   28, 'GC', F /  ! MVK
      DATA CGRID_INDEX(  33 ), SPECIES_TYPE(  33 ), CONVERT_CONC(  33 ) /   29, 'GC', F /  ! GLY
      DATA CGRID_INDEX(  34 ), SPECIES_TYPE(  34 ), CONVERT_CONC(  34 ) /   30, 'GC', F /  ! MGLY
      DATA CGRID_INDEX(  35 ), SPECIES_TYPE(  35 ), CONVERT_CONC(  35 ) /   31, 'GC', F /  ! DCB1
      DATA CGRID_INDEX(  36 ), SPECIES_TYPE(  36 ), CONVERT_CONC(  36 ) /   32, 'GC', F /  ! DCB2
      DATA CGRID_INDEX(  37 ), SPECIES_TYPE(  37 ), CONVERT_CONC(  37 ) /   33, 'GC', F /  ! BALD
      DATA CGRID_INDEX(  38 ), SPECIES_TYPE(  38 ), CONVERT_CONC(  38 ) /   62, 'GC', F /  ! BEN
      DATA CGRID_INDEX(  39 ), SPECIES_TYPE(  39 ), CONVERT_CONC(  39 ) /  126, 'GC', F /  ! BAL1
      DATA CGRID_INDEX(  40 ), SPECIES_TYPE(  40 ), CONVERT_CONC(  40 ) /   35, 'GC', F /  ! OP1
      DATA CGRID_INDEX(  41 ), SPECIES_TYPE(  41 ), CONVERT_CONC(  41 ) /   36, 'GC', F /  ! OP2
      DATA CGRID_INDEX(  42 ), SPECIES_TYPE(  42 ), CONVERT_CONC(  42 ) /   37, 'GC', F /  ! OPB
      DATA CGRID_INDEX(  43 ), SPECIES_TYPE(  43 ), CONVERT_CONC(  43 ) /   38, 'GC', F /  ! VOP3
      DATA CGRID_INDEX(  44 ), SPECIES_TYPE(  44 ), CONVERT_CONC(  44 ) /   39, 'GC', F /  ! PAA
      DATA CGRID_INDEX(  45 ), SPECIES_TYPE(  45 ), CONVERT_CONC(  45 ) /   40, 'GC', F /  ! ONIT
      DATA CGRID_INDEX(  46 ), SPECIES_TYPE(  46 ), CONVERT_CONC(  46 ) /   41, 'GC', F /  ! PAN
      DATA CGRID_INDEX(  47 ), SPECIES_TYPE(  47 ), CONVERT_CONC(  47 ) /  183, 'GC', F /  ! CO2
      DATA CGRID_INDEX(  48 ), SPECIES_TYPE(  48 ), CONVERT_CONC(  48 ) /  114, 'GC', F /  ! PPN
      DATA CGRID_INDEX(  49 ), SPECIES_TYPE(  49 ), CONVERT_CONC(  49 ) /   92, 'GC', F /  ! RCO3
      DATA CGRID_INDEX(  50 ), SPECIES_TYPE(  50 ), CONVERT_CONC(  50 ) /   48, 'GC', F /  ! HC3P
      DATA CGRID_INDEX(  51 ), SPECIES_TYPE(  51 ), CONVERT_CONC(  51 ) /  120, 'GC', F /  ! VTRPN
      DATA CGRID_INDEX(  52 ), SPECIES_TYPE(  52 ), CONVERT_CONC(  52 ) /  121, 'GC', F /  ! VHONIT
      DATA CGRID_INDEX(  53 ), SPECIES_TYPE(  53 ), CONVERT_CONC(  53 ) /   42, 'GC', F /  ! N2O5
      DATA CGRID_INDEX(  54 ), SPECIES_TYPE(  54 ), CONVERT_CONC(  54 ) /   43, 'GC', F /  ! SO2
      DATA CGRID_INDEX(  55 ), SPECIES_TYPE(  55 ), CONVERT_CONC(  55 ) /   44, 'GC', F /  ! SULF
      DATA CGRID_INDEX(  56 ), SPECIES_TYPE(  56 ), CONVERT_CONC(  56 ) /   45, 'GC', F /  ! SULRXN
      DATA CGRID_INDEX(  57 ), SPECIES_TYPE(  57 ), CONVERT_CONC(  57 ) /   46, 'GC', F /  ! ETH
      DATA CGRID_INDEX(  58 ), SPECIES_TYPE(  58 ), CONVERT_CONC(  58 ) /   47, 'GC', F /  ! HC3
      DATA CGRID_INDEX(  59 ), SPECIES_TYPE(  59 ), CONVERT_CONC(  59 ) /  243, 'AE', T /  ! ASOATJ
      DATA CGRID_INDEX(  60 ), SPECIES_TYPE(  60 ), CONVERT_CONC(  60 ) /   49, 'GC', F /  ! HC5
      DATA CGRID_INDEX(  61 ), SPECIES_TYPE(  61 ), CONVERT_CONC(  61 ) /   50, 'GC', F /  ! HC5P
      DATA CGRID_INDEX(  62 ), SPECIES_TYPE(  62 ), CONVERT_CONC(  62 ) /   54, 'GC', F /  ! ETE
      DATA CGRID_INDEX(  63 ), SPECIES_TYPE(  63 ), CONVERT_CONC(  63 ) /   55, 'GC', F /  ! ETEP
      DATA CGRID_INDEX(  64 ), SPECIES_TYPE(  64 ), CONVERT_CONC(  64 ) /   56, 'GC', F /  ! OLT
      DATA CGRID_INDEX(  65 ), SPECIES_TYPE(  65 ), CONVERT_CONC(  65 ) /   57, 'GC', F /  ! OLTP
      DATA CGRID_INDEX(  66 ), SPECIES_TYPE(  66 ), CONVERT_CONC(  66 ) /   58, 'GC', F /  ! OLI
      DATA CGRID_INDEX(  67 ), SPECIES_TYPE(  67 ), CONVERT_CONC(  67 ) /   59, 'GC', F /  ! OLIP
      DATA CGRID_INDEX(  68 ), SPECIES_TYPE(  68 ), CONVERT_CONC(  68 ) /   60, 'GC', F /  ! ACE
      DATA CGRID_INDEX(  69 ), SPECIES_TYPE(  69 ), CONVERT_CONC(  69 ) /   61, 'GC', F /  ! ORA1
      DATA CGRID_INDEX(  70 ), SPECIES_TYPE(  70 ), CONVERT_CONC(  70 ) /   63, 'GC', F /  ! BENP
      DATA CGRID_INDEX(  71 ), SPECIES_TYPE(  71 ), CONVERT_CONC(  71 ) /   64, 'GC', F /  ! PHEN
      DATA CGRID_INDEX(  72 ), SPECIES_TYPE(  72 ), CONVERT_CONC(  72 ) /   65, 'GC', F /  ! TOL
      DATA CGRID_INDEX(  73 ), SPECIES_TYPE(  73 ), CONVERT_CONC(  73 ) /   70, 'GC', F /  ! TOLP
      DATA CGRID_INDEX(  74 ), SPECIES_TYPE(  74 ), CONVERT_CONC(  74 ) /   66, 'GC', F /  ! CSL
      DATA CGRID_INDEX(  75 ), SPECIES_TYPE(  75 ), CONVERT_CONC(  75 ) /   67, 'GC', F /  ! XYL
      DATA CGRID_INDEX(  76 ), SPECIES_TYPE(  76 ), CONVERT_CONC(  76 ) /   71, 'GC', F /  ! XYLP
      DATA CGRID_INDEX(  77 ), SPECIES_TYPE(  77 ), CONVERT_CONC(  77 ) /   68, 'GC', F /  ! EBZ
      DATA CGRID_INDEX(  78 ), SPECIES_TYPE(  78 ), CONVERT_CONC(  78 ) /   72, 'GC', F /  ! EBZP
      DATA CGRID_INDEX(  79 ), SPECIES_TYPE(  79 ), CONVERT_CONC(  79 ) /   74, 'GC', F /  ! ISO
      DATA CGRID_INDEX(  80 ), SPECIES_TYPE(  80 ), CONVERT_CONC(  80 ) /  118, 'GC', F /  ! ISON
      DATA CGRID_INDEX(  81 ), SPECIES_TYPE(  81 ), CONVERT_CONC(  81 ) /  117, 'GC', F /  ! ISONP
      DATA CGRID_INDEX(  82 ), SPECIES_TYPE(  82 ), CONVERT_CONC(  82 ) /   75, 'GC', F /  ! ISOP
      DATA CGRID_INDEX(  83 ), SPECIES_TYPE(  83 ), CONVERT_CONC(  83 ) /  108, 'GC', F /  ! ISHP
      DATA CGRID_INDEX(  84 ), SPECIES_TYPE(  84 ), CONVERT_CONC(  84 ) /  109, 'GC', F /  ! IEPOX
      DATA CGRID_INDEX(  85 ), SPECIES_TYPE(  85 ), CONVERT_CONC(  85 ) /  119, 'GC', F /  ! IPX
      DATA CGRID_INDEX(  86 ), SPECIES_TYPE(  86 ), CONVERT_CONC(  86 ) /  116, 'GC', F /  ! INALD
      DATA CGRID_INDEX(  87 ), SPECIES_TYPE(  87 ), CONVERT_CONC(  87 ) /  106, 'GC', F /  ! ROH
      DATA CGRID_INDEX(  88 ), SPECIES_TYPE(  88 ), CONVERT_CONC(  88 ) /   76, 'GC', F /  ! API
      DATA CGRID_INDEX(  89 ), SPECIES_TYPE(  89 ), CONVERT_CONC(  89 ) /   77, 'GC', F /  ! APIP1
      DATA CGRID_INDEX(  90 ), SPECIES_TYPE(  90 ), CONVERT_CONC(  90 ) /   78, 'GC', F /  ! APIP2
      DATA CGRID_INDEX(  91 ), SPECIES_TYPE(  91 ), CONVERT_CONC(  91 ) /   83, 'GC', F /  ! LIM
      DATA CGRID_INDEX(  92 ), SPECIES_TYPE(  92 ), CONVERT_CONC(  92 ) /   84, 'GC', F /  ! LIMP1
      DATA CGRID_INDEX(  93 ), SPECIES_TYPE(  93 ), CONVERT_CONC(  93 ) /   85, 'GC', F /  ! LIMP2
      DATA CGRID_INDEX(  94 ), SPECIES_TYPE(  94 ), CONVERT_CONC(  94 ) /   93, 'GC', F /  ! ACTP
      DATA CGRID_INDEX(  95 ), SPECIES_TYPE(  95 ), CONVERT_CONC(  95 ) /   94, 'GC', F /  ! MEKP
      DATA CGRID_INDEX(  96 ), SPECIES_TYPE(  96 ), CONVERT_CONC(  96 ) /   95, 'GC', F /  ! KETP
      DATA CGRID_INDEX(  97 ), SPECIES_TYPE(  97 ), CONVERT_CONC(  97 ) /   96, 'GC', F /  ! MCP
      DATA CGRID_INDEX(  98 ), SPECIES_TYPE(  98 ), CONVERT_CONC(  98 ) /   97, 'GC', F /  ! MVKP
      DATA CGRID_INDEX(  99 ), SPECIES_TYPE(  99 ), CONVERT_CONC(  99 ) /   98, 'GC', F /  ! UALP
      DATA CGRID_INDEX( 100 ), SPECIES_TYPE( 100 ), CONVERT_CONC( 100 ) /   99, 'GC', F /  ! DCB3
      DATA CGRID_INDEX( 101 ), SPECIES_TYPE( 101 ), CONVERT_CONC( 101 ) /  100, 'GC', F /  ! BALP
      DATA CGRID_INDEX( 102 ), SPECIES_TYPE( 102 ), CONVERT_CONC( 102 ) /  101, 'GC', F /  ! ADDC
      DATA CGRID_INDEX( 103 ), SPECIES_TYPE( 103 ), CONVERT_CONC( 103 ) /   34, 'GC', F /  ! CHO
      DATA CGRID_INDEX( 104 ), SPECIES_TYPE( 104 ), CONVERT_CONC( 104 ) /  102, 'GC', F /  ! MCT
      DATA CGRID_INDEX( 105 ), SPECIES_TYPE( 105 ), CONVERT_CONC( 105 ) /  103, 'GC', F /  ! MCTO
      DATA CGRID_INDEX( 106 ), SPECIES_TYPE( 106 ), CONVERT_CONC( 106 ) /  104, 'GC', F /  ! MOH
      DATA CGRID_INDEX( 107 ), SPECIES_TYPE( 107 ), CONVERT_CONC( 107 ) /  105, 'GC', F /  ! EOH
      DATA CGRID_INDEX( 108 ), SPECIES_TYPE( 108 ), CONVERT_CONC( 108 ) /  107, 'GC', F /  ! ETEG
      DATA CGRID_INDEX( 109 ), SPECIES_TYPE( 109 ), CONVERT_CONC( 109 ) /   52, 'GC', F /  ! HC10P
      DATA CGRID_INDEX( 110 ), SPECIES_TYPE( 110 ), CONVERT_CONC( 110 ) /  111, 'GC', F /  ! MAHP
      DATA CGRID_INDEX( 111 ), SPECIES_TYPE( 111 ), CONVERT_CONC( 111 ) /  112, 'GC', F /  ! ORA2
      DATA CGRID_INDEX( 112 ), SPECIES_TYPE( 112 ), CONVERT_CONC( 112 ) /  113, 'GC', F /  ! ORAP
      DATA CGRID_INDEX( 113 ), SPECIES_TYPE( 113 ), CONVERT_CONC( 113 ) /  115, 'GC', F /  ! MPAN
      DATA CGRID_INDEX( 114 ), SPECIES_TYPE( 114 ), CONVERT_CONC( 114 ) /  122, 'GC', F /  ! MCTP
      DATA CGRID_INDEX( 115 ), SPECIES_TYPE( 115 ), CONVERT_CONC( 115 ) /  123, 'GC', F /  ! OLNN
      DATA CGRID_INDEX( 116 ), SPECIES_TYPE( 116 ), CONVERT_CONC( 116 ) /  124, 'GC', F /  ! OLND
      DATA CGRID_INDEX( 117 ), SPECIES_TYPE( 117 ), CONVERT_CONC( 117 ) /   79, 'GC', F /  ! APINP1
      DATA CGRID_INDEX( 118 ), SPECIES_TYPE( 118 ), CONVERT_CONC( 118 ) /   80, 'GC', F /  ! APINP2
      DATA CGRID_INDEX( 119 ), SPECIES_TYPE( 119 ), CONVERT_CONC( 119 ) /   86, 'GC', F /  ! LIMNP1
      DATA CGRID_INDEX( 120 ), SPECIES_TYPE( 120 ), CONVERT_CONC( 120 ) /   87, 'GC', F /  ! LIMNP2
      DATA CGRID_INDEX( 121 ), SPECIES_TYPE( 121 ), CONVERT_CONC( 121 ) /  125, 'GC', F /  ! ADCN
      DATA CGRID_INDEX( 122 ), SPECIES_TYPE( 122 ), CONVERT_CONC( 122 ) /   90, 'GC', F /  ! VHOM
      DATA CGRID_INDEX( 123 ), SPECIES_TYPE( 123 ), CONVERT_CONC( 123 ) /  179, 'GC', F /  ! VROCP4OXY2
      DATA CGRID_INDEX( 124 ), SPECIES_TYPE( 124 ), CONVERT_CONC( 124 ) /  172, 'GC', F /  ! VROCN1OXY6
      DATA CGRID_INDEX( 125 ), SPECIES_TYPE( 125 ), CONVERT_CONC( 125 ) /  134, 'GC', F /  ! FURANONE
      DATA CGRID_INDEX( 126 ), SPECIES_TYPE( 126 ), CONVERT_CONC( 126 ) /  178, 'GC', F /  ! VROCP3OXY2
      DATA CGRID_INDEX( 127 ), SPECIES_TYPE( 127 ), CONVERT_CONC( 127 ) /  174, 'GC', F /  ! VROCP0OXY4
      DATA CGRID_INDEX( 128 ), SPECIES_TYPE( 128 ), CONVERT_CONC( 128 ) /  127, 'GC', F /  ! BAL2
      DATA CGRID_INDEX( 129 ), SPECIES_TYPE( 129 ), CONVERT_CONC( 129 ) /   91, 'GC', F /  ! VELHOM
      DATA CGRID_INDEX( 130 ), SPECIES_TYPE( 130 ), CONVERT_CONC( 130 ) /  135, 'GC', F /  ! VROCIOXY
      DATA CGRID_INDEX( 131 ), SPECIES_TYPE( 131 ), CONVERT_CONC( 131 ) /  136, 'GC', F /  ! SLOWROC
      DATA CGRID_INDEX( 132 ), SPECIES_TYPE( 132 ), CONVERT_CONC( 132 ) /  128, 'GC', F /  ! ACRO
      DATA CGRID_INDEX( 133 ), SPECIES_TYPE( 133 ), CONVERT_CONC( 133 ) /  129, 'GC', F /  ! BDE13
      DATA CGRID_INDEX( 134 ), SPECIES_TYPE( 134 ), CONVERT_CONC( 134 ) /  130, 'GC', F /  ! BDE13P
      DATA CGRID_INDEX( 135 ), SPECIES_TYPE( 135 ), CONVERT_CONC( 135 ) /  132, 'GC', F /  ! FURAN
      DATA CGRID_INDEX( 136 ), SPECIES_TYPE( 136 ), CONVERT_CONC( 136 ) /  133, 'GC', F /  ! FURANO2
      DATA CGRID_INDEX( 137 ), SPECIES_TYPE( 137 ), CONVERT_CONC( 137 ) /  131, 'GC', F /  ! PROG
      DATA CGRID_INDEX( 138 ), SPECIES_TYPE( 138 ), CONVERT_CONC( 138 ) /  137, 'GC', F /  ! SESQ
      DATA CGRID_INDEX( 139 ), SPECIES_TYPE( 139 ), CONVERT_CONC( 139 ) /  139, 'GC', F /  ! SESQNRO2
      DATA CGRID_INDEX( 140 ), SPECIES_TYPE( 140 ), CONVERT_CONC( 140 ) /  167, 'GC', F /  ! VROCN2OXY2
      DATA CGRID_INDEX( 141 ), SPECIES_TYPE( 141 ), CONVERT_CONC( 141 ) /  138, 'GC', F /  ! SESQRO2
      DATA CGRID_INDEX( 142 ), SPECIES_TYPE( 142 ), CONVERT_CONC( 142 ) /  173, 'GC', F /  ! VROCP0OXY2
      DATA CGRID_INDEX( 143 ), SPECIES_TYPE( 143 ), CONVERT_CONC( 143 ) /  176, 'GC', F /  ! VROCP1OXY3
      DATA CGRID_INDEX( 144 ), SPECIES_TYPE( 144 ), CONVERT_CONC( 144 ) /  244, 'AE', T /  ! AGLYOLIGJ
      DATA CGRID_INDEX( 145 ), SPECIES_TYPE( 145 ), CONVERT_CONC( 145 ) /  110, 'GC', F /  ! IEPOXP
      DATA CGRID_INDEX( 146 ), SPECIES_TYPE( 146 ), CONVERT_CONC( 146 ) /  318, 'AE', T /  ! AISO3NOSJ
      DATA CGRID_INDEX( 147 ), SPECIES_TYPE( 147 ), CONVERT_CONC( 147 ) /  235, 'AE', T /  ! ASO4J
      DATA CGRID_INDEX( 148 ), SPECIES_TYPE( 148 ), CONVERT_CONC( 148 ) /  319, 'AE', T /  ! AISO3OSJ
      DATA CGRID_INDEX( 149 ), SPECIES_TYPE( 149 ), CONVERT_CONC( 149 ) /  320, 'AE', T /  ! AISO4J
      DATA CGRID_INDEX( 150 ), SPECIES_TYPE( 150 ), CONVERT_CONC( 150 ) /  321, 'AE', T /  ! AISO5J
      DATA CGRID_INDEX( 151 ), SPECIES_TYPE( 151 ), CONVERT_CONC( 151 ) /  154, 'GC', F /  ! VROCP6ALK
      DATA CGRID_INDEX( 152 ), SPECIES_TYPE( 152 ), CONVERT_CONC( 152 ) /  160, 'GC', F /  ! VROCP6ALKP
      DATA CGRID_INDEX( 153 ), SPECIES_TYPE( 153 ), CONVERT_CONC( 153 ) /  153, 'GC', F /  ! VROCP5ALK
      DATA CGRID_INDEX( 154 ), SPECIES_TYPE( 154 ), CONVERT_CONC( 154 ) /  159, 'GC', F /  ! VROCP5ALKP
      DATA CGRID_INDEX( 155 ), SPECIES_TYPE( 155 ), CONVERT_CONC( 155 ) /  152, 'GC', F /  ! VROCP4ALK
      DATA CGRID_INDEX( 156 ), SPECIES_TYPE( 156 ), CONVERT_CONC( 156 ) /  158, 'GC', F /  ! VROCP4ALKP
      DATA CGRID_INDEX( 157 ), SPECIES_TYPE( 157 ), CONVERT_CONC( 157 ) /  151, 'GC', F /  ! VROCP3ALK
      DATA CGRID_INDEX( 158 ), SPECIES_TYPE( 158 ), CONVERT_CONC( 158 ) /  157, 'GC', F /  ! VROCP3ALKP
      DATA CGRID_INDEX( 159 ), SPECIES_TYPE( 159 ), CONVERT_CONC( 159 ) /  150, 'GC', F /  ! VROCP2ALK
      DATA CGRID_INDEX( 160 ), SPECIES_TYPE( 160 ), CONVERT_CONC( 160 ) /  156, 'GC', F /  ! VROCP2ALKP
      DATA CGRID_INDEX( 161 ), SPECIES_TYPE( 161 ), CONVERT_CONC( 161 ) /  149, 'GC', F /  ! VROCP1ALK
      DATA CGRID_INDEX( 162 ), SPECIES_TYPE( 162 ), CONVERT_CONC( 162 ) /  155, 'GC', F /  ! VROCP1ALKP
      DATA CGRID_INDEX( 163 ), SPECIES_TYPE( 163 ), CONVERT_CONC( 163 ) /   51, 'GC', F /  ! HC10
      DATA CGRID_INDEX( 164 ), SPECIES_TYPE( 164 ), CONVERT_CONC( 164 ) /  166, 'GC', F /  ! VROCP6ALKP2
      DATA CGRID_INDEX( 165 ), SPECIES_TYPE( 165 ), CONVERT_CONC( 165 ) /  165, 'GC', F /  ! VROCP5ALKP2
      DATA CGRID_INDEX( 166 ), SPECIES_TYPE( 166 ), CONVERT_CONC( 166 ) /  164, 'GC', F /  ! VROCP4ALKP2
      DATA CGRID_INDEX( 167 ), SPECIES_TYPE( 167 ), CONVERT_CONC( 167 ) /  177, 'GC', F /  ! VROCP2OXY2
      DATA CGRID_INDEX( 168 ), SPECIES_TYPE( 168 ), CONVERT_CONC( 168 ) /  163, 'GC', F /  ! VROCP3ALKP2
      DATA CGRID_INDEX( 169 ), SPECIES_TYPE( 169 ), CONVERT_CONC( 169 ) /  175, 'GC', F /  ! VROCP1OXY1
      DATA CGRID_INDEX( 170 ), SPECIES_TYPE( 170 ), CONVERT_CONC( 170 ) /  162, 'GC', F /  ! VROCP2ALKP2
      DATA CGRID_INDEX( 171 ), SPECIES_TYPE( 171 ), CONVERT_CONC( 171 ) /  161, 'GC', F /  ! VROCP1ALKP2
      DATA CGRID_INDEX( 172 ), SPECIES_TYPE( 172 ), CONVERT_CONC( 172 ) /  170, 'GC', F /  ! VROCN1OXY1
      DATA CGRID_INDEX( 173 ), SPECIES_TYPE( 173 ), CONVERT_CONC( 173 ) /   53, 'GC', F /  ! HC10P2
      DATA CGRID_INDEX( 174 ), SPECIES_TYPE( 174 ), CONVERT_CONC( 174 ) /  143, 'GC', F /  ! VROCP6ARO
      DATA CGRID_INDEX( 175 ), SPECIES_TYPE( 175 ), CONVERT_CONC( 175 ) /  145, 'GC', F /  ! VROCP6AROP
      DATA CGRID_INDEX( 176 ), SPECIES_TYPE( 176 ), CONVERT_CONC( 176 ) /  168, 'GC', F /  ! VROCN2OXY4
      DATA CGRID_INDEX( 177 ), SPECIES_TYPE( 177 ), CONVERT_CONC( 177 ) /  171, 'GC', F /  ! VROCN1OXY3
      DATA CGRID_INDEX( 178 ), SPECIES_TYPE( 178 ), CONVERT_CONC( 178 ) /  142, 'GC', F /  ! VROCP5ARO
      DATA CGRID_INDEX( 179 ), SPECIES_TYPE( 179 ), CONVERT_CONC( 179 ) /  144, 'GC', F /  ! VROCP5AROP
      DATA CGRID_INDEX( 180 ), SPECIES_TYPE( 180 ), CONVERT_CONC( 180 ) /  140, 'GC', F /  ! NAPH
      DATA CGRID_INDEX( 181 ), SPECIES_TYPE( 181 ), CONVERT_CONC( 181 ) /  141, 'GC', F /  ! NAPHP
      DATA CGRID_INDEX( 182 ), SPECIES_TYPE( 182 ), CONVERT_CONC( 182 ) /  169, 'GC', F /  ! VROCN2OXY8
      DATA CGRID_INDEX( 183 ), SPECIES_TYPE( 183 ), CONVERT_CONC( 183 ) /  180, 'GC', F /  ! VROCP5OXY1
      DATA CGRID_INDEX( 184 ), SPECIES_TYPE( 184 ), CONVERT_CONC( 184 ) /  181, 'GC', F /  ! VROCP6OXY1
      DATA CGRID_INDEX( 185 ), SPECIES_TYPE( 185 ), CONVERT_CONC( 185 ) /  182, 'GC', F /  ! ECH4
      DATA CGRID_INDEX( 186 ), SPECIES_TYPE( 186 ), CONVERT_CONC( 186 ) /  322, 'AE', T /  ! ATRPNJ
      DATA CGRID_INDEX( 187 ), SPECIES_TYPE( 187 ), CONVERT_CONC( 187 ) /  245, 'AE', T /  ! AHOMJ
      DATA CGRID_INDEX( 188 ), SPECIES_TYPE( 188 ), CONVERT_CONC( 188 ) /  323, 'AE', T /  ! AHONITJ
      DATA CGRID_INDEX( 189 ), SPECIES_TYPE( 189 ), CONVERT_CONC( 189 ) /   69, 'GC', F /  ! STY
      DATA CGRID_INDEX( 190 ), SPECIES_TYPE( 190 ), CONVERT_CONC( 190 ) /   73, 'GC', F /  ! STYP
      DATA CGRID_INDEX( 191 ), SPECIES_TYPE( 191 ), CONVERT_CONC( 191 ) /  240, 'AE', T /  ! ANO3I
      DATA CGRID_INDEX( 192 ), SPECIES_TYPE( 192 ), CONVERT_CONC( 192 ) /  241, 'AE', T /  ! ANO3J
      DATA CGRID_INDEX( 193 ), SPECIES_TYPE( 193 ), CONVERT_CONC( 193 ) /  325, 'AE', T /  ! AMTN1J
      DATA CGRID_INDEX( 194 ), SPECIES_TYPE( 194 ), CONVERT_CONC( 194 ) /  185, 'GC', F /  ! CL2
      DATA CGRID_INDEX( 195 ), SPECIES_TYPE( 195 ), CONVERT_CONC( 195 ) /  186, 'GC', F /  ! CL
      DATA CGRID_INDEX( 196 ), SPECIES_TYPE( 196 ), CONVERT_CONC( 196 ) /  187, 'GC', F /  ! CLO
      DATA CGRID_INDEX( 197 ), SPECIES_TYPE( 197 ), CONVERT_CONC( 197 ) /  188, 'GC', F /  ! OCLO
      DATA CGRID_INDEX( 198 ), SPECIES_TYPE( 198 ), CONVERT_CONC( 198 ) /  189, 'GC', F /  ! CL2O2
      DATA CGRID_INDEX( 199 ), SPECIES_TYPE( 199 ), CONVERT_CONC( 199 ) /  196, 'GC', F /  ! CLOO
      DATA CGRID_INDEX( 200 ), SPECIES_TYPE( 200 ), CONVERT_CONC( 200 ) /  190, 'GC', F /  ! HOCL
      DATA CGRID_INDEX( 201 ), SPECIES_TYPE( 201 ), CONVERT_CONC( 201 ) /  191, 'GC', F /  ! CLNO
      DATA CGRID_INDEX( 202 ), SPECIES_TYPE( 202 ), CONVERT_CONC( 202 ) /  192, 'GC', F /  ! CLNO2
      DATA CGRID_INDEX( 203 ), SPECIES_TYPE( 203 ), CONVERT_CONC( 203 ) /  193, 'GC', F /  ! CLNO3
      DATA CGRID_INDEX( 204 ), SPECIES_TYPE( 204 ), CONVERT_CONC( 204 ) /  194, 'GC', F /  ! HCOCL
      DATA CGRID_INDEX( 205 ), SPECIES_TYPE( 205 ), CONVERT_CONC( 205 ) /  195, 'GC', F /  ! HCL
      DATA CGRID_INDEX( 206 ), SPECIES_TYPE( 206 ), CONVERT_CONC( 206 ) /  197, 'GC', F /  ! BR2
      DATA CGRID_INDEX( 207 ), SPECIES_TYPE( 207 ), CONVERT_CONC( 207 ) /  198, 'GC', F /  ! BR
      DATA CGRID_INDEX( 208 ), SPECIES_TYPE( 208 ), CONVERT_CONC( 208 ) /  199, 'GC', F /  ! BRO
      DATA CGRID_INDEX( 209 ), SPECIES_TYPE( 209 ), CONVERT_CONC( 209 ) /  200, 'GC', F /  ! OBRO
      DATA CGRID_INDEX( 210 ), SPECIES_TYPE( 210 ), CONVERT_CONC( 210 ) /  201, 'GC', F /  ! HOBR
      DATA CGRID_INDEX( 211 ), SPECIES_TYPE( 211 ), CONVERT_CONC( 211 ) /  202, 'GC', F /  ! BRNO
      DATA CGRID_INDEX( 212 ), SPECIES_TYPE( 212 ), CONVERT_CONC( 212 ) /  203, 'GC', F /  ! BRNO2
      DATA CGRID_INDEX( 213 ), SPECIES_TYPE( 213 ), CONVERT_CONC( 213 ) /  204, 'GC', F /  ! BRNO3
      DATA CGRID_INDEX( 214 ), SPECIES_TYPE( 214 ), CONVERT_CONC( 214 ) /  205, 'GC', F /  ! CH2BR2
      DATA CGRID_INDEX( 215 ), SPECIES_TYPE( 215 ), CONVERT_CONC( 215 ) /  206, 'GC', F /  ! CHBR3
      DATA CGRID_INDEX( 216 ), SPECIES_TYPE( 216 ), CONVERT_CONC( 216 ) /  208, 'GC', F /  ! HCOBR
      DATA CGRID_INDEX( 217 ), SPECIES_TYPE( 217 ), CONVERT_CONC( 217 ) /  207, 'GC', F /  ! HBR
      DATA CGRID_INDEX( 218 ), SPECIES_TYPE( 218 ), CONVERT_CONC( 218 ) /  209, 'GC', F /  ! I2
      DATA CGRID_INDEX( 219 ), SPECIES_TYPE( 219 ), CONVERT_CONC( 219 ) /  210, 'GC', F /  ! I
      DATA CGRID_INDEX( 220 ), SPECIES_TYPE( 220 ), CONVERT_CONC( 220 ) /  211, 'GC', F /  ! IO
      DATA CGRID_INDEX( 221 ), SPECIES_TYPE( 221 ), CONVERT_CONC( 221 ) /  212, 'GC', F /  ! OIO
      DATA CGRID_INDEX( 222 ), SPECIES_TYPE( 222 ), CONVERT_CONC( 222 ) /  213, 'GC', F /  ! I2O2
      DATA CGRID_INDEX( 223 ), SPECIES_TYPE( 223 ), CONVERT_CONC( 223 ) /  214, 'GC', F /  ! HOI
      DATA CGRID_INDEX( 224 ), SPECIES_TYPE( 224 ), CONVERT_CONC( 224 ) /  215, 'GC', F /  ! HI
      DATA CGRID_INDEX( 225 ), SPECIES_TYPE( 225 ), CONVERT_CONC( 225 ) /  216, 'GC', F /  ! INO
      DATA CGRID_INDEX( 226 ), SPECIES_TYPE( 226 ), CONVERT_CONC( 226 ) /  217, 'GC', F /  ! INO2
      DATA CGRID_INDEX( 227 ), SPECIES_TYPE( 227 ), CONVERT_CONC( 227 ) /  218, 'GC', F /  ! INO3
      DATA CGRID_INDEX( 228 ), SPECIES_TYPE( 228 ), CONVERT_CONC( 228 ) /  219, 'GC', F /  ! CH3I
      DATA CGRID_INDEX( 229 ), SPECIES_TYPE( 229 ), CONVERT_CONC( 229 ) /  220, 'GC', F /  ! CH2I2
      DATA CGRID_INDEX( 230 ), SPECIES_TYPE( 230 ), CONVERT_CONC( 230 ) /  231, 'GC', F /  ! I2O3
      DATA CGRID_INDEX( 231 ), SPECIES_TYPE( 231 ), CONVERT_CONC( 231 ) /  232, 'GC', F /  ! I2O4
      DATA CGRID_INDEX( 232 ), SPECIES_TYPE( 232 ), CONVERT_CONC( 232 ) /  221, 'GC', F /  ! BRCL
      DATA CGRID_INDEX( 233 ), SPECIES_TYPE( 233 ), CONVERT_CONC( 233 ) /  222, 'GC', F /  ! ICL
      DATA CGRID_INDEX( 234 ), SPECIES_TYPE( 234 ), CONVERT_CONC( 234 ) /  223, 'GC', F /  ! IBR
      DATA CGRID_INDEX( 235 ), SPECIES_TYPE( 235 ), CONVERT_CONC( 235 ) /  224, 'GC', F /  ! CH2IBR
      DATA CGRID_INDEX( 236 ), SPECIES_TYPE( 236 ), CONVERT_CONC( 236 ) /  225, 'GC', F /  ! CH2ICL
      DATA CGRID_INDEX( 237 ), SPECIES_TYPE( 237 ), CONVERT_CONC( 237 ) /  226, 'GC', F /  ! CHBR2CL
      DATA CGRID_INDEX( 238 ), SPECIES_TYPE( 238 ), CONVERT_CONC( 238 ) /  227, 'GC', F /  ! CHBRCL2
      DATA CGRID_INDEX( 239 ), SPECIES_TYPE( 239 ), CONVERT_CONC( 239 ) /  228, 'GC', F /  ! CH2BRCL
      DATA CGRID_INDEX( 240 ), SPECIES_TYPE( 240 ), CONVERT_CONC( 240 ) /  229, 'GC', F /  ! NO2PIJ
      DATA CGRID_INDEX( 241 ), SPECIES_TYPE( 241 ), CONVERT_CONC( 241 ) /  230, 'GC', F /  ! NO2PK
      DATA CGRID_INDEX( 242 ), SPECIES_TYPE( 242 ), CONVERT_CONC( 242 ) /  277, 'AE', T /  ! ACLI
      DATA CGRID_INDEX( 243 ), SPECIES_TYPE( 243 ), CONVERT_CONC( 243 ) /  278, 'AE', T /  ! ACLJ
      DATA CGRID_INDEX( 244 ), SPECIES_TYPE( 244 ), CONVERT_CONC( 244 ) /  279, 'AE', T /  ! ACLK
      DATA CGRID_INDEX( 245 ), SPECIES_TYPE( 245 ), CONVERT_CONC( 245 ) /  326, 'AE', T /  ! ABRJ
      DATA CGRID_INDEX( 246 ), SPECIES_TYPE( 246 ), CONVERT_CONC( 246 ) /  252, 'AE', T /  ! AFEJ
      DATA CGRID_INDEX( 247 ), SPECIES_TYPE( 247 ), CONVERT_CONC( 247 ) /  259, 'AE', T /  ! AMNJ
      DATA CGRID_INDEX( 248 ), SPECIES_TYPE( 248 ), CONVERT_CONC( 248 ) /  328, 'AE', T /  ! AHMSJ

! The below integers define the locations of mechanism species in the solver
! concentration array.

      INTEGER :: INDEX_O3          =    1
      INTEGER :: INDEX_O3P         =    2
      INTEGER :: INDEX_O1D         =    3
      INTEGER :: INDEX_H2O2        =    4
      INTEGER :: INDEX_HO          =    5
      INTEGER :: INDEX_NO2         =    6
      INTEGER :: INDEX_NO          =    7
      INTEGER :: INDEX_NO3         =    8
      INTEGER :: INDEX_HONO        =    9
      INTEGER :: INDEX_HNO3        =   10
      INTEGER :: INDEX_HNO4        =   11
      INTEGER :: INDEX_HO2         =   12
      INTEGER :: INDEX_HCHO        =   13
      INTEGER :: INDEX_CO          =   14
      INTEGER :: INDEX_ACD         =   15
      INTEGER :: INDEX_MO2         =   16
      INTEGER :: INDEX_ALD         =   17
      INTEGER :: INDEX_ETHP        =   18
      INTEGER :: INDEX_ACT         =   19
      INTEGER :: INDEX_ACO3        =   20
      INTEGER :: INDEX_UALD        =   21
      INTEGER :: INDEX_KET         =   22
      INTEGER :: INDEX_PINAL       =   23
      INTEGER :: INDEX_PINALP      =   24
      INTEGER :: INDEX_LIMAL       =   25
      INTEGER :: INDEX_LIMALP      =   26
      INTEGER :: INDEX_MEK         =   27
      INTEGER :: INDEX_HKET        =   28
      INTEGER :: INDEX_MACR        =   29
      INTEGER :: INDEX_MACP        =   30
      INTEGER :: INDEX_XO2         =   31
      INTEGER :: INDEX_MVK         =   32
      INTEGER :: INDEX_GLY         =   33
      INTEGER :: INDEX_MGLY        =   34
      INTEGER :: INDEX_DCB1        =   35
      INTEGER :: INDEX_DCB2        =   36
      INTEGER :: INDEX_BALD        =   37
      INTEGER :: INDEX_BEN         =   38
      INTEGER :: INDEX_BAL1        =   39
      INTEGER :: INDEX_OP1         =   40
      INTEGER :: INDEX_OP2         =   41
      INTEGER :: INDEX_OPB         =   42
      INTEGER :: INDEX_VOP3        =   43
      INTEGER :: INDEX_PAA         =   44
      INTEGER :: INDEX_ONIT        =   45
      INTEGER :: INDEX_PAN         =   46
      INTEGER :: INDEX_CO2         =   47
      INTEGER :: INDEX_PPN         =   48
      INTEGER :: INDEX_RCO3        =   49
      INTEGER :: INDEX_HC3P        =   50
      INTEGER :: INDEX_VTRPN       =   51
      INTEGER :: INDEX_VHONIT      =   52
      INTEGER :: INDEX_N2O5        =   53
      INTEGER :: INDEX_SO2         =   54
      INTEGER :: INDEX_SULF        =   55
      INTEGER :: INDEX_SULRXN      =   56
      INTEGER :: INDEX_ETH         =   57
      INTEGER :: INDEX_HC3         =   58
      INTEGER :: INDEX_ASOATJ      =   59
      INTEGER :: INDEX_HC5         =   60
      INTEGER :: INDEX_HC5P        =   61
      INTEGER :: INDEX_ETE         =   62
      INTEGER :: INDEX_ETEP        =   63
      INTEGER :: INDEX_OLT         =   64
      INTEGER :: INDEX_OLTP        =   65
      INTEGER :: INDEX_OLI         =   66
      INTEGER :: INDEX_OLIP        =   67
      INTEGER :: INDEX_ACE         =   68
      INTEGER :: INDEX_ORA1        =   69
      INTEGER :: INDEX_BENP        =   70
      INTEGER :: INDEX_PHEN        =   71
      INTEGER :: INDEX_TOL         =   72
      INTEGER :: INDEX_TOLP        =   73
      INTEGER :: INDEX_CSL         =   74
      INTEGER :: INDEX_XYL         =   75
      INTEGER :: INDEX_XYLP        =   76
      INTEGER :: INDEX_EBZ         =   77
      INTEGER :: INDEX_EBZP        =   78
      INTEGER :: INDEX_ISO         =   79
      INTEGER :: INDEX_ISON        =   80
      INTEGER :: INDEX_ISONP       =   81
      INTEGER :: INDEX_ISOP        =   82
      INTEGER :: INDEX_ISHP        =   83
      INTEGER :: INDEX_IEPOX       =   84
      INTEGER :: INDEX_IPX         =   85
      INTEGER :: INDEX_INALD       =   86
      INTEGER :: INDEX_ROH         =   87
      INTEGER :: INDEX_API         =   88
      INTEGER :: INDEX_APIP1       =   89
      INTEGER :: INDEX_APIP2       =   90
      INTEGER :: INDEX_LIM         =   91
      INTEGER :: INDEX_LIMP1       =   92
      INTEGER :: INDEX_LIMP2       =   93
      INTEGER :: INDEX_ACTP        =   94
      INTEGER :: INDEX_MEKP        =   95
      INTEGER :: INDEX_KETP        =   96
      INTEGER :: INDEX_MCP         =   97
      INTEGER :: INDEX_MVKP        =   98
      INTEGER :: INDEX_UALP        =   99
      INTEGER :: INDEX_DCB3        =  100
      INTEGER :: INDEX_BALP        =  101
      INTEGER :: INDEX_ADDC        =  102
      INTEGER :: INDEX_CHO         =  103
      INTEGER :: INDEX_MCT         =  104
      INTEGER :: INDEX_MCTO        =  105
      INTEGER :: INDEX_MOH         =  106
      INTEGER :: INDEX_EOH         =  107
      INTEGER :: INDEX_ETEG        =  108
      INTEGER :: INDEX_HC10P       =  109
      INTEGER :: INDEX_MAHP        =  110
      INTEGER :: INDEX_ORA2        =  111
      INTEGER :: INDEX_ORAP        =  112
      INTEGER :: INDEX_MPAN        =  113
      INTEGER :: INDEX_MCTP        =  114
      INTEGER :: INDEX_OLNN        =  115
      INTEGER :: INDEX_OLND        =  116
      INTEGER :: INDEX_APINP1      =  117
      INTEGER :: INDEX_APINP2      =  118
      INTEGER :: INDEX_LIMNP1      =  119
      INTEGER :: INDEX_LIMNP2      =  120
      INTEGER :: INDEX_ADCN        =  121
      INTEGER :: INDEX_VHOM        =  122
      INTEGER :: INDEX_VROCP4OXY2  =  123
      INTEGER :: INDEX_VROCN1OXY6  =  124
      INTEGER :: INDEX_FURANONE    =  125
      INTEGER :: INDEX_VROCP3OXY2  =  126
      INTEGER :: INDEX_VROCP0OXY4  =  127
      INTEGER :: INDEX_BAL2        =  128
      INTEGER :: INDEX_VELHOM      =  129
      INTEGER :: INDEX_VROCIOXY    =  130
      INTEGER :: INDEX_SLOWROC     =  131
      INTEGER :: INDEX_ACRO        =  132
      INTEGER :: INDEX_BDE13       =  133
      INTEGER :: INDEX_BDE13P      =  134
      INTEGER :: INDEX_FURAN       =  135
      INTEGER :: INDEX_FURANO2     =  136
      INTEGER :: INDEX_PROG        =  137
      INTEGER :: INDEX_SESQ        =  138
      INTEGER :: INDEX_SESQNRO2    =  139
      INTEGER :: INDEX_VROCN2OXY2  =  140
      INTEGER :: INDEX_SESQRO2     =  141
      INTEGER :: INDEX_VROCP0OXY2  =  142
      INTEGER :: INDEX_VROCP1OXY3  =  143
      INTEGER :: INDEX_AGLYOLIGJ   =  144
      INTEGER :: INDEX_IEPOXP      =  145
      INTEGER :: INDEX_AISO3NOSJ   =  146
      INTEGER :: INDEX_ASO4J       =  147
      INTEGER :: INDEX_AISO3OSJ    =  148
      INTEGER :: INDEX_AISO4J      =  149
      INTEGER :: INDEX_AISO5J      =  150
      INTEGER :: INDEX_VROCP6ALK   =  151
      INTEGER :: INDEX_VROCP6ALKP  =  152
      INTEGER :: INDEX_VROCP5ALK   =  153
      INTEGER :: INDEX_VROCP5ALKP  =  154
      INTEGER :: INDEX_VROCP4ALK   =  155
      INTEGER :: INDEX_VROCP4ALKP  =  156
      INTEGER :: INDEX_VROCP3ALK   =  157
      INTEGER :: INDEX_VROCP3ALKP  =  158
      INTEGER :: INDEX_VROCP2ALK   =  159
      INTEGER :: INDEX_VROCP2ALKP  =  160
      INTEGER :: INDEX_VROCP1ALK   =  161
      INTEGER :: INDEX_VROCP1ALKP  =  162
      INTEGER :: INDEX_HC10        =  163
      INTEGER :: INDEX_VROCP6ALKP2 =  164
      INTEGER :: INDEX_VROCP5ALKP2 =  165
      INTEGER :: INDEX_VROCP4ALKP2 =  166
      INTEGER :: INDEX_VROCP2OXY2  =  167
      INTEGER :: INDEX_VROCP3ALKP2 =  168
      INTEGER :: INDEX_VROCP1OXY1  =  169
      INTEGER :: INDEX_VROCP2ALKP2 =  170
      INTEGER :: INDEX_VROCP1ALKP2 =  171
      INTEGER :: INDEX_VROCN1OXY1  =  172
      INTEGER :: INDEX_HC10P2      =  173
      INTEGER :: INDEX_VROCP6ARO   =  174
      INTEGER :: INDEX_VROCP6AROP  =  175
      INTEGER :: INDEX_VROCN2OXY4  =  176
      INTEGER :: INDEX_VROCN1OXY3  =  177
      INTEGER :: INDEX_VROCP5ARO   =  178
      INTEGER :: INDEX_VROCP5AROP  =  179
      INTEGER :: INDEX_NAPH        =  180
      INTEGER :: INDEX_NAPHP       =  181
      INTEGER :: INDEX_VROCN2OXY8  =  182
      INTEGER :: INDEX_VROCP5OXY1  =  183
      INTEGER :: INDEX_VROCP6OXY1  =  184
      INTEGER :: INDEX_ECH4        =  185
      INTEGER :: INDEX_ATRPNJ      =  186
      INTEGER :: INDEX_AHOMJ       =  187
      INTEGER :: INDEX_AHONITJ     =  188
      INTEGER :: INDEX_STY         =  189
      INTEGER :: INDEX_STYP        =  190
      INTEGER :: INDEX_ANO3I       =  191
      INTEGER :: INDEX_ANO3J       =  192
      INTEGER :: INDEX_AMTN1J      =  193
      INTEGER :: INDEX_CL2         =  194
      INTEGER :: INDEX_CL          =  195
      INTEGER :: INDEX_CLO         =  196
      INTEGER :: INDEX_OCLO        =  197
      INTEGER :: INDEX_CL2O2       =  198
      INTEGER :: INDEX_CLOO        =  199
      INTEGER :: INDEX_HOCL        =  200
      INTEGER :: INDEX_CLNO        =  201
      INTEGER :: INDEX_CLNO2       =  202
      INTEGER :: INDEX_CLNO3       =  203
      INTEGER :: INDEX_HCOCL       =  204
      INTEGER :: INDEX_HCL         =  205
      INTEGER :: INDEX_BR2         =  206
      INTEGER :: INDEX_BR          =  207
      INTEGER :: INDEX_BRO         =  208
      INTEGER :: INDEX_OBRO        =  209
      INTEGER :: INDEX_HOBR        =  210
      INTEGER :: INDEX_BRNO        =  211
      INTEGER :: INDEX_BRNO2       =  212
      INTEGER :: INDEX_BRNO3       =  213
      INTEGER :: INDEX_CH2BR2      =  214
      INTEGER :: INDEX_CHBR3       =  215
      INTEGER :: INDEX_HCOBR       =  216
      INTEGER :: INDEX_HBR         =  217
      INTEGER :: INDEX_I2          =  218
      INTEGER :: INDEX_I           =  219
      INTEGER :: INDEX_IO          =  220
      INTEGER :: INDEX_OIO         =  221
      INTEGER :: INDEX_I2O2        =  222
      INTEGER :: INDEX_HOI         =  223
      INTEGER :: INDEX_HI          =  224
      INTEGER :: INDEX_INO         =  225
      INTEGER :: INDEX_INO2        =  226
      INTEGER :: INDEX_INO3        =  227
      INTEGER :: INDEX_CH3I        =  228
      INTEGER :: INDEX_CH2I2       =  229
      INTEGER :: INDEX_I2O3        =  230
      INTEGER :: INDEX_I2O4        =  231
      INTEGER :: INDEX_BRCL        =  232
      INTEGER :: INDEX_ICL         =  233
      INTEGER :: INDEX_IBR         =  234
      INTEGER :: INDEX_CH2IBR      =  235
      INTEGER :: INDEX_CH2ICL      =  236
      INTEGER :: INDEX_CHBR2CL     =  237
      INTEGER :: INDEX_CHBRCL2     =  238
      INTEGER :: INDEX_CH2BRCL     =  239
      INTEGER :: INDEX_NO2PIJ      =  240
      INTEGER :: INDEX_NO2PK       =  241
      INTEGER :: INDEX_ACLI        =  242
      INTEGER :: INDEX_ACLJ        =  243
      INTEGER :: INDEX_ACLK        =  244
      INTEGER :: INDEX_ABRJ        =  245
      INTEGER :: INDEX_AFEJ        =  246
      INTEGER :: INDEX_AMNJ        =  247
      INTEGER :: INDEX_AHMSJ       =  248

      INTEGER, PARAMETER :: N_ACT_SP = 248

      INTEGER, PARAMETER :: NRXNS = 921

      INTEGER, PARAMETER ::     ONE_REACT_REACTIONS =    0

      INTEGER, PARAMETER ::     TWO_REACT_REACTIONS =    0

      INTEGER, PARAMETER ::   THREE_REACT_REACTIONS =    0

      INTEGER, PARAMETER ::    ZERO_REACT_REACTIONS =    0

      LOGICAL, PARAMETER ::       UNITARY_REACTIONS = .FALSE.

      INTEGER, PARAMETER ::         ONE_REACT_START =    0

      INTEGER, PARAMETER ::         ONE_REACT_STOP  =   -1

      LOGICAL, PARAMETER ::       BINARY_REACTIONS  = .FALSE.

      INTEGER, PARAMETER ::         TWO_REACT_START =    0

      INTEGER, PARAMETER ::         TWO_REACT_STOP  =   -1

      LOGICAL, PARAMETER ::       TERNARY_REACTIONS = .FALSE.

      INTEGER, PARAMETER ::       THREE_REACT_START =    0

      INTEGER, PARAMETER ::       THREE_REACT_STOP  =   -1

      LOGICAL, PARAMETER ::       NULL_REACTIONS    = .FALSE.

      INTEGER, PARAMETER ::        ZERO_REACT_START =    0

      INTEGER, PARAMETER ::        ZERO_REACT_STOP  =   -1

      INTEGER, PARAMETER ::        NSUNLIGHT_RXNS   =   86

      INTEGER, PARAMETER ::        NTHERMAL_RXNS    =  835

      INTEGER, PARAMETER ::        KUNITS           =    2

      INTEGER  :: IRXXN

      INTEGER, PARAMETER :: NMPHOT =  86
      INTEGER            :: IPH( NMPHOT,3 )

      DATA ( IPH( IRXXN,1 ), IRXXN = 1, NMPHOT ) / & 
     &      1,    2,    3,    4,    5,    6,    7,    8,    9,   10, & ! O   
     &     11,   12,   13,   14,   15,   16,   17,   18,   19,   20, & ! 1   
     &     21,   22,   23,   24,   25,   26,   27,   28,   29,   30, & ! 2   
     &     31,   32,   33,   34,   35,   36,   37,   38,   39,   40, & ! 3   
     &     41,   42,   43,  404,  542,  543,  544,  545,  546,  547, & ! 4   
     &    548,  549,  550,  551,  552,  740,  741,  742,  743,  744, & ! 5   
     &    745,  746,  747,  748,  749,  750,  811,  812,  813,  814, & ! 6   
     &    815,  816,  817,  818,  819,  820,  821,  822,  823,  851, & ! 7   
     &    852,  853,  854,  855,  856,  857/     !  8   

      DATA ( IPH( IRXXN,2 ), IRXXN = 1, NMPHOT ) / & 
     &      1,    2,    3,    4,    5,    6,    7,    8,    9,   10, & ! O   
     &     11,   12,   13,   14,   15,   16,   13,   13,   17,   18, & ! 1   
     &     19,   20,   21,   22,   23,   24,   25,   25,   25,   26, & ! 2   
     &     27,   28,   28,   28,   28,   29,   30,   31,   32,   33, & ! 3   
     &     34,   35,   35,   36,    4,   37,   38,   39,   40,   41, & ! 4   
     &     42,   43,   44,   45,   46,   47,   48,   49,   50,   51, & ! 5   
     &     52,   53,   54,   55,   56,   57,   58,   59,   60,   61, & ! 6   
     &     62,   63,   64,   65,   61,   66,   67,   61,   61,   68, & ! 7   
     &     69,   70,   71,   72,   73,   74/     !  8   

      DATA ( IPH( IRXXN,3 ), IRXXN = 1, NMPHOT ) / & 
     &      1,    2,    3,    4,    5,    6,    7,    8,    9,   10, & ! O   
     &     11,   12,   13,   14,   15,   16,   17,   18,   19,   20, & ! 1   
     &     21,   22,   23,   24,   25,   26,   27,   28,   29,   30, & ! 2   
     &     31,   32,   33,   34,   35,   36,   37,   38,   39,   40, & ! 3   
     &     41,   42,   43,   44,   45,   46,   47,   48,   49,   50, & ! 4   
     &     51,   52,   53,   54,   55,   56,   57,   58,   59,   60, & ! 5   
     &     61,   62,   63,   64,   65,   66,   67,   68,   69,   70, & ! 6   
     &     71,   72,   73,   74,   75,   76,   77,   78,   79,   80, & ! 7   
     &     81,   82,   83,   84,   85,   86/     !  8   

      REAL( 8 )          :: RTDAT( 3,NRXNS )

      INTEGER, PARAMETER :: NFALLOFF =  26
      REAL( 8 )          :: RFDAT( 5,NFALLOFF )

      INTEGER            :: KTYPE( NRXNS )

      DATA ( KTYPE( IRXXN ), IRXXN = 1, NRXNS ) /  & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    3,    3,    3,    3,    2,    3,    3, & ! 4   
     &      3,    3,    3,    3,    9,    9,    3,   10,   10,    3, & ! 5   
     &      9,    3,    3,    3,   10,   10,    8,    1,    1,    3, & ! 6   
     &      3,    3,   10,    5,    1,   10,    5,    3,   10,    9, & ! 7   
     &      3,    3,    3,    3,   10,    3,    3,   10,    3,    3, & ! 8   
     &      1,    1,    1,    3,    3,    3,    3,    3,    3,    3, & ! 9   
     &      3,    3,    3,    1,    1,    1,    3,    3,    3,    1, & ! O   
     &      3,    3,    3,    4,    3,    3,    1,    3,    3,    3, & ! 1   
     &      1,    3,    3,    3,    1,    3,    3,    3,    3,    3, & ! 2   
     &      3,    3,    1,    3,    3,    3,    1,    1,    3,    3, & ! 3   
     &      1,    1,    1,    3,    1,    8,    3,    3,    3,    3, & ! 4   
     &      3,    1,    1,    3,    3,    1,    1,    1,    1,    1, & ! 5   
     &      4,    3,    3,    3,    1,    3,    3,    3,    3,    1, & ! 6   
     &      3,    3,    3,    1,    1,    1,    1,    3,   10,    1, & ! 7   
     &      1,   10,    5,   10,    5,    3,    3,    3,    3,    3, & ! 8   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 9   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! O   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 1   
     &      3,    3,    1,    3,    3,    3,    1,    1,    1,    3, & ! 2   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 3   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 4   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 5   
     &      3,    1,    3,    1,    3,    3,    3,    3,    3,    3, & ! 6   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 7   
     &      3,    3,    1,    1,    1,    1,    1,    1,    1,    1, & ! 8   
     &      3,    3,    3,    3,    3,    3,    3,    1,    3,    3, & ! 9   
     &      3,    3,    3,    3,    3,    3,    1,    3,    3,    3, & ! O   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 1   
     &      1,    3,    1,    3,    1,    3,    1,    3,    3,    3, & ! 2   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 3   
     &      3,    3,    3,    3,    3,    3,    1,    1,    1,    1, & ! 4   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 5   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 6   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    3, & ! 7   
     &      3,    3,    1,    3,    3,    1,    1,    1,    1,    1, & ! 8   
     &      1,    1,    1,    1,    1,    1,    1,    3,    1,    1, & ! 9   
     &      3,    1,    1,    0,    3,    3,    1,    3,    1,    1, & ! O   
     &      3,    1,    1,    3,    3,    1,    1,    1,    1,    1, & ! 1   
     &      3,    3,    1,    1,    1,    3,    1,    3,   -1,   -1, & ! 2   
     &     -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,    1,    1, & ! 3   
     &      1,    1,    1,    1,    1,    3,    3,    3,    3,    3, & ! 4   
     &      3,    3,    1,    1,    1,    1,    1,    1,    1,    3, & ! 5   
     &      3,    3,    3,    3,    3,    3,    1,    1,    1,    1, & ! 6   
     &      1,    1,    1,    3,    3,    3,    3,    3,    3,    3, & ! 7   
     &      1,    1,    1,    1,    1,    1,    1,    3,    3,    3, & ! 8   
     &      3,    3,    3,    3,    1,    3,    3,    1,    3,    3, & ! 9   
     &      1,    3,    3,    1,    3,    3,    1,    3,    3,    1, & ! O   
     &      3,    3,    1,    1,    1,    1,    1,    1,    1,    1, & ! 1   
     &      1,    1,    1,    1,    1,    1,    1,    1,    3,    1, & ! 2   
     &      1,    3,    1,    1,    3,    3,    1,    1,    1,   -1, & ! 3   
     &     -1,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    3,    3,    3,    3,    3,    3,    3,    3, & ! 5   
     &      3,    3,    2,   10,    1,    3,    3,    3,    3,   10, & ! 6   
     &      3,    3,   10,    5,    3,    3,    3,    3,    3,    3, & ! 7   
     &      3,    3,    3,    3,    3,    3,    3,    3,    1,    3, & ! 8   
     &      3,    3,    3,    3,    3,    3,    1,    1,    1,    1, & ! 9   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! O   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 1   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 2   
     &      1,    3,    1,    1,    1,    1,    1,    3,    1,    1, & ! 3   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 4   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 5   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 6   
     &      1,    1,    1,    3,    3,    1,    1,    3,    3,    3, & ! 7   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 8   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 9   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! O   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 1   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 2   
     &      3,    3,    3,    3,    3,    3,    3,    3,    3,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      3,    3,    3,   10,    3,    3,    3,    3,    3,    3, & ! 5   
     &     10,    1,    3,    3,    1,    3,    3,    3,    1,    1, & ! 6   
     &      1,    3,    3,    1,    3,    3,    1,    1,    1,    1, & ! 7   
     &      1,    1,    1,    1,    1,    1,    3,    1,    1,    1, & ! 8   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 9   
     &      1,    3,    3,    1,    1,    3,    1,    3,    1,    1, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    1,    1,    3,    3,   10,   10,    3, & ! 2   
     &      3,    3,    3,    3,    1,    1,    3,    3,    3,   10, & ! 3   
     &      3,    3,    2,    3,    3,    3,    3,    3,    1,    3, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    3,    3,    1, & ! 5   
     &      1,    3,    3,    3,    3,    3,    3,    3,    3,    3, & ! 6   
     &      3,    3,    3,   -1,   -1,   -1,   -1,   -1,   -1,   -1, & ! 7   
     &     -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, & ! 8   
     &     -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, & ! 9   
     &     -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, & ! O   
     &     -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1, & ! 1   
     &     -1/     !  2   

      INTEGER            :: IRXBITS( NRXNS )

      DATA ( IRXBITS( IRXXN ), IRXXN = 1, NRXNS ) / & 
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    0,    0,    0,    0,   20,    0,   16, & ! 4   
     &     32,    8,  128,    0,    0,    8,    0,    1,    1,    0, & ! 5   
     &      0,   16,    0,    0,    1,    1,    0,    0,    0,    0, & ! 6   
     &      0,    0,    1,    0,    8,    1,    0,    0,    1,    0, & ! 7   
     &     64,    0,    0,    0,    1,    0,    0,    1,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    1,    0, & ! 7   
     &      0,    1,    0,    1,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    2,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    1,    1, & ! 2   
     &      1,    1,    1,    1,    1,    1,    1,    1,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    1, & ! 3   
     &      1,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    0,  128,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    4,    1,    0,    0,    0,    0,    0,    1, & ! 6   
     &      0,    0,    1,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,   64,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      0,    0,    0,    1,    0,    0,    0,    0,    0,    0, & ! 5   
     &      1,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    0,    0,    0,    0,    1,    1,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    1, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    1,    1,    1,    1,    1,    1,    1, & ! 7   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 8   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 9   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! O   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 1   
     &      1/     !  2   

      INTEGER, PARAMETER :: NTERMS_JACOB =    61504

      INTEGER, PARAMETER :: NSTEPS_JACOB =     1842

      INTEGER            :: IORDER( NRXNS )

      DATA ( IORDER( IRXXN ), IRXXN = 1, NRXNS ) / & 
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! O   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 1   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 2   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 3   
     &      1,    1,    1,    2,    2,    2,    2,    3,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    3,    2,    2,    2,    2, & ! 5   
     &      2,    3,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    1,    2,    2,    1,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    1,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 7   
     &      1,    2,    1,    2,    1,    2,    1,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    1,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    1,    1, & ! 2   
     &      1,    1,    1,    1,    1,    2,    1,    1,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    1,    1,    1,    1, & ! 6   
     &      1,    1,    1,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 2   
     &      1,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 3   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 4   
     &      1,    1,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    3,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    1,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 3   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 1   
     &      1,    1,    1,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      1,    1,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    1,    2,    1,    2,    2,    1,    2,    2,    2, & ! 4   
     &      1,    1,    1,    1,    1,    1,    1,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    1,    1,    1,    1,    2,    2,    2, & ! 7   
     &      2,    2,    1,    1,    2,    2,    1,    1,    2,    2, & ! 8   
     &      2,    2,    2,    2,    1,    1,    2,    2,    1,    2, & ! 9   
     &      2,    2,    2,    2,    2,    1,    1,    1,    1,    1, & ! O   
     &      1,    2,    2,    2,    2,    2,    3,    2,    2,    2, & ! 1   
     &      1/     !  2   

      INTEGER, PARAMETER :: KTN1 = 297
      INTEGER            :: KRX1( KTN1 )

      DATA ( KRX1( IRXXN ), IRXXN = 1, KTN1 ) / & 
     &     68,   69,   75,   91,   92,   93,  104,  105,  106,  110, & ! O   
     &    117,  121,  125,  133,  137,  138,  141,  142,  143,  145, & ! 1   
     &    152,  153,  156,  157,  158,  159,  160,  165,  170,  174, & ! 2   
     &    175,  176,  177,  180,  181,  223,  227,  228,  229,  262, & ! 3   
     &    264,  283,  284,  285,  286,  287,  288,  289,  290,  298, & ! 4   
     &    307,  321,  323,  325,  327,  347,  348,  349,  350,  351, & ! 5   
     &    352,  353,  354,  355,  356,  357,  358,  359,  360,  361, & ! 6   
     &    362,  363,  364,  365,  366,  367,  368,  369,  370,  371, & ! 7   
     &    372,  373,  374,  375,  376,  377,  378,  379,  383,  386, & ! 8   
     &    387,  388,  389,  390,  391,  392,  393,  394,  395,  396, & ! 9   
     &    397,  399,  400,  402,  403,  407,  409,  410,  412,  413, & ! O   
     &    416,  417,  418,  419,  420,  423,  424,  425,  427,  439, & ! 1   
     &    440,  441,  442,  443,  444,  445,  453,  454,  455,  456, & ! 2   
     &    457,  458,  459,  467,  468,  469,  470,  471,  472,  473, & ! 3   
     &    481,  482,  483,  484,  485,  486,  487,  495,  498,  501, & ! 4   
     &    504,  507,  510,  513,  514,  515,  516,  517,  518,  519, & ! 5   
     &    520,  521,  522,  523,  524,  525,  526,  527,  528,  530, & ! 6   
     &    531,  533,  534,  537,  538,  539,  565,  589,  597,  598, & ! 7   
     &    599,  600,  601,  602,  603,  604,  605,  606,  607,  608, & ! 8   
     &    609,  610,  611,  612,  613,  614,  615,  616,  617,  618, & ! 9   
     &    619,  620,  621,  622,  623,  624,  625,  626,  627,  628, & ! O   
     &    629,  630,  631,  633,  634,  635,  636,  637,  639,  640, & ! 1   
     &    641,  642,  643,  644,  645,  646,  647,  648,  649,  650, & ! 2   
     &    651,  652,  653,  654,  655,  656,  657,  658,  659,  660, & ! 3   
     &    661,  662,  663,  664,  665,  666,  667,  668,  669,  670, & ! 4   
     &    671,  672,  673,  676,  677,  762,  765,  769,  770,  771, & ! 5   
     &    774,  777,  778,  779,  780,  781,  782,  783,  784,  785, & ! 6   
     &    786,  788,  789,  790,  791,  792,  793,  794,  795,  796, & ! 7   
     &    797,  798,  799,  800,  801,  804,  805,  807,  809,  810, & ! 8   
     &    824,  825,  835,  836,  849,  860,  861/     !  9   

      INTEGER, PARAMETER :: KTN2 =   3
      INTEGER            :: KRX2( KTN2 )

      DATA ( KRX2( IRXXN ), IRXXN = 1, KTN2 ) / & 
     &     48,  563,  843/

      INTEGER, PARAMETER :: KTN3 = 442
      INTEGER            :: KRX3( KTN3 )

      DATA ( KRX3( IRXXN ), IRXXN = 1, KTN3 ) / & 
     &     44,   45,   46,   47,   49,   50,   51,   52,   53,   54, & ! O   
     &     57,   60,   62,   63,   64,   70,   71,   72,   78,   81, & ! 1   
     &     82,   83,   84,   86,   87,   89,   90,   94,   95,   96, & ! 2   
     &     97,   98,   99,  100,  101,  102,  103,  107,  108,  109, & ! 3   
     &    111,  112,  113,  115,  116,  118,  119,  120,  122,  123, & ! 4   
     &    124,  126,  127,  128,  129,  130,  131,  132,  134,  135, & ! 5   
     &    136,  139,  140,  144,  147,  148,  149,  150,  151,  154, & ! 6   
     &    155,  162,  163,  164,  166,  167,  168,  169,  171,  172, & ! 7   
     &    173,  178,  186,  187,  188,  189,  190,  191,  192,  193, & ! 8   
     &    194,  195,  196,  197,  198,  199,  200,  201,  202,  203, & ! 9   
     &    204,  205,  206,  207,  208,  209,  210,  211,  212,  213, & ! O   
     &    214,  215,  216,  217,  218,  219,  220,  221,  222,  224, & ! 1   
     &    225,  226,  230,  231,  232,  233,  234,  235,  236,  237, & ! 2   
     &    238,  239,  240,  241,  242,  243,  244,  245,  246,  247, & ! 3   
     &    248,  249,  250,  251,  252,  253,  254,  255,  256,  257, & ! 4   
     &    258,  259,  260,  261,  263,  265,  266,  267,  268,  269, & ! 5   
     &    270,  271,  272,  273,  274,  275,  276,  277,  278,  279, & ! 6   
     &    280,  281,  282,  291,  292,  293,  294,  295,  296,  297, & ! 7   
     &    299,  300,  301,  302,  303,  304,  305,  306,  308,  309, & ! 8   
     &    310,  311,  312,  313,  314,  315,  316,  317,  318,  319, & ! 9   
     &    320,  322,  324,  326,  328,  329,  330,  331,  332,  333, & ! O   
     &    334,  335,  336,  337,  338,  339,  340,  341,  342,  343, & ! 1   
     &    344,  345,  346,  380,  381,  382,  384,  385,  398,  401, & ! 2   
     &    405,  406,  408,  411,  414,  415,  421,  422,  426,  428, & ! 3   
     &    446,  447,  448,  449,  450,  451,  452,  460,  461,  462, & ! 4   
     &    463,  464,  465,  466,  474,  475,  476,  477,  478,  479, & ! 5   
     &    480,  488,  489,  490,  491,  492,  493,  494,  496,  497, & ! 6   
     &    499,  500,  502,  503,  505,  506,  508,  509,  511,  512, & ! 7   
     &    529,  532,  535,  536,  553,  554,  555,  556,  557,  558, & ! 8   
     &    559,  560,  561,  562,  566,  567,  568,  569,  571,  572, & ! 9   
     &    575,  576,  577,  578,  579,  580,  581,  582,  583,  584, & ! O   
     &    585,  586,  587,  588,  590,  591,  592,  593,  594,  595, & ! 1   
     &    596,  632,  638,  674,  675,  678,  679,  680,  681,  682, & ! 2   
     &    683,  684,  685,  686,  687,  688,  689,  690,  691,  692, & ! 3   
     &    693,  694,  695,  696,  697,  698,  699,  700,  701,  702, & ! 4   
     &    703,  704,  705,  706,  707,  708,  709,  710,  711,  712, & ! 5   
     &    713,  714,  715,  716,  717,  718,  719,  720,  721,  722, & ! 6   
     &    723,  724,  725,  726,  727,  728,  729,  730,  731,  732, & ! 7   
     &    733,  734,  735,  736,  737,  738,  739,  751,  752,  753, & ! 8   
     &    755,  756,  757,  758,  759,  760,  763,  764,  766,  767, & ! 9   
     &    768,  772,  773,  775,  776,  787,  802,  803,  806,  808, & ! O   
     &    826,  827,  830,  831,  832,  833,  834,  837,  838,  839, & ! 1   
     &    841,  842,  844,  845,  846,  847,  848,  850,  858,  859, & ! 2   
     &    862,  863,  864,  865,  866,  867,  868,  869,  870,  871, & ! 3   
     &    872,  873/     !  4   

      INTEGER, PARAMETER :: KTN4 =   2
      INTEGER            :: KRX4( KTN4 )

      DATA ( KRX4( IRXXN ), IRXXN = 1, KTN4 ) / & 
     &    114,  161/

      INTEGER, PARAMETER :: KTN5 =   5
      INTEGER            :: KRX5( KTN5 )

      DATA ( KRX5( IRXXN ), IRXXN = 1, KTN5 ) / & 
     &     74,   77,  183,  185,  574/

      INTEGER, PARAMETER :: KTN6 =   0
      INTEGER            :: KRX6( 1 )

      DATA   KRX6( 1 ) / 0 /

      INTEGER, PARAMETER :: KTN7 =   0
      INTEGER            :: KRX7( 1 )

      DATA   KRX7( 1 ) / 0 /

      INTEGER, PARAMETER :: NWM =   2
      INTEGER            :: NRXWM( NWM )

      DATA ( NRXWM( IRXXN ), IRXXN = 1, NWM ) /  & 
     &     48,  563/
      REAL( 8 ),    PARAMETER :: ATM_AIR = 1.00000D+06

      INTEGER, PARAMETER :: NWW =   3
      INTEGER            :: NRXWW( NWW )

      DATA ( NRXWW( IRXXN ), IRXXN = 1, NWW ) / & 
     &     52,   56,   75/

      INTEGER, PARAMETER :: NWO2 =   3
      INTEGER            :: NRXWO2( NWO2 )

      DATA ( NRXWO2( IRXXN ), IRXXN = 1, NWO2 ) / & 
     &     48,   50,   62/
      REAL( 8 ),    PARAMETER :: ATM_O2 = 2.09500D+05

      INTEGER, PARAMETER :: NWN2 =   1
      INTEGER            :: NRXWN2( NWN2 )

      DATA ( NRXWN2( IRXXN ), IRXXN = 1, NWN2 ) / & 
     &     51/
      REAL( 8 ),    PARAMETER :: ATM_N2 = 7.80800D+05

      INTEGER, PARAMETER :: NWCH4 =   2
      INTEGER            :: NRXWCH4( NWCH4 )

      DATA ( NRXWCH4( IRXXN ), IRXXN = 1, NWCH4 ) / & 
     &     81,  585/
      REAL( 8 ),    PARAMETER :: ATM_CH4 = 1.85000D+00

      INTEGER, PARAMETER :: NWH2 =   2
      INTEGER            :: NRXWH2( NWH2 )

      DATA ( NRXWH2( IRXXN ), IRXXN = 1, NWH2 ) / & 
     &     53,  554/
      REAL( 8 ),    PARAMETER :: ATM_H2 = 5.60000D-01

      INTEGER, PARAMETER :: N_ATOMS =  14
      CHARACTER( 16 )    :: ATOMS( N_ATOMS )

      DATA ( ATOMS( IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &   'CA', 'MN', 'CL', 'HG', 'BR', 'NA', & 
     &   'SI', 'S ', 'TI', 'FE', 'K ', 'I ', & 
     &   'N ', 'C '/
      REAL( 8 )            :: CHEM_SPC_ATOMS( NUMB_MECH_SPC,N_ATOMS )
! O3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   1,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! O3P found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   2,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! O1D found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   3,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! H2O2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   4,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! HO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   5,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! NO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   6,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! NO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   7,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! NO3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   8,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! HONO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(   9,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! HNO3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  10,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! HNO4 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  11,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! HO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  12,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! HCHO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  13,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! CO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  14,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! ACD found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  15,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! MO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  16,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! ALD found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  17,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! ETHP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  18,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! ACT found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  19,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! ACO3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  20,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! UALD found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  21,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! KET found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  22,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! PINAL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  23,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! PINALP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  24,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! LIMAL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  25,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! LIMALP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  26,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! MEK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  27,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! HKET found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  28,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! MACR found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  29,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! MACP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  30,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! XO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  31,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! MVK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  32,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! GLY found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  33,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! MGLY found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  34,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! DCB1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  35,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! DCB2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  36,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! BALD found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  37,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! BEN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  38,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! BAL1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  39,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! OP1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  40,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! OP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  41,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! OPB found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  42,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! VOP3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  43,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! PAA found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  44,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! ONIT found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  45,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 4.0000D+00/
! PAN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  46,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 2.0000D+00/
! CO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  47,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! PPN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  48,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 3.0000D+00/
! RCO3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  49,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! HC3P found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  50,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! VTRPN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  51,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+01/
! VHONIT found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  52,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+01/
! N2O5 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  53,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 2.0000D+00, 0.0000D+00/
! SO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  54,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! SULF found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  55,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! SULRXN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  56,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! ETH found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  57,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! HC3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  58,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! ASOATJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS(  59,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! HC5 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  60,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! HC5P found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  61,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! ETE found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  62,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! ETEP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  63,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! OLT found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  64,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! OLTP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  65,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! OLI found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  66,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! OLIP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  67,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! ACE found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  68,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! ORA1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  69,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! BENP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  70,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! PHEN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  71,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! TOL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  72,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! TOLP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  73,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! CSL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  74,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! XYL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  75,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! XYLP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  76,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! EBZ found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  77,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! EBZP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  78,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! ISO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  79,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! ISON found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  80,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 5.0000D+00/
! ISONP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  81,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 5.0000D+00/
! ISOP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  82,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! ISHP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  83,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! IEPOX found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  84,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! IPX found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  85,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! INALD found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  86,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 5.0000D+00/
! ROH found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  87,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! API found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  88,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! APIP1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  89,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! APIP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  90,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! LIM found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  91,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! LIMP1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  92,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! LIMP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  93,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! ACTP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  94,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! MEKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  95,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! KETP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  96,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! MCP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  97,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! MVKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  98,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! UALP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS(  99,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! DCB3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 100,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! BALP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 101,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! ADDC found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 102,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! CHO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 103,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! MCT found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 104,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! MCTO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 105,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! MOH found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 106,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! EOH found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 107,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! ETEG found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 108,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! HC10P found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 109,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! MAHP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 110,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! ORA2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 111,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! ORAP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 112,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00/
! MPAN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 113,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 4.0000D+00/
! MCTP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 114,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! OLNN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 115,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 3.0000D+00/
! OLND found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 116,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 3.0000D+00/
! APINP1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 117,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+01/
! APINP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 118,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+01/
! LIMNP1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 119,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+01/
! LIMNP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 120,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+01/
! ADCN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 121,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 6.0000D+00/
! VHOM found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 122,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! VROCP4OXY2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 123,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 9.0000D+00/
! VROCN1OXY6 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 124,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! FURANONE found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 125,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! VROCP3OXY2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 126,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.1000D+01/
! VROCP0OXY4 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 127,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! BAL2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 128,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! VELHOM found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 129,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+01/
! VROCIOXY found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 130,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 5.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! SLOWROC found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 131,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00/
! ACRO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 132,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! BDE13 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 133,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! BDE13P found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 134,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.0000D+00/
! FURAN found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 135,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! FURANO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 136,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! PROG found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 137,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00/
! SESQ found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 138,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.5000D+01/
! SESQNRO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 139,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.5000D+01/
! VROCN2OXY2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 140,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.7000D+01/
! SESQRO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 141,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.5000D+01/
! VROCP0OXY2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 142,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D+01/
! VROCP1OXY3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 143,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.1000D+01/
! AGLYOLIGJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 144,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.0000D+00/
! IEPOXP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 145,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! AISO3NOSJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 146,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! ASO4J found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 147,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! AISO3OSJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 148,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! AISO4J found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 149,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! AISO5J found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 150,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D+00/
! VROCP6ALK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 151,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3000D+01/
! VROCP6ALKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 152,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3000D+01/
! VROCP5ALK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 153,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D+01/
! VROCP5ALKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 154,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D+01/
! VROCP4ALK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 155,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.8000D+01/
! VROCP4ALKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 156,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.8000D+01/
! VROCP3ALK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 157,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.1000D+01/
! VROCP3ALKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 158,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.1000D+01/
! VROCP2ALK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 159,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4000D+01/
! VROCP2ALKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 160,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4000D+01/
! VROCP1ALK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 161,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.7000D+01/
! VROCP1ALKP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 162,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.7000D+01/
! HC10 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 163,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! VROCP6ALKP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 164,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3000D+01/
! VROCP5ALKP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 165,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D+01/
! VROCP4ALKP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 166,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.8000D+01/
! VROCP2OXY2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 167,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.2000D+01/
! VROCP3ALKP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 168,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.1000D+01/
! VROCP1OXY1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 169,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.7000D+01/
! VROCP2ALKP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 170,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4000D+01/
! VROCP1ALKP2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 171,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.7000D+01/
! VROCN1OXY1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 172,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+01/
! HC10P2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 173,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! VROCP6ARO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 174,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3000D+01/
! VROCP6AROP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 175,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3000D+01/
! VROCN2OXY4 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 176,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.1000D+01/
! VROCN1OXY3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 177,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.2000D+01/
! VROCP5ARO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 178,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D+01/
! VROCP5AROP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 179,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D+01/
! NAPH found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 180,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! NAPHP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 181,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+01/
! VROCN2OXY8 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 182,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D+00/
! VROCP5OXY1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 183,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.1000D+01/
! VROCP6OXY1 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 184,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 9.0000D+00/
! ECH4 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 185,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! ATRPNJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 186,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! AHOMJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 187,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! AHONITJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 188,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! STY found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 189,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! STYP found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 190,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D+00/
! ANO3I found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 191,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! ANO3J found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 192,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! AMTN1J found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 193,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! CL2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 194,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 2.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! CL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 195,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! CLO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 196,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! OCLO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 197,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! CL2O2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 198,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 2.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! CLOO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 199,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! HOCL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 200,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! CLNO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 201,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! CLNO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 202,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! CLNO3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 203,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! HCOCL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 204,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! HCL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 205,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! BR2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 206,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! BR found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 207,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! BRO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 208,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! OBRO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 209,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! HOBR found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 210,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! BRNO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 211,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! BRNO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 212,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! BRNO3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 213,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! CH2BR2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 214,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! CHBR3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 215,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! HCOBR found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 216,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! HBR found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 217,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! I2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 218,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 2.0000D+00, 0.0000D+00, 0.0000D+00/
! I found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 219,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00/
! IO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 220,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00/
! OIO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 221,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00/
! I2O2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 222,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 2.0000D+00, 0.0000D+00, 0.0000D+00/
! HOI found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 223,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00/
! HI found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 224,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00/
! INO found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 225,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00/
! INO2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 226,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00/
! INO3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 227,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00/
! CH3I found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 228,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00/
! CH2I2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 229,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 2.0000D+00, 0.0000D+00, 1.0000D+00/
! I2O3 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 230,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 2.0000D+00, 0.0000D+00, 0.0000D+00/
! I2O4 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 231,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 2.0000D+00, 0.0000D+00, 0.0000D+00/
! BRCL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 232,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! ICL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 233,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00/
! IBR found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 234,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00/
! CH2IBR found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 235,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00/
! CH2ICL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 236,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00/
! CHBR2CL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 237,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 2.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! CHBRCL2 found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 238,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 2.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! CH2BRCL found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 239,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/
! NO2PIJ found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 240,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! NO2PK found in GC namelist 
      DATA ( CHEM_SPC_ATOMS( 241,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00/
! ACLI found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 242,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! ACLJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 243,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! ACLK found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 244,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! ABRJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 245,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! AFEJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 246,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! AMNJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 247,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00/
! AHMSJ found in AE namelist 
      DATA ( CHEM_SPC_ATOMS( 248,IRXXN ), IRXXN = 1, N_ATOMS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00/

      INTEGER, PARAMETER :: MXPRD =  27
      INTEGER            :: IRR( NRXNS,MXPRD+3 )

      DATA ( IRR( IRXXN,  1 ), IRXXN = 1, NRXNS ) / & 
     &      1,    1,    4,    6,    8,    8,    9,   10,   11,   13, & ! O   
     &     13,   15,   17,   19,   19,   21,   23,   25,   27,   22, & ! 1   
     &     28,   29,   32,   33,   33,   33,   34,   35,   36,   37, & ! 2   
     &     37,   40,   41,   42,   43,   44,   45,   46,   46,   48, & ! 3   
     &     48,   51,   52,    1,    1,    1,    1,    2,    2,    3, & ! 4   
     &      3,    3,    5,    5,   12,   12,    4,    7,    7,    7, & ! 5   
     &      7,    7,    9,    6,    6,    6,   10,    8,    8,    8, & ! 6   
     &      8,    8,    8,   53,   53,    6,   11,   11,   54,   14, & ! 7   
     &      5,   57,   58,   60,   62,   64,   66,   68,   38,   72, & ! 8   
     &     75,   77,   77,   79,   79,   79,   82,   82,   83,   83, & ! 9   
     &     81,   81,   80,    5,   80,   85,   88,   91,   23,   25, & ! O   
     &     13,   15,   17,   19,   27,   22,   28,   29,   32,   21, & ! 1   
     &     33,   34,   35,   36,  100,   37,   71,   74,  104,  106, & ! 2   
     &    107,   87,  108,   40,   41,   42,  110,   69,  111,   44, & ! 3   
     &     46,   48,  113,   45,   51,   52,   62,   64,   66,   88, & ! 4   
     &     91,   25,   51,   29,   32,   21,   35,   36,  100,  105, & ! 5   
     &     62,   64,   66,   88,   91,   51,   13,   15,   17,   29, & ! 6   
     &     21,   33,   34,   37,   71,   74,  104,  113,  101,   24, & ! 7   
     &     26,   20,   46,   49,   48,   30,  113,   16,   18,   50, & ! 8   
     &     61,   63,   65,   67,   70,   73,   76,   78,   89,   90, & ! 9   
     &    117,  118,   92,   93,  119,  120,   24,   26,   20,   49, & ! O   
     &     94,   95,   96,   30,   97,   98,   99,  101,   39,  102, & ! 1   
     &    114,  112,  115,  116,  121,   31,  128,  103,  105,   16, & ! 2   
     &     18,   50,   61,   63,   65,   67,   70,   73,   76,   78, & ! 3   
     &     89,   90,  117,  118,   92,   93,  119,  120,   24,   26, & ! 4   
     &     20,   49,   94,   95,   96,   30,   97,   98,   99,  101, & ! 5   
     &     39,  128,  102,  103,  114,  112,  115,  116,  121,   31, & ! 6   
     &     16,   18,   50,   61,   63,   65,   67,   70,   73,   76, & ! 7   
     &     78,   82,   89,   90,  117,  118,   92,   93,  119,  120, & ! 8   
     &     20,   49,   94,   95,   96,   30,   97,   98,   99,  101, & ! 9   
     &     39,  102,  114,  112,  115,  116,  121,   31,   18,   50, & ! O   
     &     61,   63,   65,   67,   70,   73,   76,   78,   82,   89, & ! 1   
     &     90,  117,  118,   92,   93,  119,  120,   20,   49,   94, & ! 2   
     &     95,   96,   30,   97,   98,   99,  101,   39,  102,  114, & ! 3   
     &    112,  115,  116,  121,   31,   49,   16,   18,   50,   61, & ! 4   
     &     63,   65,   67,   70,   73,   76,   78,   82,   89,   92, & ! 5   
     &    117,  119,   20,   49,   94,   95,   96,   30,   97,   98, & ! 6   
     &     99,  101,   39,  102,  114,  112,  115,  116,  121,  115, & ! 7   
     &    115,  116,   31,   31,   31,   90,   90,   90,   93,   93, & ! 8   
     &     93,  118,  118,  118,  120,  120,  120,   84,  130,  131, & ! 9   
     &    132,  132,  132,  132,  133,  134,  134,  134,  134,  134, & ! O   
     &    133,  133,  135,  136,  136,  125,  135,  135,  137,  138, & ! 1   
     &    139,  139,  139,  138,  138,  141,  141,  141,   33,   34, & ! 2   
     &      6,   12,    8,   84,  145,  145,   85,   86,  151,  153, & ! 3   
     &    155,  157,  159,  161,  163,  152,  154,  156,  158,  160, & ! 4   
     &    162,  109,  152,  154,  156,  158,  160,  162,  109,  152, & ! 5   
     &    154,  156,  158,  160,  162,  109,  164,  165,  166,  168, & ! 6   
     &    170,  171,  173,  164,  165,  166,  168,  170,  171,  173, & ! 7   
     &    164,  165,  166,  168,  170,  171,  173,  164,  165,  166, & ! 8   
     &    168,  170,  171,  173,  174,  175,  175,  175,  175,  175, & ! 9   
     &    178,  179,  179,  179,  179,  179,  180,  181,  181,  181, & ! O   
     &    181,  181,  182,  176,  140,  124,  177,  172,  127,  142, & ! 1   
     &    143,  169,  167,  126,  123,  183,  184,   43,  185,  186, & ! 2   
     &    188,  189,  189,  189,  190,  190,  190,  190,  190,  191, & ! 3   
     &    192,  187,  194,  196,  197,  198,  200,  201,  202,  203, & ! 4   
     &    203,  204,  195,  195,  195,  195,  195,  195,  195,  195, & ! 5   
     &    195,  195,  195,  195,  195,  195,  195,  196,  196,  196, & ! 6   
     &    196,  196,  196,  198,  197,    5,    5,    5,    5,    5, & ! 7   
     &      5,    5,    5,    5,  195,  195,  195,  195,  195,  195, & ! 8   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! 9   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! O   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! 1   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! 2   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! 3   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! 4   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! 5   
     &    195,  195,  195,  195,  195,  195,  195,  195,  195,  195, & ! 6   
     &    195,  195,  195,    5,  195,  195,  195,  196,  196,  196, & ! 7   
     &    196,  196,  196,  196,  196,  196,  196,  196,  196,  196, & ! 8   
     &    196,  196,  196,  196,  196,  196,  196,  196,  196,  196, & ! 9   
     &    196,  196,  196,  196,  196,  196,  196,  196,  196,  196, & ! O   
     &    196,  196,  196,  196,  196,  196,  196,  196,  196,  196, & ! 1   
     &    196,  196,  196,  196,  196,  196,  196,  196,  196,  196, & ! 2   
     &    196,  196,  196,  196,  196,  196,  196,  196,  196,  206, & ! 3   
     &    208,  209,  210,  211,  212,  213,  213,  214,  215,  216, & ! 4   
     &    206,  207,  207,  207,  207,  207,    5,  208,  208,  208, & ! 5   
     &    208,  208,  208,  208,  209,  207,  207,  207,  207,  207, & ! 6   
     &    207,  207,  207,  207,  207,  207,  207,  207,  207,  207, & ! 7   
     &    207,  207,  207,  207,  207,  207,  207,  207,  207,  207, & ! 8   
     &    207,  207,  207,  207,  207,  207,  207,  207,  207,  207, & ! 9   
     &    207,  207,  207,  207,  207,    5,    5,    5,  207,  207, & ! O   
     &    218,  220,  221,  222,  223,  224,  225,  226,  227,  228, & ! 1   
     &    229,  230,  231,  218,  218,  219,  219,  219,  219,  219, & ! 2   
     &    226,  227,  225,  226,    5,    5,  220,  220,  221,  220, & ! 3   
     &    220,  222,  223,  222,  221,  220,  231,  224,  218,  228, & ! 4   
     &    232,  233,  234,  235,  236,  237,  238,  207,  207,  219, & ! 5   
     &    207,  196,  196,  196,  208,  195,  195,  195,    5,  195, & ! 6   
     &    195,    5,    5,   53,   53,  240,  241,  240,  240,  241, & ! 7   
     &    200,  200,  202,  202,  202,  202,  203,  203,  203,  203, & ! 8   
     &    210,  210,  212,  212,  213,  213,  213,  213,  217,  223, & ! 9   
     &    223,  226,  226,  227,  227,  222,  222,  230,  230,  231, & ! O   
     &    231,   54,   54,   54,   54,   54,   54,   54,   54,  248, & ! 1   
     &    248/     !  2   

      DATA ( IRR( IRXXN,  2 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    5,   12,    7,    6,    0,    1,    0, & ! 4   
     &      0,    0,    0,   12,   12,   12,    5,    2,    5,   12, & ! 5   
     &     12,    7,    5,    2,    2,    5,    5,    5,   12,    7, & ! 6   
     &      6,    8,    6,    0,    0,   12,    0,    5,    5,    5, & ! 7   
     &      0,    5,    5,    5,    5,    5,    5,    5,    5,    5, & ! 8   
     &      5,    5,    8,    1,    8,    5,   12,    7,    5,    5, & ! 9   
     &     12,    7,    5,   86,    0,    5,    5,    5,    5,    5, & ! O   
     &      5,    5,    5,    5,    5,    5,    5,    5,    5,    5, & ! 1   
     &      5,    5,    5,    5,    5,    5,    5,    5,    5,    5, & ! 2   
     &      5,    5,    5,    5,    5,    5,    5,    5,    5,    5, & ! 3   
     &      5,    5,    5,    5,    5,    5,    1,    1,    1,    1, & ! 4   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 5   
     &      8,    8,    8,    8,    8,    8,    8,    8,    8,    8, & ! 6   
     &      8,    8,    8,    8,    8,    8,    8,    8,    6,    0, & ! 7   
     &      0,    6,    0,    6,    0,    6,    0,    7,    7,    7, & ! 8   
     &      7,    7,    7,    7,    7,    7,    7,    7,    7,    7, & ! 9   
     &      7,    7,    7,    7,    7,    7,    7,    7,    7,    7, & ! O   
     &      7,    7,    7,    7,    7,    7,    7,    7,    7,    7, & ! 1   
     &      7,    7,    7,    7,    7,    7,    6,    6,    6,   12, & ! 2   
     &     12,   12,   12,   12,   12,   12,   12,   12,   12,   12, & ! 3   
     &     12,   12,   12,   12,   12,   12,   12,   12,   12,   12, & ! 4   
     &     12,   12,   12,   12,   12,   12,   12,   12,   12,   12, & ! 5   
     &     12,    1,   12,   12,   12,   12,   12,   12,   12,   12, & ! 6   
     &     16,   16,   16,   16,   16,   16,   16,   16,   16,   16, & ! 7   
     &     16,   16,   16,   16,   16,   16,   16,   16,   16,   16, & ! 8   
     &     16,   16,   16,   16,   16,   16,   16,   16,   16,   16, & ! 9   
     &     16,   16,   16,   16,   16,   16,   16,   16,   20,   20, & ! O   
     &     20,   20,   20,   20,   20,   20,   20,   20,   20,   20, & ! 1   
     &     20,   20,   20,   20,   20,   20,   20,   20,   20,   20, & ! 2   
     &     20,   20,   20,   20,   20,   20,   20,   20,   20,   20, & ! 3   
     &     20,   20,   20,   20,   20,   49,    8,    8,    8,    8, & ! 4   
     &      8,    8,    8,    8,    8,    8,    8,    8,    8,    8, & ! 5   
     &      8,    8,    8,    8,    8,    8,    8,    8,    8,    8, & ! 6   
     &      8,    8,    8,    8,    8,    8,    8,    8,    8,  115, & ! 7   
     &    116,  116,    8,   49,   31,   89,   92,   82,   89,   92, & ! 8   
     &     82,   89,   92,   82,   89,   92,   82,    5,    5,    5, & ! 9   
     &      5,    1,    8,    0,    5,    7,    8,   12,   16,   20, & ! O   
     &      1,    8,    5,    7,   12,    5,    1,    8,    5,    8, & ! 1   
     &     12,    7,    8,    1,    5,   12,    8,    7,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,  147,    0,    0,    5,    5, & ! 3   
     &      5,    5,    5,    5,    5,    7,    7,    7,    7,    7, & ! 4   
     &      7,    7,    8,    8,    8,    8,    8,    8,    8,   12, & ! 5   
     &     12,   12,   12,   12,   12,   12,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    7,    7,    7,    7,    7,    7,    7, & ! 7   
     &      8,    8,    8,    8,    8,    8,    8,   12,   12,   12, & ! 8   
     &     12,   12,   12,   12,    5,   12,    7,    8,   16,   20, & ! 9   
     &      5,   12,    7,    8,   16,   20,    5,   12,    7,    8, & ! O   
     &     16,   20,    5,    5,    5,    5,    5,    5,    5,    5, & ! 1   
     &      5,    5,    5,    5,    5,    5,    5,    5,    5,    0, & ! 2   
     &      0,    5,    1,    8,   12,    7,    8,   16,   20,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    1,    0,   12,   12,    4,  199,  199,  197, & ! 5   
     &    198,  200,    7,    6,    8,  201,  203,   12,    7,    6, & ! 6   
     &    196,  196,  196,    0,    7,  194,  196,  196,  205,  200, & ! 7   
     &    202,  203,  197,  198,    0,   57,   58,   60,  163,   13, & ! 8   
     &     15,   17,   19,  106,  107,   87,   69,  111,   23,   25, & ! 9   
     &     27,   22,   28,   29,   32,   33,   34,  108,   40,   41, & ! O   
     &     43,   42,   83,   83,  110,   44,   46,   48,  113,   45, & ! 1   
     &     51,   86,   80,  125,  137,   62,   64,   66,   79,   88, & ! 2   
     &     91,   68,   21,   35,   36,  100,  138,  132,  133,   38, & ! 3   
     &     72,   75,   37,   71,   74,  104,  180,  135,  161,  159, & ! 4   
     &    157,  155,  153,  151,  131,  178,  174,  172,  177,  140, & ! 5   
     &    176,  124,  182,  142,  127,  169,  143,  167,  126,  123, & ! 6   
     &    183,  184,  130,  204,  204,   77,  189,   16,   18,   50, & ! 7   
     &     61,   63,   65,   67,   70,   73,   76,   82,   89,   90, & ! 8   
     &    117,  118,   92,   93,  119,  120,   24,   26,   20,   49, & ! 9   
     &     94,   95,   96,   30,   97,   98,   99,  101,   39,  102, & ! O   
     &    114,  112,  115,  116,  121,   31,  134,  136,  139,  141, & ! 1   
     &    152,  154,  156,  158,  160,  162,  109,  164,  165,  166, & ! 2   
     &    168,  170,  171,  173,  175,  179,  181,   78,  190,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      5,    1,   12,    6,    8,  213,  217,    5,   12,    7, & ! 5   
     &      6,    8,  208,  208,    1,   13,   15,   17,  106,  107, & ! 6   
     &     87,   23,   25,   28,   33,   34,  108,   83,   83,  110, & ! 7   
     &     62,   64,   66,   79,   88,   91,   68,   32,   21,   35, & ! 8   
     &     36,  100,  138,  132,   29,  133,   38,   72,   75,   37, & ! 9   
     &     71,   74,  104,  180,  135,  214,  216,  215,   77,  189, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    5,    8,    1,   12,    7,    6,  227, & ! 2   
     &      0,    0,  225,  226,  223,  224,   12,    7,    7,    6, & ! 3   
     &    220,    0,    8,    0,  221,  221,    0,    8,    2,    5, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  197,  198,  208, & ! 5   
     &    220,  208,  208,  220,  220,  214,  215,  239,  239,  237, & ! 6   
     &    238,  237,  238,    0,    0,    0,    0,  242,  243,  244, & ! 7   
     &    243,  245,    0,    0,  243,  245,    0,    0,  243,  245, & ! 8   
     &    243,  245,  243,  245,    0,    0,  243,  245,    0,  243, & ! 9   
     &    245,  243,  245,  243,  245,    0,    0,    0,    0,    0, & ! O   
     &      0,    4,    1,   40,   44,    6,  246,   13,   11,    5, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN,  3 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,  247,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN,  4 ), IRXXN = 1, NRXNS ) / & 
     &      2,    3,    5,    2,    7,    2,    5,    5,    5,   14, & ! O   
     &     12,   12,   12,   16,   16,   12,   12,   12,   16,   18, & ! 1   
     &     12,    5,   16,   14,   13,   12,   12,   12,   12,   38, & ! 2   
     &     39,    5,    5,    5,    5,    5,   12,   20,   16,   49, & ! 3   
     &     50,    6,   28,   12,    5,    6,    8,    1,    0,    2, & ! 4   
     &      2,    5,   12,    0,    4,    4,   12,    6,    9,    6, & ! 5   
     &     10,    6,    6,    7,    8,   10,    8,   12,    5,    6, & ! 6   
     &      7,    6,   53,    6,   10,   11,   12,    6,   12,   12, & ! 7   
     &     16,   18,   50,   61,   63,   65,   67,    5,   70,   73, & ! 8   
     &     76,   78,   78,    5,    6,   82,   83,   80,   82,   34, & ! 9   
     &      5,   32,    5,   14,   10,   29,   89,   92,   24,   26, & ! O   
     &     12,   20,   49,   94,   95,   96,   12,   30,   98,   20, & ! 1   
     &     12,   20,   12,   12,   12,  101,   59,   59,  105,   12, & ! 2   
     &     12,   12,   12,    5,    5,    5,   30,   12,   16,    5, & ! 3   
     &     31,   31,    6,   50,   52,   28,    5,    5,    5,    5, & ! 4   
     &      5,    5,   52,    5,    5,    5,    5,    5,    5,  114, & ! 5   
     &    115,  115,  115,  117,  119,   52,   12,   20,   49,   13, & ! 6   
     &     12,   12,   20,  101,   59,   59,  105,   30,   46,  122, & ! 7   
     &    122,   46,   20,   48,   49,  113,   30,   12,   12,   12, & ! 8   
     &     12,   12,   12,   12,   45,   45,   45,   45,   51,    5, & ! 9   
     &      6,    6,   51,    5,   51,    6,   51,   51,   16,   18, & ! O   
     &     20,   12,   12,   16,    6,   12,   12,   39,  128,   12, & ! 1   
     &    105,    6,    6,    6,    6,    6,   45,   45,   45,   40, & ! 2   
     &     41,   41,   41,   41,   41,   41,   41,   41,   41,   41, & ! 3   
     &     19,  122,   51,  122,   42,  122,   51,  122,   42,   42, & ! 4   
     &      5,    5,    5,   41,   41,  110,  110,   41,   41,   43, & ! 5   
     &    111,   39,   41,   74,   41,   41,   45,   45,   45,   41, & ! 6   
     &     12,   12,   12,   12,   12,   12,   12,   13,   13,   13, & ! 7   
     &     13,   12,   13,   12,   51,   12,  106,   12,   51,   12, & ! 8   
     &     12,   12,   12,   12,   12,   12,   12,   12,   12,   12, & ! 9   
     &     12,   12,   12,   13,   12,   12,   12,   12,   12,   12, & ! O   
     &     12,   12,   12,   12,   16,   16,   16,   16,   12,   23, & ! 1   
     &      5,    6,    6,   25,    5,   16,   16,   16,   16,   16, & ! 2   
     &     12,   12,  111,   12,   12,   12,   16,   16,   12,   12, & ! 3   
     &     16,   12,   16,   12,   16,   18,   12,   12,   12,   12, & ! 4   
     &     12,   17,   12,    6,    6,    6,    6,   12,    6,   12, & ! 5   
     &      6,    6,   16,   18,   20,   12,   12,   13,    6,   12, & ! 6   
     &     12,   39,  128,   12,    6,    6,   12,    6,    6,   12, & ! 7   
     &     12,    6,    6,   18,    0,  122,  122,  122,  122,  122, & ! 8   
     &    122,  122,  122,  122,  122,  122,  122,    5,   18,   18, & ! 9   
     &     30,   14,   13,   14,  134,   12,   12,   41,  106,   16, & ! O   
     &    132,  115,   35,   45,   41,   22,    5,    6,   28,  139, & ! 1   
     &    127,  126,  126,  126,  141,  142,  126,  143,  144,  144, & ! 2   
     &      9,    0,   10,  145,  146,  148,  149,  150,  152,  154, & ! 3   
     &    156,  158,  160,  162,  109,  164,  165,  166,  168,  170, & ! 4   
     &    171,  173,  164,  165,  166,  168,  170,  171,  173,  126, & ! 5   
     &    167,  169,  142,  172,  140,   41,   12,   12,   12,   12, & ! 6   
     &     12,   12,   12,  167,  143,  142,  172,  140,  140,   45, & ! 7   
     &      6,    6,    6,    6,    6,    6,    6,  143,  142,  172, & ! 8   
     &    140,  140,  140,  167,  175,  123,  123,    6,   13,   16, & ! 9   
     &    179,  126,  126,    6,   13,   16,  181,  126,  123,    6, & ! O   
     &     13,   16,  182,  182,  182,  182,  182,  182,  182,  182, & ! 1   
     &    182,  182,  182,  182,  182,  176,  124,  182,   16,  187, & ! 2   
     &    187,  190,   37,    6,  126,    6,    6,   12,   12,    9, & ! 3   
     &      9,  193,  195,  195,  196,  195,    5,  195,  195,  195, & ! 4   
     &    196,  195,  196,  205,  196,  205,  205,  196,  194,  196, & ! 5   
     &    194,  194,  201,  202,  196,  194,  194,  200,  195,  203, & ! 6   
     &    194,  197,  198,  196,  196,  195,  205,  195,  195,  196, & ! 7   
     &    200,  200,  200,  200,  205,  205,  205,  205,  205,  205, & ! 8   
     &    205,  205,  205,  205,  205,  205,  205,  205,  205,  205, & ! 9   
     &    205,  205,  205,  205,  205,  205,  205,  205,  205,  205, & ! O   
     &    205,  205,  205,  205,  205,  205,  205,  205,  205,  205, & ! 1   
     &    205,  205,  205,  205,  205,  205,  205,  205,  204,  205, & ! 2   
     &    205,  205,  205,  205,  205,  205,  205,  205,  205,  205, & ! 3   
     &    205,  205,  205,  205,  205,  205,  205,  205,  205,  205, & ! 4   
     &    205,  205,  205,  205,  205,  205,  205,  205,  205,  205, & ! 5   
     &    205,  205,  205,  205,  205,  205,  205,  205,  205,  205, & ! 6   
     &    205,  205,  205,  195,  195,  205,  205,  199,  199,   12, & ! 7   
     &     12,   12,   12,   12,  123,  123,  126,   13,   12,    5, & ! 8   
     &    199,  199,   12,    5,  199,  199,  199,  199,   16,   18, & ! 9   
     &     20,   12,   12,   16,  199,   12,   12,   39,  128,   12, & ! O   
     &    105,  199,  199,  199,  199,  199,   12,  199,  199,  143, & ! 1   
     &    164,  165,  166,  168,  170,  171,  173,  167,  143,  142, & ! 2   
     &    172,  140,  140,  199,  123,  126,  123,  199,  199,  207, & ! 3   
     &    207,  208,  207,  207,  207,  207,  208,  207,  207,  207, & ! 4   
     &    210,  208,  217,  212,  208,  206,  207,  207,  210,  207, & ! 5   
     &    213,  209,  207,  206,  208,  217,  217,  217,  217,  217, & ! 6   
     &    217,  217,  217,  217,  217,  217,  217,  217,  217,  217, & ! 7   
     &    217,  217,  217,  216,  217,  217,  217,  217,  217,  217, & ! 8   
     &    217,  217,  217,  217,  217,  217,  217,  217,  217,  217, & ! 9   
     &    217,  217,  217,  217,  217,  216,  207,  207,  217,  217, & ! O   
     &    219,  219,  219,  219,  219,  219,  219,  219,  219,  219, & ! 1   
     &    219,  220,  221,  223,  219,  220,  224,  225,  226,  218, & ! 2   
     &    219,  220,  218,  218,  220,  219,  223,  219,  220,  227, & ! 3   
     &    219,  221,  220,  220,  231,  230,  221,  219,  220,  219, & ! 4   
     &    207,  219,  219,  219,  219,  207,  207,  208,  232,  220, & ! 5   
     &    219,  207,  195,  219,  207,  205,  205,  205,   12,  205, & ! 6   
     &    205,  207,  207,   10,   10,   10,   10,  202,  202,  202, & ! 7   
     &    194,  232,  242,  243,  194,  232,  200,  200,  194,  232, & ! 8   
     &    232,  206,  232,  206,  210,  210,  232,  206,  245,  233, & ! 9   
     &    234,  233,  234,  233,  234,    0,    0,    0,    0,    0, & ! O   
     &      0,  147,  147,  147,  147,  147,  147,  248,  147,   12, & ! 1   
     &     54/     !  2   

      DATA ( IRR( IRXXN,  5 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    7,    0,    6,    7,    6,   12,    0, & ! O   
     &     14,   16,   18,   20,   14,   20,   24,   26,   18,   20, & ! 1   
     &     20,   12,   30,    0,   14,   14,   20,   20,   20,   14, & ! 2   
     &     14,   12,   12,   12,   12,   16,    6,    6,    8,    6, & ! 3   
     &      8,   22,    6,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    5, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    6,    6,    0, & ! 6   
     &      6,    0,    0,    8,    0,    0,    6,    0,   55,   47, & ! 7   
     &      0,    0,   59,   59,    0,    0,    0,   12,   71,   74, & ! 8   
     &     74,   74,   10,   12,   80,   13,   12,   13,    0,   33, & ! 9   
     &     86,   86,   86,    6,   87,   32,   90,   93,   49,   49, & ! O   
     &     14,    0,    0,    0,    0,    0,   34,   97,    0,   99, & ! 1   
     &     14,   14,   14,   14,   30,    0,   12,   12,    0,   13, & ! 2   
     &     15,   17,   15,   16,   50,  109,    0,   47,  112,   20, & ! 3   
     &      8,    8,   28,    6,    6,    8,   12,   12,   12,   24, & ! 4   
     &     25,   17,    6,   12,   12,   12,   12,   12,   12,    0, & ! 5   
     &    116,  116,  116,  118,  120,    6,   14,   10,   10,   30, & ! 6   
     &     31,   14,   14,   10,  103,  103,   10,    6,    0,    0, & ! 7   
     &      0,    0,    6,    0,    6,    0,    6,    6,    6,   16, & ! 8   
     &     16,    6,    6,    6,  123,  123,  126,  126,   12,    6, & ! 9   
     &     51,    5,    6,    6,   45,    5,  122,    6,    6,    6, & ! O   
     &      6,    6,   20,   20,   12,   20,    6,    6,    6,    6, & ! 1   
     &      6,   33,   12,   13,   33,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,  124,  124,   43,   43, & ! 3   
     &     25,    0,   23,    0,   25,    0,    5,    0,    5,    5, & ! 4   
     &     16,   18,   20,    0,    0,    0,    0,    0,    0,  111, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &     13,   13,   16,   16,   13,   13,   13,   12,   12,   12, & ! 7   
     &     12,   13,   25,   13,   13,    6,   25,   13,   13,   13, & ! 8   
     &     16,   16,   20,   13,   13,   20,   13,   20,   14,   39, & ! 9   
     &    128,   13,  105,   12,   13,    6,    6,   13,   16,   16, & ! O   
     &     16,   16,   16,   16,   12,   12,   12,   12,   16,   17, & ! 1   
     &     16,   23,   16,   13,   16,    6,    6,   47,   18,   20, & ! 2   
     &     16,   16,   16,   13,   16,   16,   39,  128,   16,   16, & ! 3   
     &     33,   16,    6,   16,   47,   47,   13,    6,   16,   16, & ! 4   
     &      6,   13,   17,   12,   12,   12,   12,    6,   12,    6, & ! 5   
     &     23,   25,    6,    6,    6,    6,    6,   20,   12,   20, & ! 6   
     &      6,    6,    6,    6,  105,   33,    6,   13,   33,   45, & ! 7   
     &      6,   13,    0,   47,    0,   87,   87,   87,   87,   87, & ! 8   
     &     87,   87,   87,   87,   87,   87,   87,    0,   59,   59, & ! 9   
     &     97,   12,   30,   12,   21,    6,    6,    0,   13,  111, & ! O   
     &     14,  116,   12,    6,  125,   33,   17,   35,   17,    0, & ! 1   
     &      0,    6,    6,  140,    0,    0,    0,  126,    0,    0, & ! 2   
     &     10,    0,    0,    0,    0,    0,    0,   10,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,  123,  126,  167,  169,  142, & ! 4   
     &    172,   45,    6,    6,    6,    6,    6,    6,    6,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,  126,  167,  169,  142, & ! 6   
     &    172,  140,  123,    6,    6,    6,    6,    6,    6,    6, & ! 7   
     &    126,  167,  169,  142,  172,  140,   22,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,   12,  143,  167,   12,   12,   12, & ! 9   
     &     12,  142,  143,   12,   12,   12,   12,  143,  167,   12, & ! O   
     &     12,   12,   27,  176,  176,  176,  176,  176,  176,  176, & ! 1   
     &    176,  176,  176,  176,  176,  140,  177,  176,    0,   10, & ! 2   
     &     10,    0,   13,  190,    0,   12,   12,   13,   13,    6, & ! 3   
     &      6,   69,    0,    2,    2,  199,  195,    7,    6,    8, & ! 4   
     &      6,   14,    0,   12,    5,    0,   12,    0,    0,    0, & ! 5   
     &    199,    5,    0,    0,    6,    7,    8,    0,    6,    0, & ! 6   
     &      0,  195,    0,    0,    6,  200,    0,   12,    0,    0, & ! 7   
     &      6,    8,    0,  199,   16,   18,   50,   61,  109,   12, & ! 8   
     &     20,   49,   94,   12,   12,   12,   12,   16,   24,   26, & ! 9   
     &     95,   96,   12,   30,   98,   12,   20,   12,    5,    5, & ! O   
     &      5,    5,   82,   34,   30,    5,   31,   31,    6,   50, & ! 1   
     &     52,   14,    5,   22,   28,   63,   65,   67,   82,   89, & ! 2   
     &     92,    5,   20,   12,   12,   12,  141,   30,  134,   70, & ! 3   
     &     73,   76,  101,   59,   59,  105,  181,   35,  162,  160, & ! 4   
     &    158,  156,  154,  152,   18,  179,  175,  182,  182,  182, & ! 5   
     &    182,  182,  182,  182,  182,  182,  182,  182,  182,  182, & ! 6   
     &    176,  124,   18,   14,   14,   78,  190,   13,   15,   16, & ! 7   
     &     16,  199,  199,  199,  124,  124,  127,   12,  199,  199, & ! 8   
     &      6,    5,  199,  199,   51,    5,  122,   12,  199,  199, & ! 9   
     &    199,  199,   20,   20,   12,   20,  199,  199,  199,  199, & ! O   
     &    199,   33,   12,    6,    6,    0,  199,  125,  126,  126, & ! 1   
     &    123,  126,  167,  169,  142,  172,  199,  199,  199,  199, & ! 2   
     &    199,  199,  199,   22,  167,  143,  167,  126,   12,    0, & ! 3   
     &      2,    2,    5,    7,    6,    8,    6,   12,   12,   14, & ! 4   
     &    207,    0,    0,    0,    6,    8,    0,   12,    0,    6, & ! 5   
     &      0,    6,    0,    0,    0,   12,   20,   49,   13,   15, & ! 6   
     &     17,   24,   26,   34,   12,   20,   17,   82,   34,   30, & ! 7   
     &     63,   65,   67,   82,   89,   92,    5,   98,   20,   12, & ! 8   
     &     12,   12,  141,   30,   30,  134,   70,   73,   76,  101, & ! 9   
     &     59,   59,  105,  181,   35,  207,   14,    0,   78,  190, & ! O   
     &      0,    2,    0,  221,    5,   12,    7,    6,    8,   16, & ! 1   
     &     12,  221,    0,  219,  227,    0,    0,    0,    0,    8, & ! 2   
     &      6,    6,    7,    6,    0,    0,    0,    6,    6,    0, & ! 3   
     &    221,  219,   10,    0,    0,    0,    0,   10,  219,   13, & ! 4   
     &    195,  195,  207,  207,  195,  195,  195,  196,  199,  207, & ! 5   
     &    208,  197,  207,  197,  221,   12,  207,   12,  195,  207, & ! 6   
     &    207,  195,  195,  240,  241,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,  191,  192,    9,    9,   10,   10,  192,  192, & ! 8   
     &      0,    0,    9,    9,   10,   10,   10,   10,    0,    0, & ! 9   
     &      0,    9,    9,   10,   10,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,  106,  111,    9,  246,    0,  192,   54, & ! 1   
     &     13/     !  2   

      DATA ( IRR( IRXXN,  6 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    6,    0, & ! O   
     &      0,   14,   14,    0,    0,   14,   14,   14,   20,   14, & ! 1   
     &     13,   20,   14,    0,    0,    0,   14,   31,   31,    0, & ! 2   
     &     12,   13,   17,   17,   17,    0,   17,    0,   47,    0, & ! 3   
     &      0,   21,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,   10,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,   56,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,   14,   12,   12, & ! 8   
     &     12,   12,    0,   16,   13,   29,    5,   12,    0,   32, & ! 9   
     &     85,   80,   84,   12,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,   17,   27,   14,    0,  102,  102,    0,    0, & ! 2   
     &      0,   15,    0,   13,   31,   31,    0,    0,   47,   31, & ! 3   
     &     13,   13,   14,    0,   23,    0,   14,   16,   16,    4, & ! 4   
     &     20,   13,   23,   20,   20,   16,   49,   49,   14,    0, & ! 5   
     &      0,    0,    0,    0,    0,   23,   10,    0,    0,   31, & ! 6   
     &     14,   10,   10,    0,  102,  102,    0,   10,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,   13,   15,   18, & ! 8   
     &     18,   13,   13,   15,  124,  124,  127,  127,    6,  122, & ! 9   
     &     45,  122,   25,  122,    6,  122,    6,   12,   47,   47, & ! O   
     &     13,   13,   31,    6,   13,   31,   14,   47,    0,   28, & ! 1   
     &      0,   12,   45,   17,   41,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,  127,  127, & ! 3   
     &     13,    0,    6,    0,    5,    0,   25,    0,  122,   12, & ! 4   
     &     47,   47,   13,    0,    0,    0,    0,    0,    0,    1, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &    106,   15,   18,   18,   17,   17,   17,  106,  106,  106, & ! 7   
     &    106,   29,   23,  106,  106,  106,   13,  106,  106,    6, & ! 8   
     &     47,   47,   13,   35,   35,   14,   28,   31,   13,   13, & ! 9   
     &     13,   28,   13,   33,   45,   13,   13,    0,   15,   18, & ! O   
     &     18,   13,   13,   17,  111,  111,  111,  111,   13,   42, & ! 1   
     &    111,   16,  111,   12,  111,   25,  111,    0,   47,   13, & ! 2   
     &     13,   35,   20,   28,   20,   14,   47,   47,   28,  105, & ! 3   
     &     47,   45,   13,    6,    0,    0,    6,   15,   31,   18, & ! 4   
     &     13,   12,   22,   37,   37,   37,   37,   13,   23,   25, & ! 5   
     &     51,   51,   47,    0,   13,   13,   35,   14,   13,   31, & ! 6   
     &     14,   47,    0,   28,    0,   12,   45,   17,   41,    0, & ! 7   
     &     13,   17,    0,    0,    0,   23,   25,   13,   23,   25, & ! 8   
     &     13,   23,   25,   13,   23,   25,   13,    0,    0,    0, & ! 9   
     &      0,    5,   31,   62,   12,  132,  132,    0,   12,   12, & ! O   
     &     12,  132,  136,  125,    5,   12,    0,  100,   12,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    6,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    6,    6,    6,    6,    6, & ! 4   
     &      6,    6,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,  126,  167,  169,  142,  172,   12,   22, & ! 7   
     &     12,   12,   12,   12,   12,   12,   12,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,  123,  176,  177,   37,  106,  111, & ! 9   
     &    126,  176,  176,  123,  106,  111,  126,  182,  182,  123, & ! O   
     &    106,  111,   22,  124,  140,   43,  140,  140,  124,  140, & ! 1   
     &    140,  140,  140,  140,  124,  124,  127,  124,    0,    0, & ! 2   
     &      0,    0,   14,    0,    0,   13,   13,   37,   37,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,   12,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,   59,   59,    0,   14, & ! 8   
     &      0,    0,    0,   13,   15,   17,    0,  112,   49,   49, & ! 9   
     &      0,    0,   34,   97,    0,   14,   14,   17,   16,   50, & ! O   
     &    109,  109,    0,   33,    0,   20,    8,    8,   28,    6, & ! 1   
     &      6,    6,   86,   33,   17,    0,    0,    0,   13,   90, & ! 2   
     &     93,   12,   99,   14,   14,   30,    0,   97,   21,   71, & ! 3   
     &     74,   74,    0,   12,   12,    0,   12,   12,    0,    0, & ! 4   
     &      0,    0,    0,    0,   59,   12,   12,  176,  176,  176, & ! 5   
     &    176,  176,   27,  176,  176,  176,  176,  176,  176,  176, & ! 6   
     &    140,  177,   59,    0,  205,   74,    0,   12,   12,   18, & ! 7   
     &     18,   13,   13,   15,  199,  199,  199,  199,   23,  205, & ! 8   
     &     51,  205,   25,  205,   45,  205,    0,   13,    0,    0, & ! 9   
     &     13,   13,   31,  199,   13,   31,   14,    0,    0,   28, & ! O   
     &      0,   12,   45,   13,   33,    0,  132,   12,    6,  199, & ! 1   
     &    199,  199,  199,  199,  199,  199,    0,  126,  167,  169, & ! 2   
     &    142,  172,   12,   12,  177,  176,  182,  127,   13,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   12, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,   14,    0,    0,   12,   12, & ! 6   
     &     15,   49,   49,   12,   14,   14,   12,    0,   33,    0, & ! 7   
     &      0,    0,    0,   13,   90,   93,   12,    0,   99,   14, & ! 8   
     &     14,   30,    0,   97,   97,   21,   71,   74,   74,    0, & ! 9   
     &     12,   12,    0,   12,   12,   12,    0,    0,   74,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &    222,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,   12,   12,   12,   12,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,  207,    0,  195,  207,  195, & ! 6   
     &    195,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,  247,    0,    0,   69, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN,  7 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    8,    0, & ! O   
     &      0,    0,    0,    0,    0,   13,    0,    0,   14,    0, & ! 1   
     &      0,   30,   21,    0,    0,    0,    0,   14,   14,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,   22,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,   33,    0,    0, & ! 8   
     &      0,    0,    0,   20,   81,   34,    0,    6,    0,   84, & ! 9   
     &      0,   13,    6,   13,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,   22,   33,   33,    0,  103,  103,    0,    0, & ! 2   
     &      0,    0,    0,    0,   17,   17,    0,    0,    0,   13, & ! 3   
     &      0,    0,    0,    0,   22,    0,   13,   18,   18,   23, & ! 4   
     &     13,   40,   22,   14,   31,   20,   31,   31,   47,    0, & ! 5   
     &      0,    0,    0,    0,    0,   22,    0,    0,    0,   34, & ! 6   
     &     13,    0,    0,    0,  121,  121,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   31, & ! 8   
     &     31,   15,   15,   17,    6,    6,    6,    6,   23,    0, & ! 9   
     &     23,    0,   12,    0,   25,    0,    0,   13,    0,    0, & ! O   
     &      0,   35,    6,   14,   28,    6,   13,    0,    0,   33, & ! 1   
     &      0,    0,    0,   22,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &     42,    0,    5,    0,   12,    0,    6,    0,    0,   13, & ! 4   
     &    111,  111,   41,    0,    0,    0,    0,    0,    0,   39, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,  106,   31,   31,  106,   22,   22,   37,   37,   37, & ! 7   
     &     37,   32,   17,  122,   12,   13,   12,  122,   12,  106, & ! 8   
     &     13,   13,  106,  106,  106,   13,  106,   13,   17,   47, & ! 9   
     &      0,   33,    0,    0,    0,   17,   33,    0,  111,   31, & ! O   
     &     31,   15,   17,   22,   37,   37,   37,   37,   29,   19, & ! 1   
     &    122,   51,  122,   16,  122,   51,  122,    0,    0,  111, & ! 2   
     &     35,  111,   14,   16,   31,   13,    0,    0,   33,   47, & ! 3   
     &      0,    0,   17,   33,    0,    0,    0,    0,   18,   31, & ! 4   
     &     17,    6,    6,   33,   33,   33,   33,   29,   17,   13, & ! 5   
     &     17,   12,    0,    0,    0,   35,    0,    6,   28,    6, & ! 6   
     &     13,    0,    0,   33,    0,    0,    0,   22,    0,    0, & ! 7   
     &     17,   22,    0,    0,    0,    5,    5,   32,    5,    5, & ! 8   
     &     32,    6,    6,   32,    6,    6,   32,    0,    0,    0, & ! 9   
     &      0,   13,   34,   20,    0,   13,   13,    0,  132,  132, & ! O   
     &      5,    0,    0,   12,   12,   16,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,   12,   12,   12,   12,   12,    0,   12, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    6,   33,   37,   37, & ! 9   
     &      0,    0,    6,   33,  123,  123,    0,    0,    6,   33, & ! O   
     &    123,  123,   13,   27,  124,   27,  124,  124,  177,  124, & ! 1   
     &    124,  124,  124,  124,  177,  177,  143,  127,    0,    0, & ! 2   
     &      0,    0,   12,    0,    0,   37,   37,  106,   16,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,   15,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,   13,   31, & ! O   
     &     31,   31,    0,   32,    0,   31,   13,   13,    0,    0, & ! 1   
     &     23,   12,   84,   12,   12,    0,    0,    0,   29,    0, & ! 2   
     &      0,   14,    0,   17,   27,   14,    0,    0,   12,   12, & ! 3   
     &     12,   12,    0,  102,  102,    0,  126,  136,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,  126,  123,  140,  140,  140, & ! 5   
     &    124,   43,   22,  140,  124,  140,  140,  140,  140,  124, & ! 6   
     &    124,  127,    0,    0,    0,   12,    0,    0,    0,   31, & ! 7   
     &     31,   17,   15,   17,   12,   12,   12,   29,   17,  122, & ! 8   
     &     45,  122,   13,  122,    6,  122,    0,   44,    0,    0, & ! 9   
     &      0,   35,  199,   14,   28,  199,   13,    0,    0,   33, & ! O   
     &      0,    0,    0,   17,   41,    0,   13,   16,    0,  205, & ! 1   
     &    205,  205,  205,  205,  205,  205,    0,   12,   12,   12, & ! 2   
     &     12,   12,  205,    0,  199,  199,  199,   12,   37,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &     12,    0,    0,    0,    0,   12,    0,    0,   32,    0, & ! 7   
     &      0,    0,    0,   29,    0,    0,   14,    0,    0,   17, & ! 8   
     &     27,   14,    0,    0,    0,   12,   12,   12,   12,    0, & ! 9   
     &    102,  102,    0,  126,  136,    0,    0,    0,   12,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,  207,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN,  8 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,   17,    0,    0,    0,    0, & ! 1   
     &      0,   31,    0,    0,    0,    0,    0,   33,   33,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,   69,    0,    0, & ! 8   
     &      0,    0,    0,   30,   32,    0,    0,   29,    0,   20, & ! 9   
     &      0,    6,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,   33,   34,   34,    0,  104,  104,    0,    0, & ! 2   
     &      0,    0,    0,    0,   22,   22,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,   13,    0,   69,   50,   50,  111, & ! 4   
     &     92,    4,   13,   34,   14,   31,   14,   14,   33,    0, & ! 5   
     &      0,    0,    0,    0,    0,   13,    0,    0,    0,   10, & ! 6   
     &     17,    0,    0,    0,   10,   10,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    6, & ! 8   
     &      6,    0,   17,   19,   12,   12,   12,   12,   17,    0, & ! 9   
     &     17,    0,   13,    0,   12,    0,    0,   44,    0,    0, & ! O   
     &      0,    0,   17,   13,    0,   13,   17,    0,    0,   41, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &     12,    0,   17,    0,   13,    0,   12,    0,    0,   44, & ! 4   
     &     44,   44,    0,    0,    0,    0,    0,    0,    0,    5, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,  107,   13,   13,  108,  106,  106,   33,   33,   33, & ! 7   
     &     33,  106,   42,    0,    6,  122,    0,    0,   25,  122, & ! 8   
     &    111,  111,   87,   87,   87,  111,   87,   17,   22,    0, & ! 9   
     &      0,   41,    0,    0,    0,   22,   41,    0,   47,   13, & ! O   
     &     13,  111,   22,  111,   33,   33,   33,   33,   32,   25, & ! 1   
     &     47,   17,   47,   47,   47,   12,   47,    0,    0,   47, & ! 2   
     &    111,   47,   13,  111,   13,   17,    0,    0,   41,    0, & ! 3   
     &      0,    0,   22,   41,    0,    0,    0,    0,    6,    6, & ! 4   
     &      0,   27,   19,  125,   34,   34,   34,   32,   19,    0, & ! 5   
     &     22,    0,    0,    0,    0,    0,    0,    0,    0,   13, & ! 6   
     &     17,    0,    0,   41,    0,    0,    0,    0,    0,    0, & ! 7   
     &     22,   45,    0,    0,    0,   12,   12,    5,   12,   12, & ! 8   
     &      5,   12,   12,    6,   12,   12,    6,    0,    0,    0, & ! 9   
     &      0,   33,   10,    5,    0,  135,  135,    0,   64,   13, & ! O   
     &     13,    0,    0,   16,   16,   14,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,   12,   34,   33,   33, & ! 9   
     &      0,    0,   12,   34,   33,   33,    0,    0,   12,   34, & ! O   
     &     33,   33,   15,   22,  177,   22,  177,  177,  127,  177, & ! 1   
     &    177,  177,  177,  177,  127,  127,  169,   43,    0,    0, & ! 2   
     &      0,    0,    5,    0,    0,    0,    0,    0,  111,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   17, & ! O   
     &     17,   17,    0,   84,    0,   13,    0,    0,    0,    0, & ! 1   
     &     22,   13,    6,   16,    0,    0,    0,    0,   34,    0, & ! 2   
     &      0,   33,    0,   22,   33,   33,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,  103,  103,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  124,  124,  124, & ! 5   
     &     27,   27,   13,  124,  177,  124,  124,  124,  124,  177, & ! 6   
     &    177,  143,    0,    0,    0,    0,    0,    0,    0,  199, & ! 7   
     &    199,    0,   17,   19,   37,   37,   37,   32,   19,    0, & ! 8   
     &     23,    0,    0,    0,   25,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,   17,   13,    0,   13,   17,    0,    0,   41, & ! O   
     &      0,    0,    0,   22,    0,    0,  135,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,  205,  205,  205, & ! 2   
     &    205,  205,    0,    0,   12,   12,   12,   37,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,   84,    0, & ! 7   
     &      0,    0,    0,   34,    0,    0,   33,    0,    0,   22, & ! 8   
     &     33,   33,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &    103,  103,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN,  9 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,   22,    0,    0,    0,    0, & ! 1   
     &      0,   14,    0,    0,    0,    0,    0,   34,   34,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    4,   29,    0,    0,   32,    0,   29, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,   34,   41,   41,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,   17,    0,    0,   61,   20,   12, & ! 4   
     &     93,   69,   17,   69,   13,   14,   47,   47,   35,    0, & ! 5   
     &      0,    0,    0,    0,    0,   17,    0,    0,    0,    6, & ! 6   
     &     45,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   15, & ! 8   
     &     13,    0,   19,   22,   37,   37,   37,   37,   19,    0, & ! 9   
     &     22,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,   34,    0,    0,   17,   33,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &     23,    0,   22,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   47, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,   17,   17,    0,   87,   87,  125,   34,   34, & ! 7   
     &     34,   87,  106,    0,   23,    0,    0,    0,    6,    0, & ! 8   
     &      0,    0,  111,    0,    0,   16,    0,   34,   33,    0, & ! 9   
     &      0,    0,    0,    0,    0,  106,   45,    0,    0,   17, & ! O   
     &     17,   47,  111,   47,  125,   34,   34,   34,  111,   13, & ! 1   
     &      0,   22,    0,    0,    0,   47,    0,    0,    0,    0, & ! 2   
     &     47,    0,   47,   47,   17,   22,    0,    0,   47,    0, & ! 3   
     &      0,    0,  111,   45,    0,    0,    0,    0,   15,   13, & ! 4   
     &      0,   15,   15,   36,  125,  125,  125,   33,   25,    0, & ! 5   
     &     13,    0,    0,    0,    0,    0,    0,    0,    0,   17, & ! 6   
     &     22,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &     45,    0,    0,    0,    0,  129,  129,   12,  129,  129, & ! 8   
     &     12,  129,  129,   12,  129,  129,   12,    0,    0,    0, & ! 9   
     &      0,   47,    6,   13,    0,   45,    0,    0,   32,   32, & ! O   
     &     62,    0,    0,    0,    0,   59,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,   37,  125,   34,   34, & ! 9   
     &      0,    0,  123,  125,   34,   34,    0,    0,   33,  125, & ! O   
     &     34,   34,   17,   13,  127,   13,  127,  172,   27,  172, & ! 1   
     &    127,  172,  172,  127,  143,  142,  167,   27,    0,    0, & ! 2   
     &      0,    0,  128,    0,    0,    0,    0,    0,   47,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   22, & ! O   
     &     22,   22,    0,   20,    0,    0,    0,    0,    0,    0, & ! 1   
     &     13,    0,    0,   14,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,   69,    0,   33,   34,   34,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,  104,  104,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  177,  177,  177, & ! 5   
     &     22,   22,   15,  177,  127,  177,  177,  177,  177,  127, & ! 6   
     &    127,  169,    0,    0,    0,    0,    0,    0,    0,   15, & ! 7   
     &     13,    0,   19,   22,   33,   33,   33,  205,   25,    0, & ! 8   
     &     17,    0,    0,    0,   12,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,   34,    0,    0,   17,   33,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,   37,  123,   33,   33,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,   20,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,   69,    0,    0,   33, & ! 8   
     &     34,   34,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &    104,  104,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 10 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,   13,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,   14,    0,    0,    0,    0,    0,    5, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    4,   14,   14, & ! 4   
     &      4,    0,    0,    0,   34,   13,   13,   13,   69,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   17, & ! 8   
     &     15,    0,   27,   28,   33,   33,   33,   33,   25,    0, & ! 9   
     &     13,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,   34,   22,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      5,    0,   13,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,   22,   22,    0,    0,    0,   36,  125,  125, & ! 7   
     &    125,   17,   22,    0,   22,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,  106,    0,  106,   34,    0, & ! 9   
     &      0,    0,    0,    0,    0,   87,    0,    0,    0,   22, & ! O   
     &     22,    0,   47,    0,   36,  125,  125,  125,   47,   12, & ! 1   
     &      0,   13,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,   34,   33,    0,    0,    0,    0, & ! 3   
     &      0,    0,   47,   47,    0,    0,    0,    0,   17,   17, & ! 4   
     &      0,   19,   28,  100,   35,   36,   36,   28,   13,    0, & ! 5   
     &     12,    0,    0,    0,    0,    0,    0,    0,    0,   34, & ! 6   
     &     33,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,  129,    0,    0, & ! 8   
     &    129,    0,    0,  129,    0,    0,  129,    0,    0,    0, & ! 9   
     &      0,    0,    0,   16,    0,    0,    0,    0,  135,  135, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,   33,   36,  125,  125, & ! 9   
     &      0,    0,   33,   36,  125,  125,    0,    0,   34,   36, & ! O   
     &    125,  125,   16,   15,   43,   15,   27,  127,   22,  127, & ! 1   
     &    142,  127,  127,  142,  169,  143,  126,   22,    0,    0, & ! 2   
     &      0,    0,   38,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,   29,    0,    0,    0,    0,    0,    0, & ! 1   
     &     17,    0,    0,   59,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,   34,   41,   41,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  172,  127,  127, & ! 5   
     &     13,   13,   17,  172,   27,  172,  127,  172,  127,  143, & ! 6   
     &    142,  167,    0,    0,    0,    0,    0,    0,    0,   17, & ! 7   
     &     15,    0,   27,   28,  125,   34,   34,    0,   13,    0, & ! 8   
     &     22,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,   34,   22,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,   33,   33,   34,   34,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,   29,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   34, & ! 8   
     &     41,   41,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 11 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,   13,    0,    0,    0,    0,    0,   12, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   14,    4,   13, & ! 4   
     &    111,    0,    0,    0,   69,   15,   33,   33,  111,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   19, & ! 8   
     &     17,    0,   45,   45,  125,   34,   34,   34,   13,    0, & ! 9   
     &     12,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,   34,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &     22,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,   33,  106,    0,    0,    0,  100,   35,   36, & ! 7   
     &     36,   33,   19,    0,   17,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,   87,    0,   87,  106,    0, & ! 9   
     &      0,    0,    0,    0,    0,   45,    0,    0,    0,   33, & ! O   
     &    111,    0,    0,    0,  100,   35,   36,   36,    0,   16, & ! 1   
     &      0,   12,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,  111,   34,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,   27,   22, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   17,   42,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &     34,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,   21,   21, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,   34,    0,   36,   36, & ! 9   
     &      0,    0,   34,    0,   36,   36,    0,    0,  125,    0, & ! O   
     &     36,   36,   18,   17,   27,   17,   22,  142,   13,  142, & ! 1   
     &     43,  142,  142,  143,  167,  169,  123,   13,    0,    0, & ! 2   
     &      0,    0,   47,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    5,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  127,   27,   43, & ! 5   
     &     15,   15,   16,  127,   22,  127,  142,  127,  142,  169, & ! 6   
     &    143,  126,    0,    0,    0,    0,    0,    0,    0,   19, & ! 7   
     &     17,    0,    0,    0,   36,  125,  125,    0,   42,    0, & ! 8   
     &     13,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,   34,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,   34,   34,  125,  125,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    5,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 12 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,   29,    0,    0,    0,    0,    0,   85, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   57,   57,   89, & ! 4   
     &      0,    0,    0,    0,  111,   22,   34,   34,   44,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   27, & ! 8   
     &     27,    0,    0,    0,   36,  125,  125,  125,   42,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,  106,   87,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,   28,   12,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,  111,   87,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,  111, & ! O   
     &     47,    0,    0,    0,   47,   47,   47,   47,    0,   22, & ! 1   
     &      0,   47,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,   47,  111,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,   19,   27, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,   22,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   47, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,  125,    0,    0,   47, & ! 9   
     &      0,    0,  125,    0,    0,   47,    0,    0,   36,    0, & ! O   
     &      0,   47,   50,   16,   22,   16,   13,  143,   15,  143, & ! 1   
     &     27,  143,  143,  169,  126,  167,  183,   15,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,   12,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  142,   22,   27, & ! 5   
     &     17,   17,   18,  142,   13,  142,   43,  142,  143,  167, & ! 6   
     &    169,  123,    0,    0,    0,    0,    0,    0,    0,   27, & ! 7   
     &     27,    0,    0,    0,  100,   35,   36,    0,   22,    0, & ! 8   
     &     12,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,  125,  125,   36,   36,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,   12,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 13 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,   32,    0,    0,    0,    0,    0,   14, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   58,   58,   90, & ! 4   
     &      0,    0,    0,    0,   17,   33,   41,   35,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   45, & ! 8   
     &     19,    0,    0,    0,  100,   35,   36,   36,   22,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,   87,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   47, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   47, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,   47,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   19, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,   36,    0,    0,    0, & ! 9   
     &      0,    0,   36,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,   95,   18,   13,   18,   15,  169,   17,  169, & ! 1   
     &     22,  169,  169,  167,  123,  126,  184,   17,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,   85,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  143,   13,   22, & ! 5   
     &     16,   16,   50,  143,   15,  143,   27,  143,  169,  126, & ! 6   
     &    167,  183,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &     19,    0,    0,    0,    0,    0,    0,    0,  205,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,   36,   36,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,   85,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 14 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,   69,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   60,   13,   22, & ! 4   
     &      0,    0,    0,    0,    0,   34,    0,   41,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &     22,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   15, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,   50,   15,   50,   17,  167,   16,  167, & ! 1   
     &     13,  167,  167,  126,   43,  123,   43,   16,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,   14,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  169,   15,   13, & ! 5   
     &     18,   18,   95,  169,   17,  169,   22,  169,  167,  123, & ! 6   
     &    126,  184,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &     22,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,   14,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 15 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,   64,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   38,   15,   58, & ! 4   
     &      0,    0,    0,    0,    0,   69,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &     28,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,   28, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,   95,   17,   95,   16,  126,   18,   43, & ! 1   
     &     15,  126,  126,  123,   41,  183,   41,   18,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  167,   17,   15, & ! 5   
     &     50,   50,    0,  167,   16,  167,   13,  167,  126,   43, & ! 6   
     &    123,   43,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &     28,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 16 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   13,   17,   40, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &     45,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,   96,   16,   96,   18,  123,   50,   27, & ! 1   
     &     17,  123,   43,   43,   27,  184,   27,   50,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  126,   16,   17, & ! 5   
     &     95,   95,    0,   43,   18,  126,   15,  126,  123,   41, & ! 6   
     &    183,   41,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 17 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   15,   19,   69, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,   18,    0,   50,   43,   95,   22, & ! 1   
     &     16,   43,   41,   41,   22,   43,   22,   95,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,  123,   18,   16, & ! 5   
     &     96,   96,    0,   27,   50,  123,   17,   43,   43,   27, & ! 6   
     &    184,   27,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 18 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   17,   22,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,   50,    0,   95,   41,   96,   13, & ! 1   
     &     18,   41,   27,   27,   13,   41,   13,   96,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   43,   50,   18, & ! 5   
     &      0,    0,    0,   22,   95,   43,   16,   41,   41,   22, & ! 6   
     &     43,   22,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 19 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   19,   28,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,   95,    0,   96,   27,    0,   15, & ! 1   
     &     50,   27,   22,   22,   15,   27,   15,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   41,   95,   50, & ! 5   
     &      0,    0,    0,   13,   96,   41,   18,   27,   27,   13, & ! 6   
     &     41,   13,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 20 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   37,  111,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,   96,    0,    0,   22,    0,   17, & ! 1   
     &     95,   22,   13,   13,   17,   22,   17,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   27,   96,   95, & ! 5   
     &      0,    0,    0,   15,    0,   27,   50,   22,   22,   15, & ! 6   
     &     27,   15,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 21 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   27,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   13,    0,   16, & ! 1   
     &     96,   13,   15,   15,   16,   13,   16,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   22,    0,   96, & ! 5   
     &      0,    0,    0,   17,    0,   22,   95,   13,   13,   17, & ! 6   
     &     22,   17,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 22 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   28,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   15,    0,   18, & ! 1   
     &      0,   15,   17,   17,   18,   15,   18,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   13,    0,    0, & ! 5   
     &      0,    0,    0,   16,    0,   13,   96,   15,   15,   16, & ! 6   
     &     13,   16,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 23 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,   69,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   17,    0,   50, & ! 1   
     &      0,   17,   16,   16,   50,   17,   50,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   15,    0,    0, & ! 5   
     &      0,    0,    0,   18,    0,   15,    0,   17,   17,   18, & ! 6   
     &     15,   18,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 24 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,  111,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   16,    0,   95, & ! 1   
     &      0,   16,   18,   18,   95,   16,   95,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   17,    0,    0, & ! 5   
     &      0,    0,    0,   50,    0,   17,    0,   16,   16,   50, & ! 6   
     &     17,   50,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 25 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   18,    0,   96, & ! 1   
     &      0,   18,   50,   50,   96,   18,   61,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   16,    0,    0, & ! 5   
     &      0,    0,    0,   95,    0,   16,    0,   18,   18,   95, & ! 6   
     &     16,   95,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 26 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   50,    0,    0, & ! 1   
     &      0,   50,   95,   95,    0,   50,   96,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   18,    0,    0, & ! 5   
     &      0,    0,    0,   96,    0,   18,    0,   50,   50,   96, & ! 6   
     &     18,   61,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 27 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   95,    0,    0, & ! 1   
     &      0,   95,   96,   96,    0,   95,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   50,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,   50,    0,   95,   95,    0, & ! 6   
     &     50,   96,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 28 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,   96,    0,    0, & ! 1   
     &      0,   96,    0,    0,    0,   61,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   95,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,   95,    0,   96,   96,    0, & ! 6   
     &     95,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 29 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,   96,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,   96,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,   96,    0,    0,    0,    0, & ! 6   
     &     61,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( IRR( IRXXN, 30 ), IRXXN = 1, NRXNS ) / & 
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &     96,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 2   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 3   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 4   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 5   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 6   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 7   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 8   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 9   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! O   
     &      0,    0,    0,    0,    0,    0,    0,    0,    0,    0, & ! 1   
     &      0/     !  2   

      DATA ( RTDAT( 1,IRXXN ), IRXXN = 1, NRXNS ) / & 
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 2   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.7000D-12, 1.0000D-14, & ! 4   
     &     3.0000D-12, 1.2000D-13, 6.1000D-34, 8.0000D-12, 3.3000D-11, & ! +   
     &     2.1500D-11, 1.6300D-10, 2.8000D-12, 4.8000D-11, 3.0000D-13, & ! 5   
     &     4.2000D-34, 1.8000D-12, 9.1000D-32, 7.1000D-31, 3.4400D-12, & ! +   
     &     6.0950D-14, 4.2500D-39, 3.0000D-12, 5.3000D-12, 3.4000D-31, & ! 6   
     &     1.8000D-30, 2.4000D-14, 2.0000D-11, 3.5000D-12, 1.7000D-11, & ! +   
     &     4.3500D-14, 8.5000D-13, 2.4000D-30, 5.8000D-27, 1.0000D-22, & ! 7   
     &     1.9000D-31, 2.1000D-27, 4.5000D-13, 2.9000D-31, 1.4400D-13, & ! +   
     &     2.4500D-12, 7.6600D-12, 7.6800D-12, 1.0100D-11, 1.0000D-28, & ! 8   
     &     5.7200D-12, 1.3300D-11, 5.5000D-30, 2.3300D-12, 1.8100D-12, & ! +   
     &     2.3300D-11, 7.0000D-12, 5.7000D-16, 1.5800D-14, 2.9500D-12, & ! 9   
     &     2.6900D-11, 4.5000D-13, 6.0000D-12, 4.6000D-12, 2.9700D-11, & ! +   
     &     3.1400D-11, 9.4200D-12, 2.4000D-11, 1.5000D-11, 4.0000D-05, & ! O   
     &     3.0000D-12, 1.2100D-11, 4.2000D-11, 5.2000D-12, 1.1000D-10, & ! +   
     &     5.5000D-12, 4.7000D-12, 4.9000D-12, 4.5600D-14, 1.5000D-12, & ! 1   
     &     2.8000D-12, 3.0000D-12, 8.0000D-12, 2.6000D-12, 5.7700D-12, & ! +   
     &     1.1000D-11, 9.2600D-13, 2.8000D-11, 2.8000D-11, 1.0000D-11, & ! 2   
     &     5.3200D-12, 6.7500D-12, 4.6500D-11, 2.0500D-10, 2.8500D-12, & ! +   
     &     3.0000D-12, 2.6000D-12, 1.4700D-11, 2.9000D-12, 3.4000D-12, & ! 3   
     &     3.4000D-12, 3.0000D-11, 4.5000D-13, 4.0000D-14, 2.9300D-12, & ! +   
     &     4.0000D-14, 4.0000D-14, 3.2000D-11, 5.3100D-12, 4.8000D-12, & ! 4   
     &     2.4000D-14, 9.1400D-15, 4.3300D-15, 4.4000D-15, 8.0500D-16, & ! +   
     &     2.8000D-15, 8.3000D-18, 1.6700D-16, 1.3600D-15, 8.5000D-16, & ! 5   
     &     1.6600D-18, 2.0000D-16, 2.0000D-16, 9.0000D-17, 2.8600D-13, & ! +   
     &     4.3920D-13, 1.7900D-13, 8.6400D-13, 1.1900D-12, 1.2200D-11, & ! 6   
     &     3.1500D-13, 2.0000D-12, 1.4000D-12, 3.7600D-12, 3.4000D-15, & ! +   
     &     5.0200D-13, 2.9000D-12, 3.7600D-12, 2.4000D-15, 3.7800D-12, & ! 7   
     &     1.0600D-12, 2.0100D-10, 2.2000D-14, 9.7000D-29, 2.9000D-02, & ! +   
     &     2.4000D-02, 9.7000D-29, 9.0000D-29, 9.7000D-29, 9.0000D-29, & ! 8   
     &     2.8000D-12, 1.6000D+16, 2.8000D-12, 2.6000D-12, 2.7000D-12, & ! +   
     &     2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, & ! 9   
     &     2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, & ! +   
     &     2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, & ! O   
     &     2.7000D-12, 2.7000D-12, 2.7000D-12, 8.1000D-12, 8.1000D-12, & ! +   
     &     2.9000D-12, 2.7000D-12, 2.7000D-12, 2.5400D-12, 2.5400D-12, & ! 1   
     &     2.5400D-12, 2.5400D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, & ! +   
     &     2.7000D-12, 2.7000D-12, 4.0000D-12, 2.7000D-12, 2.7000D-12, & ! 2   
     &     2.7000D-12, 2.0000D-11, 2.0000D-11, 2.0800D-12, 4.1000D-13, & ! +   
     &     7.5000D-13, 1.6600D-13, 1.6600D-13, 1.9000D-13, 1.6600D-13, & ! 3   
     &     1.6600D-13, 2.9100D-13, 2.9100D-13, 2.9100D-13, 2.9100D-13, & ! +   
     &     2.6000D-13, 2.7300D-13, 2.7100D-13, 2.7500D-13, 2.6000D-13, & ! 4   
     &     2.7300D-13, 2.7100D-13, 2.7500D-13, 2.7100D-13, 2.7300D-13, & ! +   
     &     4.3000D-13, 4.3000D-13, 1.1500D-13, 1.1500D-13, 1.1500D-13, & ! 5   
     &     1.8200D-13, 1.8200D-13, 2.9100D-13, 2.9100D-13, 2.3800D-13, & ! +   
     &     2.1100D-13, 2.8600D-13, 3.7500D-13, 1.0000D-11, 3.7500D-13, & ! 6   
     &     1.1500D-13, 1.6600D-13, 1.6600D-13, 3.7500D-13, 1.6600D-13, & ! +   
     &     9.5000D-14, 1.1800D-13, 9.4600D-14, 1.0000D-13, 1.7100D-13, & ! 7   
     &     1.4600D-13, 9.1800D-14, 3.5600D-14, 3.5600D-14, 3.5600D-14, & ! +   
     &     3.5600D-14, 3.4000D-14, 2.0000D-12, 1.0000D-10, 2.0000D-12, & ! 8   
     &     1.0000D-10, 2.0000D-12, 1.0000D-10, 2.0000D-12, 1.0000D-10, & ! +   
     &     2.0000D-11, 2.0000D-11, 7.5000D-13, 6.9100D-13, 6.9100D-13, & ! 9   
     &     3.4000D-14, 3.4000D-14, 8.3700D-14, 3.4000D-14, 3.5600D-14, & ! +   
     &     3.5600D-14, 3.5600D-14, 3.5600D-14, 7.5000D-13, 1.6000D-13, & ! O   
     &     9.6800D-14, 3.5600D-14, 5.9900D-15, 1.0300D-12, 6.9000D-13, & ! +   
     &     5.5900D-13, 9.4800D-13, 8.1100D-13, 5.0900D-13, 7.4000D-13, & ! 1   
     &     7.4000D-13, 7.4000D-13, 7.4000D-13, 8.4000D-14, 2.0000D-12, & ! +   
     &     1.0000D-10, 2.0000D-12, 1.0000D-10, 2.0000D-12, 1.0000D-10, & ! 2   
     &     2.0000D-12, 1.0000D-10, 2.5000D-12, 2.5000D-12, 7.5100D-13, & ! +   
     &     7.5100D-13, 7.5100D-13, 8.4000D-14, 8.4000D-14, 1.6800D-12, & ! 3   
     &     1.6800D-12, 7.4000D-13, 7.4000D-13, 7.4000D-13, 7.4000D-13, & ! +   
     &     7.5100D-13, 8.8500D-13, 5.3700D-13, 7.4000D-13, 3.4000D-14, & ! 4   
     &     2.5000D-12, 1.2000D-12, 1.2000D-12, 1.2000D-12, 1.2000D-12, & ! +   
     &     1.2000D-12, 1.2000D-12, 1.2000D-12, 2.3000D-12, 2.3000D-12, & ! 5   
     &     2.3000D-12, 2.3000D-12, 1.2000D-12, 2.3000D-12, 2.3000D-12, & ! +   
     &     2.3000D-12, 2.3000D-12, 4.0000D-12, 4.0000D-12, 1.2000D-12, & ! 6   
     &     1.2000D-12, 1.2000D-12, 1.2000D-12, 1.2000D-12, 2.5000D-12, & ! +   
     &     2.5000D-12, 2.5000D-12, 2.5000D-12, 1.2000D-12, 1.2000D-12, & ! 7   
     &     1.2000D-12, 1.2000D-12, 1.2000D-12, 1.2000D-12, 7.0000D-14, & ! +   
     &     4.2500D-14, 2.9600D-14, 1.2000D-12, 2.5000D-12, 7.1300D-17, & ! 8   
     &     1.0000D-10, 1.0000D-10, 1.0000D-10, 1.0000D-10, 1.0000D-10, & ! +   
     &     1.0000D-10, 1.0000D-10, 1.0000D-10, 1.0000D-10, 1.0000D-10, & ! 9   
     &     1.0000D-10, 1.0000D-10, 5.7800D-11, 6.8900D-12, 6.5500D-14, & ! +   
     &     8.0000D-12, 2.9000D-19, 3.4000D-15, 1.0000D+00, 1.4800D-11, & ! O   
     &     2.7000D-12, 2.3000D-12, 1.9300D-13, 2.3900D-12, 1.3700D-11, & ! +   
     &     1.3400D-14, 1.0000D-13, 5.0100D-11, 2.7000D-12, 3.7500D-13, & ! 1   
     &     4.4000D-11, 3.4300D-17, 8.9900D-12, 1.2000D-11, 1.9000D-11, & ! +   
     &     2.8400D-13, 2.7000D-12, 2.3000D-12, 1.2000D-14, 1.9700D-10, & ! 2   
     &     2.8400D-13, 2.3000D-12, 2.7000D-12, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 2.0000D+00, 5.0000D-01, 1.5300D-11, 1.6800D-11, & ! +   
     &     2.2400D-11, 2.6700D-11, 3.0900D-11, 3.3800D-11, 1.1000D-11, & ! 4   
     &     2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, & ! +   
     &     2.7000D-12, 2.7000D-12, 2.3000D-12, 2.3000D-12, 2.3000D-12, & ! 5   
     &     2.3000D-12, 2.3000D-12, 2.3000D-12, 2.3000D-12, 2.6800D-13, & ! +   
     &     2.7100D-13, 2.7800D-13, 2.8000D-13, 2.8100D-13, 2.8100D-13, & ! 6   
     &     2.6600D-13, 1.8800D-01, 1.8800D-01, 1.8800D-01, 1.8800D-01, & ! +   
     &     1.8800D-01, 1.8800D-01, 1.8800D-01, 2.7000D-12, 2.7000D-12, & ! 7   
     &     2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, 2.7000D-12, & ! +   
     &     2.3000D-12, 2.3000D-12, 2.3000D-12, 2.3000D-12, 2.3000D-12, & ! 8   
     &     2.3000D-12, 2.3000D-12, 2.7100D-13, 2.7300D-13, 2.7800D-13, & ! +   
     &     2.8000D-13, 2.8100D-13, 2.8200D-13, 2.6600D-13, 1.8100D-11, & ! 9   
     &     2.9100D-13, 2.7000D-12, 2.3000D-12, 3.5600D-14, 7.4000D-13, & ! +   
     &     1.8100D-11, 2.9100D-13, 2.7000D-12, 2.3000D-12, 3.5600D-14, & ! O   
     &     7.4000D-13, 2.3100D-11, 2.9100D-13, 2.7000D-12, 2.3000D-12, & ! +   
     &     3.5600D-14, 7.4000D-13, 5.9000D-11, 6.0700D-11, 5.5400D-11, & ! 1   
     &     5.6300D-11, 5.4600D-11, 4.5000D-11, 5.1700D-11, 4.7300D-11, & ! +   
     &     4.6000D-11, 3.8000D-11, 3.9300D-11, 3.5200D-11, 3.1200D-11, & ! 2   
     &     2.4000D-11, 2.0500D-11, 4.7400D-11, 2.4500D-12, 9.2600D-05, & ! +   
     &     9.2600D-05, 1.0200D-11, 1.4000D-17, 1.5100D-13, 2.9100D-12, & ! 3   
     &     2.7000D-12, 2.3000D-12, 2.5000D-13, 2.5000D-13, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D-02, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.3000D-11, 3.0500D-11, 3.6000D-11, & ! 5   
     &     1.4000D-11, 1.1000D-11, 1.2000D-11, 2.3000D-10, 3.4000D-11, & ! +   
     &     7.6000D-11, 3.4000D-12, 7.7000D-32, 1.8000D-31, 2.4000D-11, & ! 6   
     &     5.8000D-11, 6.2000D-12, 2.6000D-12, 6.4000D-12, 1.8000D-31, & ! +   
     &     1.0000D-12, 3.5000D-13, 1.9000D-32, 2.1600D-27, 1.1000D-13, & ! 7   
     &     2.6000D-12, 6.0000D-13, 7.4000D-12, 1.8000D-12, 3.0000D-12, & ! +   
     &     2.4000D-12, 1.2000D-12, 1.4000D-12, 6.0000D-13, 7.1000D-12, & ! 8   
     &     7.2000D-11, 1.4000D-10, 2.7000D-10, 5.2700D-10, 8.1000D-11, & ! +   
     &     8.0000D-11, 2.8000D-11, 2.6000D-11, 5.5000D-11, 6.0000D-11, & ! 9   
     &     1.7000D-10, 2.0000D-13, 2.7000D-14, 2.3000D-10, 2.3000D-10, & ! +   
     &     4.1000D-11, 1.1600D-10, 5.4000D-11, 2.6000D-10, 2.0000D-10, & ! O   
     &     3.6000D-11, 5.0000D-11, 4.6000D-10, 5.9000D-11, 1.1000D-10, & ! +   
     &     1.1000D-10, 1.1000D-10, 1.3600D-09, 1.5600D-09, 5.1000D-10, & ! 1   
     &     2.2000D-13, 1.0000D-14, 1.1300D-12, 9.1000D-10, 2.2000D-11, & ! +   
     &     4.8000D-11, 9.1000D-11, 1.3000D-10, 1.8000D-09, 2.1000D-10, & ! 2   
     &     9.5000D-11, 3.0000D-10, 3.5000D-10, 4.6000D-10, 2.4000D-10, & ! +   
     &     7.4000D-10, 1.0600D-11, 2.6000D-10, 3.7000D-10, 3.7000D-10, & ! 3   
     &     7.4000D-11, 5.3000D-10, 2.2300D-10, 3.3000D-10, 1.5000D-16, & ! +   
     &     6.0000D-11, 1.4500D-10, 9.9000D-11, 2.0000D-10, 3.5000D-11, & ! 4   
     &     1.6000D-09, 2.8000D-15, 2.0000D-10, 1.8200D-09, 1.6700D-09, & ! +   
     &     1.4500D-09, 1.2240D-09, 9.3000D-10, 8.5100D-10, 3.4300D-11, & ! 5   
     &     1.3700D-10, 1.3700D-10, 3.0900D-10, 3.7700D-10, 3.8300D-10, & ! +   
     &     4.2000D-10, 3.8900D-10, 4.0800D-10, 3.2500D-10, 3.5600D-10, & ! 6   
     &     2.5900D-10, 3.1600D-10, 2.6800D-10, 2.3900D-10, 2.1100D-10, & ! +   
     &     1.6000D-10, 1.3500D-10, 3.8800D-11, 3.5600D-11, 8.1000D-11, & ! 7   
     &     1.2900D-10, 3.6000D-10, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! +   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! 8   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! +   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! 9   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! +   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! O   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! +   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! 1   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! +   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! 2   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! +   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, & ! 3   
     &     3.2500D-12, 3.2500D-12, 3.2500D-12, 3.2500D-12, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     2.1000D-11, 1.6000D-11, 4.8000D-12, 4.3000D-31, 1.6000D-11, & ! 5   
     &     4.9000D-11, 5.5000D-12, 1.7000D-11, 4.5000D-12, 8.8000D-12, & ! +   
     &     5.5000D-31, 1.0000D-12, 2.4000D-12, 2.8000D-14, 1.5000D-15, & ! 6   
     &     1.7000D-11, 1.8000D-11, 4.3000D-11, 8.0000D-15, 2.0000D-14, & ! +   
     &     3.3000D-14, 5.7000D-11, 5.7000D-11, 8.0000D-15, 1.7000D-11, & ! 7   
     &     1.7000D-11, 2.0000D-14, 6.4000D-16, 7.3600D-15, 8.0000D-15, & ! +   
     &     1.3000D-13, 3.7000D-12, 1.0000D-11, 8.1000D-11, 2.4000D-11, & ! 8   
     &     7.4000D-11, 6.4000D-15, 2.0000D-11, 2.6000D-11, 3.7000D-11, & ! +   
     &     3.7000D-11, 7.4000D-12, 5.3000D-11, 4.3000D-12, 2.5000D-11, & ! 9   
     &     6.3000D-11, 5.0000D-16, 1.4500D-14, 7.2500D-14, 1.0000D-13, & ! +   
     &     7.2500D-14, 7.2500D-14, 7.2500D-14, 3.0000D-16, 5.0100D-12, & ! O   
     &     2.0000D-12, 8.5000D-12, 9.0000D-12, 3.1400D-14, 5.6000D-12, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.8000D-10, 1.5000D-12, & ! 2   
     &     2.1000D-11, 1.5000D-11, 1.8000D-32, 3.0000D-31, 9.1000D-11, & ! +   
     &     9.9400D+17, 2.1000D+15, 8.4000D-11, 2.9000D-11, 5.0000D-12, & ! 3   
     &     3.0000D-11, 1.4000D-11, 8.6000D-12, 1.1000D-12, 7.7000D-31, & ! +   
     &     1.5000D-11, 2.5000D+14, 2.7000D-12, 1.0000D+12, 1.5000D-10, & ! 4   
     &     1.5000D-10, 3.8000D-02, 1.3000D-12, 1.2500D-10, 4.3000D-12, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 2.6000D-11, 5.9000D-12, 1.2000D-11, & ! +   
     &     2.7000D-11, 9.5000D-13, 2.9000D-12, 2.7000D-12, 4.4000D-12, & ! 6   
     &     6.8000D-12, 4.3500D-12, 6.8000D-12, 2.1000D-13, 4.8500D-12, & ! +   
     &     2.0300D-11, 9.0000D-13, 9.4000D-13, 1.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 9   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00/           !        2   

      DATA ( RTDAT( 2,IRXXN ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00,-2.4000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.6000D+02, & ! 5   
     &     2.6600D+03, 0.0000D+00,-1.5000D+00,-2.6000D+00, 0.0000D+00, & ! +   
     &     2.7000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00,-1.6000D+00, & ! 6   
     &    -3.0000D+00, 4.6000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00,-3.0000D+00, 1.0840D+04, 0.0000D+00, & ! 7   
     &    -3.4000D+00, 1.0900D+04, 0.0000D+00,-4.1000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,-4.5000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.6500D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     4.6000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-5.6000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00,-5.6000D+00, 1.4000D+04,-5.6000D+00, 1.4000D+04, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00,-1.8000D+00,-2.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,-3.4000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00,-3.6000D+00, 8.5370D+03, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-2.4000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &    -3.1000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00,-1.0000D+00,-1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00,-3.5000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00,-2.6600D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( RTDAT( 3,IRXXN ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-9.4000D+02,-4.9000D+02, & ! 4   
     &    -1.5000D+03,-2.4500D+03, 0.0000D+00,-2.0600D+03, 5.5000D+01, & ! +   
     &     1.1000D+02, 6.0000D+01,-1.8000D+03, 2.5000D+02, 2.1000D-33, & ! 5   
     &     2.9400D-54, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.6000D+02, & ! +   
     &     6.8570D-34, 6.6350D+02, 2.5000D+02, 2.0000D+02, 0.0000D+00, & ! 6   
     &     0.0000D+00, 2.7000D-17, 0.0000D+00, 0.0000D+00, 1.2500D+02, & ! +   
     &    -1.3350D+03,-2.4500D+03, 0.0000D+00, 7.3000D+01, 0.0000D+00, & ! 7   
     &     0.0000D+00, 7.6000D+01, 6.1000D+02, 0.0000D+00, 2.7400D-33, & ! +   
     &    -1.7750D+03,-1.0200D+03,-3.7000D+02,-2.4500D+02, 0.0000D+00, & ! 8   
     &     5.0000D+02, 5.0000D+02, 0.0000D+00,-1.9300D+02, 3.5400D+02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-2.0000D+03,-4.5000D+02, & ! 9   
     &     3.9000D+02, 1.3000D+03, 3.5000D+02, 2.0000D+02, 3.9000D+02, & ! +   
     &     5.8000D+02, 5.8000D+02, 3.9000D+02, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 4.4000D+02, 4.0100D+02, 6.0000D+02, 0.0000D+00, & ! +   
     &     1.2500D+02, 3.4500D+02, 4.0500D+02,-4.2700D+02,-9.0000D+01, & ! 1   
     &     1.0000D+01, 0.0000D+00, 3.8000D+02, 6.1000D+02, 5.3300D+02, & ! +   
     &     0.0000D+00, 8.3000D+02, 1.7500D+02, 1.7500D+02, 0.0000D+00, & ! 2   
     &     2.4300D+02, 4.0500D+02, 0.0000D+00, 0.0000D+00,-3.4500D+02, & ! +   
     &     2.0000D+01, 2.0000D+02, 0.0000D+00, 1.9000D+02, 1.9000D+02, & ! 3   
     &     1.9000D+02, 0.0000D+00, 0.0000D+00, 8.5000D+02, 1.9000D+02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-2.6000D+02, 0.0000D+00, & ! 4   
     &     2.7000D-17,-2.5800D+03,-1.8000D+03,-8.4500D+02,-6.4000D+02, & ! +   
     &    -7.7000D+02, 0.0000D+00, 0.0000D+00,-2.1120D+03,-1.5200D+03, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &    -2.2820D+03,-4.5000D+02, 4.5000D+02, 4.9000D+02, 0.0000D+00, & ! 6   
     &    -4.4800D+02,-2.4400D+03,-1.9000D+03,-1.9000D+03, 0.0000D+00, & ! +   
     &    -1.0760D+03,-1.9000D+03,-1.9000D+03, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00,-5.0000D+02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.8200D+02, 0.0000D+00, 1.8400D+02, & ! 8   
     &     1.8100D+02,-1.3486D+04, 3.0000D+02, 3.6500D+02, 3.6000D+02, & ! +   
     &     3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, & ! 9   
     &     3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, & ! +   
     &     3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, & ! O   
     &     3.6000D+02, 3.6000D+02, 3.6000D+02, 2.7000D+02, 2.7000D+02, & ! +   
     &     3.0000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, & ! 1   
     &     3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, & ! +   
     &     3.6000D+02, 3.6000D+02, 0.0000D+00, 3.6000D+02, 3.6000D+02, & ! 2   
     &     3.6000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 7.5000D+02, & ! +   
     &     7.0000D+02, 1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! 3   
     &     1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! +   
     &     1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! 4   
     &     1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! +   
     &     1.0400D+03, 1.0400D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! 5   
     &     1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! +   
     &     1.3000D+03, 0.0000D+00, 9.8000D+02, 0.0000D+00, 9.8000D+02, & ! 6   
     &     1.3000D+03, 1.3000D+03, 1.3000D+03, 9.8000D+02, 1.3000D+03, & ! +   
     &     3.9000D+02, 1.5800D+02, 4.3100D+02, 4.6700D+02, 7.0800D+02, & ! 7   
     &     7.0800D+02, 7.0800D+02, 7.0800D+02, 7.0800D+02, 7.0800D+02, & ! +   
     &     7.0800D+02, 2.2100D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     5.0000D+02, 5.0000D+02, 5.0000D+02, 5.0800D+02, 5.0800D+02, & ! 9   
     &     2.2100D+02, 2.2100D+02, 0.0000D+00, 2.2100D+02, 7.0800D+02, & ! +   
     &     7.0800D+02, 7.0800D+02, 7.0800D+02, 5.0000D+02, 7.0800D+02, & ! O   
     &     7.0800D+02, 0.0000D+00, 1.5100D+03, 2.1100D+02, 4.6000D+02, & ! +   
     &     5.2200D+02, 7.6500D+02, 7.6500D+02, 7.6500D+02, 7.6500D+02, & ! 1   
     &     7.6500D+02, 7.6500D+02, 7.6500D+02, 2.2100D+02, 5.0000D+02, & ! +   
     &     0.0000D+00, 5.0000D+02, 0.0000D+00, 5.0000D+02, 0.0000D+00, & ! 2   
     &     5.0000D+02, 0.0000D+00, 5.0000D+02, 5.0000D+02, 5.6500D+02, & ! +   
     &     5.6500D+02, 5.6500D+02, 2.2100D+02, 2.2100D+02, 5.0000D+02, & ! 3   
     &     5.0000D+02, 7.6500D+02, 7.6500D+02, 7.0800D+02, 7.0800D+02, & ! +   
     &     5.6500D+02, 7.6500D+02, 7.6500D+02, 7.0800D+02, 1.5600D+03, & ! 4   
     &     5.0000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+03, & ! +   
     &     1.0000D+03, 1.0000D+03, 0.0000D+00, 5.0000D+02, 2.9500D+03, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00,-4.0000D+02, 0.0000D+00, 0.0000D+00, & ! +   
     &     3.8000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4800D+02, & ! O   
     &     3.6000D+02, 0.0000D+00, 1.3000D+03, 0.0000D+00, 0.0000D+00, & ! +   
     &    -2.2830D+03, 0.0000D+00, 0.0000D+00, 3.6000D+02, 9.8000D+02, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.3000D+03, 3.6000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     1.3000D+03, 0.0000D+00, 3.6000D+02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, & ! +   
     &     3.6000D+02, 3.6000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3000D+03, & ! +   
     &     1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! 6   
     &     1.3000D+03, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.6000D+02, 3.6000D+02, & ! 7   
     &     3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, 3.6000D+02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 1.3000D+03, 1.3000D+03, 1.3000D+03, & ! +   
     &     1.3000D+03, 1.3000D+03, 1.3000D+03, 1.3000D+03, 0.0000D+00, & ! 9   
     &     1.3000D+03, 3.6000D+02, 0.0000D+00, 7.0800D+02, 7.6500D+02, & ! +   
     &     0.0000D+00, 1.3000D+03, 3.6000D+02, 0.0000D+00, 7.0800D+02, & ! O   
     &     7.6500D+02, 0.0000D+00, 1.3000D+03, 3.6000D+02, 0.0000D+00, & ! +   
     &     7.0800D+02, 7.6500D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-1.7750D+03, 0.0000D+00, & ! +   
     &     0.0000D+00, 5.3200D+02, 0.0000D+00, 0.0000D+00, 1.3000D+03, & ! 3   
     &     3.6000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00,-2.0000D+02,-2.2700D+03,-3.7500D+02, & ! 5   
     &     2.7000D+02,-9.8000D+02, 0.0000D+00, 0.0000D+00, 1.6000D+02, & ! +   
     &     6.5000D+01,-1.3000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     1.0000D+02, 1.4500D+02, 2.9000D+02, 2.9000D+02, 0.0000D+00, & ! +   
     &    -1.5900D+03,-1.3700D+03, 0.0000D+00, 5.7300D+02, 3.5000D+02, & ! 7   
     &    -1.1000D+03, 2.3000D+02, 2.7000D+02,-2.5000D+02,-5.0000D+02, & ! +   
     &    -1.2500D+03,-3.3000D+02, 6.0000D+02, 6.7000D+02,-1.2700D+03, & ! 8   
     &    -7.0000D+01, 0.0000D+00, 0.0000D+00, 0.0000D+00,-3.0000D+01, & ! +   
     &     0.0000D+00, 4.5300D+02,-7.4500D+02, 0.0000D+00, 1.5500D+02, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 4.8000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-1.4190D+03,-3.0000D+01, & ! 7   
     &     0.0000D+00, 0.0000D+00,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! +   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! 8   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! +   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! 9   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! +   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! O   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! +   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! 1   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! +   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! 2   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! +   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, & ! 3   
     &    -1.1400D+02,-1.1400D+02,-1.1400D+02,-1.1400D+02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.4000D+02,-7.8000D+02,-3.1000D+02, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 2.0000D+02, 2.5000D+02, 4.6000D+02, 2.6000D+02, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.0000D+01, 8.6000D+02, 0.0000D+00, & ! 6   
     &    -8.0000D+02,-4.6000D+02,-5.2600D+02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00,-5.3600D+02,-5.3600D+02, 0.0000D+00,-8.0000D+02, & ! 7   
     &    -8.0000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 4.4000D+02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &    -8.4000D+02, 0.0000D+00,-3.6000D+02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &    -8.3000D+02,-1.0900D+03, 0.0000D+00, 0.0000D+00,-1.4600D+02, & ! +   
     &    -1.1859D+04,-1.3670D+04,-2.6200D+03,-2.6000D+03, 0.0000D+00, & ! 3   
     &     0.0000D+00, 5.4000D+02, 2.3000D+02, 5.4200D+02, 0.0000D+00, & ! +   
     &     5.0000D+02,-9.7700D+03, 0.0000D+00,-9.7700D+03, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00,-1.8300D+03, 0.0000D+00,-1.1200D+03, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00,-1.3000D+03,-1.7000D+02, 0.0000D+00, & ! +   
     &     0.0000D+00, 5.5000D+02, 2.2000D+02, 2.8000D+02, 7.6000D+02, & ! 6   
     &    -8.3000D+02,-8.5000D+02,-8.7000D+02,-8.8000D+02,-8.5000D+02, & ! +   
     &    -1.1100D+03,-4.2000D+02,-5.1000D+02, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   
      INTEGER            :: IRRFALL( NFALLOFF )

      DATA ( IRRFALL( IRXXN ), IRXXN = 1, NFALLOFF ) / & 
     &     55,   56,   58,   59,   61,   65,   66,   67,   73,   76, & 
     &     79,   80,   85,   88,  146,  179,  182,  184,  564,  570, & 
     &    573,  754,  761,  828,  829,  840/

      DATA ( RFDAT( 1,IRXXN ), IRXXN = 1, NFALLOFF ) / & 
     &     9.2000D+02, 3.1200D+03, 3.0000D-11, 3.6000D-11, 2.7000D+02, & 
     &     2.3000D-11, 2.8000D-11, 2.1990D+03, 1.6000D-12, 4.0000D-12, & 
     &     1.7000D-12, 0.0000D+00, 8.8000D-12, 8.3000D-13, 2.1990D+03, & 
     &     9.3000D-12, 9.3000D-12, 9.3000D-12, 1.0000D-10, 1.5000D-11, & 
     &     3.7000D-12, 2.7000D-11, 6.6000D-12, 1.7000D-11, 6.6000D-11, & 
     &     7.7000D-12/

      DATA ( RFDAT( 2,IRXXN ), IRXXN = 1, NFALLOFF ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00,-1.0000D-01,-1.0000D+00, & 
     &    -2.0000D-01, 0.0000D+00, 6.5000D-34, 1.0000D-01,-3.0000D-01, & 
     &     2.0000D-01, 0.0000D+00,-8.5000D-01, 2.0000D+00, 6.5000D-34, & 
     &    -1.5000D+00,-1.5000D+00,-1.5000D+00,-1.0000D+00,-1.9000D+00, & 
     &    -1.6000D+00, 0.0000D+00,-2.9000D+00, 0.0000D+00, 0.0000D+00, & 
     &    -1.5000D+00/

      DATA ( RFDAT( 3,IRXXN ), IRXXN = 1, NFALLOFF ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 1.3350D+03, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3350D+03, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & 
     &     0.0000D+00/

      DATA ( RFDAT( 4,IRXXN ), IRXXN = 1, NFALLOFF ) / & 
     &     0.0000D+00, 0.0000D+00, 6.0000D-01, 6.0000D-01,-5.9680D-14, & 
     &     6.0000D-01, 6.0000D-01, 0.0000D+00, 6.0000D-01, 6.0000D-01, & 
     &     6.0000D-01, 0.0000D+00, 6.0000D-01, 6.0000D-01, 0.0000D+00, & 
     &     1.4200D+00, 6.0000D-01, 6.0000D-01, 6.0000D-01, 6.0000D-01, & 
     &     6.0000D-01, 6.0000D-01, 6.0000D-01, 6.0000D-01, 6.0000D-01, & 
     &     6.0000D-01/

      DATA ( RFDAT( 5,IRXXN ), IRXXN = 1, NFALLOFF ) / & 
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 2.7000D+02, & 
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & 
     &     1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, & 
     &     3.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & 
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & 
     &     1.0000D+00/

      REAL( 8 )               :: SC( NRXNS,MXPRD )

      DATA ( SC( IRXXN,  1 ), IRXXN = 1, NRXNS ) / & 
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 2.0000D-01, 1.0000D+00, & ! +   
     &     2.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 2.0000D+00, & ! 1   
     &     1.2200D+00, 1.0000D+00, 1.0000D+00, 1.0000D-01, 1.5000D+00, & ! +   
     &     1.0000D+00, 3.4000D-01, 3.0000D-01, 2.0000D+00, 1.0000D+00, & ! 2   
     &     2.0000D+00, 1.0000D+00, 1.5000D+00, 1.5000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 7.0000D-01, 2.0000D+00, & ! +   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, 2.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 6.5000D-01, 4.7000D-01, 8.2000D-01, & ! +   
     &     8.3000D-01, 8.2000D-01, 1.0000D+00, 2.5000D-01, 4.0000D-01, & ! 9   
     &     1.0000D+00, 1.0000D+00, 1.3000D-01, 1.0000D+00, 4.0000D-02, & ! +   
     &     4.5000D-01, 1.5000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     5.7000D-01, 9.7500D-01, 9.4500D-01, 2.3000D-01, 8.3000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 5.7000D-01, 1.0000D+00, 3.1300D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 5.2000D-01, 5.2000D-01, 5.6000D-01, & ! 2   
     &     1.0000D+00, 1.5200D-01, 2.0000D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 3.5000D-01, 1.0000D-02, & ! 3   
     &     1.0000D-02, 1.0000D+00, 1.0000D+00, 6.4000D-01, 3.5000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 3.3000D-01, & ! 4   
     &     1.0000D+00, 8.0000D-02, 2.2000D-01, 4.6000D-01, 6.0550D-01, & ! +   
     &     6.6000D-01, 9.0000D-02, 3.3000D-01, 1.9000D-01, 1.6000D-01, & ! 5   
     &     1.0000D-01, 5.0000D-02, 5.0000D-02, 5.0000D-02, 1.0000D+00, & ! +   
     &     8.0000D-01, 4.3000D-01, 1.1000D-01, 9.7500D-01, 9.4500D-01, & ! 6   
     &     3.3000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 6.8000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.5200D-01, & ! 7   
     &     2.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 6.6000D-01, & ! +   
     &     2.0000D-01, 1.0000D+00, 7.8000D-01, 8.3000D-01, 0.0000D+00, & ! 9   
     &     2.0000D-04, 1.0000D-04, 2.0000D-04, 2.3700D-01, 8.2000D-01, & ! +   
     &     1.6955D+00, 8.2000D-01, 2.3000D-01, 7.7000D-01, 5.7000D-01, & ! O   
     &     7.7000D-01, 3.6000D-01, 6.4000D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 6.7000D-01, 7.7000D-01, 6.5000D-01, 1.0000D+00, & ! 1   
     &     3.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 2.0000D+00, 2.0000D+00, & ! 2   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 6.0210D-01, 7.1950D-01, 4.8200D-02, 8.5400D-02, & ! +   
     &     3.9000D-02, 1.0000D+00, 4.6100D-01, 1.0000D+00, 9.0000D-01, & ! 4   
     &     1.0000D+00, 7.7000D-01, 1.0000D+00, 7.5000D-01, 9.0000D-01, & ! +   
     &     4.4000D-01, 4.4000D-01, 1.5000D-01, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 4.1000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     7.4000D-01, 1.0000D+00, 8.9400D-01, 8.4200D-01, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 6.8000D-01, 6.8000D-01, 6.8000D-01, & ! +   
     &     6.8000D-01, 1.0000D+00, 1.0295D+00, 1.0000D+00, 3.4100D-01, & ! 8   
     &     7.5000D-01, 2.5000D-01, 1.0000D+00, 6.9000D-01, 7.5000D-01, & ! +   
     &     9.0000D-01, 9.0000D-01, 5.0000D-01, 8.3400D-01, 1.0000D+00, & ! 9   
     &     5.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, 2.0000D+00, & ! O   
     &     5.0000D-01, 1.0000D+00, 1.0000D+00, 5.0000D-01, 3.9400D-01, & ! +   
     &     3.4200D-01, 5.0000D-01, 5.0000D-01, 5.0000D-01, 7.0000D-01, & ! 1   
     &     7.0000D-01, 7.0000D-01, 7.0000D-01, 5.0000D-01, 2.5350D-01, & ! +   
     &     5.0000D-01, 8.2500D-01, 5.0000D-01, 1.0000D+00, 5.0000D-01, & ! 2   
     &     1.0000D+00, 5.0000D-01, 2.0000D+00, 1.0000D+00, 5.0000D-01, & ! +   
     &     3.3000D-01, 5.0000D-01, 6.3500D-01, 5.0000D-01, 5.0000D-01, & ! 3   
     &     5.0000D-01, 1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 5.0000D-01, 1.0000D+00, 1.0000D+00, & ! 4   
     &     2.0000D+00, 1.0000D+00, 1.0000D+00, 2.5400D-01, 4.8800D-01, & ! +   
     &     1.0000D+00, 4.7000D-01, 8.6000D-01, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.8250D+00, 1.4600D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     6.7000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 3.0000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 2.0000D+00, 1.0000D+00, & ! +   
     &     5.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, & ! 8   
     &     9.6000D-01, 9.6000D-01, 9.6000D-01, 9.6000D-01, 9.6000D-01, & ! +   
     &     9.6000D-01, 9.6000D-01, 9.6000D-01, 9.6000D-01, 9.6000D-01, & ! 9   
     &     9.6000D-01, 9.6000D-01, 1.0000D+00, 8.5200D-01, 1.0000D+00, & ! +   
     &     5.7000D-01, 8.4000D-01, 6.8000D-01, 1.0000D+00, 6.6700D-01, & ! O   
     &     9.6800D-01, 1.0000D+00, 1.0000D+00, 3.2000D-01, 7.0000D-01, & ! +   
     &     6.2000D-01, 9.0000D-01, 4.9000D-01, 8.0000D-02, 6.0000D-01, & ! 1   
     &     6.5000D-01, 2.0000D-02, 1.0000D+00, 6.1300D-01, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 9.8200D-01, 1.0000D+00, & ! 2   
     &     1.0000D+00, 1.0000D+00, 2.4700D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     5.0000D-01, 0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, & ! +   
     &     7.2000D-01, 7.4000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.4000D-01, 1.4000D-01, & ! 7   
     &     1.4000D-01, 1.4000D-01, 1.4000D-01, 1.0000D+00, 1.2000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 8.4000D-01, & ! 9   
     &     5.9500D-02, 1.0000D-04, 1.0000D+00, 6.8000D-01, 7.0000D-01, & ! +   
     &     8.4000D-01, 5.9500D-02, 1.0000D-04, 1.0000D+00, 6.8000D-01, & ! O   
     &     7.0000D-01, 8.4000D-01, 5.9500D-02, 5.9500D-02, 1.0000D+00, & ! +   
     &     6.8000D-01, 7.0000D-01, 8.5400D-02, 8.0830D-01, 1.1560D-01, & ! 1   
     &     2.0360D-01, 3.6090D-01, 7.4000D-03, 3.6950D-01, 7.3800D-02, & ! +   
     &     1.9620D-01, 2.3000D-03, 4.7900D-02, 3.2400D-02, 1.1700D-02, & ! 2   
     &     1.0300D-02, 6.1000D-03, 1.6280D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 6.2000D-01, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.3700D+00, 1.0000D+00, 6.7000D-01, & ! +   
     &     6.7000D-01, 8.0000D-01, 2.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, 2.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 9   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 2   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 6.6000D-01, & ! +   
     &     2.0000D-01, 1.0000D+00, 7.8000D-01, 8.3000D-01, 1.2000D-03, & ! 8   
     &     1.3000D-03, 1.3000D-03, 4.0000D-01, 7.6300D-01, 8.2000D-01, & ! +   
     &     1.0000D+00, 8.2000D-01, 7.7000D-01, 7.7000D-01, 1.0000D+00, & ! 9   
     &     7.7000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 6.7000D-01, 7.7000D-01, 6.5000D-01, 1.0000D+00, & ! O   
     &     3.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 9.6800D-01, 1.0000D+00, 1.0000D+00, 2.4700D-01, & ! +   
     &     7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, & ! 2   
     &     7.2000D-01, 7.4000D-01, 1.4000D-01, 1.4000D-01, 1.4000D-01, & ! +   
     &     1.4000D-01, 1.4000D-01, 1.0000D+00, 1.0000D+00, 1.0000D-04, & ! 3   
     &     1.0000D-04, 5.9500D-02, 1.0000D+00, 1.0000D+00, 2.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 3.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 9   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     1.0000D+00, 1.0000D+00, 3.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     2.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     2.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, & ! 2   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     6.4000D-01, 1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     2.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 9   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00/           !        2   

      DATA ( SC( IRXXN,  2 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! O   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 8.0000D-01, 0.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     7.8400D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 5.0000D-01, & ! +   
     &     1.0000D+00, 6.6000D-01, 3.0000D-01, 0.0000D+00, 1.0000D+00, & ! 2   
     &     2.0000D+00, 1.0000D+00, 2.5000D-01, 2.5000D-01, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 6.7000D-01, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 7.0000D-01, 0.0000D+00, & ! +   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 2.8100D-05, 1.3400D-03, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 3.5000D-01, 5.3000D-01, 1.8000D-01, & ! +   
     &     1.7000D-01, 1.8000D-01, 1.0000D+00, 2.5000D-01, 4.5000D-02, & ! 9   
     &     2.5000D-01, 7.0000D-02, 4.0000D-01, 0.0000D+00, 2.0000D-02, & ! +   
     &     9.5000D-01, 6.5000D-01, 3.5000D-01, 1.0000D+00, 1.0000D+00, & ! O   
     &     4.3000D-01, 2.5000D-02, 5.5000D-02, 7.7000D-01, 1.7000D-01, & ! +   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 1.0000D+00, 4.3000D-01, 0.0000D+00, 6.8700D-01, & ! +   
     &     2.0000D+00, 1.0000D+00, 3.3000D-01, 3.3000D-01, 2.1000D-01, & ! 2   
     &     0.0000D+00, 6.1900D-01, 5.8400D-01, 0.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 7.1900D-01, 1.0000D+00, 6.5000D-01, 4.4000D-01, & ! 3   
     &     4.4000D-01, 0.0000D+00, 1.0000D+00, 3.6000D-01, 6.5000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 6.7000D-01, & ! 4   
     &     1.0000D+00, 1.5000D-01, 3.2000D-01, 7.0000D-02, 2.1450D-01, & ! +   
     &     6.6000D-01, 1.0000D+00, 6.7000D-01, 1.4000D-01, 1.1000D-01, & ! 5   
     &     7.2000D-02, 1.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     2.0000D-01, 5.7000D-01, 8.9000D-01, 2.5000D-02, 5.5000D-02, & ! 6   
     &     6.7000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 3.2000D-01, & ! +   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, 3.3900D-01, & ! 7   
     &     3.2000D-01, 1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, & ! 8   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.3100D-01, & ! +   
     &     5.1000D-02, 1.0000D+00, 9.7000D-01, 9.5000D-01, 1.2000D-03, & ! 9   
     &     1.3000D-03, 1.3000D-03, 1.3000D-03, 7.6300D-01, 8.2000D-01, & ! +   
     &     2.3450D-01, 8.2000D-01, 7.7000D-01, 7.7000D-01, 7.0000D-02, & ! O   
     &     7.7000D-01, 6.4000D-01, 3.6000D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.3000D-01, 3.5000D-01, 5.0000D-01, & ! 1   
     &     7.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 2.8700D-01, 1.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 3.9790D-01, 2.8050D-01, 6.7470D-01, 6.3410D-01, & ! +   
     &     3.9000D-02, 0.0000D+00, 4.5500D-01, 0.0000D+00, 1.0000D-01, & ! 4   
     &     0.0000D+00, 5.0000D-01, 0.0000D+00, 2.5000D-01, 1.0000D-01, & ! +   
     &     4.4000D-01, 4.4000D-01, 1.5000D-01, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.5000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.3700D+00, 7.5000D-01, 8.0000D-02, 1.8000D-02, 1.9500D+00, & ! 7   
     &     1.5000D+00, 7.5000D-01, 1.3700D+00, 1.2846D+00, 1.3218D+00, & ! +   
     &     1.2846D+00, 1.3100D+00, 2.2050D-01, 7.5000D-01, 9.6400D-01, & ! 8   
     &     7.5000D-01, 1.0000D+00, 7.5000D-01, 9.1000D-01, 7.5000D-01, & ! +   
     &     9.0000D-01, 9.0000D-01, 5.0000D-01, 1.0000D+00, 7.5000D-01, & ! 9   
     &     2.6900D-01, 1.5000D+00, 1.1600D+00, 3.0500D-01, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     5.0000D-01, 7.0000D-01, 1.0000D+00, 5.0000D-01, 5.8000D-01, & ! +   
     &     5.1800D-01, 5.0000D-01, 5.0000D-01, 5.0000D-01, 1.0000D+00, & ! 1   
     &     9.1460D-01, 9.5180D-01, 9.1460D-01, 5.0000D-01, 3.7100D-01, & ! +   
     &     5.0000D-01, 6.5000D-01, 5.0000D-01, 5.6000D-01, 5.0000D-01, & ! 2   
     &     4.6000D-01, 5.0000D-01, 2.0000D+00, 1.0000D+00, 5.0000D-01, & ! +   
     &     5.0000D-01, 5.0000D-01, 5.0000D-01, 1.0000D+00, 5.0000D-01, & ! 3   
     &     5.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     2.0000D+00, 1.0000D+00, 1.0000D+00, 1.4000D-01, 5.5000D-02, & ! +   
     &     1.0000D+00, 7.9000D-01, 7.2000D-01, 1.0000D+00, 9.1460D-01, & ! 5   
     &     9.5180D-01, 9.1460D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     6.5000D-01, 4.6000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 5.3800D-01, 1.0000D+00, 7.0000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 2.8700D-01, 1.0000D+00, 2.0000D+00, & ! +   
     &     5.0000D-01, 5.0400D-01, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! 8   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! +   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! 9   
     &     4.8000D-01, 4.8000D-01, 0.0000D+00, 1.4900D-01, 1.0100D-03, & ! +   
     &     4.3000D-01, 5.6000D-01, 3.2000D-01, 4.7700D-01, 3.3300D-01, & ! O   
     &     9.6800D-01, 1.0000D+00, 0.0000D+00, 1.1430D+00, 3.0000D-01, & ! +   
     &     6.3000D-01, 1.0000D-01, 4.9000D-01, 9.2000D-01, 4.0000D-01, & ! 1   
     &     3.1000D-01, 1.0000D+00, 8.0000D-01, 3.8700D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.0000D+00, 2.0000D+00, 1.8000D-02, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 7.5300D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     5.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     2.8000D-01, 2.8000D-01, 2.8000D-01, 2.8000D-01, 2.8000D-01, & ! +   
     &     2.8000D-01, 2.6000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 8.6000D-01, 8.6000D-01, & ! 7   
     &     8.6000D-01, 8.6000D-01, 8.6000D-01, 8.6000D-01, 8.8000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.6000D-01, & ! 9   
     &     9.0480D-01, 1.8000D-03, 9.4050D-01, 1.3105D+00, 9.4050D-01, & ! +   
     &     1.6000D-01, 9.0480D-01, 1.8000D-03, 9.4050D-01, 1.3105D+00, & ! O   
     &     9.4050D-01, 1.6000D-01, 9.0480D-01, 1.8000D-03, 9.4050D-01, & ! +   
     &     1.3105D+00, 9.4050D-01, 2.5810D-01, 3.1300D-01, 6.8950D-01, & ! 1   
     &     7.1000D-03, 5.2000D-01, 1.1900D-01, 1.2910D-01, 4.6670D-01, & ! +   
     &     1.9470D-01, 1.3400D-01, 1.7460D-01, 7.6900D-02, 1.6700D-02, & ! 2   
     &     6.0000D-04, 4.9000D-03, 1.0300D-02, 0.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 0.0000D+00, 6.2000D-01, 1.0000D+00, 0.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.6800D+00, 1.0000D+00, 3.3000D-01, & ! +   
     &     3.3000D-01, 2.0000D-01, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 9   
     &     1.0000D+00, 1.0000D+00, 6.4000D-01, 2.3000D-01, 8.3000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 5.7000D-01, 1.0000D+00, & ! O   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 3.5000D-01, 1.0000D-02, & ! +   
     &     1.0000D-02, 1.0000D-02, 1.0000D+00, 4.0000D-02, 1.0000D+00, & ! 1   
     &     3.5000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     3.3000D-01, 1.0000D+00, 1.0000D+00, 6.5000D-01, 6.1300D-01, & ! 2   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 9.7500D-01, & ! +   
     &     9.4500D-01, 6.5000D-01, 3.1300D-01, 5.2000D-01, 5.2000D-01, & ! 3   
     &     5.6000D-01, 1.0000D+00, 5.7000D-01, 6.6700D-01, 4.7000D-01, & ! +   
     &     8.2000D-01, 8.3000D-01, 1.0000D+00, 1.5200D-01, 2.0000D-01, & ! 4   
     &     1.0000D+00, 8.4000D-01, 4.9000D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     8.4000D-01, 8.4000D-01, 7.4000D-03, 3.6090D-01, 1.1560D-01, & ! +   
     &     8.0830D-01, 2.0360D-01, 8.5400D-02, 7.3800D-02, 3.6950D-01, & ! 6   
     &     2.3000D-03, 1.9620D-01, 4.7900D-02, 3.2400D-02, 1.1700D-02, & ! +   
     &     1.0300D-02, 6.1000D-03, 8.5200D-01, 1.0000D+00, 1.0000D+00, & ! 7   
     &     8.2000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.3100D-01, & ! +   
     &     5.1000D-02, 1.0000D+00, 1.0000D+00, 1.0000D+00, 8.0000D-04, & ! 8   
     &     6.0000D-04, 6.0000D-04, 8.8000D-01, 7.6300D-01, 8.2000D-01, & ! +   
     &     6.9550D-01, 8.2000D-01, 1.0000D+00, 7.7000D-01, 5.7000D-01, & ! 9   
     &     7.7000D-01, 6.4000D-01, 3.6000D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.3000D-01, 3.5000D-01, 5.0000D-01, & ! O   
     &     7.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     0.0000D+00, 1.0000D+00, 9.2000D-01, 1.0000D+00, 7.5300D-01, & ! +   
     &     2.8000D-01, 2.8000D-01, 2.8000D-01, 2.8000D-01, 2.8000D-01, & ! 2   
     &     2.8000D-01, 1.0000D+00, 8.6000D-01, 8.6000D-01, 8.6000D-01, & ! +   
     &     8.6000D-01, 8.6000D-01, 8.6000D-01, 8.8000D-01, 1.8000D-03, & ! 3   
     &     1.8000D-03, 1.8000D-03, 1.3000D-03, 1.0000D+00, 0.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 4   
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     7.1900D-01, 2.3000D-01, 7.0000D-01, 1.0000D+00, 1.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 4.0000D-02, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 9.7500D-01, & ! 8   
     &     9.4500D-01, 6.5000D-01, 1.0000D+00, 3.1300D-01, 5.2000D-01, & ! +   
     &     5.2000D-01, 5.6000D-01, 1.0000D+00, 5.7000D-01, 5.7000D-01, & ! 9   
     &     6.6700D-01, 4.7000D-01, 8.2000D-01, 8.3000D-01, 1.0000D+00, & ! +   
     &     1.5200D-01, 2.0000D-01, 1.0000D+00, 8.4000D-01, 4.9000D-01, & ! O   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 8.2000D-01, 1.0000D+00, & ! +   
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     2.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 2.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     4.0000D-01, 1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 6   
     &     1.0000D+00, 3.0000D+00, 1.0000D+00, 1.0000D+00, 2.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.0000D+00, 1.0000D+00, 1.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 9   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! 1   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00/           !        2   

      DATA ( SC( IRXXN,  3 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     1.2200D+00, 1.0000D+00, 1.0000D+00, 9.0000D-01, 5.0000D-01, & ! +   
     &     1.0000D+00, 6.7000D-01, 7.0000D-01, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 1.0000D+00, 2.0000D-01, 2.0000D-01, 0.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 3   
     &     0.0000D+00, 2.0000D-01, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 3.3000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 3.5000D-01, 5.3000D-01, 1.8000D-01, & ! +   
     &     1.7000D-01, 1.8000D-01, 0.0000D+00, 4.0000D-01, 3.5000D-01, & ! 9   
     &     3.0000D-02, 5.0000D-01, 8.8000D-01, 0.0000D+00, 1.3000D-01, & ! +   
     &     2.0000D-02, 5.0000D-02, 1.7000D-01, 3.0000D-01, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.0000D-01, 1.3000D-01, 1.1000D-01, & ! 2   
     &     0.0000D+00, 1.7000D-01, 1.6000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 1.8400D-01, 0.0000D+00, 3.5000D-01, 7.0000D-02, & ! 3   
     &     7.0000D-02, 0.0000D+00, 0.0000D+00, 6.4000D-01, 3.5000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, 2.7000D-01, & ! 4   
     &     0.0000D+00, 6.3000D-01, 8.0000D-02, 3.2000D-01, 2.5500D-01, & ! +   
     &     3.3000D-01, 6.2000D-01, 2.7000D-01, 1.0000D-01, 2.8000D-01, & ! 5   
     &     8.0000D-03, 6.0000D-01, 6.0000D-01, 1.5000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     2.7000D-01, 1.0000D+00, 0.0000D+00, 0.0000D+00, 6.8000D-01, & ! +   
     &     6.6800D-01, 1.0000D+00, 1.0000D+00, 0.0000D+00, 8.5000D-01, & ! 7   
     &     8.0000D-02, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 4.8000D-02, & ! +   
     &     2.3100D-01, 1.6000D+00, 7.8000D-01, 8.1000D-01, 8.0000D-04, & ! 9   
     &     6.0000D-04, 6.0000D-04, 6.0000D-04, 7.6300D-01, 1.0000D+00, & ! +   
     &     7.0000D-02, 1.0000D+00, 7.7000D-01, 1.0000D+00, 1.3600D+00, & ! O   
     &     1.0000D+00, 6.4000D-01, 3.6000D-01, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 3.3000D-01, 1.6000D-01, 1.0000D+00, 5.0000D-01, & ! 1   
     &     7.0000D-01, 6.1000D-01, 1.0000D+00, 0.0000D+00, 3.2000D-01, & ! +   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 1.2400D+00, 1.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.7710D-01, 2.8050D-01, & ! +   
     &     6.2500D-02, 0.0000D+00, 5.3900D-01, 0.0000D+00, 1.0000D-01, & ! 4   
     &     0.0000D+00, 2.3000D-01, 0.0000D+00, 2.5000D-01, 1.0000D-01, & ! +   
     &     4.4000D-01, 4.4000D-01, 1.5000D-01, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.5000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     6.3000D-01, 7.5000D-01, 2.6000D-02, 1.4000D-01, 1.5000D-01, & ! 7   
     &     7.0500D-01, 1.2800D+00, 3.2000D-01, 3.2000D-01, 3.2000D-01, & ! +   
     &     3.2000D-01, 1.5900D-01, 2.7300D-01, 2.5000D-01, 4.3000D-02, & ! 8   
     &     2.5000D-01, 1.0300D+00, 2.5000D-01, 9.0000D-02, 7.5000D-01, & ! +   
     &     9.0000D-01, 4.0000D-01, 1.5000D+00, 3.3400D-01, 5.0000D-01, & ! 9   
     &     5.0000D-01, 5.0000D-01, 1.1600D+00, 7.7300D-01, 1.0000D+00, & ! +   
     &     1.0000D+00, 3.2000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! O   
     &     9.6500D-01, 1.0000D+00, 0.0000D+00, 1.0000D+00, 2.6000D-02, & ! +   
     &     1.4000D-01, 1.6000D+00, 1.0000D+00, 1.7100D+00, 3.0000D-01, & ! 1   
     &     3.0000D-01, 3.0000D-01, 3.0000D-01, 1.0480D+00, 9.1000D-02, & ! +   
     &     5.0000D-01, 1.0000D+00, 5.0000D-01, 1.0000D+00, 5.0000D-01, & ! 2   
     &     4.6000D-01, 5.0000D-01, 0.0000D+00, 2.0000D+00, 1.0000D+00, & ! +   
     &     3.3000D-01, 5.0000D-01, 2.6900D-01, 5.0000D-01, 1.1600D+00, & ! 3   
     &     5.0000D-01, 2.0000D+00, 1.0000D+00, 3.2000D-01, 1.0000D+00, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.8700D-01, 7.0000D-01, 0.0000D+00, & ! 4   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 9.2000D-02, 2.8000D-01, & ! +   
     &     1.6000D+00, 7.9000D-01, 1.1000D-01, 0.0000D+00, 8.5400D-02, & ! 5   
     &     4.8200D-02, 8.5400D-02, 7.5000D-01, 2.5350D-01, 1.0000D+00, & ! +   
     &     1.7500D-01, 5.4000D-01, 1.0000D+00, 0.0000D+00, 1.0000D+00, & ! 6   
     &     3.3000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 7.0000D-01, & ! +   
     &     6.1000D-01, 1.0000D+00, 0.0000D+00, 3.2000D-01, 0.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 1.2400D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     2.0200D-01, 1.2100D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! +   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! 9   
     &     4.8000D-01, 4.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.8000D-01, 6.8000D-01, 2.5000D-01, 3.3300D-01, & ! O   
     &     8.9500D-01, 9.2500D-01, 0.0000D+00, 8.7000D-01, 8.0000D-01, & ! +   
     &     4.2000D-01, 9.0000D-01, 5.1000D-01, 9.2000D-01, 4.0000D-01, & ! 1   
     &     6.6000D-01, 0.0000D+00, 2.0000D-01, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 7.5300D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, & ! +   
     &     7.2000D-01, 7.4000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.6000D-01, 8.6000D-01, & ! 7   
     &     8.6000D-01, 8.6000D-01, 8.6000D-01, 8.6000D-01, 8.8000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.6000D-01, & ! 9   
     &     3.5700D-02, 1.0000D-04, 5.9500D-02, 3.2000D-01, 3.0000D-01, & ! +   
     &     1.6000D-01, 3.5700D-02, 1.0000D-04, 5.9500D-02, 3.2000D-01, & ! O   
     &     3.0000D-01, 1.6000D-01, 3.5700D-02, 1.0000D-04, 5.9500D-02, & ! +   
     &     3.2000D-01, 3.0000D-01, 2.5810D-01, 1.6000D-03, 2.6420D-01, & ! 1   
     &     1.4460D-01, 8.8000D-03, 7.6730D-01, 1.5470D-01, 1.2220D-01, & ! +   
     &     4.0000D-04, 3.3490D-01, 1.0400D-02, 9.0000D-04, 5.0400D-02, & ! 2   
     &     9.0000D-03, 2.2400D-02, 5.8500D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 5.6000D-01, 0.0000D+00, 0.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 2.8100D-05, 1.3400D-03, 0.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! 9   
     &     7.1900D-01, 0.0000D+00, 3.6000D-01, 7.7000D-01, 1.7000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 4.3000D-01, 0.0000D+00, & ! O   
     &     2.0000D+00, 1.0000D+00, 1.0000D+00, 6.5000D-01, 4.4000D-01, & ! +   
     &     4.4000D-01, 4.4000D-01, 0.0000D+00, 2.0000D-02, 0.0000D+00, & ! 1   
     &     6.5000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     6.7000D-01, 1.0000D+00, 3.5000D-01, 3.1000D-01, 3.8700D-01, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.5000D-01, 2.5000D-02, & ! +   
     &     5.5000D-02, 3.5000D-01, 6.8700D-01, 3.3000D-01, 3.3000D-01, & ! 3   
     &     2.1000D-01, 0.0000D+00, 4.3000D-01, 3.3300D-01, 5.3000D-01, & ! +   
     &     1.8000D-01, 1.7000D-01, 0.0000D+00, 6.1900D-01, 5.8400D-01, & ! 4   
     &     0.0000D+00, 1.6000D-01, 4.9000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0100D-03, & ! 5   
     &     1.6000D-01, 1.6000D-01, 1.1900D-01, 5.2000D-01, 6.8950D-01, & ! +   
     &     3.1300D-01, 7.1000D-03, 2.5810D-01, 4.6670D-01, 1.2910D-01, & ! 6   
     &     1.3400D-01, 1.9470D-01, 1.7460D-01, 7.6900D-02, 1.6700D-02, & ! +   
     &     6.0000D-04, 4.9000D-03, 1.4900D-01, 0.0000D+00, 1.0000D+00, & ! 7   
     &     1.8000D-01, 0.0000D+00, 1.0000D+00, 1.0000D+00, 4.8000D-02, & ! +   
     &     2.3100D-01, 1.6000D+00, 7.8000D-01, 8.1000D-01, 1.0000D+00, & ! 8   
     &     1.0000D+00, 1.0000D+00, 8.8000D-01, 1.9500D-01, 1.8000D-01, & ! +   
     &     2.3450D-01, 1.8000D-01, 7.7000D-01, 2.3000D-01, 7.0000D-02, & ! 9   
     &     2.3000D-01, 0.0000D+00, 3.6000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.0000D+00, 3.3000D-01, 1.6000D-01, 1.0000D+00, 5.0000D-01, & ! O   
     &     7.0000D-01, 6.1000D-01, 0.0000D+00, 0.0000D+00, 3.2000D-01, & ! +   
     &     0.0000D+00, 1.0000D+00, 1.0000D+00, 2.8700D-01, 1.0000D+00, & ! 1   
     &     0.0000D+00, 8.9500D-01, 7.5000D-01, 1.0000D+00, 7.5300D-01, & ! +   
     &     7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, 7.2000D-01, & ! 2   
     &     7.2000D-01, 0.0000D+00, 8.6000D-01, 8.6000D-01, 8.6000D-01, & ! +   
     &     8.6000D-01, 8.6000D-01, 8.6000D-01, 8.8000D-01, 1.0000D-04, & ! 3   
     &     1.0000D-04, 1.0000D-04, 6.0000D-04, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     1.8400D-01, 7.7000D-01, 3.0000D-01, 1.0000D+00, 2.0000D+00, & ! 7   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 2.0000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.5000D-01, 2.5000D-02, & ! 8   
     &     5.5000D-02, 3.5000D-01, 0.0000D+00, 6.8700D-01, 3.3000D-01, & ! +   
     &     3.3000D-01, 2.1000D-01, 0.0000D+00, 4.3000D-01, 4.3000D-01, & ! 9   
     &     3.3300D-01, 5.3000D-01, 1.8000D-01, 1.7000D-01, 0.0000D+00, & ! +   
     &     6.1900D-01, 5.8400D-01, 0.0000D+00, 1.6000D-01, 4.9000D-01, & ! O   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 1.8000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     4.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00, 1.0000D+00, & ! 5   
     &     1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     2.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     2.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN,  4 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     3.5000D-01, 0.0000D+00, 0.0000D+00, 1.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 3.3000D-01, 7.0000D-01, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 8.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 6.5000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.8000D-02, 5.5500D-01, & ! 9   
     &     5.0000D-02, 0.0000D+00, 8.7000D-01, 0.0000D+00, 4.4000D-01, & ! +   
     &     0.0000D+00, 2.0000D-01, 6.5000D-01, 1.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 7.8000D-01, 1.0000D-01, 2.7000D-01, & ! 2   
     &     0.0000D+00, 5.9000D-02, 5.6000D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D-02, & ! 3   
     &     8.0000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.5000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.8000D-01, & ! 4   
     &     0.0000D+00, 1.0000D+00, 6.0000D-02, 7.0000D-02, 2.5350D-01, & ! +   
     &     3.3000D-01, 2.3000D-01, 3.8000D-01, 2.2000D-01, 1.0000D-02, & ! 5   
     &     2.0000D-03, 6.0000D-01, 6.0000D-01, 5.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     3.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 6.8000D-01, & ! +   
     &     3.3200D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.2400D-01, & ! 7   
     &     4.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 8.9000D-02, & ! +   
     &     2.3500D-01, 2.0000D-01, 1.2000D-02, 6.8000D-01, 9.9800D-01, & ! 9   
     &     9.9800D-01, 9.9800D-01, 9.9800D-01, 1.9500D-01, 0.0000D+00, & ! +   
     &     6.0450D-01, 0.0000D+00, 7.7000D-01, 0.0000D+00, 4.3000D-01, & ! O   
     &     0.0000D+00, 0.0000D+00, 3.6000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 6.7000D-01, 1.0000D+00, 6.5000D-01, 1.0000D+00, & ! 1   
     &     1.0000D+00, 3.0000D-02, 0.0000D+00, 0.0000D+00, 6.8000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.6400D-01, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     7.6200D-01, 0.0000D+00, 6.1950D-01, 0.0000D+00, 1.0000D-01, & ! 4   
     &     0.0000D+00, 2.3000D-01, 0.0000D+00, 0.0000D+00, 1.0000D-01, & ! +   
     &     1.5000D-01, 1.5000D-01, 8.5000D-01, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.5000D-01, 2.6000D-02, 1.9100D-01, 2.5000D-01, & ! 7   
     &     4.5000D-02, 2.1800D-01, 0.0000D+00, 8.5400D-02, 4.8200D-02, & ! +   
     &     8.5400D-02, 2.5000D-01, 2.3850D-01, 1.0000D+00, 9.1800D-01, & ! 8   
     &     7.5000D-01, 1.0000D+00, 1.0000D+00, 1.0100D+00, 2.5000D-01, & ! +   
     &     1.0000D+00, 1.0000D+00, 2.5000D-01, 2.5000D-01, 2.5000D-01, & ! 9   
     &     1.6600D+00, 2.5000D-01, 1.5000D+00, 2.0300D-01, 1.0000D+00, & ! +   
     &     0.0000D+00, 6.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     9.3000D-01, 7.0000D-01, 0.0000D+00, 5.0000D-01, 2.6000D-02, & ! +   
     &     1.9100D-01, 2.0000D-01, 9.4000D-01, 2.9000D-01, 0.0000D+00, & ! 1   
     &     8.5400D-02, 4.8200D-02, 8.5400D-02, 2.1900D-01, 1.1650D-01, & ! +   
     &     1.0000D+00, 1.7500D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! 2   
     &     5.4000D-01, 1.0000D+00, 0.0000D+00, 0.0000D+00, 7.5000D-01, & ! +   
     &     3.3400D-01, 5.0000D-01, 5.0000D-01, 5.0000D-01, 1.1600D+00, & ! 3   
     &     3.0000D-02, 0.0000D+00, 0.0000D+00, 6.8000D-01, 1.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.2400D+00, 7.0000D-01, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0300D-01, 4.8500D-01, & ! +   
     &     2.0000D-01, 1.0000D+00, 1.0000D+00, 1.0000D+00, 5.4880D-01, & ! 5   
     &     7.0430D-01, 5.4880D-01, 3.1800D-01, 3.7100D-01, 5.6000D-01, & ! +   
     &     1.6800D-01, 5.4000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     6.7000D-01, 0.0000D+00, 1.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     3.0000D-02, 0.0000D+00, 0.0000D+00, 6.8000D-01, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 4.6400D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     6.4000D-01, 2.8500D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! +   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! 9   
     &     4.8000D-01, 4.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 7.2000D-01, 6.8000D-01, 3.5400D-01, 0.0000D+00, & ! O   
     &     8.9500D-01, 9.2500D-01, 0.0000D+00, 4.6300D-01, 7.4000D-01, & ! +   
     &     8.0000D-02, 0.0000D+00, 0.0000D+00, 7.5000D-01, 3.2000D-01, & ! 1   
     &     3.4000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 8.6000D-01, 8.6000D-01, & ! 7   
     &     8.6000D-01, 8.6000D-01, 8.6000D-01, 0.0000D+00, 8.8000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 9.9800D-01, 4.7020D-01, 5.9500D-02, 5.9500D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 9.9800D-01, 4.7020D-01, 5.9500D-02, & ! O   
     &     5.9500D-02, 0.0000D+00, 0.0000D+00, 9.9800D-01, 4.7020D-01, & ! +   
     &     5.9500D-02, 5.9500D-02, 2.5810D-01, 7.9400D-02, 6.6300D-02, & ! 1   
     &     1.8400D-01, 1.4520D-01, 1.3200D-02, 8.2200D-02, 4.5300D-02, & ! +   
     &     1.0760D-01, 8.2000D-03, 5.5000D-02, 6.2100D-02, 2.7600D-02, & ! 2   
     &     1.4600D-02, 4.2000D-02, 1.7100D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 2.8000D-01, 0.0000D+00, 0.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 3.2000D-01, 7.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     1.8400D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.5000D-01, 7.0000D-02, & ! +   
     &     7.0000D-02, 7.0000D-02, 0.0000D+00, 1.3000D-01, 0.0000D+00, & ! 1   
     &     3.5000D-01, 1.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.7000D-01, 3.0000D-01, 1.7000D-01, 6.6000D-01, 1.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 3.5000D-01, 0.0000D+00, 4.0000D-01, 1.3000D-01, & ! 3   
     &     1.1000D-01, 0.0000D+00, 0.0000D+00, 3.3300D-01, 5.3000D-01, & ! +   
     &     1.8000D-01, 1.7000D-01, 0.0000D+00, 1.7000D-01, 1.6000D-01, & ! 4   
     &     0.0000D+00, 1.6000D-01, 5.1000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     1.6000D-01, 1.6000D-01, 7.6730D-01, 8.8000D-03, 2.6420D-01, & ! +   
     &     1.6000D-03, 1.4460D-01, 2.5810D-01, 1.2220D-01, 1.5470D-01, & ! 6   
     &     3.3490D-01, 4.0000D-04, 1.0400D-02, 9.0000D-04, 5.0400D-02, & ! +   
     &     9.0000D-03, 2.2400D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     1.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 8.9000D-02, & ! +   
     &     2.3500D-01, 2.0000D-01, 1.2000D-02, 6.8000D-01, 9.9800D-01, & ! 8   
     &     9.9800D-01, 9.9800D-01, 1.8000D-01, 2.8400D-01, 1.0000D+00, & ! +   
     &     7.0000D-02, 1.0000D+00, 4.3000D-01, 1.0000D+00, 3.6000D-01, & ! 9   
     &     1.0000D+00, 0.0000D+00, 3.6000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 6.7000D-01, 1.0000D+00, 6.5000D-01, 1.0000D+00, & ! O   
     &     1.0000D+00, 3.0000D-02, 0.0000D+00, 0.0000D+00, 6.8000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.2400D+00, 1.0000D+00, & ! 1   
     &     0.0000D+00, 8.9500D-01, 1.7000D-01, 0.0000D+00, 2.4700D-01, & ! +   
     &     2.8000D-01, 2.8000D-01, 2.8000D-01, 2.8000D-01, 2.8000D-01, & ! 2   
     &     2.8000D-01, 0.0000D+00, 8.6000D-01, 8.6000D-01, 8.6000D-01, & ! +   
     &     8.6000D-01, 8.6000D-01, 1.4000D-01, 0.0000D+00, 1.0000D+00, & ! 3   
     &     1.0000D+00, 1.0000D+00, 9.9800D-01, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 1.3000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-02, 0.0000D+00, & ! 8   
     &     0.0000D+00, 3.5000D-01, 0.0000D+00, 0.0000D+00, 4.0000D-01, & ! +   
     &     1.3000D-01, 1.1000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     3.3300D-01, 5.3000D-01, 1.8000D-01, 1.7000D-01, 0.0000D+00, & ! +   
     &     1.7000D-01, 1.6000D-01, 0.0000D+00, 1.6000D-01, 5.1000D-01, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.8000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN,  5 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     4.3400D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 3.4000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 5.0000D-01, 5.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 3.5000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D-01, 2.6000D-01, & ! 9   
     &     0.0000D+00, 0.0000D+00, 1.8000D-01, 0.0000D+00, 1.1000D-01, & ! +   
     &     0.0000D+00, 1.3000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.0000D-01, 1.0000D-02, 1.0000D-02, & ! 2   
     &     0.0000D+00, 6.1900D-01, 5.8400D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.1000D-01, & ! 3   
     &     4.1000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.1000D-01, & ! 4   
     &     0.0000D+00, 3.7000D-01, 4.0000D-02, 4.0000D-02, 6.5000D-03, & ! +   
     &     2.2000D-01, 2.0000D-02, 2.1000D-01, 5.0000D-01, 5.6000D-01, & ! 5   
     &     1.0000D-01, 1.5000D+00, 1.5000D+00, 4.8000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     2.1000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.2000D-01, & ! +   
     &     3.3200D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.2400D-01, & ! 7   
     &     4.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 9.3500D-01, & ! +   
     &     8.6400D-01, 0.0000D+00, 4.4000D-01, 2.0000D-01, 9.9800D-01, & ! 9   
     &     9.9800D-01, 9.9800D-01, 9.9800D-01, 2.8400D-01, 0.0000D+00, & ! +   
     &     1.5400D-01, 0.0000D+00, 4.3000D-01, 0.0000D+00, 5.0000D-01, & ! O   
     &     0.0000D+00, 0.0000D+00, 3.6000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.6000D-01, 6.5000D-01, 0.0000D+00, & ! 1   
     &     3.0000D-01, 2.7000D-01, 0.0000D+00, 0.0000D+00, 6.8000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     3.2250D-01, 0.0000D+00, 7.7000D-02, 0.0000D+00, 6.0000D-02, & ! 4   
     &     0.0000D+00, 2.7000D-01, 0.0000D+00, 0.0000D+00, 1.0000D-01, & ! +   
     &     4.1000D-01, 4.1000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.5000D-01, 8.2700D-01, 7.7700D-01, 2.5000D-01, & ! 7   
     &     2.5000D-01, 2.5000D-01, 1.0000D+00, 5.4880D-01, 7.0430D-01, & ! +   
     &     5.4880D-01, 2.5000D-01, 8.4500D-02, 0.0000D+00, 6.5900D-01, & ! 8   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 3.1000D-01, 1.0000D+00, & ! +   
     &     1.0000D-01, 1.0000D-01, 2.5000D-01, 2.5000D-01, 2.5000D-01, & ! 9   
     &     6.7000D-02, 2.5000D-01, 1.7500D+00, 5.2500D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 6.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     3.4800D-01, 7.0000D-01, 0.0000D+00, 5.0000D-01, 1.3000D-01, & ! +   
     &     4.2000D-02, 5.0000D-01, 6.0000D-02, 5.0000D-01, 1.0000D+00, & ! 1   
     &     5.4880D-01, 7.0430D-01, 5.4880D-01, 3.0500D-01, 1.7250D-01, & ! +   
     &     5.0000D-01, 1.6800D-01, 5.0000D-01, 1.0000D+00, 5.0000D-01, & ! 2   
     &     5.4000D-01, 5.0000D-01, 0.0000D+00, 0.0000D+00, 5.0000D-01, & ! +   
     &     5.0000D-01, 5.0000D-01, 1.0000D+00, 5.0000D-01, 1.0000D+00, & ! 3   
     &     2.7000D-01, 0.0000D+00, 0.0000D+00, 6.8000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.6400D-01, 7.0000D-01, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 1.8000D-01, 2.0000D-01, 5.0000D-01, 3.6590D-01, & ! 5   
     &     2.4750D-01, 3.6590D-01, 5.0000D-01, 1.1650D-01, 0.0000D+00, & ! +   
     &     7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-01, & ! +   
     &     2.7000D-01, 0.0000D+00, 0.0000D+00, 6.8000D-01, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.4900D-01, 1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! +   
     &     4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, 4.8000D-01, & ! 9   
     &     4.8000D-01, 4.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 6.2000D-01, 3.2000D-01, 2.0400D-01, 0.0000D+00, & ! O   
     &     7.2000D-02, 7.5000D-02, 0.0000D+00, 2.5000D-01, 7.4000D-01, & ! +   
     &     8.3000D-01, 0.0000D+00, 0.0000D+00, 1.7000D-01, 8.0000D-02, & ! 1   
     &     4.3000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 9.9800D-01, 4.7020D-01, 4.7020D-01, 4.7020D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 9.9800D-01, 4.7020D-01, 4.7020D-01, & ! O   
     &     4.7020D-01, 0.0000D+00, 0.0000D+00, 9.9800D-01, 4.7020D-01, & ! +   
     &     4.7020D-01, 4.7020D-01, 1.2910D-01, 7.9400D-02, 6.3500D-02, & ! 1   
     &     1.8400D-01, 4.2000D-02, 6.1400D-02, 4.0000D-02, 1.5160D-01, & ! +   
     &     5.6800D-02, 1.2020D-01, 1.1180D-01, 4.8900D-02, 9.2000D-02, & ! 2   
     &     7.0200D-02, 3.3000D-03, 1.2990D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 3.6000D-01, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 8.0000D-02, & ! +   
     &     8.0000D-02, 8.0000D-02, 0.0000D+00, 4.4000D-01, 0.0000D+00, & ! 1   
     &     3.5000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     3.8000D-01, 1.0000D+00, 6.5000D-01, 3.4000D-01, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 6.5000D-01, 0.0000D+00, 7.8000D-01, 1.0000D-01, & ! 3   
     &     2.7000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.9000D-02, 5.6000D-02, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 1.3200D-02, 1.4520D-01, 6.6300D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 2.5810D-01, 4.5300D-02, 8.2200D-02, & ! 6   
     &     8.2000D-03, 1.0760D-01, 5.5000D-02, 6.2100D-02, 2.7600D-02, & ! +   
     &     1.4600D-02, 4.2000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     1.0000D+00, 0.0000D+00, 4.4000D-01, 2.0000D-01, 0.0000D+00, & ! 8   
     &     8.5200D-02, 4.8100D-02, 5.1000D-01, 8.6500D-02, 0.0000D+00, & ! +   
     &     6.0450D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.3000D-01, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.6000D-01, 6.5000D-01, 0.0000D+00, & ! O   
     &     3.0000D-01, 2.7000D-01, 0.0000D+00, 0.0000D+00, 6.8000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.6400D-01, 0.0000D+00, & ! 1   
     &     0.0000D+00, 7.2000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 1.4000D-01, 1.4000D-01, 1.4000D-01, & ! +   
     &     1.4000D-01, 1.4000D-01, 0.0000D+00, 0.0000D+00, 9.9800D-01, & ! 3   
     &     9.9800D-01, 9.9800D-01, 8.5200D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-02, 0.0000D+00, & ! 8   
     &     0.0000D+00, 6.5000D-01, 0.0000D+00, 0.0000D+00, 7.8000D-01, & ! +   
     &     1.0000D-01, 2.7000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     5.9000D-02, 5.6000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN,  6 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     2.1600D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 6.7000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 5.0000D-01, 5.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 9.0000D-02, 2.8000D-02, & ! 9   
     &     0.0000D+00, 0.0000D+00, 5.1000D-01, 0.0000D+00, 3.0000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.0000D-02, 7.8000D-01, 7.9000D-01, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D-02, & ! 4   
     &     0.0000D+00, 0.0000D+00, 2.0000D-02, 9.0000D-02, 1.1050D-01, & ! +   
     &     1.1000D-01, 1.5000D-01, 2.0000D-02, 4.5000D-01, 1.0000D-01, & ! 5   
     &     2.4300D-01, 5.0000D-01, 5.0000D-01, 7.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     2.0000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 6.8000D-01, & ! +   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0400D-01, & ! +   
     &     1.8000D-02, 0.0000D+00, 6.0000D-02, 9.0000D-02, 0.0000D+00, & ! 9   
     &     8.5200D-02, 4.8100D-02, 8.5200D-02, 8.6500D-02, 0.0000D+00, & ! +   
     &     7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 5.4000D-01, 0.0000D+00, 0.0000D+00, & ! 1   
     &     7.0000D-01, 1.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.8850D-01, 0.0000D+00, 7.0000D-03, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.9800D-01, 2.5100D-01, 0.0000D+00, & ! 7   
     &     2.5000D-01, 2.5000D-01, 5.0000D-01, 3.6590D-01, 2.4750D-01, & ! +   
     &     3.6590D-01, 2.5000D-01, 1.1050D-01, 0.0000D+00, 5.3300D-01, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.1000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.2500D-01, 0.0000D+00, 0.0000D+00, & ! 9   
     &     2.5000D-01, 0.0000D+00, 5.0000D-01, 1.3500D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     2.5000D-01, 3.0000D-01, 0.0000D+00, 0.0000D+00, 2.7300D-01, & ! +   
     &     3.8100D-01, 5.0000D-01, 5.0000D-01, 5.0000D-01, 5.0000D-01, & ! 1   
     &     3.6590D-01, 2.4750D-01, 3.6590D-01, 5.0000D-01, 4.0300D-01, & ! +   
     &     0.0000D+00, 7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     1.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     5.0000D-01, 0.0000D+00, 5.0000D-01, 5.0000D-01, 2.3000D+00, & ! 3   
     &     7.0000D-01, 0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 5.0000D-01, 3.0000D-01, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 5.1900D-01, 2.4000D-02, & ! +   
     &     0.0000D+00, 2.0000D-02, 8.5000D-01, 2.5000D-01, 3.6590D-01, & ! 5   
     &     3.5220D-01, 4.5730D-01, 2.4000D-02, 1.7250D-01, 0.0000D+00, & ! +   
     &     7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D-01, & ! +   
     &     7.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.5000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     4.0000D-02, 4.0000D-02, 4.8000D-01, 4.0000D-02, 4.0000D-02, & ! +   
     &     4.8000D-01, 4.0000D-02, 4.0000D-02, 4.8000D-01, 4.0000D-02, & ! 9   
     &     4.0000D-02, 4.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.0000D-01, 6.8000D-01, 1.5000D-01, 0.0000D+00, & ! O   
     &     3.2000D-02, 0.0000D+00, 0.0000D+00, 2.3100D-01, 1.8500D-01, & ! +   
     &     1.7000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     4.0000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 5.9400D-02, 4.7020D-01, 4.7020D-01, 4.7020D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 5.9400D-02, 4.7020D-01, 4.7020D-01, & ! O   
     &     4.7020D-01, 0.0000D+00, 0.0000D+00, 4.6930D-01, 4.7020D-01, & ! +   
     &     4.7020D-01, 4.7020D-01, 2.5810D-01, 7.9400D-02, 9.5000D-02, & ! 1   
     &     1.8400D-01, 5.2400D-02, 3.0500D-02, 1.0680D-01, 5.0000D-03, & ! +   
     &     2.1220D-01, 7.6100D-02, 1.3000D-03, 1.5500D-01, 7.8800D-02, & ! 2   
     &     1.5300D-02, 9.0600D-02, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 2.8000D-01, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.1000D-01, & ! +   
     &     4.1000D-01, 4.1000D-01, 0.0000D+00, 1.1000D-01, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.1000D-01, 0.0000D+00, 0.0000D+00, 4.3000D-01, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 3.5000D-01, 0.0000D+00, 1.0000D-01, 1.0000D-02, & ! 3   
     &     1.0000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 6.1900D-01, 5.8400D-01, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 6.1400D-02, 4.2000D-02, 6.3500D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 1.2910D-01, 1.5160D-01, 4.0000D-02, & ! 6   
     &     1.2020D-01, 5.6800D-02, 1.1180D-01, 4.8900D-02, 9.2000D-02, & ! +   
     &     7.0200D-02, 3.3000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0400D-01, & ! +   
     &     1.8000D-02, 0.0000D+00, 6.0000D-02, 9.0000D-02, 9.9800D-01, & ! 8   
     &     5.4770D-01, 7.0290D-01, 1.2000D-01, 1.2850D-01, 0.0000D+00, & ! +   
     &     1.5400D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-01, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 5.4000D-01, 0.0000D+00, 0.0000D+00, & ! O   
     &     7.0000D-01, 1.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.9400D-02, & ! 3   
     &     5.9400D-02, 4.6930D-01, 5.4770D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.1000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 3.5000D-01, 0.0000D+00, 0.0000D+00, 1.0000D-01, & ! +   
     &     1.0000D-02, 1.0000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     6.1900D-01, 5.8400D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN,  7 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 6.7000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.2000D-01, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 6.8000D-02, 3.7000D-01, 1.1050D-01, & ! +   
     &     3.3000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.4000D-01, & ! 5   
     &     8.0000D-02, 5.0000D-02, 5.0000D-02, 2.5000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3200D-01, & ! +   
     &     4.5000D-02, 0.0000D+00, 1.3000D-01, 2.0000D-02, 9.9800D-01, & ! 9   
     &     5.4770D-01, 7.0290D-01, 5.4770D-01, 1.2850D-01, 0.0000D+00, & ! +   
     &     7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     3.0000D-01, 7.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.3800D-01, 0.0000D+00, 7.0000D-03, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.9700D-01, 6.1800D-01, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 2.5000D-01, 3.6590D-01, 3.5220D-01, & ! +   
     &     4.5730D-01, 2.3000D-02, 1.8350D-01, 0.0000D+00, 7.0000D-03, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     2.5000D-01, 0.0000D+00, 2.5000D-01, 1.0500D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     2.5000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 6.6200D-01, & ! +   
     &     8.2400D-01, 0.0000D+00, 5.0000D-01, 0.0000D+00, 2.5000D-01, & ! 1   
     &     3.6590D-01, 3.5220D-01, 4.5730D-01, 5.0000D-01, 1.0000D+00, & ! +   
     &     0.0000D+00, 7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-01, & ! 3   
     &     1.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 5.0000D-01, 1.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4700D-01, 2.4100D-01, & ! +   
     &     0.0000D+00, 9.0000D-02, 4.0000D-02, 2.5000D-01, 5.4880D-01, & ! 5   
     &     5.9960D-01, 4.5730D-01, 3.3000D-02, 4.0300D-01, 0.0000D+00, & ! +   
     &     1.7500D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-01, & ! +   
     &     1.8000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 4.0000D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     4.0000D-02, 0.0000D+00, 0.0000D+00, 4.0000D-02, 0.0000D+00, & ! 9   
     &     0.0000D+00, 4.0000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.7000D-02, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.7000D-02, 6.0000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 4.6930D-01, 4.7020D-01, 4.7020D-01, 4.7020D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.6930D-01, 4.7020D-01, 4.7020D-01, & ! O   
     &     4.7020D-01, 0.0000D+00, 0.0000D+00, 4.6930D-01, 4.7020D-01, & ! +   
     &     4.7020D-01, 4.7020D-01, 2.5810D-01, 3.9700D-02, 3.3200D-02, & ! 1   
     &     9.2000D-02, 6.4600D-02, 3.6700D-02, 1.0680D-01, 1.3320D-01, & ! +   
     &     7.0000D-04, 3.2500D-02, 1.4800D-01, 1.5500D-02, 1.0800D-02, & ! 2   
     &     1.0170D-01, 1.4890D-01, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 1.0000D-01, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-02, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.0000D-02, 0.0000D+00, 0.0000D+00, 4.0000D-02, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D-02, 7.8000D-01, & ! 3   
     &     7.9000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 3.0500D-02, 5.2400D-02, 9.5000D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 2.5810D-01, 5.0000D-03, 1.0680D-01, & ! 6   
     &     7.6100D-02, 2.1220D-01, 1.3000D-03, 1.5500D-01, 7.8800D-02, & ! +   
     &     1.5300D-02, 9.0600D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.3200D-01, & ! +   
     &     4.5000D-02, 0.0000D+00, 1.3000D-01, 2.0000D-02, 4.9900D-01, & ! 8   
     &     3.6510D-01, 2.4700D-01, 0.0000D+00, 3.0800D-01, 0.0000D+00, & ! +   
     &     7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     3.0000D-01, 7.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.6930D-01, & ! 3   
     &     4.6930D-01, 4.6930D-01, 3.6510D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D-02, & ! +   
     &     7.8000D-01, 7.9000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN,  8 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.4000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 4.3000D-01, 2.6000D-02, 4.5900D-01, & ! +   
     &     1.0000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D-02, & ! 5   
     &     4.2000D-01, 5.0000D-02, 5.0000D-02, 2.5000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.6500D-01, & ! +   
     &     2.0300D-01, 0.0000D+00, 3.0000D-02, 5.0000D-02, 4.9900D-01, & ! 9   
     &     3.6510D-01, 2.4700D-01, 3.6510D-01, 3.0800D-01, 0.0000D+00, & ! +   
     &     1.6450D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 2.1000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.0500D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 5.0000D-02, 2.5000D-01, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 2.5000D-01, 5.4880D-01, 5.9960D-01, & ! +   
     &     4.5730D-01, 1.8000D-02, 6.7000D-02, 0.0000D+00, 1.1900D-01, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     2.5000D-01, 0.0000D+00, 2.5000D-01, 2.5000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     5.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 6.7000D-02, & ! +   
     &     5.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.5000D-01, & ! 1   
     &     5.4880D-01, 5.9960D-01, 4.5730D-01, 0.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 1.7500D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0830D+00, & ! 3   
     &     1.0500D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.5000D-02, 6.0000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 3.1000D-02, 9.1000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.1000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.9000D-02, 1.5000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 4.6930D-01, 0.0000D+00, 4.7020D-01, 4.7020D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.6930D-01, 0.0000D+00, 4.7020D-01, & ! O   
     &     4.7020D-01, 0.0000D+00, 0.0000D+00, 4.6930D-01, 0.0000D+00, & ! +   
     &     4.7020D-01, 4.7020D-01, 2.5810D-01, 7.9400D-02, 3.9000D-02, & ! 1   
     &     1.8400D-01, 6.4600D-02, 6.0500D-02, 1.0680D-01, 1.0400D-02, & ! +   
     &     1.3690D-01, 1.2670D-01, 4.1700D-02, 1.1170D-01, 1.2270D-01, & ! 2   
     &     4.7000D-03, 1.7420D-01, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 2.0000D-01, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 3.6700D-02, 6.4600D-02, 3.3200D-02, & ! +   
     &     3.9700D-02, 9.2000D-02, 2.5810D-01, 1.3320D-01, 1.0680D-01, & ! 6   
     &     3.2500D-02, 7.0000D-04, 1.4800D-01, 1.5500D-02, 1.0800D-02, & ! +   
     &     1.0170D-01, 1.4890D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.6500D-01, & ! +   
     &     2.0300D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4950D-01, & ! 8   
     &     3.6510D-01, 3.5150D-01, 0.0000D+00, 7.1500D-02, 0.0000D+00, & ! +   
     &     7.0000D-03, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 2.1000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.6930D-01, & ! 3   
     &     4.6930D-01, 4.6930D-01, 4.5640D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN,  9 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.0000D-01, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 2.0000D-02, 1.0000D-02, 2.3050D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D-02, & ! 5   
     &     2.8000D-02, 8.0000D-02, 8.0000D-02, 1.1000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.2000D-02, & ! +   
     &     3.3000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4950D-01, & ! 9   
     &     3.6510D-01, 3.5150D-01, 4.5640D-01, 7.1500D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 2.5000D-01, 2.5000D-01, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 1.6000D-02, 1.2790D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 2.9200D-01, 2.5000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-01, & ! +   
     &     5.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D-01, & ! 1   
     &     7.0000D-01, 7.0000D-01, 7.0000D-01, 0.0000D+00, 1.1200D-01, & ! +   
     &     0.0000D+00, 5.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-01, & ! 3   
     &     5.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 9.5000D-02, 6.3000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.1200D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 7.0000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 4.6930D-01, 0.0000D+00, 0.0000D+00, 7.0000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.6930D-01, 0.0000D+00, 0.0000D+00, & ! O   
     &     7.0000D-01, 0.0000D+00, 0.0000D+00, 4.6930D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 7.0000D-01, 2.5810D-01, 7.9400D-02, 3.9000D-02, & ! 1   
     &     1.8400D-01, 6.4600D-02, 4.8800D-02, 5.3400D-02, 3.2100D-02, & ! +   
     &     8.1800D-02, 6.4700D-02, 1.1810D-01, 1.1300D-02, 1.5010D-01, & ! 2   
     &     1.6580D-01, 1.1570D-01, 8.2600D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.4000D-01, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 6.0500D-02, 6.4600D-02, 3.9000D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 2.5810D-01, 1.0400D-02, 1.0680D-01, & ! 6   
     &     1.2670D-01, 1.3690D-01, 4.1700D-02, 1.1170D-01, 1.2270D-01, & ! +   
     &     4.7000D-03, 1.7420D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.2000D-02, & ! +   
     &     3.3000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4950D-01, & ! 8   
     &     5.4770D-01, 5.9840D-01, 0.0000D+00, 8.4000D-02, 0.0000D+00, & ! +   
     &     1.6450D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.6930D-01, & ! 3   
     &     4.6930D-01, 4.6930D-01, 4.5640D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 3.4000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 10 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D-01, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 1.5000D-02, 1.0000D-02, 5.0000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D-01, & ! 5   
     &     4.9100D-01, 6.5000D-01, 7.0000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 6.5000D-02, & ! +   
     &     2.1700D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4950D-01, & ! 9   
     &     5.4770D-01, 5.9840D-01, 4.5640D-01, 8.4000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 2.5000D-01, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 5.0000D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     5.0000D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.4700D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 4.6930D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 4.6930D-01, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 2.5810D-01, 7.9400D-02, 3.9000D-02, & ! 1   
     &     1.8400D-01, 3.2300D-02, 2.0200D-02, 1.0680D-01, 1.4000D-02, & ! +   
     &     8.1800D-02, 3.6400D-02, 2.1800D-02, 1.4650D-01, 4.6600D-02, & ! 2   
     &     1.7570D-01, 2.9400D-02, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D-01, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 4.8800D-02, 6.4600D-02, 3.9000D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 2.5810D-01, 3.2100D-02, 5.3400D-02, & ! 6   
     &     6.4700D-02, 8.1800D-02, 1.1810D-01, 1.1300D-02, 1.5010D-01, & ! +   
     &     1.6580D-01, 1.1570D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.1700D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.3700D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.6930D-01, & ! 3   
     &     4.6930D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.4000D-01, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 11 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.8000D-01, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 6.0000D-03, 9.0000D-02, 1.7850D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     3.0000D-03, 0.0000D+00, 6.5000D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     3.3000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 4.8000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.9400D-02, 1.9500D-02, & ! 1   
     &     1.8400D-01, 6.4600D-02, 5.8400D-02, 1.0680D-01, 2.4000D-02, & ! +   
     &     8.1800D-02, 8.9000D-02, 1.1950D-01, 8.7100D-02, 7.8500D-02, & ! 2   
     &     1.5870D-01, 4.1500D-02, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.0200D-02, 3.2300D-02, 3.9000D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 2.5810D-01, 1.4000D-02, 1.0680D-01, & ! 6   
     &     3.6400D-02, 8.1800D-02, 2.1800D-02, 1.4650D-01, 4.6600D-02, & ! +   
     &     1.7570D-01, 2.9400D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     3.3000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 12 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 1.5300D-01, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 3.2000D-02, 4.5700D-01, 6.6500D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     4.4000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.7200D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 2.7500D-01, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.9400D-02, 3.9000D-02, & ! 1   
     &     1.8400D-01, 6.4600D-02, 3.7500D-02, 1.0680D-01, 5.9500D-02, & ! +   
     &     4.0900D-02, 6.4600D-02, 2.0700D-02, 4.0000D-03, 3.9900D-02, & ! 2   
     &     5.7300D-02, 5.6000D-02, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 5.8400D-02, 6.4600D-02, 1.9500D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 0.0000D+00, 2.4000D-02, 1.0680D-01, & ! 6   
     &     8.9000D-02, 8.1800D-02, 1.1950D-01, 8.7100D-02, 7.8500D-02, & ! +   
     &     1.5870D-01, 4.1500D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     2.7200D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 13 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 5.6000D-01, 7.3000D-01, 3.8500D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     1.3600D-01, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 7.9400D-02, 3.9000D-02, & ! 1   
     &     1.8400D-01, 6.4600D-02, 7.4000D-03, 1.0680D-01, 4.4500D-02, & ! +   
     &     8.1800D-02, 4.4800D-02, 2.9400D-02, 5.6900D-02, 7.8600D-02, & ! 2   
     &     9.7000D-03, 5.7100D-02, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 3.7500D-02, 6.4600D-02, 3.9000D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 0.0000D+00, 5.9500D-02, 1.0680D-01, & ! 6   
     &     6.4600D-02, 4.0900D-02, 2.0700D-02, 4.0000D-03, 3.9900D-02, & ! +   
     &     5.7300D-02, 5.6000D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 14 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 1.0000D-02, 1.1000D-01, 2.8000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.9000D-02, & ! 1   
     &     0.0000D+00, 6.4600D-02, 1.7900D-02, 1.0680D-01, 4.4500D-02, & ! +   
     &     8.1800D-02, 8.7000D-03, 4.5300D-02, 4.7200D-02, 7.8600D-02, & ! 2   
     &     2.5300D-02, 5.7100D-02, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 7.4000D-03, 6.4600D-02, 3.9000D-02, & ! +   
     &     7.9400D-02, 1.8400D-01, 0.0000D+00, 4.4500D-02, 1.0680D-01, & ! 6   
     &     4.4800D-02, 8.1800D-02, 2.9400D-02, 5.6900D-02, 7.8600D-02, & ! +   
     &     9.7000D-03, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 15 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 4.4000D-01, 1.7000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.9000D-02, & ! 1   
     &     0.0000D+00, 6.4600D-02, 2.5800D-02, 1.0680D-01, 4.4500D-02, & ! +   
     &     8.1800D-02, 1.0300D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! 2   
     &     3.6800D-02, 5.7100D-02, 1.6520D-01, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 1.7900D-02, 6.4600D-02, 3.9000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 1.0680D-01, & ! 6   
     &     8.7000D-03, 8.1800D-02, 4.5300D-02, 4.7200D-02, 7.8600D-02, & ! +   
     &     2.5300D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 16 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 3.0000D-02, 4.4000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.9000D-02, & ! 1   
     &     0.0000D+00, 6.4600D-02, 2.4400D-02, 0.0000D+00, 2.2200D-02, & ! +   
     &     8.1800D-02, 2.9700D-02, 6.2800D-02, 6.9800D-02, 3.9300D-02, & ! 2   
     &     5.2600D-02, 2.8500D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.5800D-02, 6.4600D-02, 3.9000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 1.0680D-01, & ! 6   
     &     1.0300D-02, 8.1800D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! +   
     &     3.6800D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 17 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 2.0000D-02, 1.7000D-02, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 3.9000D-02, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 4.4500D-02, & ! +   
     &     8.1800D-02, 2.9700D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! 2   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 6.4600D-02, 3.9000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 2.2200D-02, 0.0000D+00, & ! 6   
     &     2.9700D-02, 8.1800D-02, 6.2800D-02, 6.9800D-02, 3.9300D-02, & ! +   
     &     5.2600D-02, 2.8500D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 18 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 6.0000D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 4.4500D-02, & ! +   
     &     8.1800D-02, 2.9700D-02, 3.1400D-02, 3.4900D-02, 7.8600D-02, & ! 2   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 3.9000D-02, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 0.0000D+00, & ! 6   
     &     2.9700D-02, 8.1800D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! +   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 19 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 1.0000D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 1.2200D-02, 0.0000D+00, 4.4500D-02, & ! +   
     &     0.0000D+00, 1.4900D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! 2   
     &     2.6300D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 0.0000D+00, & ! 6   
     &     2.9700D-02, 8.1800D-02, 3.1400D-02, 3.4900D-02, 7.8600D-02, & ! +   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 20 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 3.0000D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 4.4500D-02, & ! +   
     &     0.0000D+00, 2.9700D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! 2   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 1.2200D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 0.0000D+00, & ! 6   
     &     1.4900D-02, 0.0000D+00, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! +   
     &     2.6300D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 21 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 6.0000D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 4.4500D-02, & ! +   
     &     0.0000D+00, 2.9700D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! 2   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 0.0000D+00, & ! 6   
     &     2.9700D-02, 0.0000D+00, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! +   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 22 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 4.4500D-02, & ! +   
     &     0.0000D+00, 2.9700D-02, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! 2   
     &     5.2600D-02, 4.9700D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 0.0000D+00, & ! 6   
     &     2.9700D-02, 0.0000D+00, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! +   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 23 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.9700D-02, 6.2800D-02, 6.9800D-02, 0.0000D+00, & ! 2   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 4.4500D-02, 0.0000D+00, & ! 6   
     &     2.9700D-02, 0.0000D+00, 6.2800D-02, 6.9800D-02, 7.8600D-02, & ! +   
     &     5.2600D-02, 4.9700D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 24 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.9700D-02, 6.2800D-02, 6.9800D-02, 0.0000D+00, & ! 2   
     &     5.2600D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     2.9700D-02, 0.0000D+00, 6.2800D-02, 6.9800D-02, 0.0000D+00, & ! +   
     &     5.2600D-02, 5.7100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 25 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 2.9700D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     5.3100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     2.9700D-02, 0.0000D+00, 6.2800D-02, 6.9800D-02, 0.0000D+00, & ! +   
     &     5.2600D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 26 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     5.2600D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 2.4400D-02, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     2.9700D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     5.3100D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      DATA ( SC( IRXXN, 27 ), IRXXN = 1, NRXNS ) / & 
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     5.2600D-02, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 2   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 3   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 4   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 5   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 6   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 7   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 8   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 9   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! O   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! 1   
     &     0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, 0.0000D+00, & ! +   
     &     0.0000D+00/           !        2   

      INTEGER            :: NREACT( NRXNS )

      DATA ( NREACT( IRXXN ), IRXXN = 1, NRXNS ) / & 
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! O   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 1   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 2   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 3   
     &      1,    1,    1,    2,    2,    2,    2,    1,    2,    1, & ! 4   
     &      1,    1,    1,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    1,    1,    2,    1,    2,    2,    2, & ! 7   
     &      1,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    1,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 7   
     &      1,    2,    1,    2,    1,    2,    1,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    1,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    1,    1, & ! 2   
     &      1,    1,    1,    1,    1,    2,    1,    1,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    1,    1,    1,    1, & ! 6   
     &      1,    1,    1,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 2   
     &      1,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 3   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 4   
     &      1,    1,    2,    1,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    1,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    1,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    1, & ! 3   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 4   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 7   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 9   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! O   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 1   
     &      1,    1,    1,    2,    2,    2,    2,    2,    2,    2, & ! 2   
     &      1,    1,    2,    2,    2,    2,    2,    2,    2,    2, & ! 3   
     &      2,    1,    2,    1,    2,    2,    1,    2,    2,    2, & ! 4   
     &      1,    1,    1,    1,    1,    1,    1,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    1,    1,    1,    1,    2,    2,    2, & ! 7   
     &      2,    2,    1,    1,    2,    2,    1,    1,    2,    2, & ! 8   
     &      2,    2,    2,    2,    1,    1,    2,    2,    1,    2, & ! 9   
     &      2,    2,    2,    2,    2,    1,    1,    1,    1,    1, & ! O   
     &      1,    2,    2,    2,    2,    2,    3,    2,    2,    2, & ! 1   
     &      1/     !  2   
      INTEGER            :: NPRDCT( NRXNS )

      DATA ( NPRDCT( IRXXN ), IRXXN = 1, NRXNS ) / & 
     &      1,    1,    1,    2,    1,    2,    2,    2,    4,    1, & ! O   
     &      2,    3,    3,    2,    2,    6,    3,    3,    4,    3, & ! 1   
     &      3,    7,    4,    1,    2,    2,    3,    6,    6,    2, & ! 2   
     &      3,    3,    3,    3,    3,    2,    4,    2,    3,    2, & ! 3   
     &      2,    3,    2,    1,    1,    1,    1,    1,    0,    1, & ! 4   
     &      1,    1,    1,    0,    1,    1,    1,    1,    1,    2, & ! 5   
     &      1,    1,    1,    1,    1,    1,    1,    2,    3,    1, & ! 6   
     &      2,    1,    1,    2,    1,    1,    2,    1,    3,    2, & ! 7   
     &      1,    1,    2,    2,    1,    1,    1,    5,    3,    3, & ! 8   
     &      3,    3,    2,   12,    6,    4,    3,    6,    1,   10, & ! 9   
     &      3,    5,    4,    4,    2,    2,    2,    2,    2,    2, & ! O   
     &      2,    1,    1,    1,    1,    1,    2,    2,    1,    2, & ! 1   
     &      2,    2,    6,    6,    6,    1,    5,    5,    1,    2, & ! 2   
     &      2,    3,    2,    3,    5,    5,    1,    2,    3,    4, & ! 3   
     &      3,    3,    3,    2,    6,    2,    5,   21,   17,   14, & ! 4   
     &      8,    6,    6,    6,   10,   12,   10,   11,    9,    1, & ! 5   
     &      2,    2,    2,    2,    2,    6,    3,    2,    2,    6, & ! 6   
     &      6,    3,    3,    2,    5,    5,    2,    3,    1,    1, & ! 7   
     &      1,    1,    2,    1,    2,    1,    2,    3,    3,   10, & ! 8   
     &     13,    4,    8,    8,   10,   10,   10,   10,   10,    3, & ! 9   
     &      8,    3,    5,    3,    5,    3,    3,    5,    3,    3, & ! O   
     &      3,    4,    6,    5,    4,    7,    8,    3,    2,    5, & ! 1   
     &      2,    3,    3,    4,    3,    1,    1,    1,    1,    1, & ! 2   
     &      1,    1,    1,    1,    1,    1,    2,    2,    3,    3, & ! 3   
     &      8,    1,    7,    1,    5,    1,    5,    1,    3,    5, & ! 4   
     &      5,    5,    4,    1,    1,    1,    1,    1,    1,    6, & ! 5   
     &      1,    1,    1,    1,    1,    1,    1,    1,    1,    1, & ! 6   
     &      3,    5,   10,    9,    5,    6,    6,    8,    8,    8, & ! 7   
     &      8,    9,    9,    4,    8,    5,    4,    4,    6,    5, & ! 8   
     &      5,    5,    6,    5,    5,    8,    5,    9,    9,    4, & ! 9   
     &      3,    5,    3,    3,    3,    8,    6,    2,    5,   10, & ! O   
     &      9,    6,    7,    6,    9,    9,    9,    9,    7,   10, & ! 1   
     &      5,    9,    5,    5,    5,    6,    5,    2,    3,    5, & ! 2   
     &      6,    5,    6,    6,    9,   10,    3,    3,    6,    4, & ! 3   
     &      3,    3,    7,    7,    2,    2,    3,    3,    9,   12, & ! 4   
     &      4,    7,    7,    7,    7,    7,    7,    8,    9,    4, & ! 5   
     &      7,    4,    3,    2,    3,    4,    3,    4,    4,    7, & ! 6   
     &      8,    3,    2,    5,    2,    3,    3,    4,    3,    2, & ! 7   
     &      6,    5,    1,    2,    0,    6,    6,    7,    6,    6, & ! 8   
     &      7,    6,    6,    7,    6,    6,    7,    1,    2,    2, & ! 9   
     &      2,    6,    6,    7,    3,    6,    5,    1,    8,    9, & ! O   
     &      6,    3,    3,    5,    5,    6,    2,    3,    3,    1, & ! 1   
     &      1,    2,    2,    2,    1,    1,    1,    3,    1,    1, & ! 2   
     &      2,    0,    1,    1,    1,    1,    1,    2,    1,    1, & ! 3   
     &      1,    1,    1,    1,    1,    3,    3,    3,    3,    3, & ! 4   
     &      3,    3,    2,    2,    2,    2,    2,    2,    2,    1, & ! 5   
     &      1,    1,    1,    1,    1,    1,    2,    2,    2,    2, & ! 6   
     &      2,    2,    2,    4,    4,    4,    4,    4,    3,    4, & ! 7   
     &      3,    3,    3,    3,    3,    3,    3,    1,    1,    1, & ! 8   
     &      1,    1,    1,    1,    3,    3,   10,    7,    8,    9, & ! 9   
     &      3,    3,   10,    7,    8,    9,    3,    3,    9,    7, & ! O   
     &      8,    9,   10,   13,   17,   13,   16,   25,   15,   22, & ! 1   
     &     18,   25,   24,   24,   22,   26,   23,   15,    1,    2, & ! 2   
     &      2,    1,    8,    2,    1,    4,    4,    4,    6,    2, & ! 3   
     &      2,    2,    1,    2,    2,    2,    2,    2,    2,    2, & ! 4   
     &      2,    3,    1,    2,    2,    1,    2,    1,    1,    1, & ! 5   
     &      2,    2,    1,    1,    2,    2,    2,    1,    2,    1, & ! 6   
     &      1,    2,    1,    1,    2,    2,    1,    2,    1,    1, & ! 7   
     &      2,    2,    1,    2,    2,    2,    3,    3,    2,    3, & ! 8   
     &      2,    2,    2,    3,    3,    4,    2,    3,    3,    3, & ! 9   
     &      2,    2,    3,    3,    2,    3,    3,    3,    4,    6, & ! O   
     &      6,    6,    2,   11,    2,    5,    4,    4,    3,    3, & ! 1   
     &      7,    5,    5,    7,    4,    2,    2,    2,    5,    3, & ! 2   
     &      3,    6,    3,    7,    7,    7,    2,    3,    4,    4, & ! 3   
     &      4,    4,    2,    6,    6,    2,    4,    4,    2,    2, & ! 4   
     &      2,    2,    2,    2,    3,    4,    4,   26,   17,   18, & ! 5   
     &     14,   14,   11,   23,   16,   26,   19,   25,   25,   23, & ! 6   
     &     27,   24,    3,    2,    3,    4,    2,    3,    3,    9, & ! 7   
     &     12,    4,    7,    7,    9,    9,    9,    6,   10,    4, & ! 8   
     &      9,    4,    4,    4,    6,    4,    2,    4,    2,    2, & ! 9   
     &      3,    4,    6,    5,    4,    7,    8,    2,    2,    5, & ! O   
     &      2,    3,    3,    5,    4,    1,    5,    4,    3,    4, & ! 1   
     &      4,    4,    4,    4,    4,    4,    2,    5,    5,    5, & ! 2   
     &      5,    5,    4,    3,   10,   10,    9,    9,    4,    1, & ! 3   
     &      2,    2,    2,    2,    2,    2,    2,    2,    2,    3, & ! 4   
     &      2,    1,    1,    1,    2,    2,    1,    2,    1,    2, & ! 5   
     &      1,    2,    1,    1,    1,    3,    2,    2,    3,    3, & ! 6   
     &      4,    3,    3,    3,    3,    4,    3,    2,   11,    2, & ! 7   
     &      2,    2,    2,    5,    3,    3,    6,    2,    3,    7, & ! 8   
     &      7,    7,    2,    3,    3,    4,    4,    4,    4,    2, & ! 9   
     &      6,    6,    2,    4,    4,    3,    2,    1,    4,    2, & ! O   
     &      1,    2,    1,    2,    2,    2,    2,    2,    2,    2, & ! 1   
     &      2,    2,    1,    2,    2,    1,    1,    1,    1,    2, & ! 2   
     &      2,    2,    2,    2,    1,    1,    1,    2,    2,    1, & ! 3   
     &      3,    2,    2,    1,    1,    1,    1,    2,    2,    2, & ! 4   
     &      2,    2,    2,    3,    3,    3,    3,    2,    2,    2, & ! 5   
     &      2,    2,    2,    2,    2,    3,    2,    4,    3,    3, & ! 6   
     &      3,    2,    2,    2,    2,    1,    1,    1,    1,    1, & ! 7   
     &      1,    1,    2,    2,    2,    2,    2,    2,    2,    2, & ! 8   
     &      1,    1,    2,    2,    2,    2,    2,    2,    1,    1, & ! 9   
     &      1,    2,    2,    2,    2,    0,    0,    0,    0,    0, & ! O   
     &      0,    1,    1,    2,    2,    2,    3,    1,    2,    3, & ! 1   
     &      2/     !  2   

      INTEGER, PARAMETER :: MHETERO =  60
      INTEGER            :: IHETERO( MHETERO,2 )

      DATA ( IHETERO( IRXXN,1 ), IRXXN = 1, MHETERO ) / & 
     &    429,  430,  431,  432,  433,  434,  435,  436,  437,  438, & ! O   
     &    540,  541,  874,  875,  876,  877,  878,  879,  880,  881, & ! 1   
     &    882,  883,  884,  885,  886,  887,  888,  889,  890,  891, & ! 2   
     &    892,  893,  894,  895,  896,  897,  898,  899,  900,  901, & ! 3   
     &    902,  903,  904,  905,  906,  907,  908,  909,  910,  911, & ! 4   
     &    912,  913,  914,  915,  916,  917,  918,  919,  920,  921/!5   

      DATA ( IHETERO( IRXXN,2 ), IRXXN = 1, MHETERO ) / & 
     &      1,    2,    3,    4,    5,    6,    7,    8,    6,    6, & ! O   
     &      9,    9,   10,   11,   12,   13,   14,   14,   15,   16, & ! 1   
     &     17,   18,   19,   20,   21,   22,   23,   24,   25,   26, & ! 2   
     &     27,   28,   29,   30,   31,   32,   33,   34,   35,   36, & ! 3   
     &     37,   38,   39,   40,   41,   42,   43,   44,   45,   46, & ! 4   
     &     47,   48,   49,   50,   51,   52,   53,   54,   55,   56/!5   

      INTEGER, PARAMETER :: NPHOTAB =  74
      CHARACTER( 16 )    :: PHOTAB( NPHOTAB )

      DATA ( PHOTAB( IRXXN ), IRXXN = 1, NPHOTAB ) / & 
     &   'O3O3P_NASA06    ', 'O3O1D_NASA06    ', 'H2O2_RACM2      ', & 
     &   'NO2_RACM2       ', 'NO3NO_RACM2     ', 'NO3NO2_RACM2    ', & 
     &   'HONO_RACM2      ', 'HNO3_RACM2      ', 'HNO4_RACM2      ', & 
     &   'HCHO_MOL_JPL19  ', 'HCHO_RAD_JPL19  ', 'CH3CHO_RACM2    ', & 
     &   'ALD_JPL19       ', 'CH3COCH3A_JPL19 ', 'CH3COCH3B_JPL19 ', & 
     &   'UALD_RACM2      ', 'MEK_JGR19       ', 'KET_JGR19       ', & 
     &   'HKET_RACM2      ', 'MACR_RACM2      ', 'MVK_JPL19       ', & 
     &   'GLYH2_JPL19     ', 'GLYF_JPL19      ', 'GLYHX_JPL19     ', & 
     &   'MGLY_RACM2      ', 'BALD1_CALVERT11 ', 'BALD2_CALVERT11 ', & 
     &   'OP1_RACM2       ', 'PAA_RACM2       ', 'ONIT_CALVERT08  ', & 
     &   'PAN1_JPL19      ', 'PAN2_JPL19      ', 'PPN1_JPL19      ', & 
     &   'PPN2_JPL19      ', 'TRPN_WANG2023   ', 'ACRO_09         ', & 
     &   'CL2_JPL19       ', 'CLO_JPL19       ', 'OCLO_JPL19      ', & 
     &   'CL2O2_JPL19     ', 'HOCL_JPL19      ', 'CLNO_JPL19      ', & 
     &   'CLNO2_JPL19     ', 'CLNO3_R_JPL19   ', 'CLNO3_M_JPL19   ', & 
     &   'HCOCL_JPL19     ', 'BR2_JPL19       ', 'BRO_JPL19       ', & 
     &   'OBRO_JPL19      ', 'HOBR_JPL19      ', 'BRNO_JPL19      ', & 
     &   'BRNO2_JPL19     ', 'BRNO3_R_JPL19   ', 'BRNO3_M_JPL19   ', & 
     &   'CH2BR2_JPL19    ', 'CHBR3_JPL19     ', 'HCOBR_JPL19     ', & 
     &   'I2_JPL19        ', 'IO_JPL19        ', 'OIO_JPL19       ', & 
     &   'INO3_06         ', 'HOI_JPL19       ', 'HI_JPL19        ', & 
     &   'INO_JPL19       ', 'INO2_JPL19      ', 'CH3I_JPL19      ', & 
     &   'CH2I2_JPL19     ', 'BRCL_JPL19      ', 'ICL_JPL19       ', & 
     &   'IBR_JPL19       ', 'CH2IBR_JPL19    ', 'CH2ICL_JPL19    ', & 
     &   'CHBR2CL_JPL19   ', 'CHBRCL2_JPL19   '/

      INTEGER, PARAMETER :: NHETERO =  56
      CHARACTER( 16 )    :: HETERO( NHETERO )

      DATA ( HETERO( IRXXN ), IRXXN = 1, NHETERO ) / & 
     &   'HETERO_GLY      ', 'HETERO_MGLY     ', 'HETERO_NO2      ', &
     &   'HETERO_HO2      ', 'HETERO_NO3      ', 'HETERO_IEPOX    ', &
     &   'HETERO_ISO3NOSJ ', 'HETERO_ISO3OSJ  ', 'HETERO_ANO3     ', &
     &   'HETERO_N2O5IJ   ', 'HETERO_N2O5K    ', 'HETERO_H2NO3PAIJ', &
     &   'HETERO_H2NO3PAK ', 'HETERO_H2NO3PBIJ', 'HETERO_H2NO3PBK ', &
     &   'HETERO_HOCL_ACLJ', 'HETERO_HOCL_ABRJ', 'HETERO_CLN2_WAI ', &
     &   'HETERO_CLN2_WAJ ', 'HETERO_CLN2_ACLJ', 'HETERO_CLN2_ABRJ', &
     &   'HETERO_CLN3_WAI ', 'HETERO_CLN3_WAJ ', 'HETERO_CLN3_ACLJ', &
     &   'HETERO_CLN3_ABRJ', 'HETERO_HOBR_ACLJ', 'HETERO_HOBR_ABRJ', &
     &   'HETERO_BRN2_ACLJ', 'HETERO_BRN2_ABRJ', 'HETERO_BRN3_WAI ', &
     &   'HETERO_BRN3_WAJ ', 'HETERO_BRN3_ACLJ', 'HETERO_BRN3_ABRJ', &
     &   'HETERO_HBR_ABRJ ', 'HETERO_HOI_ACLJ ', 'HETERO_HOI_ABRJ ', &
     &   'HETERO_INO2_ACLJ', 'HETERO_INO2_ABRJ', 'HETERO_INO3_ACLJ', &
     &   'HETERO_INO3_ABRJ', 'HETERO_I2O2_AI  ', 'HETERO_I2O2_AJ  ', &
     &   'HETERO_I2O3_AI  ', 'HETERO_I2O3_AJ  ', 'HETERO_I2O4_AI  ', &
     &   'HETERO_I2O4_AJ  ', 'HETERO_SO2H     ', 'HETERO_SO2O     ', &
     &   'HETERO_SO2M     ', 'HETERO_SO2P     ', 'HETERO_SO2N     ', &
     &   'HETERO_SO2T     ', 'HETERO_HMSP     ', 'HETERO_HNO4     ', &
     &   'HETERO_HMSL1    ', 'HETERO_HMSL2    '/

      CHARACTER( 16 )    :: RXLABEL( NRXNS )

      DATA ( RXLABEL( IRXXN ), IRXXN = 1, NRXNS ) / & 
     &    'R001            ', 'R002            ', 'R003            ', & ! 0   
     &    'R004            ', 'R005            ', 'R006            ', & ! 1   
     &    'R007            ', 'R008            ', 'R009            ', & ! 2   
     &    'R010            ', 'R011            ', 'R012            ', & ! 3   
     &    'R013            ', 'R014            ', 'R014a           ', & ! 4   
     &    'R015            ', 'TRP01           ', 'TRP02           ', & ! 5   
     &    'R016            ', 'R017            ', 'R018            ', & ! 6   
     &    'R019            ', 'R020            ', 'R021            ', & ! 7   
     &    'R022            ', 'R023            ', 'R024            ', & ! 8   
     &    'R025            ', 'R026            ', 'R027a           ', & ! 9   
     &    'R027b           ', 'R028            ', 'R029            ', & ! 0   
     &    'TRP03           ', 'R029a           ', 'R030            ', & ! 1   
     &    'R031            ', 'R032            ', 'R033            ', & ! 2   
     &    'R033a           ', 'R033b           ', 'TRP55           ', & ! 3   
     &    'TRP56           ', 'R034            ', 'R035            ', & ! 4   
     &    'R036            ', 'R037            ', 'R038            ', & ! 5   
     &    'R039            ', 'R040            ', 'R041            ', & ! 6   
     &    'R042            ', 'R043            ', 'R044            ', & ! 7   
     &    'R045            ', 'R046            ', 'R047            ', & ! 8   
     &    'R048            ', 'R049            ', 'R050            ', & ! 9   
     &    'R051            ', 'R052            ', 'R053            ', & ! 0   
     &    'R054            ', 'R055            ', 'R056            ', & ! 1   
     &    'R057            ', 'R058            ', 'R059            ', & ! 2   
     &    'R060            ', 'R061            ', 'R062            ', & ! 3   
     &    'R063            ', 'R064            ', 'R065            ', & ! 4   
     &    'R066            ', 'R067            ', 'R068            ', & ! 5   
     &    'R069            ', 'R070            ', 'R071            ', & ! 6   
     &    'R072            ', 'R073            ', 'R074            ', & ! 7   
     &    'R076            ', 'R077            ', 'R078            ', & ! 8   
     &    'R080            ', 'ROCARO31        ', 'ROCARO41        ', & ! 9   
     &    'ROCARO51        ', 'ROCARO61        ', 'ROCARO61a       ', & ! 0   
     &    'RAM01           ', 'RAM02           ', 'RAM03           ', & ! 1   
     &    'RAM04           ', 'RAM05           ', 'RAM06           ', & ! 2   
     &    'RAM07           ', 'RAM08           ', 'RAM09           ', & ! 3   
     &    'RAM10           ', 'RAM11           ', 'RAM12           ', & ! 4   
     &    'RAM13           ', 'R087            ', 'R088            ', & ! 5   
     &    'TRP04           ', 'TRP05           ', 'R089            ', & ! 6   
     &    'R090            ', 'R091            ', 'R092            ', & ! 7   
     &    'R093            ', 'R094            ', 'R095            ', & ! 8   
     &    'R096            ', 'R097            ', 'R098            ', & ! 9   
     &    'R099            ', 'R100            ', 'R101            ', & ! 0   
     &    'R102            ', 'R103            ', 'R104            ', & ! 1   
     &    'R105            ', 'R106            ', 'R108            ', & ! 2   
     &    'R109            ', 'R110            ', 'R111            ', & ! 3   
     &    'R112            ', 'R113            ', 'R114            ', & ! 4   
     &    'TRP06           ', 'R116            ', 'R117            ', & ! 5   
     &    'R118            ', 'R119            ', 'R120            ', & ! 6   
     &    'R121            ', 'R122            ', 'R123            ', & ! 7   
     &    'TRP07           ', 'TRP57           ', 'R126            ', & ! 8   
     &    'R127            ', 'R128            ', 'R131            ', & ! 9   
     &    'R132            ', 'TRP08           ', 'TRP09           ', & ! 0   
     &    'R133            ', 'R134            ', 'R135            ', & ! 1   
     &    'R136            ', 'R137            ', 'R138            ', & ! 2   
     &    'R140            ', 'R141            ', 'R142            ', & ! 3   
     &    'R143            ', 'R146            ', 'R147            ', & ! 4   
     &    'TRP10           ', 'R148            ', 'R149            ', & ! 5   
     &    'R150            ', 'R151            ', 'R152            ', & ! 6   
     &    'R153            ', 'R154            ', 'RBAL1           ', & ! 7   
     &    'R155            ', 'R156            ', 'R158            ', & ! 8   
     &    'R159            ', 'RBAL2           ', 'TRP11           ', & ! 9   
     &    'TRP12           ', 'R166            ', 'R167            ', & ! 0   
     &    'R168            ', 'R169            ', 'R170            ', & ! 1   
     &    'R171            ', 'R172            ', 'R173            ', & ! 2   
     &    'R174            ', 'R175            ', 'R177            ', & ! 3   
     &    'R178            ', 'R179            ', 'ROCARO33        ', & ! 4   
     &    'ROCARO43        ', 'ROCARO53        ', 'ROCARO63        ', & ! 5   
     &    'R189            ', 'TRP13           ', 'TRP14           ', & ! 6   
     &    'TRP15           ', 'R190            ', 'TRP16           ', & ! 7   
     &    'TRP17           ', 'TRP18           ', 'TRP19           ', & ! 8   
     &    'TRP20           ', 'R191            ', 'R192            ', & ! 9   
     &    'R193            ', 'R194            ', 'R195            ', & ! 0   
     &    'R196            ', 'R197            ', 'R198            ', & ! 1   
     &    'R199            ', 'R200            ', 'R201            ', & ! 2   
     &    'R202            ', 'R203            ', 'R204            ', & ! 3   
     &    'R205            ', 'R206            ', 'R207            ', & ! 4   
     &    'R208            ', 'R209            ', 'R210            ', & ! 5   
     &    'R211            ', 'R212            ', 'R213            ', & ! 6   
     &    'R214            ', 'R215            ', 'R217            ', & ! 7   
     &    'R218            ', 'R219            ', 'ROCARO32        ', & ! 8   
     &    'ROCARO42        ', 'ROCARO52        ', 'ROCARO62        ', & ! 9   
     &    'R229            ', 'TRP21           ', 'TRP22           ', & ! 0   
     &    'TRP23           ', 'R230            ', 'TRP24           ', & ! 1   
     &    'TRP25           ', 'TRP26           ', 'TRP27           ', & ! 2   
     &    'TRP28           ', 'R231            ', 'R232            ', & ! 3   
     &    'R233            ', 'R234            ', 'R235            ', & ! 4   
     &    'R236            ', 'R237            ', 'R238            ', & ! 5   
     &    'R239            ', 'RBAL3           ', 'RBAL4           ', & ! 6   
     &    'RBAL5           ', 'R240            ', 'R241            ', & ! 7   
     &    'R242            ', 'R243            ', 'R244            ', & ! 8   
     &    'R245            ', 'R246            ', 'R247            ', & ! 9   
     &    'R248            ', 'R249            ', 'R250            ', & ! 0   
     &    'R251            ', 'R253            ', 'R254            ', & ! 1   
     &    'R255            ', 'ROCARO35        ', 'ROCARO45        ', & ! 2   
     &    'ROCARO55        ', 'ROCARO65        ', 'R264            ', & ! 3   
     &    'R265            ', 'TRP29           ', 'TRP30           ', & ! 4   
     &    'TRP31           ', 'R266            ', 'TRP32           ', & ! 5   
     &    'TRP33           ', 'TRP34           ', 'R267            ', & ! 6   
     &    'R268            ', 'R269            ', 'R270            ', & ! 7   
     &    'R271            ', 'R272            ', 'R273            ', & ! 8   
     &    'R274            ', 'R275            ', 'R276            ', & ! 9   
     &    'R277            ', 'R278            ', 'R279            ', & ! 0   
     &    'R280            ', 'R281            ', 'R282            ', & ! 1   
     &    'R283            ', 'R284            ', 'R285            ', & ! 2   
     &    'R286            ', 'R287            ', 'R289            ', & ! 3   
     &    'R290            ', 'R291            ', 'ROCARO36        ', & ! 4   
     &    'ROCARO46        ', 'ROCARO56        ', 'ROCARO66        ', & ! 5   
     &    'R300            ', 'R301            ', 'TRP35           ', & ! 6   
     &    'TRP36           ', 'TRP37           ', 'R302            ', & ! 7   
     &    'TRP38           ', 'TRP39           ', 'TRP40           ', & ! 8   
     &    'R303            ', 'R304            ', 'R305            ', & ! 9   
     &    'R306            ', 'R307            ', 'R308            ', & ! 0   
     &    'R309            ', 'R310            ', 'R311            ', & ! 1   
     &    'R312            ', 'R313            ', 'R314            ', & ! 2   
     &    'R315            ', 'R316            ', 'R317            ', & ! 3   
     &    'R318            ', 'R319            ', 'R320            ', & ! 4   
     &    'R321            ', 'R322            ', 'R323            ', & ! 5   
     &    'R324            ', 'R325            ', 'R327            ', & ! 6   
     &    'R328            ', 'R329            ', 'ROCARO34        ', & ! 7   
     &    'ROCARO44        ', 'ROCARO54        ', 'ROCARO64        ', & ! 8   
     &    'R338            ', 'R339            ', 'R340            ', & ! 9   
     &    'TRP53           ', 'TRP54           ', 'R341            ', & ! 0   
     &    'R342            ', 'R343            ', 'R344            ', & ! 1   
     &    'R345            ', 'R346            ', 'R347            ', & ! 2   
     &    'R348            ', 'R349            ', 'R350            ', & ! 3   
     &    'R351            ', 'R352            ', 'R353            ', & ! 4   
     &    'R354            ', 'R355            ', 'R356            ', & ! 5   
     &    'R357            ', 'R358            ', 'R359            ', & ! 6   
     &    'R360            ', 'R361            ', 'R362            ', & ! 7   
     &    'R363            ', 'TRP41           ', 'TRP42           ', & ! 8   
     &    'TRP43           ', 'TRP44           ', 'TRP45           ', & ! 9   
     &    'TRP46           ', 'TRP47           ', 'TRP48           ', & ! 0   
     &    'TRP49           ', 'TRP50           ', 'TRP51           ', & ! 1   
     &    'TRP52           ', 'RAM17           ', 'R001c           ', & ! 2   
     &    'R002c           ', 'T17             ', 'T18             ', & ! 3   
     &    'T19             ', 'T20             ', 'T10             ', & ! 4   
     &    'T10a            ', 'T10b            ', 'T10c            ', & ! 5   
     &    'T10d            ', 'T10e            ', 'T11             ', & ! 6   
     &    'T12             ', 'R003c           ', 'R004c           ', & ! 7   
     &    'R005c           ', 'R006c           ', 'R007c           ', & ! 8   
     &    'R008c           ', 'R010c           ', 'R011c           ', & ! 9   
     &    'R012c           ', 'R013c           ', 'R014c           ', & ! 0   
     &    'R015c           ', 'R016c           ', 'R017c           ', & ! 1   
     &    'R019c           ', 'R020c           ', 'HET_GLY         ', & ! 2   
     &    'HET_MGLY        ', 'HET_NO2         ', 'HET_HO2         ', & ! 3   
     &    'HET_NO3         ', 'HET_IEPOX       ', 'HET_ISO3TET     ', & ! 4   
     &    'HET_IEPOXOS     ', 'HET_IPX         ', 'HET_INALD       ', & ! 5   
     &    'ROCALK1c        ', 'ROCALK2c        ', 'ROCALK3c        ', & ! 6   
     &    'ROCALK4c        ', 'ROCALK5c        ', 'ROCALK6c        ', & ! 7   
     &    'HC1001          ', 'ROCALK7c        ', 'ROCALK8c        ', & ! 8   
     &    'ROCALK9c        ', 'ROCALK10c       ', 'ROCALK11c       ', & ! 9   
     &    'ROCALK12c       ', 'HC1002          ', 'ROCALK13c       ', & ! 0   
     &    'ROCALK14c       ', 'ROCALK15c       ', 'ROCALK16c       ', & ! 1   
     &    'ROCALK17c       ', 'ROCALK18c       ', 'HC1003          ', & ! 2   
     &    'ROCALK19c       ', 'ROCALK20c       ', 'ROCALK21c       ', & ! 3   
     &    'ROCALK22c       ', 'ROCALK23c       ', 'ROCALK24c       ', & ! 4   
     &    'HC1004          ', 'ROCALK25c       ', 'ROCALK26c       ', & ! 5   
     &    'ROCALK27c       ', 'ROCALK28c       ', 'ROCALK29c       ', & ! 6   
     &    'ROCALK30c       ', 'HC1005          ', 'ROCALK31c       ', & ! 7   
     &    'ROCALK32c       ', 'ROCALK33c       ', 'ROCALK34c       ', & ! 8   
     &    'ROCALK35c       ', 'ROCALK36c       ', 'HC1006          ', & ! 9   
     &    'ROCALK37c       ', 'ROCALK38c       ', 'ROCALK39c       ', & ! 0   
     &    'ROCALK40c       ', 'ROCALK41c       ', 'ROCALK42c       ', & ! 1   
     &    'HC1007          ', 'ROCALK43c       ', 'ROCALK44c       ', & ! 2   
     &    'ROCALK45c       ', 'ROCALK46c       ', 'ROCALK47c       ', & ! 3   
     &    'ROCALK48c       ', 'HC1008          ', 'ROCARO01        ', & ! 4   
     &    'ROCARO02        ', 'ROCARO03        ', 'ROCARO04        ', & ! 5   
     &    'ROCARO05        ', 'ROCARO06        ', 'ROCARO11        ', & ! 6   
     &    'ROCARO12        ', 'ROCARO13        ', 'ROCARO14        ', & ! 7   
     &    'ROCARO15        ', 'ROCARO16        ', 'ROCARO21        ', & ! 8   
     &    'ROCARO22        ', 'ROCARO23        ', 'ROCARO24        ', & ! 9   
     &    'ROCARO25        ', 'ROCARO26        ', 'ROCOXY1c        ', & ! 0   
     &    'ROCOXY2c        ', 'ROCOXY3c        ', 'ROCOXY4c        ', & ! 1   
     &    'ROCOXY5c        ', 'ROCOXY6c        ', 'ROCOXY7c        ', & ! 2   
     &    'ROCOXY8c        ', 'ROCOXY9c        ', 'ROCOXY10c       ', & ! 3   
     &    'ROCOXY11c       ', 'ROCOXY12c       ', 'ROCOXY13c       ', & ! 4   
     &    'ROCOXY14c       ', 'ROCOXY15c       ', 'ROCOXY16c       ', & ! 5   
     &    'R364            ', 'TRP58           ', 'TRP59           ', & ! 6   
     &    'ROCARO71        ', 'ROCARO71a       ', 'ROCARO71b       ', & ! 7   
     &    'ROCARO72        ', 'ROCARO73        ', 'ROCARO74        ', & ! 8   
     &    'ROCARO75        ', 'ROCARO76        ', 'HET_ANO3I       ', & ! 9   
     &    'HET_ANO3J       ', 'SOAphot01       ', 'CLP01           ', & ! 0   
     &    'CLP02           ', 'CLP03           ', 'CLP07           ', & ! 1   
     &    'CLP10           ', 'CLP12           ', 'CLP14           ', & ! 2   
     &    'CLP17           ', 'CLP18           ', 'CLP19           ', & ! 3   
     &    'CLT01           ', 'CLT02           ', 'CLT03           ', & ! 4   
     &    'CLT04           ', 'CLT05           ', 'CLT06           ', & ! 5   
     &    'CLT07           ', 'CLT08           ', 'CLT09           ', & ! 6   
     &    'CLT10           ', 'CLT11           ', 'CLT12           ', & ! 7   
     &    'CLT14           ', 'CLT15           ', 'CLT17           ', & ! 8   
     &    'CLT22           ', 'CLT23           ', 'CLT25           ', & ! 9   
     &    'CLT27           ', 'CLT28           ', 'CLT31           ', & ! 0   
     &    'CLT32           ', 'CLT35           ', 'CLT36           ', & ! 1   
     &    'CLT37           ', 'CLT38           ', 'CLT39           ', & ! 2   
     &    'CLT40           ', 'CLT41           ', 'CLT43           ', & ! 3   
     &    'CLT44           ', 'CLT45           ', 'CLTO01          ', & ! 4   
     &    'CLTO02          ', 'CLTO03          ', 'CLTO04          ', & ! 5   
     &    'CLTO05          ', 'CLTO06          ', 'CLTO07          ', & ! 6   
     &    'CLTO08          ', 'CLTO09          ', 'CLTO10          ', & ! 7   
     &    'CLTO11          ', 'CLTO12          ', 'CLTO13          ', & ! 8   
     &    'CLTO14          ', 'CLTO15          ', 'CTLO16          ', & ! 9   
     &    'CLTO17          ', 'CLTO18          ', 'CLTO19          ', & ! 0   
     &    'CLTO20          ', 'CLTO21          ', 'CLTO22          ', & ! 1   
     &    'CLTO23          ', 'CLTO24          ', 'CLTO25          ', & ! 2   
     &    'CLTO26          ', 'CLTO27          ', 'CLTO28          ', & ! 3   
     &    'CLTO29          ', 'CLTO29a         ', 'CLTO30          ', & ! 4   
     &    'CLTO31          ', 'CLTO32          ', 'CLTO33          ', & ! 5   
     &    'CLTO34          ', 'CLTO35          ', 'CLTO36          ', & ! 6   
     &    'CLTO37          ', 'CLTO38          ', 'CLTO39          ', & ! 7   
     &    'CLTO40          ', 'CLTO41          ', 'CLTO42          ', & ! 8   
     &    'CLTO43          ', 'CLTO44          ', 'CLTO45          ', & ! 9   
     &    'CLTO46          ', 'CLTO47          ', 'CLTO48          ', & ! 0   
     &    'CLTO49          ', 'CLTO50          ', 'CLTO51          ', & ! 1   
     &    'CLTO52          ', 'CLTO53          ', 'CLTO54          ', & ! 2   
     &    'CLTO55          ', 'CLTO56          ', 'CLTO58          ', & ! 3   
     &    'CLTO59          ', 'CLTO60          ', 'CLTO61          ', & ! 4   
     &    'CLTO62          ', 'CLTO63          ', 'CLTO64          ', & ! 5   
     &    'CLTO65          ', 'CLTO66          ', 'CLTO67          ', & ! 6   
     &    'CLTO68          ', 'CLTO69          ', 'CLTO70          ', & ! 7   
     &    'CLTO71          ', 'CLTO72          ', 'CLTO73          ', & ! 8   
     &    'CLTO74          ', 'CLTO75          ', 'CLTO76          ', & ! 9   
     &    'CLTO77          ', 'CLTO78          ', 'CLTO79          ', & ! 0   
     &    'CLTO80          ', 'CLTO81          ', 'CLTO82          ', & ! 1   
     &    'CLTO83          ', 'CLTO84          ', 'CLTO85          ', & ! 2   
     &    'CLTO86          ', 'CLTO87          ', 'CLTO88          ', & ! 3   
     &    'CLTO89          ', 'CLTO90          ', 'CLTO91          ', & ! 4   
     &    'CLT092          ', 'CLT093          ', 'CLOR01          ', & ! 5   
     &    'CLOR02          ', 'CLOR03          ', 'CLOR04          ', & ! 6   
     &    'CLOR05          ', 'CLOR06          ', 'CLOR07          ', & ! 7   
     &    'CLOR08          ', 'CLOR09          ', 'CLOR10          ', & ! 8   
     &    'CLOR12          ', 'CLOR13          ', 'CLOR14          ', & ! 9   
     &    'CLOR15          ', 'CLOR16          ', 'CLOR17          ', & ! 0   
     &    'CLOR18          ', 'CLOR19          ', 'CLOR20          ', & ! 1   
     &    'CLOR21          ', 'CLOR22          ', 'CLOR23          ', & ! 2   
     &    'CLOR24          ', 'CLOR25          ', 'CLOR26          ', & ! 3   
     &    'CLOR27          ', 'CLOR28          ', 'CLOR29          ', & ! 4   
     &    'CLOR30          ', 'CLOR31          ', 'CLOR32          ', & ! 5   
     &    'CLOR33          ', 'CLOR34          ', 'CLOR35          ', & ! 6   
     &    'CLOR36          ', 'CLOR37          ', 'CLOR38          ', & ! 7   
     &    'CLOR39          ', 'CLOR40          ', 'CLOR41          ', & ! 8   
     &    'CLOR42          ', 'CLOR43          ', 'CLOR44          ', & ! 9   
     &    'CLOR45          ', 'CLOR46          ', 'CLOR47          ', & ! 0   
     &    'CLOR48          ', 'CLOR49          ', 'CLOR50          ', & ! 1   
     &    'CLOR51          ', 'CLOR52          ', 'CLOR53          ', & ! 2   
     &    'CLOR54          ', 'CLOR55          ', 'CLOR56          ', & ! 3   
     &    'CLOR57          ', 'CLOR58          ', 'CLOR59          ', & ! 4   
     &    'CLOR60          ', 'CLOR61          ', 'CLOR62          ', & ! 5   
     &    'CLOR63          ', 'BRP01           ', 'BRP02           ', & ! 6   
     &    'BRP03           ', 'BRP04           ', 'BRP06           ', & ! 7   
     &    'BRP08           ', 'BRP12           ', 'BRP13           ', & ! 8   
     &    'BRP15           ', 'BRP16           ', 'BRP17           ', & ! 9   
     &    'BRT01           ', 'BRT02           ', 'BRT03           ', & ! 0   
     &    'BRT06           ', 'BRT07           ', 'BRT08           ', & ! 1   
     &    'BRT10           ', 'BRT11           ', 'BRT12           ', & ! 2   
     &    'BRT13           ', 'BRT15           ', 'BRT16           ', & ! 3   
     &    'BRT17           ', 'BRT18           ', 'BRT19           ', & ! 4   
     &    'BRTO01          ', 'BRTO02          ', 'BRTO03          ', & ! 5   
     &    'BRTO04          ', 'BRTO05          ', 'BRTO06          ', & ! 6   
     &    'BRTO07          ', 'BRTO08          ', 'BRTO09          ', & ! 7   
     &    'BRTO10          ', 'BRTO11          ', 'BRTO12          ', & ! 8   
     &    'BRTO13          ', 'BRTO13a         ', 'BRTO14          ', & ! 9   
     &    'BRTO15          ', 'BRTO16          ', 'BRTO17          ', & ! 0   
     &    'BRTO18          ', 'BRTO19          ', 'BRTO20          ', & ! 1   
     &    'BRTO21          ', 'BRTO22          ', 'BRTO23          ', & ! 2   
     &    'BRTO24          ', 'BRTO25          ', 'BRTO26          ', & ! 3   
     &    'BRTO27          ', 'BRTO28          ', 'BRTO29          ', & ! 4   
     &    'BRTO30          ', 'BRTO31          ', 'BRTO32          ', & ! 5   
     &    'BRTO34          ', 'BRTO35          ', 'BRT036          ', & ! 6   
     &    'BRTO37          ', 'BRTO38          ', 'BRTO39          ', & ! 7   
     &    'BRTO40          ', 'BRTO41          ', 'BRTO42          ', & ! 8   
     &    'BRTO43          ', 'BRT044          ', 'BRT045          ', & ! 9   
     &    'IP01            ', 'IP02            ', 'IP03            ', & ! 0   
     &    'IP05            ', 'IP09            ', 'IP10            ', & ! 1   
     &    'IP11            ', 'IP12            ', 'IP16            ', & ! 2   
     &    'IP19            ', 'IP20            ', 'IO21            ', & ! 3   
     &    'IO22            ', 'IT02            ', 'IT03            ', & ! 4   
     &    'IT04            ', 'IT05            ', 'IT06            ', & ! 5   
     &    'IT07            ', 'IT09            ', 'IT10            ', & ! 6   
     &    'IT11            ', 'IT12            ', 'IT13            ', & ! 7   
     &    'IT15            ', 'IT16            ', 'IT20            ', & ! 8   
     &    'IT21            ', 'IT22            ', 'IT23            ', & ! 9   
     &    'IT25            ', 'IT30            ', 'IT33            ', & ! 0   
     &    'IT34            ', 'IT35            ', 'IT36            ', & ! 1   
     &    'IT37            ', 'IT38            ', 'IT39            ', & ! 2   
     &    'ITO01           ', 'MP01            ', 'MP03            ', & ! 3   
     &    'MP04            ', 'MP06            ', 'MP07            ', & ! 4   
     &    'MP08            ', 'MP09            ', 'MT01            ', & ! 5   
     &    'MT02            ', 'MT03            ', 'MT04            ', & ! 6   
     &    'MT05            ', 'MT06            ', 'MT07            ', & ! 7   
     &    'MT08            ', 'MTO01           ', 'MTO02           ', & ! 8   
     &    'MTO03           ', 'MTO06           ', 'MTO07           ', & ! 9   
     &    'MTO08           ', 'MTO09           ', 'MTO10           ', & ! 0   
     &    'HET_N2O5IJ      ', 'HET_N2O5K       ', 'HET_H2NO3PIJA   ', & ! 1   
     &    'HET_H2NO3PKA    ', 'HET_H2NO3PIB    ', 'HET_H2NO3PJB    ', & ! 2   
     &    'HET_H2NO3PKB    ', 'HET_HOCL_ACLJ   ', 'HET_HOCL_ABRJ   ', & ! 3   
     &    'HET_CLN2_WAI    ', 'HET_CLN2_WAJ    ', 'HET_CLN2_ACLJ   ', & ! 4   
     &    'HET_CLN2_ABRJ   ', 'HET_CLN3_WAI    ', 'HET_CLN3_WAJ    ', & ! 5   
     &    'HET_CLN3_ACLJ   ', 'HET_CLN3_ABRJ   ', 'HET_HOBR_ACLJ   ', & ! 6   
     &    'HET_HOBR_ABRJ   ', 'HET_BRN2_ACLJ   ', 'HET_BRN2_ABRJ   ', & ! 7   
     &    'HET_BRN3_WAI    ', 'HET_BRN3_WAJ    ', 'HET_BRN3_ACLJ   ', & ! 8   
     &    'HET_BRN3_ABRJ   ', 'HET_HBR_ABRJ    ', 'HET_HOI_ACLJ    ', & ! 9   
     &    'HET_HOI_ABRJ    ', 'HET_INO2_ACLJ   ', 'HET_INO2_ABRJ   ', & ! 0   
     &    'HET_INO3_ACLJ   ', 'HET_INO3_ABRJ   ', 'HET_I2O2_AI     ', & ! 1   
     &    'HET_I2O2_AJ     ', 'HET_I2O3_AI     ', 'HET_I2O3_AJ     ', & ! 2   
     &    'HET_I2O4_AI     ', 'HET_I2O4_AJ     ', 'HET_SO2_H2O2    ', & ! 3   
     &    'HET_SO2_O3      ', 'HET_SO2_MEPX    ', 'HET_SO2_PACD    ', & ! 4   
     &    'HET_SO2_NO2     ', 'HET_SO2_TMI     ', 'HET_SO2_HCHO    ', & ! 5   
     &    'HET_SO2_HNO4    ', 'HET_HMS_L1      ', 'HET_HMS_L2      '/! 6  

!    NSPECIAL     = Number of special rate coefficients
!    SPECIAL      = Names of special rate coefficients
!    NSPECIAL_RXN = Number of reactions with special rates
!    ISPECIAL     = Pointers to reactions using special rates and their special rate coefficients
!    MAXSPECTERMS = Max Number of terms type used by special rate coefficients
!    KC_COEFFS    = Coefficients of standard rate coefficients  times concentration terms 
!    INDEX_KTERMS  = Pointers to standard rate coefficients in  special rate coefficients
!    INDEX_CTERMS  = Pointers to species concentrations in  special rate coefficients
!    OPERATOR_COEFFS = Coefficients of preceeding special  rate coefficients used in special coefficient 
!    OPERATORS       = Pointers to preceeding special  rate coefficients used in special coefficient 

! Special Rate information not available ..
      INTEGER, PARAMETER :: NSPECIAL_RXN = 0
      INTEGER            :: ISPECIAL( 1, 2 )

! Special Rate information not available ...
      INTEGER, PARAMETER :: NSPECIAL = 0

! Special Rate information not available ...
      CHARACTER( 16 )    :: SPECIAL( 1 )

      INTEGER, PARAMETER :: MAXSPECTERMS =   1
      REAL( 8 )          :: KC_COEFFS( NSPECIAL + 1, MAXSPECTERMS)
      INTEGER            :: INDEX_KTERMS( NSPECIAL + 1, MAXSPECTERMS)
      INTEGER            :: INDEX_CTERMS( NSPECIAL + 1, MAXSPECTERMS)
      REAL( 8 )          :: OPERATOR_COEFFS( NSPECIAL + 1, MAXSPECTERMS)
      INTEGER            :: OPERATORS( NSPECIAL + 1, MAXSPECTERMS)


!    Steady-state species section
!    N_SS_SPC     = Number of species assumed to be in steady-state
!    SS_SPC_DIM   = Dimension paramete for steady-state species
!    SS_SPC       = Names of species assumed to be in steady-state
!    MAX_SS_LOSS  = Max no. of SS loss rxns for any SS species
!    MAX_SS_PROD  = Max no. of SS prod rxns for any SS species
!    N_LOSS_RXNS  = No. of SS loss rxns for each SS species
!    N_PROD_RXNS  = No. of SS prod rxns for each SS species
!    SS_LOSS_RXNS = List of SS loss rxns for each SS species
!    SS_PROD_RXNS = List of SS prod rxns for each SS species
!    SS_PROD_COEF = List of SS prod yields for each SS species
!    SS_RCT_IND   = SS species index if it is a rxn reactant

      INTEGER, PARAMETER :: N_SS_SPC =   0

      INTEGER, PARAMETER :: SS_SPC_DIM =   1

      INTEGER, PARAMETER :: MAX_SS_LOSS =   0

      INTEGER, PARAMETER :: MAX_SS_PROD =   0

      CHARACTER( 16 )    :: SS_SPC( 1 )

      INTEGER            :: N_LOSS_RXNS( 1 )
      INTEGER            :: N_PROD_RXNS( 1 )
      INTEGER            :: SS_LOSS_RXNS( 1, 1 )
      INTEGER            :: SS_PROD_RXNS( 1, 1 )
      INTEGER            :: SS_RCT_IND( 1 )

      REAL               :: SS_PROD_COEF( 1,1 ) 
       LOGICAL,  PARAMETER :: USE_SPECIAL_RATES = .FALSE.
! pointers and names to specific photolysis rates
       INTEGER, PARAMETER  :: IJ_O3O3P_NASA06     =   1
       INTEGER, PARAMETER  :: IJ_O3O1D_NASA06     =   2
       INTEGER, PARAMETER  :: IJ_H2O2_RACM2       =   3
       INTEGER, PARAMETER  :: IJ_NO2_RACM2        =   4
       INTEGER, PARAMETER  :: IJ_NO3NO_RACM2      =   5
       INTEGER, PARAMETER  :: IJ_NO3NO2_RACM2     =   6
       INTEGER, PARAMETER  :: IJ_HONO_RACM2       =   7
       INTEGER, PARAMETER  :: IJ_HNO3_RACM2       =   8
       INTEGER, PARAMETER  :: IJ_HNO4_RACM2       =   9
       INTEGER, PARAMETER  :: IJ_HCHO_MOL_JPL19   =  10
       INTEGER, PARAMETER  :: IJ_HCHO_RAD_JPL19   =  11
       INTEGER, PARAMETER  :: IJ_CH3CHO_RACM2     =  12
       INTEGER, PARAMETER  :: IJ_ALD_JPL19        =  13
       INTEGER, PARAMETER  :: IJ_CH3COCH3A_JPL19  =  14
       INTEGER, PARAMETER  :: IJ_CH3COCH3B_JPL19  =  15
       INTEGER, PARAMETER  :: IJ_UALD_RACM2       =  16
       INTEGER, PARAMETER  :: IJ_MEK_JGR19        =  17
       INTEGER, PARAMETER  :: IJ_KET_JGR19        =  18
       INTEGER, PARAMETER  :: IJ_HKET_RACM2       =  19
       INTEGER, PARAMETER  :: IJ_MACR_RACM2       =  20
       INTEGER, PARAMETER  :: IJ_MVK_JPL19        =  21
       INTEGER, PARAMETER  :: IJ_GLYH2_JPL19      =  22
       INTEGER, PARAMETER  :: IJ_GLYF_JPL19       =  23
       INTEGER, PARAMETER  :: IJ_GLYHX_JPL19      =  24
       INTEGER, PARAMETER  :: IJ_MGLY_RACM2       =  25
       INTEGER, PARAMETER  :: IJ_BALD1_CALVERT11  =  26
       INTEGER, PARAMETER  :: IJ_BALD2_CALVERT11  =  27
       INTEGER, PARAMETER  :: IJ_OP1_RACM2        =  28
       INTEGER, PARAMETER  :: IJ_PAA_RACM2        =  29
       INTEGER, PARAMETER  :: IJ_ONIT_CALVERT08   =  30
       INTEGER, PARAMETER  :: IJ_PAN1_JPL19       =  31
       INTEGER, PARAMETER  :: IJ_PAN2_JPL19       =  32
       INTEGER, PARAMETER  :: IJ_PPN1_JPL19       =  33
       INTEGER, PARAMETER  :: IJ_PPN2_JPL19       =  34
       INTEGER, PARAMETER  :: IJ_TRPN_WANG2023    =  35
       INTEGER, PARAMETER  :: IJ_ACRO_09          =  36
       INTEGER, PARAMETER  :: IJ_CL2_JPL19        =  37
       INTEGER, PARAMETER  :: IJ_CLO_JPL19        =  38
       INTEGER, PARAMETER  :: IJ_OCLO_JPL19       =  39
       INTEGER, PARAMETER  :: IJ_CL2O2_JPL19      =  40
       INTEGER, PARAMETER  :: IJ_HOCL_JPL19       =  41
       INTEGER, PARAMETER  :: IJ_CLNO_JPL19       =  42
       INTEGER, PARAMETER  :: IJ_CLNO2_JPL19      =  43
       INTEGER, PARAMETER  :: IJ_CLNO3_R_JPL19    =  44
       INTEGER, PARAMETER  :: IJ_CLNO3_M_JPL19    =  45
       INTEGER, PARAMETER  :: IJ_HCOCL_JPL19      =  46
       INTEGER, PARAMETER  :: IJ_BR2_JPL19        =  47
       INTEGER, PARAMETER  :: IJ_BRO_JPL19        =  48
       INTEGER, PARAMETER  :: IJ_OBRO_JPL19       =  49
       INTEGER, PARAMETER  :: IJ_HOBR_JPL19       =  50
       INTEGER, PARAMETER  :: IJ_BRNO_JPL19       =  51
       INTEGER, PARAMETER  :: IJ_BRNO2_JPL19      =  52
       INTEGER, PARAMETER  :: IJ_BRNO3_R_JPL19    =  53
       INTEGER, PARAMETER  :: IJ_BRNO3_M_JPL19    =  54
       INTEGER, PARAMETER  :: IJ_CH2BR2_JPL19     =  55
       INTEGER, PARAMETER  :: IJ_CHBR3_JPL19      =  56
       INTEGER, PARAMETER  :: IJ_HCOBR_JPL19      =  57
       INTEGER, PARAMETER  :: IJ_I2_JPL19         =  58
       INTEGER, PARAMETER  :: IJ_IO_JPL19         =  59
       INTEGER, PARAMETER  :: IJ_OIO_JPL19        =  60
       INTEGER, PARAMETER  :: IJ_INO3_06          =  61
       INTEGER, PARAMETER  :: IJ_HOI_JPL19        =  62
       INTEGER, PARAMETER  :: IJ_HI_JPL19         =  63
       INTEGER, PARAMETER  :: IJ_INO_JPL19        =  64
       INTEGER, PARAMETER  :: IJ_INO2_JPL19       =  65
       INTEGER, PARAMETER  :: IJ_CH3I_JPL19       =  66
       INTEGER, PARAMETER  :: IJ_CH2I2_JPL19      =  67
       INTEGER, PARAMETER  :: IJ_BRCL_JPL19       =  68
       INTEGER, PARAMETER  :: IJ_ICL_JPL19        =  69
       INTEGER, PARAMETER  :: IJ_IBR_JPL19        =  70
       INTEGER, PARAMETER  :: IJ_CH2IBR_JPL19     =  71
       INTEGER, PARAMETER  :: IJ_CH2ICL_JPL19     =  72
       INTEGER, PARAMETER  :: IJ_CHBR2CL_JPL19    =  73
       INTEGER, PARAMETER  :: IJ_CHBRCL2_JPL19    =  74
       INTEGER, PARAMETER  :: IK_HETERO_GLY       =   1
       INTEGER, PARAMETER  :: IK_HETERO_MGLY      =   2
       INTEGER, PARAMETER  :: IK_HETERO_NO2       =   3
       INTEGER, PARAMETER  :: IK_HETERO_HO2       =   4
       INTEGER, PARAMETER  :: IK_HETERO_NO3       =   5
       INTEGER, PARAMETER  :: IK_HETERO_IEPOX     =   6
       INTEGER, PARAMETER  :: IK_HETERO_ISO3NOSJ  =   7
       INTEGER, PARAMETER  :: IK_HETERO_ISO3OSJ   =   8
       INTEGER, PARAMETER  :: IK_HETERO_ANO3      =   9
       INTEGER, PARAMETER  :: IK_HETERO_N2O5IJ    =  10
       INTEGER, PARAMETER  :: IK_HETERO_N2O5K     =  11
       INTEGER, PARAMETER  :: IK_HETERO_H2NO3PAIJ =  12
       INTEGER, PARAMETER  :: IK_HETERO_H2NO3PAK  =  13
       INTEGER, PARAMETER  :: IK_HETERO_H2NO3PBIJ =  14
       INTEGER, PARAMETER  :: IK_HETERO_H2NO3PBK  =  15
       INTEGER, PARAMETER  :: IK_HETERO_HOCL_ACLJ =  16
       INTEGER, PARAMETER  :: IK_HETERO_HOCL_ABRJ =  17
       INTEGER, PARAMETER  :: IK_HETERO_CLN2_WAI  =  18
       INTEGER, PARAMETER  :: IK_HETERO_CLN2_WAJ  =  19
       INTEGER, PARAMETER  :: IK_HETERO_CLN2_ACLJ =  20
       INTEGER, PARAMETER  :: IK_HETERO_CLN2_ABRJ =  21
       INTEGER, PARAMETER  :: IK_HETERO_CLN3_WAI  =  22
       INTEGER, PARAMETER  :: IK_HETERO_CLN3_WAJ  =  23
       INTEGER, PARAMETER  :: IK_HETERO_CLN3_ACLJ =  24
       INTEGER, PARAMETER  :: IK_HETERO_CLN3_ABRJ =  25
       INTEGER, PARAMETER  :: IK_HETERO_HOBR_ACLJ =  26
       INTEGER, PARAMETER  :: IK_HETERO_HOBR_ABRJ =  27
       INTEGER, PARAMETER  :: IK_HETERO_BRN2_ACLJ =  28
       INTEGER, PARAMETER  :: IK_HETERO_BRN2_ABRJ =  29
       INTEGER, PARAMETER  :: IK_HETERO_BRN3_WAI  =  30
       INTEGER, PARAMETER  :: IK_HETERO_BRN3_WAJ  =  31
       INTEGER, PARAMETER  :: IK_HETERO_BRN3_ACLJ =  32
       INTEGER, PARAMETER  :: IK_HETERO_BRN3_ABRJ =  33
       INTEGER, PARAMETER  :: IK_HETERO_HBR_ABRJ  =  34
       INTEGER, PARAMETER  :: IK_HETERO_HOI_ACLJ  =  35
       INTEGER, PARAMETER  :: IK_HETERO_HOI_ABRJ  =  36
       INTEGER, PARAMETER  :: IK_HETERO_INO2_ACLJ =  37
       INTEGER, PARAMETER  :: IK_HETERO_INO2_ABRJ =  38
       INTEGER, PARAMETER  :: IK_HETERO_INO3_ACLJ =  39
       INTEGER, PARAMETER  :: IK_HETERO_INO3_ABRJ =  40
       INTEGER, PARAMETER  :: IK_HETERO_I2O2_AI   =  41
       INTEGER, PARAMETER  :: IK_HETERO_I2O2_AJ   =  42
       INTEGER, PARAMETER  :: IK_HETERO_I2O3_AI   =  43
       INTEGER, PARAMETER  :: IK_HETERO_I2O3_AJ   =  44
       INTEGER, PARAMETER  :: IK_HETERO_I2O4_AI   =  45
       INTEGER, PARAMETER  :: IK_HETERO_I2O4_AJ   =  46
       INTEGER, PARAMETER  :: IK_HETERO_SO2H      =  47
       INTEGER, PARAMETER  :: IK_HETERO_SO2O      =  48
       INTEGER, PARAMETER  :: IK_HETERO_SO2M      =  49
       INTEGER, PARAMETER  :: IK_HETERO_SO2P      =  50
       INTEGER, PARAMETER  :: IK_HETERO_SO2N      =  51
       INTEGER, PARAMETER  :: IK_HETERO_SO2T      =  52
       INTEGER, PARAMETER  :: IK_HETERO_HMSP      =  53
       INTEGER, PARAMETER  :: IK_HETERO_HNO4      =  54
       INTEGER, PARAMETER  :: IK_HETERO_HMSL1     =  55
       INTEGER, PARAMETER  :: IK_HETERO_HMSL2     =  56
       END MODULE RXNS_DATA
