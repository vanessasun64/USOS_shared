function J = CRACMM3M_J(Met,Jmethod)
% Calculates photolysis frequencies for the cracmm3m mechanism in the CMAQ model
% Met: structure containing required meteorological constraints. Required vars depend on Jmethod.
%       Met.SZA: solar zenith angle in degrees
%       Met.ALT: altitude, meters
%       Met.O3col: overhead ozone column, DU
%       Met.albedo: surface reflectance, 0-1 (unitless)
%       Met.T: temperature, T
%       Met.P: pressure, mbar
%       Met.LFlux: name of a text file containing an actinic flux spectrum
%
% Jmethod: numeric flag or string specifying how to calculate J-values. Default is 'MCM'.
%       0 or 'MCM':      use MCMv3.3.1 parameterization.
%                         Some reactions are not included in MCM. For these, 'HYBRID' values are used.
%                         Required Met fields: SZA
%       1 or 'BOTTOMUP': bottom-up integration of cross sections/quantum yields.
%                         See J_BottomUp.m for more info.
%                         Required Met fields: LFlux, T, P
%       2 or 'HYBRID':   Interpolation of hybrid J-values from TUV solar spectra.
%                         See J_TUVhybrid.m for more info.
%                         Required Met fields: SZA, ALT, O3col, albedo
%
% OUTPUTS:
% J: structure of J-values.
%
% INPUTS
struct2var(Met)

if nargin<2
    Jmethod = 'MCM';
elseif ischar(Jmethod)
    Jmethod = upper(Jmethod);
end

% J-Values
switch Jmethod
    case {0,'MCM'}
        error(['MCM option not functional for cracmm3m mechanism.'])

    case {1,'BOTTOMUP'}
        Jmcm = J_BottomUp(LFlux,T,P);

    case {2,'HYBRID'}
        Jmcm = J_Hybrid(SZA,ALT,O3col,albedo);

    otherwise
        fprintf('Jmethod = %f\n',Jmethod);
        error(['MCMv331_J: invalid Jmethod option selected'])

end
%rename
J=struct;
J.JO3O3P_NASA06 = Jmcm.J_O3O3P_NASA06;
J.JO3O1D_NASA06 = Jmcm.J_O3O1D_NASA06;
J.JH2O2_RACM2 = Jmcm.J_H2O2_RACM2;
J.JNO2_RACM2 = Jmcm.J_NO2_RACM2;
J.JNO3NO_RACM2 = Jmcm.J_NO3NO_RACM2;
J.JNO3NO2_RACM2 = Jmcm.J_NO3NO2_RACM2;
J.JHONO_RACM2 = Jmcm.J_HONO_RACM2;
J.JHNO3_RACM2 = Jmcm.J_HNO3_RACM2;
J.JHNO4_RACM2 = Jmcm.J_HNO4_RACM2;
J.JHCHO_MOL_JPL19 = Jmcm.J_HCHO_MOL_JPL19;
J.JHCHO_RAD_JPL19 = Jmcm.J_HCHO_RAD_JPL19;
J.JCH3CHO_RACM2 = Jmcm.J_CH3CHO_RACM2;
J.JALD_JPL19 = Jmcm.J_ALD_JPL19;
J.JCH3COCH3A_JPL19 = Jmcm.J_CH3COCH3A_JPL19;
J.JCH3COCH3B_JPL19 = Jmcm.J_CH3COCH3B_JPL19;
J.JUALD_RACM2 = Jmcm.J_UALD_RACM2;
J.JMEK_JGR19 = Jmcm.J_MEK_JGR19;
J.JKET_JGR19 = Jmcm.J_KET_JGR19;
J.JHKET_RACM2 = Jmcm.J_HKET_RACM2;
J.JMACR_RACM2 = Jmcm.J_MACR_RACM2;
J.JMVK_JPL19 = Jmcm.J_MVK_JPL19;
J.JGLYH2_JPL19 = Jmcm.J_GLYH2_JPL19;
J.JGLYF_JPL19 = Jmcm.J_GLYF_JPL19;
J.JGLYHX_JPL19 = Jmcm.J_GLYHX_JPL19;
J.JMGLY_RACM2 = Jmcm.J_MGLY_RACM2;
J.JBALD1_CALVERT11 = Jmcm.J_BALD1_CALVERT11;
J.JBALD2_CALVERT11 = Jmcm.J_BALD2_CALVERT11;
J.JOP1_RACM2 = Jmcm.J_OP1_RACM2;
J.JPAA_RACM2 = Jmcm.J_PAA_RACM2;
J.JONIT_CALVERT08 = Jmcm.J_ONIT_CALVERT08;
J.JPAN1_JPL19 = Jmcm.J_PAN1_JPL19;
J.JPAN2_JPL19 = Jmcm.J_PAN2_JPL19;
J.JPPN1_JPL19 = Jmcm.J_PPN1_JPL19;
J.JPPN2_JPL19 = Jmcm.J_PPN2_JPL19;
J.JTRPN_WANG2023 = Jmcm.J_TRPN_WANG2023;
J.JACRO_09 = Jmcm.J_ACRO_09;
J.JCL2_JPL19 = Jmcm.J_CL2_JPL19;
J.JCLO_JPL19 = Jmcm.J_CLO_JPL19;
J.JOCLO_JPL19 = Jmcm.J_OCLO_JPL19;
J.JCL2O2_JPL19 = Jmcm.J_CL2O2_JPL19;
J.JHOCL_JPL19 = Jmcm.J_HOCL_JPL19;
J.JCLNO_JPL19 = Jmcm.J_CLNO_JPL19;
J.JCLNO2_JPL19 = Jmcm.J_CLNO2_JPL19;
J.JCLNO3_R_JPL19 = Jmcm.J_CLNO3_R_JPL19;
J.JCLNO3_M_JPL19 = Jmcm.J_CLNO3_M_JPL19;
J.JHCOCL_JPL19 = Jmcm.J_HCOCL_JPL19;
J.JBR2_JPL19 = Jmcm.J_BR2_JPL19;
J.JBRO_JPL19 = Jmcm.J_BRO_JPL19;
J.JOBRO_JPL19 = Jmcm.J_OBRO_JPL19;
J.JHOBR_JPL19 = Jmcm.J_HOBR_JPL19;
J.JBRNO_JPL19 = Jmcm.J_BRNO_JPL19;
J.JBRNO2_JPL19 = Jmcm.J_BRNO2_JPL19;
J.JBRNO3_R_JPL19 = Jmcm.J_BRNO3_R_JPL19;
J.JBRNO3_M_JPL19 = Jmcm.J_BRNO3_M_JPL19;
J.JCH2BR2_JPL19 = Jmcm.J_CH2BR2_JPL19;
J.JCHBR3_JPL19 = Jmcm.J_CHBR3_JPL19;
J.JHCOBR_JPL19 = Jmcm.J_HCOBR_JPL19;
J.JI2_JPL19 = Jmcm.J_I2_JPL19;
J.JIO_JPL19 = Jmcm.J_IO_JPL19;
J.JOIO_JPL19 = Jmcm.J_OIO_JPL19;
J.JINO3_06 = Jmcm.J_INO3_06;
J.JHOI_JPL19 = Jmcm.J_HOI_JPL19;
J.JHI_JPL19 = Jmcm.J_HI_JPL19;
J.JINO_JPL19 = Jmcm.J_INO_JPL19;
J.JINO2_JPL19 = Jmcm.J_INO2_JPL19;
J.JCH3I_JPL19 = Jmcm.J_CH3I_JPL19;
J.JCH2I2_JPL19 = Jmcm.J_CH2I2_JPL19;
J.JBRCL_JPL19 = Jmcm.J_BRCL_JPL19;
J.JICL_JPL19 = Jmcm.J_ICL_JPL19;
J.JIBR_JPL19 = Jmcm.J_IBR_JPL19;
J.JCH2IBR_JPL19 = Jmcm.J_CH2IBR_JPL19;
J.JCH2ICL_JPL19 = Jmcm.J_CH2ICL_JPL19;
J.JCHBR2CL_JPL19 = Jmcm.J_CHBR2CL_JPL19;
J.JCHBRCL2_JPL19 = Jmcm.J_CHBRCL2_JPL19;
