
% # of species =125
% # of reactions =329


SpeciesToAdd = {'AACD';'ACET';'ALD2';'ALDX';'BENZ';'BR';'BR2';'BRCL';'BRN2';'BRN3';'BRO';'BZO2';...
'C2O3';'CAT1';'CGLY';'CH3I';'CH4';'CL';'CL2';'CLAD';'CLAO';'CLN2';'CLN3';'CLO';...
'CO';'CRES';'CRO';'CRON';'CXO3';'DMS';'ECH4';'EPOX';'EPX2';'ETH';'ETHA';'ETHY';...
'ETOH';'FACD';'FMBR';'FMCL';'FORM';'GLY';'GLYD';'H2';'H2O';'H2O2';'HBR';'HCL';...
'HCO3';'HI';'HIO3';'HNO3';'HO2';'HOBR';'HOCL';'HOI';'HONO';'HPLD';'I';'I2';'I2O2';...
'IBR';'ICL';'INO2';'INO3';'INTR';'IO';'IOLE';'ISO2';'ISOP';'ISPD';'ISPX';'IXOY';...
'KET';'MB2';'MB2C';'MB3';'MBC';'MBC2';'MEO2';'MEOH';'MEPX';'MGLY';'MI2';'MIB';...
'MIC';'N2O5';'NO';'NO2';'NO3';'NTR1';'NTR2';'O';'O1D';'O2';'O3';'O3P';'OH';'OIO';...
'OLE';'OPAN';'OPEN';'OPO3';'PACD';'PAN';'PANX';'PAR';'PNA';'PRPA';'ROOH';...
'ROR';'SO2';'SULF';'TERP';'TO2';'TOL';'XLO2';'XO2';'XO2H';'XO2N';'XOPN';'XPAR';
'XPRP'; 'XYL'};

RO2ToAdd = {'RO2'};

AddSpecies


%1
i=i+1;
Rnames{i} = 'NO2 = NO+ O3P';
k(:,i)=JNO2;
Gstr{i,1}='NO2';
fNO2(i)=fNO2(i)-1; fNO(i)=fNO(i)+1; fO3P(i)=fO3P(i)+1;

%2
i=i+1;
Rnames{i} = 'O3P = O3';
k(:,i)=6.00e-34.*((T./300).^-2.6).*M.^2.*0.21;
Gstr{i,1}='O3P';
fO3P(i)=fO3P(i)-1; fO3(i)=fO3(i)+1;

%3
i=i+1;
Rnames{i} = 'NO + O3 = NO2';
k(:,i)=2.07e-12.*exp(-1400./T);
Gstr{i,1}='NO'; Gstr{i,2}='O3';
fNO(i)=fNO(i)-1; fO3(i)=fO3(i)-1; fNO2(i)=fNO2(i)+1;

%4
i=i+1;
Rnames{i} = 'NO + O3P = NO2';
k(:,i)=K_O3P_NO;
Gstr{i,1}='NO'; Gstr{i,2}='O3P';
fNO(i)=fNO(i)-1; fO3P(i)=fO3P(i)-1; fNO2(i)=fNO2(i)+1;

%5
i=i+1;
Rnames{i} = 'O3P + NO2 = NO';
k(:,i)=5.10e-12.*exp(198./T);
Gstr{i,1}='O3P'; Gstr{i,2}='NO2';
fO3P(i)=fO3P(i)-1; fNO2(i)=fNO2(i)-1; fNO(i)=fNO(i)+1;

%6
i=i+1;
Rnames{i} = 'O3P + NO2 = NO3';
k(:,i)=K_O3P_NO2;
Gstr{i,1}='O3P'; Gstr{i,2}='NO2';
fO3P(i)=fO3P(i)-1; fNO2(i)=fNO2(i)-1; fNO3(i)=fNO3(i)+1;

%7
i=i+1;
Rnames{i} = 'O3P + O3 = ';
k(:,i)=8.00e-12.*exp(-2060./T);
Gstr{i,1}='O3P'; Gstr{i,2}='O3';
fO3P(i)=fO3P(i)-1; fO3(i)=fO3(i)-1;

%8
i=i+1;
Rnames{i} = 'O3 = O3P';
k(:,i)=JO3P;
Gstr{i,1}='O3';
fO3(i)=fO3(i)-1; fO3P(i)=fO3P(i)+1;

%9
i=i+1;
Rnames{i} = 'O3 = O1D';
k(:,i)=JO1D;
Gstr{i,1}='O3';
fO3(i)=fO3(i)-1; fO1D(i)=fO1D(i)+1;

%10
i=i+1;
Rnames{i} = 'O1D = O3P';
k(:,i)=2.23e-11.*exp(115./T).*M;
Gstr{i,1}='O1D';
fO1D(i)=fO1D(i)-1; fO3P(i)=fO3P(i)+1;

%11
i=i+1;
Rnames{i} = 'O1D = OH+ OH';
k(:,i)=2.140e-10.*H2O;
Gstr{i,1}='O1D';
fO1D(i)=fO1D(i)-1; fOH(i)=fOH(i)+1; fOH(i)=fOH(i)+1;

%12
i=i+1;
Rnames{i} = 'O3 + OH = HO2';
k(:,i)=1.70e-12.*exp(-940./T);
Gstr{i,1}='O3'; Gstr{i,2}='OH';
fO3(i)=fO3(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+1;

%13
i=i+1;
Rnames{i} = 'O3 + HO2 = OH';
k(:,i)=2.03e-16.*((T./300).^4.57).*exp(693./T);
Gstr{i,1}='O3'; Gstr{i,2}='HO2';
fO3(i)=fO3(i)-1; fHO2(i)=fHO2(i)-1; fOH(i)=fOH(i)+1;

%14
i=i+1;
Rnames{i} = 'O3P + OH = HO2';
k(:,i)=2.40e-11.*exp(110./T);
Gstr{i,1}='O3P'; Gstr{i,2}='OH';
fO3P(i)=fO3P(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+1;

%15
i=i+1;
Rnames{i} = 'O3P + HO2 = OH';
k(:,i)=3.00e-11.*exp(200./T);
Gstr{i,1}='O3P'; Gstr{i,2}='HO2';
fO3P(i)=fO3P(i)-1; fHO2(i)=fHO2(i)-1; fOH(i)=fOH(i)+1;

%16
i=i+1;
Rnames{i} = 'OH + OH = O3P';
k(:,i)=6.20e-14.*((T./298).^2.6).*exp(945./T);
Gstr{i,1}='OH'; Gstr{i,2}='OH';
fOH(i)=fOH(i)-1; fOH(i)=fOH(i)-1; fO3P(i)=fO3P(i)+1;

%17
i=i+1;
Rnames{i} = 'OH + OH = H2O2';
k(:,i)=K_OH_OH;
Gstr{i,1}='OH'; Gstr{i,2}='OH';
fOH(i)=fOH(i)-1; fOH(i)=fOH(i)-1; fH2O2(i)=fH2O2(i)+1;

%18
i=i+1;
Rnames{i} = 'OH + HO2 = ';
k(:,i)=4.80e-11.*exp(250./T);
Gstr{i,1}='OH'; Gstr{i,2}='HO2';
fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)-1;

%19
i=i+1;
Rnames{i} = 'HO2 + HO2 = H2O2';
k(:,i)=2.20e-13.*exp(600./T)+1.90e-33.*exp(980./T).*M;
Gstr{i,1}='HO2'; Gstr{i,2}='HO2';
fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)-1; fH2O2(i)=fH2O2(i)+1;

%20
i=i+1;
Rnames{i} = 'HO2 + HO2 = H2O2';
k(:,i)=3.08e-34.*exp(2800./T)+2.66e-54.*exp(3180./T).*M.*H2O;
Gstr{i,1}='HO2'; Gstr{i,2}='HO2';
fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)-1; fH2O2(i)=fH2O2(i)+1;

%21
i=i+1;
Rnames{i} = 'H2O2 = OH+ OH';
k(:,i)=JH2O2;
Gstr{i,1}='H2O2';
fH2O2(i)=fH2O2(i)-1; fOH(i)=fOH(i)+1; fOH(i)=fOH(i)+1;

%22
i=i+1;
Rnames{i} = 'H2O2 + OH = HO2';
k(:,i)=1.8E-12;
Gstr{i,1}='H2O2'; Gstr{i,2}='OH';
fH2O2(i)=fH2O2(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+1;

%23
i=i+1;
Rnames{i} = 'H2O2 + O3P = OH+ HO2';
k(:,i)=1.40e-12.*exp(-2000./T);
Gstr{i,1}='H2O2'; Gstr{i,2}='O3P';
fH2O2(i)=fH2O2(i)-1; fO3P(i)=fO3P(i)-1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1;

%24
i=i+1;
Rnames{i} = 'NO + NO = NO2+ NO2';
k(:,i)=4.25e-39.*exp(664./T).*M.*0.21;
Gstr{i,1}='NO'; Gstr{i,2}='NO';
fNO(i)=fNO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fNO2(i)=fNO2(i)+1;

%25
i=i+1;
Rnames{i} = 'NO + HO2 = NO2+ OH';
k(:,i)=3.45e-12.*exp(270./T);
Gstr{i,1}='NO'; Gstr{i,2}='HO2';
fNO(i)=fNO(i)-1; fHO2(i)=fHO2(i)-1; fNO2(i)=fNO2(i)+1; fOH(i)=fOH(i)+1;

%26
i=i+1;
Rnames{i} = 'NO2 + O3 = NO3';
k(:,i)=1.40e-13.*exp(-2470./T);
Gstr{i,1}='NO2'; Gstr{i,2}='O3';
fNO2(i)=fNO2(i)-1; fO3(i)=fO3(i)-1; fNO3(i)=fNO3(i)+1;

%27
i=i+1;
Rnames{i} = 'NO3 = O3P+ NO2';
k(:,i)=JNO3_NO2;
Gstr{i,1}='NO3';
fNO3(i)=fNO3(i)-1; fO3P(i)=fO3P(i)+1; fNO2(i)=fNO2(i)+1;

%28
i=i+1;
Rnames{i} = 'NO3 = NO';
k(:,i)=JNO3_NO;
Gstr{i,1}='NO3';
fNO3(i)=fNO3(i)-1; fNO(i)=fNO(i)+1;

%29
i=i+1;
Rnames{i} = 'NO + NO3 = NO2+ NO2';
k(:,i)=1.80e-11.*exp(110./T);
Gstr{i,1}='NO'; Gstr{i,2}='NO3';
fNO(i)=fNO(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fNO2(i)=fNO2(i)+1;

%30
i=i+1;
Rnames{i} = 'NO2 + NO3 = NO+ NO2';
k(:,i)=4.50e-14.*exp(-1260./T);
Gstr{i,1}='NO2'; Gstr{i,2}='NO3';
fNO2(i)=fNO2(i)-1; fNO3(i)=fNO3(i)-1; fNO(i)=fNO(i)+1; fNO2(i)=fNO2(i)+1;

%31
i=i+1;
Rnames{i} = 'O3P + NO3 = NO2';
k(:,i)=1.70e-11;
Gstr{i,1}='O3P'; Gstr{i,2}='NO3';
fO3P(i)=fO3P(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1;

%32
i=i+1;
Rnames{i} = 'NO3 + OH = NO2+ HO2';
k(:,i)=2.0e-11;
Gstr{i,1}='NO3'; Gstr{i,2}='OH';
fNO3(i)=fNO3(i)-1; fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1;

%33
i=i+1;
Rnames{i} = 'NO3 + HO2 = NO2+ OH';
k(:,i)=4.0e-12;
Gstr{i,1}='NO3'; Gstr{i,2}='HO2';
fNO3(i)=fNO3(i)-1; fHO2(i)=fHO2(i)-1; fNO2(i)=fNO2(i)+1; fOH(i)=fOH(i)+1;

%34
i=i+1;
Rnames{i} = 'O3 + NO3 = NO2';
k(:,i)=1.0e-17;
Gstr{i,1}='O3'; Gstr{i,2}='NO3';
fO3(i)=fO3(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1;

%35
i=i+1;
Rnames{i} = 'NO3 + NO3 = NO2+ NO2';
k(:,i)=8.50e-13.*exp(-2450./T);
Gstr{i,1}='NO3'; Gstr{i,2}='NO3';
fNO3(i)=fNO3(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fNO2(i)=fNO2(i)+1;

%36
i=i+1;
Rnames{i} = 'NO2 + NO3 = N2O5';
k(:,i)=K_NO3_NO2;
Gstr{i,1}='NO2'; Gstr{i,2}='NO3';
fNO2(i)=fNO2(i)-1; fNO3(i)=fNO3(i)-1; fN2O5(i)=fN2O5(i)+1;

%37
i=i+1;
Rnames{i} = 'N2O5 = NO2+ NO3';
k(:,i)=K_N2O5;
Gstr{i,1}='N2O5';
fN2O5(i)=fN2O5(i)-1; fNO2(i)=fNO2(i)+1; fNO3(i)=fNO3(i)+1;

%38
i=i+1;
Rnames{i} = 'N2O5 = NO2+ NO3';
k(:,i)=JN2O5;
Gstr{i,1}='N2O5';
fN2O5(i)=fN2O5(i)-1; fNO2(i)=fNO2(i)+1; fNO3(i)=fNO3(i)+1;

%39
i=i+1;
Rnames{i} = 'N2O5 = HNO3+ HNO3';
k(:,i)=1.0e-22.*H2O;
Gstr{i,1}='N2O5';
fN2O5(i)=fN2O5(i)-1; fHNO3(i)=fHNO3(i)+1; fHNO3(i)=fHNO3(i)+1;

%40
i=i+1;
Rnames{i} = 'NO + OH = HONO';
k(:,i)=K_OH_NO;
Gstr{i,1}='NO'; Gstr{i,2}='OH';
fNO(i)=fNO(i)-1; fOH(i)=fOH(i)-1; fHONO(i)=fHONO(i)+1;

%41
i=i+1;
Rnames{i} = 'NO + NO2 = HONO+ HONO';
k(:,i)=5e-40.*H2O;
Gstr{i,1}='NO'; Gstr{i,2}='NO2';
fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)-1; fHONO(i)=fHONO(i)+1; fHONO(i)=fHONO(i)+1;

%42
i=i+1;
Rnames{i} = 'HONO + HONO = NO+ NO2';
k(:,i)=1e-20;
Gstr{i,1}='HONO'; Gstr{i,2}='HONO';
fHONO(i)=fHONO(i)-1; fHONO(i)=fHONO(i)-1; fNO(i)=fNO(i)+1; fNO2(i)=fNO2(i)+1;

%43
i=i+1;
Rnames{i} = 'HONO = NO+ OH';
k(:,i)=JHONO;
Gstr{i,1}='HONO';
fHONO(i)=fHONO(i)-1; fNO(i)=fNO(i)+1; fOH(i)=fOH(i)+1;

%44
i=i+1;
Rnames{i} = 'HONO + OH = NO2';
k(:,i)=2.50e-12.*exp(260./T);
Gstr{i,1}='HONO'; Gstr{i,2}='OH';
fHONO(i)=fHONO(i)-1; fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)+1;

%45
i=i+1;
Rnames{i} = 'NO2 + OH = HNO3';
k(:,i)=K_OH_NO2;
Gstr{i,1}='NO2'; Gstr{i,2}='OH';
fNO2(i)=fNO2(i)-1; fOH(i)=fOH(i)-1; fHNO3(i)=fHNO3(i)+1;

%46
i=i+1;
Rnames{i} = 'HNO3 + OH = NO3';
k(:,i)=K_OH_HNO3;
Gstr{i,1}='HNO3'; Gstr{i,2}='OH';
fHNO3(i)=fHNO3(i)-1; fOH(i)=fOH(i)-1; fNO3(i)=fNO3(i)+1;

%47
i=i+1;
Rnames{i} = 'HNO3 = NO2+ OH';
k(:,i)=JHNO3;
Gstr{i,1}='HNO3';
fHNO3(i)=fHNO3(i)-1; fNO2(i)=fNO2(i)+1; fOH(i)=fOH(i)+1;

%48
i=i+1;
Rnames{i} = 'NO2 + HO2 = PNA';
k(:,i)=K_HO2_NO2;
Gstr{i,1}='NO2'; Gstr{i,2}='HO2';
fNO2(i)=fNO2(i)-1; fHO2(i)=fHO2(i)-1; fPNA(i)=fPNA(i)+1;

%49
i=i+1;
Rnames{i} = 'PNA = NO2+ HO2';
k(:,i)=K_PNA;
Gstr{i,1}='PNA';
fPNA(i)=fPNA(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1;

%50
i=i+1;
Rnames{i} = 'PNA = 0.59NO2+ 0.41NO3+ 0.41OH+ 0.59HO2';
k(:,i)=JPNA;
Gstr{i,1}='PNA';
fPNA(i)=fPNA(i)-1; fNO2(i)=fNO2(i)+0.59; fNO3(i)=fNO3(i)+0.41; fOH(i)=fOH(i)+0.41; fHO2(i)=fHO2(i)+0.59;

%51
i=i+1;
Rnames{i} = 'PNA + OH = NO2';
k(:,i)=3.20e-13.*exp(690./T);
Gstr{i,1}='PNA'; Gstr{i,2}='OH';
fPNA(i)=fPNA(i)-1; fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)+1;

%52
i=i+1;
Rnames{i} = 'SO2 + OH = SULF+ HO2';
k(:,i)=K_SO2_OH;
Gstr{i,1}='SO2'; Gstr{i,2}='OH';
fSO2(i)=fSO2(i)-1; fOH(i)=fOH(i)-1; fSULF(i)=fSULF(i)+1; fHO2(i)=fHO2(i)+1;

%53
i=i+1;
Rnames{i} = 'C2O3 + NO = MEO2+ NO2+ RO2';
k(:,i)=7.50e-12.*exp(290./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='NO';
fC2O3(i)=fC2O3(i)-1; fNO(i)=fNO(i)-1; fMEO2(i)=fMEO2(i)+1; fNO2(i)=fNO2(i)+1; fRO2(i)=fRO2(i)+1;

%54
i=i+1;
Rnames{i} = 'C2O3 + NO2 = PAN';
k(:,i)=K_C2O3_NO2;
Gstr{i,1}='C2O3'; Gstr{i,2}='NO2';
fC2O3(i)=fC2O3(i)-1; fNO2(i)=fNO2(i)-1; fPAN(i)=fPAN(i)+1;

%55
i=i+1;
Rnames{i} = 'PAN = C2O3+ NO2';
k(:,i)=K_PAN;
Gstr{i,1}='PAN';
fPAN(i)=fPAN(i)-1; fC2O3(i)=fC2O3(i)+1; fNO2(i)=fNO2(i)+1;

%56
i=i+1;
Rnames{i} = 'PAN = 0.6C2O3+ 0.4MEO2+ 0.6NO2+ 0.4NO3+ 0.4RO2';
k(:,i)=JPAN;
Gstr{i,1}='PAN';
fPAN(i)=fPAN(i)-1; fC2O3(i)=fC2O3(i)+0.6; fMEO2(i)=fMEO2(i)+0.4; fNO2(i)=fNO2(i)+0.6; fNO3(i)=fNO3(i)+0.4; fRO2(i)=fRO2(i)+0.4;

%57
i=i+1;
Rnames{i} = 'C2O3 + HO2 = 0.13AACD+ 0.5MEO2+ 0.37PACD+ 0.13O3+ 0.5RO2+ 0.5OH';
k(:,i)=3.14e-12.*exp(580./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='HO2';
fC2O3(i)=fC2O3(i)-1; fHO2(i)=fHO2(i)-1; fAACD(i)=fAACD(i)+0.13; fMEO2(i)=fMEO2(i)+0.5; fPACD(i)=fPACD(i)+0.37; fO3(i)=fO3(i)+0.13; fRO2(i)=fRO2(i)+0.5; fOH(i)=fOH(i)+0.5;

%58
i=i+1;
Rnames{i} = 'C2O3 + RO2 = MEO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='RO2';
fC2O3(i)=fC2O3(i)-1; fRO2(i)=fRO2(i)-1; fMEO2(i)=fMEO2(i)+1;

%59
i=i+1;
Rnames{i} = 'C2O3 + C2O3 = MEO2+ MEO2+ RO2+ RO2';
k(:,i)=2.90e-12.*exp(500./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='C2O3';
fC2O3(i)=fC2O3(i)-1; fC2O3(i)=fC2O3(i)-1; fMEO2(i)=fMEO2(i)+1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1; fRO2(i)=fRO2(i)+1;

%60
i=i+1;
Rnames{i} = 'C2O3 + CXO3 = ALD2+ MEO2+ RO2+ XO2H+ RO2';
k(:,i)=2.90e-12.*exp(500./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='CXO3';
fC2O3(i)=fC2O3(i)-1; fCXO3(i)=fCXO3(i)-1; fALD2(i)=fALD2(i)+1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%61
i=i+1;
Rnames{i} = 'CXO3 + NO = ALD2+ XO2H+ NO2+ RO2';
k(:,i)=6.70e-12.*exp(340./T);
Gstr{i,1}='CXO3'; Gstr{i,2}='NO';
fCXO3(i)=fCXO3(i)-1; fNO(i)=fNO(i)-1; fALD2(i)=fALD2(i)+1; fXO2H(i)=fXO2H(i)+1; fNO2(i)=fNO2(i)+1; fRO2(i)=fRO2(i)+1;

%62
i=i+1;
Rnames{i} = 'CXO3 + NO2 = PANX';
k(:,i)=K_C2O3_NO2./1.19e0;
Gstr{i,1}='CXO3'; Gstr{i,2}='NO2';
fCXO3(i)=fCXO3(i)-1; fNO2(i)=fNO2(i)-1; fPANX(i)=fPANX(i)+1;

%63
i=i+1;
Rnames{i} = 'PANX = CXO3+ NO2';
k(:,i)=K_PAN./1.19e0;
Gstr{i,1}='PANX';
fPANX(i)=fPANX(i)-1; fCXO3(i)=fCXO3(i)+1; fNO2(i)=fNO2(i)+1;

%64
i=i+1;
Rnames{i} = 'PANX = 0.4ALD2+ 0.6CXO3+ 0.4XO2H+ 0.6NO2+ 0.4NO3+ 0.4RO2';
k(:,i)=JPAN;
Gstr{i,1}='PANX';
fPANX(i)=fPANX(i)-1; fALD2(i)=fALD2(i)+0.4; fCXO3(i)=fCXO3(i)+0.6; fXO2H(i)=fXO2H(i)+0.4; fNO2(i)=fNO2(i)+0.6; fNO3(i)=fNO3(i)+0.4; fRO2(i)=fRO2(i)+0.4;

%65
i=i+1;
Rnames{i} = 'CXO3 + HO2 = 0.13AACD+ 0.5MEO2+ 0.37PACD+ 0.13O3+ 0.5RO2+ 0.5OH';
k(:,i)=3.14e-12.*exp(580./T);
Gstr{i,1}='CXO3'; Gstr{i,2}='HO2';
fCXO3(i)=fCXO3(i)-1; fHO2(i)=fHO2(i)-1; fAACD(i)=fAACD(i)+0.13; fMEO2(i)=fMEO2(i)+0.5; fPACD(i)=fPACD(i)+0.37; fO3(i)=fO3(i)+0.13; fRO2(i)=fRO2(i)+0.5; fOH(i)=fOH(i)+0.5;

%66
i=i+1;
Rnames{i} = 'CXO3 + RO2 = MEO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='CXO3'; Gstr{i,2}='RO2';
fCXO3(i)=fCXO3(i)-1; fRO2(i)=fRO2(i)-1; fMEO2(i)=fMEO2(i)+1;

%67
i=i+1;
Rnames{i} = 'CXO3 + CXO3 = MEO2+ MEO2+ RO2+ RO2';
k(:,i)=2.90e-12.*exp(500./T);
Gstr{i,1}='CXO3'; Gstr{i,2}='CXO3';
fCXO3(i)=fCXO3(i)-1; fCXO3(i)=fCXO3(i)-1; fMEO2(i)=fMEO2(i)+1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1; fRO2(i)=fRO2(i)+1;

%68
i=i+1;
Rnames{i} = 'NO + RO2 = NO';
k(:,i)=2.40e-12.*exp(360./T);
Gstr{i,1}='NO'; Gstr{i,2}='RO2';
fNO(i)=fNO(i)-1; fRO2(i)=fRO2(i)-1; fNO(i)=fNO(i)+1;

%69
i=i+1;
Rnames{i} = 'RO2 + HO2 = HO2';
k(:,i)=4.80e-13.*exp(800./T) ;
Gstr{i,1}='RO2'; Gstr{i,2}='HO2';
fRO2(i)=fRO2(i)-1; fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)+1;

%70
i=i+1;
Rnames{i} = 'RO2 + RO2 = ';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='RO2'; Gstr{i,2}='RO2';
fRO2(i)=fRO2(i)-1; fRO2(i)=fRO2(i)-1;

%71
i=i+1;
Rnames{i} = 'MEO2 + NO = FORM+ NO2+ HO2';
k(:,i)=2.30e-12.*exp(360./T) ;
Gstr{i,1}='MEO2'; Gstr{i,2}='NO';
fMEO2(i)=fMEO2(i)-1; fNO(i)=fNO(i)-1; fFORM(i)=fFORM(i)+1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1;

%72
i=i+1;
Rnames{i} = 'MEO2 + HO2 = 0.1FORM+ 0.9MEPX';
k(:,i)=3.80e-13.*exp(780./T) ;
Gstr{i,1}='MEO2'; Gstr{i,2}='HO2';
fMEO2(i)=fMEO2(i)-1; fHO2(i)=fHO2(i)-1; fFORM(i)=fFORM(i)+0.1; fMEPX(i)=fMEPX(i)+0.9;

%73
i=i+1;
Rnames{i} = 'C2O3 + MEO2 = 0.1AACD+ FORM+ 0.9MEO2+ 0.9RO2+ 0.9HO2';
k(:,i)=2.00e-12.*exp(500./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='MEO2';
fC2O3(i)=fC2O3(i)-1; fMEO2(i)=fMEO2(i)-1; fAACD(i)=fAACD(i)+0.1; fFORM(i)=fFORM(i)+1; fMEO2(i)=fMEO2(i)+0.9; fRO2(i)=fRO2(i)+0.9; fHO2(i)=fHO2(i)+0.9;

%74
i=i+1;
Rnames{i} = 'MEO2 + RO2 = 0.685FORM+ 0.315MEOH+ RO2+ 0.37HO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='MEO2'; Gstr{i,2}='RO2';
fMEO2(i)=fMEO2(i)-1; fRO2(i)=fRO2(i)-1; fFORM(i)=fFORM(i)+0.685; fMEOH(i)=fMEOH(i)+0.315; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+0.37;

%75
i=i+1;
Rnames{i} = 'NO + XO2H = NO2+ HO2';
k(:,i)=2.70e-12.*exp(360./T);
Gstr{i,1}='NO'; Gstr{i,2}='XO2H';
fNO(i)=fNO(i)-1; fXO2H(i)=fXO2H(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1;

%76
i=i+1;
Rnames{i} = 'XO2H + HO2 = ROOH';
k(:,i)=6.80e-13.*exp(800./T);
Gstr{i,1}='XO2H'; Gstr{i,2}='HO2';
fXO2H(i)=fXO2H(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1;

%77
i=i+1;
Rnames{i} = 'C2O3 + XO2H = 0.2AACD+ 0.8MEO2+ 0.8RO2+ 0.8HO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='XO2H';
fC2O3(i)=fC2O3(i)-1; fXO2H(i)=fXO2H(i)-1; fAACD(i)=fAACD(i)+0.2; fMEO2(i)=fMEO2(i)+0.8; fRO2(i)=fRO2(i)+0.8; fHO2(i)=fHO2(i)+0.8;

%78
i=i+1;
Rnames{i} = 'XO2H + RO2 = RO2+ 0.6HO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='XO2H'; Gstr{i,2}='RO2';
fXO2H(i)=fXO2H(i)-1; fRO2(i)=fRO2(i)-1; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+0.6;

%79
i=i+1;
Rnames{i} = 'NO + XO2 = NO2';
k(:,i)=2.70e-12.*exp(360./T);
Gstr{i,1}='NO'; Gstr{i,2}='XO2';
fNO(i)=fNO(i)-1; fXO2(i)=fXO2(i)-1; fNO2(i)=fNO2(i)+1;

%80
i=i+1;
Rnames{i} = 'XO2 + HO2 = ROOH';
k(:,i)=6.80e-13.*exp(800./T);
Gstr{i,1}='XO2'; Gstr{i,2}='HO2';
fXO2(i)=fXO2(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1;

%81
i=i+1;
Rnames{i} = 'C2O3 + XO2 = 0.2AACD+ 0.8MEO2+ 0.8RO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='XO2';
fC2O3(i)=fC2O3(i)-1; fXO2(i)=fXO2(i)-1; fAACD(i)=fAACD(i)+0.2; fMEO2(i)=fMEO2(i)+0.8; fRO2(i)=fRO2(i)+0.8;

%82
i=i+1;
Rnames{i} = 'XO2 + RO2 = RO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='XO2'; Gstr{i,2}='RO2';
fXO2(i)=fXO2(i)-1; fRO2(i)=fRO2(i)-1; fRO2(i)=fRO2(i)+1;

%83
i=i+1;
Rnames{i} = 'NO + XO2N = 0.5NTR1+ 0.5NTR2';
k(:,i)=2.70e-12.*exp(360./T);
Gstr{i,1}='NO'; Gstr{i,2}='XO2N';
fNO(i)=fNO(i)-1; fXO2N(i)=fXO2N(i)-1; fNTR1(i)=fNTR1(i)+0.5; fNTR2(i)=fNTR2(i)+0.5;

%84
i=i+1;
Rnames{i} = 'XO2N + HO2 = ROOH';
k(:,i)=6.80e-13.*exp(800./T);
Gstr{i,1}='XO2N'; Gstr{i,2}='HO2';
fXO2N(i)=fXO2N(i)-1; fHO2(i)=fHO2(i)-1; fROOH(i)=fROOH(i)+1;

%85
i=i+1;
Rnames{i} = 'C2O3 + XO2N = 0.2AACD+ 0.8MEO2+ 0.8RO2+ 0.8HO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='XO2N';
fC2O3(i)=fC2O3(i)-1; fXO2N(i)=fXO2N(i)-1; fAACD(i)=fAACD(i)+0.2; fMEO2(i)=fMEO2(i)+0.8; fRO2(i)=fRO2(i)+0.8; fHO2(i)=fHO2(i)+0.8;

%86
i=i+1;
Rnames{i} = 'XO2N + RO2 = RO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='XO2N'; Gstr{i,2}='RO2';
fXO2N(i)=fXO2N(i)-1; fRO2(i)=fRO2(i)-1; fRO2(i)=fRO2(i)+1;

%87
i=i+1;
Rnames{i} = 'MEPX + OH = 0.4FORM+ 0.6MEO2+ 0.6RO2+ 0.4OH';
k(:,i)=5.30e-12.*exp(190./T) ;
Gstr{i,1}='MEPX'; Gstr{i,2}='OH';
fMEPX(i)=fMEPX(i)-1; fOH(i)=fOH(i)-1; fFORM(i)=fFORM(i)+0.4; fMEO2(i)=fMEO2(i)+0.6; fRO2(i)=fRO2(i)+0.6; fOH(i)=fOH(i)+0.4;

%88
i=i+1;
Rnames{i} = 'MEPX = MEO2+ RO2+ OH';
k(:,i)=JCOOH;
Gstr{i,1}='MEPX';
fMEPX(i)=fMEPX(i)-1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1; fOH(i)=fOH(i)+1;

%89
i=i+1;
Rnames{i} = 'ROOH + OH = 0.54XO2H+ 0.06XO2N+ 0.6RO2+ 0.4OH';
k(:,i)=5.30e-12.*exp(190./T);
Gstr{i,1}='ROOH'; Gstr{i,2}='OH';
fROOH(i)=fROOH(i)-1; fOH(i)=fOH(i)-1; fXO2H(i)=fXO2H(i)+0.54; fXO2N(i)=fXO2N(i)+0.06; fRO2(i)=fRO2(i)+0.6; fOH(i)=fOH(i)+0.4;

%90
i=i+1;
Rnames{i} = 'ROOH = OH+ HO2';
k(:,i)=JCOOH;
Gstr{i,1}='ROOH';
fROOH(i)=fROOH(i)-1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1;

%91
i=i+1;
Rnames{i} = 'NTR1 + OH = NTR2';
k(:,i)=2e-12;
Gstr{i,1}='NTR1'; Gstr{i,2}='OH';
fNTR1(i)=fNTR1(i)-1; fOH(i)=fOH(i)-1; fNTR2(i)=fNTR2(i)+1;

%92
i=i+1;
Rnames{i} = 'NTR1 = NO2';
k(:,i)=JNTR;
Gstr{i,1}='NTR1';
fNTR1(i)=fNTR1(i)-1; fNO2(i)=fNO2(i)+1;

%93
i=i+1;
Rnames{i} = 'FACD + OH = HO2';
k(:,i)=4.5e-13;
Gstr{i,1}='FACD'; Gstr{i,2}='OH';
fFACD(i)=fFACD(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+1;

%94
i=i+1;
Rnames{i} = 'AACD + OH = MEO2+ RO2';
k(:,i)=4.00e-14.*exp(850./T);
Gstr{i,1}='AACD'; Gstr{i,2}='OH';
fAACD(i)=fAACD(i)-1; fOH(i)=fOH(i)-1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1;

%95
i=i+1;
Rnames{i} = 'PACD + OH = C2O3';
k(:,i)=1.00e-14;
Gstr{i,1}='PACD'; Gstr{i,2}='OH';
fPACD(i)=fPACD(i)-1; fOH(i)=fOH(i)-1; fC2O3(i)=fC2O3(i)+1;

%96
i=i+1;
Rnames{i} = 'FORM + OH = CO+ HO2';
k(:,i)=5.40e-12.*exp(135./T);
Gstr{i,1}='FORM'; Gstr{i,2}='OH';
fFORM(i)=fFORM(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%97
i=i+1;
Rnames{i} = 'FORM = HO2+ CO+ HO2';
k(:,i)=JHCHO_R;
Gstr{i,1}='FORM';
fFORM(i)=fFORM(i)-1; fHO2(i)=fHO2(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%98
i=i+1;
Rnames{i} = 'FORM = CO+ H2';
k(:,i)=JHCHO_M;
Gstr{i,1}='FORM';
fFORM(i)=fFORM(i)-1; fCO(i)=fCO(i)+1; fH2(i)=fH2(i)+1;

%99
i=i+1;
Rnames{i} = 'FORM + NO3 = HNO3+ CO+ HO2';
k(:,i)=5.50e-16;
Gstr{i,1}='FORM'; Gstr{i,2}='NO3';
fFORM(i)=fFORM(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%100
i=i+1;
Rnames{i} = 'FORM + HO2 = HCO3';
k(:,i)=9.70e-15.*exp(625./T);
Gstr{i,1}='FORM'; Gstr{i,2}='HO2';
fFORM(i)=fFORM(i)-1; fHO2(i)=fHO2(i)-1; fHCO3(i)=fHCO3(i)+1;

%101
i=i+1;
Rnames{i} = 'HCO3 = FORM+ HO2';
k(:,i)=2.40e+12.*exp(-7000./T) ;
Gstr{i,1}='HCO3';
fHCO3(i)=fHCO3(i)-1; fFORM(i)=fFORM(i)+1; fHO2(i)=fHO2(i)+1;

%102
i=i+1;
Rnames{i} = 'HCO3 + NO = FACD+ NO2+ HO2';
k(:,i)=5.60e-12;
Gstr{i,1}='HCO3'; Gstr{i,2}='NO';
fHCO3(i)=fHCO3(i)-1; fNO(i)=fNO(i)-1; fFACD(i)=fFACD(i)+1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1;

%103
i=i+1;
Rnames{i} = 'HCO3 + HO2 = 0.5FACD+ 0.5MEPX+ 0.2OH+ 0.2HO2';
k(:,i)=5.60e-15.*exp(2300./T);
Gstr{i,1}='HCO3'; Gstr{i,2}='HO2';
fHCO3(i)=fHCO3(i)-1; fHO2(i)=fHO2(i)-1; fFACD(i)=fFACD(i)+0.5; fMEPX(i)=fMEPX(i)+0.5; fOH(i)=fOH(i)+0.2; fHO2(i)=fHO2(i)+0.2;

%104
i=i+1;
Rnames{i} = 'ALD2 + OH = C2O3';
k(:,i)=4.70e-12.*exp(345./T);
Gstr{i,1}='ALD2'; Gstr{i,2}='OH';
fALD2(i)=fALD2(i)-1; fOH(i)=fOH(i)-1; fC2O3(i)=fC2O3(i)+1;

%105
i=i+1;
Rnames{i} = 'ALD2 + NO3 = C2O3+ HNO3';
k(:,i)=1.40e-12.*exp(-1860./T);
Gstr{i,1}='ALD2'; Gstr{i,2}='NO3';
fALD2(i)=fALD2(i)-1; fNO3(i)=fNO3(i)-1; fC2O3(i)=fC2O3(i)+1; fHNO3(i)=fHNO3(i)+1;

%106
i=i+1;
Rnames{i} = 'ALD2 = MEO2+ RO2+ CO+ HO2';
k(:,i)=JCCHO_R;
Gstr{i,1}='ALD2';
fALD2(i)=fALD2(i)-1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%107
i=i+1;
Rnames{i} = 'ALDX + OH = CXO3';
k(:,i)=4.90e-12.*exp(405./T);
Gstr{i,1}='ALDX'; Gstr{i,2}='OH';
fALDX(i)=fALDX(i)-1; fOH(i)=fOH(i)-1; fCXO3(i)=fCXO3(i)+1;

%108
i=i+1;
Rnames{i} = 'ALDX + NO3 = CXO3+ HNO3';
k(:,i)=6.3e-15;
Gstr{i,1}='ALDX'; Gstr{i,2}='NO3';
fALDX(i)=fALDX(i)-1; fNO3(i)=fNO3(i)-1; fCXO3(i)=fCXO3(i)+1; fHNO3(i)=fHNO3(i)+1;

%109
i=i+1;
Rnames{i} = 'ALDX = ALD2+ XO2H+ RO2+ CO+ HO2';
k(:,i)=JC2CHO;
Gstr{i,1}='ALDX';
fALDX(i)=fALDX(i)-1; fALD2(i)=fALD2(i)+1; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%110
i=i+1;
Rnames{i} = 'GLYD + OH = 0.8C2O3+ 0.2GLY+ 0.2HO2';
k(:,i)=8e-12;
Gstr{i,1}='GLYD'; Gstr{i,2}='OH';
fGLYD(i)=fGLYD(i)-1; fOH(i)=fOH(i)-1; fC2O3(i)=fC2O3(i)+0.8; fGLY(i)=fGLY(i)+0.2; fHO2(i)=fHO2(i)+0.2;

%111
i=i+1;
Rnames{i} = 'GLYD = 0.74FORM+ 0.11GLY+ 0.15MEOH+ 0.11XO2H+ 0.11RO2+ 0.89CO+ 0.19OH+ 1.4HO2';
k(:,i)=JGLYD;
Gstr{i,1}='GLYD';
fGLYD(i)=fGLYD(i)-1; fFORM(i)=fFORM(i)+0.74; fGLY(i)=fGLY(i)+0.11; fMEOH(i)=fMEOH(i)+0.15; fXO2H(i)=fXO2H(i)+0.11; fRO2(i)=fRO2(i)+0.11; fCO(i)=fCO(i)+0.89; fOH(i)=fOH(i)+0.19; fHO2(i)=fHO2(i)+1.4;

%112
i=i+1;
Rnames{i} = 'GLYD + NO3 = C2O3+ HNO3';
k(:,i)=1.40e-12.*exp(-1860./T);
Gstr{i,1}='GLYD'; Gstr{i,2}='NO3';
fGLYD(i)=fGLYD(i)-1; fNO3(i)=fNO3(i)-1; fC2O3(i)=fC2O3(i)+1; fHNO3(i)=fHNO3(i)+1;

%113
i=i+1;
Rnames{i} = 'GLY + OH = 0.2XO2+ 0.2RO2+ 1.8CO+ HO2';
k(:,i)=3.10e-12.*exp(340./T);
Gstr{i,1}='GLY'; Gstr{i,2}='OH';
fGLY(i)=fGLY(i)-1; fOH(i)=fOH(i)-1; fXO2(i)=fXO2(i)+0.2; fRO2(i)=fRO2(i)+0.2; fCO(i)=fCO(i)+1.8; fHO2(i)=fHO2(i)+1;

%114
i=i+1;
Rnames{i} = 'GLY = CO+ HO2+ CO+ HO2';
k(:,i)=JGLY;
Gstr{i,1}='GLY';
fGLY(i)=fGLY(i)-1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%115
i=i+1;
Rnames{i} = 'GLY + NO3 = HNO3+ 0.5XO2+ 0.5RO2+ 1.5CO+ HO2';
k(:,i)=4e-16;
Gstr{i,1}='GLY'; Gstr{i,2}='NO3';
fGLY(i)=fGLY(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+1; fXO2(i)=fXO2(i)+0.5; fRO2(i)=fRO2(i)+0.5; fCO(i)=fCO(i)+1.5; fHO2(i)=fHO2(i)+1;

%116
i=i+1;
Rnames{i} = 'MGLY = C2O3+ CO+ HO2';
k(:,i)=JMGLY;
Gstr{i,1}='MGLY';
fMGLY(i)=fMGLY(i)-1; fC2O3(i)=fC2O3(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%117
i=i+1;
Rnames{i} = 'MGLY + NO3 = C2O3+ HNO3+ XO2+ RO2';
k(:,i)=5e-16;
Gstr{i,1}='MGLY'; Gstr{i,2}='NO3';
fMGLY(i)=fMGLY(i)-1; fNO3(i)=fNO3(i)-1; fC2O3(i)=fC2O3(i)+1; fHNO3(i)=fHNO3(i)+1; fXO2(i)=fXO2(i)+1; fRO2(i)=fRO2(i)+1;

%118
i=i+1;
Rnames{i} = 'MGLY + OH = C2O3+ CO';
k(:,i)=1.90e-12.*exp(575./T);
Gstr{i,1}='MGLY'; Gstr{i,2}='OH';
fMGLY(i)=fMGLY(i)-1; fOH(i)=fOH(i)-1; fC2O3(i)=fC2O3(i)+1; fCO(i)=fCO(i)+1;

%119
i=i+1;
Rnames{i} = 'OH + H2 = HO2';
k(:,i)=7.70e-12.*exp(-2100./T) ;
Gstr{i,1}='OH'; Gstr{i,2}='H2';
fOH(i)=fOH(i)-1; fH2(i)=fH2(i)-1; fHO2(i)=fHO2(i)+1;

%120
i=i+1;
Rnames{i} = 'CO + OH = HO2';
k(:,i)=K_CO_OH;
Gstr{i,1}='CO'; Gstr{i,2}='OH';
fCO(i)=fCO(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+1;

%121
i=i+1;
Rnames{i} = 'CH4 + OH = MEO2+ RO2';
k(:,i)=1.85e-12.*exp(-1690./T);
Gstr{i,1}='CH4'; Gstr{i,2}='OH';
fCH4(i)=fCH4(i)-1; fOH(i)=fOH(i)-1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1;

%122
i=i+1;
Rnames{i} = 'ETHA + OH = 0.991ALD2+ 0.991XO2H+ 0.009XO2N+ RO2';
k(:,i)=6.90e-12.*exp(-1000./T);
Gstr{i,1}='ETHA'; Gstr{i,2}='OH';
fETHA(i)=fETHA(i)-1; fOH(i)=fOH(i)-1; fALD2(i)=fALD2(i)+0.991; fXO2H(i)=fXO2H(i)+0.991; fXO2N(i)=fXO2N(i)+0.009; fRO2(i)=fRO2(i)+1;

%123
i=i+1;
Rnames{i} = 'MEOH + OH = FORM+ HO2';
k(:,i)=2.85e-12.*exp(-345./T);
Gstr{i,1}='MEOH'; Gstr{i,2}='OH';
fMEOH(i)=fMEOH(i)-1; fOH(i)=fOH(i)-1; fFORM(i)=fFORM(i)+1; fHO2(i)=fHO2(i)+1;

%124
i=i+1;
Rnames{i} = 'ETOH + OH = 0.95ALD2+ 0.078FORM+ 0.011GLYD+ 0.1XO2H+ 0.1RO2+ 0.9HO2';
k(:,i)=3.00e-12.*exp(20./T);
Gstr{i,1}='ETOH'; Gstr{i,2}='OH';
fETOH(i)=fETOH(i)-1; fOH(i)=fOH(i)-1; fALD2(i)=fALD2(i)+0.95; fFORM(i)=fFORM(i)+0.078; fGLYD(i)=fGLYD(i)+0.011; fXO2H(i)=fXO2H(i)+0.1; fRO2(i)=fRO2(i)+0.1; fHO2(i)=fHO2(i)+0.9;

%125
i=i+1;
Rnames{i} = 'KET = 0.5ALD2+ 0.5C2O3+ 0.5CXO3+ 0.5MEO2- 2.5PAR+ 0.5XO2H+ RO2';
k(:,i)=JMEK;
Gstr{i,1}='KET';
fKET(i)=fKET(i)-1; fALD2(i)=fALD2(i)+0.5; fC2O3(i)=fC2O3(i)+0.5; fCXO3(i)=fCXO3(i)+0.5; fMEO2(i)=fMEO2(i)+0.5; fPAR(i)=fPAR(i)-2.5; fXO2H(i)=fXO2H(i)+0.5; fRO2(i)=fRO2(i)+1;

%126
i=i+1;
Rnames{i} = 'ACET = 0.62C2O3+ 1.38MEO2+ 1.38RO2+ 0.38CO';
k(:,i)=JACET;
Gstr{i,1}='ACET';
fACET(i)=fACET(i)-1; fC2O3(i)=fC2O3(i)+0.62; fMEO2(i)=fMEO2(i)+1.38; fRO2(i)=fRO2(i)+1.38; fCO(i)=fCO(i)+0.38;

%127
i=i+1;
Rnames{i} = 'ACET + OH = C2O3+ FORM+ XO2+ RO2';
k(:,i)=1.41e-12.*exp(-620.6./T);
Gstr{i,1}='ACET'; Gstr{i,2}='OH';
fACET(i)=fACET(i)-1; fOH(i)=fOH(i)-1; fC2O3(i)=fC2O3(i)+1; fFORM(i)=fFORM(i)+1; fXO2(i)=fXO2(i)+1; fRO2(i)=fRO2(i)+1;

%128
i=i+1;
Rnames{i} = 'PRPA + OH = XPRP';
k(:,i)=7.60e-12.*exp(-585./T);
Gstr{i,1}='PRPA'; Gstr{i,2}='OH';
fPRPA(i)=fPRPA(i)-1; fOH(i)=fOH(i)-1; fXPRP(i)=fXPRP(i)+1;

%129
i=i+1;
Rnames{i} = 'PAR + OH = XPAR';
k(:,i)=8.10e-13;
Gstr{i,1}='PAR'; Gstr{i,2}='OH';
fPAR(i)=fPAR(i)-1; fOH(i)=fOH(i)-1; fXPAR(i)=fXPAR(i)+1;

%130
i=i+1;
Rnames{i} = 'ROR = 0.42ACET+ 0.74ALD2+ 0.37ALDX+ 0.2KET- 2.7PAR+ 0.02ROR+ 0.94XO2H+ 0.04XO2N+ 0.98RO2';
k(:,i)=5.70e+12.*exp(-5780./T);
Gstr{i,1}='ROR';
fROR(i)=fROR(i)-1; fACET(i)=fACET(i)+0.42; fALD2(i)=fALD2(i)+0.74; fALDX(i)=fALDX(i)+0.37; fKET(i)=fKET(i)+0.2; fPAR(i)=fPAR(i)-2.7; fROR(i)=fROR(i)+0.02; fXO2H(i)=fXO2H(i)+0.94; fXO2N(i)=fXO2N(i)+0.04; fRO2(i)=fRO2(i)+0.98;

%131
i=i+1;
Rnames{i} = 'ROR = KET+ HO2';
k(:,i)=1.50e-14.*exp(-200./T).*M.*0.21;
Gstr{i,1}='ROR';
fROR(i)=fROR(i)-1; fKET(i)=fKET(i)+1; fHO2(i)=fHO2(i)+1;

%132
i=i+1;
Rnames{i} = 'ROR + NO2 = NTR1';
k(:,i)=8.60e-12.*exp(400./T);
Gstr{i,1}='ROR'; Gstr{i,2}='NO2';
fROR(i)=fROR(i)-1; fNO2(i)=fNO2(i)-1; fNTR1(i)=fNTR1(i)+1;

%133
i=i+1;
Rnames{i} = 'ETHY + OH = 0.3FACD+ 0.7GLY+ 0.3CO+ 0.7OH+ 0.3HO2';
k(:,i)=K_OH_ETHY;
Gstr{i,1}='ETHY'; Gstr{i,2}='OH';
fETHY(i)=fETHY(i)-1; fOH(i)=fOH(i)-1; fFACD(i)=fFACD(i)+0.3; fGLY(i)=fGLY(i)+0.7; fCO(i)=fCO(i)+0.3; fOH(i)=fOH(i)+0.7; fHO2(i)=fHO2(i)+0.3;

%134
i=i+1;
Rnames{i} = 'ETH + OH = 1.56FORM+ 0.22GLYD+ XO2H+ RO2';
k(:,i)=K_OH_ETH;
Gstr{i,1}='ETH'; Gstr{i,2}='OH';
fETH(i)=fETH(i)-1; fOH(i)=fOH(i)-1; fFORM(i)=fFORM(i)+1.56; fGLYD(i)=fGLYD(i)+0.22; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%135
i=i+1;
Rnames{i} = 'ETH + O3 = 0.42FACD+ FORM+ 0.35CO+ 0.17OH+ 0.27HO2';
k(:,i)=6.82e-15.*exp(-2500./T);
Gstr{i,1}='ETH'; Gstr{i,2}='O3';
fETH(i)=fETH(i)-1; fO3(i)=fO3(i)-1; fFACD(i)=fFACD(i)+0.42; fFORM(i)=fFORM(i)+1; fCO(i)=fCO(i)+0.35; fOH(i)=fOH(i)+0.17; fHO2(i)=fHO2(i)+0.27;

%136
i=i+1;
Rnames{i} = 'ETH + NO3 = 1.125FORM+ 0.5NTR1+ 0.5XO2+ 0.5XO2H+ 0.5NO2+ RO2';
k(:,i)=3.30e-12.*exp(-2880./T);
Gstr{i,1}='ETH'; Gstr{i,2}='NO3';
fETH(i)=fETH(i)-1; fNO3(i)=fNO3(i)-1; fFORM(i)=fFORM(i)+1.125; fNTR1(i)=fNTR1(i)+0.5; fXO2(i)=fXO2(i)+0.5; fXO2H(i)=fXO2H(i)+0.5; fNO2(i)=fNO2(i)+0.5; fRO2(i)=fRO2(i)+1;

%137
i=i+1;
Rnames{i} = 'OLE + OH = 0.488ALD2+ 0.488ALDX+ 0.781FORM- 0.73PAR+ 0.195XO2+ 0.976XO2H+ 0.024XO2N+ 1.195RO2';
k(:,i)=K_OH_OLE;
Gstr{i,1}='OLE'; Gstr{i,2}='OH';
fOLE(i)=fOLE(i)-1; fOH(i)=fOH(i)-1; fALD2(i)=fALD2(i)+0.488; fALDX(i)=fALDX(i)+0.488; fFORM(i)=fFORM(i)+0.781; fPAR(i)=fPAR(i)-0.73; fXO2(i)=fXO2(i)+0.195; fXO2H(i)=fXO2H(i)+0.976; fXO2N(i)=fXO2N(i)+0.024; fRO2(i)=fRO2(i)+1.195;

%138
i=i+1;
Rnames{i} = 'OLE + O3 = 0.13AACD+ 0.295ALD2+ 0.27ALDX+ 0.09FACD+ 0.555FORM+ 0.075GLY+ 0.04H2O2+ 0.075MGLY- 0.79PAR+ 0.15XO2H+ 0.15RO2+ 0.378CO+ 0.334OH+ 0.08HO2';
k(:,i)=5.50e-15.*exp(-1880./T);
Gstr{i,1}='OLE'; Gstr{i,2}='O3';
fOLE(i)=fOLE(i)-1; fO3(i)=fO3(i)-1; fAACD(i)=fAACD(i)+0.13; fALD2(i)=fALD2(i)+0.295; fALDX(i)=fALDX(i)+0.27; fFACD(i)=fFACD(i)+0.09; fFORM(i)=fFORM(i)+0.555; fGLY(i)=fGLY(i)+0.075; fH2O2(i)=fH2O2(i)+0.04; fMGLY(i)=fMGLY(i)+0.075; fPAR(i)=fPAR(i)-0.79; fXO2H(i)=fXO2H(i)+0.15; fRO2(i)=fRO2(i)+0.15; fCO(i)=fCO(i)+0.378; fOH(i)=fOH(i)+0.334; fHO2(i)=fHO2(i)+0.08;

%139
i=i+1;
Rnames{i} = 'OLE + NO3 = 0.25ALD2+ 0.375ALDX+ 0.5FORM+ 0.5NTR1- PAR+ 0.48XO2+ 0.48XO2H+ 0.04XO2N+ 0.5NO2+ RO2';
k(:,i)=4.60e-13.*exp(-1155./T);
Gstr{i,1}='OLE'; Gstr{i,2}='NO3';
fOLE(i)=fOLE(i)-1; fNO3(i)=fNO3(i)-1; fALD2(i)=fALD2(i)+0.25; fALDX(i)=fALDX(i)+0.375; fFORM(i)=fFORM(i)+0.5; fNTR1(i)=fNTR1(i)+0.5; fPAR(i)=fPAR(i)-1.0; fXO2(i)=fXO2(i)+0.48; fXO2H(i)=fXO2H(i)+0.48; fXO2N(i)=fXO2N(i)+0.04; fNO2(i)=fNO2(i)+0.5; fRO2(i)=fRO2(i)+1;

%140
i=i+1;
Rnames{i} = 'IOLE + OH = 1.3ALD2+ 0.7ALDX+ XO2H+ RO2';
k(:,i)=1.05e-11.*exp(519./T);
Gstr{i,1}='IOLE'; Gstr{i,2}='OH';
fIOLE(i)=fIOLE(i)-1; fOH(i)=fOH(i)-1; fALD2(i)=fALD2(i)+1.3; fALDX(i)=fALDX(i)+0.7; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%141
i=i+1;
Rnames{i} = 'IOLE + O3 = 0.08AACD+ 0.732ALD2+ 0.442ALDX+ 0.128FORM+ 0.24GLY+ 0.08H2O2+ 0.06MGLY+ 0.29PAR+ 0.3XO2H+ 0.3RO2+ 0.245CO+ 0.5OH';
k(:,i)=4.70e-15.*exp(-1013./T);
Gstr{i,1}='IOLE'; Gstr{i,2}='O3';
fIOLE(i)=fIOLE(i)-1; fO3(i)=fO3(i)-1; fAACD(i)=fAACD(i)+0.08; fALD2(i)=fALD2(i)+0.732; fALDX(i)=fALDX(i)+0.442; fFORM(i)=fFORM(i)+0.128; fGLY(i)=fGLY(i)+0.24; fH2O2(i)=fH2O2(i)+0.08; fMGLY(i)=fMGLY(i)+0.06; fPAR(i)=fPAR(i)+0.29; fXO2H(i)=fXO2H(i)+0.3; fRO2(i)=fRO2(i)+0.3; fCO(i)=fCO(i)+0.245; fOH(i)=fOH(i)+0.5;

%142
i=i+1;
Rnames{i} = 'IOLE + NO3 = 0.5ALD2+ 0.625ALDX+ 0.5NTR1+ PAR+ 0.48XO2+ 0.48XO2H+ 0.04XO2N+ 0.5NO2+ RO2';
k(:,i)=3.7e-13;
Gstr{i,1}='IOLE'; Gstr{i,2}='NO3';
fIOLE(i)=fIOLE(i)-1; fNO3(i)=fNO3(i)-1; fALD2(i)=fALD2(i)+0.5; fALDX(i)=fALDX(i)+0.625; fNTR1(i)=fNTR1(i)+0.5; fPAR(i)=fPAR(i)+1; fXO2(i)=fXO2(i)+0.48; fXO2H(i)=fXO2H(i)+0.48; fXO2N(i)=fXO2N(i)+0.04; fNO2(i)=fNO2(i)+0.5; fRO2(i)=fRO2(i)+1;

%143
i=i+1;
Rnames{i} = 'ISOP + OH = ISO2+ RO2';
k(:,i)=2.70e-11.*exp(390./T);
Gstr{i,1}='ISOP'; Gstr{i,2}='OH';
fISOP(i)=fISOP(i)-1; fOH(i)=fOH(i)-1; fISO2(i)=fISO2(i)+1; fRO2(i)=fRO2(i)+1;

%144
i=i+1;
Rnames{i} = 'ISO2 + NO = 0.673FORM+ 0.1INTR+ 0.9ISPD+ 0.082XO2H+ 0.9NO2+ 0.082RO2+ 0.818HO2';
k(:,i)=2.39e-12.*exp(365./T);
Gstr{i,1}='ISO2'; Gstr{i,2}='NO';
fISO2(i)=fISO2(i)-1; fNO(i)=fNO(i)-1; fFORM(i)=fFORM(i)+0.673; fINTR(i)=fINTR(i)+0.1; fISPD(i)=fISPD(i)+0.9; fXO2H(i)=fXO2H(i)+0.082; fNO2(i)=fNO2(i)+0.9; fRO2(i)=fRO2(i)+0.082; fHO2(i)=fHO2(i)+0.818;

%145
i=i+1;
Rnames{i} = 'ISO2 + HO2 = 0.12FORM+ 0.12ISPD+ 0.88ISPX+ 0.12OH+ 0.12HO2';
k(:,i)=7.43e-13.*exp(700./T);
Gstr{i,1}='ISO2'; Gstr{i,2}='HO2';
fISO2(i)=fISO2(i)-1; fHO2(i)=fHO2(i)-1; fFORM(i)=fFORM(i)+0.12; fISPD(i)=fISPD(i)+0.12; fISPX(i)=fISPX(i)+0.88; fOH(i)=fOH(i)+0.12; fHO2(i)=fHO2(i)+0.12;

%146
i=i+1;
Rnames{i} = 'C2O3 + ISO2 = 0.2AACD+ 0.598FORM+ ISPD+ 0.8MEO2+ 0.072XO2H+ 0.872RO2+ 0.728HO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='ISO2';
fC2O3(i)=fC2O3(i)-1; fISO2(i)=fISO2(i)-1; fAACD(i)=fAACD(i)+0.2; fFORM(i)=fFORM(i)+0.598; fISPD(i)=fISPD(i)+1.0; fMEO2(i)=fMEO2(i)+0.8; fXO2H(i)=fXO2H(i)+0.072; fRO2(i)=fRO2(i)+0.872; fHO2(i)=fHO2(i)+0.728;

%147
i=i+1;
Rnames{i} = 'ISO2 + RO2 = 0.598FORM+ ISPD+ 0.072XO2H+ 1.072RO2+ 0.728HO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='ISO2'; Gstr{i,2}='RO2';
fISO2(i)=fISO2(i)-1; fRO2(i)=fRO2(i)-1; fFORM(i)=fFORM(i)+0.598; fISPD(i)=fISPD(i)+1.0; fXO2H(i)=fXO2H(i)+0.072; fRO2(i)=fRO2(i)+1.072; fHO2(i)=fHO2(i)+0.728;

%148
i=i+1;
Rnames{i} = 'ISO2 = HPLD+ HO2';
k(:,i)=3.30e+9.*exp(-8300./T);
Gstr{i,1}='ISO2';
fISO2(i)=fISO2(i)-1; fHPLD(i)=fHPLD(i)+1; fHO2(i)=fHO2(i)+1;

%149
i=i+1;
Rnames{i} = 'ISOP + O3 = 0.15ALDX+ 0.2CXO3+ 0.6FORM+ 0.65ISPD+ 0.35PAR+ 0.2XO2+ 0.2RO2+ 0.066CO+ 0.266OH+ 0.066HO2';
k(:,i)=1.03e-14.*exp(-1995./T);
Gstr{i,1}='ISOP'; Gstr{i,2}='O3';
fISOP(i)=fISOP(i)-1; fO3(i)=fO3(i)-1; fALDX(i)=fALDX(i)+0.15; fCXO3(i)=fCXO3(i)+0.2; fFORM(i)=fFORM(i)+0.6; fISPD(i)=fISPD(i)+0.65; fPAR(i)=fPAR(i)+0.35; fXO2(i)=fXO2(i)+0.2; fRO2(i)=fRO2(i)+0.2; fCO(i)=fCO(i)+0.066; fOH(i)=fOH(i)+0.266; fHO2(i)=fHO2(i)+0.066;

%150
i=i+1;
Rnames{i} = 'ISOP + NO3 = 0.35FORM+ 0.35ISPD+ 0.65NTR2+ 0.33XO2+ 0.64XO2H+ 0.03XO2N+ 0.35NO2+ RO2';
k(:,i)=3.03e-12.*exp(-448./T);
Gstr{i,1}='ISOP'; Gstr{i,2}='NO3';
fISOP(i)=fISOP(i)-1; fNO3(i)=fNO3(i)-1; fFORM(i)=fFORM(i)+0.35; fISPD(i)=fISPD(i)+0.35; fNTR2(i)=fNTR2(i)+0.65; fXO2(i)=fXO2(i)+0.33; fXO2H(i)=fXO2H(i)+0.64; fXO2N(i)=fXO2N(i)+0.03; fNO2(i)=fNO2(i)+0.35; fRO2(i)=fRO2(i)+1;

%151
i=i+1;
Rnames{i} = 'ISPD + OH = 0.137ACET+ 0.269C2O3+ 0.269GLYD+ 0.115MEO2+ 0.115MGLY+ 0.457OPO3+ 0.117PAR+ 0.521XO2+ 0.022XO2N+ 0.658RO2+ 0.137CO+ 0.137HO2';
k(:,i)=5.58e-12.*exp(511./T);
Gstr{i,1}='ISPD'; Gstr{i,2}='OH';
fISPD(i)=fISPD(i)-1; fOH(i)=fOH(i)-1; fACET(i)=fACET(i)+0.137; fC2O3(i)=fC2O3(i)+0.269; fGLYD(i)=fGLYD(i)+0.269; fMEO2(i)=fMEO2(i)+0.115; fMGLY(i)=fMGLY(i)+0.115; fOPO3(i)=fOPO3(i)+0.457; fPAR(i)=fPAR(i)+0.117; fXO2(i)=fXO2(i)+0.521; fXO2N(i)=fXO2N(i)+0.022; fRO2(i)=fRO2(i)+0.658; fCO(i)=fCO(i)+0.137; fHO2(i)=fHO2(i)+0.137;

%152
i=i+1;
Rnames{i} = 'ISPD + O3 = 0.17ACET+ 0.04ALD2+ 0.143C2O3+ 0.15FACD+ 0.231FORM+ 0.17GLY+ 0.531MGLY+ 0.543CO+ 0.461OH+ 0.398HO2';
k(:,i)=3.88e-15.*exp(-1770./T);
Gstr{i,1}='ISPD'; Gstr{i,2}='O3';
fISPD(i)=fISPD(i)-1; fO3(i)=fO3(i)-1; fACET(i)=fACET(i)+0.17; fALD2(i)=fALD2(i)+0.04; fC2O3(i)=fC2O3(i)+0.143; fFACD(i)=fFACD(i)+0.15; fFORM(i)=fFORM(i)+0.231; fGLY(i)=fGLY(i)+0.17; fMGLY(i)=fMGLY(i)+0.531; fCO(i)=fCO(i)+0.543; fOH(i)=fOH(i)+0.461; fHO2(i)=fHO2(i)+0.398;

%153
i=i+1;
Rnames{i} = 'ISPD + NO3 = 0.717CXO3+ 0.113GLYD+ 0.717HNO3+ 0.113MGLY+ 0.142NTR2+ 0.717PAR+ 0.142XO2+ 0.142XO2H+ 0.142NO2+ 0.284RO2';
k(:,i)=4.10e-12.*exp(-1860./T);
Gstr{i,1}='ISPD'; Gstr{i,2}='NO3';
fISPD(i)=fISPD(i)-1; fNO3(i)=fNO3(i)-1; fCXO3(i)=fCXO3(i)+0.717; fGLYD(i)=fGLYD(i)+0.113; fHNO3(i)=fHNO3(i)+0.717; fMGLY(i)=fMGLY(i)+0.113; fNTR2(i)=fNTR2(i)+0.142; fPAR(i)=fPAR(i)+0.717; fXO2(i)=fXO2(i)+0.142; fXO2H(i)=fXO2H(i)+0.142; fNO2(i)=fNO2(i)+0.142; fRO2(i)=fRO2(i)+0.284;

%154
i=i+1;
Rnames{i} = 'ISPD = 0.17ACET+ 0.208C2O3+ 0.26FORM+ 0.128GLYD+ 0.34MEO2+ 0.24OLE+ 0.24PAR+ 0.16XO2+ 0.34XO2H+ 0.84RO2+ 0.76HO2';
k(:,i)=JACRO .* 0.0036./0.0065;
Gstr{i,1}='ISPD';
fISPD(i)=fISPD(i)-1; fACET(i)=fACET(i)+0.17; fC2O3(i)=fC2O3(i)+0.208; fFORM(i)=fFORM(i)+0.26; fGLYD(i)=fGLYD(i)+0.128; fMEO2(i)=fMEO2(i)+0.34; fOLE(i)=fOLE(i)+0.24; fPAR(i)=fPAR(i)+0.24; fXO2(i)=fXO2(i)+0.16; fXO2H(i)=fXO2H(i)+0.34; fRO2(i)=fRO2(i)+0.84; fHO2(i)=fHO2(i)+0.76;

%155
i=i+1;
Rnames{i} = 'ISPX + OH = 0.029ALDX+ 0.904EPOX+ 0.029IOLE+ 0.067ISO2+ 0.067RO2+ 0.933OH';
k(:,i)=2.23e-11.*exp(372./T);
Gstr{i,1}='ISPX'; Gstr{i,2}='OH';
fISPX(i)=fISPX(i)-1; fOH(i)=fOH(i)-1; fALDX(i)=fALDX(i)+0.029; fEPOX(i)=fEPOX(i)+0.904; fIOLE(i)=fIOLE(i)+0.029; fISO2(i)=fISO2(i)+0.067; fRO2(i)=fRO2(i)+0.067; fOH(i)=fOH(i)+0.933;

%156
i=i+1;
Rnames{i} = 'HPLD = ISPD+ OH';
k(:,i)=JHPLD;
Gstr{i,1}='HPLD';
fHPLD(i)=fHPLD(i)-1; fISPD(i)=fISPD(i)+1; fOH(i)=fOH(i)+1;

%157
i=i+1;
Rnames{i} = 'HPLD + NO3 = HNO3+ ISPD';
k(:,i)=6.00e-12.*exp(-1860./T);
Gstr{i,1}='HPLD'; Gstr{i,2}='NO3';
fHPLD(i)=fHPLD(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+1; fISPD(i)=fISPD(i)+1;

%158
i=i+1;
Rnames{i} = 'EPOX + OH = EPX2+ RO2';
k(:,i)=5.78e-11.*exp(-400./T);
Gstr{i,1}='EPOX'; Gstr{i,2}='OH';
fEPOX(i)=fEPOX(i)-1; fOH(i)=fOH(i)-1; fEPX2(i)=fEPX2(i)+1; fRO2(i)=fRO2(i)+1;

%159
i=i+1;
Rnames{i} = 'EPX2 + HO2 = 0.074FACD+ 0.375FORM+ 0.275GLY+ 0.275GLYD+ 0.275MGLY+ 2.175PAR+ 0.251CO+ 1.125OH+ 0.825HO2';
k(:,i)=7.43e-13.*exp(700./T);
Gstr{i,1}='EPX2'; Gstr{i,2}='HO2';
fEPX2(i)=fEPX2(i)-1; fHO2(i)=fHO2(i)-1; fFACD(i)=fFACD(i)+0.074; fFORM(i)=fFORM(i)+0.375; fGLY(i)=fGLY(i)+0.275; fGLYD(i)=fGLYD(i)+0.275; fMGLY(i)=fMGLY(i)+0.275; fPAR(i)=fPAR(i)+2.175; fCO(i)=fCO(i)+0.251; fOH(i)=fOH(i)+1.125; fHO2(i)=fHO2(i)+0.825;

%160
i=i+1;
Rnames{i} = 'EPX2 + NO = 0.375FORM+ 0.275GLY+ 0.275GLYD+ 0.275MGLY+ 2.175PAR+ NO2+ 0.251CO+ 0.125OH+ 0.825HO2';
k(:,i)=2.39e-12.*exp(365./T);
Gstr{i,1}='EPX2'; Gstr{i,2}='NO';
fEPX2(i)=fEPX2(i)-1; fNO(i)=fNO(i)-1; fFORM(i)=fFORM(i)+0.375; fGLY(i)=fGLY(i)+0.275; fGLYD(i)=fGLYD(i)+0.275; fMGLY(i)=fMGLY(i)+0.275; fPAR(i)=fPAR(i)+2.175; fNO2(i)=fNO2(i)+1; fCO(i)=fCO(i)+0.251; fOH(i)=fOH(i)+0.125; fHO2(i)=fHO2(i)+0.825;

%161
i=i+1;
Rnames{i} = 'C2O3 + EPX2 = 0.2AACD+ 0.3FORM+ 0.22GLY+ 0.22GLYD+ 0.8MEO2+ 0.22MGLY+ 1.74PAR+ 0.8RO2+ 0.2CO+ 0.1OH+ 0.66HO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='EPX2';
fC2O3(i)=fC2O3(i)-1; fEPX2(i)=fEPX2(i)-1; fAACD(i)=fAACD(i)+0.2; fFORM(i)=fFORM(i)+0.3; fGLY(i)=fGLY(i)+0.22; fGLYD(i)=fGLYD(i)+0.22; fMEO2(i)=fMEO2(i)+0.8; fMGLY(i)=fMGLY(i)+0.22; fPAR(i)=fPAR(i)+1.74; fRO2(i)=fRO2(i)+0.8; fCO(i)=fCO(i)+0.2; fOH(i)=fOH(i)+0.1; fHO2(i)=fHO2(i)+0.66;

%162
i=i+1;
Rnames{i} = 'EPX2 + RO2 = 0.375FORM+ 0.275GLY+ 0.275GLYD+ 0.275MGLY+ 2.175PAR+ RO2+ 0.251CO+ 0.125OH+ 0.825HO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='EPX2'; Gstr{i,2}='RO2';
fEPX2(i)=fEPX2(i)-1; fRO2(i)=fRO2(i)-1; fFORM(i)=fFORM(i)+0.375; fGLY(i)=fGLY(i)+0.275; fGLYD(i)=fGLYD(i)+0.275; fMGLY(i)=fMGLY(i)+0.275; fPAR(i)=fPAR(i)+2.175; fRO2(i)=fRO2(i)+1; fCO(i)=fCO(i)+0.251; fOH(i)=fOH(i)+0.125; fHO2(i)=fHO2(i)+0.825;

%163
i=i+1;
Rnames{i} = 'INTR + OH = 0.078ALDX+ 0.185FACD+ 0.592FORM+ 0.331GLYD+ 0.104INTR+ 0.266NTR2+ 0.098OLE+ 2.7PAR+ 0.63XO2+ 0.37XO2H+ 0.444NO2+ 0.185NO3+ RO2';
k(:,i)=3.1e-11;
Gstr{i,1}='INTR'; Gstr{i,2}='OH';
fINTR(i)=fINTR(i)-1; fOH(i)=fOH(i)-1; fALDX(i)=fALDX(i)+0.078; fFACD(i)=fFACD(i)+0.185; fFORM(i)=fFORM(i)+0.592; fGLYD(i)=fGLYD(i)+0.331; fINTR(i)=fINTR(i)+0.104; fNTR2(i)=fNTR2(i)+0.266; fOLE(i)=fOLE(i)+0.098; fPAR(i)=fPAR(i)+2.7; fXO2(i)=fXO2(i)+0.63; fXO2H(i)=fXO2H(i)+0.37; fNO2(i)=fNO2(i)+0.444; fNO3(i)=fNO3(i)+0.185; fRO2(i)=fRO2(i)+1;

%164
i=i+1;
Rnames{i} = 'TERP + OH = 0.47ALDX+ 0.28FORM+ 1.66PAR+ 0.5XO2+ 0.75XO2H+ 0.25XO2N+ 1.5RO2';
k(:,i)=1.50e-11.*exp(449./T);
Gstr{i,1}='TERP'; Gstr{i,2}='OH';
fTERP(i)=fTERP(i)-1; fOH(i)=fOH(i)-1; fALDX(i)=fALDX(i)+0.47; fFORM(i)=fFORM(i)+0.28; fPAR(i)=fPAR(i)+1.66; fXO2(i)=fXO2(i)+0.5; fXO2H(i)=fXO2H(i)+0.75; fXO2N(i)=fXO2N(i)+0.25; fRO2(i)=fRO2(i)+1.5;

%165
i=i+1;
Rnames{i} = 'TERP + O3 = 0.21ALDX+ 0.39CXO3+ 0.24FORM+ 7.0PAR+ 0.69XO2+ 0.07XO2H+ 0.18XO2N+ 0.94RO2+ 0.001CO+ 0.57OH';
k(:,i)=1.20e-15.*exp(-821./T);
Gstr{i,1}='TERP'; Gstr{i,2}='O3';
fTERP(i)=fTERP(i)-1; fO3(i)=fO3(i)-1; fALDX(i)=fALDX(i)+0.21; fCXO3(i)=fCXO3(i)+0.39; fFORM(i)=fFORM(i)+0.24; fPAR(i)=fPAR(i)+7.0; fXO2(i)=fXO2(i)+0.69; fXO2H(i)=fXO2H(i)+0.07; fXO2N(i)=fXO2N(i)+0.18; fRO2(i)=fRO2(i)+0.94; fCO(i)=fCO(i)+0.001; fOH(i)=fOH(i)+0.57;

%166
i=i+1;
Rnames{i} = 'TERP + NO3 = 0.47ALDX+ 0.53NTR2+ 0.75XO2+ 0.28XO2H+ 0.25XO2N+ 0.47NO2+ 1.28RO2';
k(:,i)=3.70e-12.*exp(175./T);
Gstr{i,1}='TERP'; Gstr{i,2}='NO3';
fTERP(i)=fTERP(i)-1; fNO3(i)=fNO3(i)-1; fALDX(i)=fALDX(i)+0.47; fNTR2(i)=fNTR2(i)+0.53; fXO2(i)=fXO2(i)+0.75; fXO2H(i)=fXO2H(i)+0.28; fXO2N(i)=fXO2N(i)+0.25; fNO2(i)=fNO2(i)+0.47; fRO2(i)=fRO2(i)+1.28;

%167
i=i+1;
Rnames{i} = 'BENZ + OH = 0.352BZO2+ 0.53CRES+ 0.118OPEN+ 0.352RO2+ 0.118OH+ 0.53HO2';
k(:,i)=2.30e-12.*exp(-190./T);
Gstr{i,1}='BENZ'; Gstr{i,2}='OH';
fBENZ(i)=fBENZ(i)-1; fOH(i)=fOH(i)-1; fBZO2(i)=fBZO2(i)+0.352; fCRES(i)=fCRES(i)+0.53; fOPEN(i)=fOPEN(i)+0.118; fRO2(i)=fRO2(i)+0.352; fOH(i)=fOH(i)+0.118; fHO2(i)=fHO2(i)+0.53;

%168
i=i+1;
Rnames{i} = 'BZO2 + NO = 0.918GLY+ 0.082NTR2+ 0.918OPEN+ 0.918NO2+ 0.918HO2';
k(:,i)=2.70e-12.*exp(360./T);
Gstr{i,1}='BZO2'; Gstr{i,2}='NO';
fBZO2(i)=fBZO2(i)-1; fNO(i)=fNO(i)-1; fGLY(i)=fGLY(i)+0.918; fNTR2(i)=fNTR2(i)+0.082; fOPEN(i)=fOPEN(i)+0.918; fNO2(i)=fNO2(i)+0.918; fHO2(i)=fHO2(i)+0.918;

%169
i=i+1;
Rnames{i} = 'BZO2 + C2O3 = GLY+ MEO2+ OPEN+ RO2+ HO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='BZO2'; Gstr{i,2}='C2O3';
fBZO2(i)=fBZO2(i)-1; fC2O3(i)=fC2O3(i)-1; fGLY(i)=fGLY(i)+1; fMEO2(i)=fMEO2(i)+1; fOPEN(i)=fOPEN(i)+1; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+1;

%170
i=i+1;
Rnames{i} = 'BZO2 + HO2 = ';
k(:,i)=1.90e-13.*exp(1300./T);
Gstr{i,1}='BZO2'; Gstr{i,2}='HO2';
fBZO2(i)=fBZO2(i)-1; fHO2(i)=fHO2(i)-1;

%171
i=i+1;
Rnames{i} = 'BZO2 + RO2 = GLY+ OPEN+ RO2+ HO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='BZO2'; Gstr{i,2}='RO2';
fBZO2(i)=fBZO2(i)-1; fRO2(i)=fRO2(i)-1; fGLY(i)=fGLY(i)+1; fOPEN(i)=fOPEN(i)+1; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+1;

%172
i=i+1;
Rnames{i} = 'TOL + OH = 0.18CRES+ 0.1OPEN+ 0.65TO2+ 0.07XO2H+ 0.72RO2+ 0.1OH+ 0.18HO2';
k(:,i)=1.80e-12.*exp(340./T);
Gstr{i,1}='TOL'; Gstr{i,2}='OH';
fTOL(i)=fTOL(i)-1; fOH(i)=fOH(i)-1; fCRES(i)=fCRES(i)+0.18; fOPEN(i)=fOPEN(i)+0.1; fTO2(i)=fTO2(i)+0.65; fXO2H(i)=fXO2H(i)+0.07; fRO2(i)=fRO2(i)+0.72; fOH(i)=fOH(i)+0.1; fHO2(i)=fHO2(i)+0.18;

%173
i=i+1;
Rnames{i} = 'NO + TO2 = 0.417GLY+ 0.443MGLY+ 0.14NTR2+ 0.66OPEN+ 0.2XOPN+ 0.86NO2+ 0.86HO2';
k(:,i)=2.70e-12.*exp(360./T);
Gstr{i,1}='NO'; Gstr{i,2}='TO2';
fNO(i)=fNO(i)-1; fTO2(i)=fTO2(i)-1; fGLY(i)=fGLY(i)+0.417; fMGLY(i)=fMGLY(i)+0.443; fNTR2(i)=fNTR2(i)+0.14; fOPEN(i)=fOPEN(i)+0.66; fXOPN(i)=fXOPN(i)+0.2; fNO2(i)=fNO2(i)+0.86; fHO2(i)=fHO2(i)+0.86;

%174
i=i+1;
Rnames{i} = 'C2O3 + TO2 = 0.48GLY+ MEO2+ 0.52MGLY+ 0.77OPEN+ 0.23XOPN+ RO2+ HO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='TO2';
fC2O3(i)=fC2O3(i)-1; fTO2(i)=fTO2(i)-1; fGLY(i)=fGLY(i)+0.48; fMEO2(i)=fMEO2(i)+1; fMGLY(i)=fMGLY(i)+0.52; fOPEN(i)=fOPEN(i)+0.77; fXOPN(i)=fXOPN(i)+0.23; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+1;

%175
i=i+1;
Rnames{i} = 'TO2 + HO2 = ';
k(:,i)=1.90e-13.*exp(1300./T);
Gstr{i,1}='TO2'; Gstr{i,2}='HO2';
fTO2(i)=fTO2(i)-1; fHO2(i)=fHO2(i)-1;

%176
i=i+1;
Rnames{i} = 'TO2 + RO2 = 0.48GLY+ 0.52MGLY+ 0.77OPEN+ 0.23XOPN+ RO2+ HO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='TO2'; Gstr{i,2}='RO2';
fTO2(i)=fTO2(i)-1; fRO2(i)=fRO2(i)-1; fGLY(i)=fGLY(i)+0.48; fMGLY(i)=fMGLY(i)+0.52; fOPEN(i)=fOPEN(i)+0.77; fXOPN(i)=fXOPN(i)+0.23; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+1;

%177
i=i+1;
Rnames{i} = 'XYL + OH = 0.155CRES+ 0.544XLO2+ 0.058XO2H+ 0.244XOPN+ 0.602RO2+ 0.244OH+ 0.155HO2';
k(:,i)=1.85e-11;
Gstr{i,1}='XYL'; Gstr{i,2}='OH';
fXYL(i)=fXYL(i)-1; fOH(i)=fOH(i)-1; fCRES(i)=fCRES(i)+0.155; fXLO2(i)=fXLO2(i)+0.544; fXO2H(i)=fXO2H(i)+0.058; fXOPN(i)=fXOPN(i)+0.244; fRO2(i)=fRO2(i)+0.602; fOH(i)=fOH(i)+0.244; fHO2(i)=fHO2(i)+0.155;

%178
i=i+1;
Rnames{i} = 'NO + XLO2 = 0.221GLY+ 0.675MGLY+ 0.14NTR2+ 0.3OPEN+ 0.56XOPN+ 0.86NO2+ 0.86HO2';
k(:,i)=2.70e-12.*exp(360./T);
Gstr{i,1}='NO'; Gstr{i,2}='XLO2';
fNO(i)=fNO(i)-1; fXLO2(i)=fXLO2(i)-1; fGLY(i)=fGLY(i)+0.221; fMGLY(i)=fMGLY(i)+0.675; fNTR2(i)=fNTR2(i)+0.14; fOPEN(i)=fOPEN(i)+0.3; fXOPN(i)=fXOPN(i)+0.56; fNO2(i)=fNO2(i)+0.86; fHO2(i)=fHO2(i)+0.86;

%179
i=i+1;
Rnames{i} = 'XLO2 + HO2 = ';
k(:,i)=1.90e-13.*exp(1300./T);
Gstr{i,1}='XLO2'; Gstr{i,2}='HO2';
fXLO2(i)=fXLO2(i)-1; fHO2(i)=fHO2(i)-1;

%180
i=i+1;
Rnames{i} = 'C2O3 + XLO2 = 0.26GLY+ MEO2+ 0.77MGLY+ 0.35OPEN+ 0.65XOPN+ RO2+ HO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='XLO2';
fC2O3(i)=fC2O3(i)-1; fXLO2(i)=fXLO2(i)-1; fGLY(i)=fGLY(i)+0.26; fMEO2(i)=fMEO2(i)+1; fMGLY(i)=fMGLY(i)+0.77; fOPEN(i)=fOPEN(i)+0.35; fXOPN(i)=fXOPN(i)+0.65; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+1;

%181
i=i+1;
Rnames{i} = 'XLO2 + RO2 = 0.26GLY+ 0.77MGLY+ 0.35OPEN+ 0.65XOPN+ RO2+ HO2';
k(:,i)=6.50e-14.*exp(500./T) ;
Gstr{i,1}='XLO2'; Gstr{i,2}='RO2';
fXLO2(i)=fXLO2(i)-1; fRO2(i)=fRO2(i)-1; fGLY(i)=fGLY(i)+0.26; fMGLY(i)=fMGLY(i)+0.77; fOPEN(i)=fOPEN(i)+0.35; fXOPN(i)=fXOPN(i)+0.65; fRO2(i)=fRO2(i)+1; fHO2(i)=fHO2(i)+1;

%182
i=i+1;
Rnames{i} = 'CRES + OH = 0.732CAT1+ 0.2CRO+ 0.025GLY+ 0.025OPEN+ 0.02XO2N+ 0.02RO2+ HO2';
k(:,i)=1.70e-12.*exp(950./T);
Gstr{i,1}='CRES'; Gstr{i,2}='OH';
fCRES(i)=fCRES(i)-1; fOH(i)=fOH(i)-1; fCAT1(i)=fCAT1(i)+0.732; fCRO(i)=fCRO(i)+0.2; fGLY(i)=fGLY(i)+0.025; fOPEN(i)=fOPEN(i)+0.025; fXO2N(i)=fXO2N(i)+0.02; fRO2(i)=fRO2(i)+0.02; fHO2(i)=fHO2(i)+1;

%183
i=i+1;
Rnames{i} = 'CRES + NO3 = 0.3CRO+ 0.24GLY+ HNO3+ 0.24MGLY+ 0.48OPO3+ 0.48XO2+ 0.12XO2H+ 0.1XO2N+ 0.7RO2';
k(:,i)=1.4e-11;
Gstr{i,1}='CRES'; Gstr{i,2}='NO3';
fCRES(i)=fCRES(i)-1; fNO3(i)=fNO3(i)-1; fCRO(i)=fCRO(i)+0.3; fGLY(i)=fGLY(i)+0.24; fHNO3(i)=fHNO3(i)+1; fMGLY(i)=fMGLY(i)+0.24; fOPO3(i)=fOPO3(i)+0.48; fXO2(i)=fXO2(i)+0.48; fXO2H(i)=fXO2H(i)+0.12; fXO2N(i)=fXO2N(i)+0.1; fRO2(i)=fRO2(i)+0.7;

%184
i=i+1;
Rnames{i} = 'CRO + NO2 = CRON';
k(:,i)=2.1e-12;
Gstr{i,1}='CRO'; Gstr{i,2}='NO2';
fCRO(i)=fCRO(i)-1; fNO2(i)=fNO2(i)-1; fCRON(i)=fCRON(i)+1;

%185
i=i+1;
Rnames{i} = 'CRO + HO2 = CRES';
k(:,i)=5.5e-12;
Gstr{i,1}='CRO'; Gstr{i,2}='HO2';
fCRO(i)=fCRO(i)-1; fHO2(i)=fHO2(i)-1; fCRES(i)=fCRES(i)+1;

%186
i=i+1;
Rnames{i} = 'CRON + OH = 0.5CRO+ NTR2';
k(:,i)=1.53e-12;
Gstr{i,1}='CRON'; Gstr{i,2}='OH';
fCRON(i)=fCRON(i)-1; fOH(i)=fOH(i)-1; fCRO(i)=fCRO(i)+0.5; fNTR2(i)=fNTR2(i)+1;

%187
i=i+1;
Rnames{i} = 'CRON + NO3 = 0.5CRO+ HNO3+ NTR2';
k(:,i)=3.8e-12;
Gstr{i,1}='CRON'; Gstr{i,2}='NO3';
fCRON(i)=fCRON(i)-1; fNO3(i)=fNO3(i)-1; fCRO(i)=fCRO(i)+0.5; fHNO3(i)=fHNO3(i)+1; fNTR2(i)=fNTR2(i)+1;

%188
i=i+1;
Rnames{i} = 'CRON = FORM+ HONO+ OPEN+ HO2';
k(:,i)=JCRON;
Gstr{i,1}='CRON';
fCRON(i)=fCRON(i)-1; fFORM(i)=fFORM(i)+1; fHONO(i)=fHONO(i)+1; fOPEN(i)=fOPEN(i)+1; fHO2(i)=fHO2(i)+1;

%189
i=i+1;
Rnames{i} = 'XOPN = 0.3C2O3+ 0.4GLY+ XO2H+ 0.7CO+ 0.7HO2';
k(:,i)=0.05.*JNO2;
Gstr{i,1}='XOPN';
fXOPN(i)=fXOPN(i)-1; fC2O3(i)=fC2O3(i)+0.3; fGLY(i)=fGLY(i)+0.4; fXO2H(i)=fXO2H(i)+1; fCO(i)=fCO(i)+0.7; fHO2(i)=fHO2(i)+0.7;

%190
i=i+1;
Rnames{i} = 'XOPN + OH = 0.4GLY+ MGLY+ RO2+ XO2H+ XO2H+ RO2';
k(:,i)=9e-11;
Gstr{i,1}='XOPN'; Gstr{i,2}='OH';
fXOPN(i)=fXOPN(i)-1; fOH(i)=fOH(i)-1; fGLY(i)=fGLY(i)+0.4; fMGLY(i)=fMGLY(i)+1; fRO2(i)=fRO2(i)+1; fXO2H(i)=fXO2H(i)+1; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%191
i=i+1;
Rnames{i} = 'XOPN + O3 = 0.1ALD2+ 0.6C2O3+ 1.2MGLY+ 0.3XO2H+ 0.3RO2+ 0.5CO+ 0.5OH';
k(:,i)=1.08e-16.*exp(-500./T);
Gstr{i,1}='XOPN'; Gstr{i,2}='O3';
fXOPN(i)=fXOPN(i)-1; fO3(i)=fO3(i)-1; fALD2(i)=fALD2(i)+0.1; fC2O3(i)=fC2O3(i)+0.6; fMGLY(i)=fMGLY(i)+1.2; fXO2H(i)=fXO2H(i)+0.3; fRO2(i)=fRO2(i)+0.3; fCO(i)=fCO(i)+0.5; fOH(i)=fOH(i)+0.5;

%192
i=i+1;
Rnames{i} = 'XOPN + NO3 = 0.25MGLY+ 0.5NTR2+ 0.25OPEN+ 0.45XO2+ 0.45XO2H+ 0.1XO2N+ 0.5NO2+ RO2';
k(:,i)=3e-12;
Gstr{i,1}='XOPN'; Gstr{i,2}='NO3';
fXOPN(i)=fXOPN(i)-1; fNO3(i)=fNO3(i)-1; fMGLY(i)=fMGLY(i)+0.25; fNTR2(i)=fNTR2(i)+0.5; fOPEN(i)=fOPEN(i)+0.25; fXO2(i)=fXO2(i)+0.45; fXO2H(i)=fXO2H(i)+0.45; fXO2N(i)=fXO2N(i)+0.1; fNO2(i)=fNO2(i)+0.5; fRO2(i)=fRO2(i)+1;

%193
i=i+1;
Rnames{i} = 'OPEN = OPO3+ CO+ HO2';
k(:,i)=0.028.*JNO2;
Gstr{i,1}='OPEN';
fOPEN(i)=fOPEN(i)-1; fOPO3(i)=fOPO3(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%194
i=i+1;
Rnames{i} = 'OPEN + OH = 0.4GLY+ 0.6OPO3+ 0.4XO2H+ 0.4RO2';
k(:,i)=4.4e-11;
Gstr{i,1}='OPEN'; Gstr{i,2}='OH';
fOPEN(i)=fOPEN(i)-1; fOH(i)=fOH(i)-1; fGLY(i)=fGLY(i)+0.4; fOPO3(i)=fOPO3(i)+0.6; fXO2H(i)=fXO2H(i)+0.4; fRO2(i)=fRO2(i)+0.4;

%195
i=i+1;
Rnames{i} = 'OPEN + O3 = 0.02ALD2+ 0.12C2O3+ 0.08FORM+ 1.4GLY+ 0.24MGLY+ 1.98CO+ 0.5OH+ 0.56HO2';
k(:,i)=5.40e-17.*exp(-500./T);
Gstr{i,1}='OPEN'; Gstr{i,2}='O3';
fOPEN(i)=fOPEN(i)-1; fO3(i)=fO3(i)-1; fALD2(i)=fALD2(i)+0.02; fC2O3(i)=fC2O3(i)+0.12; fFORM(i)=fFORM(i)+0.08; fGLY(i)=fGLY(i)+1.4; fMGLY(i)=fMGLY(i)+0.24; fCO(i)=fCO(i)+1.98; fOH(i)=fOH(i)+0.5; fHO2(i)=fHO2(i)+0.56;

%196
i=i+1;
Rnames{i} = 'OPEN + NO3 = HNO3+ OPO3';
k(:,i)=3.8e-12;
Gstr{i,1}='OPEN'; Gstr{i,2}='NO3';
fOPEN(i)=fOPEN(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+1; fOPO3(i)=fOPO3(i)+1;

%197
i=i+1;
Rnames{i} = 'CAT1 + OH = 0.5CRO+ 0.14FORM+ 0.2HO2';
k(:,i)=5e-11;
Gstr{i,1}='CAT1'; Gstr{i,2}='OH';
fCAT1(i)=fCAT1(i)-1; fOH(i)=fOH(i)-1; fCRO(i)=fCRO(i)+0.5; fFORM(i)=fFORM(i)+0.14; fHO2(i)=fHO2(i)+0.2;

%198
i=i+1;
Rnames{i} = 'CAT1 + NO3 = CRO+ HNO3';
k(:,i)=1.7e-10;
Gstr{i,1}='CAT1'; Gstr{i,2}='NO3';
fCAT1(i)=fCAT1(i)-1; fNO3(i)=fNO3(i)-1; fCRO(i)=fCRO(i)+1; fHNO3(i)=fHNO3(i)+1;

%199
i=i+1;
Rnames{i} = 'NO + OPO3 = 0.2CXO3+ 0.5GLY+ NO2+ 0.5CO+ 0.8HO2';
k(:,i)=6.70e-12.*exp(340./T);
Gstr{i,1}='NO'; Gstr{i,2}='OPO3';
fNO(i)=fNO(i)-1; fOPO3(i)=fOPO3(i)-1; fCXO3(i)=fCXO3(i)+0.2; fGLY(i)=fGLY(i)+0.5; fNO2(i)=fNO2(i)+1; fCO(i)=fCO(i)+0.5; fHO2(i)=fHO2(i)+0.8;

%200
i=i+1;
Rnames{i} = 'OPO3 + NO2 = OPAN';
k(:,i)=K_C2O3_NO2./1.19e0;
Gstr{i,1}='OPO3'; Gstr{i,2}='NO2';
fOPO3(i)=fOPO3(i)-1; fNO2(i)=fNO2(i)-1; fOPAN(i)=fOPAN(i)+1;

%201
i=i+1;
Rnames{i} = 'OPAN = OPO3+ NO2';
k(:,i)=K_PAN./1.19e0;
Gstr{i,1}='OPAN';
fOPAN(i)=fOPAN(i)-1; fOPO3(i)=fOPO3(i)+1; fNO2(i)=fNO2(i)+1;

%202
i=i+1;
Rnames{i} = 'OPO3 + HO2 = 0.13AACD+ 0.5MEO2+ 0.37PACD+ 0.13O3+ 0.5RO2+ 0.5OH';
k(:,i)=3.14e-12.*exp(580./T);
Gstr{i,1}='OPO3'; Gstr{i,2}='HO2';
fOPO3(i)=fOPO3(i)-1; fHO2(i)=fHO2(i)-1; fAACD(i)=fAACD(i)+0.13; fMEO2(i)=fMEO2(i)+0.5; fPACD(i)=fPACD(i)+0.37; fO3(i)=fO3(i)+0.13; fRO2(i)=fRO2(i)+0.5; fOH(i)=fOH(i)+0.5;

%203
i=i+1;
Rnames{i} = 'C2O3 + OPO3 = ALDX+ MEO2+ RO2+ XO2+ RO2';
k(:,i)=2.90e-12.*exp(500./T);
Gstr{i,1}='C2O3'; Gstr{i,2}='OPO3';
fC2O3(i)=fC2O3(i)-1; fOPO3(i)=fOPO3(i)-1; fALDX(i)=fALDX(i)+1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1; fXO2(i)=fXO2(i)+1; fRO2(i)=fRO2(i)+1;

%204
i=i+1;
Rnames{i} = 'OPO3 + RO2 = 0.2AACD+ 0.8ALDX+ 0.8XO2H+ 1.8RO2';
k(:,i)=4.40e-13.*exp(1070./T);
Gstr{i,1}='OPO3'; Gstr{i,2}='RO2';
fOPO3(i)=fOPO3(i)-1; fRO2(i)=fRO2(i)-1; fAACD(i)=fAACD(i)+0.2; fALDX(i)=fALDX(i)+0.8; fXO2H(i)=fXO2H(i)+0.8; fRO2(i)=fRO2(i)+1.8;

%205
i=i+1;
Rnames{i} = 'OPAN + OH = 0.5GLY+ 0.5NTR2+ 0.5NO2+ CO';
k(:,i)=3.6e-11;
Gstr{i,1}='OPAN'; Gstr{i,2}='OH';
fOPAN(i)=fOPAN(i)-1; fOH(i)=fOH(i)-1; fGLY(i)=fGLY(i)+0.5; fNTR2(i)=fNTR2(i)+0.5; fNO2(i)=fNO2(i)+0.5; fCO(i)=fCO(i)+1;

%206
i=i+1;
Rnames{i} = 'PANX + OH = ALD2+ NO2';
k(:,i)=3e-12;
Gstr{i,1}='PANX'; Gstr{i,2}='OH';
fPANX(i)=fPANX(i)-1; fOH(i)=fOH(i)-1; fALD2(i)=fALD2(i)+1; fNO2(i)=fNO2(i)+1;

%207
i=i+1;
Rnames{i} = 'NTR2 = HNO3';
k(:,i)=2.3e-5;
Gstr{i,1}='NTR2';
fNTR2(i)=fNTR2(i)-1; fHNO3(i)=fHNO3(i)+1;

%208
i=i+1;
Rnames{i} = 'ECH4 + OH = MEO2+ RO2';
k(:,i)=1.85e-12.*exp(-1690./T);
Gstr{i,1}='ECH4'; Gstr{i,2}='OH';
fECH4(i)=fECH4(i)-1; fOH(i)=fOH(i)-1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1;

%209
i=i+1;
Rnames{i} = 'I2 = I+ I';
k(:,i)=JI2;
Gstr{i,1}='I2';
fI2(i)=fI2(i)-1; fI(i)=fI(i)+1; fI(i)=fI(i)+1;

%210
i=i+1;
Rnames{i} = 'HOI = I+ OH';
k(:,i)=JHOI;
Gstr{i,1}='HOI';
fHOI(i)=fHOI(i)-1; fI(i)=fI(i)+1; fOH(i)=fOH(i)+1;

%211
i=i+1;
Rnames{i} = 'I + O3 = IO';
k(:,i)=2.10e-11.*exp(-830./T);
Gstr{i,1}='I'; Gstr{i,2}='O3';
fI(i)=fI(i)-1; fO3(i)=fO3(i)-1; fIO(i)=fIO(i)+1;

%212
i=i+1;
Rnames{i} = 'IO = I+ O';
k(:,i)=JIO ;
Gstr{i,1}='IO';
fIO(i)=fIO(i)-1; fI(i)=fI(i)+1; fO(i)=fO(i)+1;

%213
i=i+1;
Rnames{i} = 'IO + IO = 0.4I+ 0.6I2O2+ 0.4OIO';
k(:,i)=5.40e-11.*exp(180./T);
Gstr{i,1}='IO'; Gstr{i,2}='IO';
fIO(i)=fIO(i)-1; fIO(i)=fIO(i)-1; fI(i)=fI(i)+0.4; fI2O2(i)=fI2O2(i)+0.6; fOIO(i)=fOIO(i)+0.4;

%214
i=i+1;
Rnames{i} = 'IO + HO2 = HOI';
k(:,i)=1.40e-11.*exp(540./T);
Gstr{i,1}='IO'; Gstr{i,2}='HO2';
fIO(i)=fIO(i)-1; fHO2(i)=fHO2(i)-1; fHOI(i)=fHOI(i)+1;

%215
i=i+1;
Rnames{i} = 'IO + NO = I+ NO2';
k(:,i)=7.15e-12.*exp(300./T);
Gstr{i,1}='IO'; Gstr{i,2}='NO';
fIO(i)=fIO(i)-1; fNO(i)=fNO(i)-1; fI(i)=fI(i)+1; fNO2(i)=fNO2(i)+1;

%216
i=i+1;
Rnames{i} = 'IO + NO2 = INO3';
k(:,i)=K_IO_NO2;
Gstr{i,1}='IO'; Gstr{i,2}='NO2';
fIO(i)=fIO(i)-1; fNO2(i)=fNO2(i)-1; fINO3(i)=fINO3(i)+1;

%217
i=i+1;
Rnames{i} = 'OIO = I';
k(:,i)=JOIO;
Gstr{i,1}='OIO';
fOIO(i)=fOIO(i)-1; fI(i)=fI(i)+1;

%218
i=i+1;
Rnames{i} = 'OIO + OH = HIO3';
k(:,i)=K_OIO_OH;
Gstr{i,1}='OIO'; Gstr{i,2}='OH';
fOIO(i)=fOIO(i)-1; fOH(i)=fOH(i)-1; fHIO3(i)=fHIO3(i)+1;

%219
i=i+1;
Rnames{i} = 'IO + OIO = IXOY';
k(:,i)=1e-10;
Gstr{i,1}='IO'; Gstr{i,2}='OIO';
fIO(i)=fIO(i)-1; fOIO(i)=fOIO(i)-1; fIXOY(i)=fIXOY(i)+1;

%220
i=i+1;
Rnames{i} = 'NO + OIO = IO+ NO2';
k(:,i)=1.10e-12.*exp(542./T);
Gstr{i,1}='NO'; Gstr{i,2}='OIO';
fNO(i)=fNO(i)-1; fOIO(i)=fOIO(i)-1; fIO(i)=fIO(i)+1; fNO2(i)=fNO2(i)+1;

%221
i=i+1;
Rnames{i} = 'I2O2 = I+ OIO';
k(:,i)=1.0e1;
Gstr{i,1}='I2O2';
fI2O2(i)=fI2O2(i)-1; fI(i)=fI(i)+1; fOIO(i)=fOIO(i)+1;

%222
i=i+1;
Rnames{i} = 'I2O2 + O3 = IXOY';
k(:,i)=1.0e-12;
Gstr{i,1}='I2O2'; Gstr{i,2}='O3';
fI2O2(i)=fI2O2(i)-1; fO3(i)=fO3(i)-1; fIXOY(i)=fIXOY(i)+1;

%223
i=i+1;
Rnames{i} = 'INO3 = I+ NO3';
k(:,i)=JINO3;
Gstr{i,1}='INO3';
fINO3(i)=fINO3(i)-1; fI(i)=fI(i)+1; fNO3(i)=fNO3(i)+1;

%224
i=i+1;
Rnames{i} = 'INO3 = HNO3+ HOI';
k(:,i)=2.5e-22.*H2O;
Gstr{i,1}='INO3';
fINO3(i)=fINO3(i)-1; fHNO3(i)=fHNO3(i)+1; fHOI(i)=fHOI(i)+1;

%225
i=i+1;
Rnames{i} = 'XPRP = XO2N+ RO2';
k(:,i)=K_XPRP;
Gstr{i,1}='XPRP';
fXPRP(i)=fXPRP(i)-1; fXO2N(i)=fXO2N(i)+1; fRO2(i)=fRO2(i)+1;

%226
i=i+1;
Rnames{i} = 'XPRP = 0.732ACET+ 0.268ALDX+ 0.268PAR+ XO2H+ RO2';
k(:,i)=1.0e0;
Gstr{i,1}='XPRP';
fXPRP(i)=fXPRP(i)-1; fACET(i)=fACET(i)+0.732; fALDX(i)=fALDX(i)+0.268; fPAR(i)=fPAR(i)+0.268; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%227
i=i+1;
Rnames{i} = 'XPAR = XO2N+ RO2';
k(:,i)=K_XPAR;
Gstr{i,1}='XPAR';
fXPAR(i)=fXPAR(i)-1; fXO2N(i)=fXO2N(i)+1; fRO2(i)=fRO2(i)+1;

%228
i=i+1;
Rnames{i} = 'XPAR = 0.126ALDX- 0.126PAR+ 0.874ROR+ 0.874XO2+ 0.126XO2H+ RO2';
k(:,i)=1.0e0;
Gstr{i,1}='XPAR';
fXPAR(i)=fXPAR(i)-1; fALDX(i)=fALDX(i)+0.126; fPAR(i)=fPAR(i)-0.126; fROR(i)=fROR(i)+0.874; fXO2(i)=fXO2(i)+0.874; fXO2H(i)=fXO2H(i)+0.126; fRO2(i)=fRO2(i)+1;

%229
i=i+1;
Rnames{i} = 'INTR = HNO3';
k(:,i)=1.40e-4;
Gstr{i,1}='INTR';
fINTR(i)=fINTR(i)-1; fHNO3(i)=fHNO3(i)+1;

%230
i=i+1;
Rnames{i} = 'SO2 = SULF';
k(:,i)=0.0e0;
Gstr{i,1}='SO2';
fSO2(i)=fSO2(i)-1; fSULF(i)=fSULF(i)+1;

%231
i=i+1;
Rnames{i} = 'DMS + OH = FORM+ MEO2+ SO2+ RO2';
k(:,i)=1.12e-11.*exp(-250./T);
Gstr{i,1}='DMS'; Gstr{i,2}='OH';
fDMS(i)=fDMS(i)-1; fOH(i)=fOH(i)-1; fFORM(i)=fFORM(i)+1; fMEO2(i)=fMEO2(i)+1; fSO2(i)=fSO2(i)+1; fRO2(i)=fRO2(i)+1;

%232
i=i+1;
Rnames{i} = 'DMS + OH + O2 = MEO2+ SULF+ RO2';
k(:,i)=1.28e-37.*exp(4480./T);
Gstr{i,1}='DMS'; Gstr{i,2}='OH'; Gstr{i,3}='O2';
fDMS(i)=fDMS(i)-1; fOH(i)=fOH(i)-1; fO2(i)=fO2(i)-1; fMEO2(i)=fMEO2(i)+1; fSULF(i)=fSULF(i)+1; fRO2(i)=fRO2(i)+1;

%233
i=i+1;
Rnames{i} = 'DMS + NO3 = FORM+ HNO3+ MEO2+ SO2+ RO2';
k(:,i)=1.90e-13.*exp(520./T);
Gstr{i,1}='DMS'; Gstr{i,2}='NO3';
fDMS(i)=fDMS(i)-1; fNO3(i)=fNO3(i)-1; fFORM(i)=fFORM(i)+1; fHNO3(i)=fHNO3(i)+1; fMEO2(i)=fMEO2(i)+1; fSO2(i)=fSO2(i)+1; fRO2(i)=fRO2(i)+1;

%234
i=i+1;
Rnames{i} = 'NO2 + OH = H2O+ HNO3';
k(:,i)=1.1e-30.*H2O;
Gstr{i,1}='NO2'; Gstr{i,2}='OH';
fNO2(i)=fNO2(i)-1; fOH(i)=fOH(i)-1; fH2O(i)=fH2O(i)+1; fHNO3(i)=fHNO3(i)+1;

%235
i=i+1;
Rnames{i} = 'CL2 = CL+ CL';
k(:,i)=JCL2;
Gstr{i,1}='CL2';
fCL2(i)=fCL2(i)-1; fCL(i)=fCL(i)+1; fCL(i)=fCL(i)+1;

%236
i=i+1;
Rnames{i} = 'ICL = CL+ I';
k(:,i)=JICL;
Gstr{i,1}='ICL';
fICL(i)=fICL(i)-1; fCL(i)=fCL(i)+1; fI(i)=fI(i)+1;

%237
i=i+1;
Rnames{i} = 'HOCL = CL+ OH';
k(:,i)=JHOCL;
Gstr{i,1}='HOCL';
fHOCL(i)=fHOCL(i)-1; fCL(i)=fCL(i)+1; fOH(i)=fOH(i)+1;

%238
i=i+1;
Rnames{i} = 'CL + O3 = CLO';
k(:,i)=2.30e-11.*exp(-200./T);
Gstr{i,1}='CL'; Gstr{i,2}='O3';
fCL(i)=fCL(i)-1; fO3(i)=fO3(i)-1; fCLO(i)=fCLO(i)+1;

%239
i=i+1;
Rnames{i} = 'CL + HO2 = 0.22CLO+ 0.78HCL+ 0.22OH';
k(:,i)=3.00e-11.*exp(120./T);
Gstr{i,1}='CL'; Gstr{i,2}='HO2';
fCL(i)=fCL(i)-1; fHO2(i)=fHO2(i)-1; fCLO(i)=fCLO(i)+0.22; fHCL(i)=fHCL(i)+0.78; fOH(i)=fOH(i)+0.22;

%240
i=i+1;
Rnames{i} = 'CL + H2 = HCL+ HO2';
k(:,i)=3.05e-11.*exp(-2270./T);
Gstr{i,1}='CL'; Gstr{i,2}='H2';
fCL(i)=fCL(i)-1; fH2(i)=fH2(i)-1; fHCL(i)=fHCL(i)+1; fHO2(i)=fHO2(i)+1;

%241
i=i+1;
Rnames{i} = 'CLO + CLO = 1.4CL+ 0.3CL2';
k(:,i)=1.63e-14;
Gstr{i,1}='CLO'; Gstr{i,2}='CLO';
fCLO(i)=fCLO(i)-1; fCLO(i)=fCLO(i)-1; fCL(i)=fCL(i)+1.4; fCL2(i)=fCL2(i)+0.3;

%242
i=i+1;
Rnames{i} = 'CLO + IO = CL+ I';
k(:,i)=5.00e-13.*exp(300./T);
Gstr{i,1}='CLO'; Gstr{i,2}='IO';
fCLO(i)=fCLO(i)-1; fIO(i)=fIO(i)-1; fCL(i)=fCL(i)+1; fI(i)=fI(i)+1;

%243
i=i+1;
Rnames{i} = 'CLO + HO2 = HOCL';
k(:,i)=2.60e-12.*exp(290./T);
Gstr{i,1}='CLO'; Gstr{i,2}='HO2';
fCLO(i)=fCLO(i)-1; fHO2(i)=fHO2(i)-1; fHOCL(i)=fHOCL(i)+1;

%244
i=i+1;
Rnames{i} = 'CLO + MEO2 = CL+ FORM+ HO2';
k(:,i)=1.80e-11.*exp(-600./T);
Gstr{i,1}='CLO'; Gstr{i,2}='MEO2';
fCLO(i)=fCLO(i)-1; fMEO2(i)=fMEO2(i)-1; fCL(i)=fCL(i)+1; fFORM(i)=fFORM(i)+1; fHO2(i)=fHO2(i)+1;

%245
i=i+1;
Rnames{i} = 'CLO + NO = CL+ NO2';
k(:,i)=6.40e-12.*exp(290./T);
Gstr{i,1}='CLO'; Gstr{i,2}='NO';
fCLO(i)=fCLO(i)-1; fNO(i)=fNO(i)-1; fCL(i)=fCL(i)+1; fNO2(i)=fNO2(i)+1;

%246
i=i+1;
Rnames{i} = 'CLO + NO2 = CLN3';
k(:,i)=K_CLO_NO2;
Gstr{i,1}='CLO'; Gstr{i,2}='NO2';
fCLO(i)=fCLO(i)-1; fNO2(i)=fNO2(i)-1; fCLN3(i)=fCLN3(i)+1;

%247
i=i+1;
Rnames{i} = 'CLN3 = CLO+ NO2';
k(:,i)=K_CLO_NO2./(2.98e-28.*exp(13264./T));
Gstr{i,1}='CLN3';
fCLN3(i)=fCLN3(i)-1; fCLO(i)=fCLO(i)+1; fNO2(i)=fNO2(i)+1;

%248
i=i+1;
Rnames{i} = 'CLN3 = CLO+ NO2';
k(:,i)=JCLNO3_CLO_NO2;
Gstr{i,1}='CLN3';
fCLN3(i)=fCLN3(i)-1; fCLO(i)=fCLO(i)+1; fNO2(i)=fNO2(i)+1;

%249
i=i+1;
Rnames{i} = 'CLN3 = CL+ NO3';
k(:,i)=JCLNO3_Cl_NO3;
Gstr{i,1}='CLN3';
fCLN3(i)=fCLN3(i)-1; fCL(i)=fCL(i)+1; fNO3(i)=fNO3(i)+1;

%250
i=i+1;
Rnames{i} = 'CLN2 = CL+ NO2';
k(:,i)=JCLNO2;
Gstr{i,1}='CLN2';
fCLN2(i)=fCLN2(i)-1; fCL(i)=fCL(i)+1; fNO2(i)=fNO2(i)+1;

%251
i=i+1;
Rnames{i} = 'HCL + N2O5 = CLN2+ HNO3';
k(:,i)=6e-13;
Gstr{i,1}='HCL'; Gstr{i,2}='N2O5';
fHCL(i)=fHCL(i)-1; fN2O5(i)=fN2O5(i)-1; fCLN2(i)=fCLN2(i)+1; fHNO3(i)=fHNO3(i)+1;

%252
i=i+1;
Rnames{i} = 'CLN3 = HNO3+ HOCL';
k(:,i)=2.5e-22.*H2O;
Gstr{i,1}='CLN3';
fCLN3(i)=fCLN3(i)-1; fHNO3(i)=fHNO3(i)+1; fHOCL(i)=fHOCL(i)+1;

%253
i=i+1;
Rnames{i} = 'CL + FORM = HCL+ CO+ HO2';
k(:,i)=8.10e-11.*exp(-30./T);
Gstr{i,1}='CL'; Gstr{i,2}='FORM';
fCL(i)=fCL(i)-1; fFORM(i)=fFORM(i)-1; fHCL(i)=fHCL(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%254
i=i+1;
Rnames{i} = 'ALD2 + CL = C2O3+ HCL';
k(:,i)=7.3e-11;
Gstr{i,1}='ALD2'; Gstr{i,2}='CL';
fALD2(i)=fALD2(i)-1; fCL(i)=fCL(i)-1; fC2O3(i)=fC2O3(i)+1; fHCL(i)=fHCL(i)+1;

%255
i=i+1;
Rnames{i} = 'ALDX + CL = CXO3+ HCL';
k(:,i)=1.4e-10;
Gstr{i,1}='ALDX'; Gstr{i,2}='CL';
fALDX(i)=fALDX(i)-1; fCL(i)=fCL(i)-1; fCXO3(i)=fCXO3(i)+1; fHCL(i)=fHCL(i)+1;

%256
i=i+1;
Rnames{i} = 'CL + GLY = HCL+ 0.2XO2+ 0.2RO2+ 1.8CO+ HO2';
k(:,i)=3.8e-11;
Gstr{i,1}='CL'; Gstr{i,2}='GLY';
fCL(i)=fCL(i)-1; fGLY(i)=fGLY(i)-1; fHCL(i)=fHCL(i)+1; fXO2(i)=fXO2(i)+0.2; fRO2(i)=fRO2(i)+0.2; fCO(i)=fCO(i)+1.8; fHO2(i)=fHO2(i)+1;

%257
i=i+1;
Rnames{i} = 'CL + GLYD = 0.8C2O3+ 0.2GLY+ HCL+ 0.2HO2';
k(:,i)=6.6e-11;
Gstr{i,1}='CL'; Gstr{i,2}='GLYD';
fCL(i)=fCL(i)-1; fGLYD(i)=fGLYD(i)-1; fC2O3(i)=fC2O3(i)+0.8; fGLY(i)=fGLY(i)+0.2; fHCL(i)=fHCL(i)+1; fHO2(i)=fHO2(i)+0.2;

%258
i=i+1;
Rnames{i} = 'CL + MGLY = C2O3+ HCL+ CO';
k(:,i)=4.8e-11;
Gstr{i,1}='CL'; Gstr{i,2}='MGLY';
fCL(i)=fCL(i)-1; fMGLY(i)=fMGLY(i)-1; fC2O3(i)=fC2O3(i)+1; fHCL(i)=fHCL(i)+1; fCO(i)=fCO(i)+1;

%259
i=i+1;
Rnames{i} = 'ACET + CL = C2O3+ FORM+ HCL+ XO2+ RO2';
k(:,i)=1.63e-11.*exp(-610./T);
Gstr{i,1}='ACET'; Gstr{i,2}='CL';
fACET(i)=fACET(i)-1; fCL(i)=fCL(i)-1; fC2O3(i)=fC2O3(i)+1; fFORM(i)=fFORM(i)+1; fHCL(i)=fHCL(i)+1; fXO2(i)=fXO2(i)+1; fRO2(i)=fRO2(i)+1;

%260
i=i+1;
Rnames{i} = 'CL + KET = 0.5ALD2+ 0.5C2O3+ 0.5CXO3+ HCL+ 0.5MEO2- 2.5PAR+ 0.5XO2H+ RO2';
k(:,i)=2.77e-11.*exp(76./T);
Gstr{i,1}='CL'; Gstr{i,2}='KET';
fCL(i)=fCL(i)-1; fKET(i)=fKET(i)-1; fALD2(i)=fALD2(i)+0.5; fC2O3(i)=fC2O3(i)+0.5; fCXO3(i)=fCXO3(i)+0.5; fHCL(i)=fHCL(i)+1; fMEO2(i)=fMEO2(i)+0.5; fPAR(i)=fPAR(i)-2.5; fXO2H(i)=fXO2H(i)+0.5; fRO2(i)=fRO2(i)+1;

%261
i=i+1;
Rnames{i} = 'CL + MEOH = FORM+ HCL+ HO2';
k(:,i)=5.5e-11;
Gstr{i,1}='CL'; Gstr{i,2}='MEOH';
fCL(i)=fCL(i)-1; fMEOH(i)=fMEOH(i)-1; fFORM(i)=fFORM(i)+1; fHCL(i)=fHCL(i)+1; fHO2(i)=fHO2(i)+1;

%262
i=i+1;
Rnames{i} = 'CL + ETOH = 0.95ALD2+ 0.078FORM+ 0.011GLYD+ 0.1XO2H+ 0.1RO2+ 0.9HO2';
k(:,i)=9.6e-11;
Gstr{i,1}='CL'; Gstr{i,2}='ETOH';
fCL(i)=fCL(i)-1; fETOH(i)=fETOH(i)-1; fALD2(i)=fALD2(i)+0.95; fFORM(i)=fFORM(i)+0.078; fGLYD(i)=fGLYD(i)+0.011; fXO2H(i)=fXO2H(i)+0.1; fRO2(i)=fRO2(i)+0.1; fHO2(i)=fHO2(i)+0.9;

%263
i=i+1;
Rnames{i} = 'CL + ISPD = 0.48C2O3+ 0.5CLAD+ 0.34CLAO+ 0.17HCL+ 0.17OPO3+ 0.48XO2+ 0.32XO2H+ 0.04XO2N+ 0.84RO2+ 0.32CO';
k(:,i)=2.20e-10;
Gstr{i,1}='CL'; Gstr{i,2}='ISPD';
fCL(i)=fCL(i)-1; fISPD(i)=fISPD(i)-1; fC2O3(i)=fC2O3(i)+0.48; fCLAD(i)=fCLAD(i)+0.5; fCLAO(i)=fCLAO(i)+0.34; fHCL(i)=fHCL(i)+0.17; fOPO3(i)=fOPO3(i)+0.17; fXO2(i)=fXO2(i)+0.48; fXO2H(i)=fXO2H(i)+0.32; fXO2N(i)=fXO2N(i)+0.04; fRO2(i)=fRO2(i)+0.84; fCO(i)=fCO(i)+0.32;

%264
i=i+1;
Rnames{i} = 'FMCL = HCL+ CO';
k(:,i)=6.94e-5;
Gstr{i,1}='FMCL';
fFMCL(i)=fFMCL(i)-1; fHCL(i)=fHCL(i)+1; fCO(i)=fCO(i)+1;

%265
i=i+1;
Rnames{i} = 'CLAD = CL+ MEO2+ RO2+ XO2+ RO2+ CO';
k(:,i)=JCLAD;
Gstr{i,1}='CLAD';
fCLAD(i)=fCLAD(i)-1; fCL(i)=fCL(i)+1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1; fXO2(i)=fXO2(i)+1; fRO2(i)=fRO2(i)+1; fCO(i)=fCO(i)+1;

%266
i=i+1;
Rnames{i} = 'CLAD + OH = FMCL+ RO2+ XO2+ XO2H+ RO2';
k(:,i)=3.10e-12;
Gstr{i,1}='CLAD'; Gstr{i,2}='OH';
fCLAD(i)=fCLAD(i)-1; fOH(i)=fOH(i)-1; fFMCL(i)=fFMCL(i)+1; fRO2(i)=fRO2(i)+1; fXO2(i)=fXO2(i)+1; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%267
i=i+1;
Rnames{i} = 'CLAO + OH = C2O3+ FMCL+ RO2+ XO2+ RO2';
k(:,i)=4.20e-13;
Gstr{i,1}='CLAO'; Gstr{i,2}='OH';
fCLAO(i)=fCLAO(i)-1; fOH(i)=fOH(i)-1; fC2O3(i)=fC2O3(i)+1; fFMCL(i)=fFMCL(i)+1; fRO2(i)=fRO2(i)+1; fXO2(i)=fXO2(i)+1; fRO2(i)=fRO2(i)+1;

%268
i=i+1;
Rnames{i} = 'CH4 + CL = HCL+ MEO2+ RO2';
k(:,i)=7.10e-12.*exp(-1270./T);
Gstr{i,1}='CH4'; Gstr{i,2}='CL';
fCH4(i)=fCH4(i)-1; fCL(i)=fCL(i)-1; fHCL(i)=fHCL(i)+1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1;

%269
i=i+1;
Rnames{i} = 'CL + ECH4 = HCL+ MEO2+ RO2';
k(:,i)=7.10e-12.*exp(-1270./T);
Gstr{i,1}='CL'; Gstr{i,2}='ECH4';
fCL(i)=fCL(i)-1; fECH4(i)=fECH4(i)-1; fHCL(i)=fHCL(i)+1; fMEO2(i)=fMEO2(i)+1; fRO2(i)=fRO2(i)+1;

%270
i=i+1;
Rnames{i} = 'CL + ETHA = 0.991ALD2+ HCL+ 0.991XO2H+ 0.009XO2N+ RO2';
k(:,i)=7.20e-11.*exp(-70./T);
Gstr{i,1}='CL'; Gstr{i,2}='ETHA';
fCL(i)=fCL(i)-1; fETHA(i)=fETHA(i)-1; fALD2(i)=fALD2(i)+0.991; fHCL(i)=fHCL(i)+1; fXO2H(i)=fXO2H(i)+0.991; fXO2N(i)=fXO2N(i)+0.009; fRO2(i)=fRO2(i)+1;

%271
i=i+1;
Rnames{i} = 'CL + PRPA = HCL+ XPRP';
k(:,i)=1.4e-10;
Gstr{i,1}='CL'; Gstr{i,2}='PRPA';
fCL(i)=fCL(i)-1; fPRPA(i)=fPRPA(i)-1; fHCL(i)=fHCL(i)+1; fXPRP(i)=fXPRP(i)+1;

%272
i=i+1;
Rnames{i} = 'CL + PAR = HCL+ XPAR';
k(:,i)=4.5e-11;
Gstr{i,1}='CL'; Gstr{i,2}='PAR';
fCL(i)=fCL(i)-1; fPAR(i)=fPAR(i)-1; fHCL(i)=fHCL(i)+1; fXPAR(i)=fXPAR(i)+1;

%273
i=i+1;
Rnames{i} = 'CL + ETHY = 0.21CL+ 0.26FMCL+ 0.21GLY+ 0.53HCL+ 1.32CO+ 0.79HO2';
k(:,i)=K_ETHY_CL;
Gstr{i,1}='CL'; Gstr{i,2}='ETHY';
fCL(i)=fCL(i)-1; fETHY(i)=fETHY(i)-1; fCL(i)=fCL(i)+0.21; fFMCL(i)=fFMCL(i)+0.26; fGLY(i)=fGLY(i)+0.21; fHCL(i)=fHCL(i)+0.53; fCO(i)=fCO(i)+1.32; fHO2(i)=fHO2(i)+0.79;

%274
i=i+1;
Rnames{i} = 'CL + ETH = CLAD+ XO2H+ RO2';
k(:,i)=K_ETH_CL;
Gstr{i,1}='CL'; Gstr{i,2}='ETH';
fCL(i)=fCL(i)-1; fETH(i)=fETH(i)-1; fCLAD(i)=fCLAD(i)+1; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%275
i=i+1;
Rnames{i} = 'CL + OLE = 0.18ALDX+ 0.3CLAD+ 0.5CLAO+ 0.2HCL- PAR+ 0.92XO2H+ 0.08XO2N+ RO2';
k(:,i)=2.15e-10;
Gstr{i,1}='CL'; Gstr{i,2}='OLE';
fCL(i)=fCL(i)-1; fOLE(i)=fOLE(i)-1; fALDX(i)=fALDX(i)+0.18; fCLAD(i)=fCLAD(i)+0.3; fCLAO(i)=fCLAO(i)+0.5; fHCL(i)=fHCL(i)+0.2; fPAR(i)=fPAR(i)-1.0; fXO2H(i)=fXO2H(i)+0.92; fXO2N(i)=fXO2N(i)+0.08; fRO2(i)=fRO2(i)+1;

%276
i=i+1;
Rnames{i} = 'CL + IOLE = 0.17ALD2+ 0.23ALDX+ 0.17C2O3+ 0.56CLAO+ 0.44HCL+ 0.73XO2H+ 0.1XO2N+ 0.83RO2';
k(:,i)=3.50e-10;
Gstr{i,1}='CL'; Gstr{i,2}='IOLE';
fCL(i)=fCL(i)-1; fIOLE(i)=fIOLE(i)-1; fALD2(i)=fALD2(i)+0.17; fALDX(i)=fALDX(i)+0.23; fC2O3(i)=fC2O3(i)+0.17; fCLAO(i)=fCLAO(i)+0.56; fHCL(i)=fHCL(i)+0.44; fXO2H(i)=fXO2H(i)+0.73; fXO2N(i)=fXO2N(i)+0.1; fRO2(i)=fRO2(i)+0.83;

%277
i=i+1;
Rnames{i} = 'CL + ISOP = 0.58CLAD+ 0.22CLAO+ 0.05FMCL+ 0.1FORM+ 0.15HCL+ 0.15ISPD+ 0.88XO2+ 0.88XO2H+ 0.12XO2N+ 1.88RO2';
k(:,i)=7.60e-11.*exp(500./T);
Gstr{i,1}='CL'; Gstr{i,2}='ISOP';
fCL(i)=fCL(i)-1; fISOP(i)=fISOP(i)-1; fCLAD(i)=fCLAD(i)+0.58; fCLAO(i)=fCLAO(i)+0.22; fFMCL(i)=fFMCL(i)+0.05; fFORM(i)=fFORM(i)+0.1; fHCL(i)=fHCL(i)+0.15; fISPD(i)=fISPD(i)+0.15; fXO2(i)=fXO2(i)+0.88; fXO2H(i)=fXO2H(i)+0.88; fXO2N(i)=fXO2N(i)+0.12; fRO2(i)=fRO2(i)+1.88;

%278
i=i+1;
Rnames{i} = 'CL + TERP = 0.15CLAD+ 0.15CLAO+ 0.15FMCL+ 0.55HCL+ 0.7XO2H+ 0.3XO2N+ RO2';
k(:,i)=5.3e-10;
Gstr{i,1}='CL'; Gstr{i,2}='TERP';
fCL(i)=fCL(i)-1; fTERP(i)=fTERP(i)-1; fCLAD(i)=fCLAD(i)+0.15; fCLAO(i)=fCLAO(i)+0.15; fFMCL(i)=fFMCL(i)+0.15; fHCL(i)=fHCL(i)+0.55; fXO2H(i)=fXO2H(i)+0.7; fXO2N(i)=fXO2N(i)+0.3; fRO2(i)=fRO2(i)+1;

%279
i=i+1;
Rnames{i} = 'CL + TOL = 0.9CRES+ HCL+ 0.9XO2H+ 0.1XO2N+ RO2';
k(:,i)=5.6e-11;
Gstr{i,1}='CL'; Gstr{i,2}='TOL';
fCL(i)=fCL(i)-1; fTOL(i)=fTOL(i)-1; fCRES(i)=fCRES(i)+0.9; fHCL(i)=fHCL(i)+1; fXO2H(i)=fXO2H(i)+0.9; fXO2N(i)=fXO2N(i)+0.1; fRO2(i)=fRO2(i)+1;

%280
i=i+1;
Rnames{i} = 'CL + XYL = 0.9CRES+ HCL+ 0.9XO2H+ 0.1XO2N+ RO2';
k(:,i)=1.4e-10;
Gstr{i,1}='CL'; Gstr{i,2}='XYL';
fCL(i)=fCL(i)-1; fXYL(i)=fXYL(i)-1; fCRES(i)=fCRES(i)+0.9; fHCL(i)=fHCL(i)+1; fXO2H(i)=fXO2H(i)+0.9; fXO2N(i)=fXO2N(i)+0.1; fRO2(i)=fRO2(i)+1;

%281
i=i+1;
Rnames{i} = 'CL + CRES = 0.732CAT1+ 0.2CRO+ 0.025GLY+ HCL+ 0.025OPEN+ 0.02XO2N+ 0.02RO2+ HO2';
k(:,i)=1.9e-10;
Gstr{i,1}='CL'; Gstr{i,2}='CRES';
fCL(i)=fCL(i)-1; fCRES(i)=fCRES(i)-1; fCAT1(i)=fCAT1(i)+0.732; fCRO(i)=fCRO(i)+0.2; fGLY(i)=fGLY(i)+0.025; fHCL(i)=fHCL(i)+1; fOPEN(i)=fOPEN(i)+0.025; fXO2N(i)=fXO2N(i)+0.02; fRO2(i)=fRO2(i)+0.02; fHO2(i)=fHO2(i)+1;

%282
i=i+1;
Rnames{i} = 'CL + DMS = FORM+ HCL+ MEO2+ SO2+ RO2';
k(:,i)=1.8e-10;
Gstr{i,1}='CL'; Gstr{i,2}='DMS';
fCL(i)=fCL(i)-1; fDMS(i)=fDMS(i)-1; fFORM(i)=fFORM(i)+1; fHCL(i)=fHCL(i)+1; fMEO2(i)=fMEO2(i)+1; fSO2(i)=fSO2(i)+1; fRO2(i)=fRO2(i)+1;

%283
i=i+1;
Rnames{i} = 'BR2 = BR+ BR';
k(:,i)=JBR2;
Gstr{i,1}='BR2';
fBR2(i)=fBR2(i)-1; fBR(i)=fBR(i)+1; fBR(i)=fBR(i)+1;

%284
i=i+1;
Rnames{i} = 'IBR = BR+ I';
k(:,i)=JIBR;
Gstr{i,1}='IBR';
fIBR(i)=fIBR(i)-1; fBR(i)=fBR(i)+1; fI(i)=fI(i)+1;

%285
i=i+1;
Rnames{i} = 'BRCL = BR+ CL';
k(:,i)=JBRCL;
Gstr{i,1}='BRCL';
fBRCL(i)=fBRCL(i)-1; fBR(i)=fBR(i)+1; fCL(i)=fCL(i)+1;

%286
i=i+1;
Rnames{i} = 'HOBR = BR+ OH';
k(:,i)=JHOBR;
Gstr{i,1}='HOBR';
fHOBR(i)=fHOBR(i)-1; fBR(i)=fBR(i)+1; fOH(i)=fOH(i)+1;

%287
i=i+1;
Rnames{i} = 'BR + O3 = BRO';
k(:,i)=1.60e-11.*exp(-780./T);
Gstr{i,1}='BR'; Gstr{i,2}='O3';
fBR(i)=fBR(i)-1; fO3(i)=fO3(i)-1; fBRO(i)=fBRO(i)+1;

%288
i=i+1;
Rnames{i} = 'BR + HO2 = HBR';
k(:,i)=4.80e-12.*exp(-310./T);
Gstr{i,1}='BR'; Gstr{i,2}='HO2';
fBR(i)=fBR(i)-1; fHO2(i)=fHO2(i)-1; fHBR(i)=fHBR(i)+1;

%289
i=i+1;
Rnames{i} = 'BR + NO2 = BRN2';
k(:,i)=K_BR_NO2;
Gstr{i,1}='BR'; Gstr{i,2}='NO2';
fBR(i)=fBR(i)-1; fNO2(i)=fNO2(i)-1; fBRN2(i)=fBRN2(i)+1;

%290
i=i+1;
Rnames{i} = 'BR + NO3 = BRO+ NO2';
k(:,i)=1.6e-11;
Gstr{i,1}='BR'; Gstr{i,2}='NO3';
fBR(i)=fBR(i)-1; fNO3(i)=fNO3(i)-1; fBRO(i)=fBRO(i)+1; fNO2(i)=fNO2(i)+1;

%291
i=i+1;
Rnames{i} = 'BR2 + OH = BR+ HOBR';
k(:,i)=2.10e-11.*exp(240./T);
Gstr{i,1}='BR2'; Gstr{i,2}='OH';
fBR2(i)=fBR2(i)-1; fOH(i)=fOH(i)-1; fBR(i)=fBR(i)+1; fHOBR(i)=fHOBR(i)+1;

%292
i=i+1;
Rnames{i} = 'HBR + OH = BR';
k(:,i)=5.50e-12.*exp(200./T);
Gstr{i,1}='HBR'; Gstr{i,2}='OH';
fHBR(i)=fHBR(i)-1; fOH(i)=fOH(i)-1; fBR(i)=fBR(i)+1;

%293
i=i+1;
Rnames{i} = 'BRO = BR+ O';
k(:,i)=JBRO;
Gstr{i,1}='BRO';
fBRO(i)=fBRO(i)-1; fBR(i)=fBR(i)+1; fO(i)=fO(i)+1;

%294
i=i+1;
Rnames{i} = 'BRO + BRO = 1.7BR+ 0.15BR2';
k(:,i)=1.50e-12.*exp(230./T);
Gstr{i,1}='BRO'; Gstr{i,2}='BRO';
fBRO(i)=fBRO(i)-1; fBRO(i)=fBRO(i)-1; fBR(i)=fBR(i)+1.7; fBR2(i)=fBR2(i)+0.15;

%295
i=i+1;
Rnames{i} = 'BRO + CLO = BR+ CL';
k(:,i)=3.10e-12.*exp(420./T);
Gstr{i,1}='BRO'; Gstr{i,2}='CLO';
fBRO(i)=fBRO(i)-1; fCLO(i)=fCLO(i)-1; fBR(i)=fBR(i)+1; fCL(i)=fCL(i)+1;

%296
i=i+1;
Rnames{i} = 'BRO + IO = BR+ I';
k(:,i)=5.50e-12.*exp(760./T);
Gstr{i,1}='BRO'; Gstr{i,2}='IO';
fBRO(i)=fBRO(i)-1; fIO(i)=fIO(i)-1; fBR(i)=fBR(i)+1; fI(i)=fI(i)+1;

%297
i=i+1;
Rnames{i} = 'BRO + HO2 = HOBR';
k(:,i)=4.50e-12.*exp(460./T);
Gstr{i,1}='BRO'; Gstr{i,2}='HO2';
fBRO(i)=fBRO(i)-1; fHO2(i)=fHO2(i)-1; fHOBR(i)=fHOBR(i)+1;

%298
i=i+1;
Rnames{i} = 'BRO + NO = BR+ NO2';
k(:,i)=8.80e-12.*exp(260./T);
Gstr{i,1}='BRO'; Gstr{i,2}='NO';
fBRO(i)=fBRO(i)-1; fNO(i)=fNO(i)-1; fBR(i)=fBR(i)+1; fNO2(i)=fNO2(i)+1;

%299
i=i+1;
Rnames{i} = 'BRO + NO2 = BRN3';
k(:,i)=K_BRO_NO2;
Gstr{i,1}='BRO'; Gstr{i,2}='NO2';
fBRO(i)=fBRO(i)-1; fNO2(i)=fNO2(i)-1; fBRN3(i)=fBRN3(i)+1;

%300
i=i+1;
Rnames{i} = 'BRN2 = BR+ NO2';
k(:,i)=JBRNO2;
Gstr{i,1}='BRN2';
fBRN2(i)=fBRN2(i)-1; fBR(i)=fBR(i)+1; fNO2(i)=fNO2(i)+1;

%301
i=i+1;
Rnames{i} = 'BRN3 = 0.85BR+ 0.15BRO+ 0.15NO2+ 0.85NO3';
k(:,i)=JBRNO3;
Gstr{i,1}='BRN3';
fBRN3(i)=fBRN3(i)-1; fBR(i)=fBR(i)+0.85; fBRO(i)=fBRO(i)+0.15; fNO2(i)=fNO2(i)+0.15; fNO3(i)=fNO3(i)+0.85;

%302
i=i+1;
Rnames{i} = 'BRN3 = HNO3+ HOBR';
k(:,i)=2.5e-22.*H2O;
Gstr{i,1}='BRN3';
fBRN3(i)=fBRN3(i)-1; fHNO3(i)=fHNO3(i)+1; fHOBR(i)=fHOBR(i)+1;

%303
i=i+1;
Rnames{i} = 'BR + FORM = HBR+ CO+ HO2';
k(:,i)=1.70e-11.*exp(-800./T);
Gstr{i,1}='BR'; Gstr{i,2}='FORM';
fBR(i)=fBR(i)-1; fFORM(i)=fFORM(i)-1; fHBR(i)=fHBR(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%304
i=i+1;
Rnames{i} = 'ALD2 + BR = C2O3+ HBR';
k(:,i)=1.80e-11.*exp(-460./T);
Gstr{i,1}='ALD2'; Gstr{i,2}='BR';
fALD2(i)=fALD2(i)-1; fBR(i)=fBR(i)-1; fC2O3(i)=fC2O3(i)+1; fHBR(i)=fHBR(i)+1;

%305
i=i+1;
Rnames{i} = 'ALDX + BR = CXO3+ HBR';
k(:,i)=5.75e-11.*exp(-575./T);
Gstr{i,1}='ALDX'; Gstr{i,2}='BR';
fALDX(i)=fALDX(i)-1; fBR(i)=fBR(i)-1; fCXO3(i)=fCXO3(i)+1; fHBR(i)=fHBR(i)+1;

%306
i=i+1;
Rnames{i} = 'BR + ETH = FMBR+ FORM+ XO2H+ RO2';
k(:,i)=6.35e-15.*exp(-440./T);
Gstr{i,1}='BR'; Gstr{i,2}='ETH';
fBR(i)=fBR(i)-1; fETH(i)=fETH(i)-1; fFMBR(i)=fFMBR(i)+1; fFORM(i)=fFORM(i)+1; fXO2H(i)=fXO2H(i)+1; fRO2(i)=fRO2(i)+1;

%307
i=i+1;
Rnames{i} = 'BR + OLE = ALD2+ FMBR- PAR+ 0.92XO2H+ 0.08XO2N+ RO2';
k(:,i)=3.6e-12;
Gstr{i,1}='BR'; Gstr{i,2}='OLE';
fBR(i)=fBR(i)-1; fOLE(i)=fOLE(i)-1; fALD2(i)=fALD2(i)+1; fFMBR(i)=fFMBR(i)+1; fPAR(i)=fPAR(i)-1.0; fXO2H(i)=fXO2H(i)+0.92; fXO2N(i)=fXO2N(i)+0.08; fRO2(i)=fRO2(i)+1;

%308
i=i+1;
Rnames{i} = 'BR + IOLE = 0.5ACET+ 0.5ALDX+ FMBR+ 0.9XO2H+ 0.1XO2N+ RO2';
k(:,i)=9.3e-12;
Gstr{i,1}='BR'; Gstr{i,2}='IOLE';
fBR(i)=fBR(i)-1; fIOLE(i)=fIOLE(i)-1; fACET(i)=fACET(i)+0.5; fALDX(i)=fALDX(i)+0.5; fFMBR(i)=fFMBR(i)+1; fXO2H(i)=fXO2H(i)+0.9; fXO2N(i)=fXO2N(i)+0.1; fRO2(i)=fRO2(i)+1;

%309
i=i+1;
Rnames{i} = 'BR + ISOP = FMBR+ ISPD+ 0.88XO2H+ 0.12XO2N+ RO2';
k(:,i)=7.4e-11;
Gstr{i,1}='BR'; Gstr{i,2}='ISOP';
fBR(i)=fBR(i)-1; fISOP(i)=fISOP(i)-1; fFMBR(i)=fFMBR(i)+1; fISPD(i)=fISPD(i)+1; fXO2H(i)=fXO2H(i)+0.88; fXO2N(i)=fXO2N(i)+0.12; fRO2(i)=fRO2(i)+1;

%310
i=i+1;
Rnames{i} = 'BR + TERP = FMBR+ 0.7XO2H+ 0.3XO2N+ RO2';
k(:,i)=2.5e-11;
Gstr{i,1}='BR'; Gstr{i,2}='TERP';
fBR(i)=fBR(i)-1; fTERP(i)=fTERP(i)-1; fFMBR(i)=fFMBR(i)+1; fXO2H(i)=fXO2H(i)+0.7; fXO2N(i)=fXO2N(i)+0.3; fRO2(i)=fRO2(i)+1;

%311
i=i+1;
Rnames{i} = 'FMBR = HBR+ CO';
k(:,i)=2.78e-4;
Gstr{i,1}='FMBR';
fFMBR(i)=fFMBR(i)-1; fHBR(i)=fHBR(i)+1; fCO(i)=fCO(i)+1;

%312
i=i+1;
Rnames{i} = 'CH3I = I+ MEO2';
k(:,i)=JCH3I;
Gstr{i,1}='CH3I';
fCH3I(i)=fCH3I(i)-1; fI(i)=fI(i)+1; fMEO2(i)=fMEO2(i)+1;

%313
i=i+1;
Rnames{i} = 'MI2 = FORM+ I+ I';
k(:,i)=JMI2;
Gstr{i,1}='MI2';
fMI2(i)=fMI2(i)-1; fFORM(i)=fFORM(i)+1; fI(i)=fI(i)+1; fI(i)=fI(i)+1;

%314
i=i+1;
Rnames{i} = 'MIB = BR+ FORM+ I';
k(:,i)=JMIB;
Gstr{i,1}='MIB';
fMIB(i)=fMIB(i)-1; fBR(i)=fBR(i)+1; fFORM(i)=fFORM(i)+1; fI(i)=fI(i)+1;

%315
i=i+1;
Rnames{i} = 'MIC = CL+ FORM+ I';
k(:,i)=JMIC;
Gstr{i,1}='MIC';
fMIC(i)=fMIC(i)-1; fCL(i)=fCL(i)+1; fFORM(i)=fFORM(i)+1; fI(i)=fI(i)+1;

%316
i=i+1;
Rnames{i} = 'MB3 = BR+ BR+ BR+ CO+ HO2';
k(:,i)=JCHBR3;
Gstr{i,1}='MB3';
fMB3(i)=fMB3(i)-1; fBR(i)=fBR(i)+1; fBR(i)=fBR(i)+1; fBR(i)=fBR(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%317
i=i+1;
Rnames{i} = 'MB3 + OH = BR+ BR+ BR+ CO';
k(:,i)=9.00e-13.*exp(-360./T);
Gstr{i,1}='MB3'; Gstr{i,2}='OH';
fMB3(i)=fMB3(i)-1; fOH(i)=fOH(i)-1; fBR(i)=fBR(i)+1; fBR(i)=fBR(i)+1; fBR(i)=fBR(i)+1; fCO(i)=fCO(i)+1;

%318
i=i+1;
Rnames{i} = 'MB2 + OH = BR+ BR+ CO+ HO2';
k(:,i)=2.00e-12.*exp(-840./T);
Gstr{i,1}='MB2'; Gstr{i,2}='OH';
fMB2(i)=fMB2(i)-1; fOH(i)=fOH(i)-1; fBR(i)=fBR(i)+1; fBR(i)=fBR(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1;

%319
i=i+1;
Rnames{i} = 'MBC + OH = BR+ MEO2';
k(:,i)=2.10e-12.*exp(-880./T);
Gstr{i,1}='MBC'; Gstr{i,2}='OH';
fMBC(i)=fMBC(i)-1; fOH(i)=fOH(i)-1; fBR(i)=fBR(i)+1; fMEO2(i)=fMEO2(i)+1;

%320
i=i+1;
Rnames{i} = 'MBC2 + OH = BR+ MEO2';
k(:,i)=9.40e-13.*exp(-510./T);
Gstr{i,1}='MBC2'; Gstr{i,2}='OH';
fMBC2(i)=fMBC2(i)-1; fOH(i)=fOH(i)-1; fBR(i)=fBR(i)+1; fMEO2(i)=fMEO2(i)+1;

%321
i=i+1;
Rnames{i} = 'MB2C + OH = BR+ MEO2';
k(:,i)=9.00e-13.*exp(-420./T);
Gstr{i,1}='MB2C'; Gstr{i,2}='OH';
fMB2C(i)=fMB2C(i)-1; fOH(i)=fOH(i)-1; fBR(i)=fBR(i)+1; fMEO2(i)=fMEO2(i)+1;

%322
i=i+1;
Rnames{i} = 'I + HO2 = HI';
k(:,i)=1.50e-11.*exp(-1090./T);
Gstr{i,1}='I'; Gstr{i,2}='HO2';
fI(i)=fI(i)-1; fHO2(i)=fHO2(i)-1; fHI(i)=fHI(i)+1;

%323
i=i+1;
Rnames{i} = 'HI + OH = I';
k(:,i)=3e-11;
Gstr{i,1}='HI'; Gstr{i,2}='OH';
fHI(i)=fHI(i)-1; fOH(i)=fOH(i)-1; fI(i)=fI(i)+1;

%324
i=i+1;
Rnames{i} = 'I + NO2 = INO2';
k(:,i)=K_I_NO2;
Gstr{i,1}='I'; Gstr{i,2}='NO2';
fI(i)=fI(i)-1; fNO2(i)=fNO2(i)-1; fINO2(i)=fINO2(i)+1;

%325
i=i+1;
Rnames{i} = 'INO2 = I+ NO2';
k(:,i)=JINO2;
Gstr{i,1}='INO2';
fINO2(i)=fINO2(i)-1; fI(i)=fI(i)+1; fNO2(i)=fNO2(i)+1;

%326
i=i+1;
Rnames{i} = 'INO2 + INO2 = I2+ NO2+ NO2';
k(:,i)=4.70e-12.*exp(-1670./T);
Gstr{i,1}='INO2'; Gstr{i,2}='INO2';
fINO2(i)=fINO2(i)-1; fINO2(i)=fINO2(i)-1; fI2(i)=fI2(i)+1; fNO2(i)=fNO2(i)+1; fNO2(i)=fNO2(i)+1;

%327
i=i+1;
Rnames{i} = 'BR + BRN2 = BR2+ NO2';
k(:,i)=5e-11;
Gstr{i,1}='BR'; Gstr{i,2}='BRN2';
fBR(i)=fBR(i)-1; fBRN2(i)=fBRN2(i)-1; fBR2(i)=fBR2(i)+1; fNO2(i)=fNO2(i)+1;

%328
i=i+1;
Rnames{i} = 'GLY = CGLY';
k(:,i)=1.0e-6;
Gstr{i,1}='GLY';
fGLY(i)=fGLY(i)-1; fCGLY(i)=fCGLY(i)+1;

%329
i=i+1;
Rnames{i} = 'MGLY = CGLY';
k(:,i)=1.0e-6;
Gstr{i,1}='MGLY';
fMGLY(i)=fMGLY(i)-1; fCGLY(i)=fCGLY(i)+1;


