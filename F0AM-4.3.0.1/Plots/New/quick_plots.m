clear all; close all; 

pth_win="C:\Users\u1545774\Documents\GitHub\USOS_shared\F0AM-4.3.0.1\Runs\USOS_8_6_2024\Run24\USOS_8_6_2024_run.mat"; 
load(pth_win)
%% 

% Set start/end date in MST time... 
foam_start=datetime('08/06/2024','InputFormat','MM/dd/yyyy');
foam_end=datetime('08/07/2024','InputFormat','MM/dd/yyyy');
yr=year(foam_start);mon=month(foam_start);dy=day(foam_start); 
t_start=datetime(yr,mon,dy,0,0,0);
t_end=datetime(yr,mon,dy,23,30,0);
[USOS, sun]= get_subset_USOS(t_start, t_end);
%% 
% PlotConc('CH4', S); hold on 
% plot(S.Time, USOS.CH4_Piccaro, 'k')
% legend('Model','Data')
% purtyPlot

% PlotConc('PAN', S); hold on 
% plot(S.Time, USOS.PAN_CIMS, 'k')
% legend('Model','Data')
% purtyPlot

% PlotConc('PANX', S); hold on 
% plot(S.Time, USOS.PPN_CIMS, 'k')
% legend('Model','Data')
% purtyPlot
% 

% figure
% plot(S.Time, S.Met.jcorr)
% ylabel('jcorr')
%% 

% 
% S.Conc.NOx = S.Conc.NO+S.Conc.NO2;
% PlotConcGroup({'NO','NO2','NOx'},S,3,'sortem',0,'ptype','line')
% hold on
% plot(S.Time,USOS.NO_LIF,'b--')
% plot(S.Time,USOS.NO2_LIF,'r--')
% plot(S.Time,USOS.NO_LIF+USOS.NO2_LIF,'y--')
% text(0.55,0.7,'solid: model')
% text(0.55,0.6,'dash: observed')
% legend('NO','NO2','NOx')



PlotConc('O3', S); hold on 
plot(S.Time, USOS.O3_ppbv, 'k')
legend('Model','Data')
purtyPlot
% 
% PlotConc('ISPD', S); hold on 
% plot(S.Time, USOS.MVK_MACR_PTR+USOS.C5H10O3_CIMS, 'k')
% legend('Model','Data')
% purtyPlot
% 
% PlotConc('FACD', S); hold on 
% plot(S.Time, USOS.HCOOH_CIMS, 'k')
% legend('Model','Data')
% purtyPlot
% 
% PlotRates('FACD', S,5 );


%% 

% 
% 
% NO2rates = PlotRates('NO2',S,15,'unit','ppb_h','sumEq',1, 'plotme',0);
% 
% area(S.Time, NO2rates.Prod(:,3:width(NO2rates.Prod)));
% legend(NO2rates.Pnames{3:length(NO2rates.Pnames)})
% 
% %Or, we could look at OH reactivity to see which VOC are contributing most.
% % Let's group measured species together to make the plot easier to interpret.
% Inorg = {'Inorganic';'CO';'H2';'O3';'HO2';'H2O2';'NO2'};
% Terp  = {'MTerp','APINENE','BPINENE','LIMONENE'};
% Alk   = {'Alk','C2H4','C2H6','C3H8','IC4H10','IC5H12','NC5H12','NC6H14','NC10H22'};
% Arom  = {'Arom','BENZENE','TOLUENE','EBENZ','TM124B','TM135B','MXYL','OXYL','PXYL','BENZAL'};
% oVOC  = {'oVOC','CH3CHO','C2H5CHO','C3H7CHO','HOCH2CHO','GLYOX','CH3OH','C2H5OH','ACETOL','BIACET'};
% MVKMACR = {'MVKMACR','MVK','MACR'};
% Reactants = {Inorg,Alk,Arom,oVOC,'HCHO',MVKMACR,Terp,'C5H8'};
% PlotReactivity('OH',S3,Reactants,'ptype','bar');
% hold on
% plot(SOAS.Time,SOAS.kOH,'k*-')
% text(0.05,0.9,'black line: observed')
% 

%% 

% % Get the MODELED net O3 production rate: 
O3rates = PlotRates('O3',S,5,'unit','ppb_h','sumEq',1);
O3netRate = sum(sum(O3rates.Prod,2) + sum(O3rates.Loss,2),2);

figure
plot(S.Time,O3rates.Prod,'k-',S.Time,O3rates.Loss,'r--')
xlabel('Hour of Day')
ylabel('Ozone (ppb h^-^1)')
legend('P(O3)','L(O3)')
purtyPlot
%% 


O3netRate(O3netRate<0)=0; 
OPE=O3netRate./LNOxrate;
figure
plot(S.Time, OPE)
ylabel('OPE')
% 
% 
% 
% % The net ozone rate calculated above should be roughly equal to the rate of NO + XO2 = NO2.
% % Let's find reactions that destroy NO and a peroxy radical and produce NO2.
% XO2names = {'BZO2'; 'C2O3';'CXO3';'EPX2'; 'ISO2';'MEO2';'OPO3';'TO2';'XLO2';'HO2';};
% [~,iNO]  = ismember('NO',S.Cnames); %index location of species
% [~,iNO2] = ismember('NO2',S.Cnames);
% [~,iXO2] = ismember(XO2names,S.Cnames);
% iNOtoNO2 = S.Chem.f(:,iNO)==-1 & ... %use stoichiometric coefficients to ID reactions
%            S.Chem.f(:,iNO2)==1 & ...
%            sum(S.Chem.f(:,iXO2)==-1,2);
% PNO2fromXO2 = sum(S.Chem.Rates(:,iNOtoNO2),2)*3600; %total NO2 production from XO2+NO, ppb/h
% figure
% plot(S.Time,O3netRate,'k-',S.Time,PNO2fromXO2,'r--')
% xlabel('Hour of Day')
% ylabel('Ozone Production (ppb h^-^1)')
% legend('O_3 Net','XO_2 + NO')
% purtyPlot
% They don't quite add up. Why? The assumption that P(O3)net = sum(k[NO][XO2]) does not account for:
% 1) NO2 production from PANs
% 2) O3 losses that do not make NO2, like reaction with VOC or HO2, or dilution


