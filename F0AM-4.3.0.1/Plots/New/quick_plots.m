clear all; close all; 

pth="/Users/vanessasun/Documents/phd/utah/research/USOS_shared/F0AM-4.3.0.1/Runs/USOS_8_6_2024/Run23/USOS_8_6_2024_run.mat"; 
load(pth)

% Set start/end date in MST time... 
foam_start=datetime('08/06/2024','InputFormat','MM/dd/yyyy');
foam_end=datetime('08/07/2024','InputFormat','MM/dd/yyyy');
yr=year(foam_start);mon=month(foam_start);dy=day(foam_start); 
t_start=datetime(yr,mon,dy,0,0,0);
t_end=datetime(yr,mon,dy,23,30,0);
[USOS, sun]= get_subset_USOS(t_start, t_end);


% PlotConc('PAN', S); hold on 
% plot(S.Time, USOS.PAN_CIMS, 'k')
% legend('Model','Data')
% purtyPlot
% 
% PlotConc('PANX', S); hold on 
% plot(S.Time, USOS.PPN_CIMS, 'k')
% legend('Model','Data')
% purtyPlot
% 
% 
% 
% PlotConc('O3', S); hold on 
% plot(S.Time, USOS.O3_ppbv, 'k')
% legend('Model','Data')
% purtyPlot
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

% Total loss rate of nOX: 
NOx={'NO','NO2','NO3','HONO','CRON'};
NOx_res={'CLN2','CLN3','INO2','INO3','BRN2','BRN3','N2O5', 'PAN', 'PANX','HNO3', 'INTR','NTR1','NTR2','OPAN','PNA'}; 
[~,iNOx] = ismember(NOx,S.Cnames);
[~,i_res] = ismember(NOx_res,S.Cnames);
iLostNOx = sum(S.Chem.f(:,iNOx)==-1,2) & ... %use stoichiometric coefficients to ID reactions
           sum(S.Chem.f(:,i_res)==1,2);
LNOx_res = sum(S.Chem.Rates(:,iLostNOx),2)*3600; %total NOx lost to make res species

rxns=S.Chem.Rnames{iLostNOx};
NOxrates = PlotRates(['Nox_res',NOx_res],S,10,'unit','ppb_h','sumEq',1);
LNOxrate = sum(sum(NOxrates.Prod,2) + sum(NOxrates.Loss,2),2);

figure
plot(S.Time,LNOxrate,'k-',S.Time,LNOx_res,'r--')
xlabel('Hour of Day')
ylabel('NOx Loss Rate (ppb h^-^1)')
legend('P(NOx res)','Loses NOx to Res')
purtyPlot


NO2rates = PlotRates('NO2',S,15,'unit','ppb_h','sumEq',1, 'plotme',0);

area(S.Time, NO2rates.Prod(:,3:width(NO2rates.Prod)));
legend(NO2rates.Pnames{3:length(NO2rates.Pnames)})

% Get the MODELED net O3 production rate: 
O3rates = PlotRates('O3',S,5,'unit','ppb_h','sumEq',1);
O3netRate = sum(sum(O3rates.Prod,2) + sum(O3rates.Loss,2),2);

O3netRate(O3netRate<0)=0; 
OPE=O3netRate./LNOxrate;
figure
plot(S.Time, OPE)



% The net ozone rate calculated above should be roughly equal to the rate of NO + XO2 = NO2.
% Let's find reactions that destroy NO and a peroxy radical and produce NO2.
XO2names = {'BZO2'; 'C2O3';'CXO3';'EPX2'; 'ISO2';'MEO2';'OPO3';'TO2';'XLO2';'HO2';};
[~,iNO]  = ismember('NO',S.Cnames); %index location of species
[~,iNO2] = ismember('NO2',S.Cnames);
[~,iXO2] = ismember(XO2names,S.Cnames);
iNOtoNO2 = S.Chem.f(:,iNO)==-1 & ... %use stoichiometric coefficients to ID reactions
           S.Chem.f(:,iNO2)==1 & ...
           sum(S.Chem.f(:,iXO2)==-1,2);
PNO2fromXO2 = sum(S.Chem.Rates(:,iNOtoNO2),2)*3600; %total NO2 production from XO2+NO, ppb/h
figure
plot(S.Time,O3netRate,'k-',S.Time,PNO2fromXO2,'r--')
xlabel('Hour of Day')
ylabel('Ozone Production (ppb h^-^1)')
legend('O_3 Net','XO_2 + NO')
purtyPlot
% They don't quite add up. Why? The assumption that P(O3)net = sum(k[NO][XO2]) does not account for:
% 1) NO2 production from PANs
% 2) O3 losses that do not make NO2, like reaction with VOC or HO2, or dilution


