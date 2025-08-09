clear all; close all; 

%Load run
pth_win="C:\Users\u1545774\Documents\GitHub\USOS_shared\F0AM-4.3.0.1\Runs\USOS_8_6_2024\Run26\USOS_8_6_2024_run.mat"; 
load(pth_win)

%Load USOS observational data
% Set start/end date in MST time... 
foam_start=datetime('08/06/2024','InputFormat','MM/dd/yyyy');
foam_end=datetime('08/07/2024','InputFormat','MM/dd/yyyy');
yr=year(foam_start);mon=month(foam_start);dy=day(foam_start); 
t_start=datetime(yr,mon,dy,0,0,0);
t_end=datetime(yr,mon,dy,23,30,0);
[USOS, sun]= get_subset_USOS(t_start, t_end);
%% 
PlotConc('O3', S); hold on 
plot(S.Time, USOS.O3_ppbv, 'k')
legend('Model','Data')
purtyPlot

S.Conc.NOx = S.Conc.NO+S.Conc.NO2;
PlotConcGroup({'NO','NO2','NOx'},S,3,'sortem',0,'ptype','line'); hold on
legend('NO','NO2','NOx')


%%
%Get the Ozone Production and Loss Rates
figure
O3rates = PlotRates('O3',S,10,'unit','ppb_h','sumEq',1);
O3netRate = sum(sum(O3rates.Prod,2) + sum(O3rates.Loss,2),2); %Net Ozone

o3_prod_rxns = S.Chem.Rnames(O3rates.iRx_Prod);

%% 
figure
plot(S.Time,O3rates.Prod,'k-',S.Time,O3rates.Loss,'r--')
xlabel('Hour of Day')
ylabel('Ozone (ppb h^-^1)')
legend('P(O3)','L(O3)')
purtyPlot

plot(S.Time,O3netRate)

%%

% % Total loss rate of NOx: 
% NOx={'NO','NO2','NO3','HONO','CRON'}; % HONO: nitrous acid, CRON: nitro-cresols
% 
% %NOx reservoir species
% NOx_res={'N2O5', 'PAN','HNO3', 'PANX','CLN2','CLN3','INO2','INO3','BRN2','BRN3', 'INTR','NTR1','NTR2','OPAN','PNA'}; 
% [~,iNOx] = ismember(NOx,S.Cnames);
% [~,i_res] = ismember(NOx_res,S.Cnames);
% 
% %identify where you're losing NOx and also where you're making a NOx reservoir species
% iLostNOx = sum(S.Chem.f(:,iNOx)==-1,2) & ... %use stoichiometric coefficients to ID reactions; identifies where you're losing NOx
%            sum(S.Chem.f(:,i_res)==1,2); %identifies where you're making NOx reservoir species
% 
% LNOx_res = sum(S.Chem.Rates(:,iLostNOx),2)*3600; %net loss rate of NOx to make res species; conversion from ppb/s to ppb/h
% 
% loss_rxns=S.Chem.Rnames(iLostNOx); %identify which reactions occurred where you lost NOx and are making res species
% disp(loss_rxns)



%%
%Compare total ozone concentration to total reservoir species concentration
LNOx_res_conc = ExtractSpecies(NOx_res,S)
cell_convert = struct2cell(LNOx_res_conc);
cell_to_mat = cell2mat(transpose(cell_convert));
LNOx_res_conc_totalSum = sum(cell_to_mat,2);
disp(LNOx_res_conc_totalSum)

figure
scatter(LNOx_res_conc_totalSum,S.Conc.O3);
xlabel('NOz conc (ppb)');
ylabel('O3 conc (ppb)')

%%

NOxrates = PlotRates(['NOx',NOx],S,10,'unit','ppb_h','sumEq',1, 'plotme',1); %Plot of overall NOx production and loss
purtyPlot

NOx_res_rates = PlotRates(['NOx_{res}',NOx_res],S,10,'unit','ppb_h','sumEq',1, 'plotme',1); %Plot of overall NOx res species production and loss
purtyPlot

%%
NOx_net_rate = sum(sum(NOxrates.Prod,2) + sum(NOxrates.Loss,2),2); %Net NOx 
NOx_res_net_rate = sum(sum(NOx_res_rates.Prod,2) + sum(NOx_res_rates.Loss,2),2); %Net NOx reservoir 

figure
plot(S.Time, NOx_res_net_rate,'y', S.Time, LNOx_res, 'r')
legend('NOx_res_net_rate', 'LNOx_res')

figure
plot(S.Time, NOx_net_rate,'g',S.Time,NOx_res_net_rate, 'b')
ylabel('NOx Net Rate')
legend('NOx', 'NOx res')


%%
figure
plot(S.Time,NOx_res_net_rate,'k-',S.Time,LNOx_res,'r--')
xlabel('Hour of Day')
ylabel('NOx Loss Rate (ppb h^-^1)')
legend('P(NOx res)','Total NOx lost to make res')
purtyPlot

%%
%Plot of OPE
% figure
% O3netRate(O3netRate<0)=0; 
% OPE = O3netRate./NOx_res_net_rate;
% plot(S.Time,OPE, 'b')
% ylabel('OPE')
% purtyPlot
%%
figure
plot(S.Time, O3netRate, 'm', S.Time, LNOx_res, 'c');
legend('Net P(O3)','NOx Lost to make NOx res')
xlabel('Time (Hours)')
ylabel('Rate (ppb/hr)')
purtyPlot

%%
figure
O3_frac = O3netRate./LNOx_res
plot(S.Time, O3_frac,'b')
xlabel('Time (Hours)')
ylabel('OPE (ppb/hr)')
purtyPlot

%%
% figure
% rate_RO2_NO_rxn = (S.Chem.Rates(:,68) .* 3600) .* S.Conc.RO2 .* S.Conc.NO;
% frac_contr_from_RO2_NO = rate_RO2_NO_rxn ./LNOx_res
% plot(S.Time, frac_contr_from_RO2_NO, 'r', S.Time, O3_frac, 'b')
% legend('Contribution of HO2 + NO to OPE', 'Total OPE')
% plot(S.Time, frac_contr_from_RO2_NO, 'r')
% xlabel('Time (hours)')
% ylabel('Fractional contribution of RO2 + NO rxn')
% title('Contribution of RO2 + NO to OPE')

%%
figure
PO3_from_HO2_NO = (S.Chem.Rates(:,25) .* 3600) %ppb/s to ppb/h conversion
frac_contr_from_HO2_NO = PO3_from_HO2_NO./LNOx_res
% plot(S.Time, frac_contr_from_HO2_NO, 'r', S.Time,OPE, 'b')
% legend('Contribution of HO2 + NO to OPE', 'Total OPE')
plot(S.Time, frac_contr_from_HO2_NO, 'r')
legend('Contribution of HO2 + NO to OPE')
%%
contribution_from_NOx_res = NOx_res_net_rate./NOx_net_rate;
plot(S.Time, contribution_from_NOx_res,'m',S.Time, frac_contr_from_HO2_NO,'y',S.Time,OPE, 'b')
legend('Contribution of RO2 rxns', 'Contribution of HO2 + NO to OPE', 'Total OPE')

%%
%Numerator of OPE: Production of O3 from
%HO2 + NO -> NO2 + OH and RO2 + NO -> RO + NO2 rxns

RO2_list = {'C2O3','CXO3','MEO2','HCO3','ISO2','EPX2','BZO2','TO2','XLO2','OPO3'};
[~,iRO2] = ismember(RO2_list,S.Cnames);
[~,iNO]  = ismember('NO',S.Cnames); %index location of species
[~,iNO2] = ismember('NO2',S.Cnames)

%identify where you're destroying NO and RO2 to make NO2
iRO2_NO_conversion = S.Chem.f(:,iNO)==-1 & ... %identifies where you're losing NO
                     S.Chem.f(:,iNO2)==1 & ...%identifies where you're making NO2
                     sum(S.Chem.f(:,iRO2)==-1,2); %identifies where you're losing RO2; 

iRO2_NO_conversion_net = sum(S.Chem.Rates(:,iRO2_NO_conversion),2)*3600; %total NO2 production from RO2+NO in ppb/h
conversion_rxns=S.Chem.Rnames(iRO2_NO_conversion);

figure
plot(S.Time, iRO2_NO_conversion_net);
xlabel('Time (hours)')
ylabel('Net NO2 Production Rate (ppb/hr)')
title('Total NO2 produced from RO2+NO rxn')

%With HO2+NO rxn, we should have approximately the ozone production rate
figure
PO3_from_HO2_NO = (S.Chem.Rates(:,25) .* 3600); %ppb/s to ppb/h conversion
PO3_from_HO2_NO_and_RO2_NO = iRO2_NO_conversion_net + PO3_from_HO2_NO;
plot(S.Time, PO3_from_HO2_NO_and_RO2_NO);
xlabel('Time (hrs)');
ylabel('Rate (ppb/hr)');
title('Net NO to NO2 conversion from HO2+NO and RO2+NO');

figure
RO2_NO_approx = (S.Chem.Rates(:,68) .* 3600);
plot(S.Time, iRO2_NO_conversion_net,'r', S.Time,RO2_NO_approx,'k');
legend('Calculated conversion rate', 'RO2 + NO rxn approximation')

%%
%Denominator for OPE, loss of NOx to make NOx reservoirs, organic nitrates
% Total loss rate of NOx: 
NOx={'NO','NO2','NO3'} 

%NOx reservoir species
NOx_res={'N2O5', 'PAN','HNO3', 'PANX','CLN2','CLN3','INO2','INO3','BRN2','BRN3', 'INTR','NTR1','NTR2','OPAN','PNA','HONO','CRON'}; 
[~,iNOx] = ismember(NOx,S.Cnames);
[~,i_res] = ismember(NOx_res,S.Cnames);

%identify where you're losing NOx and also where you're making a NOx reservoir species
iLostNOx = sum(S.Chem.f(:,iNOx)==-1,2) & ... %use stoichiometric coefficients to ID reactions; identifies where you're losing NOx
           sum(S.Chem.f(:,i_res)==1,2); %identifies where you're making NOx reservoir species

LNOx_res = sum(S.Chem.Rates(:,iLostNOx),2)*3600; %net loss rate of NOx to make res species; conversion from ppb/s to ppb/h

loss_rxns=S.Chem.Rnames(iLostNOx); %identify which reactions occurred where you lost NOx and are making res species

%NEW OPE Estimate
estimated_OPE = PO3_from_HO2_NO_and_RO2_NO / LNOx_res;

figure
plot(S.Time,estimated_OPE);
xlabel('Time (hrs)');
ylabel('OPE (ppb/hr)');
title('Estimated OPE with only whole yields');


%% IF COEFFICIENTS > 1 ARE ALLOWED:
%Numerator of OPE: Production of O3 from
%HO2 + NO -> NO2 + OH and RO2 + NO -> RO + NO2 rxns

RO2_list = {'C2O3','CXO3','MEO2','HCO3','ISO2','EPX2','BZO2','TO2','XLO2','OPO3'};
[~,iRO2] = ismember(RO2_list,S.Cnames);
[~,iNO]  = ismember('NO',S.Cnames); %index location of species
[~,iNO2] = ismember('NO2',S.Cnames)

%identify where you're destroying NO and RO2 to make NO2
iRO2_NO_conversion_with_partial_yield = S.Chem.f(:,iNO)==-1 & ... %identifies where you're losing NO
                     S.Chem.f(:,iNO2)>0 & ...%identifies where you're making NO2
                     sum(S.Chem.f(:,iRO2)==-1,2); %identifies where you're losing RO2; 

iRO2_NO_conversion_net_with_partial_yield = sum(S.Chem.Rates(:,iRO2_NO_conversion_with_partial_yield),2)*3600; %total NO2 production from RO2+NO in ppb/h
conversion_rxns_with_partial_yield=S.Chem.Rnames(iRO2_NO_conversion_with_partial_yield);

figure
plot(S.Time, iRO2_NO_conversion_net_with_partial_yield);
xlabel('Time (hours)')
ylabel('Net NO2 Production Rate (ppb/hr)')
title('Total NO2 produced from RO2+NO rxn, including partial yield')

PO3_from_HO2_NO_and_RO2_NO_partial_yield = PO3_from_HO2_NO + iRO2_NO_conversion_net_with_partial_yield;

%%
%Denominator for OPE, loss of NOx to make NOx reservoirs and organic nitrates, RO2+NO -> RONO2

iLostNOx_with_partial_yield = sum(S.Chem.f(:,iNOx)==-1,2) & ... %identifies where you're losing NOx
           sum(S.Chem.f(:,i_res)>0,2); %identifies where you're making NOx reservoir species

lost_NOx_rates = S.Chem.Rates(:,iLostNOx_with_partial_yield)
loss_rxns_with_partial_yield= S.Chem.f(:,i_res)


%%

LNOx_res_with_partial_yield = sum(S.Chem.Rates(:,iLostNOx_with_partial_yield),2)*3600; %net loss rate of NOx to make res species; conversion from ppb/s to ppb/h

loss_rxns_with_partial_yield=S.Chem.Rnames(iLostNOx_with_partial_yield); %identify which reactions occurred where you lost NOx and are making res species

estimated_OPE_partial_yield = PO3_from_HO2_NO_and_RO2_NO_partial_yield / LNOx_res_with_partial_yield;

%%
%Compare numerators to RO2 + NO rxn in mech
figure
plot(S.Time, iRO2_NO_conversion_net,'r', S.Time,RO2_NO_approx,'k', S.Time,iRO2_NO_conversion_net_with_partial_yield,'m');
legend('Each RO2 + NO, with coeff = 1', 'RO2+NO rxn in mech', 'Each RO2 + NO, with partial yields')
title('Conversion of NO to NO2 from')


%%


%%
% RO2_NO_rxns=S.Chem.Rnames(iRO2_NO_conversion);
% [~,iRO2_NO_rxns] = ismember(RO2_NO_rxns,S.Chem.Rnames);
% %%
% RO2_conc_extract = ExtractSpecies(RO2_list,S)
% RO2_cell_convert = struct2cell(RO2_conc_extract);
% RO2_cell_to_mat = cell2mat(transpose(RO2_cell_convert));
% 
% RO2_NO_rates_total = S.Chem.Rates(:,iRO2_NO_rxns) .* 3600 
% PO3_from_RO2_NO = RO2_NO_rates_total .* RO2_cell_to_mat .* S.Conc.NO;
% sum_PO3_from_RO2_NO = sum(PO3_from_RO2_NO,2);
% 
% figure
% plot(S.Time, sum_PO3_from_RO2_NO)
% xlabel("Time (hrs)")
% ylabel('Summed kRO2+NO[RO2][NO] (ppb/hr)')
% title('Ozone Production from RO2+NO rxn')
% 
% %%
% figure
% frac_contr_from_RO2_NO_list = sum_PO3_from_RO2_NO ./ LNOx_res;
% plot(S.Time, frac_contr_from_RO2_NO_list, 'm')
% title('Contribution of RO2 + NO to OPE')
% 
% figure
% sum_HO2_and_RO2_rxns = frac_contr_from_HO2_NO + frac_contr_from_RO2_NO_list;
% plot(S.Time,sum_HO2_and_RO2_rxns,'c')
% title('Contribution to O3 from HO2+NO and RO2+NO rxns')

% PlotConc('PAN', S); hold on 
% plot(S.Time, USOS.PAN_CIMS, 'k')
% legend('Model','Data')
% purtyPlot

