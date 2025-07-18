%% Only to get the observational values

% Set start/end date in MST time... 
foam_start=datetime('08/06/2024','InputFormat','MM/dd/yyyy');
foam_end=datetime('08/07/2024','InputFormat','MM/dd/yyyy');

% Extract the year, month day: 
yr=year(foam_start);mon=month(foam_start);dy=day(foam_start); 
RUNNAME= strcat('USOS','_',num2str(mon), '_', num2str(dy),'_', num2str(yr)) ;

% Use to set t_start and t_end to pass to get_subset to select only data
% during this time period : 
t_start=datetime(yr,mon,dy,0,0,0); 
t_end=datetime(yr,mon,dy,23,30,0);

% Pass to our function to load USOS data and select the appropriate time
% period. Also include sun position calculation. Outputs are: 
%   USOS: struct with USOS data from only this period 
%   utc_time: MATLAB datetime object added to USOS.utc_time
%   Time_MST: MATLAB datetime object added to USOS.Time_MST
%   sun: struct with zenith and azimuth estimates for sun position
[USOS, sun]= get_subset_USOS(t_start, t_end);
%% 

%%%%% Actual Plotting
Svarload = strcat(RUNNAME,'_rep_save.mat')
Svar_path = fullfile(Svarload)
load(Svar_path)

FieldList = fieldnames(Svar_path);
for iField = 1:numel(FieldList)
    Field    = FieldList{iField};
        out.(Field) = USOS.(Field);
end
%% 


% First, let's separate the three days using SplitRun.
% The first day is effectively "spin-up" for secondary and intermediate species.
SplitRun(Svarload.S,'custom',Svarload.S.repIndex)
%% 

% Now let's see how well we simulated NO and NO2, since only total NOx was "fixed".
S3.Conc.NOx = S3.Conc.NO+S3.Conc.NO2;
PlotConcGroup({'NO','NO2','NOx'},S3,3,'sortem',0,'ptype','line')
hold on
plot(USOS.timehr_output,USOS.NO_LIF,'b--')
hold off

hold on
plot(USOS.timehr_output,USOS.NO2_LIF,'r--')
hold off
% hold on
% plot(timehr_output,USOS.NO_LIF+USOS.NO2_LIF,'g--')
% hold off

text(0.55,0.7,'solid: model')
text(0.55,0.5,'dash: observed')
legend('NO','NO2','NOx')

figure
plot(USOS.timehr_output,USOS.NO_LIF./USOS.NO2_LIF,'k-')
hold on
plot(S3.Time,S3.Conc.NO./S3.Conc.NO2,'k--')
xlabel('Model Time')
ylabel('NO/NO2')
legend('Obs','Model')
% 
% Now, let's see how ozone did over the three days.
PlotConc('O3',{S1,S2,S3})
hold on
plot(USOS.timehr_output,USOS.O3_ppbv,'k-')
legend('Day 1','Day 2','Day 3','Obs')

%% PLOTTING AND ANALYSIS

% First, let's separate the three days using SplitRun.
% The first day is effectively "spin-up" for secondary and intermediate species.
%SplitRun(S,'custom',repIndex)
% 
% % % Now let's see how well we simulated NO and NO2, since only total NOx was "fixed".
% S3.Conc.NOx = S3.Conc.NO+S3.Conc.NO2;
% PlotConcGroup({'NO','NO2','NOx'},S3,3,'sortem',0,'ptype','line')
% hold on
% plot(USOS.time_local,USOS.NO,'b--')
% plot(USOS.time_local,USOS.NO2,'r--')
% plot(USOS.time_local,USOS.NO+USOS.NO2,'y--')
% text(0.55,0.7,'solid: model')
% text(0.55,0.6,'dash: observed')
% legend('NO','NO2','NOx')
% 
% figure
% plot(USOS.time_local,USOS.NO./USOS.NO2,'k-')
% hold on
% plot(S3.time_local,S3.Conc.NO./S3.Conc.NO2,'k--')
% xlabel('Model Time')
% ylabel('NO/NO2')
% legend('Obs','Model')
% % 
% % Now, let's see how ozone did over the three days.
% PlotConc('O3',{S1,S2,S3})
% hold on
% plot(USOS.time_local,USOS.O3,'k-')
% legend('Day 1','Day 2','Day 3','Obs')
% 
% % Next, let's look at ozone production on the last day.
% O3rates = PlotRates('O3',S3,5,'unit','ppb_h','sumEq',1);
% O3netRate = sum(O3rates.Prod + O3rates.Loss,2);
% 
% % The net ozone rate calculated above should be roughly equal to the rate of NO + XO2 = NO2.
% % Let's find reactions that destroy NO and a peroxy radical and produce NO2.
% XO2names = [S3.Cnames(S3.iRO2);'HO2'];
% [~,iNO]  = ismember('NO',S3.Cnames); %index location of species
% [~,iNO2] = ismember('NO2',S3.Cnames);
% [~,iXO2] = ismember(XO2names,S3.Cnames);
% iNOtoNO2 = S3.Chem.f(:,iNO)==-1 & ... %use stoichiometric coefficients to ID reactions
%            S3.Chem.f(:,iNO2)==1 & ...
%            sum(S3.Chem.f(:,iXO2)==-1,2);
% PNO2fromXO2 = sum(S3.Chem.Rates(:,iNOtoNO2),2)*3600; %total NO2 production from XO2+NO, ppb/h
% 
% figure
% plot(S3.time_local,O3netRate,'k-',S3.time_local,PNO2fromXO2,'r--')
% xlabel('Hour of Day')
% ylabel('Ozone Production (ppb h^-^1)')
% legend('O_3 Net','XO_2 + NO')
% purtyPlot
% % They don't quite add up. Why? The assumption that P(O3)net = sum(k[NO][XO2]) does not account for:
% % 1) NO2 production from PANs
% % 2) O3 losses that do not make NO2, like reaction with VOC or HO2, or dilution
% 
% % The O3 rates plot does not give us much info on what RO2 species are driving O3 production.
% % Let's look at NO2 production instead, averaged over the day.
% iday = S3.time_local>=8 & S3.time_local<=18; %peak production hours
% PlotRatesAvg('NO2',S3,10,'pts2avg',iday,'unit','ppb_h','sumEq',1) %sumEq=1 to take net rates for PAN equilibria
% 
% %Or, we could look at OH reactivity to see which VOC are contributing most.
% % Let's group measured species together to make the plot easier to interpret.
% %Inorg = {'Inorganic';'CO';'H2';'O3';'HO2';'H2O2';'NO2'};
% %Terp  = {'MTerp','APINENE','BPINENE','LIMONENE'};
% %Alk   = {'Alk','C2H4','C2H6','C3H8','IC4H10','IC5H12','NC5H12','NC6H14','NC10H22'};
% Arom  = {'Arom','BENZ','TOL'};
% %oVOC  = {'oVOC','CH3CHO','C2H5CHO','C3H7CHO','HOCH2CHO','GLYOX','CH3OH','C2H5OH','ACETOL','BIACET'};
% % MVKMACR = {'MVKMACR','MVK','MACR'};
% % Reactants = {Inorg,Alk,Arom,oVOC,'HCHO',MVKMACR,Terp,'C5H8'};
% % PlotReactivity('OH',S3,Reactants,'ptype','bar');
% % hold on
% % plot(SOAS.Time,SOAS.kOH,'k*-')
% % text(0.05,0.9,'black line: observed')


