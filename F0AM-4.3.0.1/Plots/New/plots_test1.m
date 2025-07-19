%% Get the observational values
% Change this date when you plot a new date
foam_start=datetime('08/06/2024','InputFormat','MM/dd/yyyy');
foam_end=datetime('08/07/2024','InputFormat','MM/dd/yyyy');

% Extract the year, month day: 
yr=year(foam_start);mon=month(foam_start);dy=day(foam_start); 

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

%% Get Model Run values
savedir = '/Users/vanessasun/Documents/phd/utah/research/USOS_shared/F0AM-4.3.0.1/Runs/';
save_runname = strcat('USOS','_',num2str(mon), '_', num2str(dy),'_', num2str(yr));
full_savepath = strcat(savedir,save_runname,'_rep_save.mat');

Svar_path = fullfile(full_savepath);
svar = load(Svar_path);

%
%Separate the three days of Model Run values using SplitRun.
%The first day is effectively "spin-up" for secondary and intermediate species.
SplitRun(svar.S,'custom',svar.S.repInd)

%% 

%Test to see if we reached steady-state for a long-lived species
PlotConc('GLYD',{S1,S2,S3, S4})
legend('Day 1','Day 2','Day 3','Day 4')
% 
% % Now let's see how well we simulated NO and NO2, since only total NOx was "fixed".
% S3.Conc.NOx = S3.Conc.NO+S3.Conc.NO2;
% PlotConcGroup({'NO','NO2','NOx'},S3,3,'sortem',0,'ptype','line')
% hold on
% plot(USOS.timehr_output,USOS.NO_LIF,'b--')
% hold off
% 
% hold on
% plot(USOS.timehr_output,USOS.NO2_LIF,'r--')
% hold off
% % hold on
% % plot(timehr_output,USOS.NO_LIF+USOS.NO2_LIF,'g--')
% % hold off
% 
% text(0.55,0.7,'solid: model')
% text(0.55,0.5,'dash: observed')
% legend('NO','NO2','NOx')
% 
% figure
% plot(USOS.timehr_output,USOS.NO_LIF./USOS.NO2_LIF,'k-')
% hold on
% plot(S3.Time,S3.Conc.NO./S3.Conc.NO2,'k--')
% xlabel('Model Time')
% ylabel('NO/NO2')
% legend('Obs','Model')
% 
% Now, let's see how ozone did over the three days.
% PlotConc('O3',{S1,S2,S3})
% hold on
% plot(USOS.timehr_output,USOS.O3_ppbv,'k-')
% legend('Day 1','Day 2','Day 3','Obs')