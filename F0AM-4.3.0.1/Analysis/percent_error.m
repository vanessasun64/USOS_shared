clear all; close all; 

pth_win="C:\Users\u1545774\Documents\GitHub\USOS_shared\F0AM-4.3.0.1\Runs\USOS_8_6_2024\Run23\USOS_8_6_2024_run.mat"; 
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

half_day_time = S.Time(24:48);
obs_o3_half_day = USOS.O3_ppbv(24:48);
model_o3_half_day = S.Conc.O3(24:48);

percent_error_num = abs(model_o3_half_day - obs_o3_half_day)
percent_error_den = obs_o3_half_day
percent_error_calc = (percent_error_num ./ percent_error_den) .* 100
mean(percent_error_calc)