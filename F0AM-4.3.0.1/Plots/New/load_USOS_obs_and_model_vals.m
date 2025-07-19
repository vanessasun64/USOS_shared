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

%% Get Model Runs values
savedir = '/Users/vanessasun/Documents/phd/utah/research/USOS_shared/F0AM-4.3.0.1/Runs/';
save_runname = strcat('USOS','_',num2str(mon), '_', num2str(dy),'_', num2str(yr));
full_savepath = strcat(savedir,save_runname);

Svarload = strcat(full_savepath,'_rep_save.mat')
Svar_path = fullfile(Svarload)
load(Svar_path)

FieldList = fieldnames(Svar_path);
for iField = 1:numel(FieldList)
    Field    = FieldList{iField};
        out.(Field) = USOS.(Field);
end