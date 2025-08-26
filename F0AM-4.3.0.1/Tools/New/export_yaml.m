clearvars

%% Load the Runs
% Change this date when you plot a new date
foam_start=datetime('08/06/2024','InputFormat','MM/dd/yyyy');
foam_end=datetime('08/07/2024','InputFormat','MM/dd/yyyy');

% Extract the year, month day: 
yr=year(foam_start);mon=month(foam_start);dy=day(foam_start); 

% Use to set t_start and t_end to pass to get_subset to select only data
% during this time period : 
t_start=datetime(yr,mon,dy,0,0,0); 
t_end=datetime(yr,mon,dy,23,30,0);

savedir = '/Users/vanessasun/Documents/phd/utah/research/USOS_shared/F0AM-4.3.0.1/Runs/';
save_runname = strcat('USOS','_',num2str(mon), '_', num2str(dy),'_', num2str(yr));
full_savepath = strcat(savedir,save_runname,'/Run26/',save_runname);
Svarload = strcat(full_savepath,'_run.mat')
Svar_path = fullfile(Svarload)
load(Svar_path)

%% Load the observational values
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

%%



%%
%WriteYaml will not accept any datetimes, so we convert the datetimes in
%USOS to cells
convert_datetime_usos_time_mst = cellstr(USOS.time_MST);
convert_datetime_usos_time_utc = cellstr(USOS.time_UTC);
USOS.time_MST = convert_datetime_usos_time_mst;
USOS.time_UTC = convert_datetime_usos_time_utc;

obs_savedir_mac = '/Users/vanessasun/Documents/phd/utah/research/USOS_shared/F0AM-4.3.0.1/Campaign_Data/obs_yaml/parked/';
obs_date = strcat('USOS','_',num2str(mon), '_', num2str(dy),'_', num2str(yr));

full_obs_savepath = strcat(obs_savedir_mac,obs_date,'.yaml');
WriteYaml(full_obs_savepath,USOS)


yaml_savepath_run_mac = strcat(full_savepath,'_run.yaml')
WriteYaml(yaml_savepath_run_mac,S)


%%
 % % S.Obs=struct; 
 % % for f=1:length(fields)
 % %     field_i=fields{f}; 
 % %     S.Obs.(field)=USOS.(field); 
 % % end
 % % Overwrite with the obs appended to it too!  
 % %save(savepath, savename)
 % 
 %  % fields=fieldnames(USOS);
 % % S.Obs=struct; 
 % % for f=1:length(fields)
 % %     field_i=fields{f}; 
 % %     S.Obs.(field)=USOS.(field); 
 % % end
 % 
 %         % Loop over all vars in the .mat file and pull out a struc only containing data 
 %        % within our preferred time window... 
 %        FieldList = fieldnames(USOS);
 %        for iField = 1:numel(FieldList)
 %           Field    = FieldList{iField};
 %           if length(USOS.(Field))==length(USOS.time_MST)
 %              out.(Field) = USOS.(Field)(cond);
 %           else 
 %              out.(Field) = USOS.(Field);
 %           end
 %        end
 % 
 %        out=orderfields(out); % Sort the fieldnames alphabetically.