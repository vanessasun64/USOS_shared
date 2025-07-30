function [out, sun]= get_subset_USOS(t_start, t_end)
% Function to get a subset of the data structure based on time, date,
    
    % Load in MATLAB version of the filled merge of campaign's data! 
    dirpath2USOS='/Users/vanessasun/Documents/phd/utah/research/F0AM-4.3.0.1/Campaign_Data/matlab_merge/parked/corrected/';
    USOSfilename = '20240804_20240808_30min_CSL_mobile_lab_parked_with_interp_pan_interp_struct_for_MATLAB.mat';
    
    dirpath2USOS_win= '\Users\u1545774\Documents\GitHub\USOS_shared\F0AM-4_3_0_1\Campaign_Data\matlab_merge\parked\corrected\';
    USOSfilename_win = '20240804_20240808_30min_CSL_mobile_lab_parked_with_interp_struct_for_MATLAB.mat';
    
    fullpath2USOS = fullfile('C:\',dirpath2USOS_win,USOSfilename_win);
    load(fullpath2USOS); % loads struct named USOS. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Get time to make sense in MATLAB and add a local time field to the struct 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % USOS time_UTC was a numpy timstamp object... It loads here in
    % nanoseconds since some epoch.... 
    
    % Convert USOS.Time_UTC from nano Seconds to Seconds then convert to a 
    % MATLAB datetime object assuming its seconds since that epoch
    % Add that BACK into the USOS struct so we can use it nicely in MATLAB
    % & not do this again: 
    USOS.time_UTC= datetime(USOS.time_UTC./1e9, 'ConvertFrom', 'posixtime');
    USOS.time_MST= datetime(USOS.time_local./1e9, 'ConvertFrom', 'posixtime');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Subselect USOS data ONLY when we want it: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create a 30min step array between the start and stop time, to know
    % how big out data SHOULD be, so we can check if we have a full day's worth or not! 
    ideal=(t_start:minutes(30):t_end);

    % Get Indices in data of the time window we've sub-selected data for. 
    cond=find(USOS.time_MST>= t_start & USOS.time_MST <= t_end); 

    if length(cond) == length(ideal)
        % Loop over all vars in the .mat file and pull out a struc only containing data 
        % within our preferred time window... 
        FieldList = fieldnames(USOS);
        for iField = 1:numel(FieldList)
           Field    = FieldList{iField};
           if length(USOS.(Field))==length(USOS.time_MST)
              out.(Field) = USOS.(Field)(cond);
           else 
              out.(Field) = USOS.(Field);
           end
        end

        out=orderfields(out); % Sort the fieldnames alphabetically.

       %% Set up F0AM Sun Meteorology Parameters

       %{
       P, T and RH were measured at the site and will be updated every step of the simulation.
       SZA was not measured, so we can use a function to calculate it.
       kdil is a physical loss constant for all species; 1 per day is a typical value.
       %}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Get the sun position solar zenith angle for this subselection of the data 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get local time as nums for y m d h m s ... 
        [y,mn,d]= ymd(out.time_MST); [h,m,s]= hms(out.time_MST); 

        time.year=y; 
        time.month=mn;
        time.day=d;
        time.hour=h;
        time.min=m;
        time.sec=s;
        time.UTC=-6; % offset hour from UTC. Local time = Greenwich time + time.UTC
        location.latitude   = out.Lat; % (in degrees, north of equator is positive)
        location.longitude  = out.Lon; %  (in degrees, positive for east of Greenwich)
        location.altitude   = out.Altitude_m; % Needs altitude above mean sea level (in meters) 
        sun = sun_position(time,location); %fields zenith and azimuth
        
        timestamp_collection = 0:0.5:23.5;
        out.timehr_output = transpose(timestamp_collection); %index assigning each timestamp to what hour of day it is
        %disp(out.timehr_output)
    else % If you didn't have data for the whole day, then return empties!
        out=[]; sun=[]; 
    end

    
end

    



