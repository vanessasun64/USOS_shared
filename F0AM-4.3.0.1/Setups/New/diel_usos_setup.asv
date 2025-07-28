% ExampleSetup_DielCycle.m
% This example shows a model setup for simulation of an "average" diurnal cycle at a ground location.
% In particular, we will try to simulate ozone production.
% Read comments in each section for a guided tour.
%
% 20151126 GMW

clear

%% OBSERVATIONS & METEOROLOGY
%{
Constraints are taken from observations during the USOS field campaign.
The file loaded below contains observations (resolution: 30 minutes) for a subset
of multiple days for the campaign. For our F0AM run, we only want one day
to be analyzed at a time so we are importing our data as a variable and
limiting to the observations for the one day.

Note that constraints CANNOT contains NaNs or negative numbers (data in this file has already been filtered).
%}

% load 20240804_20240808_30min_CSL_mobile_lab_parked_with_interp_struct_for_MATLAB.mat %structure "USOS"
% Set start/end date in MST time... 
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


%Set where to store the total file
savedir = '/Users/vanessasun/Documents/phd/utah/research/USOS_shared/F0AM-4.3.0.1/Runs/';
runname_str = strcat('USOS','_',num2str(mon), '_', num2str(dy),'_', num2str(yr));
dir_path = strcat(savedir,runname_str,'/','Run9/');
mkdir(dir_path);
full_savepath = strcat(dir_path,runname_str);

%% Set Model Options

for spn=0:1
        % Set the date of the simulation. As start for 1st dat, and as start+NDays for the rest. 
       
    load_spinup = spn;    % Set=1 to load data in file "spin_up" as InitConc for F0AM. 
            
    % Now use that to configure F0AM to repeat runs and such! 
    if load_spinup ==1
        type='_run';
        spin_up= strcat(full_savepath,'_spinup.mat'); 
        nRep=0;
        hold=1; % Hold these species constant throughout model step.   Always set to 1. 
        evolve=0; % Don't hold these species constant throughout model step, allow to evolve on regular run  
    else
        nRep=4; %Set number of days to repeat on spin up.
        hold=1; % Hold these species constant throughout model step. Always set to 1. 
        evolve=1; % Also hold these species constant throughout model step on spin up! 
        type='_spinup'; 
    end

    %% 
    %Save USOS output as a CSV file for use in Python
    USOS_observed_table=struct2table(USOS);
    USOS_observed_save_path = strcat(full_savepath,'_observed_conc','.csv');
    writetable(USOS_observed_table, USOS_observed_save_path)
    
    %% 
    
    
    %Set other Met parameters besides sun position
    %{
    P, T and RH were measured at the site and will be updated every step of the simulation.
    SZA was not measured, so we can use a function to calculate it.
    kdil is a physical loss constant for all species; 1 per day is a typical value.
    %}
    
    Met = {...
    %   names       %values
        'P'          USOS.Pressure_mb; %Pressure, mbar
        'T'          USOS.Temp_K; %Temperature, K
        'RH'         USOS.RH_percent; %Relative Humidity, %
        'SZA'        sun.zenith; %solar zenith angle, degrees
        'kdil'       1/(32*60*60); %dilution constant, /s
        'jcorr'      0.3; %optimizes comparison b/w model and observed NO/NO2

        %example: 'J3' USOS.jBr2*0.3 (jcorr)
        
        'J4' USOS.jNO2_meas
        'Jn31' USOS.jBrCl
        'Jn24' USOS.jBr2
        USOS.jCCl4
        USOS.jCH2Oa
        USOS.jCH2Ob 
        USOS.jClNO2
        USOS.jClOa
        USOS.jClOb
        USOS.jCl2
        USOS.jHNO2
        USOS.jHNO3
        USOS.jI2
        USOS.jNO2 
        USOS.jNO3a
        USOS.jNO3b 
        USOS.jN2O5 
        USOS.jO3
        %ozone column info here
       get a time varying jcorr ()
        };
    
    %% CHEMICAL CONCENTRATIONS
    %{
    Concentrations are initialized using observations or fixed values.
    Species with HoldMe = 1 will be held constant throughout each step.
    Species with HoldMe = 0 are only initialized at the start of the run, because
     ModelOptions.LinkSteps=1 (see below). For this particular case, NO2 and O3 are
     unconstrained because we are investigating ozone production.
    When many species are used, it helps to organize alphabetically or by functional group.
    %}

    InitConc = {...
        % names           conc(ppb)                HoldMe
        
        %Inorganics
        'O3'                USOS.O3_ppbv            evolve;
        'CO'                USOS.CO_Piccaro         hold;
        'DMS'               USOS.DMS_PTR            hold;
    
        %NOy
        'NO'                USOS.NO_LIF             0;
        'NO2'               USOS.NO2_LIF            0;
        'HONO'              USOS.HONO_CIMS          hold;
        'HNO3'              USOS.HNO3_CIMS          hold;
        %'PAN'               USOS.PAN_CIMS           1;
        %'PANX'              USOS.PPN_CIMS           1;
        'NOx'               {'NO2','NO'}        []; %family conservation
    
        %Biogenics
        'ISOP'              USOS.Isoprene_PTR       hold;
        'TERP'              USOS.Monoterpenes_PTR   hold;
    
        %CxHy (Hydrocarbons)
        'CH4'               USOS.CH4_Piccaro        hold;
    
        %Aromatics
        'BENZ'              USOS.Benzene_PTR        hold;
        'TOL'               USOS.Toluene_PTR        hold;
    
        %Oxygenates
        'MEOH'              USOS.Methanol_PTR       hold;
        'ALD2'              USOS.Acetaldehyde_PTR   hold;
        'ETOH'              USOS.Ethanol_PTR        hold;
        'FORM'              USOS.HCHO_CRDS          hold; %Formaldehyde 
        'FACD'              USOS.HCOOH_CIMS         hold; %Formic Acid
        
        %Br and Cl
        'BR2'               USOS.Br2_CIMS           hold;
        'CL2'               USOS.Cl2_CIMS           hold;
        'CLN2'              USOS.ClNO2_CIMS         hold;
        'BRCL'              USOS.BrCl_CIMS          hold;
        'BRO'               USOS.BrO_CIMS           hold;
        };

    %% Decide if you want to load in spin - up concentrations or not. 
        if load_spinup ==1 % If this is a spin up, then load the spin up file, and 
            % Use contents to initialize the Met input instead. 
            load(spin_up); 
            % Pull the last point from that run, and spin it into InitConc & Met ready
            % for the "real" model run ... 
            [InitConc,Met_old] = Run2Init(S,length(S.Conc.NO2), InitConc);    
        end

    %% CHEMISTRY
    %{
    ChemFiles is a cell array of strings specifying functions and scripts for the chemical mechanism.
    THE FIRST CELL is always a function for generic K-values.
    THE SECOND CELL is always a function for J-values (photolysis frequencies).
    All other inputs are scripts for mechanisms and sub-mechanisms.
    Here we give example using MCMv3.3.1. Note that this mechanism was extracted from the MCM website for
    the specific set of initial species included above.
    %}
    ChemFiles = {...
       'CB6r5h_K(Met)'; %contains all k-rates
       'CB6r5h_J(Met,2)'; %photolysis file, flag of 2 means we're using the Hybrid Method
       'CB6r5h_AllRxns'; %Lists all reactions and species to add
       };
    
    %% DILUTION CONCENTRATIONS
    %{
    Background concentrations, along with the value of kdil in Met, determine the dilution rate for chemical species.
    Here we stick with the default value of 0 for all species, which effectively makes dilution a first-order loss.
    %}
    BkgdConc = {'DEFAULT'       0
    'PAN' 0.2   
    };
    
    %% OPTIONS
    %{
    "Verbose" can be set from 0-3; this just affects the level of detail printed to the command
      window regarding model progress.
    "EndPointsOnly" is set to 1 because we only want the last point of each step.
    "LinkSteps" is set to 1 so that non-constrained species are carried over between steps.
    "IntTime" is the integration time for each step, equal to the spacing of the data (60 minutes).
    "TimeStamp" is set to the hour-of-day for observations.
    "SavePath" give the filename only (in this example); the default save directory is the UWCMv3\Runs folder.
    "FixNOx" forces total NOx to be reset to constrained values at the beginning of every step.
    %}
    
    ModelOptions.Verbose        = 2;
    ModelOptions.EndPointsOnly  = 1;
    ModelOptions.LinkSteps      = 1;
    ModelOptions.IntTime        = 1800; %change to timestep to match frequency of data
    ModelOptions.TimeStamp      = USOS.timehr_output; 
    %ModelOptions.SavePath       = full_savepath;
    % ModelOptions.FixNOx         = 1; %if you use this, disable family conservation above.
    
    if spn == 1
       ModelOptions.SavePath = strcat(full_savepath,'_run.mat');
    else
       ModelOptions.SavePath = strcat(full_savepath,'_spinup.mat');
    end

%% INPUT REPLICATION AND INTERPOLATION
% For this particular scenario, it might be desirable to modify the inputs in a few ways.
% This sections demonstrates how to do so.

% INTERPOLATION
% Inputs currently have a time resolution of 60 minutes, but this is pretty coarse (the sun can move
% a lot in 60 minutes). The InputInterp function allows you to interpolate all inputs to a finer
% time resolution. NOTES:
%   - If your native data is fast (e.g., 1 Hz), it is generally better practice to bin-average that 
%       data to your desired resolution rather than average down to 60 minutes and then interpolate as done here.
%   - Make sure you adjust ModelOptions.IntTime too!
% To turn this on, set the "0" to "1" below.
% if 0
%     dt = 1800; %time spacing, seconds
% 
%     Time_interp = (0:dt:(86400-dt))'/3600; %interpolation timebase, fractional hours (to match SOAS.Time)
%     circularFlag = 1; % time wraps around at midnight
%     [Met,InitConc,BkgdConc] = ...
%         InputInterp(Met,InitConc,BkgdConc,USOS.time_local,Time_interp,circularFlag);
%     ModelOptions.TimeStamp = Time_interp;
%     ModelOptions.IntTime = dt;
% end

%% REPLICATION
% Sometimes you may want to run the same inputs for multiple times. Typically, this scenario would
% be ground-site observations over one or more days, and you need a "spin-up" for non-measured
% species. The InputReplicate function lets you do this. Note, this only makes sense to use if
% ModelOptions.LinkSteps = 1. This replaces the "ModelOptions.Repeat" functionality in model
% versions prior to F0AMv4.
% Here, we run the same contraints for 3 days.
% The output "repIndex" is used to separate the days with SplitRun later.

    if nRep > 1
        [Met,InitConc,BkgdConc,repIndex] = InputReplicate(Met,InitConc,BkgdConc,nRep);
        ModelOptions.TimeStamp = repmat(ModelOptions.TimeStamp,nRep,1);
    end

%% MODEL RUN
% Now we call the model. Note this may take several minutes to run, depending on your system.
% Output will be saved in the "SavePath" above and will also be written to the structure S.

    S = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);
end
%% Corrected file:
%F0AM is set up to save without the variable repIndex. We re-save a new
%file with its inclusion with ending '_rep_save'

%Save missing variable repIndex as new MATLAB file
S(1).repInd = repIndex;
save_struct_name = strcat(full_savepath,'_rep_save','.mat');
save(save_struct_name, 'S')


% %Split up the run into each day of data
% SplitRun(S,'custom',S.repInd);
% 
% %Save each day's data as its own file
% day_split = {'S1', 'S2', 'S3', 'S4'};
% for i = 1:length(day_split)
%     save_runday_name = strcat(full_savepath,'_spinup_',day_split{i},'.mat');
%     save(save_runday_name,day_split{i})
% end

%%
%% 
%Save species concentrations as csv file for use in plotting via Python
speciesConc_table=struct2table(S.Conc);
species_conc_save_path = strcat(full_savepath,'_model_conc','.csv');
writetable(speciesConc_table, species_conc_save_path)
%%
chem_struct = S.Chem;
photolysis_rates_save_path = strcat(full_savepath,'_model_photolysis_rates','.csv');
writematrix(chem_struct.Rates,photolysis_rates_save_path);
%% 
transpose_rnames = transpose(chem_struct.Rnames);
rnames_table = cell2table(transpose_rnames);
photolysis_names_save_path = strcat(full_savepath,'_model_photolysis_rxn_names','.csv');
writetable(rnames_table,photolysis_names_save_path)

%%
spinup_path_load = strcat(full_savepath,'_spinup.mat');
spinup_var = importdata(spinup_path_load);
spinup_speciesConc_table=struct2table(spinup_var.Conc);
spinup_conc_save_path = strcat(full_savepath,'_spinup_conc','.csv');
writetable(spinup_speciesConc_table, spinup_conc_save_path)