clear all 

datapath_mac='/Users/vanessasun/Documents/phd/utah/research/USOS_shared/F0AM-4.3.0.1/Campaign_Data/matlab_merge/';
datapath_win = 'C:\Users\u1545774\Documents\GitHub\USOS_shared\F0AM-4.3.0.1\Campaign_Data\matlab_merge\';
parkeddir_mac = 'parked/';
parkeddir_win = 'parked\';
drivingdir_mac = 'driving/';
drivingdir_win = 'driving\';
originaldir_mac = 'original/';
originaldir_win = 'original\';
correcteddir_mac = 'corrected/';
correcteddir_win = 'corrected\';
% originalmerge = 'original_merge/';
matlabstructs_mac = 'matlab_structs/';
matlabstructs_win = 'matlab_structs\';

d=dir(fullfile(datapath_mac,parkeddir_mac,originaldir_mac,'*.mat')); %select parkeddir or drivingdir as needed
for i=1:numel(d)
  mat_file = fullfile(datapath_mac,parkeddir_mac,originaldir_mac,d(i).name); %select parkeddir or drivingdir as needed
  filename_with_extension = d(i).name;
  filename_without_extension = extractBefore(filename_with_extension, ".");
  
  load(mat_file);
  disp(['Now loading file:', mat_file]);
  cols=fieldnames(USOS); 
  for nm=1:length(cols)
    c=cols(nm);
    if isa(USOS.(c{1}),'char') ==1 % Turn char arrays into cell arrays. 
        USOS.(c{1})=cellstr(USOS.(c{1}));
    end
    if size(USOS.(c{1}))==[1,length(USOS.(c{1}))] % Make everything Nx1 not 1xN 
        USOS.(c{1})=USOS.(c{1})'; % So that all our indicies line up later... 
    end
  end
  add_closing = '_struct_for_MATLAB.mat';
  save(strcat(datapath_mac,parkeddir_mac,correcteddir_mac, filename_without_extension,add_closing),'USOS'); %select parkeddir or drivingdir as needed
end

clear all 