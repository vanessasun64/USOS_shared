%USOS Conversion
savedir = 'C:\Users\u1545774\Documents\GitHub\USOS_shared\F0AM-4.3.0.1\Runs\'
runname_str = strcat('USOS','_',num2str(mon), '_', num2str(dy),'_', num2str(yr));

usos_yaml_save_path = strcat(savedir,runname_str,'\',runnumber,'.yaml');
USOS.time_UTC=cellstr(USOS.time_UTC);
USOS.time_MST=cellstr(USOS.time_MST);

WriteYaml(yaml_save_path,USOS);

%%
% savedir = 'C:\Users\u1545774\Documents\GitHub\USOS_shared\F0AM-4.3.0.1\Runs\'
% runname_str = strcat('USOS_8_6_2024');
% 
% usos_yaml_save_path = strcat(savedir,runname_str,'\',runname,'.yaml');
% USOS.time_UTC=cellstr(USOS.time_UTC);
% USOS.time_MST=cellstr(USOS.time_MST);
% 
% WriteYaml(yaml_save_path,USOS)
