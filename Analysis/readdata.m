%% MT_GS x Biom
% Authors: MS, GS
% 28.10.2024
clc
clear all

dirData = uigetdir; % open the folder where the data of the participants is in
myFiles = dir(fullfile(dirData)); % creates

for idxParticipant = 8: 19 %length(myFiles) % starting at 3 because there are additional files in (at 1 and 2) we want to skip
    stopath = fullfile(myFiles(idxParticipant).folder,myFiles(idxParticipant).name); % this line creates the path (directory) to current participant
    directory = dir(stopath);
    for idxConditions = 1:length(dir(stopath))
            matFiles = dir(fullfile(stopath, directory(idxConditions).name, '*.mat')); % get all the files of that participant with extension _pos_global.sto
            txtFiles = dir(fullfile(stopath, directory(idxConditions).name, '*.txt'));
            for idxFile = 1 : length(matFiles)
               matfile = fullfile(matFiles(idxFile).folder, matFiles(idxFile).name); % get current file
               eval('a = importdata (matfile)');
               list_of_names{idxFile, idxParticipant-2} = matFiles(idxFile).name;
               out{idxFile, idxParticipant-2} = a.Trajectories.Labeled.Data;
               skelt_out{idxFile, idxParticipant-2} = a.Skeletons.PositionData;
            end
            for idxTxtFile = 1:length(txtFiles)
               txtFile = fullfile(txtFiles(idxTxtFile).folder, txtFiles(idxTxtFile).name); % get current file
               aa = readmatrix (txtFile);
               list_of_txt_names{idxTxtFile, idxParticipant-2} = txtFiles(idxTxtFile).name;
               out_fd{idxTxtFile, idxParticipant-2} = aa(18:end,:);
            end
    end
end
out_fd_off = out_fd; 