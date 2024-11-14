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
%%
save_CoM =[];
save_velocity = [];
save_acc = [];
save_UpArm_l_acc = [];
save_UpArm_r_acc = [];
save_LoArm_l_acc = [];
save_LoArm_r_acc = [];
save_VE_max = [];
save_VE_total = [];
save_LE_total = [];
save_LE_max = [];
save_Knee_angle_l = [];
save_Knee_angle_r = [];
save_Elbow_angle_l = [];
save_Elbow_angle_r = [];
save_distance = [];

%% Values for CoM

                % Head Shoul UpArm LoArm Hand  Trunk PelvW Thigh Shank Foot
Segment_Length = [0.13 0.129 0.186 0.146 0.108 0.288 0.045 0.245 0.246 0.039 0.055 0.152]; % könnte man noch indivduaksiieren

                % Hand LoArm UpArm Foot Shank Thigh Pelvis Trunk
CoM_Segment_p = [0.506 0.430 0.436 0.50 0.433 0.433 0.105 0.50];
CoM_Segment_d = [0.494 0.570 0.564 0.50 0.567 0.567 0.895 0.50];

Radius_of_Marker = 1.2;
sizeArray = size(out); 

list_of_markers = {'C7'	'LSHO'	'RSHO'	'RBAK'	'CLAV'	'STRN'	'T10'	'SAC'	'RUPA'	'RELB'	'RFRM'	'RWRB'	'RWRA'	'RFIN'	'LUPA'	'LELB'	'LFRM'	'LWRB'	'LWRA'	'LFIN'	'LASI'	'RASI'	'LTHI'	'RTHI'	'RKNE'	'LKNE'	'RTIB'	'LTIB'	'RANK'	'LANK'	'RTOE'	'LTOE'	'RHEE'	'LHEE'};

%demo = readtable('zzz_demographische_daten_SocCogxBiom.xlsx');
Height = [0,0,0,0,0,0,0, 195, 180 , 175];
%Height = readtable("heights.csv");
Weight =  [0,0,0,0,0,0,0,68,72,72];
%Weight = readtable("weights.csv");
fc = 1.0;  
Fs = 200;                             % Sampling Frequency (Hz)
%%
for idxPart =  9: 9%sizeArray(2)
    % for idxTrial = 1:length(out{idxPart}) 
    for idxTrial = 2:sizeArray(1)      
        % if isempty(out{idxTrial, idxPart}) == 1 %|| idxPart == 1 && idxTrial == 13 || idxPart == 2 && idxTrial == 11 || idxPart == 5 && idxTrial == 2 || idxPart == 5 && idxTrial == 7 || idxPart == 5 && idxTrial == 9 || idxPart == 5 && idxTrial == 10 || idxPart == 5 && idxTrial == 13 || idxPart == 5 && idxTrial == 15 || idxPart == 11 && idxTrial == 3 || idxPart == 11 && idxTrial == 15 || idxPart == 22 && idxTrial == 16 || idxPart == 23 && idxTrial == 1 || idxPart == 23 && idxTrial == 3 || idxPart == 23 && idxTrial == 8 || idxPart == 28 && idxTrial == 2 || idxPart == 41 && idxTrial == 2|| idxPart == 41 && idxTrial == 6|| idxPart == 41 && idxTrial == 7 ||idxPart == 2 && idxTrial == 7 || idxPart == 7 && idxTrial == 9 || idxPart == 8 && idxTrial == 4 || idxPart == 11 && idxTrial == 7 || idxPart == 22 && idxTrial == 16  || idxPart == 22 && idxTrial == 17 || idxPart == 34 && idxTrial == 3 || idxPart == 5 && idxTrial == 4
            % continue
        %else
            markerset = [];
            for idxMarker = 1:42
                marker_tmp = out{idxTrial, idxPart}(idxMarker,:,:);
                marker = squeeze(marker_tmp);
          
                %% Velocity
                marker_tmp = marker;
                marker_tmp(isnan(marker_tmp)) = 0;
                [b,ba] = butter(4,fc/(Fs/2));
                marker_filtered = filtfilt(b,ba,marker_tmp');
                %marker_filtered(iszero(marker_tmp)) = ;
                markerset = [markerset; marker_filtered'];
            end  
 
            %calculate offset
               % for j = 8: 8
                   for i = 1:size(out_fd,1)
                    if isempty(out_fd{i,idxPart})
                        continue
                    end
                    M = mean(out_fd{i,idxPart});
                   
                    out_fd_off{i,idxPart}(:,4) = out_fd{i,idxPart}(:,4)  - (out_fd{i,idxPart}(7500,4)); %Weight(i-2)*9.81;
                    out_fd_off{i,idxPart}(:,9) = out_fd{i,idxPart}(:,9)  - (out_fd{i,idxPart}(5,9)); %Weight(i-2)*9.81; 5
                    out_fd_off{i,idxPart}(:,14) = out_fd{i,idxPart}(:,14)  - (out_fd{i,idxPart}(5,14)); %Weight(i-2)*9.81; 5
                    out_fd_off{i,idxPart}(:,19) = out_fd{i,idxPart}(:,19)  - (out_fd{i,idxPart}(7500,19)); %Weight(i-2)*9.81;
                    %figure()
                    %title([num2str(i),num2str(idxPart)])
                    
                   
                    end
           % end
%%
Height_of_Participant = Height(idxPart);
for idxSegment = 1: length(Segment_Length)
    Segment_Length_adj(idxSegment) = Segment_Length(idxSegment)*Height_of_Participant(1,1);
end
%% 

CoM_UpArm = Segment_Length_adj(3) * CoM_Segment_p(3);
CoM_LoArm = Segment_Length_adj(4) * CoM_Segment_p(2);
CoM_Thigh = Segment_Length_adj(8) * CoM_Segment_p(6);
CoM_Shank = Segment_Length_adj(9) * CoM_Segment_p(5);

CoM_Hand = Segment_Length_adj(5) * CoM_Segment_p(1);
CoM_Foot = Segment_Length_adj(10) * CoM_Segment_p(4);
CoM_Pelvis = Segment_Length_adj(7) * CoM_Segment_p(7);
CoM_Trunk = Segment_Length_adj(6) * CoM_Segment_p(8);

%%
CoM_Trunk_xyz = [];
CoM_Trunk_xyz(1,:) = (markerset(45,:)+markerset(17,:)+markerset(101,:)+markerset(97,:))/4; % muss ich anpassen
CoM_Trunk_xyz(2,:) = (markerset(46,:)+markerset(18,:)+markerset(102,:)+markerset(98,:))/4;
CoM_Trunk_xyz(3,:) = (markerset(47,:)+markerset(19,:)+markerset(103,:)+markerset(99,:))/4;
% plot3(CoM_Trunk_xyz(1,:), CoM_Trunk_xyz(2,:), CoM_Trunk_xyz(3,:))
% axis equal
%% CoM UpperArm_r
RSHO = [];
RELB = [];
UpArm_r = [];
RSHO(1,:) = (markerset(45,:));%Top marker
RSHO(2,:) = (markerset(46,:));
RSHO(3,:) = (markerset(47,:));
RELB(1,:) = (markerset(57,:));
RELB(2,:) = (markerset(58,:));
RELB(3,:) = (markerset(59,:));
%
UpArm_r(1,:) = RELB(1,:)-RSHO(1,:);
UpArm_r(2,:) = RELB(2,:)-RSHO(2,:);
UpArm_r(3,:) = RELB(3,:)-RSHO(3,:);

CoM_UpArm_r = CoM_Segment_p(3) * UpArm_r;
%% CoM UpperArm_l
% CoM_LowArm_xyz_r(1,:) = (markerset(5,:)+markerset(9,:)+markerset(85,:)+markerset(89,:))/4;
% CoM_LowArm_xyz_r(2,:) = (markerset(6,:)+markerset(10,:)+markerset(86,:)+markerset(90,:))/4;
% CoM_LowArm_xyz_r(3,:) = (markerset(7,:)+markerset(11,:)+markerset(87,:)+markerset(91,:))/4;
LSHO = [];
LELB = [];
UpArm_l = [];
LSHO(1,:) = (markerset(17,:));
LSHO(2,:) = (markerset(18,:));
LSHO(3,:) = (markerset(19,:));
LELB(1,:) = (markerset(29,:));
LELB(2,:) = (markerset(30,:));
LELB(3,:) = (markerset(31,:));
%
UpArm_l(1,:) = LELB(1,:)-LSHO(1,:);
UpArm_l(2,:) = LELB(2,:)-LSHO(2,:);
UpArm_l(3,:) = LELB(3,:)-LSHO(3,:);

CoM_UpArm_l = CoM_Segment_p(3) * UpArm_l;

%% CoM LowerArm_l
% CoM_UpArm_xyz_l(1,:) = (markerset(5,:)+markerset(9,:)+markerset(85,:)+markerset(89,:))/4;
% CoM_UpArm_xyz_l(2,:) = (markerset(6,:)+markerset(10,:)+markerset(86,:)+markerset(90,:))/4;
% CoM_UpArm_xyz_l(3,:) = (markerset(7,:)+markerset(11,:)+markerset(87,:)+markerset(91,:))/4;
LWRI = [];
LoArm_l = [];

LWRI(1,:) = (markerset(33,:)); % wrist outside
LWRI(2,:) = (markerset(34,:));
LWRI(3,:) = (markerset(35,:));

%
LoArm_l(1,:) = LWRI(1,:)-LELB(1,:);
LoArm_l(2,:) = LWRI(2,:)-LELB(2,:);
LoArm_l(3,:) = LWRI(3,:)-LELB(3,:);

CoM_LoArm_l = CoM_Segment_p(2) * LoArm_l;
%% CoM LowerArm_r
% CoM_LowArm_xyz_l(1,:) = (markerset(5,:)+markerset(9,:)+markerset(85,:)+markerset(89,:))/4;
% CoM_LowArm_xyz_l(2,:) = (markerset(6,:)+markerset(10,:)+markerset(86,:)+markerset(90,:))/4;
% CoM_LowArm_xyz_l(3,:) = (markerset(7,:)+markerset(11,:)+markerset(87,:)+markerset(91,:))/4;
RWRI = [];
LoArm_r = [];

RWRI(1,:) = (markerset(61,:));
RWRI(2,:) = (markerset(62,:));
RWRI(3,:) = (markerset(63,:));

%
LoArm_r(1,:) = RWRI(1,:)-RELB(1,:);
LoArm_r(2,:) = RWRI(2,:)-RELB(2,:);
LoArm_r(3,:) = RWRI(3,:)-RELB(3,:);

CoM_LoArm_r = CoM_Segment_p(2) * LoArm_r;
%% CoM Thigh_r
% CoM_Thigh_xyz_r(1,:) = (markerset(5,:)+markerset(9,:)+markerset(85,:)+markerset(89,:))/4;
% CoM_Thigh_xyz_r(2,:) = (markerset(6,:)+markerset(10,:)+markerset(86,:)+markerset(90,:))/4;
% CoM_Thigh_xyz_r(3,:) = (markerset(7,:)+markerset(11,:)+markerset(87,:)+markerset(91,:))/4;
RASI = [];
RKNE = [];
Thigh_r = [];
RASI(1,:) = (markerset(101,:));
RASI(2,:) = (markerset(102,:));
RASI(3,:) = (markerset(103,:));
RKNE(1,:) = (markerset(133,:));
RKNE(2,:) = (markerset(134,:));
RKNE(3,:) = (markerset(135,:));
%
Thigh_r(1,:) = RKNE(1,:)-RASI(1,:);
Thigh_r(2,:) = RKNE(2,:)-RASI(2,:);
Thigh_r(3,:) = RKNE(3,:)-RASI(3,:);

CoM_Thigh_r = CoM_Segment_p(6) * Thigh_r;
%% CoM Shank_r
% CoM_Shank_xyz_r(1,:) = (markerset(5,:)+markerset(9,:)+markerset(85,:)+markerset(89,:))/4;
% CoM_Shank_xyz_r(2,:) = (markerset(6,:)+markerset(10,:)+markerset(86,:)+markerset(90,:))/4;
% CoM_Shank_xyz_r(3,:) = (markerset(7,:)+markerset(11,:)+markerset(87,:)+markerset(91,:))/4;
RANK = [];
Shank_r = [];

RANK(1,:) = (markerset(141,:));
RANK(2,:) = (markerset(142,:));
RANK(3,:) = (markerset(143,:));

%
Shank_r(1,:) = RANK(1,:)-RKNE(1,:);
Shank_r(2,:) = RANK(2,:)-RKNE(2,:);
Shank_r(3,:) = RANK(3,:)-RKNE(3,:);

CoM_Shank_r = CoM_Segment_p(5) * Shank_r;
%% CoM Thigh_l
% CoM_Thigh_xyz_l(1,:) = (markerset(5,:)+markerset(9,:)+markerset(85,:)+markerset(89,:))/4;
% CoM_Thigh_xyz_l(2,:) = (markerset(6,:)+markerset(10,:)+markerset(86,:)+markerset(90,:))/4;
% CoM_Thigh_xyz_l(3,:) = (markerset(7,:)+markerset(11,:)+markerset(87,:)+markerset(91,:))/4;
LASI = [];
LKNE = [];
Thigh_l = [];

LASI(1,:) = (markerset(97,:));
LASI(2,:) = (markerset(98,:));
LASI(3,:) = (markerset(99,:));
LKNE(1,:) = (markerset(109,:));
LKNE(2,:) = (markerset(110,:));
LKNE(3,:) = (markerset(111,:));
%
Thigh_l(1,:) = LKNE(1,:)-LASI(1,:);
Thigh_l(2,:) = LKNE(2,:)-LASI(2,:);
Thigh_l(3,:) = LKNE(3,:)-LASI(3,:);

CoM_Thigh_l = CoM_Segment_p(6) * Thigh_l;
%% CoM Shank_l
% CoM_Shank_xyz_l(1,:) = (markerset(5,:)+markerset(9,:)+markerset(85,:)+markerset(89,:))/4;
% CoM_Shank_xyz_l(2,:) = (markerset(6,:)+markerset(10,:)+markerset(86,:)+markerset(90,:))/4;
% CoM_Shank_xyz_l(3,:) = (markerset(7,:)+markerset(11,:)+markerset(87,:)+markerset(91,:))/4;
LANK = [];
Shank_l = [];

LANK(1,:) = (markerset(117,:));
LANK(2,:) = (markerset(118,:));
LANK(3,:) = (markerset(119,:));

%
Shank_l(1,:) = LANK(1,:)-LKNE(1,:);
Shank_l(2,:) = LANK(2,:)-LKNE(2,:);
Shank_l(3,:) = LANK(3,:)-LKNE(3,:);

CoM_Shank_l = CoM_Segment_p(5) * Shank_l;


%%
rel_CoM_UpArm_r = RSHO + CoM_UpArm_r;
rel_CoM_UpArm_l = LSHO + CoM_UpArm_l;
rel_CoM_LoArm_r = RELB + CoM_UpArm_r;
rel_CoM_LoArm_l = LELB + CoM_UpArm_l;

rel_CoM_Thigh_r = RASI + CoM_Thigh_r;
rel_CoM_Thigh_l = LASI + CoM_Thigh_l;
rel_CoM_Shank_r = RKNE + CoM_Shank_r;
rel_CoM_Shank_l = LKNE + CoM_Shank_l;

%% werte anpassen, da Kopf bei mir noch da (was ist Headset ?)
m = Weight(idxPart);
CoM = (rel_CoM_UpArm_r * 0.0318906606 *m + rel_CoM_UpArm_l * 0.0318906606 *m + rel_CoM_LoArm_l * 0.0182232346 *m+ rel_CoM_LoArm_r * 0.0182232346 *m+ rel_CoM_Shank_l * 0.0529612756 *m+ rel_CoM_Shank_r * 0.0529612756 *m+ rel_CoM_Thigh_l * 0.113895216 *m+ rel_CoM_Thigh_r * 0.113895216 *m+ CoM_Trunk_xyz * 0.566059226*m)/m;
vel_tmp = diff(CoM(1,:))*Fs/1000;

%% Cutting of longer trials
cut_velocity = 0.19;
min_walk_time = 1;
disp("Hallo");
foot_middle = [];
if contains(list_of_names(idxTrial, idxPart), "baseline")
    continue
    % eventuell "leere" function einfügen
else
    disp(list_of_names(idxTrial, idxPart));
    foot_middle(1,:) = ((squeeze(skelt_out{idxTrial, idxPart}(1,20,:)) + squeeze(skelt_out{idxTrial, idxPart}(1,24,:)))/2)';
    foot_middle(2,:) = ((squeeze(skelt_out{idxTrial, idxPart}(2,20,:)) + squeeze(skelt_out{idxTrial, idxPart}(2,24,:)))/2)';    
    foot_middle(3,:) = ((squeeze(skelt_out{idxTrial, idxPart}(3,20,:)) + squeeze(skelt_out{idxTrial, idxPart}(3,24,:)))/2)';

end

    end
end
%%
clf
plot(out_fd_off{33, 8}(:,4))
hold on
plot(out_fd_off{33, 8}(:,9))
plot(out_fd_off{33, 8}(:,14))
plot(out_fd_off{33, 8}(:,19))
shg
%% test trails schneiden
ms = CoM(1,:);
Df = diff(ms)
plot(ms)
hold on
plot(Df*200)
%%
plot(CoM(1,:),CoM(2,:))
hold on 
plot(foot_middle(1,:), foot_middle(2,:))
%%
diff_C_FM1 = CoM(1,:)- foot_middle(1,:);
diff_C_FM2 = CoM(2,:)- foot_middle(2,:);

plot(diff_C_FM1)
hold on
plot(diff_C_FM2)
%%
euclid = sqrt(diff_C_FM1.^2 + diff_C_FM2.^2);
first = true;
trials = [];
if (first)
    walking = euclid > 30;
    start = find(walking);
    disp(start(1,1));
   % e = CoM(1,start(1,1):end);
    ends = CoM(1,:) < CoM(1, start(1,1));
    e = find(ends);
    disp(e(start(1,1)+1));
    ttt = CoM(1,start(1,1):e(start(1,1)+1));
    trials{1} = ttt;
    first = false;
%else?
    cut_e = euclid(1, e(start(1,1)+1):end) > 50;
    cut_c = CoM(:, e(start(1,1)+1):end);
    start2 = find(cut_e);
    disp(start2(1,1));
    ends2 = cut_c(1,:) < cut_c(1, start2(1,1));
    e2 = find(ends2);
    disp(e2(start2(1,1)+1));
    ttt2 = cut_c(1, start2(1,1):e2(start2(1,1)+1));
    trials{2} = ttt2;
end

