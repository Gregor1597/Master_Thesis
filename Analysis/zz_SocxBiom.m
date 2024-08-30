%% SocCog x Biom
% Authors: MS, NeS, PaF, JaK
% 07.08.2024

clc
clear all

dirData = uigetdir; % open the folder where the data of the participants is in
myFiles = dir(fullfile(dirData)); % creates

for idxParticipant = 3 : length(myFiles) % starting at 3 because there are additional files in (at 1 and 2) we want to skip
    stopath = fullfile(myFiles(idxParticipant).folder,myFiles(idxParticipant).name); % this line creates the path (directory) to current participant
    directory = dir(stopath);
    for idxConditions = 1:length(dir(stopath))
            matFiles = dir(fullfile(stopath, directory(idxConditions).name, '*.mat')); % get all the files of that participant with extension _pos_global.sto
            for idxFile = 1 : length(matFiles)
                matfile = fullfile(matFiles(idxFile).folder, matFiles(idxFile).name); % get current file
                eval('a = importdata (matfile)');
                list_of_names{idxFile, idxParticipant-2} = matFiles(idxFile).name;
                out{idxFile, idxParticipant-2} = a.Trajectories.Labeled.Data;                 
            end
    end
end
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

sizeArray = size(out); 
list_of_markers = {'C7'	'LSHO'	'RSHO'	'RBAK'	'CLAV'	'STRN'	'T10'	'SAC'	'RUPA'	'RELB'	'RFRM'	'RWRB'	'RWRA'	'RFIN'	'LUPA'	'LELB'	'LFRM'	'LWRB'	'LWRA'	'LFIN'	'LASI'	'RASI'	'LTHI'	'RTHI'	'RKNE'	'LKNE'	'RTIB'	'LTIB'	'RANK'	'LANK'	'RTOE'	'LTOE'	'RHEE'	'LHEE'};

%demo = readtable('zzz_demographische_daten_SocCogxBiom.xlsx');
%Height = [demo.DD05_01];
Height = readtable("heights.csv");
%Weight = [demo.DD04_01];
Weight = readtable("weights.csv");
for idxPart =  1: sizeArray(2)
    % for idxTrial = 1:length(out{idxPart}) 
    for idxTrial = 2:sizeArray(1)      
        if isempty(out{idxTrial, idxPart}) == 1 %|| idxPart == 1 && idxTrial == 13 || idxPart == 2 && idxTrial == 11 || idxPart == 5 && idxTrial == 2 || idxPart == 5 && idxTrial == 7 || idxPart == 5 && idxTrial == 9 || idxPart == 5 && idxTrial == 10 || idxPart == 5 && idxTrial == 13 || idxPart == 5 && idxTrial == 15 || idxPart == 11 && idxTrial == 3 || idxPart == 11 && idxTrial == 15 || idxPart == 22 && idxTrial == 16 || idxPart == 23 && idxTrial == 1 || idxPart == 23 && idxTrial == 3 || idxPart == 23 && idxTrial == 8 || idxPart == 28 && idxTrial == 2 || idxPart == 41 && idxTrial == 2|| idxPart == 41 && idxTrial == 6|| idxPart == 41 && idxTrial == 7 ||idxPart == 2 && idxTrial == 7 || idxPart == 7 && idxTrial == 9 || idxPart == 8 && idxTrial == 4 || idxPart == 11 && idxTrial == 7 || idxPart == 22 && idxTrial == 16  || idxPart == 22 && idxTrial == 17 || idxPart == 34 && idxTrial == 3 || idxPart == 5 && idxTrial == 4
            continue
        else
            markerset = [];
            for idxMarker = 1:34
                marker_tmp = out{idxTrial, idxPart}(idxMarker,:,:);
                marker = squeeze(marker_tmp);
                markerset = [markerset; marker];
            end         

%% Plots for testing
% hold on
% for x = 1:4:136
%   plot(markerset(x+1,:), markerset(x,:))
%   axis equal
% end
% shg
%% Values for CoM

                % Head Shoul UpArm LoArm Hand  Trunk PelvW Thigh Shank Foot
Segment_Length = [0.13 0.129 0.186 0.146 0.108 0.288 0.045 0.245 0.246 0.039 0.055 0.152]; % könnte man noch indivduaksiieren

                % Hand LoArm UpArm Foot Shank Thigh Pelvis Trunk
CoM_Segment_p = [0.506 0.430 0.436 0.50 0.433 0.433 0.105 0.50];
CoM_Segment_d = [0.494 0.570 0.564 0.50 0.567 0.567 0.895 0.50];

Radius_of_Marker = 1.2;
%% Segment Lengths


Height_of_Participant = cell2mat(Height(idxPart));
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
% LASI(1) = nanmean(markerset(81,:));
% LASI(2) = nanmean(markerset(82,:));
% LASI(3) = nanmean(markerset(83,:));
% RASI(1) = nanmean(markerset(85,:));
% RASI(2) = nanmean(markerset(86,:));
% RASI(3) = nanmean(markerset(87,:));
%%
%ASI_distance = abs(norm(LASI - RASI));

%% CoM Trunk
CoM_Trunk_xyz = [];
CoM_Trunk_xyz(1,:) = (markerset(5,:)+markerset(9,:)+markerset(81,:)+markerset(85,:))/4; % muss ich anpassen
CoM_Trunk_xyz(2,:) = (markerset(6,:)+markerset(10,:)+markerset(82,:)+markerset(86,:))/4;
CoM_Trunk_xyz(3,:) = (markerset(7,:)+markerset(11,:)+markerset(83,:)+markerset(87,:))/4;
% plot3(CoM_Trunk_xyz(1,:), CoM_Trunk_xyz(2,:), CoM_Trunk_xyz(3,:))
% axis equal
%% CoM UpperArm_r
RSHO = [];
RELB = [];
UpArm_r = [];
RSHO(1,:) = (markerset(9,:));
RSHO(2,:) = (markerset(10,:));
RSHO(3,:) = (markerset(11,:));
RELB(1,:) = (markerset(37,:));
RELB(2,:) = (markerset(38,:));
RELB(3,:) = (markerset(39,:));
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
LSHO(1,:) = (markerset(5,:));
LSHO(2,:) = (markerset(6,:));
LSHO(3,:) = (markerset(7,:));
LELB(1,:) = (markerset(61,:));
LELB(2,:) = (markerset(62,:));
LELB(3,:) = (markerset(63,:));
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

LWRI(1,:) = (markerset(69,:));
LWRI(2,:) = (markerset(70,:));
LWRI(3,:) = (markerset(71,:));

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

RWRI(1,:) = (markerset(45,:));
RWRI(2,:) = (markerset(46,:));
RWRI(3,:) = (markerset(47,:));

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
RASI(1,:) = (markerset(85,:));
RASI(2,:) = (markerset(86,:));
RASI(3,:) = (markerset(87,:));
RKNE(1,:) = (markerset(97,:));
RKNE(2,:) = (markerset(98,:));
RKNE(3,:) = (markerset(99,:));
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

RANK(1,:) = (markerset(113,:));
RANK(2,:) = (markerset(114,:));
RANK(3,:) = (markerset(115,:));

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

LASI(1,:) = (markerset(81,:));
LASI(2,:) = (markerset(82,:));
LASI(3,:) = (markerset(83,:));
LKNE(1,:) = (markerset(101,:));
LKNE(2,:) = (markerset(102,:));
LKNE(3,:) = (markerset(103,:));
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

%% CoM
% 028/100
% /87.8
% t = rel_CoM_UpArm_l *0.0318906606 *m;
% plot3(t(1,:), t(2,:), t(3,:))
%%
t=rel_CoM_UpArm_l*0.05;

% plot(rel_CoM_UpArm_l(1,:)) 
% hold on
% plot(t(1,:), 'LineStyle','-', 'Color','r', 'LineWidth',5)

%plot(rel_CoM_UpArm_l(1,:)*0.5)


%%

m_tmp = cell2mat(Weight(idxPart));
m = str2num(m_tmp);
%%
CoM = [];
%% werte anpassen, da Kopf bei mir noch da (was ist Headset ?)
CoM = (rel_CoM_UpArm_r * 0.0318906606 *m + rel_CoM_UpArm_l * 0.0318906606 *m + rel_CoM_LoArm_l * 0.0182232346 *m+ rel_CoM_LoArm_r * 0.0182232346 *m+ rel_CoM_Shank_l * 0.0529612756 *m+ rel_CoM_Shank_r * 0.0529612756 *m+ rel_CoM_Thigh_l * 0.113895216 *m+ rel_CoM_Thigh_r * 0.113895216 *m+ CoM_Trunk_xyz * 0.566059226*m)/m;

%% Plot of CoM per segment
% figure(1)
% hold on
% plot3(rel_CoM_LoArm_l(1,:), rel_CoM_LoArm_l(2,:), rel_CoM_LoArm_l(3,:))
% plot3(rel_CoM_UpArm_l(1,:), rel_CoM_UpArm_l(2,:), rel_CoM_UpArm_l(3,:))
% plot3(rel_CoM_Thigh_l(1,:), rel_CoM_Thigh_l(2,:), rel_CoM_Thigh_l(3,:))
% plot3(rel_CoM_Shank_l(1,:), rel_CoM_Shank_l(2,:), rel_CoM_Shank_l(3,:))
% 
% plot3(rel_CoM_LoArm_r(1,:), rel_CoM_LoArm_r(2,:), rel_CoM_LoArm_r(3,:))
% plot3(rel_CoM_UpArm_r(1,:), rel_CoM_UpArm_r(2,:), rel_CoM_UpArm_r(3,:))
% plot3(rel_CoM_Thigh_r(1,:), rel_CoM_Thigh_r(2,:), rel_CoM_Thigh_r(3,:))
% plot3(rel_CoM_Shank_r(1,:), rel_CoM_Shank_r(2,:), rel_CoM_Shank_r(3,:))
% plot3(CoM_Trunk_xyz(1,:), CoM_Trunk_xyz(2,:), CoM_Trunk_xyz(3,:))
% plot3(CoM(1,:), CoM(2,:), CoM(3,:), 'LineWidth', 2, 'Color', "r", "MarkerSize",8, "MarkerFaceColor", "r")
% axis equal
% shg

%% Velocity
CoM_time = CoM(1,:) * (1/460); 
% plot(CoM_time)
%%
samples = 0:length(CoM(1,:));         % Sample Indices Vector
Fs = 460;                             % Sampling Frequency (Hz)
t_tmp = samples/Fs; 
t = t_tmp(1,2:end);
% plot(t,CoM(1,:))
%%
%%
fc = 1.0;                            % check!
fs = 460;
CoM_tmp = CoM(1,:);
CoM_tmp(isnan(CoM_tmp)) = 0;
[b,ba] = butter(4,fc/(fs/2));            
CoM_filtered = filtfilt(b,ba,CoM_tmp);
%%
vel_tmp = diff(-CoM_filtered(1,:))*460/1000;
%eventuell drehen abhängig von Ausrichtung des Koordinatenssystems

vel_filtered= lowpass(vel_tmp, 4, 460);
% plot(t(2:end),vel_filtered)
%% Acceleration
acc = diff(diff(-CoM_filtered(1,:)));
% plot(t(3:end),acc_lp)

%freqz(b,a,[],fs)

%% Starting and Stopping
vel_with_zeroes = diff(-CoM_tmp);
Start = [];
Stop = [];
startpoint = 137;
for idxVel = startpoint+1:length(vel_with_zeroes)
    if vel_with_zeroes(idxVel-1) < 0.19 && vel_with_zeroes(idxVel) >= 0.19
        Start = [t(idxVel) vel_with_zeroes(idxVel)];
        for idxStart = (idxVel + 920) : length(vel_with_zeroes)
            if vel_with_zeroes(idxStart-1) > 0.19 && vel_with_zeroes(idxStart) <= 0.19
                Stop = [t(idxStart) vel_with_zeroes(idxStart)];
            end
        end
        break
    end
end
save_CoM{idxTrial-1,idxPart} = CoM(:,Start(1)*460:Stop(1)*460);

%%


for idxVel = 2:length(vel_filtered)
    if vel_filtered(idxVel-1) < 0.19 && vel_filtered(idxVel) >= 0.19
        Start = [t(idxVel) vel_filtered(idxVel)];
    end
    if vel_filtered(idxVel-1) > 0.19 && vel_filtered(idxVel) <= 0.19
        Stop = [t(idxVel) vel_filtered(idxVel)];
    end
end

%%
% %plot(acc, 'Color', [0.7 0.7 0.7])
% hold on
% %plot(acc*1000)
% plot(t(3:end), acc_lp)
% plot(Start(1), Start(2), 'Marker','diamond', 'Color','r')
% plot(Stop(1), Stop(2), 'Marker','diamond', 'Color','m')
% 
% %%plot(acc*100)
%% Distance to group
Distance_to_Group_Center = CoM(1,round(Stop(1)*460));
save_distance{idxTrial-1, idxPart} = Distance_to_Group_Center;
%% Lateral CoM excursion
Lateral_Excursion_total = CoM(2,round(Stop(1)*460)) - CoM(2,round(Start(1)*460));
Lateral_Excursion_max = max(abs(CoM(2,round(Start(1)*460):round(Stop(1)*460))-CoM(2,round(Start(1)*460))));
save_LE_total{idxTrial-1, idxPart} = Lateral_Excursion_total;
save_LE_max{idxTrial-1, idxPart} = Lateral_Excursion_max;
% plot(CoM(2,:))
% hold on
% plot(Start(1)*460, Start(2), 'Marker','diamond', 'Color','r')
% plot(Stop(1)*460, Stop(2), 'Marker','diamond', 'Color','m')
%% Vertical CoM excursion
Vertical_Excursion_total = CoM(3,round(Stop(1)*460)) - CoM(3,round(Start(1)*460));
Vertical_Excursion_max = max(abs(CoM(3,round(Start(1)*460):round(Stop(1)*460))-CoM(3,round(Start(1)*460))));
save_VE_total{idxTrial-1, idxPart} = Vertical_Excursion_total;
save_VE_max{idxTrial-1, idxPart} = Vertical_Excursion_max;
% plot(CoM(3,:))
% hold on
% plot(Start(1)*460, Start(2), 'Marker','diamond', 'Color','r')
% plot(Stop(1)*460, Stop(2), 'Marker','diamond', 'Color','m')
%% Touchdown / Takeoff

for idxTO = round(Start(1)*460):round(Stop(1)*460)
    RHEE(idxTO+1-round(Start(1)*460)) = -markerset(131,idxTO);
    LHEE(idxTO+1-round(Start(1)*460)) = -markerset(135,idxTO);
    RTOE(idxTO+1-round(Start(1)*460)) = -markerset(123,idxTO);
    LTOE(idxTO+1-round(Start(1)*460)) = -markerset(127,idxTO);
end
if length(RHEE) > round(Stop(1)*460)-1
    RHEE = RHEE(1:round(Stop(1)*460)-1);
end    
if length(LHEE) > round(Stop(1)*460)-1
    LHEE = LHEE(1:round(Stop(1)*460)-1);
end
if length(RTOE) > round(Stop(1)*460)-1
    RTOE = RTOE(1:round(Stop(1)*460)-1);
end    
if length(LTOE) > round(Stop(1)*460)-1
    LTOE = LTOE(1:round(Stop(1)*460)-1);
end    
%%

LFOO =(LTOE + LHEE) / 2;
RFOO =(RTOE + RHEE) / 2;

cutoff_freq_foot = 7.0;
[foot_b,foot_a] = butter(4,cutoff_freq_foot/(fs/2));            
LFOO_filt = filtfilt(foot_b,foot_a,LFOO(~isnan(LFOO)));
RFOO_filt = filtfilt(foot_b,foot_a,RFOO(~isnan(RFOO)));

LFOO_vel = -diff(LFOO_filt);
RFOO_vel = -diff(RFOO_filt);

%%

c = [];
d = [];
e = [];
f = [];
HS_dummye = [];
HS_dummyf = [];
HS_dummyc = [];
HS_dummyd = [];
[HS_rtmpx HS_rtmpy] = findpeaks(-RFOO_vel, 'MinPeakProminence', 0.5);
[NHS_rtmpx NHS_rtmpy] =  findpeaks(-RFOO_vel, 'MinPeakProminence', 1, 'MinPeakDistance', 350);

[HS_ltmpx HS_ltmpy] = findpeaks(-LFOO_vel, 'MinPeakProminence', 0.5);
[NHS_ltmpx NHS_ltmpy] =  findpeaks(-LFOO_vel, 'MinPeakProminence', 1, 'MinPeakDistance', 350);

for idxHS = 1:length(HS_ltmpx)
    for idxNHS = 1:length(NHS_ltmpx)
        if isequal(HS_ltmpx(idxHS), NHS_ltmpx(idxNHS))
            HS_ltmpx(idxHS) = 0;
            HS_ltmpy(idxHS) = 0;
        end
    end
end
for idxe = 1:length(HS_ltmpx)
    if HS_ltmpx(idxe) ~= 0
        HS_dummye(end+1) = HS_ltmpx(idxe); %
        HS_dummyf(end+1) = HS_ltmpy(idxe); %
    end   
end
for idxHSr = 1:length(HS_rtmpx)
    for idxNHSr = 1:length(NHS_rtmpx)
        if isequal(HS_rtmpx(idxHSr), NHS_rtmpx(idxNHSr))
            HS_rtmpx(idxHSr) = 0;
            HS_rtmpy(idxHSr) = 0;
        end
    end
end
for idxer = 1:length(HS_rtmpx)
    if HS_rtmpx(idxer) ~= 0
        HS_dummyc(end+1) = HS_rtmpx(idxer); % HS_dummy ehem c
        HS_dummyd(end+1) = HS_rtmpy(idxer); % HS_dummy ehem d
    end  

end

[c d] = findpeaks(RFOO_vel, 'MinPeakProminence', 1, 'MinPeakDistance', 350); % TO_ry TO_rx
[e f] = findpeaks(LFOO_vel, 'MinPeakProminence', 1, 'MinPeakDistance', 350); % TO_ly TO_lx
%%
if length(d) > length(f)
for idxD = 2:length(d)
    if idxD > length(f)+1
        break
    end
    if idxD > length(d)
        break
    end
   
    if d(idxD-1) < f(idxD-1) && d(idxD) < f(idxD-1)
        d(idxD-1) = [];
    end
    if length(d) > length(f)+1
        d(end) = [];
    end
end
else
  for idxD = 2:length(f)
  if idxD > length(d)+1
        break
    end
    if idxD > length(f)
        break
    end
    if f(idxD-1) < d(idxD-1) && f(idxD) < d(idxD-1)
        f(idxD-1) = [];
    end
    if length(f) > length(d)+1
        f(end) = [];
    end
  end
end
%%
StepLRight = [];
StepLLeft = [];
StepWRight = [];
StepWLeft = [];
if d(1) < f(1)

    for idxStep = 1:length(d)-1
        if idxStep > length(f)
            break
        end
        StepLRight(idxStep) = markerset(129, d(idxStep+1))-markerset(133, f(idxStep));
        StepWRight(idxStep) = markerset(130, d(idxStep+1))-markerset(134, f(idxStep));
    end
    for idxStep = 1:length(f)-1
        StepLLeft(idxStep) = markerset(133, f(idxStep))-markerset(129, d(idxStep));
        StepWLeft(idxStep) = markerset(134, f(idxStep))-markerset(130, d(idxStep));
    end
else 
    for idxStep = 1:length(f)-1
        if idxStep > length(d)
            break
        end
        StepLRight(idxStep) = markerset(133, f(idxStep+1))-markerset(129, d(idxStep));
        StepWRight(idxStep) = markerset(134, f(idxStep+1))-markerset(130, d(idxStep));
    end
    for idxStep = 1:length(d)-1
        StepLLeft(idxStep) = markerset(129, d(idxStep))-markerset(133, f(idxStep));
        StepWLeft(idxStep) = markerset(130, d(idxStep))-markerset(134, f(idxStep));
    end
end

save_StepL_Right{idxTrial-1, idxPart} = StepLRight;
save_StepL_Left{idxTrial-1, idxPart} = StepLLeft; 
save_Width_Right{idxTrial-1, idxPart} = StepWRight;
save_Wigth_Left{idxTrial-1, idxPart} = StepWLeft; 
%%
% figure(2)
% plot(RHEE)
% hold on
 %plot(LHEE)
 % hold on
 % plot(LFOO_vel)
 % plot(f(1), -e(1), 'Marker','diamond', 'Color','r')
 % plot(f(2), -e(2), 'Marker','diamond', 'Color','r')
 % plot(TO_lx(1), TO_ly(1), 'Marker', 'diamond', 'Color', 'b')
 % plot(TO_lx(2), TO_ly(2), 'Marker', 'diamond', 'Color', 'b')
 % plot(TO_lx(3), TO_ly(3), 'Marker', 'diamond', 'Color', 'b')

 % plot(f(3), e(3), 'Marker','diamond', 'Color','r')
 % plot(f(4), e(4), 'Marker','diamond', 'Color','r')
 % plot(f(5), e(5), 'Marker','diamond', 'Color','r')
% plot(f(6), e(6), 'Marker','diamond', 'Color','r')

%% Joint Angles (Knee, Elbow)
V0 = CoM_Shank_r;
V1 = CoM_Shank_l;
V2 = CoM_Thigh_r;
V3 = CoM_Thigh_l;

W0 = CoM_LoArm_r;
W1 = CoM_LoArm_l;
W2 = CoM_UpArm_r;
W3 = CoM_UpArm_l;

%plot(t,atan2d(norm(cross(u,v)), dot(u,v)))
for idx = 1:length(V0)
    Elbow_angle_r(idx) = acosd(dot(W0(:,idx),W2(:,idx)) / (norm(W0(:,idx)) * norm(W2(:,idx))) );
end
for idx = 1:length(V0)
    Elbow_angle_l(idx) = acosd(dot(W1(:,idx),W3(:,idx)) / (norm(W1(:,idx)) * norm(W3(:,idx))) );
end

for idx = 1:length(V0)
    Knee_angle_r(idx) = acosd(dot(V0(:,idx),V2(:,idx)) / (norm(V0(:,idx)) * norm(V2(:,idx))) );
end
for idx = 1:length(V0)
    Knee_angle_l(idx) = acosd(dot(V1(:,idx),V3(:,idx)) / (norm(V1(:,idx)) * norm(V3(:,idx))) );
end
save_velocity{idxTrial-1,idxPart} = vel_tmp(:,Start(1)*460:Stop(1)*460);
save_Elbow_angle_l{idxTrial-1, idxPart} = Elbow_angle_l(:,Start(1)*460:Stop(1)*460);
save_Elbow_angle_r{idxTrial-1, idxPart} = Elbow_angle_r(:,Start(1)*460:Stop(1)*460);
save_Knee_angle_l{idxTrial-1, idxPart} = Knee_angle_l(:,Start(1)*460:Stop(1)*460);
save_Knee_angle_r{idxTrial-1, idxPart} = Knee_angle_r(:,Start(1)*460:Stop(1)*460);

%%
% clf
% plot(t, Elbow_angle_r)
% hold on
% plot(t, Elbow_angle_l)
%% Acc Arm

UpArm_l_acc = diff(diff(-CoM_UpArm_l(1,:)));
UpArm_r_acc = diff(diff(-CoM_UpArm_r(1,:)));
LoArm_l_acc = diff(diff(-CoM_LoArm_l(1,:)));
LoArm_r_acc = diff(diff(-CoM_LoArm_r(1,:)));

save_acc{idxTrial-1, idxPart} = acc(:,Start(1)*460:Stop(1)*460-1);

save_UpArm_l_acc{idxTrial-1, idxPart} = UpArm_l_acc(:,Start(1)*460:Stop(1)*460-1);
save_UpArm_r_acc{idxTrial-1, idxPart} = UpArm_r_acc(:,Start(1)*460:Stop(1)*460-1);
save_LoArm_l_acc{idxTrial-1, idxPart} = LoArm_l_acc(:,Start(1)*460:Stop(1)*460-1);
save_LoArm_r_acc{idxTrial-1, idxPart} = LoArm_r_acc(:,Start(1)*460:Stop(1)*460-1);
%plot(t(3:end),UpArm_l_acc(1,:))
%% Eye Movement
        end
    end
    text = ['Participant ', num2str(idxPart), ' done'];
    disp(text)
end
%%
save zz_Results_SocCogxBiom.mat save_acc save_velocity save_LoArm_r_acc save_LoArm_l_acc save_UpArm_r_acc save_UpArm_l_acc save_Knee_angle_r save_Knee_angle_l save_Elbow_angle_r save_Elbow_angle_l save_VE_max save_VE_total save_LE_max save_LE_total save_distance save_CoM save_StepL_Left save_StepL_Right save_Wigth_Left save_Width_Right -mat
%%

