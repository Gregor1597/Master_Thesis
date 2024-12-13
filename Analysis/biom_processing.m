% %% MT_GS x Biom
% % Authors: MS, GS
% % 28.10.2024
%% Values for CoM

                % Head Shoul UpArm LoArm Hand  Trunk PelvW Thigh Shank Foot
Segment_Length = [0.13 0.129 0.186 0.146 0.108 0.288 0.045 0.245 0.246 0.039 0.055 0.152]; % könnte man noch indivduaksiieren

                % Hand LoArm UpArm Foot Shank Thigh Pelvis Trunk
CoM_Segment_p = [0.506 0.430 0.436 0.50 0.433 0.433 0.105 0.50];
CoM_Segment_d = [0.494 0.570 0.564 0.50 0.567 0.567 0.895 0.50];

Radius_of_Marker = 1.2;
sizeArray = size(out); 
results_r = struct();
results_l = struct();
list_of_markers = {'C7'	'LSHO'	'RSHO'	'RBAK'	'CLAV'	'STRN'	'T10'	'SAC'	'RUPA'	'RELB'	'RFRM'	'RWRB'	'RWRA'	'RFIN'	'LUPA'	'LELB'	'LFRM'	'LWRB'	'LWRA'	'LFIN'	'LASI'	'RASI'	'LTHI'	'RTHI'	'RKNE'	'LKNE'	'RTIB'	'LTIB'	'RANK'	'LANK'	'RTOE'	'LTOE'	'RHEE'	'LHEE'};

Height = readtable("heights.csv");
Weight = readtable("weights.csv");
fc = 10;  
Fs = 200; % Sampling Frequency (Hz)
counter = 1;
dt = 1 / 200; % Zeitintervall in Sekunden (200 Hz)
for idxPart = 8: 17%sizeArray(2)
    
    counter = counter +1;
    for idxTrial = 1:9
            markerset = [];
            skeleton = squeeze(skelt_out{idxTrial, idxPart});
            skeleton(isnan(skeleton)) = 0;
            if idxPart == 14
                m_amount = 41;
            else 
                m_amount = 42;
            end
            first_iter= true;
            for idxMarker = 1:m_amount
                marker_tmp = out{idxTrial, idxPart}(idxMarker,:,:);
                marker = squeeze(marker_tmp);
                %probleme die ich indviduell behebe
                
                if idxPart == 16 && idxTrial == 9
                    if(first_iter)
                        skeleton = skeleton(:,:,8000:end);
                      
                        first_iter=false;
                    end
                    marker = marker(:,8000:end);
                    
                end
                if idxPart == 16 && idxTrial == 7
                    if(first_iter)
                        skeleton = skeleton(:,:,850:end);
              
                        first_iter=false;
                    end
                    marker = marker(:,850:end);
                end
                %% filter
                marker_tmp = marker;
                marker_tmp(isnan(marker_tmp)) = 0;
                [b,ba] = butter(4,2/(Fs/2));
                marker_filtered = filtfilt(b,ba,marker_tmp');
                %marker_filtered(iszero(marker_tmp)) = ;
                markerset = [markerset; marker_filtered'];
            end  
 
            %calculate offset
                   for i = 1:size(out_fd,1)
                    if isempty(out_fd{i,idxPart})
                        continue
                    end
                    M = mean(out_fd{i,idxPart});
                   
                    out_fd_off{i,idxPart}(:,4) = out_fd{i,idxPart}(:,4)  - (out_fd{i,idxPart}(7500,4)); %Weight(i-2)*9.81;
                    out_fd_off{i,idxPart}(:,9) = out_fd{i,idxPart}(:,9)  - (out_fd{i,idxPart}(5,9)); %Weight(i-2)*9.81; 5
                    out_fd_off{i,idxPart}(:,14) = out_fd{i,idxPart}(:,14)  - (out_fd{i,idxPart}(5,14)); %Weight(i-2)*9.81; 5
                    out_fd_off{i,idxPart}(:,19) = out_fd{i,idxPart}(:,19)  - (out_fd{i,idxPart}(7500,19)); %Weight(i-2)*9.81;  
                   
                    end
           % end
%%
Height_of_Participant = table2array(Height(idxPart-7,1));
for idxSegment = 1: length(Segment_Length)
    Segment_Length_adj(idxSegment) = Segment_Length(idxSegment)*Height_of_Participant(1,1);
end


%%
%sonderfall participant 14
sond_cond = 0;
if(idxPart==14)
    sond_cond = 4;
end
CoM_Trunk_xyz = [];
CoM_Trunk_xyz(1,:) = (markerset(45-sond_cond,:)+markerset(17-sond_cond,:)+markerset(101-sond_cond,:)+markerset(97-sond_cond,:))/4; % muss ich anpassen
CoM_Trunk_xyz(2,:) = (markerset(46-sond_cond,:)+markerset(18-sond_cond,:)+markerset(102-sond_cond,:)+markerset(98-sond_cond,:))/4;
CoM_Trunk_xyz(3,:) = (markerset(47-sond_cond,:)+markerset(19-sond_cond,:)+markerset(103-sond_cond,:)+markerset(99-sond_cond,:))/4;
% plot3(CoM_Trunk_xyz(1,:), CoM_Trunk_xyz(2,:), CoM_Trunk_xyz(3,:))
% axis equal
%% CoM UpperArm_r
RSHO = [];
RELB = [];
UpArm_r = [];
RSHO(1,:) = (markerset(45-sond_cond,:));%Top marker
RSHO(2,:) = (markerset(46-sond_cond,:));
RSHO(3,:) = (markerset(47-sond_cond,:));
RELB(1,:) = (markerset(57-sond_cond,:));
RELB(2,:) = (markerset(58-sond_cond,:));
RELB(3,:) = (markerset(59-sond_cond,:));
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
LSHO(1,:) = (markerset(17-sond_cond,:));
LSHO(2,:) = (markerset(18-sond_cond,:));
LSHO(3,:) = (markerset(19-sond_cond,:));
LELB(1,:) = (markerset(29-sond_cond,:));
LELB(2,:) = (markerset(30-sond_cond,:));
LELB(3,:) = (markerset(31-sond_cond,:));
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

LWRI(1,:) = (markerset(33-sond_cond,:)); % wrist outside
LWRI(2,:) = (markerset(34-sond_cond,:));
LWRI(3,:) = (markerset(35-sond_cond,:));

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

RWRI(1,:) = (markerset(61-sond_cond,:));
RWRI(2,:) = (markerset(62-sond_cond,:));
RWRI(3,:) = (markerset(63-sond_cond,:));

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
RASI(1,:) = (markerset(101-sond_cond,:));
RASI(2,:) = (markerset(102-sond_cond,:));
RASI(3,:) = (markerset(103-sond_cond,:));
RKNE(1,:) = (markerset(133-sond_cond,:));
RKNE(2,:) = (markerset(134-sond_cond,:));
RKNE(3,:) = (markerset(135-sond_cond,:));
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

RANK(1,:) = (markerset(141-sond_cond,:));
RANK(2,:) = (markerset(142-sond_cond,:));
RANK(3,:) = (markerset(143-sond_cond,:));

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

LASI(1,:) = (markerset(97-sond_cond,:));
LASI(2,:) = (markerset(98-sond_cond,:));
LASI(3,:) = (markerset(99-sond_cond,:));
LKNE(1,:) = (markerset(109-sond_cond,:));
LKNE(2,:) = (markerset(110-sond_cond,:));
LKNE(3,:) = (markerset(111-sond_cond,:));
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

LANK(1,:) = (markerset(117-sond_cond,:));
LANK(2,:) = (markerset(118-sond_cond,:));
LANK(3,:) = (markerset(119-sond_cond,:));

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
%% CoM_Head_neck
CoM_Head_xyz = [];

if(idxPart==14)
    CoM_Head_xyz(1,:) = (markerset(5-sond_cond,:)+markerset(9-sond_cond,:)+markerset(13-sond_cond,:) + markerset(75-sond_cond,:))/5; % muss ich anpassen
    CoM_Head_xyz(2,:) = (markerset(6-sond_cond,:)+markerset(10-sond_cond,:)+markerset(14-sond_cond,:) + markerset(76-sond_cond,:))/5;
    CoM_Head_xyz(3,:) = (markerset(7-sond_cond,:)+markerset(11-sond_cond,:)+markerset(15-sond_cond,:) + markerset(77-sond_cond,:))/5;
else
    CoM_Head_xyz(1,:) = (markerset(1,:)+markerset(5-sond_cond,:)+markerset(9-sond_cond,:)+markerset(13-sond_cond,:) + markerset(75-sond_cond,:))/5; % muss ich anpassen
    CoM_Head_xyz(2,:) = (markerset(2,:)+markerset(6-sond_cond,:)+markerset(10-sond_cond,:)+markerset(14-sond_cond,:) + markerset(76-sond_cond,:))/5;
    CoM_Head_xyz(3,:) = (markerset(3,:)+markerset(7-sond_cond,:)+markerset(11-sond_cond,:)+markerset(15-sond_cond,:) + markerset(77-sond_cond,:))/5;
end
%% CoM Hand
Hand_l = [];
Hand_l(1,:) = skeleton(1,11,:);
Hand_l(2,:) = skeleton(2,11,:);
Hand_l(3,:) = skeleton(3,11,:);

CoM_Hand_l = Hand_l * CoM_Segment_p(1);
Hand_r = [];
Hand_r(1,:) = skeleton(1,16,:);
Hand_r(2,:) = skeleton(2,16,:);
Hand_r(3,:) = skeleton(3,16,:);

CoM_Hand_r = Hand_r * CoM_Segment_p(1);
%% CoM Foot
Foot_l = [];
Foot_l(1,:) = markerset(117-sond_cond,:)+ markerset(165-sond_cond,:)/2;
Foot_l(2,:) = markerset(118-sond_cond,:)+markerset(166-sond_cond,:)/2;
Foot_l(3,:) = markerset(119-sond_cond,:)+markerset(167-sond_cond,:)/2;

CoM_Foot_l = Foot_l * CoM_Segment_p(4);
Foot_r = [];
Foot_r(1,:) = markerset(141-sond_cond,:)+markerset(153-sond_cond,:)/2;
Foot_r(2,:) = markerset(142-sond_cond,:)+markerset(154-sond_cond,:)/2;
Foot_r(3,:) = markerset(143-sond_cond,:)+markerset(155-sond_cond,:)/2;

CoM_Foot_r = Foot_r * CoM_Segment_p(4);
%% werte anpassen, da Kopf bei mir noch da (was ist Headset ?)
m = table2array(Weight(idxPart-7,1));

CoM_Foot_r(isnan(CoM_Foot_r)) = 0;
CoM_Foot_l(isnan(CoM_Foot_l)) = 0;

CoM_Hand_r(isnan(CoM_Hand_r)) = 0;
CoM_Hand_l(isnan(CoM_Hand_l)) = 0;

CoM_Head_xyz(isnan(CoM_Head_xyz)) = 0;
CoM_Trunk_xyz(isnan(CoM_Trunk_xyz)) = 0;

rel_CoM_UpArm_r(isnan(rel_CoM_UpArm_r)) = 0;
rel_CoM_UpArm_l(isnan(rel_CoM_UpArm_l)) = 0;

rel_CoM_LoArm_r(isnan(rel_CoM_LoArm_r)) = 0;
rel_CoM_LoArm_l(isnan(rel_CoM_LoArm_l)) = 0;

rel_CoM_Shank_r(isnan(rel_CoM_Shank_r)) = 0;
rel_CoM_Shank_l(isnan(rel_CoM_Shank_l)) = 0;

rel_CoM_Thigh_r(isnan(rel_CoM_Thigh_r)) = 0;
rel_CoM_Thigh_l(isnan(rel_CoM_Thigh_l)) = 0;

CoM = (CoM_Foot_r*0.0145*m + CoM_Foot_l * 0.0145*m +CoM_Hand_r*0.006*m + CoM_Hand_l*0.006*m + 0.081*CoM_Head_xyz*m + rel_CoM_UpArm_r * 0.028 *m + rel_CoM_UpArm_l * 0.028 *m + rel_CoM_LoArm_l * 0.016 *m+ rel_CoM_LoArm_r * 0.016*m + rel_CoM_Shank_l * 0.0465 *m+ rel_CoM_Shank_r * 0.0465 *m+ rel_CoM_Thigh_l * 0.1 *m+ rel_CoM_Thigh_r * 0.1 *m+ CoM_Trunk_xyz * 0.497*m)/m;
vel_tmp = diff(CoM(1,:))*Fs/1000;

%% Cutting of longer trials
cut_velocity = 0.19;
min_walk_time = 1;

% Anzahl der zusätzlichen Frames nach dem Endpunkt
extra_frames = 0; % z. B. 50 Frames (entspricht 0.25 Sekunden bei 200 Hz)
foot_middle = [];
if contains(list_of_names(idxTrial, idxPart), "baseline")
    velocity_CoM = sqrt(diff(CoM(1,:)).^2 + diff(CoM(2,:)).^2) / dt;
    velocity_CoM(velocity_CoM> 1300) = 0;
    [c,ca] = butter(4,fc/(Fs/2));
    velocity_CoM_filt = filtfilt(c,ca,velocity_CoM');
else
    foot_middle(1,:) = ((skeleton(1,20,:) + skeleton(1,24,:))/2);
    foot_middle(2,:) = ((skeleton(2,20,:) + skeleton(2,24,:))/2);    
    foot_middle(3,:) = ((skeleton(3,20,:) + skeleton(3,24,:))/2);
    RHeel = [];
    LHeel = [];
    
    LHeel(1,:) = (markerset(121-sond_cond,:));
    LHeel(2,:) = (markerset(122-sond_cond,:));
    LHeel(3,:) = (markerset(123-sond_cond,:));
    RHeel(1,:) = (markerset(145-sond_cond,:));
    RHeel(2,:) = (markerset(146-sond_cond,:));
    RHeel(3,:) = (markerset(147-sond_cond,:));
    % Berechnung der Geschwindigkeit des CoM (in mm/s)
   
    velocity_CoM = sqrt(diff(CoM(1,:)).^2 + diff(CoM(2,:)).^2) / dt;
    velocity_CoM(velocity_CoM> 1300) =0;
    [c,ca] = butter(4,fc/(Fs/2));
    velocity_CoM_filt = filtfilt(c,ca,velocity_CoM');

    velocity_window = 21; % Fenstergröße für Geschwindigkeit (z. B. 11 Frames)
    acceleration_window = 21; % Fenstergröße für Beschleunigung (z. B. 11 Frames)

    % Geschwindigkeit berechnen und glätten
    raw_velocity_x = diff(CoM(1,:)) / dt; % Originale Geschwindigkeit
    smoothed_velocity_x = movmean(raw_velocity_x, velocity_window); % Geglättete Geschwindigkeit

    % Beschleunigung berechnen aus geglätteter Geschwindigkeit
    smoothed_acceleration_x = diff(smoothed_velocity_x) / dt;

    % Gleitender Mittelwert der Beschleunigung
    mean_acceleration_x = movmean(smoothed_acceleration_x, acceleration_window);

    % Berechnung der euklidischen Disatnz zwischen Mittelfuß (CoP) & CoM
    diff_C_FM1 = CoM(1,:)- foot_middle(1,:);
    diff_C_FM2 = CoM(2,:)- foot_middle(2,:);
    euclid = sqrt(diff_C_FM1.^2 + diff_C_FM2.^2);
    
    % Erkennung der stabilen Start- und Endphasen
    stable_walking_start = find(euclid > 20);
    stable_walking_end = find(velocity_CoM_filt(1:end-1) < 190 & diff(velocity_CoM_filt) < 0);
    
    
    % Parameter: Mindestabstand zwischen zwei Startpunkten
    min_frames_between_starts = 2500; % z. B. 500 Frames (entspricht 2.5 Sekunden bei 200 Hz)
    
    % Trials initialisieren
    trials = {};
    foot_angles_l ={};
    foot_angles_r ={};
    knee_angles_l ={};
    knee_angles_r ={};
    foot_RoM_l =[];
    foot_RoM_r =[];
    knee_RoM_l =[];
    knee_RoM_r =[];

    trials_marker = {};
    start_idx = 1;
    last_end_idx = 0; % Index des letzten Endpunktes
    
    % fürs debuging
    timings = [];
    timinge = [];
    
    while ~isempty(stable_walking_start)
        % Finde den nächsten Startpunkt
         trial_start = stable_walking_start(1);
   
        
        if  RHeel(1,trial_start) > 800 || trial_start > length(mean_acceleration_x)
            stable_walking_start(1) = []; % Entferne diesen Startpunkt
            continue;
        end
        if mean_acceleration_x(trial_start) < 0 
            stable_walking_start(1) =[];
            continue;
        end
        if(idxPart == 12 && idxTrial == 9 && RHeel(1, trial_start) > 500 ) 
            stable_walking_start(1) =[];
            continue;
        end% sonderfall
        
        if(idxPart==14 && idxTrial==7)
            min_frames_between_starts=2000;
        else
            min_frames_between_starts=2500;
        end
         % Entferne Startpunkte, die vor dem letzten Endpunkt liegen oder zu
        % nah dran sind und nicht im Bereichder start force plates liegen
        if trial_start < last_end_idx + min_frames_between_starts  && ~isempty(trials)
            stable_walking_start(1) = []; % Entferne diesen Startpunkt
            continue;
        end
        % Entferne alle verbleibenden Startpunkte bis trial_start
        stable_walking_start(stable_walking_start <= trial_start) = [];
        
        % Finde den nächsten Endpunkt nach dem Startpunkt
        trial_end = stable_walking_end(find(stable_walking_end > trial_start+400, 1)); % + 400, weil mindestens 2 sek laufzeit
        
        % Falls kein Endpunkt gefunden wird, abbrechen
        if isempty(trial_end)
            break;
        end
        
        % Sicherstellen, dass der Endpunkt nach dem Startpunkt liegt und auch
        % eine ausreichende Geschwindkeit 
        if trial_end > trial_start && max(velocity_CoM_filt(trial_start:trial_end)) > 620
            % Speichere den aktuellen Trial
  
            % **Endpunkt erweitern**:
            % Stelle sicher, dass wir nicht über die Datenlänge hinausgehen
            if(idxPart==14 &&idxTrial==7 && trial_start==31328)
                trial_start = 31073;
            end
            if(idxPart==13 &&idxTrial==8 )
                if( trial_start==19181)
                    trial_start = 19059;
                end
                if( trial_start==35746)
                    trial_start = 35605;
                end
                if( trial_start==38964)
                    trial_start = 38786;
                end
            end
            if(idxPart==12 &&idxTrial==8 && trial_start==25685)
                trial_start = 25615;
            end
            if(idxPart==10 &&idxTrial==7 && trial_start==33949)
                trial_start = 33805;
             end
            if(idxPart==10 &&idxTrial==8 && trial_start==25470)
                trial_start = 25355;
            end
            if(idxPart==8 &&idxTrial==6 && trial_start==27570)
                trial_start = 27490;
             end
            if(idxPart==8 &&idxTrial==7)
                if( trial_start==9968)
                    trial_start = 9877;
                end
                if( trial_start==32995)
                    trial_start = 32955;
                end
            end
            trial_end_extended = min(trial_end + extra_frames, size(CoM, 2));

       
            trials{start_idx} = CoM(:, trial_start:trial_end_extended);
            trials_marker{start_idx} = markerset(:, trial_start:trial_end_extended);
            timings(start_idx) = trial_start;
            timinge(start_idx) = trial_end;
            

            %Range of Motion
            V0 = CoM_Shank_r(:, trial_start:trial_end_extended);
            V1 = CoM_Shank_l(:, trial_start:trial_end_extended);
            V2 = CoM_Thigh_r(:, trial_start:trial_end_extended);
            V3 = CoM_Thigh_l(:, trial_start:trial_end_extended);
            Knee_angle_r =[];
            Knee_angle_l =[];
            V4 = CoM_Foot_r(:, trial_start:trial_end_extended);
            V5 = CoM_Foot_l(:, trial_start:trial_end_extended);
            Foot_angle_r =[];
            Foot_angle_l =[];

            for idx = 1:length(V0)
                Knee_angle_r(idx) = atan2d(norm(cross(V0(:,idx), V2(:,idx))), dot(V0(:,idx), V2(:,idx)));            
            end
            for idx = 1:length(V1)
                Knee_angle_l(idx) = atan2d(norm(cross(V1(:,idx), V3(:,idx))), dot(V1(:,idx), V3(:,idx)));
            end
            for idx = 1:length(V0)
                Foot_angle_r(idx) = atan2d(norm(cross(V4(:,idx), V0(:,idx))), dot(V4(:,idx), V0(:,idx)));
            end
            for idx = 1:length(V1)
                Foot_angle_l(idx) = atan2d(norm(cross(V5(:,idx), V1(:,idx))), dot(V5(:,idx), V1(:,idx)));
            end

            foot_angles_l{start_idx} = Foot_angle_l;
            foot_angles_r{start_idx} = Foot_angle_r;
            knee_angles_l{start_idx} = Knee_angle_l;
            knee_angles_r{start_idx} = Knee_angle_r;

            start_idx = start_idx + 1;
         % Aktualisiere den letzten Endpunkt
            last_end_idx = trial_end;
        end
        
        % Entferne alle Endpunkte bis einschließlich trial_end
        stable_walking_end(stable_walking_end <= trial_end) = [];
    end
    
    %Anzahl der gefundenen Trials ausgeben
    % disp(['Anzahl der Trials gefunden: ', num2str(length(trials))]);
    % figure
    % scatter(timings, velocity_CoM_filt(timings), "k", "filled")
    % hold on
    % scatter(timinge, velocity_CoM_filt(timinge), "r", "filled")
    % plot(velocity_CoM_filt, "b")
    % plot(RHeel(1,:), "r")
    % ylabel('Geschwindigkeit CoM (mm/s)');
    % title(["Particpant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart)])
    % shg

   
end

%% schritte schneiden
cond = list_of_names(idxTrial, idxPart);
% Mindestlänge in Frames (z. B. 50 Frames bei 200 Hz entspricht 0.25 Sekunden)
min_step_length = 50;
if(contains(cond, "baseline"))
    RHeel = [];
    LHeel = [];
    foot_angles_l =[];
    foot_angles_r =[];
    knee_angles_l =[];
    knee_angles_r =[];
    min_frame_distance = 140; % Mindestens 100 Frames Abstand zwischen Touchdown-Punkten
    LKNE(1,:) = (markerset(109-sond_cond,:));
    LKNE(2,:) = (markerset(110-sond_cond,:));
    LKNE(3,:) = (markerset(111-sond_cond,:));
    RKNE(1,:) = (markerset(133-sond_cond,:));
    RKNE(2,:) = (markerset(134-sond_cond,:));
    RKNE(3,:) = (markerset(135-sond_cond,:));
    LHeel(1,:) = (markerset(121-sond_cond,:));
    LHeel(2,:) = (markerset(122-sond_cond,:));
    LHeel(3,:) = (markerset(123-sond_cond,:));
    RHeel(1,:) = (markerset(145-sond_cond,:));
    RHeel(2,:) = (markerset(146-sond_cond,:));
    RHeel(3,:) = (markerset(147-sond_cond,:));
    % Dynamische Anpassung der MinPeakHeight
    heel_range_r = max(RHeel(3,:)) - min(RHeel(3,:));
    heel_range_l = max(LHeel(3,:)) - min(LHeel(3,:));

    % Faktor für MinPeakHeight (z. B. 20% des Signalbereichs)
    min_peak_factor = 0.66;

    % Dynamisch berechneter MinPeakHeight
    min_peak_height_r = heel_range_r * min_peak_factor;
    min_peak_height_l = heel_range_l * min_peak_factor;

    %% RHeel
    
    signal = RHeel(3,:); 
    [maxima, max_indices] = findpeaks(signal, 'MinPeakHeight', min_peak_height_r);
    
    search_range = 100;

    % Speicher für Schnittpunkte
    step_boundaries_r = []; % Start- und Endpunkte der Schritte

    for ii = 1:length(max_indices)
        max_idx = max_indices(ii);

        % Suche linkes Minimum
        left_range = max(1, max_idx - search_range):max_idx; % Bereich links vom Maximum
        [left_min, left_idx] = min(signal(left_range));
        left_idx = left_idx + left_range(1) - 1; % Absoluter Index

        % Suche rechtes Minimum
        right_range = max_idx:min(length(signal), max_idx + search_range); % Bereich rechts vom Maximum
        [right_min, right_idx] = min(signal(right_range));
        right_idx = right_idx + right_range(1) - 1; % Absoluter Index

        % Speichere Schrittgrenzen
        step_boundaries_r = [step_boundaries_r; left_idx, right_idx];
    end

    % %Minima visualisieren
    % figure;
    % plot(signal);
    % hold on;
    % plot(max_indices, maxima, 'ro'); % Maxima
    % plot(step_boundaries_r(:, 1), signal(step_boundaries_r(:, 1)), 'go'); % Linke Minima
    % plot(step_boundaries_r(:, 2), signal(step_boundaries_r(:, 2)), 'mo'); % Rechte Minima
    % title(["Right Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart), " Trial: ", i])
    % xlabel('Frame');
    % ylabel('Amplitude');
    % hold off;

    % Schritte extrahieren
    steps_r = cell(size(step_boundaries_r, 1), 1);

    for ii = 1:size(step_boundaries_r, 1)
        step_start = step_boundaries_r(ii, 1);
        step_end = step_boundaries_r(ii, 2);
        steps_r{ii} = signal(step_start:step_end);
        
    end
    valid_steps_r = step_boundaries_r((step_boundaries_r(:,2) - step_boundaries_r(:,1)) >= min_step_length, :);
    last_step_incomplete_r = false;
    if ~isempty(valid_steps_r) && valid_steps_r(end,2) >= size(RHeel, 2) - 1
        last_step_incomplete_r = true;
        valid_steps_r(end, :) = [];
    end

    num_steps_r = size(valid_steps_r, 1);
    step_lengths_r = zeros(num_steps_r, 1);

    % Berechnung der Schrittlänge für jeden Schritt rechts
    for ii = 1:num_steps_r
        start_idx = valid_steps_r(ii, 1);
        end_idx = valid_steps_r(ii, 2);   % Rechtes Minimum (Ende des Schritts)
        %Range of Motion
        V0 = CoM_Shank_r(:,start_idx:end_idx);
        V2 = CoM_Thigh_r(:,start_idx:end_idx);
        V4 = CoM_Foot_r(:,start_idx:end_idx);
        Knee_angle_r =[];
        Foot_angle_r =[];
        for idx = 1:length(V0)
            Knee_angle_r(idx) = atan2d(norm(cross(V0(:,idx), V2(:,idx))), dot(V0(:,idx), V2(:,idx)));
        end

        for idx = 1:length(V0)
            Foot_angle_r(idx) = atan2d(norm(cross(V4(:,idx), V0(:,idx))), dot(V4(:,idx), V0(:,idx)));
        end

        foot_angles_r(ii) = max(Foot_angle_r)-min(Foot_angle_r);
        knee_angles_r(ii) = max(Knee_angle_r)-min(Knee_angle_r);
        % Schrittlänge berechnen (Differenz der x-Koordinaten)
        step_lengths_r(ii) = abs(RHeel(1,end_idx) - RHeel(1,start_idx));
    end
    %% Lheel
 
    signal = LHeel(3,:); 
    [maxima, max_indices] = findpeaks(signal, 'MinPeakHeight', min_peak_height_l);
    
    search_range = 100;

    % Speicher für Schnittpunkte
    step_boundaries_l = []; % Start- und Endpunkte der Schritte

    for ii = 1:length(max_indices)
        max_idx = max_indices(ii);

        % Suche linkes Minimum
        left_range = max(1, max_idx - search_range):max_idx; % Bereich links vom Maximum
        [left_min, left_idx] = min(signal(left_range));
        left_idx = left_idx + left_range(1) - 1; % Absoluter Index

        % Suche rechtes Minimum
        right_range = max_idx:min(length(signal), max_idx + search_range); % Bereich rechts vom Maximum
        [right_min, right_idx] = min(signal(right_range));
        right_idx = right_idx + right_range(1) - 1; % Absoluter Index

        % Speichere Schrittgrenzen
        step_boundaries_l = [step_boundaries_l; left_idx, right_idx];
    end

    % Schritte extrahieren
    steps_l = cell(size(step_boundaries_l, 1), 1);

    for ii = 1:size(step_boundaries_l, 1)
        step_start = step_boundaries_l(ii, 1);
        step_end = step_boundaries_l(ii, 2);
        steps_l{ii} = signal(step_start:step_end);
    end
    valid_steps_l = step_boundaries_l((step_boundaries_l(:,2) - step_boundaries_l(:,1)) >= min_step_length, :);
    last_step_incomplete_l = false;
    if ~isempty(valid_steps_l) && valid_steps_l(end,2) >= size(LHeel, 2) - 1
        last_step_incomplete_l = true;
        valid_steps_l(end, :) = [];
    end

    num_steps_l = size(valid_steps_l, 1);
    step_lengths_l = zeros(num_steps_l, 1);

    % Berechnung der Schrittlänge für jeden Schritt rechts
    for ii = 1:num_steps_l
        start_idx = valid_steps_l(ii, 1);
        end_idx = valid_steps_l(ii, 2);   % Rechtes Minimum (Ende des Schritts)
        % Schrittlänge berechnen (Differenz der x-Koordinaten)
        step_lengths_l(ii) = abs(LHeel(1,end_idx) - LHeel(1,start_idx));
        %Range of Motion
        V0 = CoM_Shank_l(:,start_idx:end_idx);
        V2 = CoM_Thigh_l(:,start_idx:end_idx);
        V4 = CoM_Foot_l(:,start_idx:end_idx);
        Knee_angle_l =[];
        Foot_angle_l =[];

        for idx = 1:length(V0)
            Knee_angle_l(idx) = atan2d(norm(cross(V0(:,idx), V2(:,idx))), dot(V0(:,idx), V2(:,idx)));
        end

        for idx = 1:length(V0)
            Foot_angle_l(idx) = atan2d(norm(cross(V0(:,idx), V4(:,idx))), dot(V0(:,idx), V4(:,idx)));
        end
        foot_angles_l(ii) = max(Foot_angle_l)-min(Foot_angle_l);
        knee_angles_l(ii) = max(Knee_angle_l)-min(Knee_angle_l);
    end
    
    walking_distance = max(CoM(1,:)) - min(CoM(1,:));% sum(abs(diff(CoM(1,:))));
    HeelR_velocity = sqrt(diff(RHeel(1,:)).^2 + diff(RHeel(2,:)).^2 + diff(RHeel(3,:)).^2) ./ dt;
    HeelL_velocity = sqrt(diff(LHeel(1,:)).^2 + diff(LHeel(2,:)).^2 + diff(LHeel(3,:)).^2) ./ dt;

    % Speichere die Ergebnisse
    results_l(idxPart).conditions(1).trials(idxTrial).numStrides = num_steps_l;
    results_l(idxPart).conditions(1).trials(idxTrial).strideLengths = step_lengths_l;
    results_l(idxPart).conditions(1).trials(idxTrial).meanStepLength = mean(step_lengths_l);
    results_l(idxPart).conditions(1).trials(idxTrial).velocity = velocity_CoM_filt;
    results_l(idxPart).conditions(1).trials(idxTrial).walkingDistance = walking_distance;
    results_l(idxPart).conditions(1).trials(idxTrial).footAngle = foot_angles_l;
    results_l(idxPart).conditions(1).trials(idxTrial).kneeAngle = knee_angles_l;
    results_l(idxPart).conditions(1).trials(idxTrial).HeelL_velocity = HeelL_velocity;
    results_l(idxPart).conditions(1).trials(idxTrial).HeelL_z = LHeel(3, :);
    results_l(idxPart).conditions(1).trials(idxTrial).KneeL_z = LKNE(3, :);
    results_l(idxPart).conditions(1).trials(idxTrial).lastStepIncomplete = last_step_incomplete_l;

    results_r(idxPart).conditions(1).trials(idxTrial).HeelR_z = RHeel(3, :);
    results_r(idxPart).conditions(1).trials(idxTrial).KneeR_z = RKNE(3, :);
    % Speichere die Ergebnisse
    results_r(idxPart).conditions(1).trials(idxTrial).numStrides = num_steps_r;
    results_r(idxPart).conditions(1).trials(idxTrial).strideLengths = step_lengths_r;
    results_r(idxPart).conditions(1).trials(idxTrial).meanStepLength = mean(step_lengths_r);
    results_r(idxPart).conditions(1).trials(idxTrial).velocity = velocity_CoM_filt;
    results_r(idxPart).conditions(1).trials(idxTrial).walkingDistance = walking_distance;
    results_r(idxPart).conditions(1).trials(idxTrial).footAngle = foot_angles_r;
    results_r(idxPart).conditions(1).trials(idxTrial).kneeAngle = knee_angles_r;
    results_r(idxPart).conditions(1).trials(idxTrial).HeelR_velocity = HeelR_velocity;
    results_r(idxPart).conditions(1).trials(idxTrial).lastStepIncomplete = last_step_incomplete_r;

else %not baseline
    for i=1 :length(trials_marker)
        %%
        marker_set = trials_marker{i};
        RHeel = [];
        LHeel = [];
        LKNEE =[];
        RKNEE =[];
        min_frame_distance = 125; % Mindestens 100 Frames Abstand zwischen Touchdown-Punkten
        min_start_dist = 10; % mindestens abstand zum Start (soll gegen noise am Anfang helfen)
        n = 200; % Anzahl der letzten Frames, die geprüft werden
        m = 40; % Erweiterter Suchbereich bei fehlendem Minimum
        threshold = 5; % Toleranzbereich um das Minimum (z. B. in mm)
        LKNEE(1,:) = (marker_set(109-sond_cond,:));
        LKNEE(2,:) = (marker_set(110-sond_cond,:));
        LKNEE(3,:) = (marker_set(111-sond_cond,:));
        RKNEE(1,:) = (marker_set(133-sond_cond,:));
        RKNEE(2,:) = (marker_set(134-sond_cond,:));
        RKNEE(3,:) = (marker_set(135-sond_cond,:));
        LHeel(1,:) = (marker_set(121-sond_cond,:));
        LHeel(2,:) = (marker_set(122-sond_cond,:));
        LHeel(3,:) = (marker_set(123-sond_cond,:));
        RHeel(1,:) = (marker_set(145-sond_cond,:));
        RHeel(2,:) = (marker_set(146-sond_cond,:));
        RHeel(3,:) = (marker_set(147-sond_cond,:));
        % Dynamische Anpassung der MinPeakHeight
        heel_range_r = max(RHeel(3,:)) - min(RHeel(3,:));
        heel_range_l = max(LHeel(3,:)) - min(LHeel(3,:));
        
        % Faktor für MinPeakHeight (z. B. 20% des Signalbereichs)
        min_peak_factor = 0.66;
        
        % Dynamisch berechneter MinPeakHeight
        min_peak_height_r = heel_range_r * min_peak_factor;
        min_peak_height_l = heel_range_l * min_peak_factor;

        %% RHeel
        signal = RHeel(3,:); 

        [maxima, max_indices] = findpeaks(signal, 'MinPeakHeight', min_peak_height_r);

        % Reichweite für Suche nach Minima (z. B. 100 Frames nach links und rechts)
        search_range = 100;

        % Speicher für Schnittpunkte
        step_boundaries_r = []; % Start- und Endpunkte der Schritte

        for ii = 1:length(max_indices)
            max_idx = max_indices(ii);

            % Suche linkes Minimum
            left_range = max(1, max_idx - search_range):max_idx; % Bereich links vom Maximum
            [left_min, left_idx] = min(signal(left_range));
            left_idx = left_idx + left_range(1) - 1; % Absoluter Index

            % Suche rechtes Minimum
            right_range = max_idx:min(length(signal), max_idx + search_range); % Bereich rechts vom Maximum
            [right_min, right_idx] = min(signal(right_range));
            right_idx = right_idx + right_range(1) - 1; % Absoluter Index
  
            % Speichere Schrittgrenzen
            step_boundaries_r = [step_boundaries_r; left_idx, right_idx];
        end

        % %Minima visualisieren
        % figure;
        % plot(signal);
        % hold on;
        % plot(max_indices, maxima, 'ro'); % Maxima
        % plot(step_boundaries_r(:, 1), signal(step_boundaries_r(:, 1)), 'go'); % Linke Minima
        % plot(step_boundaries_r(:, 2), signal(step_boundaries_r(:, 2)), 'mo'); % Rechte Minima
        % title(["Right Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart), " Trial: ", i])
        % xlabel('Frame');
        % ylabel('Amplitude');
        % hold off;

        % Schritte extrahieren
        steps_r = cell(size(step_boundaries_r, 1), 1);

        for ii = 1:size(step_boundaries_r, 1)
            step_start = step_boundaries_r(ii, 1);
            step_end = step_boundaries_r(ii, 2);
            steps_r{ii} = signal(step_start:step_end);
           
        end

        
        %% Lheel
      
        signal = LHeel(3,:); 

        [maxima, max_indices] = findpeaks(signal, 'MinPeakHeight', min_peak_height_l);

        % Reichweite für Suche nach Minima (z. B. 100 Frames nach links und rechts)
        search_range = 100;

        % Speicher für Schnittpunkte
        step_boundaries_l = []; % Start- und Endpunkte der Schritte

        for ii = 1:length(max_indices)
            max_idx = max_indices(ii);

            % Suche linkes Minimum
            left_range = max(1, max_idx - search_range):max_idx; % Bereich links vom Maximum
            [left_min, left_idx] = min(signal(left_range));
            left_idx = left_idx + left_range(1) - 1; % Absoluter Index

            % Suche rechtes Minimum
            right_range = max_idx:min(length(signal), max_idx + search_range); % Bereich rechts vom Maximum
            [right_min, right_idx] = min(signal(right_range));
            right_idx = right_idx + right_range(1) - 1; % Absoluter Index

            % Speichere Schrittgrenzen
            step_boundaries_l = [step_boundaries_l; left_idx, right_idx];
        end

        % %Minima visualisieren
        % figure;
        % plot(signal);
        % hold on;
        % plot(max_indices, maxima, 'ro'); % Maxima
        % plot(step_boundaries_l(:, 1), signal(step_boundaries_l(:, 1)), 'go'); % Linke Minima
        % plot(step_boundaries_l(:, 2), signal(step_boundaries_l(:, 2)), 'mo'); % Rechte Minima
        % title(["Left Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart), " Trial: ", i])
        % xlabel('Frame');
        % ylabel('Amplitude');
        % hold off;

        % Schritte extrahieren
        steps_l = cell(size(step_boundaries_l, 1), 1);

        for ii = 1:size(step_boundaries_l, 1)
            step_start = step_boundaries_l(ii, 1);
            step_end = step_boundaries_l(ii, 2);
            steps_l{ii} = signal(step_start:step_end);
        end
        

%%
        % Standardwert für unvollständigen letzten Schritt
        last_step_incomplete_r = false;
        last_step_incomplete_l = false;
        valid_steps_r = step_boundaries_r((step_boundaries_r(:,2) - step_boundaries_r(:,1)) >= min_step_length, :);
        valid_steps_l = step_boundaries_l((step_boundaries_l(:,2) - step_boundaries_l(:,1)) >= min_step_length, :);
        % Überprüfung auf unvollständige Schritte rechts
        if ~isempty(valid_steps_r) && valid_steps_r(end,2) >= size(RHeel, 2) - 1
            last_step_incomplete_r = true;
            valid_steps_r(end, :) = [];
        end
        
        % Überprüfung auf unvollständige Schritte links
        if ~isempty(valid_steps_l) && valid_steps_l(end,2) >= size(LHeel, 2) - 1
            last_step_incomplete_l = true;
            valid_steps_l(end, :) = [];
        end
         % Schrittlängen berechnen
        num_steps_l = size(valid_steps_l, 1);
        step_lengths_l = zeros(num_steps_l, 1);
        % Berechnung der Schrittlänge für jeden Schritt
        for ii = 1:num_steps_l
            start_idx = valid_steps_l(ii, 1);
            end_idx = valid_steps_l(ii, 2); % Rechtes Minimum (Ende des Schritts)

            % Schrittlänge berechnen (Differenz der x-Koordinaten)
            step_lengths_l(ii) = abs(LHeel(1,end_idx) - LHeel(1,start_idx));
            foot_RoM_l(ii) = max(foot_angles_l{i}(start_idx:end_idx))-min(foot_angles_l{i}(start_idx:end_idx));
            knee_RoM_l(ii) = max(knee_angles_l{i}(start_idx:end_idx))-min(knee_angles_l{i}(start_idx:end_idx));
        end
        num_steps_r = size(valid_steps_r, 1);
        step_lengths_r = zeros(num_steps_r, 1);
        % Berechnung der Schrittlänge für jeden Schritt rechts
        for ii = 1:num_steps_r
            start_idx = valid_steps_r(ii, 1);
            end_idx = valid_steps_r(ii, 2);  % Rechtes Minimum (Ende des Schritts)

            % Schrittlänge berechnen (Differenz der x-Koordinaten)
            step_lengths_r(ii) = abs(RHeel(1,end_idx) - RHeel(1,start_idx));

             %Range of Motion
            foot_RoM_r(ii) = max(foot_angles_r{i}(start_idx:end_idx))-min(foot_angles_r{i}(start_idx:end_idx));
            knee_RoM_r(ii) = max(knee_angles_r{i}(start_idx:end_idx))-min(knee_angles_r{i}(start_idx:end_idx));
        end

        walking_distance = max(CoM(1,:)) - min(CoM(1,:));%sum(abs(diff(CoM(1,:))));
        velocity_CoM_results = sqrt(diff(trials{1,i}(1,:)).^2 + diff(trials{1,i}(2,:)).^2) / dt;
        velocity_CoM_results(velocity_CoM_results> 1300) =0;
        [d,da] = butter(4,fc/(Fs/2));
        velocity_CoM_results_f = filtfilt(d,da,velocity_CoM_results');
        HeelR_velocity = sqrt(diff(RHeel(1,:)).^2 + diff(RHeel(2,:)).^2 + diff(RHeel(3,:)).^2) ./ dt;
        HeelL_velocity = sqrt(diff(LHeel(1,:)).^2 + diff(LHeel(2,:)).^2 + diff(LHeel(3,:)).^2) ./ dt;
        
        
        %%
        % Speichere die Ergebnisse
        results_l(idxPart).conditions(idxTrial-4).trials(i).numStrides = num_steps_l;
        results_l(idxPart).conditions(idxTrial-4).trials(i).strideLengths = step_lengths_l;
        results_l(idxPart).conditions(idxTrial-4).trials(i).meanStepLength = mean(step_lengths_l);
        results_l(idxPart).conditions(idxTrial-4).trials(i).lastStepIncomplete = last_step_incomplete_l;
        results_l(idxPart).conditions(idxTrial -4).trials(i).velocity = velocity_CoM_results_f;
        results_l(idxPart).conditions(idxTrial -4).trials(i).walkingDistance = walking_distance;
        results_l(idxPart).conditions(idxTrial -4).trials(i).footAngle = foot_RoM_l;
        results_l(idxPart).conditions(idxTrial -4).trials(i).kneeAngle = knee_RoM_l;
        results_l(idxPart).conditions(idxTrial -4).trials(i).HeelL_velocity = HeelL_velocity;
        results_l(idxPart).conditions(idxTrial-4).trials(i).HeelL_z = LHeel(3, :);
        results_l(idxPart).conditions(idxTrial-4).trials(i).KneeL_z = LKNEE(3, :);
        
        % Speichere die Ergebnisse
        results_r(idxPart).conditions(idxTrial-4).trials(i).numStrides = num_steps_r;
        results_r(idxPart).conditions(idxTrial-4).trials(i).strideLengths = step_lengths_r;
        results_r(idxPart).conditions(idxTrial-4).trials(i).meanStepLength = mean(step_lengths_r);
        results_r(idxPart).conditions(idxTrial-4).trials(i).lastStepIncomplete = last_step_incomplete_r;
        results_r(idxPart).conditions(idxTrial -4).trials(i).velocity = velocity_CoM_results_f;
        results_r(idxPart).conditions(idxTrial -4).trials(i).walkingDistance = walking_distance;
        results_r(idxPart).conditions(idxTrial -4).trials(i).footAngle = foot_RoM_r;
        results_r(idxPart).conditions(idxTrial -4).trials(i).kneeAngle = knee_RoM_r;
        results_r(idxPart).conditions(idxTrial -4).trials(i).HeelR_velocity = HeelR_velocity;
        results_r(idxPart).conditions(idxTrial-4).trials(i).HeelR_z = RHeel(3, :);
        results_r(idxPart).conditions(idxTrial-4).trials(i).KneeR_z = RKNEE(3, :);
    end
end
    end
end