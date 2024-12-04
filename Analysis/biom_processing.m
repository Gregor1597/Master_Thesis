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
Fs = 200; 
%%
% Sampling Frequency (Hz)
counter = 1;
dt = 1 / 200; % Zeitintervall in Sekunden (200 Hz)
for idxPart = 8: 17%sizeArray(2)
    
    counter = counter +1;
    for idxTrial = 1:9%sizeArray(1)      
        % if isempty(out{idxTrial, idxPart}) == 1 %|| idxPart == 1 && idxTrial == 13 || idxPart == 2 && idxTrial == 11 || idxPart == 5 && idxTrial == 2 || idxPart == 5 && idxTrial == 7 || idxPart == 5 && idxTrial == 9 || idxPart == 5 && idxTrial == 10 || idxPart == 5 && idxTrial == 13 || idxPart == 5 && idxTrial == 15 || idxPart == 11 && idxTrial == 3 || idxPart == 11 && idxTrial == 15 || idxPart == 22 && idxTrial == 16 || idxPart == 23 && idxTrial == 1 || idxPart == 23 && idxTrial == 3 || idxPart == 23 && idxTrial == 8 || idxPart == 28 && idxTrial == 2 || idxPart == 41 && idxTrial == 2|| idxPart == 41 && idxTrial == 6|| idxPart == 41 && idxTrial == 7 ||idxPart == 2 && idxTrial == 7 || idxPart == 7 && idxTrial == 9 || idxPart == 8 && idxTrial == 4 || idxPart == 11 && idxTrial == 7 || idxPart == 22 && idxTrial == 16  || idxPart == 22 && idxTrial == 17 || idxPart == 34 && idxTrial == 3 || idxPart == 5 && idxTrial == 4
            % continue
        %else
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


    V0 = CoM_Shank_r;
    V1 = CoM_Shank_l;
    V2 = CoM_Thigh_r;
    V3 = CoM_Thigh_l;
    V4 = CoM_Foot_r;
    V5 = CoM_Foot_l;

    for idx = 1:length(V0)
        Knee_angle_r(idx) = acosd(dot(V0(:,idx),V2(:,idx)) / (norm(V0(:,idx)) * norm(V2(:,idx))) );
    end
    for idx = 1:length(V1)
        Knee_angle_l(idx) = acosd(dot(V1(:,idx),V3(:,idx)) / (norm(V1(:,idx)) * norm(V3(:,idx))) );
    end
    for idx = 1:length(V0)
        Foot_angle_r(idx) = acosd(dot(V4(:,idx),V0(:,idx)) / (norm(V4(:,idx)) * norm(V0(:,idx))) );
    end
    for idx = 1:length(V1)
        Foot_angle_l(idx) = acosd(dot(V5(:,idx),V1(:,idx)) / (norm(V5(:,idx)) * norm(V1(:,idx))) );
    end
    foot_angles_l = max(Foot_angle_l)-min(Foot_angle_l);
    foot_angles_r = max(Foot_angle_r)-min(Foot_angle_r);
    knee_angles_l = max(Knee_angle_l)-min(Knee_angle_l);
    knee_angles_r = max(Knee_angle_r)-min(Knee_angle_r);
else
    foot_middle(1,:) = ((skeleton(1,20,:) + skeleton(1,24,:))/2);
    foot_middle(2,:) = ((skeleton(2,20,:) + skeleton(2,24,:))/2);    
    foot_middle(3,:) = ((skeleton(3,20,:) + skeleton(3,24,:))/2);
    RHeel = [];
    LHeel = [];
    if(idxPart == 14)
        LHeel(1,:) = (markerset(117,:));
        LHeel(2,:) = (markerset(118,:));
        LHeel(3,:) = (markerset(119,:));
        RHeel(1,:) = (markerset(141,:));
        RHeel(2,:) = (markerset(142,:));
        RHeel(3,:) = (markerset(143,:));
    else
        LHeel(1,:) = (markerset(121,:));
        LHeel(2,:) = (markerset(122,:));
        LHeel(3,:) = (markerset(123,:));
        RHeel(1,:) = (markerset(145,:));
        
        RHeel(2,:) = (markerset(146,:));
        RHeel(3,:) = (markerset(147,:));
    end
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
                Knee_angle_r(idx) = acosd(dot(V0(:,idx),V2(:,idx)) / (norm(V0(:,idx)) * norm(V2(:,idx))) );
            end
            for idx = 1:length(V1)
                Knee_angle_l(idx) = acosd(dot(V1(:,idx),V3(:,idx)) / (norm(V1(:,idx)) * norm(V3(:,idx))) );
            end
            for idx = 1:length(V0)
                Foot_angle_r(idx) = acosd(dot(V4(:,idx),V0(:,idx)) / (norm(V4(:,idx)) * norm(V0(:,idx))) );
            end
            for idx = 1:length(V1)
                Foot_angle_l(idx) = acosd(dot(V5(:,idx),V1(:,idx)) / (norm(V5(:,idx)) * norm(V1(:,idx))) );
            end

            foot_angles_l{start_idx} = max(Foot_angle_l)-min(Foot_angle_l);
            foot_angles_r{start_idx} = max(Foot_angle_r)-min(Foot_angle_r);
            knee_angles_l{start_idx} = max(Knee_angle_l)-min(Knee_angle_l);
            knee_angles_r{start_idx} = max(Knee_angle_r)-min(Knee_angle_r);

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
if(contains(cond, "baseline"))
    RHeel = [];
    LHeel = [];
    min_frame_distance = 140; % Mindestens 100 Frames Abstand zwischen Touchdown-Punkten
    if(idxPart == 14)
        LHeel(1,:) = (markerset(117,:));
        LHeel(2,:) = (markerset(118,:));
        LHeel(3,:) = (markerset(119,:));
        RHeel(1,:) = (markerset(141,:));
        RHeel(2,:) = (markerset(142,:));
        RHeel(3,:) = (markerset(143,:));
    else
        LHeel(1,:) = (markerset(121,:));
        LHeel(2,:) = (markerset(122,:));
        LHeel(3,:) = (markerset(123,:));
        RHeel(1,:) = (markerset(145,:));
        RHeel(2,:) = (markerset(146,:));
        RHeel(3,:) = (markerset(147,:));
    end
    %% RHeel
    steps =[];
    touchdowns_r = [];
    start_value = RHeel(3,1);
    % Dynamische Anpassung der Prominenz
    signal = RHeel(3,:); % Beispielsignal (Höhendaten des linken Fußes)
    start_value_r = signal(1);
    
    % Dynamische Prominenz basierend auf Signalbereich
    signal_range = max(signal) - min(signal);
    % if( idxPart==14)
    %       prominence_factor = 0.0075; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % else
    prominence_factor = 0.01; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % end
    dynamic_prominence = signal_range * prominence_factor;

    % Finde lokale Minima mit dynamischer Prominenz
    %touchdowns_r = find(islocalmin(signal, "MinProminence", dynamic_prominence));
    if (idxPart == 14)
        touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif(idxPart == 11 && idxTrial == 2 || idxPart == 11 && idxTrial == 3)
        touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",0.5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif(idxPart == 13 && idxTrial == 1)
        touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif(idxPart == 15 && idxTrial == 4)
        touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    else
        touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",dynamic_prominence)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    end
    touchdowns2 = touchdowns_r( touchdowns_r > 0 );
    filtered_tds_r = touchdowns2([true, diff(touchdowns2) > min_frame_distance]); % filter alle Touchdownpunkte, die den mindestabstand voneinander haben
    %steps{1} = RHeel(3,1: filtered_tds_r(1));
    for i = 1:length(filtered_tds_r)-1
        steps{i} = RHeel(3,filtered_tds_r(i):filtered_tds_r(i+1)-1);
    end
    figure
    scatter(filtered_tds_r, RHeel(3,filtered_tds_r), "r", "filled")
    hold on
    plot(RHeel(3,:))
    plot(diff(RHeel(1,:)))
    title([" Right Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart) ])
    disp(list_of_names(idxTrial, idxPart));

    %% Lheel
    steps_l =[];
    touchdowns_l = [];
    start_value_l = LHeel(3,1);
    min_start_dist = 70;
    % Dynamische Anpassung der Prominenz
    signal = LHeel(3,:); % Beispielsignal (Höhendaten des linken Fußes)
    start_value_l = signal(1);
    
    % Dynamische Prominenz basierend auf Signalbereich
    signal_range = max(signal) - min(signal);
    % if( idxPart==14)
    %       prominence_factor = 0.0075; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % else
    prominence_factor = 0.01; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % end
    dynamic_prominence = signal_range * prominence_factor;

    % Finde lokale Minima mit dynamischer Prominenz
    %touchdowns_l = find(islocalmin(signal, "MinProminence", dynamic_prominence));
    if (idxPart == 9)
        touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence", 2.1)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif( idxPart == 14 && idxTrial == 3 || idxPart == 14 && idxTrial == 4  || idxPart == 14 && idxTrial == 5)
        touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence", 3.5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif(idxPart == 17 && idxTrial == 1 )
        touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",2)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif(idxPart == 10 && idxTrial == 3 || idxPart == 12 && idxTrial == 3 )
        touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif(idxPart == 15)
        if(idxTrial==1 || idxTrial == 3)
            touchdowns_l = find(islocalmin(LHeel(3,:))); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
        elseif(idxTrial==4)
            touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.2)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
        elseif(idxTrial==5)
            touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.2));
            helpee = find(islocalmin(LHeel(3,:),"MinProminence",0.1));% finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
            touchdowns_l(end+1)= helpee(end);
        else
            touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.1)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
        end
        min_start_dist = 245;
    elseif(idxPart == 16 && idxTrial == 4  )
        touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.95)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif (idxPart == 12 && idxTrial == 1 || idxPart == 16 && idxTrial == 5)
        touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    elseif(idxPart == 11 && idxTrial == 5 || idxPart == 12 && idxTrial == 4 )
        touchdowns_l = find(islocalmin(LHeel(3,:))); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

    else
        touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",dynamic_prominence)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    end

    touchdowns2_l = touchdowns_l( touchdowns_l > min_start_dist ); %% mindestabstand zum startpunkt
    filtered_tds_l = touchdowns2_l([true, diff(touchdowns2_l) > min_frame_distance]); % filter alle Touchdownpunkte, die den mindestabstand voneinander haben
    %steps_l{1} = LHeel(3,1: filtered_tds_l(1));
    for i = 1:length(filtered_tds_l)-1
        steps_l{i} = LHeel(3,filtered_tds_l(i):filtered_tds_l(i+1)-1);
    end

        figure;
        plot(LHeel(3,:))
        hold on
        plot(diff(LHeel(1,:)))
        scatter(filtered_tds_l, LHeel(3,filtered_tds_l), "r", "filled")
        title(["Left Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart)])
        disp(list_of_names(idxTrial, idxPart));

        % Schrittlängen berechnen
        step_lengths_l = diff(filtered_tds_l); % Differenz zwischen Touchdown-Indizes = Schrittlänge
        num_steps_l = length(step_lengths_l);   % Anzahl der Schritte im aktuellen Trial
        walking_distance = max(CoM(1,:)) - min(CoM(1,:));% sum(abs(diff(CoM(1,:))));
        % Speichere die Ergebnisse
        results_l(idxPart).conditions(1).trials(idxTrial).numStrides = num_steps_l;
        results_l(idxPart).conditions(1).trials(idxTrial).strideLengths = step_lengths_l;
        results_l(idxPart).conditions(1).trials(idxTrial).meanStepLength = mean(step_lengths_l);
        results_l(idxPart).conditions(1).trials(idxTrial).velocity = velocity_CoM_filt;
        results_l(idxPart).conditions(1).trials(idxTrial).walkingDistance = walking_distance;
        results_l(idxPart).conditions(1).trials(idxTrial).footAngle = foot_angles_l;
        results_l(idxPart).conditions(1).trials(idxTrial).kneeAngle = knee_angles_l;
       

        step_lengths_r = diff(filtered_tds_r); % DFifferenz zwischen Touchdown-Indizes = Schrittlänge
        num_steps_r = length(step_lengths_r);   % Anzahl der Schritte im aktuellen Trial
        
        % Speichere die Ergebnisse
        results_r(idxPart).conditions(1).trials(idxTrial).numStrides = num_steps_r;
        results_r(idxPart).conditions(1).trials(idxTrial).strideLengths = step_lengths_r;
        results_r(idxPart).conditions(1).trials(idxTrial).meanStepLength = mean(step_lengths_r);
        results_r(idxPart).conditions(1).trials(idxTrial).velocity = velocity_CoM_filt;
        results_r(idxPart).conditions(1).trials(idxTrial).walkingDistance = walking_distance;
        results_r(idxPart).conditions(1).trials(idxTrial).footAngle = foot_angles_r;
        results_r(idxPart).conditions(1).trials(idxTrial).kneeAngle = knee_angles_r;
        
    
else
    for i=1 :length(trials_marker)
        %%
        marker_set = trials_marker{i};
        RHeel = [];
        LHeel = [];
        min_frame_distance = 100; % Mindestens 100 Frames Abstand zwischen Touchdown-Punkten
        min_start_dist = 10; % mindestens abstand zum Start (soll gegen noise am Anfang helfen)
        n = 200; % Anzahl der letzten Frames, die geprüft werden
        m = 40; % Erweiterter Suchbereich bei fehlendem Minimum
        threshold = 5; % Toleranzbereich um das Minimum (z. B. in mm)
        if(idxPart == 14)
            LHeel(1,:) = (marker_set(117,:));
            LHeel(2,:) = (marker_set(118,:));
            LHeel(3,:) = (marker_set(119,:));
            RHeel(1,:) = (marker_set(141,:));
            RHeel(2,:) = (marker_set(142,:));
            RHeel(3,:) = (marker_set(143,:));
        else
            LHeel(1,:) = (marker_set(121,:));
            LHeel(2,:) = (marker_set(122,:));
            LHeel(3,:) = (marker_set(123,:));
            RHeel(1,:) = (marker_set(145,:));
            RHeel(2,:) = (marker_set(146,:));
            RHeel(3,:) = (marker_set(147,:));
        end
        %% RHeel
        steps_r = {};
        touchdowns_r = [];

        % Dynamische Anpassung der Prominenz
        signal = RHeel(3,:); % Beispielsignal (Höhendaten des linken Fußes)
        start_value_r = signal(1);
        
        % Dynamische Prominenz basierend auf Signalbereich
        signal_range = max(signal) - min(signal);
        if( idxPart==14)
              prominence_factor = 0.075; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
        else
            prominence_factor = 0.1; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
        end
        dynamic_prominence = signal_range * prominence_factor;

        % Finde lokale Minima mit dynamischer Prominenz
        touchdowns_r = find(islocalmin(signal, "MinProminence", dynamic_prominence));

        % Filter Touchdowns mit Mindestabstand zum Start
        touchdowns2_r = touchdowns_r(touchdowns_r > min_start_dist);

        % Filter Touchdowns mit Mindestabstand zwischen den Punkten
        filtered_tds_r = touchdowns2_r([true, diff(touchdowns2_r) > min_frame_distance]);
        
        % Überprüfen des Endpunkts
        
        % Bereich für die letzten n Frames des Signals
        last_n_frames = signal(end-n+1:end);
        
        % Minimum und Toleranzbereich in den letzten n Frames
        min_last_n = min(last_n_frames);
        lower_bound = min_last_n - threshold;
        upper_bound = min_last_n + threshold;
        
        % Finde Kandidatenpunkte in den letzten n Frames innerhalb des Bereichs
        candidates = find(last_n_frames >= lower_bound & last_n_frames <= upper_bound);
        
        if ~isempty(candidates)
            % Wähle den frühesten bzw niedrigsten Punkt (relativ zum Ende des Signals)
            candidate_index = find(last_n_frames == min(last_n_frames(candidates)));
            global_index = length(signal) - n + candidate_index; % Absoluter Index im Signal

            % Prüfe Mindestabstand zum letzten Touchdown
            if global_index - filtered_tds_r(end) > min_frame_distance
                filtered_tds_r = [filtered_tds_r, global_index]; % Hinzufügen des neuen Touchdowns
            else
                % Fallback: Erweiterten Bereich prüfen (m Frames)
                last_m_frames = signal(max(1, end-m+1):end); % Sicherstellen, dass Index nicht negativ wird
                [min_val, min_idx] = min(last_m_frames); % Minimum im erweiterten Bereich
                global_index = max(1, length(signal) - m + min_idx); % Absoluter Index im Signal
                %disp(global_index);
                % Prüfe Mindestabstand und füge den niedrigsten Punkt hinzu
                if global_index - filtered_tds_r(end) > min_frame_distance
                    filtered_tds_r = [filtered_tds_r, global_index];
                end

            end
        else
            %Fallback: Erweiterten Bereich prüfen (m Frames)
            last_m_frames = signal(max(1, end-m+1):end); % Sicherstellen, dass Index nicht negativ wird
            [min_val, min_idx] = min(last_m_frames); % Minimum im erweiterten Bereich
            global_index = max(1, length(signal) - m + min_idx); % Absoluter Index im Signal
            %disp(global_index);
            % Prüfe Mindestabstand und füge den niedrigsten Punkt hinzu
            if global_index - filtered_tds_r(end) > min_frame_distance
                filtered_tds_r = [filtered_tds_r, global_index];
            end
        end
        % Schritte segmentieren
        
       %steps_r{1} = signal(1:filtered_tds_r(1));
        for ii = 1:length(filtered_tds_r)-1
            steps_r{ii} = signal(filtered_tds_r(ii):filtered_tds_r(ii+1)-1);
        end
        % figure
        % scatter(filtered_tds_r, RHeel(3,filtered_tds_r), "r", "filled")
        % hold on
        % plot(RHeel(3,:))
        % title([" Right Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart), " Trial: ", i])
        % disp(list_of_names(idxTrial, idxPart));

        %% Lheel
        steps_l = {};
        touchdowns_l = [];

        % Dynamische Anpassung der Prominenz
        signal = LHeel(3,:); % Beispielsignal (Höhendaten des linken Fußes)
        start_value_l = signal(1);

        % Dynamische Prominenz basierend auf Signalbereich
        signal_range = max(signal) - min(signal);
        if(idxPart==14)
              prominence_factor = 0.075; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
        else
            prominence_factor = 0.1; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
        end
        dynamic_prominence = signal_range * prominence_factor;

        % Finde lokale Minima mit dynamischer Prominenz
        touchdowns_l = find(islocalmin(signal, "MinProminence", dynamic_prominence));

        % Filter Touchdowns mit Mindestabstand zum Start
        touchdowns2_l = touchdowns_l(touchdowns_l > min_start_dist);

        % Filter Touchdowns mit Mindestabstand zwischen den Punkten
        filtered_tds_l = touchdowns2_l([true, diff(touchdowns2_l) > min_frame_distance]);
        
        % Überprüfen des Endpunkts
        % Bereich für die letzten n Frames des Signals
        last_n_frames = signal(end-n+1:end);
        
        % Minimum und Toleranzbereich in den letzten n Frames
        min_last_n = min(last_n_frames);
        lower_bound = min_last_n - threshold;
        upper_bound = min_last_n + threshold;
        
        % Finde Kandidatenpunkte in den letzten n Frames innerhalb des Bereichs
        candidates = find(last_n_frames >= lower_bound & last_n_frames <= upper_bound);
        
        if ~isempty(candidates)
            % Wähle den frühesten bzw niedrigsten Punkt (relativ zum Ende des Signals)
            candidate_index = find(last_n_frames == min(last_n_frames(candidates)));
            global_index = length(signal) - n + candidate_index; % Absoluter Index im Signal

            % Prüfe Mindestabstand zum letzten Touchdown
            if global_index - filtered_tds_l(end) > min_frame_distance
                filtered_tds_l = [filtered_tds_l, global_index]; % Hinzufügen des neuen Touchdowns
            else
                % Fallback: Erweiterten Bereich prüfen (m Frames)
                last_m_frames = signal(max(1, end-m+1):end); % Sicherstellen, dass Index nicht negativ wird
                [min_val, min_idx] = min(last_m_frames); % Minimum im erweiterten Bereich
                global_index = max(1, length(signal) - m + min_idx); % Absoluter Index im Signal
                %disp(global_index);
                % Prüfe Mindestabstand und füge den niedrigsten Punkt hinzu
                if global_index - filtered_tds_l(end) > min_frame_distance
                    filtered_tds_l = [filtered_tds_l, global_index];
                end
       
            end
        else
            % Fallback: Erweiterten Bereich prüfen (m Frames)
            last_m_frames = signal(max(1, end-m+1):end); % Sicherstellen, dass Index nicht negativ wird
            [min_val, min_idx] = min(last_m_frames); % Minimum im erweiterten Bereich
            global_index = max(1, length(signal) - m + min_idx); % Absoluter Index im Signal
            %disp(global_index);
            % Prüfe Mindestabstand und füge den niedrigsten Punkt hinzu
            if global_index - filtered_tds_l(end) > min_frame_distance
                filtered_tds_l = [filtered_tds_l, global_index];
            end
        end
        
        % Schritte segmentieren
        %steps_l{1} = signal(1:filtered_tds_l(1));
        for j = 1:length(filtered_tds_l)-1
            steps_l{j} = signal(filtered_tds_l(j):filtered_tds_l(j+1)-1);
        end
        % figure;
        % plot(LHeel(3,:))
        % hold on
        % scatter(filtered_tds_l, LHeel(3,filtered_tds_l), "r", "filled")
        % title(["Left Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart), " Trial: ", i])
        % disp(list_of_names(idxTrial, idxPart));

        if(idxPart==17 && idxTrial == 6 && i == 7)
            steps_r(1) = [];
        end

        if(idxPart==17 && idxTrial == 6 && i == 10)
            steps_l(1) = [];
        end


         % Schrittlängen berechnen
        step_lengths_l = diff(filtered_tds_l); % Differenz zwischen Touchdown-Indizes = Schrittlänge
        num_steps_l = length(step_lengths_l);   % Anzahl der Schritte im aktuellen Trial
        walking_distance = max(CoM(1,:)) - min(CoM(1,:));%sum(abs(diff(CoM(1,:))));
        velocity_CoM_results = sqrt(diff(trials{1,i}(1,:)).^2 + diff(trials{1,i}(2,:)).^2) / dt;
        velocity_CoM_results(velocity_CoM_results> 1300) =0;
        [d,da] = butter(4,fc/(Fs/2));
        velocity_CoM_results_f = filtfilt(d,da,velocity_CoM_results');
        
        %disp(mean(velocity_CoM_results));
        %%
        % Speichere die Ergebnisse
        results_l(idxPart).conditions(idxTrial-4).trials(i).numStrides = num_steps_l;
        results_l(idxPart).conditions(idxTrial-4).trials(i).strideLengths = step_lengths_l;
        results_l(idxPart).conditions(idxTrial-4).trials(i).meanStepLength = mean(step_lengths_l);
        results_l(idxPart).conditions(idxTrial -4).trials(i).velocity = velocity_CoM_results_f;
        results_l(idxPart).conditions(idxTrial -4).trials(i).walkingDistance = walking_distance;
        results_l(idxPart).conditions(idxTrial -4).trials(i).footAngle = foot_angles_l{1,i};
        results_l(idxPart).conditions(idxTrial -4).trials(i).kneeAngle = knee_angles_l{1,i};
        
        step_lengths_r = diff(filtered_tds_r); % Differenz zwischen Touchdown-Indizes = Schrittlänge
        num_steps_r = length(step_lengths_r);   % Anzahl der Schritte im aktuellen Trial
        
        % Speichere die Ergebnisse
        results_r(idxPart).conditions(idxTrial-4).trials(i).numStrides = num_steps_r;
        results_r(idxPart).conditions(idxTrial-4).trials(i).strideLengths = step_lengths_r;
        results_r(idxPart).conditions(idxTrial-4).trials(i).meanStepLength = mean(step_lengths_r);
        results_r(idxPart).conditions(idxTrial -4).trials(i).velocity = velocity_CoM_results_f;
        results_r(idxPart).conditions(idxTrial -4).trials(i).walkingDistance = walking_distance;
        results_r(idxPart).conditions(idxTrial -4).trials(i).footAngle = foot_angles_r{1,i};
        results_r(idxPart).conditions(idxTrial -4).trials(i).kneeAngle = knee_angles_r{1,i};
            
    end
end
    end
end

%%
% Combine results_l and results_r into a unified structure
combined_results = struct();
num_trials = [5,10,10,10,10];
conditions = ["Baseline", "Large", "NoAvatar","Normal", "Small" ];
for participant = 8:17
    participant_name = sprintf('Participant%d', participant);
    combined_results.(participant_name) = struct();
    
    for condition = 1:5
        condition_name = sprintf('Condition%d', condition);
        combined_results.(participant_name).(condition_name) = struct();
        if participant == 16 && condition == 5
            num_trials = [5,10,10,10,9];
        elseif participant == 14 && condition ==3
            num_trials = [5,10,9,10,10];
        else
            num_trials = [5,10,10,10,10];
        end
        
        for trial = 1:num_trials(condition)
            % Left foot data
            trial_data_l = results_l(participant).conditions(condition).trials(trial);
            
            % Right foot data
            trial_data_r = results_r(participant).conditions(condition).trials(trial);
            
            % Combine data
            trial_data_combined = struct();
            trial_data_combined.numStrides_l = trial_data_l.numStrides;
            trial_data_combined.strideLengths_l = trial_data_l.strideLengths;
            trial_data_combined.walking_distance_l = trial_data_l.walkingDistance;
            trial_data_combined.footAngle_l = trial_data_l.footAngle;
            trial_data_combined.kneeAngle_l = trial_data_l.kneeAngle;
            trial_data_combined.mean_velocity = mean(trial_data_r.velocity);
            trial_data_combined.numStrides_r = trial_data_r.numStrides;
            trial_data_combined.strideLengths_r = trial_data_r.strideLengths;
            trial_data_combined.walking_distance_r = trial_data_r.walkingDistance;
            trial_data_combined.footAngle_r = trial_data_r.footAngle;
            trial_data_combined.kneeAngle_r = trial_data_r.kneeAngle;

            trial_data_combined.velocity = trial_data_r.velocity;
            
            % Add to combined results
            trial_name = sprintf('Trial%d', trial);
            combined_results.(participant_name).(condition_name).(trial_name) = trial_data_combined;
        end
    end
end

% Save as JSON
jsonText = jsonencode(combined_results, 'PrettyPrint', true);
fileID = fopen('combined_step_analysis_results_with_arrays.json', 'w');
fwrite(fileID, jsonText);
fclose(fileID);

disp('Export with arrays as JSON completed.');

%%
clf
plot(out_fd_off{33, 8}(:,4))
hold on
plot(out_fd_off{33, 8}(:,9))
plot(out_fd_off{33, 8}(:,14))
plot(out_fd_off{33, 8}(:,19))
shg

