%%
diff_C_FM1 = CoM(1,:)- foot_middle(1,:);
diff_C_FM2 = CoM(2,:)- foot_middle(2,:);
euclid = sqrt(diff_C_FM1.^2 + diff_C_FM2.^2);
first = true;
trials = [];
if (first)
    walking = euclid > 15;
    start = find(walking);
    disp(start(1,1));
   % e = CoM(1,start(1,1):end);
    %ends = CoM(1,:) < CoM(1, start(1,1));
    ends = abs(diff(CoM(1,:))) < 0.00095;
    e = find(ends);
    disp(e(start(1,1)+1));
    ttt = CoM(1,start(1,1):e(start(1,1)+1));
    trials{1} = ttt;
    first = false;
%else?
    cut_e = euclid(1, e(start(1,1)+1):end) > 15;
    cut_c = CoM(:, e(start(1,1)+1):end);
    start2 = find(cut_e);
    disp(start2(1,1));
    %ends2 = cut_c(1,:) < cut_c(1, start2(1,1));
    ends2 = abs(diff(cut_c(1,:))) < 0.00095;
    e2 = find(ends2);
    disp(e2(start2(1,1)+1));
    ttt2 = cut_c(1, start2(1,1):e2(start2(1,1)+1));
    trials{2} = ttt2;
    
    cut_e3 = euclid(1, e2(start2(1,1)+1):end) > 15;
    cut_c3 = CoM(:, e2(start2(1,1)+1):end);
    start3 = find(cut_e3);
    disp(start3(1,1));
    %ends2 = cut_c(1,:) < cut_c(1, start2(1,1));
    ends3 = abs(diff(cut_c3(1,:))) < 0.00095;
    e3 = find(ends3);
    disp(e3(start3(1,1)+1));
    ttt3 = cut_c3(1, start3(1,1):e3(start3(1,1)+1));
    trials{3} = ttt3;
end

%%
%% test trails schneiden
ms = CoM(1,:);
Df = diff(ms);
plot(ms)
hold on
plot(Df*200)
%%
plot(CoM(1,:),CoM(2,:))
hold on 
plot(foot_middle(1,:), foot_middle(2,:))

%%
%%
% Berechnung der Geschwindigkeit des CoM (in mm/s)
dt = 1 / 200; % Zeitintervall in Sekunden (200 Hz)
velocity_CoM = sqrt(diff(CoM(1,:)).^2 + diff(CoM(2,:)).^2) / dt;

% Berechnung der euklidischen Disatnz zwischen Mittelfuß (CoP) & CoM
diff_C_FM1 = CoM(1,:)- foot_middle(1,:);
diff_C_FM2 = CoM(2,:)- foot_middle(2,:);
euclid = sqrt(diff_C_FM1.^2 + diff_C_FM2.^2);

% Erkennung der stabilen Start- und Endphasen
stable_walking_start = find(euclid > 15);
stable_walking_end = find(velocity_CoM(1:end-1) < 190 & diff(velocity_CoM) < 0);


% Parameter: Mindestabstand zwischen zwei Startpunkten
min_frames_between_starts = 0; % z. B. 500 Frames (entspricht 2.5 Sekunden bei 200 Hz)

% Trials initialisieren
trials = {};
start_idx = 1;
last_end_idx = 0; % Index des letzten Endpunktes
odd_counter = 1; % uj zurücklaufen herauszufiltern

% fürs debuging
timings = [];
timinge = [];

while ~isempty(stable_walking_start)
    % Finde den nächsten Startpunkt
    trial_start = stable_walking_start(1);
    
    % Entferne Startpunkte, die vor dem letzten Endpunkt liegen oder zu nah dran sind
    if trial_start < last_end_idx + min_frames_between_starts
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
    if trial_end > trial_start && max(velocity_CoM(trial_start:trial_end)) > 700
        % Speichere den aktuellen Trial
        if mod(odd_counter,2) == 1
            trials{start_idx} = CoM(:, trial_start:trial_end);
            timings(start_idx) = trial_start;
            timinge(start_idx) = trial_end;
            start_idx = start_idx + 1;
        end
         % Aktualisiere den letzten Endpunkt
        last_end_idx = trial_end;
        odd_counter = odd_counter +1;
    end
    
    % Entferne alle Endpunkte bis einschließlich trial_end
    stable_walking_end(stable_walking_end <= trial_end) = [];
end

% Anzahl der gefundenen Trials ausgeben
disp(['Anzahl der Trials gefunden: ', num2str(length(trials))]);


%%
% Debugging: Plot der Start- und Endpunkte
figure;
subplot(2,1,1);
plot(euclid);
hold on;
scatter(stable_walking_start, euclid(stable_walking_start), 'g', 'filled');
scatter(stable_walking_end, euclid(stable_walking_end), 'r', 'filled');
ylabel('Distanz CoM-CoP (mm)');
title('Ganginitiation und Gangende');

subplot(2,1,2);
plot(velocity_CoM);
hold on;
scatter(stable_walking_start, velocity_CoM(stable_walking_start), 'g', 'filled');
scatter(stable_walking_end, velocity_CoM(stable_walking_end), 'r', 'filled');
ylabel('Geschwindigkeit CoM (mm/s)');
title('Geschwindigkeit');
%%
scatter(timings, velocity_CoM(timings), "k", "filled")
hold on
scatter(timinge, velocity_CoM(timinge), "r", "filled")
plot(velocity_CoM, "b")
ylabel('Geschwindigkeit CoM (mm/s)');
shg


%%
RHeel = [];
LHeel = [];
min_frame_distance = 150; % Mindestens 100 Frames Abstand zwischen Touchdown-Punkten
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
touchdowns = [];
start_value = RHeel(3,1);
if (idxPart == 14)
    touchdowns = find(islocalmin(RHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

elseif(idxPart == 11 && idxTrial == 2 || idxPart == 11 && idxTrial == 3)
     touchdowns = find(islocalmin(RHeel(3,:),"MinProminence",0.5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

elseif(idxPart == 13 && idxTrial == 1)
     touchdowns = find(islocalmin(RHeel(3,:),"MinProminence",5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

elseif(idxPart == 15 && idxTrial == 4)
    touchdowns = find(islocalmin(RHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz

else
    touchdowns = find(islocalmin(RHeel(3,:),"MinProminence",1)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
end
touchdowns2 = touchdowns( touchdowns > 100 );
filtered_tds = touchdowns2([true, diff(touchdowns2) > min_frame_distance]); % filter alle Touchdownpunkte, die den mindestabstand voneinander haben
steps{1} = RHeel(3,1: filtered_tds(1));
for i = 2:length(filtered_tds)-1
    steps{i} = RHeel(3,filtered_tds(i):filtered_tds(i+1)-1);
end
figure
scatter(filtered_tds, RHeel(3,filtered_tds), "r", "filled")
hold on
plot(RHeel(3,:))
title([" Right Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart)])
disp(list_of_names(idxTrial, idxPart));
%%
 if idxPart == 8 && idxTrial == 9
        deleter = 520 : 1180;
        foot_middle(:, deleter) = [];
    end          
    if idxPart == 16 && idxTrial == 9
        foot_middle = foot_middle(:,7000:end);
    end
    if idxPart == 14 && idxTrial == 9
        %deleter = [920 : 1080, 4720:4880, 8295:8468 ,10670:10840, 13695:13840, 16930:17090, 19855:20015, 20550:20750,22820:22970, 26620:26837];
        %foot_middle(:, deleter) = [];
    end
    if idxPart == 10 && idxTrial == 6
        deleter = 26700 : 27150;
        foot_middle(:, deleter) = [];
    end
    if idxPart == 9 && idxTrial == 7
        % deleter = [10000 : 10612, 25637:25969];
        % foot_middle(:, deleter) = [];
    end
    if idxPart == 14 && idxTrial == 9
        deleter = [1655 : 2245,6047 : 6683, 9571 : 10285, 13356 : 14044, 17084: 17721, 20951 : 21539, 24473: 25024, 28021: 28605, 31700: 32260];
        foot_middle(:,deleter) = [];
    end
    if idxPart == 8 && idxTrial == 8
        deleter = 2300 : 2700;
        foot_middle(:, deleter) = [];
    end         
%%
if idxPart == 8 && idxTrial == 9
                    % deleter = 520 : 1180;
                    % if(first_iter)
                    %     skeleton(:,:,deleter) = [];
                    %     first_iter=false;
                    % end
                    % marker(:, deleter) = [];
                end           
                if idxPart == 8 && idxTrial == 8
                    % deleter = 2300 : 2700;
                    % if(first_iter)
                    %     skeleton(:,:,deleter) = [];
                    %     first_iter=false;
                    % end
                    % marker(:, deleter) = [];
                end 
                % if idxPart == 15 && idxTrial == 6
                %     deleter = [18268 : 18970];
                %     if(first_iter)
                %         skeleton(:,:,deleter) = [];
                %         first_iter=false;
                %     end
                %     marker(:, deleter) = [];
                % end    
                % if idxPart == 15 && idxTrial == 7
                %     deleter = 5850 : 7000;
                %     if(first_iter)
                %         skeleton(:,:,deleter) = [];
                %         first_iter=false;
                %     end
                %     marker(:, deleter) = [];
                % end     
                 if idxPart == 10 && idxTrial == 6
                    % deleter = 26700 : 27150;
                    % if(first_iter)
                    %     skeleton(:,:,deleter) = [];
                    %     first_iter=false;
                    % end
                    % marker(:, deleter) = [];
                end
                % if idxPart == 9 && idxTrial == 7
                %     deleter = [1200:1800, 9450:11090, 13325:14535,17174:17991, 20970:22055, 24797:25945, 29003:30103, 32768:33928];
                %     if(first_iter)
                %         skeleton(:,:,deleter) = [];
                %         first_iter=false;
                %     end
                %     marker(:, deleter) = [];
                % end
                if idxPart == 9 && idxTrial == 8
                    % deleter = [1875:3375, 11303:11390, 14915: 15700,20127:20990, 29645:30178,33046:34095, 37581:38444,39750:44123,49550:49990];
                    % if(first_iter)
                    %     skeleton(:,:,deleter) = [];
                    %     first_iter=false;
                    % end
                    % marker(:, deleter) = [];
                end
                if idxPart == 9 && idxTrial == 6
                    % deleter = [11750:12010,36440:37005];
                    % if(first_iter)
                    %     skeleton(:,:,deleter) = [];
                    %     first_iter=false;
                    % end
                    % marker(:, deleter) = [];
                end
                if idxPart == 14 && idxTrial == 7
                    % deleter = [1655 : 2245,6047 : 6683, 9571 : 10285, 13356 : 14044, 17084: 17721, 20951 : 21539, 24473: 25024, 28021: 28605, 31700: 32260];
                    % delter = deleter - 639;
                    if(first_iter)
                    %     skeleton(:,:,deleter) = [];
                        skeleton = skeleton(:,:,400:end);
                        first_iter=false;
                    end
                    marker = marker(:,400:end);
                    % marker(:,deleter) = [];
                end


                %%
num_trials = [5,9,9,9,9];
save('step_analysis_results_l.mat', 'results_l');

% Optionale Ausgabe als CSV
results_table_l = [];
for participant = 8:17
    for condition = 1:5
        for trial = 1:num_trials(condition)
            trial_data = results_l(participant).conditions(condition).trials(trial);
            results_table_l = [results_table_l; participant, condition, trial, trial_data.numStrides, trial_data.meanStepLength, trial_data.velocity_mean, mean(trial_data.footAngle), mean(trial_data.kneeAngle)];
        end
    end
end

% Speichern als CSV
writematrix(results_table_l,'step_analysis_results_l.csv');

save('step_analysis_results_r.mat', 'results_r');

% Optionale Ausgabe als CSV
results_table_r = [];
for participant = 8:17
    for condition = 1:5
        for trial = 1:num_trials(condition)
            trial_data = results_r(participant).conditions(condition).trials(trial);
            results_table_r = [results_table_r; participant, condition, trial, trial_data.numStrides, trial_data.meanStepLength, trial_data.velocity_mean, mean(trial_data.footAngle), mean(trial_data.kneeAngle)];
        end
    end
end

% Speichern als CSV
writematrix(results_table_r,'step_analysis_results_r.csv');
%%
%  start_value_r = signal(1);
       % 
       %  % Dynamische Prominenz basierend auf Signalbereich
       %  signal_range = max(signal) - min(signal);
       %  % if( idxPart==14)
       %  %       prominence_factor = 0.0075; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
       %  % else
       %  %     prominence_factor = 0.01; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
       %  % end
       %  dynamic_prominence = signal_range * prominence_factor;
       % 
       %  % Finde lokale Minima mit dynamischer Prominenz
       %  touchdowns_r = find(islocalmin(signal, "MinProminence", dynamic_prominence));
       % 
       %  % Filter Touchdowns mit Mindestabstand zum Start
       %  touchdowns2_r = touchdowns_r(touchdowns_r > min_start_dist);
       % 
       %  % Filter Touchdowns mit Mindestabstand zwischen den Punkten
       %  filtered_tds_r = touchdowns2_r([true, diff(touchdowns2_r) > min_frame_distance]);
       % 
       %  % Überprüfen des Endpunkts
       % 
       %  % Bereich für die letzten n Frames des Signals
       %  last_n_frames = signal(end-n+1:end);
       % 
       %  % Minimum und Toleranzbereich in den letzten n Frames
       %  min_last_n = min(last_n_frames);
       %  lower_bound = min_last_n - threshold;
       %  upper_bound = min_last_n + threshold;
       % 
       %  % Finde Kandidatenpunkte in den letzten n Frames innerhalb des Bereichs
       %  candidates = find(last_n_frames >= lower_bound & last_n_frames <= upper_bound);
       % 
       %  if ~isempty(candidates)
       %      % Wähle den frühesten bzw niedrigsten Punkt (relativ zum Ende des Signals)
       %      candidate_index = find(last_n_frames == min(last_n_frames(candidates)));
       %      global_index = length(signal) - n + candidate_index; % Absoluter Index im Signal
       % 
       %      % Prüfe Mindestabstand zum letzten Touchdown
       %      if global_index - filtered_tds_r(end) > min_frame_distance
       %          filtered_tds_r = [filtered_tds_r, global_index]; % Hinzufügen des neuen Touchdowns
       %      else
       %          % Fallback: Erweiterten Bereich prüfen (m Frames)
       %          last_m_frames = signal(max(1, end-m+1):end); % Sicherstellen, dass Index nicht negativ wird
       %          [min_val, min_idx] = min(last_m_frames); % Minimum im erweiterten Bereich
       %          global_index = max(1, length(signal) - m + min_idx); % Absoluter Index im Signal
       %          %disp(global_index);
       %          % Prüfe Mindestabstand und füge den niedrigsten Punkt hinzu
       %          if global_index - filtered_tds_r(end) > min_frame_distance
       %              filtered_tds_r = [filtered_tds_r, global_index];
       %          end
       % 
       %      end
       %  else
       %      %Fallback: Erweiterten Bereich prüfen (m Frames)
       %      last_m_frames = signal(max(1, end-m+1):end); % Sicherstellen, dass Index nicht negativ wird
       %      [min_val, min_idx] = min(last_m_frames); % Minimum im erweiterten Bereich
       %      global_index = max(1, length(signal) - m + min_idx); % Absoluter Index im Signal
       %      %disp(global_index);
       %      % Prüfe Mindestabstand und füge den niedrigsten Punkt hinzu
       %      if global_index - filtered_tds_r(end) > min_frame_distance
       %          filtered_tds_r = [filtered_tds_r, global_index];
       %      end
       %  end
       %  % Schritte segmentieren
       % 
       % %steps_r{1} = signal(1:filtered_tds_r(1));
       %  for ii = 1:length(filtered_tds_r)-1
       %      steps_r{ii} = signal(filtered_tds_r(ii):filtered_tds_r(ii+1)-1);
       %  end

        % figure
        % scatter(filtered_tds_r, RHeel(3,filtered_tds_r), "r", "filled")
        % hold on
        % plot(RHeel(3,:))
        % title([" Right Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart), " Trial: ", i])
        % disp(list_of_names(idxTrial, idxPart));
        %%
        % start_value_r = signal(1);
    % 
    % % Dynamische Prominenz basierend auf Signalbereich
    % signal_range = max(signal) - min(signal);
    % % if( idxPart==14)
    % %       prominence_factor = 0.0075; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % % else
    % prominence_factor = 0.01; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % % end
    % dynamic_prominence = signal_range * prominence_factor;
    % 
    % % Finde lokale Minima mit dynamischer Prominenz
    % %touchdowns_r = find(islocalmin(signal, "MinProminence", dynamic_prominence));
    % if (idxPart == 14)
    %     touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    % 
    % elseif(idxPart == 11 && idxTrial == 2 || idxPart == 11 && idxTrial == 3)
    %     touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",0.5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    % 
    % elseif(idxPart == 13 && idxTrial == 1)
    %     touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    % 
    % elseif(idxPart == 15 && idxTrial == 4)
    %     touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    % 
    % else
    %     touchdowns_r = find(islocalmin(RHeel(3,:),"MinProminence",dynamic_prominence)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    % end
    % touchdowns2 = touchdowns_r( touchdowns_r > 0 );
    % filtered_tds_r = touchdowns2([true, diff(touchdowns2) > min_frame_distance]); % filter alle Touchdownpunkte, die den mindestabstand voneinander haben
    % %steps{1} = RHeel(3,1: filtered_tds_r(1));
    % for i = 1:length(filtered_tds_r)-1
    %     steps{i} = RHeel(3,filtered_tds_r(i):filtered_tds_r(i+1)-1);
    % end
    % figure
    % scatter(filtered_tds_r, RHeel(3,filtered_tds_r), "r", "filled")
    % hold on
    % plot(RHeel(3,:))
    % 
    % title([" Right Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart) ])
    % disp(list_of_names(idxTrial, idxPart));
    %%
    % start_value_l = signal(1);
    %
    % % Dynamische Prominenz basierend auf Signalbereich
    % signal_range = max(signal) - min(signal);
    % % if( idxPart==14)
    % %       prominence_factor = 0.0075; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % % else
    % prominence_factor = 0.01; % Faktor zur Anpassung (z. B. 10 % des Signalbereichs)
    % % end
    % dynamic_prominence = signal_range * prominence_factor;
    %
    % % Finde lokale Minima mit dynamischer Prominenz
    % %touchdowns_l = find(islocalmin(signal, "MinProminence", dynamic_prominence));
    % if (idxPart == 9)
    %     touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence", 2.1)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %
    % elseif( idxPart == 14 && idxTrial == 3 || idxPart == 14 && idxTrial == 4  || idxPart == 14 && idxTrial == 5)
    %     touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence", 3.5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %
    % elseif(idxPart == 17 && idxTrial == 1 )
    %     touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",2)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %
    % elseif(idxPart == 10 && idxTrial == 3 || idxPart == 12 && idxTrial == 3 )
    %     touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.5)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %
    % elseif(idxPart == 15)
    %     if(idxTrial==1 || idxTrial == 3)
    %         touchdowns_l = find(islocalmin(LHeel(3,:))); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %     elseif(idxTrial==4)
    %         touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.2)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %     elseif(idxTrial==5)
    %         touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.2));
    %         helpee = find(islocalmin(LHeel(3,:),"MinProminence",0.1));% finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %         touchdowns_l(end+1)= helpee(end);
    %     else
    %         touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.1)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %     end
    %     min_start_dist = 245;
    % elseif(idxPart == 16 && idxTrial == 4  )
    %     touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.95)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %
    % elseif (idxPart == 12 && idxTrial == 1 || idxPart == 16 && idxTrial == 5)
    %     touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",0.3)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %
    % elseif(idxPart == 11 && idxTrial == 5 || idxPart == 12 && idxTrial == 4 )
    %     touchdowns_l = find(islocalmin(LHeel(3,:))); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    %
    % else
    %     touchdowns_l = find(islocalmin(LHeel(3,:),"MinProminence",dynamic_prominence)); % finde alle lokalen minima mit notfalls indviduell angepasster Prominenz
    % end
    %
    % touchdowns2_l = touchdowns_l( touchdowns_l > min_start_dist ); %% mindestabstand zum startpunkt
    % filtered_tds_l = touchdowns2_l([true, diff(touchdowns2_l) > min_frame_distance]); % filter alle Touchdownpunkte, die den mindestabstand voneinander haben
    % %steps_l{1} = LHeel(3,1: filtered_tds_l(1));
    % for i = 1:length(filtered_tds_l)-1
    %     steps_l{i} = LHeel(3,filtered_tds_l(i):filtered_tds_l(i+1)-1);
    % end
    %
    %     figure;
    %     plot(LHeel(3,:))
    %     hold on
    %     scatter(filtered_tds_l, LHeel(3,filtered_tds_l), "r", "filled")
    %     title(["Left Heel: Participant: ",idxPart , " ;Condition: ", list_of_names(idxTrial, idxPart)])
    %     disp(list_of_names(idxTrial, idxPart));