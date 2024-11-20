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

