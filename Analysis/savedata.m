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
            trial_data_combined.meanstrideLength_l = trial_data_l.meanStepLength;
            trial_data_combined.last_step_incomplete_l = trial_data_l.lastStepIncomplete;
            trial_data_combined.walking_distance_l = trial_data_l.walkingDistance;
            trial_data_combined.footAngle_l = trial_data_l.footAngle;
            trial_data_combined.kneeAngle_l = trial_data_l.kneeAngle;
            trial_data_combined.mean_velocity = mean(trial_data_r.velocity);
            trial_data_combined.numStrides_r = trial_data_r.numStrides;
            trial_data_combined.strideLengths_r = trial_data_r.strideLengths;
            trial_data_combined.meanstrideLength_r = trial_data_r.meanStepLength;
            trial_data_combined.last_step_incomplete_r = trial_data_r.lastStepIncomplete;
            trial_data_combined.walking_distance_r = trial_data_r.walkingDistance;
            trial_data_combined.footAngle_r = trial_data_r.footAngle;
            trial_data_combined.kneeAngle_r = trial_data_r.kneeAngle;
            trial_data_combined.HeelR_velocity = mean(trial_data_r.HeelR_velocity);
            trial_data_combined.HeelL_velocity = mean(trial_data_l.HeelL_velocity);
            trial_data_combined.velocity = trial_data_r.velocity;
            trial_data_combined.HeelL_vertical = trial_data_l.HeelL_z;
            trial_data_combined.KneeL_vertical = trial_data_l.KneeL_z;
            trial_data_combined.HeelR_vertical = trial_data_r.HeelR_z;
            trial_data_combined.KneeR_vertical = trial_data_r.KneeR_z;
            
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