% Extract and save mean force plate data
mean_force_plate_data = struct();

% Initialize arrays for each condition
condition_arrays = struct('Baseline', zeros(8000,1), 'Large', zeros(8000,1), 'NoAvatar', zeros(8000,1), 'Normal', zeros(8000,1), 'Small', zeros(8000,1));

conditions = {"Baseline", "Large", "NoAvatar", "Normal", "Small"};
condition_start_indices = [1, 6, 16, 26, 36];
num_trials_per_condition = [5, 10, 10, 10, 10];

for cond_idx = 1:5
    condition_name = conditions{cond_idx};
    start_idx = condition_start_indices(cond_idx);
    end_idx = start_idx + num_trials_per_condition(cond_idx) - 1;

    % Collect data for all participants and trials
    data_fp1 = zeros(8000, 0);
    data_fp2 = zeros(8000, 0);
    data_fp3 = zeros(8000, 0);
    data_fp4 = zeros(8000, 0);
    
    for participant = 9:17
        for trial = start_idx:end_idx
            trial_data = out_fd_off{trial, participant};
            data_fp1 = [data_fp1, trial_data(:,4)];
            data_fp2 = [data_fp2, trial_data(:,19)];
            data_fp3 = [data_fp3, trial_data(:,14)];
            data_fp4 = [data_fp4, trial_data(:,9)];
        end
    end

    % Compute mean across trials for each force plate
    condition_arrays.(condition_name)(:,1) = mean(data_fp1, 2);
    condition_arrays.(condition_name)(:,2) = mean(data_fp2, 2);
    condition_arrays.(condition_name)(:,3) = mean(data_fp3, 2);
    condition_arrays.(condition_name)(:,4) = mean(data_fp4, 2);
end

% Save averaged force plate data as JSON
jsonText_fp = jsonencode(condition_arrays, 'PrettyPrint', true);
fileID_fp = fopen('mean_force_plate_data.json', 'w');
fwrite(fileID_fp, jsonText_fp);
fclose(fileID_fp);

disp('Mean force plate data export completed.');
%%
% Define consistent colors for conditions
condition_colors = lines(5);
fp_labels = {"FP1", "FP2", "FP3", "FP4"};

% Plot pairwise comparisons
pairwise_comparisons = {
    {'Baseline', 'Large'},
    {'Baseline', 'NoAvatar'},
    {'Baseline', 'Normal'},
    {'Baseline', 'Small'},
    {'Normal', 'Large'},
    {'Normal', 'NoAvatar'},
    {'Normal', 'Small'}
};

for i = 1:length(pairwise_comparisons)
    comparison = pairwise_comparisons{i};
    figure;
    hold on;
    % Plot each force plate with consistent colors
    for fp_idx = 1:4
        plot(linspace(0, 8, 8000), condition_arrays.(comparison{1})(:,fp_idx), 'DisplayName', sprintf('%s (%s)', comparison{1}, fp_labels{fp_idx}), 'Color', [condition_colors(fp_idx,:) 0.3], 'LineWidth', 1.5);
        plot(linspace(0, 8, 8000), condition_arrays.(comparison{2})(:,fp_idx), 'DisplayName', sprintf('%s (%s)', comparison{2}, fp_labels{fp_idx}), 'Color', condition_colors(fp_idx,:), 'LineWidth', 2.5);
    end
    xlabel('Time (s)');
    ylabel('Ground Reaction Force (N)');
    legend('show');
    title(sprintf('%s vs %s - Ground Reaction Force', comparison{1}, comparison{2}));
    saveas(gcf, sprintf('%s_vs_%s_Force_Plot.png', comparison{1}, comparison{2}));
    hold off;
end

% Plot all conditions together
figure;
hold on;
handles = [];
labels = {};
for cond_idx = 1:5
    condition_name = conditions{cond_idx};
    for fp_idx = 1:4
        h = plot(linspace(0, 8, 8000), condition_arrays.(condition_name)(:,fp_idx), 'Color', condition_colors(cond_idx,:), 'LineWidth', 2);
        if fp_idx == 1 % Add only one legend entry per condition
            handles = [handles, h];
            labels = [labels, condition_name];
        end
    end
end
xlabel('Time (s)');
ylabel('Ground Reaction Force (N)');
legend(handles, labels, 'Location', 'best');
title('All Conditions - Ground Reaction Force');
saveas(gcf, 'All_Conditions_Force_Plot.png');
hold off;

disp('Pairwise and combined plotting completed.');

