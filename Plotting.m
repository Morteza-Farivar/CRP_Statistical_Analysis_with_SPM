%% Combined Analysis and Plotting for Specified Segments
% Load necessary data from Workspace
crp_data_extracted = evalin('base', 'crp_data_extracted'); % Extracted CRP data
var_data_extracted = evalin('base', 'var_data_extracted'); % Extracted Variability data
crp_mean_data = evalin('base', 'crp_mean_data'); % Mean CRP data
segment_names = {'R_thigh_R_shank', 'R_shank_R_foot', 'R_thigh_pelvis'}; % Specified segments
alpha_level = 0.05; % Significance level
max_gait_cycles = 50; % Limit to 50 gait cycles per trial

% Define original trial labels and desired walking speed labels with correct order
trial_names_original = {'Trial_100', 'Trial_120', 'Trial_140', 'Trial_160', 'Trial_60', 'Trial_80'}; % Original order
trial_names_ws = {'WS 60', 'WS 80', 'WS 100', 'WS 120', 'WS 140', 'WS 160'}; % Correct order

% Mapping for reordering the legend
desired_order = [5, 6, 1, 2, 3, 4]; % Ensuring correct mapping
trial_colors = lines(length(trial_names_original)); % Keeping the correct colors

% Initialize storage for SPM results
spm_results = struct();

% Predefine labels for combined figure
labels = {'A', 'B', 'C', 'D', 'E', 'F'}; % Labels for subplots

% Create a single figure for combined plots
figure('Color', 'w');
num_segments = length(segment_names); % Number of segments

% Loop through each segment for calculations and plotting
for i = 1:num_segments
    segment_name = segment_names{i};
    disp(['Processing segment: ', segment_name]);

    % ---------------------- SPM Analysis ----------------------
    crp_segment = crp_data_extracted.(segment_name);
    var_segment = var_data_extracted.(segment_name);

    % Initialize matrices for concatenated data
    num_points = size(crp_segment{1, 1}, 2); % Number of gait cycle points
    concatenated_crp = [];
    concatenated_var = [];
    trial_ids_crp = [];
    subj_ids_crp = [];
    trial_ids_var = [];
    subj_ids_var = [];

    % Concatenate CRP and Variability data
    for trial_idx = 1:length(trial_names_original)
        for subj_idx = 1:size(crp_segment, 1)
            crp_gait_cycles = crp_segment{subj_idx, trial_idx};
            var_gait_cycles = var_segment{subj_idx, trial_idx};

            if ~isempty(crp_gait_cycles)
                num_cycles = min(size(crp_gait_cycles, 1), max_gait_cycles);
                concatenated_crp = [concatenated_crp; crp_gait_cycles(1:num_cycles, :)];
                trial_ids_crp = [trial_ids_crp; repmat(trial_idx, num_cycles, 1)];
                subj_ids_crp = [subj_ids_crp; repmat(subj_idx, num_cycles, 1)];
            end

            if ~isempty(var_gait_cycles)
                num_cycles = min(size(var_gait_cycles, 1), max_gait_cycles);
                concatenated_var = [concatenated_var; var_gait_cycles(1:num_cycles, :)];
                trial_ids_var = [trial_ids_var; repmat(trial_idx, num_cycles, 1)];
                subj_ids_var = [subj_ids_var; repmat(subj_idx, num_cycles, 1)];
            end
        end
    end

    % Perform SPM analysis for CRP
    spm_crp = spm1d.stats.anova1rm(concatenated_crp, trial_ids_crp, subj_ids_crp);
    spmi_crp = spm_crp.inference(alpha_level);

    % Perform SPM analysis for Variability
    spm_var = spm1d.stats.anova1rm(concatenated_var, trial_ids_var, subj_ids_var);
    spmi_var = spm_var.inference(alpha_level);

    % ---------------------- Combined Plotting ----------------------
    % CRP Plot
    subplot(num_segments, 2, (i-1)*2 + 1);
    segment_mean_data = crp_mean_data.(segment_name); % [Trials x Points]
    num_trials = size(segment_mean_data, 1);

    hold on;
    plot_handles = gobjects(num_trials, 1); % Store handles for legend
    for trial_idx = 1:num_trials
        plot_handles(trial_idx) = plot(1:num_points, segment_mean_data(trial_idx, :), 'LineWidth', 1.5, ...
                                       'Color', trial_colors(trial_idx, :));
    end
    hold off;

    % Fix Legend with correct order
    legend(plot_handles(desired_order), trial_names_ws, 'Location', 'northeast', 'FontSize', 4);

    title(['CRP - ', strrep(segment_name, '_', ' '), ' - FW']);
    xlabel('Gait Cycle (%)');
    ylabel('CRP (degrees)');
    box on;
    grid off;

    % Add label to subplot (e.g., A, C, E)
    annotation('textbox', [0.1, 0.9 - (i-1)*0.3, 0.02, 0.02], ...
               'String', labels{(i-1)*2 + 1}, ...
               'FontWeight', 'bold', 'FontSize', 10, 'EdgeColor', 'none');

    % ---------------------- SPM Plot ----------------------
    subplot(num_segments, 2, (i-1)*2 + 2);
    spmi_crp.plot();
    spmi_crp.plot_threshold_label();
    spmi_crp.plot_p_values();

    % Adjust p-value text size
    ax = gca;
    for t_idx = 1:length(ax.Children)
        if isa(ax.Children(t_idx), 'matlab.graphics.primitive.Text')
            ax.Children(t_idx).FontSize = 6; % Reduce font size for p-values
        end
    end

    title(['SPM Results for CRP - ', strrep(segment_name, '_', ' ')]);
    xlabel('Gait Cycle (%)');
    ylabel('SPM{F}');
    box on;

    % Add label to subplot (e.g., B, D, F)
    annotation('textbox', [0.55, 0.9 - (i-1)*0.3, 0.02, 0.02], ...
               'String', labels{(i-1)*2 + 2}, ...
               'FontWeight', 'bold', 'FontSize', 10, 'EdgeColor', 'none');
end

disp('All CRP and SPM plots generated successfully in a single figure!');