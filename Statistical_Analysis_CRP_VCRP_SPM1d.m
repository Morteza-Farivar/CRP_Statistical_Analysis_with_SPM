%% Step 1: Extract CRP Data and Trial Mapping from Workspace (All Gait Cycles)

% Define segment coupling names and trial identifiers
segment_names = {'L_shank_L_foot', 'L_thigh_L_shank', 'L_thigh_pelvis', ...
                 'R_shank_R_foot', 'R_thigh_pelvis', 'R_thigh_R_shank', 'trunk_pelvis'};
trials = {'100', '120', '140', '160', '60', '80'}; 

% Save `trials` to Workspace for later steps
assignin('base', 'trials', trials);

% Initialize storage for extracted CRP and Variability data
crp_data_extracted = struct(); % Structure to store CRP data
var_data_extracted = struct(); % Structure to store Variability data

% Loop through each segment to extract and organize data
for segment_idx = 1:length(segment_names)
    segment_name = segment_names{segment_idx};
    
    % Access CRP and Variability data for the current segment (Y-axis only)
    crp_segment = crp_all.(segment_name).Y; % CRP data
    var_segment = crp_variability.(segment_name).Y; % Variability data

    % Remove Subject 3 (outlier) from both CRP and Variability data
    crp_segment(3, :) = []; % Remove Subject 3 from CRP
    var_segment(3, :) = []; % Remove Subject 3 from Variability

    % Initialize storage for this segment
    crp_data_extracted.(segment_name) = cell(size(crp_segment)); % Cell array to handle varying lengths
    var_data_extracted.(segment_name) = cell(size(var_segment));

    % Loop through each trial and subject to store all gait cycles
    for trial_idx = 1:size(crp_segment, 2)
        for subj_idx = 1:size(crp_segment, 1)
            % Store all gait cycles for the current trial and subject
            crp_data_extracted.(segment_name){subj_idx, trial_idx} = crp_segment{subj_idx, trial_idx};
            var_data_extracted.(segment_name){subj_idx, trial_idx} = var_segment{subj_idx, trial_idx};
        end
    end
    
    % Display progress
    disp(['Data extracted for segment: ', segment_name]);
end

% Save extracted CRP and Variability data to Workspace for further analysis
assignin('base', 'crp_data_extracted', crp_data_extracted);
assignin('base', 'var_data_extracted', var_data_extracted);

disp('Step 1 completed: All gait cycles extracted and saved successfully.');

%% Step 2: Calculate Mean Across Subjects for Each Trial (Optimized and Warning-Free)
% Load extracted data from Step 1
crp_data_extracted = evalin('base', 'crp_data_extracted'); % Extracted CRP data
var_data_extracted = evalin('base', 'var_data_extracted'); % Extracted Variability data
segment_names = evalin('base', 'segment_names'); % Segment coupling names
trials = evalin('base', 'trials'); % Trials

% Initialize storage for mean data
crp_mean_data = struct();
var_mean_data = struct();

% Loop through each segment
for segment_idx = 1:length(segment_names)
    segment_name = segment_names{segment_idx};
    
    % Access CRP and Variability data for this segment
    crp_data = crp_data_extracted.(segment_name);
    var_data = var_data_extracted.(segment_name);
    
    % Check if data is empty
    if isempty(crp_data) || isempty(var_data)
        warning('Data for segment "%s" is empty. Skipping...', segment_name);
        continue;
    end
    
    % Initialize storage for mean data for this segment
    num_trials = length(trials); % Use trials directly
    num_points = size(crp_data{1, 1}, 2); % Number of points in gait cycle
    crp_mean = zeros(num_trials, num_points); % [Trials x Points]
    var_mean = zeros(num_trials, num_points); % [Trials x Points]
    
    % Preallocate for trial data
    max_subjects = size(crp_data, 1); % Maximum number of subjects
    max_gait_cycles = 200; % Estimate maximum number of gait cycles per subject
    trial_crp_data = NaN(max_subjects * max_gait_cycles, num_points); % Preallocate for all gait cycles
    trial_var_data = NaN(max_subjects * max_gait_cycles, num_points); % Preallocate for all gait cycles

    % Loop through each trial
    for trial_idx = 1:num_trials
        trial_name = trials{trial_idx}; % Optional use for clarity
        
        % Reset counters for valid entries
        valid_crp_count = 0;
        valid_var_count = 0;
        
        for subj_idx = 1:size(crp_data, 1) % Number of subjects
            % Append CRP data for the current subject (if not empty)
            if ~isempty(crp_data{subj_idx, trial_idx})
                num_cycles = size(crp_data{subj_idx, trial_idx}, 1);
                trial_crp_data(valid_crp_count + 1:valid_crp_count + num_cycles, :) = crp_data{subj_idx, trial_idx};
                valid_crp_count = valid_crp_count + num_cycles;
            end
            
            % Append Variability data for the current subject (if not empty)
            if ~isempty(var_data{subj_idx, trial_idx})
                num_cycles = size(var_data{subj_idx, trial_idx}, 1);
                trial_var_data(valid_var_count + 1:valid_var_count + num_cycles, :) = var_data{subj_idx, trial_idx};
                valid_var_count = valid_var_count + num_cycles;
            end
        end
        
        % Trim excess preallocated rows
        trial_crp_data = trial_crp_data(1:valid_crp_count, :);
        trial_var_data = trial_var_data(1:valid_var_count, :);
        
        % Calculate mean across all valid entries for the trial
        crp_mean(trial_idx, :) = mean(trial_crp_data, 1, 'omitnan');
        var_mean(trial_idx, :) = mean(trial_var_data, 1, 'omitnan');
    end
    
    % Store the mean data in the structure
    crp_mean_data.(segment_name) = crp_mean;
    var_mean_data.(segment_name) = var_mean;
    
    disp(['Mean calculated for segment: ', segment_name]);
end

% Save mean data to Workspace for further analysis
assignin('base', 'crp_mean_data', crp_mean_data);
assignin('base', 'var_mean_data', var_mean_data);

disp('Step 2 completed: Mean data calculated and saved successfully.');


%% Step 3: Advanced SPM Output Analysis and Reporting (Full Standalone)
% Statistical Parametric Mapping (SPM) without external toolbox dependency
% Author: AutoRefined for portability and compatibility

% Load required data
crp_data_extracted = evalin('base', 'crp_data_extracted');
var_data_extracted = evalin('base', 'var_data_extracted');
segment_names = evalin('base', 'segment_names');
trials = {'100', '120', '140', '160', '60', '80'};
alpha_level = 0.05;
max_gait_cycles = 50;

% Initialize
spm_results = struct();

% Loop through segments
for segment_idx = 1:length(segment_names)
    segment_name = segment_names{segment_idx};
    disp(['Running advanced SPM analysis for segment: ', segment_name]);

    crp_segment = crp_data_extracted.(segment_name);
    var_segment = var_data_extracted.(segment_name);

    % Prepare concatenated data
    total_crp_cycles = 0;
    total_var_cycles = 0;

    for trial_idx = 1:length(trials)
        for subj_idx = 1:size(crp_segment, 1)
            if ~isempty(crp_segment{subj_idx, trial_idx})
                total_crp_cycles = total_crp_cycles + min(size(crp_segment{subj_idx, trial_idx}, 1), max_gait_cycles);
            end
            if ~isempty(var_segment{subj_idx, trial_idx})
                total_var_cycles = total_var_cycles + min(size(var_segment{subj_idx, trial_idx}, 1), max_gait_cycles);
            end
        end
    end

    % Initialize matrices
    num_points = size(crp_segment{1, 1}, 2);
    concatenated_crp = zeros(total_crp_cycles, num_points);
    concatenated_var = zeros(total_var_cycles, num_points);
    trial_ids_crp = zeros(total_crp_cycles, 1);
    subj_ids_crp = zeros(total_crp_cycles, 1);
    trial_ids_var = zeros(total_var_cycles, 1);
    subj_ids_var = zeros(total_var_cycles, 1);

    start_idx_crp = 1;
    start_idx_var = 1;

    for trial_idx = 1:length(trials)
        for subj_idx = 1:size(crp_segment, 1)
            crp_data = crp_segment{subj_idx, trial_idx};
            var_data = var_segment{subj_idx, trial_idx};

            if ~isempty(crp_data)
                N = min(size(crp_data, 1), max_gait_cycles);
                concatenated_crp(start_idx_crp:start_idx_crp+N-1, :) = crp_data(1:N, :);
                trial_ids_crp(start_idx_crp:start_idx_crp+N-1) = trial_idx;
                subj_ids_crp(start_idx_crp:start_idx_crp+N-1) = subj_idx;
                start_idx_crp = start_idx_crp + N;
            end

            if ~isempty(var_data)
                N = min(size(var_data, 1), max_gait_cycles);
                concatenated_var(start_idx_var:start_idx_var+N-1, :) = var_data(1:N, :);
                trial_ids_var(start_idx_var:start_idx_var+N-1) = trial_idx;
                subj_ids_var(start_idx_var:start_idx_var+N-1) = subj_idx;
                start_idx_var = start_idx_var + N;
            end
        end
    end

    % ----------------------
    % ANOVA (Manual SPM replacement)
    % ----------------------

    % Perform ANOVA for CRP
    [p_vals_crp, F_vals_crp] = deal(nan(1, num_points));
    for t = 1:num_points
        tbl = table(trial_ids_crp, subj_ids_crp, concatenated_crp(:, t), ...
                    'VariableNames', {'Trial', 'Subject', 'Value'});
        rm = fitrm(tbl, 'Value~Trial+Subject');
        ranova_tbl = ranova(rm);
        p_vals_crp(t) = ranova_tbl.pValue(1);
        F_vals_crp(t) = ranova_tbl.F(1);
    end

    % Perform ANOVA for Variability
    [p_vals_var, F_vals_var] = deal(nan(1, num_points));
    for t = 1:num_points
        tbl = table(trial_ids_var, subj_ids_var, concatenated_var(:, t), ...
                    'VariableNames', {'Trial', 'Subject', 'Value'});
        rm = fitrm(tbl, 'Value~Trial+Subject');
        ranova_tbl = ranova(rm);
        p_vals_var(t) = ranova_tbl.pValue(1);
        F_vals_var(t) = ranova_tbl.F(1);
    end

    % Save Results
    spm_results.(segment_name).CRP.p = p_vals_crp;
    spm_results.(segment_name).CRP.F = F_vals_crp;
    spm_results.(segment_name).Variability.p = p_vals_var;
    spm_results.(segment_name).Variability.F = F_vals_var;

    % Plot Results
    figure;
    subplot(2, 1, 1);
    plot(F_vals_crp, 'LineWidth', 2);
    title(['ANOVA-F for CRP - ', segment_name]);
    xlabel('Gait Cycle (%)'); ylabel('F-statistic');
    yline(finv(1-alpha_level, length(trials)-1, total_crp_cycles-length(trials)), '--r', 'Threshold');

    subplot(2, 1, 2);
    plot(p_vals_crp, 'LineWidth', 2);
    title(['p-values for CRP - ', segment_name]);
    xlabel('Gait Cycle (%)'); ylabel('p-value');
    yline(alpha_level, '--r', 'alpha=0.05');

    figure;
    subplot(2, 1, 1);
    plot(F_vals_var, 'LineWidth', 2);
    title(['ANOVA-F for Variability - ', segment_name]);
    xlabel('Gait Cycle (%)'); ylabel('F-statistic');
    yline(finv(1-alpha_level, length(trials)-1, total_var_cycles-length(trials)), '--r', 'Threshold');

    subplot(2, 1, 2);
    plot(p_vals_var, 'LineWidth', 2);
    title(['p-values for Variability - ', segment_name]);
    xlabel('Gait Cycle (%)'); ylabel('p-value');
    yline(alpha_level, '--r', 'alpha=0.05');
end

assignin('base', 'spm_results', spm_results);
disp('✅ Step 3 completed: Advanced ANOVA-based results saved.');

%% Step 4: Pairwise Comparison using SPM T-test (Final Version)
% Load necessary data from Workspace
crp_data_extracted = evalin('base', 'crp_data_extracted'); % Extracted CRP data
var_data_extracted = evalin('base', 'var_data_extracted'); % Extracted Variability data
segment_names = evalin('base', 'segment_names'); % Segment coupling names
trials = evalin('base', 'trials'); % Trials
alpha_level = 0.05; % Significance level
max_gait_cycles = 50; % Limit to 50 gait cycles per trial

% Generate all pairwise combinations
num_trials = length(trials);
pairwise_combinations = nchoosek(1:num_trials, 2); % All pair combinations (15 pairs)

% Initialize storage for results
spm_pairwise_results = struct();

% Loop through each segment
for segment_idx = 1:length(segment_names)
    segment_name = segment_names{segment_idx};
    disp(['Running pairwise comparisons for segment: ', segment_name]);

    % Access CRP and Variability data for the current segment
    crp_segment = crp_data_extracted.(segment_name);
    var_segment = var_data_extracted.(segment_name);

    % Initialize storage for this segment
    spm_pairwise_results.(segment_name).CRP = [];
    spm_pairwise_results.(segment_name).Variability = [];

    % Loop through each pair of trials
    for pair_idx = 1:size(pairwise_combinations, 1)
        trial_a_idx = pairwise_combinations(pair_idx, 1);
        trial_b_idx = pairwise_combinations(pair_idx, 2);

        trial_a_name = trials{trial_a_idx};
        trial_b_name = trials{trial_b_idx};
        disp(['Comparing Trial ', trial_a_name, ' vs Trial ', trial_b_name]);

        % Extract and concatenate data for the two trials
        [data_a_crp, data_b_crp] = extract_gait_cycle_data(crp_segment, trial_a_idx, trial_b_idx, max_gait_cycles);
        [data_a_var, data_b_var] = extract_gait_cycle_data(var_segment, trial_a_idx, trial_b_idx, max_gait_cycles);

        % ----------------------
        % SPM T-test for CRP
        % ----------------------
        spm_crp = spm1d.stats.ttest2(data_a_crp, data_b_crp);
        spmi_crp = spm_crp.inference(alpha_level);

        % Save CRP results
        spm_pairwise_results.(segment_name).CRP(pair_idx).Trials = [trial_a_name, ' vs ', trial_b_name];
        spm_pairwise_results.(segment_name).CRP(pair_idx).Clusters = spmi_crp.clusters;
        spm_pairwise_results.(segment_name).CRP(pair_idx).p_values = spmi_crp.p_set;

        % Visualization for CRP
        figure;
        spmi_crp.plot();
        title(['SPM T-test for CRP - ', segment_name, ' (', trial_a_name, ' vs ', trial_b_name, ')']);
        xlabel('Gait Cycle (%)');
        ylabel('SPM{t}');
        spmi_crp.plot_threshold_label();
        spmi_crp.plot_p_values();

        % ----------------------
        % SPM T-test for Variability
        % ----------------------
        spm_var = spm1d.stats.ttest2(data_a_var, data_b_var);
        spmi_var = spm_var.inference(alpha_level);

        % Save Variability results
        spm_pairwise_results.(segment_name).Variability(pair_idx).Trials = [trial_a_name, ' vs ', trial_b_name];
        spm_pairwise_results.(segment_name).Variability(pair_idx).Clusters = spmi_var.clusters;
        spm_pairwise_results.(segment_name).Variability(pair_idx).p_values = spmi_var.p_set;

        % Visualization for Variability
        figure;
        spmi_var.plot();
        title(['SPM T-test for Variability - ', segment_name, ' (', trial_a_name, ' vs ', trial_b_name, ')']);
        xlabel('Gait Cycle (%)');
        ylabel('SPM{t}');
        spmi_var.plot_threshold_label();
        spmi_var.plot_p_values();

        disp(['Pairwise comparison completed for trials: ', trial_a_name, ' vs ', trial_b_name]);
    end
end

% Save results to Workspace
assignin('base', 'spm_pairwise_results', spm_pairwise_results);

disp('All pairwise comparisons completed successfully. Results saved to Workspace.');