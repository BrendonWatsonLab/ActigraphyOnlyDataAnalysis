% MATLAB Script for Generating Figure 5 Components as Separate Figures
% Date: Oct 08, 2025
% Author: Noah Muscat & ActigraphyOnly Research Assistant
%
% Description:
% REVISED per PI feedback.
% - 5A: Improved significance bar visualization.
% - 5B: Rebuilt as separate M/F plots with 4 violins each; stats added.
% - 5C: REMOVED.
% - 5D: REBUILT as a 24h difference bar plot (last 2 weeks); interaction stats added.
%
% ASSUMPTIONS: violinplot.m is on the MATLAB path.
% ---------------------------------------------------------------------
%% --- Configuration & Setup ---
clear; clc; close all;
% --- File and Path Parameters ---
filename = '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/ActivityAnalysis/ActigraphyOnly/Unified_AO1to12_Combined_Data.csv'; 
savePath = '/Users/noahmuscat/Desktop/WatsonLab/AOPaper/Figure5Plots';
% --- Key Data Column Names ---
normalizedActivityVar = 'NormalizedActivity';
pixelDiffVar = 'SelectedPixelDifference';
animalVar = 'Animal';
conditionVar = 'Condition';
relativeDayVar = 'RelativeDay';
timeVarZT = 'DateZT';
ztHourVar = 'ZT_Time'; % Integer ZT hour (0-23)
% --- Analysis Parameters ---
hoursPerDay = 24;
ztLightsOnStart = 0; ztLightsOnEnd = 11;
ztLightsOffStart = 12; ztLightsOffEnd = 23;
% --- Animal Grouping ---
allAnimals = string({'AO1', 'AO2', 'AO3', 'AO4', 'AO5', 'AO6', 'AO7', 'AO8', 'AO9', 'AO10', 'AO11', 'AO12'});
maleAnimals = string({'AO1', 'AO2', 'AO3', 'AO7', 'AO9', 'AO12'});
femaleAnimals = string({'AO4', 'AO5', 'AO6', 'AO8', 'AO10', 'AO11'});
fourCondAnimals = string({'AO5', 'AO6', 'AO7', 'AO8', 'AO9','AO10','AO11','AO12'}); 
% Define colors
maleColor = [0.2 0.5470 0.8410]; 
femaleColor = [0.9500 0.4250 0.1980]; 
grayColor = [0.6 0.6 0.6];
lightingConditionsFig5_Overall = {'300Lux', '1000Lux', 'FullDark', '300LuxEnd'};
%% --- Setup Save Directory ---
disp(['Figures will be saved to: ', savePath]);
if ~isfolder(savePath), try mkdir(savePath); disp('Save directory created.'); catch ME, error('SaveDirErr:CouldNotCreate', 'Could not create save directory "%s". Err: %s', savePath, ME.message); end, end
%% --- Load Data ---
disp(['Loading data from: ', filename]);
try
    opts = detectImportOptions(filename); opts.VariableNamingRule = 'preserve';
    expectedVars = {pixelDiffVar, normalizedActivityVar, animalVar, conditionVar, relativeDayVar, timeVarZT, ztHourVar};
    if ~all(ismember(expectedVars, opts.VariableNames)), error('LoadErr:MissingColumn', 'One or more expected columns not found in %s.', filename); end
    
    opts = setvartype(opts, {pixelDiffVar, normalizedActivityVar, relativeDayVar, ztHourVar}, 'double');
    opts = setvartype(opts, {animalVar, conditionVar}, 'string');
    opts = setvartype(opts, timeVarZT, 'datetime');
    opts = setvaropts(opts, timeVarZT, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    
    dataTable = readtable(filename, opts);
    disp('Data loaded successfully.');
    dataTable.DayOfCondition = floor(dataTable.(relativeDayVar));
catch ME_load, error('LoadErr:FileProcessing', 'Failed to load/parse CSV "%s". Err: %s', filename, ME_load.message); end
%% --- Data Preprocessing ---
disp('Starting data preprocessing...');
dataTable = sortrows(dataTable, timeVarZT);
colsToClean = {pixelDiffVar, normalizedActivityVar, timeVarZT, relativeDayVar, ztHourVar, animalVar, conditionVar, 'DayOfCondition'};
nanRowsFilter = false(height(dataTable), 1);
for c_idx = 1:length(colsToClean), colData = dataTable.(colsToClean{c_idx}); if isdatetime(colData), nanRowsFilter = nanRowsFilter | isnat(colData); elseif isnumeric(colData), nanRowsFilter = nanRowsFilter | isnan(colData); elseif isstring(colData), nanRowsFilter = nanRowsFilter | ismissing(colData); end, end
if any(nanRowsFilter), fprintf('%d rows removed due to missing values.\n', sum(nanRowsFilter)); dataTable(nanRowsFilter, :) = []; end
if isempty(dataTable), error('PreprocErr:NoData', 'No valid data after NaN removal.'); end
dataTable = dataTable(ismember(dataTable.(animalVar), allAnimals), :);
if isempty(dataTable), error('PreprocErr:NoAnimalData', 'No data found for the animals specified in the `allAnimals` list.'); end
dataTable.Condition = categorical(dataTable.Condition, lightingConditionsFig5_Overall, 'Ordinal',true);
dataTable.Sex = categorical(ismember(dataTable.(animalVar), maleAnimals), [0, 1], {'Female', 'Male'});
dataTable.LightPeriod = categorical((dataTable.(ztHourVar) >= ztLightsOnStart & dataTable.(ztHourVar) <= ztLightsOnEnd), [0, 1], {'Off', 'On'});
disp('Data preprocessing complete.');

%% =====================================================================
% GENERATING FIGURE 5 COMPONENTS AS SEPARATE FIGURES
% =====================================================================
activity_metrics = {normalizedActivityVar, pixelDiffVar};
metric_suffixes = {'Normalized', 'PixelDiff'};
metric_ylabels = {'Normalized Activity', 'Pixel Difference'};

%% FIG 5A: Sex differences in 24h profiles by Condition (WITH STATS)
disp('--- Starting Figure 5A Components ---');
fprintf('\n--- STATS: Figure 5A (300Lux vs 1000Lux Hourly) ---\n');
for met_idx = 1:length(activity_metrics)
    current_metric_var = activity_metrics{met_idx}; current_metric_suffix = metric_suffixes{met_idx}; current_metric_ylabel = metric_ylabels{met_idx};
    sex_groups_5A = {maleAnimals, femaleAnimals}; sex_labels_5A = {'Male', 'Female'};
    conditionsFor5A_plot = {'300Lux', '1000Lux', 'FullDark'};
    
    % Pre-calculation loop to find global Y-axis limits for each metric
    global_min_y = Inf; global_max_y = -Inf;
    for sex_idx_pre = 1:length(sex_groups_5A)
        current_sex_animals_pre = sex_groups_5A{sex_idx_pre};
        for i_cond_pre = 1:length(conditionsFor5A_plot)
            currentCondStr_pre = conditionsFor5A_plot{i_cond_pre};
            animalsForThisProfile_pre = current_sex_animals_pre;
            if strcmp(currentCondStr_pre, 'FullDark'), animalsForThisProfile_pre = intersect(current_sex_animals_pre, fourCondAnimals); end
            if isempty(animalsForThisProfile_pre), continue; end
            sexCondData_pre = dataTable(ismember(dataTable.Animal, animalsForThisProfile_pre) & dataTable.Condition == string(currentCondStr_pre), :);
            if isempty(sexCondData_pre), continue; end
            hourlyStats_pre = groupsummary(sexCondData_pre, ztHourVar, 'mean', current_metric_var);
            global_min_y = min(global_min_y, min(hourlyStats_pre.(['mean_',current_metric_var])));
            global_max_y = max(global_max_y, max(hourlyStats_pre.(['mean_',current_metric_var])));
        end
    end
    y_range = global_max_y - global_min_y; y_lims_final = [global_min_y - 0.1*y_range, global_max_y + 0.15*y_range]; % Increased top buffer
    
    for sex_idx_5A = 1:length(sex_groups_5A)
        current_sex_animals_5A = sex_groups_5A{sex_idx_5A}; 
        current_sex_label_5A = sex_labels_5A{sex_idx_5A};
        
        hFig5A_current = figure('Name', sprintf('Figure 5A: %s Profiles (%s)', current_sex_label_5A, current_metric_suffix),'NumberTitle', 'off', 'Color', 'w','Units', 'normalized', 'OuterPosition', [0.1+sex_idx_5A*0.05 0.1+met_idx*0.05 0.4 0.6]);
        ax_current_5A = axes('Parent', hFig5A_current); hold(ax_current_5A, 'on');
        lineColors_cond_5A = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]}; lineStyles_cond_5A = {'-', '--', ':'};
        plotHandles_current_5A = gobjects(length(conditionsFor5A_plot),1);
        for i_cond = 1:length(conditionsFor5A_plot)
            currentCondStr = conditionsFor5A_plot{i_cond}; animalsForThisProfile = current_sex_animals_5A;
            if strcmp(currentCondStr, 'FullDark'), animalsForThisProfile = intersect(current_sex_animals_5A, fourCondAnimals); end
            if isempty(animalsForThisProfile), plotHandles_current_5A(i_cond) = plot(ax_current_5A, NaN, NaN, 'DisplayName', sprintf('%s (No Data)', currentCondStr), 'HandleVisibility','off'); continue; end
            sexCondData = [];
            if strcmp(currentCondStr, '300Lux')
                tempData = dataTable(ismember(dataTable.(animalVar), animalsForThisProfile),:);
                for i_an_5a = 1:length(animalsForThisProfile), animalDataOnly = tempData(tempData.(animalVar) == animalsForThisProfile(i_an_5a),:); animalCondsOrder = unique(animalDataOnly.Condition,'stable'); if ~isempty(animalCondsOrder) && animalCondsOrder(1) == "300Lux", sexCondData = [sexCondData; animalDataOnly(animalDataOnly.Condition == "300Lux", :)]; end, end
            else, sexCondData = dataTable(ismember(dataTable.(animalVar), animalsForThisProfile) & dataTable.Condition == string(currentCondStr), :); end
            if isempty(sexCondData), plotHandles_current_5A(i_cond) = plot(ax_current_5A, NaN, NaN, 'DisplayName', sprintf('%s (No Data)', currentCondStr), 'HandleVisibility','off'); continue; end
            hourlyStats = groupsummary(sexCondData, ztHourVar, {'mean','std','nnz'}, current_metric_var); hourlyStats = sortrows(hourlyStats, ztHourVar);
            meanAct = hourlyStats.(['mean_',current_metric_var]); stdAct = hourlyStats.(['std_',current_metric_var]); countAct = hourlyStats.(['nnz_',current_metric_var]); semAct = stdAct ./ sqrt(countAct); semAct(countAct <= 1) = NaN;
            ztH = hourlyStats.(ztHourVar); fullDayZT=(0:23)'; meanProfile=NaN(hoursPerDay,1); semProfile=NaN(hoursPerDay,1);
            [lia,locb]=ismember(fullDayZT,ztH); meanProfile(lia)=meanAct(locb(lia)); if ~all(isnan(semAct)), semProfile(lia)=semAct(locb(lia)); end
            plotHandles_current_5A(i_cond) = plot(ax_current_5A, fullDayZT, meanProfile, 'LineStyle', lineStyles_cond_5A{i_cond}, 'LineWidth', 1.5, 'Color', lineColors_cond_5A{i_cond}, 'DisplayName', currentCondStr);
            if ~all(isnan(semProfile)), x_patch=[fullDayZT;flipud(fullDayZT)]; y_patch=[meanProfile+semProfile;flipud(meanProfile-semProfile)]; valid_patch=~isnan(x_patch)&~isnan(y_patch); if any(valid_patch), fill(ax_current_5A,x_patch(valid_patch),y_patch(valid_patch),lineColors_cond_5A{i_cond}, 'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off'); end, end
        end
        
        ylim(ax_current_5A, y_lims_final);
        
        yLims_curr=ylim(ax_current_5A); patch(ax_current_5A,[ztLightsOffStart,ztLightsOffEnd+1,ztLightsOffEnd+1,ztLightsOffStart],[yLims_curr(1),yLims_curr(1),yLims_curr(2),yLims_curr(2)],[0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        uistack(plotHandles_current_5A(isgraphics(plotHandles_current_5A)),'top');
        if strcmp(current_metric_var, normalizedActivityVar), plot(ax_current_5A,xlim(ax_current_5A),[0 0],'k--','LineWidth',0.75,'HandleVisibility','off'); end
        
        % Stats for 300 vs 1000 Lux by hour
        fprintf('\n-- Metric: %s, Sex: %s --\n', current_metric_suffix, current_sex_label_5A);
        p_values_5A = ones(24, 1);
        for zt = 0:23
            d300 = dataTable(ismember(dataTable.Animal, current_sex_animals_5A) & dataTable.Condition=='300Lux' & dataTable.ZT_Time==zt, :).(current_metric_var);
            d1000 = dataTable(ismember(dataTable.Animal, current_sex_animals_5A) & dataTable.Condition=='1000Lux' & dataTable.ZT_Time==zt, :).(current_metric_var);
            if ~isempty(d300) && ~isempty(d1000), [~, p] = ttest2(d300, d1000); p_values_5A(zt+1) = p; end
        end
        p_corrected_5A = p_values_5A * 24; % Bonferroni Correction
        significant_hours = find(p_corrected_5A < 0.05);
        fprintf('Significant differences (p < 0.05, Bonferroni) at ZT: %s\n', num2str(significant_hours' - 1));

        % <<< CHANGED: Plot significance bar as continuous blocks
        if ~isempty(significant_hours)
            y_pos = y_lims_final(2) - 0.05*y_range;
            block_starts = significant_hours([true; diff(significant_hours)>1]);
            block_ends = significant_hours([diff(significant_hours)>1; true]);
            for i = 1:length(block_starts)
                plot(ax_current_5A, [block_starts(i)-1-0.5, block_ends(i)-1+0.5], [y_pos, y_pos], 'k-', 'LineWidth', 2.5);
            end
        end

        hold(ax_current_5A,'off'); xlabel(ax_current_5A,'ZT Hour','FontSize',10); ylabel(ax_current_5A,['Mean ', strrep(current_metric_ylabel,'_',' ')],'FontSize',10); title(ax_current_5A,sprintf('Fig 5A: %s Profiles (%s)', current_sex_label_5A, current_metric_suffix),'FontSize',11);
        xticks(ax_current_5A,0:6:23);xlim(ax_current_5A,[-0.5,23.5]);grid(ax_current_5A,'on');set(ax_current_5A,'FontSize',9); legend(ax_current_5A, plotHandles_current_5A(isgraphics(plotHandles_current_5A)), 'Location','northwest','FontSize',8);
        fig_filename = fullfile(savePath, sprintf('Figure5A_SameY_%s_%s.png', current_sex_label_5A, current_metric_suffix));
        exportgraphics(hFig5A_current, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]); close(hFig5A_current);
    end
end
disp('--- Finished Figure 5A Components ---');

%% FIG 5B: REBUILT as separate M/F plots
disp('--- Starting Figure 5B Components (Rebuilt) ---');
if ~exist('violinplot.m','file'),error('ViolinPlotNotFound:Fig5B','violinplot.m is required.');end
fprintf('\n--- STATS: Figure 5B (On vs Off, 300 vs 1000) ---\n');
for met_idx = 1:length(activity_metrics)
    current_metric_var = activity_metrics{met_idx}; current_metric_suffix = metric_suffixes{met_idx}; current_metric_ylabel = metric_ylabels{met_idx};
    sex_groups_5B = {maleAnimals, femaleAnimals}; sex_labels_5B = {'Male', 'Female'};
    
    for sex_idx = 1:length(sex_groups_5B)
        current_sex_animals = sex_groups_5B{sex_idx};
        current_sex_label = sex_labels_5B{sex_idx};
        
        hFig5B_current = figure('Name', sprintf('Fig 5B: %s Diurnality (%s)', current_sex_label, current_metric_suffix), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized');
        ax_current_5B = axes('Parent', hFig5B_current); hold(ax_current_5B, 'on');
        
        % <<< CHANGED: Draw the patch first to act as a true background
        p = patch(ax_current_5B, [2.5, 4.5, 4.5, 2.5], [0 0 1 1], [0.9 0.9 0.9], 'FaceAlpha',0.3, 'EdgeColor','none', 'HandleVisibility','off');
        uistack(p, 'bottom'); % Send to very back, behind grid lines
        
        conds = {'300Lux', '1000Lux'};
        periods = {'On', 'Off'};
        group_order = {'300Lux-On', '1000Lux-On', '300Lux-Off', '1000Lux-Off'};
        plot_data_vec = []; plot_group_vec = []; all_animal_means = table();
        for i_c = 1:length(conds)
            for i_p = 1:length(periods)
                groupName = sprintf('%s-%s', conds{i_c}, periods{i_p});
                data_segment = dataTable(ismember(dataTable.Animal, current_sex_animals) & dataTable.Condition==conds{i_c} & dataTable.LightPeriod==periods{i_p}, :);
                if isempty(data_segment), continue; end
                
                animal_means = groupsummary(data_segment, animalVar, 'mean', current_metric_var);
                
                plot_data_vec = [plot_data_vec; animal_means.(['mean_', current_metric_var])];
                plot_group_vec = [plot_group_vec; repmat(categorical({groupName}), height(animal_means), 1)];
                animal_means.Properties.VariableNames{animalVar} = 'Animal';
                animal_means.Properties.VariableNames{end} = groupName;
                if isempty(all_animal_means)
                    all_animal_means = animal_means(:,{'Animal', groupName});
                else
                    all_animal_means = outerjoin(all_animal_means, animal_means(:,{'Animal', groupName}), 'Keys','Animal','MergeKeys',true);
                end
            end
        end
        
        if isempty(plot_data_vec), close(hFig5B_current); continue; end

        violinplot(plot_data_vec, plot_group_vec, 'GroupOrder', group_order, 'ViolinColor', grayColor, 'ShowData', false, 'Parent', ax_current_5B);
        
        % <<< CHANGED: Update patch Y-data after plotting
        new_ylim = ylim(ax_current_5B);
        p.YData = [new_ylim(1), new_ylim(1), new_ylim(2), new_ylim(2)];
        
        y_max_val = max(plot_data_vec,[],'all','omitnan');
        y_range_val = diff(ylim(ax_current_5B));
        if y_range_val == 0, y_range_val = y_max_val; end
        
        [~,p_300] = ttest(all_animal_means.('300Lux-On'), all_animal_means.('300Lux-Off'));
        plot_sig_bar(ax_current_5B, [1, 3], p_300, y_max_val + 0.1*y_range_val);
        
        [~,p_1000] = ttest(all_animal_means.('1000Lux-On'), all_animal_means.('1000Lux-Off'));
        plot_sig_bar(ax_current_5B, [2, 4], p_1000, y_max_val + 0.25*y_range_val);
        
        [~,p_on] = ttest2(all_animal_means.('300Lux-On'), all_animal_means.('1000Lux-On'));
        plot_sig_bar(ax_current_5B, [1, 2], p_on, y_max_val + 0.1*y_range_val);
        
        [~,p_off] = ttest2(all_animal_means.('300Lux-Off'), all_animal_means.('1000Lux-Off'));
        plot_sig_bar(ax_current_5B, [3, 4], p_off, y_max_val + 0.1*y_range_val);
        
        title(sprintf('Fig 5B: %s Diurnality (%s)', current_sex_label, current_metric_suffix));
        ylabel(['Mean ' strrep(current_metric_ylabel, '_', ' ')]);
        ylim([min(0, min(plot_data_vec,[],'all','omitnan') - 0.1*y_range_val) y_max_val + 0.35*y_range_val]);
        fprintf('\n-- Metric: %s, Sex: %s --\n', current_metric_suffix, current_sex_label);
        fprintf('Paired T-Test (On vs Off) 300Lux: p = %.4f\n', p_300);
        fprintf('Paired T-Test (On vs Off) 1000Lux: p = %.4f\n', p_1000);
        fprintf('Unpaired T-Test (300 vs 1000) On Period: p = %.4f\n', p_on);
        fprintf('Unpaired T-Test (300 vs 1000) Off Period: p = %.4f\n', p_off);
        
        fig_filename = fullfile(savePath, sprintf('Figure5B_Violins_Separate_%s_%s.png', current_sex_label, current_metric_suffix));
        exportgraphics(hFig5B_current, fig_filename, 'Resolution', 300); close(hFig5B_current);
    end
end
disp('--- Finished Figure 5B Components ---');

%% FIG 5C (Original Peak Shift)
disp('--- Figure 5C (Original Peak Shift) has been removed. ---');

%% FIG 5D: REBUILT to show 1000L-300L Difference (Last 2 Weeks)
disp('--- Starting Figure 5D (Rebuilt Difference Plot) ---');
fprintf('\n--- STATS: Figure 5D (Male Change vs Female Change Hourly) ---\n');
for met_idx = 1:length(activity_metrics)
    current_metric_var = activity_metrics{met_idx}; current_metric_suffix = metric_suffixes{met_idx};
    
    hFig5D_current = figure('Name', sprintf('Figure 5D: 1000L-300L Diff by Sex (%s)', current_metric_suffix), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.2 0.2 0.7 0.5]);
    ax_current_5D = axes('Parent', hFig5D_current);
    
    male_300_all = get_last_2_weeks_profile(dataTable, maleAnimals, '300Lux', ztHourVar, current_metric_var);
    male_1000_all = get_last_2_weeks_profile(dataTable, maleAnimals, '1000Lux', ztHourVar, current_metric_var);
    female_300_all = get_last_2_weeks_profile(dataTable, femaleAnimals, '300Lux', ztHourVar, current_metric_var);
    female_1000_all = get_last_2_weeks_profile(dataTable, femaleAnimals, '1000Lux', ztHourVar, current_metric_var);

    male_diffs_individual = male_1000_all - male_300_all;
    female_diffs_individual = female_1000_all - female_300_all;
    
    male_diff_mean = mean(male_diffs_individual, 1, 'omitnan');
    female_diff_mean = mean(female_diffs_individual, 1, 'omitnan');
    
    bar(ax_current_5D, 0:23, male_diff_mean, 0.8, 'FaceColor', maleColor, 'EdgeColor','k', 'DisplayName','Male');
    hold on;
    bar(ax_current_5D, 0:23, female_diff_mean, 0.5, 'FaceColor', femaleColor, 'EdgeColor','k', 'DisplayName','Female');
    plot(xlim, [0 0], 'k--');
    
    fprintf('\n-- Metric: %s --\n', current_metric_suffix);
    p_values_5D = ones(24, 1);
    for zt = 1:24
        [~, p] = ttest2(male_diffs_individual(:, zt), female_diffs_individual(:, zt));
        p_values_5D(zt) = p;
    end
    p_corrected_5D = p_values_5D * 24;
    significant_hours = find(p_corrected_5D < 0.05);
    fprintf('Significant interaction (p < 0.05, Bonferroni) at ZT: %s\n', num2str(significant_hours' - 1));

    y_lims = ylim(ax_current_5D);
    y_range_5D = diff(y_lims);
    for hr = significant_hours'
        zt_hour = hr-1;
        y_max_bar = max(male_diff_mean(hr), female_diff_mean(hr));
        y_pos = y_max_bar + 0.05 * y_range_5D;
        text(ax_current_5D, zt_hour, y_pos, '*', 'FontSize', 16, 'HorizontalAlignment', 'center');
    end

    title(sprintf('Fig 5D: Change in Activity (Last 2wks 1000L - Last 2wks 300L) (%s)', current_metric_suffix));
    ylabel('Change in Activity'); xlabel('ZT Hour');
    legend; grid on; xlim([-0.75, 23.75]);
    
    fig_filename_5D = fullfile(savePath, sprintf('Figure5D_DifferenceBars_Last2wks_%s.png', current_metric_suffix));
    exportgraphics(hFig5D_current, fig_filename_5D, 'Resolution', 300); close(hFig5D_current);
end
disp('--- Finished Figure 5D ---');
%% --- End of Script ---
disp('--- Script for Figure 5 Finished ---');

%% --- Helper Functions ---
function profile_data = get_last_2_weeks_profile(dataTable, animals, cond, ztHourVar, current_metric_var)
    all_animal_profiles = [];
    for an_idx = 1:length(animals)
        animal_data = dataTable(dataTable.Animal == animals(an_idx) & dataTable.Condition == cond, :);
        if isempty(animal_data), continue; end
        unique_days = unique(animal_data.DayOfCondition);
        if isempty(unique_days), continue; end
        last_14_days = unique_days(max(1, end-13):end);
        final_data = animal_data(ismember(animal_data.DayOfCondition, last_14_days), :);
        if isempty(final_data), continue; end
        hourly_means = groupsummary(final_data, ztHourVar, 'mean', current_metric_var);
            
        an_profile = NaN(1, 24);
        [lia, locb] = ismember(0:23, hourly_means.(ztHourVar));
        if any(lia)
            an_profile(lia) = hourly_means.(['mean_', current_metric_var])(locb(lia));
        end
        all_animal_profiles(end+1, :) = an_profile;
    end
    profile_data = all_animal_profiles;
end

function str = pval_to_asterisk(p)
    if p < 0.001, str = '***';
    elseif p < 0.01, str = '**';
    elseif p < 0.05, str = '*';
    else, str = ''; end % Return empty for non-significant
end

function plot_sig_bar(ax, x_coords, p, y_pos)
    ast_str = pval_to_asterisk(p);
    if ~isempty(ast_str)
        plot(ax, x_coords, [y_pos, y_pos], '-k', 'LineWidth', 1.2);
        text(ax, mean(x_coords), y_pos, ast_str, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16);
    end
end