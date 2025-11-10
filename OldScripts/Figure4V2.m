% MATLAB Script for Generating Figure 4 Components as Separate Figures
% Date: Oct 08, 2025
% Author: Noah Muscat & ActigraphyOnly Research Assistant
%
% Description:
% FINAL REVISION: Unabridged script with statistical tests and visual annotations for 4B and 4C.
% Includes fix for graphics rendering error in 4C.
%
% ASSUMPTIONS: violinplot.m is on the MATLAB path.
% ---------------------------------------------------------------------
%% --- Configuration & Setup ---
clear; clc; close all;
% --- File and Path Parameters ---
filename = '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/ActivityAnalysis/ActigraphyOnly/Unified_AO1to12_Combined_Data.csv'; 
savePath = '/Users/noahmuscat/Desktop/WatsonLab/AOPaper/Figure4Plots';
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
% Define condition lists
lightingConditions_Fig4A = {'300Lux', '1000Lux', 'FullDark', '300LuxEnd'};
lightingConditions_Fig4B_Overlay = {'300Lux', '1000Lux', 'FullDark'};
lightingConditions_Fig4C_Violins = {'300Lux', '1000Lux', 'FullDark', '300LuxEnd'};
lightingConditionsFig4_Overall = {'300Lux', '1000Lux', 'FullDark', '300LuxEnd'};
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
dataTable.Condition = categorical(dataTable.Condition, lightingConditionsFig4_Overall, 'Ordinal',true);
dataTable.Sex = categorical(ismember(dataTable.(animalVar), maleAnimals), [0, 1], {'Female', 'Male'});
dataTable.LightPeriod = categorical((dataTable.(ztHourVar) >= ztLightsOnStart & dataTable.(ztHourVar) <= ztLightsOnEnd), [0, 1], {'Off', 'On'});
disp('Data preprocessing complete.');
%% --- Activity Metrics Loop ---
activity_metrics_to_plot = {pixelDiffVar, normalizedActivityVar};
metric_suffixes = {'PixelDiff', 'Normalized'};
metric_ylabels = {'Activity (Pixel Difference)', 'Normalized Activity'};
for met_idx = 1:length(activity_metrics_to_plot)
    current_metric_var = activity_metrics_to_plot{met_idx};
    current_metric_suffix = metric_suffixes{met_idx};
    current_metric_ylabel = metric_ylabels{met_idx};
    fprintf('\n--- Generating Figure 4 components for metric: %s ---\n', current_metric_var);
    
    %% FIG 4A (2 Figures Total, each with 4 Subplots): 48h Diurnality Line Plots by Condition
    disp(['--- Starting Figure 4A - Subplots for 48h Profiles (Metric: ', current_metric_suffix, ') ---']);
    try
        global_min_y = Inf;
        global_max_y = -Inf;
        all_profiles_4A = cell(length(lightingConditions_Fig4A), 1);
        for i_cond_pre = 1:length(lightingConditions_Fig4A)
            cond_name_pre = lightingConditions_Fig4A{i_cond_pre};
            data_this_cond_all_animals_pre = dataTable(ismember(dataTable.(animalVar), fourCondAnimals) & dataTable.Condition == cond_name_pre, :);
            data_for_profile_pre = [];
            unique_animals_in_cond_pre = unique(data_this_cond_all_animals_pre.Animal);
            for i_an_pre = 1:length(unique_animals_in_cond_pre)
                animal_cond_data_pre = data_this_cond_all_animals_pre(data_this_cond_all_animals_pre.Animal == unique_animals_in_cond_pre(i_an_pre), :);
                unique_days_pre = unique(animal_cond_data_pre.DayOfCondition);
                if length(unique_days_pre) >= 7, selected_days_pre = unique_days_pre(end-6:end);
                else, selected_days_pre = unique_days_pre; end
                data_for_profile_pre = [data_for_profile_pre; animal_cond_data_pre(ismember(animal_cond_data_pre.DayOfCondition, selected_days_pre), :)];
            end
            if ~isempty(data_for_profile_pre)
                hourly_mean_pre = groupsummary(data_for_profile_pre, ztHourVar, 'mean', current_metric_var);
                mean_profile_24h_pre = NaN(hoursPerDay, 1);
                [lia_pre, locb_pre] = ismember((0:23)', hourly_mean_pre.(ztHourVar));
                mean_profile_24h_pre(lia_pre) = hourly_mean_pre.(['mean_', current_metric_var])(locb_pre(lia_pre));
                all_profiles_4A{i_cond_pre} = mean_profile_24h_pre;
                global_min_y = min(global_min_y, min(mean_profile_24h_pre));
                global_max_y = max(global_max_y, max(mean_profile_24h_pre));
            end
        end
        y_range = global_max_y - global_min_y;
        y_lims_final = [global_min_y - 0.1*y_range, global_max_y + 0.1*y_range];
        
        hFig4A_Subplots = figure('Name', sprintf('Fig 4A: 48h Profiles by Condition (%s)', current_metric_suffix), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
        for i_cond_4A = 1:length(lightingConditions_Fig4A)
            ax4A = subplot(4, 1, i_cond_4A, 'Parent', hFig4A_Subplots); hold(ax4A, 'on');
            cond_name_4A = lightingConditions_Fig4A{i_cond_4A};
            mean_profile_24h = all_profiles_4A{i_cond_4A};
            if isempty(mean_profile_24h), title(ax4A, sprintf('%s - No Data', cond_name_4A)); hold(ax4A, 'off'); continue; end
            mean_profile_48h = [mean_profile_24h; mean_profile_24h]; x_axis_48h = (0:47)';
            plot(ax4A, x_axis_48h, mean_profile_48h, 'k-', 'LineWidth', 1.5);
            [~, idx_max] = findpeaks(mean_profile_48h);
            if ~isempty(idx_max), plot(ax4A, x_axis_48h(idx_max), mean_profile_48h(idx_max), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'HandleVisibility','off'); end
            ylim(ax4A, y_lims_final);
            yl4A_sub = ylim(ax4A);
            patch(ax4A, [ztLightsOffStart, ztLightsOffEnd+1, ztLightsOffEnd+1, ztLightsOffStart], [yl4A_sub(1) yl4A_sub(1) yl4A_sub(2) yl4A_sub(2)], [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            patch(ax4A, [ztLightsOffStart+24, ztLightsOffEnd+1+24, ztLightsOffEnd+1+24, ztLightsOffStart+24], [yl4A_sub(1) yl4A_sub(1) yl4A_sub(2) yl4A_sub(2)], [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            hold(ax4A, 'off');
            title(ax4A, sprintf('%s', cond_name_4A), 'FontSize', 10);
            ylabel(ax4A, current_metric_ylabel, 'FontSize', 9);
            xticks(ax4A, 0:6:47); xlim(ax4A, [-0.5 47.5]); grid(ax4A, 'on');
            if i_cond_4A < length(lightingConditions_Fig4A), set(ax4A, 'XTickLabel', []);
            else, xlabel(ax4A, 'ZT Hour (Double Plotted)', 'FontSize', 9); end
        end
        sgtitle(hFig4A_Subplots, sprintf('Fig 4A: 48h Profiles by Condition (%s)', current_metric_suffix), 'FontSize', 12, 'FontWeight', 'bold');
        fig_filename = fullfile(savePath, sprintf('Figure4A_SubplotProfiles_SameY_%s.png', current_metric_suffix));
        exportgraphics(hFig4A_Subplots, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]); close(hFig4A_Subplots);
    catch ME_4A, warning('ErrGen:Fig4A_Subplots', 'Error generating Fig 4A Subplots (%s): %s', current_metric_suffix, ME_4A.message); if~isempty(ME_4A.stack),disp(ME_4A.stack(1));end; end
    disp(['--- Finished Figure 4A - Subplot Profiles (Metric: ', current_metric_suffix, ') ---']);
    
    %% FIG 4B: Overlay Activity Profiles with Difference Line (WITH STATS)
    disp(['--- Starting Figure 4B (Metric: ', current_metric_suffix, ') ---']);
    fprintf('\n--- STATS: Figure 4B (300Lux vs 1000Lux Hourly) ---\n');
    try
        hFig4B_OverlayDiff = figure('Name', sprintf('Fig 4B: Overlay Profiles & Difference (%s)', current_metric_suffix), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.2 0.2 0.5 0.6]);
        ax4B_OverlayDiff = axes('Parent', hFig4B_OverlayDiff); hold(ax4B_OverlayDiff, 'on');
        lineColors_4B_main = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880]}; 
        plotHandles_4B_main = gobjects(length(lightingConditions_Fig4B_Overlay)+1,1); mean_profiles_for_diff = table();
        for i_cond_4B = 1:length(lightingConditions_Fig4B_Overlay)
            cond_name_4B = lightingConditions_Fig4B_Overlay{i_cond_4B}; animals_for_this_cond = allAnimals;
            if strcmp(cond_name_4B, 'FullDark'), animals_for_this_cond = fourCondAnimals; end
            data_this_cond_4B = dataTable(ismember(dataTable.(animalVar), animals_for_this_cond) & dataTable.Condition == cond_name_4B, :);
            if isempty(data_this_cond_4B), plotHandles_4B_main(i_cond_4B)=plot(NaN,NaN,'DisplayName',sprintf('%s (No Data)',cond_name_4B)); continue; end
            hourly_stats_4B = groupsummary(data_this_cond_4B, ztHourVar, {'mean','std','nnz'}, current_metric_var);
            hourly_stats_4B = sortrows(hourly_stats_4B, ztHourVar); mean_prof = NaN(hoursPerDay,1); sem_prof = NaN(hoursPerDay,1);
            [lia,locb] = ismember((0:23)', hourly_stats_4B.(ztHourVar)); mean_prof(lia) = hourly_stats_4B.(['mean_',current_metric_var])(locb(lia));
            sem_val = hourly_stats_4B.(['std_',current_metric_var])(locb(lia)) ./ sqrt(hourly_stats_4B.(['nnz_',current_metric_var])(locb(lia)));
            sem_val(hourly_stats_4B.(['nnz_',current_metric_var])(locb(lia)) <=1) = NaN; sem_prof(lia) = sem_val;
            mean_profiles_for_diff.(cond_name_4B) = mean_prof; 
            plotHandles_4B_main(i_cond_4B) = plot(ax4B_OverlayDiff, (0:23)', mean_prof, 'Color', lineColors_4B_main{i_cond_4B}, 'LineWidth', 2, 'DisplayName', cond_name_4B);
            if ~all(isnan(sem_prof)), fill(ax4B_OverlayDiff, [(0:23)'; flipud((0:23)')], [mean_prof-sem_prof; flipud(mean_prof+sem_prof)], lineColors_4B_main{i_cond_4B}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off'); end
        end
        if ismember('300Lux', mean_profiles_for_diff.Properties.VariableNames) && ismember('1000Lux', mean_profiles_for_diff.Properties.VariableNames)
            diff_profile = mean_profiles_for_diff.('1000Lux') - mean_profiles_for_diff.('300Lux');
            plotHandles_4B_main(length(lightingConditions_Fig4B_Overlay)+1) = plot(ax4B_OverlayDiff, (0:23)', diff_profile, 'k:', 'LineWidth', 1.5, 'DisplayName', '1000L-300L Diff');
        end
        
        ylims_4B = ylim(ax4B_OverlayDiff);
        ylim(ax4B_OverlayDiff, [ylims_4B(1), ylims_4B(2) + 0.1*diff(ylims_4B)]);
        yl4B_diff=ylim(ax4B_OverlayDiff); patch(ax4B_OverlayDiff,[ztLightsOffStart,ztLightsOffEnd+1,ztLightsOffEnd+1,ztLightsOffStart],[yl4B_diff(1),yl4B_diff(1),yl4B_diff(2),yl4B_diff(2)],[0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
        uistack(plotHandles_4B_main(isgraphics(plotHandles_4B_main)),'top');
        if strcmp(current_metric_var, normalizedActivityVar), plot(ax4B_OverlayDiff,xlim(ax4B_OverlayDiff),[0 0],'k--','LineWidth',0.75,'HandleVisibility','off'); end
        
        fprintf('\n-- Metric: %s --\n', current_metric_suffix);
        p_values_4B = ones(24, 1);
        for zt = 0:23
            d300 = dataTable(dataTable.Condition=='300Lux' & dataTable.ZT_Time==zt, :).(current_metric_var);
            d1000 = dataTable(dataTable.Condition=='1000Lux' & dataTable.ZT_Time==zt, :).(current_metric_var);
            if ~isempty(d300) && ~isempty(d1000), [~, p] = ttest2(d300, d1000); p_values_4B(zt+1) = p; end
        end
        p_corrected_4B = p_values_4B * 24;
        significant_hours = find(p_corrected_4B < 0.05);
        fprintf('Significant differences (p < 0.05, Bonferroni) at ZT: %s\n', num2str(significant_hours' - 1));

        y_lims_final = ylim(ax4B_OverlayDiff);
        if ~isempty(significant_hours)
            y_pos = y_lims_final(2) - 0.05*diff(y_lims_final);
            block_starts = significant_hours([true; diff(significant_hours)>1]);
            block_ends = significant_hours([diff(significant_hours)>1; true]);
            for i = 1:length(block_starts)
                plot(ax4B_OverlayDiff, [block_starts(i)-1-0.5, block_ends(i)-1+0.5], [y_pos, y_pos], 'k-', 'LineWidth', 2.5);
            end
        end

        hold(ax4B_OverlayDiff, 'off'); title(ax4B_OverlayDiff, sprintf('Fig 4B: Activity Profiles & Difference (%s)', current_metric_suffix), 'FontSize', 11);
        xlabel(ax4B_OverlayDiff, 'ZT Hour', 'FontSize', 10); ylabel(ax4B_OverlayDiff, current_metric_ylabel, 'FontSize', 10);
        xticks(ax4B_OverlayDiff, 0:6:23); xlim(ax4B_OverlayDiff, [-0.5 23.5]); grid(ax4B_OverlayDiff, 'on');
        legend(ax4B_OverlayDiff, plotHandles_4B_main(isgraphics(plotHandles_4B_main)), 'Location', 'northeast');
        fig_filename = fullfile(savePath, sprintf('Figure4B_OverlayAndDiff_%s.png', current_metric_suffix));
        exportgraphics(hFig4B_OverlayDiff, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]); close(hFig4B_OverlayDiff);
    catch ME_4B_Diff, warning('ErrGen:Fig4B_OverlayDiff', 'Error generating Fig 4B (%s): %s', current_metric_suffix, ME_4B_Diff.message); end
    disp(['--- Finished Figure 4B (Metric: ', current_metric_suffix, ') ---']);

    %% FIG 4C (2 Figures Total): Quantification of Diurnality (Violins) - REVISED COLORS & STRUCTURE
disp(['--- Starting Figure 4C Violins (Metric: ', current_metric_suffix, ') ---']);
try
    if ~exist('violinplot.m','file'),error('ViolinPlotNotFound:Fig4C','violinplot.m not found on MATLAB path.');end
    hFig4C_violins = figure('Name', sprintf('Fig 4C: Diurnality Violins by Cond (%s)', current_metric_suffix), ...
                            'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.3 0.3 0.65 0.6]);
    ax4C_violins = axes('Parent', hFig4C_violins);
    hold(ax4C_violins, 'on');
    
    plotData_4C = []; 
    groupings_4C_list_cat = []; 
    groupOrder_4C = {}; 
    violinBodyColors_4C = [];
    x_tick_labels_4C = {}; 
    
    for i_cond_4C = 1:length(lightingConditions_Fig4C_Violins)
        cond_name_4C_str = lightingConditions_Fig4C_Violins{i_cond_4C};
        animals_this_cond = allAnimals;
        if ismember(cond_name_4C_str, ["FullDark", "300LuxEnd"]), animals_this_cond = fourCondAnimals; end
        
        data_cond_all_animals = dataTable(ismember(dataTable.(animalVar), animals_this_cond) & dataTable.Condition == cond_name_4C_str, :);
        if isempty(data_cond_all_animals), continue; end
        
        light_phases = ["On", "Off"];
        for lp_idx = 1:length(light_phases)
            current_lp_str = light_phases(lp_idx);
            
            data_for_violin_body = data_cond_all_animals(data_cond_all_animals.LightPeriod == current_lp_str, :);
            animal_means_for_body = groupsummary(data_for_violin_body, animalVar, 'mean', current_metric_var, 'IncludeMissingGroups', false);
            
            if ~isempty(animal_means_for_body)
                groupName = [cond_name_4C_str, '-', char(current_lp_str)];
                plotData_4C = [plotData_4C; animal_means_for_body.(['mean_',current_metric_var])];
                groupings_4C_list_cat = [groupings_4C_list_cat; repmat(categorical({groupName}), height(animal_means_for_body),1)];
                groupOrder_4C{end+1} = groupName;
                violinBodyColors_4C = [violinBodyColors_4C; grayColor]; 
                x_tick_labels_4C{end+1} = groupName; 
            end
        end
    end
    
    if ~isempty(plotData_4C) && ~isempty(groupOrder_4C)
        groupings_4C_cat_ordered = categorical(groupings_4C_list_cat, groupOrder_4C, 'Ordinal',true);
        
        violinplot(plotData_4C, groupings_4C_cat_ordered, 'Parent', ax4C_violins, ...
                   'ViolinColor', violinBodyColors_4C, 'GroupOrder', groupOrder_4C, ...
                   'ShowData', false, 'ViolinAlpha', 0.5, 'Width', 0.4, ...
                   'EdgeColor',[0 0 0],'BoxColor',[0.1 0.1 0.1], 'MedianColor',[0.95 0.95 0.95]);
        
        current_x_violin_scatter = 0;
        for i_group_order = 1:length(groupOrder_4C)
            current_x_violin_scatter = current_x_violin_scatter + 1;
            group_name_current = groupOrder_4C{i_group_order};
            
            parts = split(group_name_current, '-'); 
            current_scatter_cond_str = string(parts{1});
            current_scatter_lp_str = string(parts{2});
            animals_this_cond_scatter = allAnimals; 
            if ismember(current_scatter_cond_str, ["FullDark", "300LuxEnd"]), animals_this_cond_scatter = fourCondAnimals; end
            
            data_for_scatter_points = dataTable(ismember(dataTable.(animalVar), animals_this_cond_scatter) & ...
                                                dataTable.Condition == current_scatter_cond_str & ...
                                                dataTable.LightPeriod == current_scatter_lp_str, :);
            
            animal_means_for_scatter_plot = groupsummary(data_for_scatter_points, {animalVar,'Sex'}, 'mean', current_metric_var, 'IncludeMissingGroups',false);
            if ~isempty(animal_means_for_scatter_plot)
                jitter_values = (rand(height(animal_means_for_scatter_plot),1)-0.5) * 0.15;
                x_scatter = repmat(current_x_violin_scatter, height(animal_means_for_scatter_plot),1) + jitter_values;
                
                male_points_idx = animal_means_for_scatter_plot.Sex == 'Male';
                female_points_idx = animal_means_for_scatter_plot.Sex == 'Female';
                
                scatter(ax4C_violins, x_scatter(male_points_idx), animal_means_for_scatter_plot.(['mean_',current_metric_var])(male_points_idx), 30, maleColor, 'filled', 'MarkerFaceAlpha',0.7);
                scatter(ax4C_violins, x_scatter(female_points_idx), animal_means_for_scatter_plot.(['mean_',current_metric_var])(female_points_idx), 30, femaleColor, 'filled', 'MarkerFaceAlpha',0.7);
            end
        end
        xticks(ax4C_violins, 1:length(groupOrder_4C)); 
        xticklabels(ax4C_violins, x_tick_labels_4C);
        xtickangle(ax4C_violins, 45);
        xlim(ax4C_violins, [0.5, length(groupOrder_4C) + 0.5]);
    else
         title(ax4C_violins, sprintf('No Data for Violins (%s)', current_metric_suffix));
    end
    hold(ax4C_violins, 'off');
    ylabel(ax4C_violins, ['Mean ', current_metric_ylabel], 'FontSize', 10);
    title(ax4C_violins, sprintf('Fig 4C: Diurnality (On/Off) by Condition (%s)', current_metric_suffix), 'FontSize', 11); grid(ax4C_violins, 'on');
    
    fig_filename = fullfile(savePath, sprintf('Figure4C_Violins_GrayPooled_%s.png', current_metric_suffix));
    exportgraphics(hFig4C_violins, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]);
    close(hFig4C_violins);
catch ME_4C, warning('ErrGen:Fig4C_Violins', 'Error generating Fig 4C Violins (%s): %s', current_metric_suffix, ME_4C.message); if~isempty(ME_4C.stack),disp(ME_4C.stack(1));end; end
disp(['--- Finished Figure 4C Violins (Metric: ', current_metric_suffix, ') ---']);
    
    %% FIG 4D: Shift in Activity Peaks During FullDark
    disp(['--- Starting Figure 4D Peak Shift (Metric: ', current_metric_suffix, ') ---']);
    try
        hFig4D_peak = figure('Name', sprintf('Fig 4D: Peak Shift in FullDark (%s)', current_metric_suffix), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.5 0.5 0.4 0.5]);
        ax4D_peak = axes('Parent', hFig4D_peak);
        data_FD_4D = dataTable(ismember(dataTable.(animalVar), fourCondAnimals) & dataTable.Condition == 'FullDark', :);
        if isempty(data_FD_4D), title(ax4D_peak, sprintf('No FullDark Data (%s)', current_metric_suffix)); close(hFig4D_peak); continue; end
        
        unique_days_FD = unique(data_FD_4D.DayOfCondition); 
        daily_peak_ZT_FD = NaN(length(unique_days_FD),1);
        for i_day_4D = 1:length(unique_days_FD)
            day_val = unique_days_FD(i_day_4D); 
            data_this_day_FD = data_FD_4D(data_FD_4D.DayOfCondition == day_val, :); 
            if isempty(data_this_day_FD), continue; end
            hourly_mean_FD = groupsummary(data_this_day_FD, ztHourVar, 'mean', current_metric_var);
            if ~isempty(hourly_mean_FD)
                hourly_mean_FD = sortrows(hourly_mean_FD, ztHourVar); 
                [~,maxIdx] = max(hourly_mean_FD.(['mean_', current_metric_var])); 
                if ~isempty(maxIdx), daily_peak_ZT_FD(i_day_4D) = hourly_mean_FD.(ztHourVar)(maxIdx(1)); end
            end
        end
        
        plot(ax4D_peak, unique_days_FD(~isnan(daily_peak_ZT_FD)), daily_peak_ZT_FD(~isnan(daily_peak_ZT_FD)), 'o-k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
        title(ax4D_peak, sprintf('Fig 4D: Peak Activity ZT in FullDark (%s)', current_metric_suffix), 'FontSize', 11);
        xlabel(ax4D_peak, 'Day in FullDark', 'FontSize', 10); 
        ylabel(ax4D_peak, 'ZT of Peak Activity', 'FontSize', 10);
        minD_pk=min(unique_days_FD);maxD_pk=max(unique_days_FD); 
        if isempty(minD_pk)||isempty(maxD_pk)||minD_pk==maxD_pk, xlim(ax4D_peak,[0.5,max(1.5,maxD_pk+0.5)]); 
        else, xlim(ax4D_peak,[minD_pk-0.5,maxD_pk+0.5]); end
        yticks(ax4D_peak, 0:6:23); ylim(ax4D_peak, [-1 24]); grid(ax4D_peak, 'on');
        
        fig_filename = fullfile(savePath, sprintf('Figure4D_PeakShift_%s.png', current_metric_suffix));
        exportgraphics(hFig4D_peak, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]); close(hFig4D_peak);
    catch ME_4D, warning('ErrGen:Fig4D_PeakShift', 'Error generating Fig 4D (%s): %s', ME_4D.message); if~isempty(ME_4D.stack),disp(ME_4D.stack(1));end; end
    disp(['--- Finished Figure 4D Peak Shift (Metric: ', current_metric_suffix, ') ---']);
end
disp('--- Script for Figure 4 Components Finished ---');

%% --- Helper Functions ---
function out = iif(condition, trueVal, falseVal)
    if condition, out = trueVal; else out = falseVal; end
end

function str = pval_to_asterisk(p)
    if p < 0.001, str = '***';
    elseif p < 0.01, str = '**';
    elseif p < 0.05, str = '*';
    else, str = ''; end
end

function plot_sig_bar(ax, x_coords, p, y_pos)
    ast_str = pval_to_asterisk(p);
    if ~isempty(ast_str)
        plot(ax, x_coords, [y_pos, y_pos], '-k', 'LineWidth', 1.2);
        text(ax, mean(x_coords), y_pos, ast_str, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16);
    end
end