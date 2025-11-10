% MATLAB Script for Generating Figure 3 Components as Separate Figures
% Date: Oct 08, 2025
% Author: Noah Muscat & ActigraphyOnly Research Assistant
%
% Description:
% FINAL REVISION: Added 'drawnow' command to fix saving error in 3D.
%
% ASSUMPTIONS: violinplot.m is on the MATLAB path.
% ---------------------------------------------------------------------
%% --- Configuration & Setup ---
clear; clc; close all; % Start fresh
% --- File and Path Parameters ---
filename = '/Users/noahmuscat/University of Michigan Dropbox/Noah Muscat/ActivityAnalysis/ActigraphyOnly/Unified_AO1to12_Combined_Data.csv'; 
savePath = '/Users/noahmuscat/Desktop/'; 
% --- Key Data Column Names ---
pixelDiffVar = 'SelectedPixelDifference';
normalizedActivityVar = 'NormalizedActivity';
animalVar = 'Animal';
conditionVar = 'Condition';
relativeDayVar = 'RelativeDay';
timeVarZT = 'DateZT';
ztHourVar = 'ZT_Time'; % Integer ZT hour of the day (0-23)
% --- Analysis Parameters ---
samplingIntervalMinutes = 5;
hoursPerDay = 24;
daysForPoolingFig3AB = 14; 
targetConditionFig3AB = '300Lux'; 
% ZT definitions
ztLightsOnStart = 0; ztLightsOnEnd = 11;
ztLightsOffStart = 12; ztLightsOffEnd = 23;
% --- Animal Grouping ---
allAnimals = string({'AO1', 'AO2', 'AO3', 'AO4', 'AO5', 'AO6', 'AO7', 'AO8', 'AO9', 'AO10', 'AO11', 'AO12'});
maleAnimals = string({'AO1', 'AO2', 'AO3', 'AO7', 'AO9', 'AO12'});
femaleAnimals = string({'AO4', 'AO5', 'AO6', 'AO8', 'AO10', 'AO11'});
sexColorMap = containers.Map();
sexColorMap('Male') = [0.2 0.5470 0.8410]; 
sexColorMap('Female') = [0.9500 0.4250 0.1980];
grayColor = [0.6 0.6 0.6];
%% --- Setup Save Directory ---
disp(['Figures will be saved to: ', savePath]);
if ~isfolder(savePath), try mkdir(savePath); disp('Save directory created.'); catch ME, error('SaveDirErr:CouldNotCreate', 'Could not create save directory "%s". Err: %s', savePath, ME.message); end, end
%% --- Load Data ---
disp(['Loading data from: ', filename]);
try
    opts = detectImportOptions(filename);
    opts.VariableNamingRule = 'preserve';
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
for c_idx = 1:length(colsToClean)
    colData = dataTable.(colsToClean{c_idx});
    if isdatetime(colData), nanRowsFilter = nanRowsFilter | isnat(colData);
    elseif isnumeric(colData), nanRowsFilter = nanRowsFilter | isnan(colData);
    elseif isstring(colData), nanRowsFilter = nanRowsFilter | ismissing(colData);
    end
end
if any(nanRowsFilter), fprintf('%d rows removed due to missing values.\n', sum(nanRowsFilter)); dataTable(nanRowsFilter, :) = []; end
if isempty(dataTable), error('PreprocErr:NoData', 'No valid data after NaN removal.'); end
dataTable = dataTable(ismember(dataTable.(animalVar), allAnimals), :);
if isempty(dataTable), error('PreprocErr:NoAnimalData', 'No data found for the animals specified in the `allAnimals` list.'); end 
dataTable.Condition = categorical(dataTable.Condition); 
dataTable.Sex = categorical(ismember(dataTable.(animalVar), maleAnimals), [0, 1], {'Female', 'Male'});
disp('Data preprocessing complete.');
%% --- Activity Metrics Loop ---
activity_metrics_to_plot_fig3 = {pixelDiffVar, normalizedActivityVar};
metric_suffixes_fig3 = {'PixelDiff', 'Normalized'};
metric_ylabels_fig3 = {'Activity (Pixel Difference)', 'Normalized Activity'};

% Create a single baseline data table for all of Figure 3
initialBaselineData = [];
for i_an = 1:length(allAnimals)
    currentAnimalID = allAnimals(i_an);
    animalData = dataTable(dataTable.(animalVar) == currentAnimalID, :);
    animalConds = unique(animalData.Condition, 'stable');
    if ~isempty(animalConds) && animalConds(1) == string(targetConditionFig3AB)
        animalTargetCondData = animalData(animalData.Condition == string(targetConditionFig3AB), :);
        if ~isempty(animalTargetCondData)
            uniqueDays = unique(animalTargetCondData.DayOfCondition);
            if isempty(uniqueDays), continue; end
            daysToKeep = uniqueDays(max(1, end-daysForPoolingFig3AB+1):end);
            finalDataForPooling = animalTargetCondData(ismember(animalTargetCondData.DayOfCondition, daysToKeep), :);
            initialBaselineData = [initialBaselineData; finalDataForPooling];
        end
    end
end
if isempty(initialBaselineData), error('DataErr:Fig3_NoBaseline', 'No baseline data found for Figure 3.'); end

for met_idx_fig3 = 1:length(activity_metrics_to_plot_fig3)
    current_metric_var_fig3 = activity_metrics_to_plot_fig3{met_idx_fig3};
    current_metric_suffix_fig3 = metric_suffixes_fig3{met_idx_fig3};
    current_metric_ylabel_fig3 = metric_ylabels_fig3{met_idx_fig3};
    fprintf('\n--- Generating Figure 3 components for metric: %s ---\n', current_metric_var_fig3);
    
    %% FIG 3A: Averaged Actogram
    disp(['--- Starting Figure 3A: Averaged Actogram (Metric: ', current_metric_suffix_fig3, ') ---']);
    try
        alignedAnimalData = NaN(daysForPoolingFig3AB, hoursPerDay, length(allAnimals));
        for i_an = 1:length(allAnimals)
            currentAnimalID = allAnimals(i_an);
            animalData = dataTable(dataTable.(animalVar) == currentAnimalID, :);
            animalConds = unique(animalData.Condition, 'stable');
            if ~isempty(animalConds) && animalConds(1) == string(targetConditionFig3AB)
                animalTargetCondData = animalData(animalData.Condition == string(targetConditionFig3AB), :);
                if isempty(animalTargetCondData), continue; end
                uniqueDays = unique(animalTargetCondData.DayOfCondition);
                daysToKeep = uniqueDays(max(1, end-daysForPoolingFig3AB+1):end);
                numDaysToProcess = min(length(daysToKeep), daysForPoolingFig3AB);
                for i_day_idx = 1:numDaysToProcess
                    dayNum = daysToKeep(i_day_idx);
                    dayData = animalTargetCondData(animalTargetCondData.DayOfCondition == dayNum, :);
                    if ~isempty(dayData)
                        hourlyMeans = groupsummary(dayData, ztHourVar, 'mean', current_metric_var_fig3);
                        [lia, locb] = ismember(hourlyMeans.(ztHourVar), (0:23)');
                        alignedAnimalData(i_day_idx, locb(lia), i_an) = hourlyMeans.(['mean_', current_metric_var_fig3])(lia);
                    end
                end
            end
        end
        actogramMatrix = mean(alignedAnimalData, 3, 'omitnan');
        
        hFig3A_current = figure('Name', sprintf('Fig 3A: Averaged Actogram (%s)', current_metric_suffix_fig3), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.5 0.6]);
        ax3A = axes('Parent', hFig3A_current); cla(ax3A);
        imagesc(ax3A, (0:23), (1:daysForPoolingFig3AB), actogramMatrix); 
        colormap(ax3A, 'parula'); ax3A.YDir = 'reverse';
        xlabel(ax3A, 'ZT Hour','FontSize',10); xticks(ax3A, [0,6,12,18,23]);
        ylabel(ax3A, sprintf('Day of %s (Aligned)',targetConditionFig3AB),'FontSize',10); 
        ylim(ax3A,[0.5, daysForPoolingFig3AB + 0.5]);
        yticks(ax3A, 2:2:daysForPoolingFig3AB);
        title(ax3A,sprintf('Averaged Actogram (%s)\nLast %d days of %s',strrep(current_metric_suffix_fig3,'_',' '),daysForPoolingFig3AB,targetConditionFig3AB),'FontSize',11);
        cb3A = colorbar(ax3A); ylabel(cb3A, current_metric_ylabel_fig3,'FontSize',10);
        hold(ax3A,'on'); 
        patch(ax3A,[ztLightsOffStart-0.5, ztLightsOffEnd+0.5, ztLightsOffEnd+0.5, ztLightsOffStart-0.5],[0.5,0.5,daysForPoolingFig3AB+0.5,daysForPoolingFig3AB+0.5],[0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off'); 
        hold(ax3A,'off');
        set(ax3A,'Layer','top');
        fig_filename = fullfile(savePath, sprintf('Figure3A_Actogram_Averaged_%s.png', current_metric_suffix_fig3));
        exportgraphics(hFig3A_current, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]);
        close(hFig3A_current);
    catch ME_3A, warning('ErrGen:Fig3A','Error Fig 3A (%s): %s', current_metric_suffix_fig3, ME_3A.message); if~isempty(ME_3A.stack),disp(ME_3A.stack(1));end; end
    disp('--- Finished Figure 3A ---');
    
    %% FIG 3B: Pooled Periodogram
    disp(['--- Starting Figure 3B: Pooled Periodogram (Metric: ', current_metric_suffix_fig3, ') ---']);
    try
        hFig3B_current = figure('Name', sprintf('Fig 3B: Periodogram (%s)', current_metric_suffix_fig3), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.15 0.15 0.5 0.5]);
        ax3B = axes('Parent', hFig3B_current); cla(ax3B);
        tempDataForPerio = initialBaselineData;
        tempDataForPerio.RoundedTime5Min = dateshift(tempDataForPerio.(timeVarZT),'start','minute',samplingIntervalMinutes);
        pooled5MinStats = groupsummary(tempDataForPerio,'RoundedTime5Min','mean',current_metric_var_fig3);
        pooled5MinStats = sortrows(pooled5MinStats,'RoundedTime5Min');
        activitySignal = pooled5MinStats.(['mean_',current_metric_var_fig3]);
        timeSignalHours = hours(pooled5MinStats.RoundedTime5Min - pooled5MinStats.RoundedTime5Min(1));
        validIdx = ~isnan(activitySignal) & ~isnan(timeSignalHours);
        activitySignal = activitySignal(validIdx); timeSignalHours = timeSignalHours(validIdx);
        
        periodogramRangeHours = [0.5, 48]; oversamplingFactor = 4;
        if length(activitySignal) < (hoursPerDay * (60/samplingIntervalMinutes) / 2), error('DataErr:Fig3B_Short','Insufficient data for periodogram.'); end
        fmax_for_plomb = 1 / periodogramRangeHours(1);
        [pxx, f_hz] = plomb(activitySignal, timeSignalHours, fmax_for_plomb, oversamplingFactor);
        periodHours = 1 ./ f_hz;
        validDataIdx = isfinite(periodHours) & isfinite(pxx) & ~isnan(periodHours) & ~isnan(pxx);
        periodHours = periodHours(validDataIdx); pxx = pxx(validDataIdx);
        displayRangeIdx = periodHours >= periodogramRangeHours(1) & periodHours <= periodogramRangeHours(2);
        periodHours = periodHours(displayRangeIdx); pxx = pxx(displayRangeIdx);
        [periodHours, sortIdx] = sort(periodHours); pxx = pxx(sortIdx);
        
        plot(ax3B,periodHours,pxx,'LineWidth',1.5,'Color','k');
        xlabel(ax3B,'Period (hours)','FontSize',10); ylabel(ax3B,'Power','FontSize',10);
        title(ax3B,sprintf('Periodogram (%s)\nLast %dd %s',current_metric_suffix_fig3,daysForPoolingFig3AB,targetConditionFig3AB),'FontSize',11);
        grid(ax3B,'on'); xlim(ax3B,periodogramRangeHours); set(ax3B,'FontSize',9);
        
        fig_filename = fullfile(savePath, sprintf('Figure3B_Periodogram_Last14d_%s.png', current_metric_suffix_fig3));
        exportgraphics(hFig3B_current, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]);
        close(hFig3B_current);
    catch ME_3B, warning('ErrGen:Fig3B','Error Fig 3B (%s): %s',current_metric_suffix_fig3, ME_3B.message); if~isempty(ME_3B.stack),disp(ME_3B.stack(1));end; end
    disp('--- Finished Figure 3B ---');
    
    %% FIG 3C: Pooled 24-Hour Activity Profile
    disp(['--- Starting Figure 3C: Pooled 24-Hour Profile (Metric: ', current_metric_suffix_fig3, ') ---']);
    try
        hFig3C_current = figure('Name', sprintf('Fig 3C: 24h Profile (%s)', current_metric_suffix_fig3), 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.2 0.2 0.5 0.5]);
        ax3C = axes('Parent', hFig3C_current); cla(ax3C);
        
        hourlyStats_3C = groupsummary(initialBaselineData, ztHourVar, {'mean','std','nnz'}, current_metric_var_fig3);
        hourlyStats_3C = sortrows(hourlyStats_3C, ztHourVar);
        meanAct_3C = hourlyStats_3C.(['mean_',current_metric_var_fig3]);
        stdAct_3C  = hourlyStats_3C.(['std_',current_metric_var_fig3]);
        countAct_3C= hourlyStats_3C.(['nnz_',current_metric_var_fig3]);
        semAct_3C  = stdAct_3C ./ sqrt(countAct_3C); semAct_3C(countAct_3C <= 1) = NaN;
        ztH_3C = hourlyStats_3C.(ztHourVar);
        
        errorbar(ax3C,ztH_3C,meanAct_3C,semAct_3C,'o-','LineWidth',1.5,'MarkerSize',4,'CapSize',3,'Color',[0.2 0.2 0.2]);
        hold(ax3C,'on');
        yl3C=ylim(ax3C); patch(ax3C,[ztLightsOffStart-0.5, ztLightsOffEnd+0.5, ztLightsOffEnd+0.5, ztLightsOffStart-0.5], [yl3C(1) yl3C(1) yl3C(2) yl3C(2)], [0.85 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
        if strcmp(current_metric_var_fig3, normalizedActivityVar), plot(ax3C,xlim(ax3C),[0 0],'k--','LineWidth',0.75,'HandleVisibility','off'); end
        hold(ax3C,'off');
        
        xlabel(ax3C,'ZT Hour','FontSize',10); ylabel(ax3C,['Mean ',current_metric_ylabel_fig3],'FontSize',10);
        title(ax3C,sprintf('Pooled 24h Activity Profile (%s - Baseline)',current_metric_suffix_fig3),'FontSize',11);
        xticks(ax3C,0:3:23); xlim(ax3C,[-0.5,23.5]); grid(ax3C,'on'); set(ax3C,'FontSize',9);
        
        fig_filename = fullfile(savePath, sprintf('Figure3C_Profile_Baseline_%s.png', current_metric_suffix_fig3));
        exportgraphics(hFig3C_current, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]);
        close(hFig3C_current);
    catch ME_3C, warning('ErrGen:Fig3C','Error Fig 3C (%s): %s',current_metric_suffix_fig3, ME_3C.message); if~isempty(ME_3C.stack),disp(ME_3C.stack(1));end; end
    disp('--- Finished Figure 3C ---');
    
%% FIG 3D: Diurnality Quantification (Violin Plots) (WITH STATS)
    disp(['--- Starting Figure 3D (Metric: ', current_metric_suffix_fig3, ') ---']);
    fprintf('\n--- STATS: Figure 3D (Bright vs Dark) ---\n');
    try
        if ~exist('violinplot.m','file'),error('ViolinPlotNotFound:Fig3D','violinplot.m not found.');end
        
        hFig3D_current = figure('Name', sprintf('Fig 3D: Diurnality Violins (%s)', current_metric_suffix_fig3), ...
            'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.25 0.25 0.7 0.6], 'Visible', 'off'); 
        
        ax3D = axes('Parent', hFig3D_current); cla(ax3D); hold(ax3D, 'on');
        
        time_windows = {[0 23], [0 11], [3 9], [12 23], [15 21], [10 14], [22 2]};
        window_labels = {'Full Day', 'Bright (0-11)', 'Mid-Bright (3-9)', 'Dark (12-23)', 'Mid-Dark (15-21)', 'Bright-Dark (10-14)', 'Dark-Bright (22-2)'};
        animal_window_means_list = [];
        
        % --- MODIFICATION START ---
        % Use the main 'dataTable' which contains all conditions, instead of 'initialBaselineData'.
        dataForFig3D = dataTable;
        allAnimalsInDataset = unique(dataForFig3D.Animal);
        % --- MODIFICATION END ---
        
        for i_an = 1:length(allAnimalsInDataset)
            currentAnimalID = allAnimalsInDataset(i_an);
            % --- MODIFICATION: Use the full dataset for this animal ---
            T_animal = dataForFig3D(dataForFig3D.Animal == currentAnimalID, :); 
            if isempty(T_animal), continue; end
            currentSex = T_animal.Sex(1);
            for i_win = 1:length(time_windows)
                window = time_windows{i_win}; windowLabel = window_labels{i_win};
                if window(1) <= window(2), idx_time = T_animal.(ztHourVar) >= window(1) & T_animal.(ztHourVar) <= window(2);
                else, idx_time = T_animal.(ztHourVar) >= window(1) | T_animal.(ztHourVar) <= window(2); end
                mean_act_this_win = mean(T_animal.(current_metric_var_fig3)(idx_time), 'omitnan');
                if ~isnan(mean_act_this_win), animal_window_means_list = [animal_window_means_list; {currentAnimalID, windowLabel, currentSex, mean_act_this_win}]; end
            end
        end
        if isempty(animal_window_means_list), error('DataErr:Fig3D_NoData', 'No data calculated for Fig 3D.'); end
        
        plotTable_3D = cell2table(animal_window_means_list, 'VariableNames', {'Animal', 'Window', 'Sex', 'MeanActivity'});
        plotTable_3D.Window = categorical(plotTable_3D.Window, window_labels);
        
        violinplot(plotTable_3D.MeanActivity, plotTable_3D.Window, 'Parent', ax3D, 'GroupOrder', window_labels, 'ShowData', false, 'ViolinColor', grayColor, 'ViolinAlpha', 0.4);
        
        jitterAmount = 0.15;
        unique_windows = get(ax3D, 'XTickLabel');
        for i_win = 1:length(unique_windows)
            current_window_label = unique_windows{i_win};
            window_data = plotTable_3D(plotTable_3D.Window == current_window_label, :);
            x_base = i_win; 
            x_scatter = x_base + (rand(height(window_data), 1) - 0.5) * jitterAmount;
            male_idx = window_data.Sex == 'Male';
            female_idx = window_data.Sex == 'Female';
            scatter(ax3D, x_scatter(male_idx), window_data.MeanActivity(male_idx), 40, sexColorMap('Male'), 'filled', 'MarkerFaceAlpha', 0.9, 'LineWidth', 0.75, 'MarkerEdgeColor', 'k');
            scatter(ax3D, x_scatter(female_idx), window_data.MeanActivity(female_idx), 40, sexColorMap('Female'), 'filled', 'MarkerFaceAlpha', 0.9, 'LineWidth', 0.75, 'MarkerEdgeColor', 'k');
        end
        
        bright_data = plotTable_3D(plotTable_3D.Window == 'Bright (0-11)', :);
        dark_data = plotTable_3D(plotTable_3D.Window == 'Dark (12-23)', :);
        stats_table = outerjoin(bright_data(:,{'Animal', 'MeanActivity'}), dark_data(:,{'Animal', 'MeanActivity'}), 'Keys','Animal', 'MergeKeys', true);
        stats_table.Properties.VariableNames(2:3) = {'Bright', 'Dark'};
        [~, p] = ttest(stats_table.Bright, stats_table.Dark);
        fprintf('Paired T-Test (Bright vs Dark): p = %.4f\n', p);
        y_max = max(plotTable_3D.MeanActivity, [], 'all', 'omitnan');
        y_range = diff(ylim(ax3D));
        if y_range==0, y_range=y_max; end
        plot_sig_bar(ax3D, [2, 4], p, y_max + 0.1*y_range);
        ylim([min(0, min(plotTable_3D.MeanActivity,[],'all','omitnan')-0.1*y_range) y_max + 0.2*y_range]);
        
        xticks(ax3D, 1:length(window_labels)); xticklabels(ax3D, window_labels); xtickangle(ax3D, 45);
        xlim(ax3D, [0.5, length(window_labels) + 0.5]);
        ylabel(ax3D,['Mean ',current_metric_ylabel_fig3],'FontSize',10);
        % --- MODIFICATION: Updated title to reflect the change ---
        title(ax3D,sprintf('Diurnality Quantifications (%s) - All Conditions',current_metric_suffix_fig3),'FontSize',11);
        grid(ax3D, 'off'); set(ax3D,'FontSize',9); 
        if strcmp(current_metric_var_fig3, normalizedActivityVar), line(ax3D, xlim(ax3D), [0 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.75, 'HandleVisibility', 'off'); end
        hold(ax3D, 'off');
        
        hFig3D_current.Visible = 'on';
        drawnow;
        
        fig_filename = fullfile(savePath, sprintf('Figure3D_Violins_PooledSex_Stats_%s_AllConditions.png', current_metric_suffix_fig3));
        exportgraphics(hFig3D_current, fig_filename, 'Resolution', 300); disp(['Saved: ', fig_filename]);
        close(hFig3D_current);
    catch ME_3D, warning('ErrGen:Fig3D','Error Fig 3D (%s): %s', ME_3D.message); if~isempty(ME_3D.stack),disp(ME_3D.stack(1));end; end
    disp('--- Finished Figure 3D ---');
end
disp('--- Script for Figure 3 Components Finished ---');

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