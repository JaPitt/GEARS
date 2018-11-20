% Copyright 2018 Jake Alan Pitt and Julio R. Banga
 
% This file is part of GEARS.
 
% GEARS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
 
% GEARS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
 
% You should have received a copy of the GNU General Public License
% along with GEARS.  If not, see <http://www.gnu.org/licenses/>.
 
% File author: Jake Alan Pitt (jp00191.su@gmail.com)


function [] = Write_GEARS_results_html_summary(Results, Fitting_data, Save_loc, Validation_data, Figures_folder)
% A function that writes the results of the GEARS procedure into a html file.
% Results         - The GEARS results should be in the format all GEARS results hold. (structure)
% Fitting_data    - The data that has been fit to. Should be in the processed GEARS format. (structure)
% Save_loc        - The name of where the hmtl file will be saved. Please dont use an extension. (string)
% Validation_data - The data used to validate to. Should be in the processed GEARS format. (structure) [Optional]
% Figures_folder  - The folder containing figures to be included in the html report. [Optional]


%% Check inputs

Ins = who;

    if sum(ismember({'Results', 'Fitting_data', 'Save_loc'}, Ins)) ~= 3
        
    error('The first 3 inputs are required to use this function')   
        
    end 
    
    if ismember('Validation_data', Ins) 
        
        if ~isempty(Validation_data)
        
        Val = 1; 

        else

        Val = 0;    

        Validation_Exp_names = [];
    
        end
        
    else

    Val = 0;    

    Validation_Exp_names = [];

    end

    if ismember('Figures_folder', Ins) 
        
    Include_figs = 1;

    else

    Include_figs = 1;

    end

Fitting_params = reshape(Results.Model_information.Fitting_params, [length(Results.Model_information.Fitting_params) 1]);


%% Create folders if they dont exist already

    if exist(Save_loc, 'dir') == 0

    mkdir(Save_loc)

    end

    if Include_figs

        if exist([Save_loc filesep 'Figures'], 'dir') == 0

        mkdir([Save_loc filesep 'Figures'])

        end

    end


%% Prep figures 

Figure_locators = GEARS_prep_figures_for_html(Figures_folder, [Save_loc filesep 'Figures']); 

html_figure_code = @(Figure_name, Figure_caption) cellstr(strcat('<p> <center>', Figure_caption, '<br> <img src="./Figures/', Figure_locators.(Figure_name), '" height="500" width="750">  </br> </center> </pr>'));


%% Copy the GEARS logo to the figures folder in html summary

Logo_path = [fileparts(which('Run_GEARS')) filesep 'Documentation' filesep 'GEARS_logo.svg'];

copyfile(Logo_path, [Save_loc filesep 'Figures' filesep 'GEARS_logo.svg'])

Markup_tag = Results.Markup_tag;

GEARS_logo_path = './Figures/GEARS_logo';


%% Find out which pages are to be used

Results_fields = fieldnames(Results);

    if ismember('Regularised_estimation', fieldnames(Results.Global_parameter_estimation))

    Reg = 1;

    else

    Reg = 0;

    end

    if ismember('Statistics', Results_fields)

    Stats = 1;

    else

    Stats = 0;

    end


%% Create the page links and names for the navigation bar 

Page_names = {'Regularised solution'; 'Statistics for the regularised solution'; 'Non-regularised solution'; 'Statistics for the non-regularised solution'; 'Parameter bounding'; 'Sampling'; 'Warning log'; 'Settings used'};

Page_links = {'Regularised solution.html'; 'Statistics for the regularised solution.html'; 'Non-regularised solution.html'; 'Statistics for the non-regularised solution.html'; 'Parameter bounding.html'; 'Sampling.html'; 'Warning log.html'; 'Settings used.html'};

    if Reg 

    Reg_pages = [0 0 0 0 0 0 0 0]';

    else

    Reg_pages = [1 1 0 0 1 1 0 0]';

    end

    if Stats 

    Stat_pages = [0 0 0 0 0 0 0 0]';

    else

    Stat_pages = [0 1 0 1 0 0 0 0]';

    end

Page_names = Page_names(~(Reg_pages | Stat_pages));

Page_links = Page_links(~(Reg_pages | Stat_pages));


%% Settings used page

% GEARS options 

GEARS_options_settings = fieldnames(Results.Settings_used.GEARS_options);

GEARS_options_value = cell(size(GEARS_options_settings));

    for i = 1:length(GEARS_options_settings)

    GEARS_options_value(i) = Convert_col_cell_str(Results.Settings_used.GEARS_options.(GEARS_options_settings{i}));

    end

GEARS_options_table = table(GEARS_options_settings, GEARS_options_value, 'Variablenames', {'Setting', 'Value'});

GEARS_options_table_html = GEARS_html_table_converter(GEARS_options_table, 'GEARS options');

% AMICI options

AMICI_options_settings = fieldnames(Results.Settings_used.AMICI_options);

AMICI_options_value = cell(size(AMICI_options_settings));

    for i = 1:length(AMICI_options_settings)

    AMICI_options_value(i) = Convert_col_cell_str(Results.Settings_used.AMICI_options.(AMICI_options_settings{i}));

    end

AMICI_options_table = table(AMICI_options_settings, AMICI_options_value, 'Variablenames', {'Setting', 'Value'});

AMICI_options_table_html = GEARS_html_table_converter(AMICI_options_table, 'AMICI options');

% MEIGO options

MEIGO_options_settings = fieldnames(Results.Settings_used.MEIGO_options);

MEIGO_options_settings = MEIGO_options_settings(~ismember(MEIGO_options_settings, {'local', 'log_var'}));

MEIGO_options_value = cell(size(MEIGO_options_settings));

    for i = 1:length(MEIGO_options_settings)

    MEIGO_options_value(i) = Convert_col_cell_str(Results.Settings_used.MEIGO_options.(MEIGO_options_settings{i}));

    end

MEIGO_options_settings_local = fieldnames(Results.Settings_used.MEIGO_options.local);

MEIGO_options_value_local = cell(size(MEIGO_options_settings_local));

    for i = 1:length(MEIGO_options_settings_local)

    MEIGO_options_value_local(i) = Convert_col_cell_str(Results.Settings_used.MEIGO_options.local.(MEIGO_options_settings_local{i}));

    end

MEIGO_options_settings = [MEIGO_options_settings; strcat('local.', MEIGO_options_settings_local)];

MEIGO_options_value = [MEIGO_options_value; MEIGO_options_value_local];

MEIGO_options_table = table(MEIGO_options_settings, MEIGO_options_value, 'Variablenames', {'Setting', 'Value'});

MEIGO_options_table_html = GEARS_html_table_converter(MEIGO_options_table, 'MEIGO options (please note the option log_var has not been included as it is a vector)');

Settings_used_page_content = [GEARS_options_table_html; {'<p></p>'}; AMICI_options_table_html; {'<p></p>'}; MEIGO_options_table_html];

GEARS_html_template(Settings_used_page_content, GEARS_logo_path, Page_names, 'Settings used', Page_links, Markup_tag, Save_loc)


%% Warning log page

Warning_classifications = {'<p style = "padding: 5px"> GEARS warnings are classified as follows:'; '<ul>'; '<li> Red warnings - Warnings that denote problems users should pay direct attention to and may affect the main results of the procedure'; ...
 '<li> Orange warnings - Warnings that denote issues directly affecting some noncrucial parts of the procedure'; ...
'<li> Yellow warnings - Warnings that should be noted but in general should not directly alter the GEARS process.'; ...
'</ul>'; '</p>'};

    if ~isempty([Results.Warning_log.Red; Results.Warning_log.Orange; Results.Warning_log.Orange])

    List_warnings = @(Warning) strcat('<li>', Warning);

        if ~isempty(Results.Warning_log.Red)

        Warnings_red = [{'<p>'; '<div style="background-color: #ff6666; margin-left: 0px; margin-right: 0px; padding-bottom: 8px; padding-left: 25px; padding-right: 10px; padding-top: 10px;">'; ...
        '<h2>Red warnings</h2>'}; List_warnings(Results.Warning_log.Red) ;{'</ol>'; '</div>'; '</p>'}];

        else

        Warnings_red = [];

        end

        if ~isempty(Results.Warning_log.Orange)

        Warnings_orange = [{'<p>'; '<div style="background-color: #ffbf80; margin-left: 0px; margin-right: 0px; padding-bottom: 10px; padding-left: 25px; padding-right: 10px; padding-top: 10px;">'; ...
        '<h2>Orange warnings</h2>'}; List_warnings(Results.Warning_log.Orange) ;{'</ol>'; '</div>'; '</p>'}];

        else

        Warnings_orange = [];

        end

        if ~isempty(Results.Warning_log.Yellow)

        Warnings_yellow = [{'<p>'; '<div style="background-color: #ffff99;; margin-left: 0px; margin-right: 0px; padding-bottom: 10px; padding-left: 25px; padding-right: 10px; padding-top: 10px;">'; ...
        '<h2>Yellow warnings</h2>'}; List_warnings(Results.Warning_log.Yellow) ;{'</ol>'; '</div>'; '</p>'}];

        else

        Warnings_yellow = [];

        end

    Warnings_page_content = [Warning_classifications; Warnings_red; Warnings_orange; Warnings_yellow];

    else

    Warnings_page_content = {'<p> Congratulations, there are no warnings to report.</p>'};

    end

GEARS_html_template(Warnings_page_content, GEARS_logo_path, Page_names, 'Warning log', Page_links, Markup_tag, Save_loc)


%% Sampling page

    if Reg

    Sample_table = table(Results.Sampling.Sample_size/prod(Results.Parameter_bounding.Original_bounds.Upper - Results.Parameter_bounding.Original_bounds.Lower), ...
    median(Results.Sampling.Cost_cut_offs), Results.Sampling.Sample_size, ...
    'variablenames', {'Sample_dens_approx', 'Avg_Cost_cut_off', 'Sample_size'});

        if Results.Settings_used.GEARS_options.Plot_samples

        Sampling_contents = [GEARS_html_table_converter(Sample_table, 'A Sampling summary.'); html_figure_code('Sample_plot', 'The parameter-cost sampling distribution.')];

        else

        Sampling_contents = GEARS_html_table_converter(Sample_table, 'A Sampling summary.');

        end

    GEARS_html_template(Sampling_contents, GEARS_logo_path, Page_names, 'Sampling', Page_links, Markup_tag, Save_loc)

    end

%% Parameter bounding page

    if Reg

    Parameter_bounds_original_table = table(Results.Parameter_bounding.Original_bounds.Lower', Results.Parameter_bounding.Original_bounds.Upper', 'variablenames', {'Lower', 'Upper'}, 'rownames', Fitting_params);

    Parameter_bounds_reduced_table = table(Results.Parameter_bounding.Reduced_bounds.Lower', Results.Parameter_bounding.Reduced_bounds.Upper', 'variablenames', {'Lower', 'Upper'}, 'rownames', Fitting_params);

        if Results.Settings_used.GEARS_options.Plot_bounds_param_confidence

        Bounding_contents = [GEARS_html_table_converter(Parameter_bounds_original_table, 'The original parameter bounds.'); {'<p></p>'}; GEARS_html_table_converter(Parameter_bounds_reduced_table, 'The reduced parameter bounds.'); ...
        {'<p></p>'}; html_figure_code('Bounds_param_confidence_plot', 'Plot showing the reduced parameter bounds alongside the estimated solution.')];

        else

        Bounding_contents = [GEARS_html_table_converter(Parameter_bounds_original_table, 'The original parameter bounds.'); {'<p></p>'}; GEARS_html_table_converter(Parameter_bounds_reduced_table, 'The reduced parameter bounds.')];        

        end

    GEARS_html_template(Bounding_contents, GEARS_logo_path, Page_names, 'Parameter bounding', Page_links, Markup_tag, Save_loc)

    end


%% Regularised Statistics page

Reg_stats_contents = [];

    if Reg && Stats

        if Results.Settings_used.GEARS_options.Calculate_NRMSE

        Fitting = [Results.Statistics.NRMSE.Fitting.Regularised_estimation.All_experiments Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.All_experiments];

        Fitting_Exp_names = fieldnames(Fitting_data);

            if length(Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise) > 1

            Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

            Fitting = [Fitting; Results.Statistics.NRMSE.Fitting.Regularised_estimation.Experiment_wise  Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise];

            end

        Fitting_NRMSE_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

        Fitting_NRMSE_table_html = GEARS_html_table_converter(Fitting_NRMSE_table, 'A table showing the NRMSE for each solution in fitting.');

        Reg_stats_contents = [Reg_stats_contents; Fitting_NRMSE_table_html; {'<p></p>' }];

            if Val

            Validation = [Results.Statistics.NRMSE.Validation.Regularised_estimation.All_experiments Results.Statistics.NRMSE.Validation.Non_regularised_estimation.All_experiments];

            Validation_Exp_names = fieldnames(Validation_data);

                if length(Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise) > 1

                Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                Validation = [Validation; Results.Statistics.NRMSE.Validation.Regularised_estimation.Experiment_wise  Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise];

                end

            Validation_NRMSE_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

            Validation_NRMSE_table_html = GEARS_html_table_converter(Validation_NRMSE_table, 'A table showing the NRMSE for each solution in cross-validation.');

            Reg_stats_contents = [Reg_stats_contents; Validation_NRMSE_table_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Plot_residuals

        Regularised_fitting_residuals_html = html_figure_code('Regularised_fitting_residuals', 'A plot showing the residuals for the regularised solution in fitting.');

        Reg_stats_contents = [Reg_stats_contents; Regularised_fitting_residuals_html; {'<p></p>' }];

            if Val

            Regularised_validation_residuals_html = html_figure_code('Regularised_validation_residuals', 'A plot showing the residuals for the regularised solution in cross-validation.');

            Reg_stats_contents = [Reg_stats_contents; Regularised_validation_residuals_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Calculate_R2

        Fitting = [Results.Statistics.R2.Fitting.Regularised_estimation.All_experiments Results.Statistics.R2.Fitting.Non_regularised_estimation.All_experiments];

        Fitting_Exp_names = fieldnames(Fitting_data);

            if length(Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise) > 1

            Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

            Fitting = [Fitting; Results.Statistics.R2.Fitting.Regularised_estimation.Experiment_wise  Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise];

            end

        Fitting_R2_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

        Fitting_R2_table_html = GEARS_html_table_converter(Fitting_R2_table, 'A table showing the R2 for each solution in fitting.');

        Reg_stats_contents = [Reg_stats_contents; Fitting_R2_table_html; {'<p></p>' }];

            if Val

            Validation = [Results.Statistics.R2.Validation.Regularised_estimation.All_experiments Results.Statistics.R2.Validation.Non_regularised_estimation.All_experiments];

            Validation_Exp_names = fieldnames(Validation_data);

                if length(Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise) > 1

                Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                Validation = [Validation; Results.Statistics.R2.Validation.Regularised_estimation.Experiment_wise  Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise];

                end

            Validation_R2_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

            Validation_R2_table_html = GEARS_html_table_converter(Validation_R2_table, 'A table showing the R2 for each solution in cross-validation.');

            Reg_stats_contents = [Reg_stats_contents; Validation_R2_table_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Plot_pred_vs_meas

        Regularised_fitting_pred_vs_meas_html = html_figure_code('Regularised_fitting_pred_vs_meas', 'A plot of residuals against measurements for the regularised solution in fitting.');

        Reg_stats_contents = [Reg_stats_contents; Regularised_fitting_pred_vs_meas_html; {'<p></p>' }];

            if Val

            Regularised_validation_pred_vs_meas_html = html_figure_code('Regularised_validation_pred_vs_meas', 'A plot of residuals against measurements for the regularised solution in cross-validation.');

            Reg_stats_contents = [Reg_stats_contents; Regularised_validation_pred_vs_meas_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Plot_correlation_matrix

        Regularised_fitting_corr_matrix_html = html_figure_code('Regularised_fitting_corr_matrix', 'A heatmap of the correlation matrix for the regularised solution.');

        Reg_stats_contents = [Reg_stats_contents; Regularised_fitting_corr_matrix_html; {'<p></p>' }];

        end

        if Results.Settings_used.GEARS_options.Calculate_chi2

        Fitting = [{num2str(Results.Statistics.chi2.Fitting.Regularised_estimation.Probability)} {num2str(Results.Statistics.chi2.Fitting.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Fitting.Regularised_estimation.Conclusion}, {Results.Statistics.chi2.Fitting.Non_regularised_estimation.Conclusion}];

        Fitting_Exp_names = {'p'; 'Conclusion'};

        Fitting_chi2_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

        Fitting_chi2_table_html = GEARS_html_table_converter(Fitting_chi2_table, 'A table showing the chi2 for each solution in fitting.');

        Reg_stats_contents = [Reg_stats_contents; Fitting_chi2_table_html; {'<p></p>' }];

            if Val

            Validation = [{num2str(Results.Statistics.chi2.Validation.Regularised_estimation.Probability)} {num2str(Results.Statistics.chi2.Validation.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Validation.Regularised_estimation.Conclusion} {Results.Statistics.chi2.Validation.Non_regularised_estimation.Conclusion}];

            Validation_Exp_names = {'p'; 'Conclusion'};

            Validation_chi2_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

            Validation_chi2_table_html = GEARS_html_table_converter(Validation_chi2_table, 'A table showing the chi2 for each solution in cross-validation.');

            Reg_stats_contents = [Reg_stats_contents; Validation_chi2_table_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Plot_convergence_curves
    
        Convergence_html = html_figure_code('Convergence_curves', 'A plot showing the convergence curves for the regularised estimation.');

        Reg_stats_contents = [Reg_stats_contents; Convergence_html; {'<p></p>' }];

        end

    GEARS_html_template(Reg_stats_contents, GEARS_logo_path, Page_names, 'Statistics for the regularised solution', Page_links, Markup_tag, Save_loc)

    end


%% Non-regularised Statistics page

Non_reg_stats_contents = [];

    if Stats

        if Results.Settings_used.GEARS_options.Calculate_NRMSE

            if Reg

            Non_reg_stats_contents = [Non_reg_stats_contents; Fitting_NRMSE_table_html; {'<p></p>' }];

                if Val

                Non_reg_stats_contents = [Non_reg_stats_contents; Validation_NRMSE_table_html; {'<p></p>' }];

                end

            else 

            Fitting = Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.All_experiments;

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_NRMSE_table = table(Fitting(:, 1), 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            Fitting_NRMSE_table_html = GEARS_html_table_converter(Fitting_NRMSE_table, 'A table showing the NRMSE for the non-regularised solution in fitting.');

            Non_reg_stats_contents = [Non_reg_stats_contents; Fitting_NRMSE_table_html; {'<p></p>' }];

                if Val

                Validation = Results.Statistics.NRMSE.Validation.Non_regularised_estimation.All_experiments;

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation;  Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_NRMSE_table = table(Validation(:, 1), 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                Validation_NRMSE_table_html = GEARS_html_table_converter(Validation_NRMSE_table, 'A table showing the NRMSE for the non-regularised solution in cross-validation.');

                Non_reg_stats_contents = [Non_reg_stats_contents; Validation_NRMSE_table_html; {'<p></p>' }];
             
                end

            end

        end

        if Results.Settings_used.GEARS_options.Plot_residuals

        Non_regularised_fitting_residuals_html = html_figure_code('Non_regularised_fitting_residuals', 'A plot showing the residuals for the non-regularised solution in fitting.');

        Non_reg_stats_contents = [Non_reg_stats_contents; Non_regularised_fitting_residuals_html; {'<p></p>' }];

            if Val

            Non_regularised_validation_residuals_html = html_figure_code('Non_regularised_validation_residuals', 'A plot showing the residuals for the non-regularised solution in cross-validation.');

            Non_reg_stats_contents = [Non_reg_stats_contents; Non_regularised_validation_residuals_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Calculate_R2

            if Reg

            Non_reg_stats_contents = [Non_reg_stats_contents; Fitting_R2_table_html; {'<p></p>' }];

                if Val

                Non_reg_stats_contents = [Non_reg_stats_contents; Validation_R2_table_html; {'<p></p>' }];

                end

            else

            Fitting = Results.Statistics.R2.Fitting.Non_regularised_estimation.All_experiments;

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting;  Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_R2_table = table(Fitting(:, 1), 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            Fitting_R2_table_html = GEARS_html_table_converter(Fitting_R2_table, 'A table showing the R2 for the non-regularised solution in fitting.');

            Non_reg_stats_contents = [Non_reg_stats_contents; Fitting_R2_table_html; {'<p></p>' }];

                if Val

                Validation = Results.Statistics.R2.Validation.Non_regularised_estimation.All_experiments;

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation; Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_R2_table = table(Validation(:, 1), 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                Validation_R2_table_html = GEARS_html_table_converter(Validation_R2_table, 'A table showing the R2 for the non-regularised solution in cross-validation.');

                Non_reg_stats_contents = [Reg_stats_contents; Validation_R2_table_html; {'<p></p>' }];

                end

            end

        end

        if Results.Settings_used.GEARS_options.Plot_pred_vs_meas

        Non_regularised_fitting_pred_vs_meas_html = html_figure_code('Non_regularised_fitting_pred_vs_meas', 'A plot of residuals against measurements for the non-regularised solution in fitting.');

        Non_reg_stats_contents = [Reg_stats_contents; Non_regularised_fitting_pred_vs_meas_html; {'<p></p>' }];

            if Val

            Non_regularised_validation_pred_vs_meas_html = html_figure_code('Non_regularised_validation_pred_vs_meas', 'A plot of residuals against measurements for the non-regularised solution in cross-validation.');

            Non_reg_stats_contents = [Reg_stats_contents; Non_regularised_validation_pred_vs_meas_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Plot_correlation_matrix

        Non_regularised_fitting_corr_matrix_html = html_figure_code('Non_regularised_fitting_corr_matrix', 'A heatmap of the correlation matrix for the non-regularised solution.');

        Non_reg_stats_contents = [Reg_stats_contents; Non_regularised_fitting_corr_matrix_html; {'<p></p>' }];

        end

        if Results.Settings_used.GEARS_options.Calculate_chi2

            if Reg

            Non_reg_stats_contents = [Non_reg_stats_contents; Fitting_chi2_table_html; {'<p></p>' }];

                if Val

                Non_reg_stats_contents = [Non_reg_stats_contents; Validation_chi2_table_html; {'<p></p>' }];

                end

            else 

            Fitting = [{num2str(Results.Statistics.chi2.Fitting.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Fitting.Non_regularised_estimation.Conclusion}];

            Fitting_Exp_names = {'p'; 'Conclusion'};

            Fitting_chi2_table = table(Fitting(:, 1),  'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            Fitting_chi2_table_html = GEARS_html_table_converter(Fitting_chi2_table, 'A table showing the chi2 for the non-regularised solution in fitting.');

            Non_reg_stats_contents = [Non_reg_stats_contents; Fitting_chi2_table_html; {'<p></p>' }];

                if Val

                Validation = [{num2str(Results.Statistics.chi2.Validation.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Validation.Non_regularised_estimation.Conclusion}];

                Validation_Exp_names = {'p'; 'Conclusion'};

                Validation_chi2_table = table(Validation(:, 1), 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                Validation_chi2_table_html = GEARS_html_table_converter(Validation_chi2_table, 'A table showing the chi2 for the non-regularised solution in cross-validation.');

                Non_reg_stats_contents = [Non_reg_stats_contents; Validation_chi2_table_html; {'<p></p>' }];

                end

            end

        end

    GEARS_html_template(Non_reg_stats_contents, GEARS_logo_path, Page_names, 'Statistics for the non-regularised solution', Page_links, Markup_tag, Save_loc)

    end


%% Regularised solution page

    if Reg

    Regularised_solution_contents = [];

    Headers = {'Parameter_value', 'Confidence_95', 'Coeff_of_var_percentage', 'Bounds_status', 'P_ref'};

    Reg_param_summary = table(Results.Global_parameter_estimation.Regularised_estimation.xbest', Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(: ,3), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 4), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5), Results.Regularisation.P_ref, 'variablenames', Headers, 'rownames', Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 1));  

    Reg_param_summary_html = GEARS_html_table_converter(Reg_param_summary, ['A summary of key results from the regularised estimation (with alpha = ' num2str(Results.Regularisation.alpha) ').']);

    Regularised_solution_contents = [Regularised_solution_contents; Reg_param_summary_html; {'<p></p>' }];

        if Results.Settings_used.GEARS_options.Plot_trajectories_with_uncertainty

        Regularised_fit_with_uncertainty_html = html_figure_code('Regularised_fit_with_uncertainty', 'A plot showing the regularised fit with uncertanty.');

        Regularised_solution_contents = [Regularised_solution_contents; Regularised_fit_with_uncertainty_html; {'<p></p>' }];

            if Val

            Regularised_validation_with_uncertainty_html = html_figure_code('Regularised_validation_with_uncertainty', 'A plot showing the regularised cross-validation with uncertanty.');

            Regularised_solution_contents = [Regularised_solution_contents; Regularised_validation_with_uncertainty_html; {'<p></p>' }];

            end

        end

        if Results.Settings_used.GEARS_options.Plot_trajectories

        Fitting_figures_html = html_figure_code('Fitting_figures', 'A plot showing the regularised and non-regularised fits.');

        Regularised_solution_contents = [Regularised_solution_contents; Fitting_figures_html; {'<p></p>' }];

            if Val

            Validation_figures_html = html_figure_code('Validation_figures', 'A plot showing the regularised and non-regularised solutions in cross-validation.');

            Regularised_solution_contents = [Regularised_solution_contents; Validation_figures_html; {'<p></p>' }];

            end

        end

    GEARS_html_template(Regularised_solution_contents, GEARS_logo_path, Page_names, 'Regularised solution', Page_links, Markup_tag, Save_loc)

    end


%% Non-regularised solution page

Non_regularised_solution_contents = [];

Headers = {'Parameter_value', 'Confidence_95', 'Coeff_of_var_percentage', 'Bounds_status'};

Non_reg_param_summary = table(Results.Global_parameter_estimation.Non_regularised_estimation.xbest', Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(: ,3), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 4), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5), 'variablenames', Headers, 'rownames', Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 1));  

Non_reg_param_summary_html = GEARS_html_table_converter(Non_reg_param_summary, 'A summary of key results from the non-regularised estimation.');

Non_regularised_solution_contents = [Non_regularised_solution_contents; Non_reg_param_summary_html; {'<p></p>' }];

    if Results.Settings_used.GEARS_options.Plot_trajectories_with_uncertainty

    Non_regularised_fit_with_uncertainty_html = html_figure_code('Non_regularised_fit_with_uncertainty', 'A plot showing the non-regularised fit with uncertanty.');

    Non_regularised_solution_contents = [Non_regularised_solution_contents; Non_regularised_fit_with_uncertainty_html; {'<p></p>' }];

        if Val

        Non_regularised_validation_with_uncertainty_html = html_figure_code('Non_regularised_validation_with_uncertainty', 'A plot showing the non-regularised cross-validation with uncertanty.');
        
        Non_regularised_solution_contents = [Non_regularised_solution_contents; Non_regularised_validation_with_uncertainty_html; {'<p></p>' }];

        end

    end

    if Results.Settings_used.GEARS_options.Plot_trajectories

        if Reg

        Fitting_figures_html = html_figure_code('Fitting_figures', 'A plot showing the regularised and non-regularised fits.');

        else

        Fitting_figures_html = html_figure_code('Fitting_figures', 'A plot showing the non-regularised fit.');

        end

        Non_regularised_solution_contents = [Non_regularised_solution_contents; Fitting_figures_html; {'<p></p>' }];

        if Val

            if Reg

            Validation_figures_html = html_figure_code('Validation_figures', 'A plot showing the regularised and non-regularised solutions in cross-validation.');

            else

            Validation_figures_html = html_figure_code('Validation_figures', 'A plot showing the non-regularised solution in cross-validation.');

            end

        Non_regularised_solution_contents = [Non_regularised_solution_contents; Validation_figures_html; {'<p></p>' }];

        end

    end

GEARS_html_template(Non_regularised_solution_contents, GEARS_logo_path, Page_names, 'Non-regularised solution', Page_links, Markup_tag, Save_loc)


end


function [] = GEARS_html_template(Content, GEARS_logo_path, Page_names, This_page_name, Page_links, Markup_tag, Folder_save_path)
% A function that creates html pages in the GEARS format.
% Content - The content to be on the page. (Cell array).
% GEARS_logo_path - The page to the GEARS logo, typically in the documentation folder (string).
% Page_names - The list of page names to be included in the navigation bar. (cell array)
% This_page_name - The name of this specific page. (string)
% Page_links - The links to be used as hyperlink in the navigation bar, will be assigned in the same order as Page_names (cell array)
% Markup_tag - The markup identification tag used within GEARS. (string)
% Folder_save_path - The path to the folder that page should be saved in. (string)


%% Check inputs

Inputs = who;

    if sum(ismember({'Content', 'GEARS_logo_path', 'Page_names', 'This_page_name', 'Page_links', 'Markup_tag', 'Folder_save_path'}, Inputs)) ~= length(Inputs) 

    error('All inputs are required to run this function')

    end

Content = reshape(Content, length(Content), 1);

Page_names = reshape(Page_names, length(Page_names), 1);

Page_links = reshape(Page_links, length(Page_links), 1);



%% Format inputs

Title = [This_page_name ' for GEARS run: ' Markup_tag];

Nav_item = @(Name, Links) strcat('<li class="head"><A style="text-decoration: none" href="', Links, '">', Name, '</A></li>');

Nav_item_active = @(Name, Links) strcat('<li class="active"><A style="text-decoration: none" href="', Links, '">', Name, '</A></li>');

This_page_pos = find(ismember(Page_names, This_page_name));

    if This_page_pos == 1

    Nav_bar = [Nav_item_active(Page_names(1), Page_links(1));  Nav_item(Page_names(2:end), Page_links(2:end))];

    else

        if This_page_pos == length(Page_names)

        Nav_bar = [Nav_item(Page_names(1:end-1), Page_links(1:end-1)); Nav_item_active(Page_names(end), Page_links(end));  ];

        else

        Nav_bar = [Nav_item(Page_names(1:This_page_pos-1), Page_links(1:This_page_pos-1)); Nav_item_active(Page_names(This_page_pos), Page_links(This_page_pos)); Nav_item(Page_names(This_page_pos+1:end), Page_links(This_page_pos+1:end))];

        end

    end


%% Write the page into the template

HTML_page = [{'<!DOCTYPE html>'; '<html>'; '<head>'; '<style>'; 'body {margin:0;}'; 'ul.head {list-style-type: none;'; 'margin: 0;'; ...
'padding: 0;'; 'overflow: hidden;'; 'background-color: black;'; 'position: -webkit-sticky; /* Safari */'; 'position: sticky;'; ...
'top: 0;'; 'top: 0;'; 'width: 100%;'; 'margin-bottom:2px'; 'text-decoration: none}'; 'li.head {	border-right: 1px solid white;'; ...
'float: left;'; 'display: block;'; 'color: #EBF1DE;'; 'text-align: center;'; 'padding: 10px 10px;}'; '.active {'; 'background-color: #4CAF50;}'; ...
'li.active {'; 'border-right: 1px solid white;'; 'float: left;'; 'display: block;'; 'color: black;'; 'text-align: center;'; 'padding: 10px 10px;'; ...
'background-color: #EBF1DE;}'; 'li.head a:hover {background-color: #0099ff;'; 'color = black;}'; 'li.active a:hover {background-color: #0099ff;'; ...
'color = black;}'; 'li.head a:visited{color:#EBF1DE}'; 'li.head a:link{color:#EBF1DE}'; 'li.active a:link{color:black}'; 'li.active a:visited{color:black;'; 'background-color: #EBF1DE;}'; 'h1.head {display: block;'; ...
'font-size: 2em;'; 'margin-top: 0em;'; 'margin-bottom: 0em;'; 'margin-left: 0;'; 'margin-right: 0;'; 'font-weight: bold;'; 'vertical-align:middle;'; ...
'line-height: 4;';'padding:5px}'; '</style>'; '</head>'; '<body>'; '<font size="4">'; ...
 ['<h1 class="head"> <span> ' Title ' </span> </center> <img src="' GEARS_logo_path '.svg" alt="GEARS_logo" height="150" width="175" align="right"></h1> ']}; ...          % Title
{'<font size="4">'; '<ul class="head">'}; ... 
Nav_bar;
 {'</ul>'; '</font>';'<p></p>'}; 
Content; ...                                                                                                                                             % Actual content
{'</html>'; '</font>'; '</body>'}];


%% Write the html file

Z  = fopen([Folder_save_path filesep This_page_name '.html'], 'wt');

    for i = 1:length(HTML_page)
        
    fprintf(Z, '%s\n', [char(HTML_page(i))]);
                
    end

fclose(Z);

end


function [Figure_locators] = GEARS_prep_figures_for_html(Figures_folder, Figures_save_folder)
% A function that creates svg copies of all the figures saved by GEARS for use in the html summary.
% Figures_folder - The folder where the figures are saved (string).
% Figures_save_folder - The path to where the figures should be saved.

Dir_info = dir(Figures_folder);

Is_dir = [Dir_info.isdir];  % Find the index for directories
  
Figures_to_include = strcat(Figures_folder ,filesep ,{Dir_info(~Is_dir).name}');

File_size_mb = [Dir_info(~Is_dir).bytes]'/2^20; % Convert to mbs

Large_figures = File_size_mb > 1; % Figures larger than 1mb will be in the html in png format rather than svg format.

Format = repmat({'.svg'}, size(Large_figures));

Format(Large_figures) = {'.png'};

    for i = 1:length(Figures_to_include)

    Handle = openfig(Figures_to_include{i}, 'visible', 'invisible');

        if length(Handle) > 1

            for j = 1:length(Handle)

            [~, Name ,~] = fileparts(Figures_to_include{i});

            Sub_handle = Handle(2);

            saveas(Sub_handle, [Figures_save_folder filesep Name '_' num2str(j) Format{i}])  

            end

        else

        j = 0;

        [~, Name ,~] = fileparts(Figures_to_include{i});        

        saveas(Handle, [Figures_save_folder filesep Name Format{i}])        

        end

    close(Handle)

        if j == 0

        Figure_locators.(Name) = strcat(Name, Format{i});

        else

        Figure_locators.(Name) = strcat(strrep(cellstr(strcat(Name, '_', num2str((1:j)'))), ' ', ''), Format{i});

        end

    end

end


function [HTML_table] = GEARS_html_table_converter(Table, Caption)
% A function that writes a HTML table code containing the table in Table.
% Table - The table to be written in HTML (Table)
% Caption - A caption for the table (String)[Optional]

%% Check inputs 

Inputs = who;

    if isempty(Table)

    error('A table is required')

    end

    if ~ismember('Caption', Inputs)

    Caption = '';

    end
    


%% Extract names from table

Table_col_headers = Table.Properties.VariableNames;

Table_row_headers = Table.Properties.RowNames;


%% Find table size

Num_cols = size(Table, 2);

Num_rows = size(Table, 1);


%% Extract data from table

Table_cols = cell(size(Table_col_headers));

    for i = 1:Num_cols

    Table_cols{i} = Convert_col_cell_str(Table.(Table_col_headers{i}));

    end

% Jake insert a conversion to str here. We need all of them to be strs

Table_info = cell(size(Table));

    for i = 1:Num_cols

    Table_info(:, i) = Table_cols{i};

    end


%% HTML code definitions

HTML_header = @(x) strcat('<th class="GEARS_table">', x , '</th>');

HTML_data = @(x) strcat('<td class="GEARS_table">', x , '</td>');

HTML_row = @(x) ['<tr class="GEARS_table">'; x; '</tr>'];


%% Write HTML table

HTML_table = cell(10^5, 1); % Initialise an empty table

HTML_table(1:9) = {'<style>'; 'table.GEARS_table{border: 2px solid black; border-collapse: collapse; margin-left:10px; padding: 3px; width: 99%}'; 'th.GEARS_table{border: 2px solid black; text-align: left;  background-color: #84A246; color: black;}'; ...
'td.GEARS_table{border: 1px solid black; text-align: left}'; 'tr.GEARS_table:nth-child(even) {background-color: #EBF1DE;}'; 'tr.GEARS_table:hover {background-color: #000000; color: white}'; '</style>'; ...
'<table class="GEARS_table">'; ['<caption>' Caption '</caption>']};

    if isempty(Table_row_headers)

    Current_addition = HTML_row(HTML_header(Table_col_headers'));

    else

    Current_addition = HTML_row(HTML_header([{''} ; Table_col_headers']));

    end

HTML_table(10: 10 + length(Current_addition) - 1) = Current_addition;

Current_position = 10 + length(Current_addition);

    for i = 1:Num_rows

        if isempty(Table_row_headers)

        Current_addition = HTML_row(HTML_data(Table_info(i, :)'));

        else

        Current_addition = HTML_row([HTML_header(Table_row_headers{i}); HTML_data(Table_info(i, :)')]);

        end

    HTML_table(Current_position: Current_position + length(Current_addition) - 1) = Current_addition;

    Current_position = Current_position + length(Current_addition) - 1;

    end

HTML_table = [HTML_table(1:Current_position); '</table>'];

end


function [Cell_str_col] = Convert_col_cell_str(Col)
% A function that converts table columns into cell strings
% Col - The table column

Type = class(Col);

sprintf_v = @(x) cellfun(@(y) sprintf('%e', y), reshape(num2cell(x), [length(x), 1]), 'uni', 0);

Options = {'double'; 'logical'; 'cell'; 'char'};

Class_found = ismember(Options, Type);

    if sum(Class_found) == 0 

    error('The format type was not found')

    end

    if find(Class_found) == 1

    Cell_str_col = sprintf_v(Col);

    end

    if find(Class_found) == 2 

    Cell_str_col = cellstr(num2str(Col));

    end

    if find(Class_found) == 3 

    Not_numbers = cellfun(@isempty, cellfun(@str2num, Col, 'uni', 0));
        
    Cell_str_col = Col;

        if sum(Not_numbers) ~= length(Not_numbers)

        Cell_str_nums = Cell_str_col(~Not_numbers);

        Num_values = str2double(Cell_str_nums);

        Not_logical = Num_values ~= 1 & Num_values ~= 0;

        Cell_str_nums(Not_logical) = sprintf_v((str2num(cell2mat(Cell_str_col(Not_logical)))));

        Cell_str_col(~Not_numbers) = Cell_str_nums;

        Includes_plus_minus = ~cellfun(@isempty, (cellfun(@(x) strfind(x, '+-'), Col, 'uni', 0)));

        Num_without_pm = cellstr(num2str(cellfun(@(x) abs(str2num(x)), Col(Includes_plus_minus))));

        Num_without_pm = sprintf_v((str2num(cell2mat(Num_without_pm))));

        Num_with_pm = strcat('+-', Num_without_pm);

        Cell_str_col(Includes_plus_minus) = Num_with_pm;

        end

    end

    if find(Class_found) == 4

        if sum(size(Col) == 1) > 0

            if isempty(str2num(Col))

            Cell_str_col = {Col};

            else

            Cell_str_col = sprintf_v(str2num(Col));

            Includes_plus_minus = ~cellfun(@isempty, (cellfun(@(x) strfind(x, '+-'), Col, 'uni', 0)));

            Num_without_pm = cellstr(num2str(cellfun(@(x) abs(str2num(x)), Col(Includes_plus_minus))));

            Num_without_pm = sprintf_v((str2num(cell2mat(Num_without_pm))));

            Num_with_pm = strcat('+-', Num_without_pm);

            Cell_str_col(Includes_plus_minus) = Num_with_pm;

            end       
            
        else

        Cell_str_col = cellstr(Col);

            if isempty(str2num(Col{1}))

            Cell_str_col = Col;

            else

            Cell_str_col = sprintf_v(str2num(cell2mat(Col)));

            Includes_plus_minus = ~cellfun(@isempty, (cellfun(@(x) strfind(x, '+-'), Col, 'uni', 0)));

            Num_without_pm = cellstr(num2str(cellfun(@(x) abs(str2num(x)), Col(Includes_plus_minus))));

            Num_without_pm = sprintf_v((str2num(cell2mat(Num_without_pm))));

            Num_with_pm = strcat('+-', Num_without_pm);

            Cell_str_col(Includes_plus_minus) = Num_with_pm;

            end

        end

    end

end

