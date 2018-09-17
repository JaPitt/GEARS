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


function [] = Write_GEARS_results_xls_summary(Results, Fitting_data, Save_loc, Validation_data)
% A function that writes the results of the GEARS procedure into an xls spreadsheet. Formating will be performed if excel is installed.
% Results         - The GEARS results should be in the format all GEARS results hold. (structure)
% Fitting_data    - The data that has been fit to. Should be in the processed GEARS format. (structure)
% Save_loc        - The name of where the excel file will be saved. Please dont use an extension. (string)
% Validation_data - The data used to validate to. Should be in the processed GEARS format. (structure) [Optional]


%% Check inputs
warning('OFF', 'MATLAB:xlswrite:AddSheet') % To stop matlab outputting warnings for adding sheets (is turned back on at the end of the function)

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

Fitting_params = reshape(Results.Model_information.Fitting_params, [length(Results.Model_information.Fitting_params) 1]);

filename = [Save_loc '.xls'];

    if exist(filename, 'file') ~= 0
        
    delete(filename) % To overwrite
                
    end

   
%% Global parameter estimation

    if ismember('Regularised_estimation', fieldnames(Results.Global_parameter_estimation))

    Reg = 1;

    else

    Reg = 0;

    end

a = table({Results.Markup_tag});

    try

    writetable(a, filename, 'Sheet', 'Global parameter estimation', 'Range', 'A1', 'WriteRowNames', true, 'WriteVariablenames', false)

    Headers = {'Parameter_value', 'Confidence_95', 'Coeff_of_var_percentage', 'Bounds_status'};
    
    catch
        
    warning('There was an error exporting to xls, the xls writeup will be skipped.')    
       
    return 
        
    end

Non_reg_param_summary = table(Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 2), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(: ,3), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 4), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5), 'variablenames', Headers, 'rownames', Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 1));  

Headers =[Headers 'P_ref'];

Non_reg_param_summary.Properties.DimensionNames{1} = 'Parameter_names';

    if Reg 

    Reg_param_summary = table(Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 2), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(: ,3), Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 4), Results.Global_parameter_estimation.Non_regularised_estimation.Parameter_summary(:, 5), Results.Regularisation.P_ref, 'variablenames', Headers, 'rownames', Results.Global_parameter_estimation.Regularised_estimation.Parameter_summary(:, 1));    

    Reg_param_summary.Properties.DimensionNames{1} = 'Parameter_names';

    a = table([], 'Variablename', {'Regularised_parameter_estimation_results'});

    writetable(a, filename, 'Sheet', 'Global parameter estimation', 'Range', 'C1', 'WriteRowNames', true)

    writetable(Reg_param_summary, filename, 'Sheet', 'Global parameter estimation', 'Range', 'A3', 'WriteRowNames', true)

    a = table({['alpha: ' num2str(Results.Regularisation.alpha)]});

    writetable(a, filename, 'Sheet', 'Global parameter estimation', 'Range', 'F2', 'WriteRowNames', true, 'WriteVariablenames', false)

    a = table([], 'Variablename', {'Non_regularised_parameter_estimation_results'});

    writetable(a, filename, 'Sheet', 'Global parameter estimation', 'Range', 'K1', 'WriteRowNames', true)

    writetable(Non_reg_param_summary, filename, 'Sheet', 'Global parameter estimation', 'Range', 'I3', 'WriteRowNames', true)

    else

    a = table([], 'Variablename', {'Non_regularised_parameter_estimation_results'});

    writetable(a, filename, 'Sheet', 'Global parameter estimation', 'Range', 'C1', 'WriteRowNames', true)

    writetable(Non_reg_param_summary, filename, 'Sheet', 'Global parameter estimation', 'Range', 'A3', 'WriteRowNames', true)

    end


%% Statistics

Results_fields = fieldnames(Results);

    if ismember('Statistics', Results_fields)
        
    a = table({Results.Markup_tag});

    writetable(a, filename, 'Sheet', 'Statistics', 'Range', 'A1', 'WriteRowNames', true, 'WriteVariablenames', false)    
        
    Row_index = 1;    
                
        if Results.Settings_used.GEARS_options.Calculate_NRMSE

            if Reg
            
            Fitting = [Results.Statistics.NRMSE.Fitting.Regularised_estimation.All_experiments Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.All_experiments];

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.NRMSE.Fitting.Regularised_estimation.Experiment_wise  Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_NRMSE_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            else

            Fitting = Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.All_experiments;

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.NRMSE.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_NRMSE_table = table(Fitting, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            end
        
            if Val

                if Reg
            
                Validation = [Results.Statistics.NRMSE.Validation.Regularised_estimation.All_experiments Results.Statistics.NRMSE.Validation.Non_regularised_estimation.All_experiments];

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation; Results.Statistics.NRMSE.Validation.Regularised_estimation.Experiment_wise  Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_NRMSE_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                Validation_NRMSE_table.Properties.DimensionNames{1} = 'Experiment';

                a = table([], 'Variablename', {'NRMSEs_for_the_cross_validation'});

                writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['G' num2str(Row_index)], 'WriteRowNames', true)

                writetable(Validation_NRMSE_table, filename, 'Sheet', 'Statistics', 'Range', ['F' num2str(Row_index + 2)], 'WriteRowNames', true)   
 
                else

                Validation = Results.Statistics.NRMSE.Validation.Non_regularised_estimation.All_experiments;

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation; Results.Statistics.NRMSE.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_NRMSE_table = table(Validation, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                Validation_NRMSE_table.Properties.DimensionNames{1} = 'Experiment';

                a = table([], 'Variablename', {'NRMSEs_for_the_cross_validation'});

                writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['F' num2str(Row_index)], 'WriteRowNames', true)

                writetable(Validation_NRMSE_table, filename, 'Sheet', 'Statistics', 'Range', ['E' num2str(Row_index + 2)], 'WriteRowNames', true)
                
                end

            end
                
        Fitting_NRMSE_table.Properties.DimensionNames{1} = 'Experiment';
                
        a = table([], 'Variablename', {'NRMSEs_for_the_fitting'});
            
        writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['B' num2str(Row_index)], 'WriteRowNames', true)
             
        writetable(Fitting_NRMSE_table, filename, 'Sheet', 'Statistics', 'Range', ['A' num2str(Row_index + 2)], 'WriteRowNames', true)   
        
        NRMSE_start_row = num2str(Row_index + 2);
        
        Row_index = Row_index + max(length(Validation_Exp_names), length(Fitting_Exp_names)) + 6;
          
        end
        
        if Results.Settings_used.GEARS_options.Calculate_R2

            if Reg
            
            Fitting = [Results.Statistics.R2.Fitting.Regularised_estimation.All_experiments Results.Statistics.R2.Fitting.Non_regularised_estimation.All_experiments];

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.R2.Fitting.Regularised_estimation.Experiment_wise  Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_R2_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            else

            Fitting = Results.Statistics.R2.Fitting.Non_regularised_estimation.All_experiments;

            Fitting_Exp_names = fieldnames(Fitting_data);

                if length(Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise) > 1

                Fitting_Exp_names = ['All experiments' ;Fitting_Exp_names];

                Fitting = [Fitting; Results.Statistics.R2.Fitting.Non_regularised_estimation.Experiment_wise];

                end

            Fitting_R2_table = table(Fitting, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);

            end
        
            if Val

                if Reg
            
                Validation = [Results.Statistics.R2.Validation.Regularised_estimation.All_experiments Results.Statistics.R2.Validation.Non_regularised_estimation.All_experiments];

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation; Results.Statistics.R2.Validation.Regularised_estimation.Experiment_wise  Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_R2_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                Validation_R2_table.Properties.DimensionNames{1} = 'Experiment';

                a = table([], 'Variablename', {'R2s_for_the_cross_validation'});
                 
                writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['G' num2str(Row_index)], 'WriteRowNames', true)

                writetable(Validation_R2_table, filename, 'Sheet', 'Statistics', 'Range', ['F' num2str(Row_index + 2)], 'WriteRowNames', true) 
   
                else

                Validation = Results.Statistics.R2.Validation.Non_regularised_estimation.All_experiments;

                Validation_Exp_names = fieldnames(Validation_data);

                    if length(Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise) > 1

                    Validation_Exp_names = ['All experiments' ;Validation_Exp_names];

                    Validation = [Validation; Results.Statistics.R2.Validation.Non_regularised_estimation.Experiment_wise];

                    end

                Validation_R2_table = table(Validation, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 

                Validation_R2_table.Properties.DimensionNames{1} = 'Experiment';

                a = table([], 'Variablename', {'R2s_for_the_cross_validation'});
                
                writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['F' num2str(Row_index)], 'WriteRowNames', true)

                writetable(Validation_R2_table, filename, 'Sheet', 'Statistics', 'Range', ['E' num2str(Row_index + 2)], 'WriteRowNames', true) 

                end
                    
            end
                
        a = table([], 'Variablename', {'R2s_for_the_fitting_procedure'});
            
        writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['B' num2str(Row_index)], 'WriteRowNames', true)
        
        Fitting_R2_table.Properties.DimensionNames{1} = 'Experiment';
            
        writetable(Fitting_R2_table, filename, 'Sheet', 'Statistics', 'Range', ['A' num2str(Row_index + 2)], 'WriteRowNames', true)  
        
        R2_start_row = num2str(Row_index + 2);
        
        Row_index = Row_index + max(length(Validation_Exp_names), length(Fitting_Exp_names)) + 6;
          
        end
        
        if Results.Settings_used.GEARS_options.Calculate_chi2

            if Reg
            
            Fitting = [{num2str(Results.Statistics.chi2.Fitting.Regularised_estimation.Probability)} {num2str(Results.Statistics.chi2.Fitting.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Fitting.Regularised_estimation.Conclusion}, {Results.Statistics.chi2.Fitting.Non_regularised_estimation.Conclusion}];

            Fitting_Exp_names = {'p'; 'Conclusion'};

            Fitting_chi2_table = table(Fitting(:, 1), Fitting(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);
            
            Fitting_chi2_table.Properties.DimensionNames{1} = '.';

                if Val

                Validation = [{num2str(Results.Statistics.chi2.Validation.Regularised_estimation.Probability)} {num2str(Results.Statistics.chi2.Validation.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Validation.Regularised_estimation.Conclusion} {Results.Statistics.chi2.Validation.Non_regularised_estimation.Conclusion}];

                Validation_Exp_names = {'p'; 'Conclusion'};

                Validation_chi2_table = table(Validation(:, 1), Validation(:, 2), 'variablenames', {'Regularised_estimation', 'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 
            
                Validation_chi2_table.Properties.DimensionNames{1} = '.';

                a = table([], 'Variablename', {'chi2s_for_the_cross_validation'});

                writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['G' num2str(Row_index)], 'WriteRowNames', true)

                writetable(Validation_chi2_table, filename, 'Sheet', 'Statistics', 'Range', ['F' num2str(Row_index + 2)], 'WriteRowNames', true)    

                end

            a = table([], 'Variablename', {'chi2s_for_the_fitting_procedure'});

            writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['B' num2str(Row_index)], 'WriteRowNames', true)

            writetable(Fitting_chi2_table, filename, 'Sheet', 'Statistics', 'Range', ['A' num2str(Row_index + 2)], 'WriteRowNames', true)  

            chi2_start_row = num2str(Row_index + 2);

            else

            Fitting = [{num2str(Results.Statistics.chi2.Fitting.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Fitting.Non_regularised_estimation.Conclusion}];

            Fitting_Exp_names = {'p'; 'Conclusion'};

            Fitting_chi2_table = table(Fitting, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Fitting_Exp_names);
            
            Fitting_chi2_table.Properties.DimensionNames{1} = '.';

                if Val

                Validation = [{num2str(Results.Statistics.chi2.Validation.Non_regularised_estimation.Probability)}; {Results.Statistics.chi2.Validation.Non_regularised_estimation.Conclusion}];

                Validation_Exp_names = {'p'; 'Conclusion'};

                Validation_chi2_table = table(Validation, 'variablenames', {'Non_regularised_estimation'}, 'rownames', Validation_Exp_names); 
                
                Validation_chi2_table.Properties.DimensionNames{1} = '.';

                a = table([], 'Variablename', {'chi2s_for_the_cross_validation'});

                writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['F' num2str(Row_index)], 'WriteRowNames', true)

                writetable(Validation_chi2_table, filename, 'Sheet', 'Statistics', 'Range', ['E' num2str(Row_index + 2)], 'WriteRowNames', true)    

                end

            a = table([], 'Variablename', {'chi2s_for_the_fitting_procedure'});

            writetable(a, filename, 'Sheet', 'Statistics', 'Range', ['B' num2str(Row_index)], 'WriteRowNames', true)

            writetable(Fitting_chi2_table, filename, 'Sheet', 'Statistics', 'Range', ['A' num2str(Row_index + 2)], 'WriteRowNames', true)  

            chi2_start_row = num2str(Row_index + 2);

            end
                  
        end
          
    end


%% Parameter bounding

    if Reg

    a = table({Results.Markup_tag});

    writetable(a, filename, 'Sheet', 'Parameter bounding', 'Range', 'A1', 'WriteRowNames', true, 'WriteVariablenames', false)  

    Parameter_bounds_original = table(Results.Parameter_bounding.Original_bounds.Lower', Results.Parameter_bounding.Original_bounds.Upper', 'variablenames', {'Lower', 'Upper'}, 'rownames', Fitting_params);

    Parameter_bounds_reduced = table(Results.Parameter_bounding.Reduced_bounds.Lower', Results.Parameter_bounding.Reduced_bounds.Upper', 'variablenames', {'Lower', 'Upper'}, 'rownames', Fitting_params);

    Parameter_bounds_original.Properties.DimensionNames{1} = 'Parameter_names';

    a = table([], 'Variablename', {'Original_parameter_bounds'});

    writetable(a, filename, 'Sheet', 'Parameter bounding', 'Range', 'G1', 'WriteRowNames', true)

    writetable(Parameter_bounds_original, filename, 'Sheet', 'Parameter bounding', 'Range', 'F3', 'WriteRowNames', true)

    Parameter_bounds_reduced.Properties.DimensionNames{1} = 'Parameter_names';

    a = table([], 'Variablename', {'Reduced_parameter_bounds'});

    writetable(a, filename, 'Sheet', 'Parameter bounding', 'Range', 'B1', 'WriteRowNames', true)

    writetable(Parameter_bounds_reduced, filename, 'Sheet', 'Parameter bounding', 'Range', 'A3', 'WriteRowNames', true)


    %% Sampling

    a = table({Results.Markup_tag});

    writetable(a, filename, 'Sheet', 'Parameter sampling', 'Range', 'A1', 'WriteRowNames', true, 'WriteVariablenames', false) 

    Sample_table = table(Results.Sampling.Sample_size/prod(Results.Parameter_bounding.Original_bounds.Upper - Results.Parameter_bounding.Original_bounds.Lower), ...
    median(Results.Sampling.Cost_cut_offs), Results.Sampling.Sample_size, ...
    'variablenames', {'Sample_dens_approx', 'Avg_Cost_cut_off', 'Sample_size'});

    a = table([], 'Variablename', {'Parameter_sampling_results'});

    writetable(a, filename, 'Sheet', 'Parameter sampling', 'Range', 'C1', 'WriteRowNames', true)

    writetable(Sample_table, filename, 'Sheet', 'Parameter sampling', 'Range', 'A3', 'WriteRowNames', true)

    end


%% Full path if a partial path is used

    if exist([cd filesep filename], 'file')
        
    filename = [cd filesep filename];               
        
    end

    %% Try and format xls document
    
    try 
        
    newExcel = actxserver('excel.application');

    excelWB = newExcel.Workbooks.Open(filename, 1, false);

    newExcel.DisplayAlerts = false;

    excelWB.Sheets.Item(1).Delete; 
  
    excelWB.Save();

    excelWB.Close();

    newExcel.Quit();

    delete(newExcel);

    xlsAutoFitCol(filename, 'Global parameter estimation', 'A:M')

        if Reg

        xlsAutoFitCol(filename, 'Parameter bounding', 'A:K')

        end
    
    catch
        
    Ins = who;
    
        if ismember({'excelWB'}, Ins)
            
        % Close Workbook
        excelWB.Close();     
                        
        end
        
        if ismember({'newExcel'}, Ins)
            
        % Quit Excel
        newExcel.Quit();  
        
        end
        
    end
        
	if ismember('Statistics', Results_fields)
            
	xlsAutoFitCol(filename, 'Statistics', 'A:H')
        
	end
        
    if Reg

    xlsAutoFitCol(filename, 'Parameter sampling', 'A:D')
    
    % Colour tables
    table_start = 'I3';

    table_size = [length(Fitting_params) + 1, 5];

    Excel_colour_table(filename, 1, table_start, table_size)

    table_start = 'A3';

    table_size = [length(Fitting_params) + 1, 6];

    Excel_colour_table(filename, 1, table_start, table_size)

    else

    table_start = 'A3';

    table_size = [length(Fitting_params) + 1, 5];

    Excel_colour_table(filename, 1, table_start, table_size)

    end

	if ismember('Statistics', Results_fields)             

        if Results.Settings_used.GEARS_options.Calculate_NRMSE

        table_start = ['A' NRMSE_start_row];

            if Reg

            table_size = [length(Fitting_Exp_names) + 1, 3];

            else

            table_size = [length(Fitting_Exp_names) + 1, 2];

            end

        Excel_colour_table(filename, 2, table_start, table_size)    

            if Val           

                if Reg

                table_start = ['F' NRMSE_start_row];

                table_size = [length(Validation_Exp_names) + 1, 3];

                else

                table_start = ['E' NRMSE_start_row];

                table_size = [length(Validation_Exp_names) + 1, 2];

                end

            Excel_colour_table(filename, 2, table_start, table_size)     

            end

        end

        if Results.Settings_used.GEARS_options.Calculate_R2

        table_start = ['A' R2_start_row];

            if Reg

            table_size = [length(Fitting_Exp_names) + 1, 3];

            else

            table_size = [length(Fitting_Exp_names) + 1, 2];

            end

        Excel_colour_table(filename, 2, table_start, table_size)    

        	if Val

                if Reg

                table_start = ['F' R2_start_row];

                table_size = [length(Validation_Exp_names) + 1, 3];

                else

                table_start = ['E' R2_start_row];

                table_size = [length(Validation_Exp_names) + 1, 2];

                end

            Excel_colour_table(filename, 2, table_start, table_size)     

            end                

        end
        
        if Results.Settings_used.GEARS_options.Calculate_chi2

        table_start = ['A' chi2_start_row];

            if Reg

            table_size = [3, 3];

            else

            table_size = [3, 2];

            end

        Excel_colour_table(filename, 2, table_start, table_size)    

            if Val

            table_start = ['F' chi2_start_row];

                if Reg

                table_start = ['F' chi2_start_row];

                table_size = [3, 3];

                else
                    
                table_start = ['E' chi2_start_row];    

                table_size = [3, 2];

                end

            Excel_colour_table(filename, 2, table_start, table_size)     

            end                

        end

        if Reg

        table_start = 'A3';

        table_size = [length(Fitting_params) + 1, 3];

        Excel_colour_table(filename, 3, table_start, table_size)

        table_start = 'F3';

        table_size = [length(Fitting_params) + 1, 3];

        Excel_colour_table(filename, 3, table_start, table_size)

        table_start = 'A3';

        table_size = [2, 3];

        Excel_colour_table(filename, 4, table_start, table_size)

        end

    else

        if Reg

        table_start = 'A3';

        table_size = [length(Fitting_params) + 1, 3];

        Excel_colour_table(filename, 2, table_start, table_size)

        table_start = 'F3';

        table_size = [length(Fitting_params) + 1, 3];

        Excel_colour_table(filename, 2, table_start, table_size)

        table_start = 'A3';

        table_size = [2, 3];

        Excel_colour_table(filename, 3, table_start, table_size)    

        end

    end               

warning('ON', 'MATLAB:xlswrite:AddSheet') % Turn the warning back on 
end

function Excel_colour_table(File_name, Sheet_num, Table_start, Table_size)
% A function that formats xls tables in excel with colours.
% Filename    - The file name of the xls file. (string)
% Sheet_num   - The number of the sheet the table is in. (string)
% Table_start - The position of the top left corner of the table. (string)
% Table_size  - The dimensions of the table. (vector)

    try        
    %% find positions to colour
    
    % Connect to Excel
    newExcel = actxserver('Excel.application');

    excelWB = newExcel.Workbooks.Open(File_name, 0, false);

    ypos_start = char(strrep(regexp(Table_start,'\d*','Match'), ' ', ''));

    xpos_start = Table_start(1:end - length(ypos_start));

    Table_start_nums = [xls_column_index_finder_GEARS(xpos_start) str2double(ypos_start)];

    Column_range = {xpos_start, xls_column_index_finder_GEARS(Table_start_nums(1) + Table_size(2) - 1)};

    Row_range = [Table_start_nums(2) Table_start_nums(2) + Table_size(1)];

    title_range = [char(Column_range(1)) num2str(Row_range(1)) ':'  char(Column_range(2)) num2str(Row_range(1))];

    excelWB.Worksheets.Item(Sheet_num).Range(title_range).Interior.ColorIndex = 42;

    % WB.Worksheets.Item(sheet_num).Range(title_range).Font.color = 2;

        for Row_to_colour = Row_range(1) + 1:2:Row_range(2)

        Colour_range = [char(Column_range(1)) num2str(Row_to_colour) ':'  char(Column_range(2)) num2str(Row_to_colour)];    

        excelWB.Worksheets.Item(Sheet_num).Range(Colour_range).Interior.ColorIndex = 43;

        end        

    % Save Workbook
    excelWB.Save();

    % Close Workbook
    excelWB.Close();

    % Quit Excel
    newExcel.Quit();

    catch
        
    Ins = who;   
        
        if ismember({'excelWB'}, Ins)
            
        % Close Workbook
        excelWB.Close();     
                        
        end
        
        if ismember({'newExcel'}, Ins)
            
        % Quit Excel
        newExcel.Quit();  
        
        end    
  
    end

end

function xlsAutoFitCol(File_name, Sheet_name, Range)
% A function that sets the coloumn size of an xls spreadsheet to auto for a given range of columns
% Filename   - The file name of the xls file. (string)
% Sheet_name - The name of the sheet to be considered. (string)
% Range      - The range of columns that should be considered. (string)
    
[fpath, file, ext] = fileparts(char(File_name));

    if isempty(fpath)
    
    fpath = pwd;
    
    end

    try 
        
    newExcel = actxserver('excel.Application');    
        
    set(newExcel,'Visible',0);

    excelWB = invoke(newExcel.Workbooks, 'open', [fpath filesep file ext]);

    sheet = get(newExcel.Worksheets, 'Item',Sheet_name);

    invoke(sheet,'Activate');

    ExAct = newExcel.Activesheet;

    ExActRange = get(ExAct,'Range',Range);

    ExActRange.Select;

    invoke(newExcel.Selection.Columns,'Autofit');

    invoke(excelWB, 'Save');

    invoke(newExcel, 'Quit');

    delete(newExcel);
    
    catch
        
    Ins = who;
    
        if ismember({'Excel'}, Ins)
        
        invoke(newExcel, 'Quit');

        delete(newExcel);
    
        end
        
    end

end

