%%% Jake function for log level 3 cost interface % Start
function[Cost,g,residual] = Cost_log_level_3(Param_values, varargin)

    if strcmp(class(varargin{end}), 'function_handle') && strcmp(class(varargin{end - 1}), 'function_handle') % assume this is sufficiant to find if fjac is defined.
    
    log_var = varargin{end - 2};
        
    Param_values(log_var) = exp(Param_values(log_var));
    
    [Cost,g,residual] = feval(varargin{end-1}, Param_values, varargin{1:end - 3});

    else
        
    log_var = varargin{end - 1}; 
    
    Param_values(log_var) = exp(Param_values(log_var));
        
    [Cost,g,residual] = feval(varargin{end}, Param_values, varargin{1:end - 2});    
        
    end

end
%%% Jake function for log level 3 cost interface % End
