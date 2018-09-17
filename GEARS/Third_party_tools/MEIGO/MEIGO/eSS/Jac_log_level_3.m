%%% Jake Interface for log level 3 jacobian % Start
function[Jobj, Jg, Jres] = Jac_log_level_3(Param_values, varargin)

log_var = varargin{end - 2};
      
Param_values(log_var) = exp(Param_values(log_var));

[Jobj, Jg, Jres] = feval(varargin{end}, Param_values, varargin{1:end - 3});

end
%%% Jake Interface for log level 3 jacobian % End
