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


function [Simulate] = Generate_simulation_handle_GEARS(Model_name, Fixed_params, Param_initials)
% A function that creates the simulation handle used to call AMICI in the GEARS package. This function does not create the AMICI files. 
% Model_name   - The name of the model (string)
% Fixed_params - The names of the fixed parameters in the system. Only checks if there are any. (cell array) Optional
% Param_initials - A flag regarding if any initial conditions contain parameters

Ins = who;

Simulation_name = ['simulate_' Model_name '_model'];

    if ~ismember('Fixed_params', Ins)
   
    Fixed_params = [];
        
    end
    
    if ~ismember('Param_initials', Ins)
   
    Param_initials = false;
        
    end

    if isempty(Fixed_params)
        
    Simulate = @(Param_values, Data, Int_opts, Time_points) feval(Simulation_name, Time_points, [Param_values; Data.Initial_conditions'], [], [], Int_opts);  
    
    else 
        
    Simulate = @(Param_values, Data, Int_opts, Time_points) feval(Simulation_name, Time_points, [Param_values; Data.Initial_conditions'; Data.Fixed_params], [], [], Int_opts);
    
    end
    
    if Param_initials
        
        if isempty(Fixed_params)
        
        Simulate = @(Param_values, Data, Int_opts, Time_points) feval(Simulation_name, Time_points, [Param_values(~Data.Initial_params_only); eval(['[' num2str(Data.Initial_conditions) ']'''])], [], [], Int_opts);  
    
        else 
        
        Simulate = @(Param_values, Data, Int_opts, Time_points) feval(Simulation_name, Time_points, [Param_values(~Data.Initial_params_only); eval(['[' num2str(Data.Initial_conditions) ']''']); Data.Fixed_params], [], [], Int_opts);
    
        end   
             
    end
        
end

