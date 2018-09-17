function [model] = FHN_model()
% The model file for the FHN model as created by GEARS for use with AMICI.


%% Variables

syms t V R

model.sym.x = [V R];


%% Initial conditions

model.sym.x0 = sym('IC', (size(model.sym.x)));


%% Parameters and Inputs

syms a b g

Fitting_params = [a b g];

model.sym.p = [Fitting_params model.sym.x0];


%% Differential equations

model.sym.xdot = sym(zeros(size(model.sym.x)));

% dV
model.sym.xdot(1) = g*(V-V^3/3+R);

% dR
model.sym.xdot(2) = (-1/g)*(V-a+b*R);


%% Observables

model.sym.y = [];

end

