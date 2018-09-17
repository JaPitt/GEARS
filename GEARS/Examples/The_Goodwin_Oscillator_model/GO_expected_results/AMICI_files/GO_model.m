function [model] = GO_model()
% The model file for the GO model as created by GEARS for use with AMICI.


%% Variables

syms t x1 x2 x3

model.sym.x = [x1 x2 x3];


%% Initial conditions

model.sym.x0 = sym('IC', (size(model.sym.x)));


%% Parameters and Inputs

syms k1 k2 k3 k4 k5 k6 Ki n

Fitting_params = [k1 k2 k3 k4 k5 k6 Ki n];

model.sym.p = [Fitting_params model.sym.x0];


%% Differential equations

model.sym.xdot = sym(zeros(size(model.sym.x)));

% dx1
model.sym.xdot(1) = (k1*Ki^n)/(Ki^n+x3^n)-k2*x1;

% dx2
model.sym.xdot(2) = k3*x1-k4*x2;

% dx3
model.sym.xdot(3) = k5*x2-k6*x3;


%% Observables

model.sym.y = [];

end

