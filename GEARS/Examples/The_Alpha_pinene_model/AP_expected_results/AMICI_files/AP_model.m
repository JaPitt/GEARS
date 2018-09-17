function [model] = AP_model()
% The model file for the AP model as created by GEARS for use with AMICI.


%% Variables

syms t x1 x2 x3 x4 x5

model.sym.x = [x1 x2 x3 x4 x5];


%% Initial conditions

model.sym.x0 = sym('IC', (size(model.sym.x)));


%% Parameters and Inputs

syms p1 p2 p3 p4 p5

Fitting_params = [p1 p2 p3 p4 p5];

model.sym.p = [Fitting_params model.sym.x0];


%% Differential equations

model.sym.xdot = sym(zeros(size(model.sym.x)));

% dx1
model.sym.xdot(1) = -(p1+p2)*x1;

% dx2
model.sym.xdot(2) = p1*x1;

% dx3
model.sym.xdot(3) = p2*x1-(p3+p4)*x3+p5*x5;

% dx4
model.sym.xdot(4) = p3*x3;

% dx5
model.sym.xdot(5) = p4*x3-p5*x5;


%% Observables

model.sym.y = [];

end

