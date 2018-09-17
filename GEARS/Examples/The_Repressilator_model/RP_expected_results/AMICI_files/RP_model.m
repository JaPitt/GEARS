function [model] = RP_model()
% The model file for the RP model as created by GEARS for use with AMICI.


%% Variables

syms t p1 p2 p3 m1 m2 m3

model.sym.x = [p1 p2 p3 m1 m2 m3];


%% Initial conditions

model.sym.x0 = sym('IC', (size(model.sym.x)));


%% Parameters and Inputs

syms alpha0 alpha1 n beta1

Fitting_params = [alpha0 alpha1 n beta1];

model.sym.p = [Fitting_params model.sym.x0];


%% Differential equations

model.sym.xdot = sym(zeros(size(model.sym.x)));

% dp1
model.sym.xdot(1) = beta1*(m1-p1);

% dp2
model.sym.xdot(2) = beta1*(m2-p2);

% dp3
model.sym.xdot(3) = beta1*(m3-p3);

% dm1
model.sym.xdot(4) = alpha0+alpha1/(1+p3^n)-m1;

% dm2
model.sym.xdot(5) = alpha0+alpha1/(1+p1^n)-m2;

% dm3
model.sym.xdot(6) = alpha0+alpha1/(1+p2^n)-m3;


%% Observables

model.sym.y = [];

end

