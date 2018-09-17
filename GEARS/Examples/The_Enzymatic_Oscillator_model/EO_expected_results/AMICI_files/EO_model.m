function [model] = EO_model()
% The model file for the EO model as created by GEARS for use with AMICI.


%% Variables

syms t alpha1 beta1 gamma1

model.sym.x = [alpha1 beta1 gamma1];


%% Initial conditions

model.sym.x0 = sym('IC', (size(model.sym.x)));


%% Parameters and Inputs

syms v_Km1_r1 L1_r2 sigma1_r2 L2_r3 d_r3 sigma2_r3 ks_r4

Fitting_params = [v_Km1_r1 L1_r2 sigma1_r2 L2_r3 d_r3 sigma2_r3 ks_r4];

model.sym.p = [Fitting_params model.sym.x0];


%% Differential equations

model.sym.xdot = sym(zeros(size(model.sym.x)));

% dalpha1
model.sym.xdot(1) = v_Km1_r1-(alpha1*sigma1_r2*(alpha1+1)*(beta1+1)^2)/(L1_r2*10^6+(alpha1+1)^2*(beta1+1)^2);

% dbeta1
model.sym.xdot(2) = (50*alpha1*sigma1_r2*(alpha1+1)*(beta1+1)^2)/(L1_r2*10^6+(alpha1+1)^2*(beta1+1)^2)-(sigma2_r3*(gamma1+1)^2*(d_r3*beta1/100+1)*beta1)/(L2_r3+(gamma1+1)^2*(d_r3*beta1/100+1)^2);

% dgamma1
model.sym.xdot(3) = (sigma2_r3*(gamma1+1)^2*(d_r3*beta1/100+1)*beta1)/(50*(L2_r3+(gamma1+1)^2*(d_r3*beta1/100+1)^2))-ks_r4*gamma1;


%% Observables

model.sym.y = [];

end

