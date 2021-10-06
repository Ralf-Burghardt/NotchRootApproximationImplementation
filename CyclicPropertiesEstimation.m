function [E,K_,n_] = CyclicPropertiesEstimation(R_m)
% This function implements the FKM-Method proposed by Wächter et. al. to 
% estimate cyclic material properties for the steel material group. 
% [E, K_, n_] = CyclicPropertiesEstimation(R_m) 
% 
% Parameters: 
%    R_m - tensile strength in MPa
% Returns:
%    E   - Youngs Modulus 
%    K_  - cyclic hardening coefficient
%    n_  - cylic hardening exponent
%
% Ref: https://doi.org/10.21268/20161013-153328
% 
% Created by Ralf Burghardt 


a_sig       = 3.1148;
a_eps       = 1033;
b_sig       = 0.897;
b_eps       = -1.235;
eps_grenz   = 0.338;

E           = 206000;
n_          = 0.187;
K_          = (a_sig.*R_m.^b_sig)./(min(eps_grenz , a_eps.*R_m.^b_eps).^n_);

end

