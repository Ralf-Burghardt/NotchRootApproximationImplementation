function [ sigma, i ] = Neuber_quotient( E, K_, n_, K_p, sigma_e, mode, eps)
%Neuber_quotient performs the notch root approximation according to the
% extended version of Neuber's rule. In this case the root finding problem
% is formulated using the quotient between the elastic-plastic strain
% calculated using the Ramberg-Osgood equation (epsilon_el_pl_RO) and the
% elastic-plastic strain using Neuber's rule (epsilon_el_pl_NRA) 
% The root finding problem therefore is: 
% 0 = epsilon_el_pl_RO / epsilon_el_pl_NRA - 1
%
%[ sigma ] = Neuber_quotient( E, K_, n_, K_p, sigma_e, mode, eps)
%
%Parameter:
%   E           - Youngs Modulus 
%   K_          - cyclic hardening coefficient
%   n_          - cylic hardening exponent
%   K_p         - limit load factor
%   sigma_e     - linear-elastic determined stress in MPa
%   mode        - termination criterion for Newton's method
%                 1 --> 10 steps (Standard in FKM-Guideline nonlinear) 
%                 2 --> abs(f(x)) < tol 
%   eps         - startvalue of the tolerance to be reached
%
%Returns:
%   sigma       - notch stress in MPa
%   i           - number of iterations performed in Newton's method
%
%Created by: Ralf Burghardt, Dr.-Ing. Michael Waechter

f       = @(sig_x) (sig_x ./ (E) + 2*( sig_x / (2*K_)) .^ (1/n_)) / (sigma_e/sig_x * K_p * (((sigma_e/K_p) / (E)) + 2*((sigma_e/K_p) / (2*K_))^(1/n_))) - 1;
df      = @(sig_x) (sig_x/E + 2*(sig_x/(2*K_))^(1/n_))/(K_p*sigma_e*(2*(sigma_e/(2*K_*K_p))^(1/n_) + sigma_e/(E*K_p))) + (sig_x*(1/E + (sig_x/(2*K_))^(1/n_ - 1)/(K_*n_)))/(K_p*sigma_e*(2*(sigma_e/(2*K_*K_p))^(1/n_) + sigma_e/(E*K_p)));

[sigma, i]  = newton(f,df,sigma_e,mode, eps);