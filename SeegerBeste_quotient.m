function [ sigma, i ] = SeegerBeste_quotient( E, K_, n_, K_p, sigma_e, mode,eps)
%SeegerBeste_quotient performs the notch root approximation according to the
% the Seeger and Beste approach. In this case the root finding problem
% is formulated using the quotient between the elastic-plastic strain
% calculated using the Ramberg-Osgood equation (epsilon_el_pl_RO) and the
% elastic-plastic strain using the Seeger and Beste approach (epsilon_el_pl_NRA) 
% The root finding problem therefore is: 
% 0 = epsilon_el_pl_RO / epsilon_el_pl_NRA - 1
%
%[ sigma ] = SeegerBeste_quotient( E, K_, n_, K_p, sigma_e, mode)
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

u           = @(sigma) ((pi / 2) * (( ((sigma_e/(sigma) - 1)) / (K_p -1))));
e_Stern     =  (sigma_e/K_p)/E + 2 * ((sigma_e/K_p)/(2*K_))^(1/n_);

f           = @(sigma) ( sigma / E + 2 * (sigma/(2 * K_))^(1/n_)) / ( sigma/E * ( (sigma_e / sigma)^2 * 2/(u(sigma))^2 * log(1/cos(u(sigma)))- sigma_e/sigma + 1) * ((e_Stern * E * K_p)/sigma_e) ) - 1;
df          = @(sigma) (sigma_e*(1/E + (sigma/(2*K_))^(1/n_ - 1)/(K_*n_)))/(K_p*e_Stern*sigma*((8*sigma_e^2*log(1/cos((pi*(sigma_e/sigma - 1))/(2*(K_p - 1))))*(K_p - 1)^2)/(pi^2*sigma^2*(sigma_e/sigma - 1)^2) - sigma_e/sigma + 1)) - (sigma_e*(sigma/E + 2*(sigma/(2*K_))^(1/n_)))/(K_p*e_Stern*sigma^2*((8*sigma_e^2*log(1/cos((pi*(sigma_e/sigma - 1))/(2*(K_p - 1))))*(K_p - 1)^2)/(pi^2*sigma^2*(sigma_e/sigma - 1)^2) - sigma_e/sigma + 1)) - (sigma_e*(sigma/E + 2*(sigma/(2*K_))^(1/n_))*(sigma_e/sigma^2 - (16*sigma_e^2*log(1/cos((pi*(sigma_e/sigma - 1))/(2*(K_p - 1))))*(K_p - 1)^2)/(pi^2*sigma^3*(sigma_e/sigma - 1)^2) + (16*sigma_e^3*log(1/cos((pi*(sigma_e/sigma - 1))/(2*(K_p - 1))))*(K_p - 1)^2)/(pi^2*sigma^4*(sigma_e/sigma - 1)^3) - (4*sigma_e^3*sin((pi*(sigma_e/sigma - 1))/(2*(K_p - 1)))*(K_p - 1))/(pi*sigma^4*cos((pi*(sigma_e/sigma - 1))/(2*(K_p - 1)))*(sigma_e/sigma - 1)^2)))/(K_p*e_Stern*sigma*((8*sigma_e^2*log(1/cos((pi*(sigma_e/sigma - 1))/(2*K_p - 2)))*(K_p - 1)^2)/(pi^2*sigma^2*(sigma_e/sigma - 1)^2) - sigma_e/sigma + 1)^2);
start_value = sigma_e * ( 1 - ( (1 - (1/K_p))/1000));

[sigma, i] = newton(f,df,start_value,mode,eps);