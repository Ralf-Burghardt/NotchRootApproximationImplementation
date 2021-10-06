% This script provides examples of how the notch approximation methods
% available in this folder can be invoked and used. 
% 
% This script was used (in a much extended version) for the investigations 
% in the paper 
% "Determination of local stresses and strains within
% the notch strain approach: Efficient implementation of notch root approximations".
% 
% Created by: Ralf Burghardt, Dr.-Ing. Michael Waechter

clear
clc
close all

% Get the cyclic material properties for the steel material group by using 
% the method of Waechter using a tensile strength of 800 MPa 
R_m = 800;
[E,K_,n_] = CyclicPropertiesEstimation(R_m);

% Additionally to the cyclic material properties a limit load factor is
% necessary. In this example we consider a notched bar (K_t = 3)
% under tensile stress where the limit load factor equals the stress 
% concentration factor (K_p = K_t)
K_p = 3;


%% Notch root approximation with extended version of Neuber's rule (Variant 1 in the paper)
% using the difference as formulation for the root finding problem and 
% 10 Iteration steps in Newton's method 
% (proposed method in the FKM Guideline Nonlinear)

% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = Neuber_difference(E,...
                                          K_,...
                                          n_,...
                                          K_p,...
                                          AST(k,2),...
                                          1,...         % 10 Step Newton
                                          0);           % not required
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with extended version of Neuber's rule (Variant 2 in the paper)
% using the difference as formulation for the root finding problem and
% a tolerance eps = 10^-7 in Newton's method
eps = 10^-7;
% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = Neuber_difference(E,...
                                          K_,...
                                          n_,...
                                          K_p,...
                                          AST(k,2),...
                                          2,...        % termination if abs(f) < eps 
                                          eps);       
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with extended version of Neuber's rule (Variant 3 in the paper)
% using the quotient as formulation for the root finding problem and
% a tolerance eps = 10^-3 in Newton's method
eps = 10^-3;
% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = Neuber_quotient(E,...
                                          K_,...
                                          n_,...
                                          K_p,...
                                          AST(k,2),...
                                          2,...        % termination if abs(f) < eps 
                                          eps);       
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with the Seeger and Beste approach (Variant 4 in the paper)
% using the difference as formulation for the root finding problem and 
% 10 Iteration steps in Newton's method
% (proposed method in the FKM Guideline Nonlinear)

% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the oad-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = SeegerBeste_difference(E,...
                                               K_,...
                                               n_,...
                                               K_p,...
                                               AST(k,2),...
                                               1,...         % 10 Step Newton
                                               0);           % not required
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with the Seeger and Beste approach (Variant 5 in the paper)
% using the difference as formulation for the root finding problem and
% a tolerance eps = 10^-7 in Newton's method
eps = 10^-7;
% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = SeegerBeste_difference(E,...
                                               K_,...
                                               n_,...
                                               K_p,...
                                               AST(k,2),...
                                               2,...        % termination if abs(f) < eps 
                                               eps);       
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with the Seeger and Beste approach (Variant 6 in the paper)
% using the quotient as formulation for the root finding problem and
% a tolerance eps = 10^-3 in Newton's method
eps = 10^-3;
% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = SeegerBeste_quotient(E,...
                                             K_,...
                                             n_,...
                                             K_p,...
                                             AST(k,2),...
                                             2,...        % termination if abs(f) < eps 
                                             eps);       
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with the Seeger and Beste approach (Variant 7 in the paper)
% using the difference as formulation for the root finding problem,  
% 10 Iteration steps in Newton's method and an approximation of the
% derivative

% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = SeegerBeste_difference_numdiff(E,...
                                                       K_,...
                                                       n_,...
                                                       K_p,...
                                                       AST(k,2),...
                                                       1,...         % 10 Step Newton
                                                       0);           % not required
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with the Seeger and Beste approach (Variant 8 in the paper)
% using the difference as formulation for the root finding problem,
% a tolerance eps = 10^-7 in Newton's method and an approximation of the
% derivative
eps = 10^-7;
% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = SeegerBeste_difference_numdiff(E,...
                                                       K_,...
                                                       n_,...
                                                       K_p,...
                                                       AST(k,2),...
                                                       2,...        % termination if abs(f) < eps 
                                                       eps);       
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

%% Notch root approximation with the Seeger and Beste approach (Variant 9 in the paper)
% using the quotient as formulation for the root finding problem,
% a tolerance eps = 10^-3 in Newton's method and an approximation of the
% derivative
eps = 10^-3;
% Setup the discretised load-notch-strain-curve for the hysteresis branch
% (200 steps as proposed in the FKM Guideline Nonlinear)
AST = zeros(200,4);
AST(2:200,1) = (1:1:200-1)';

% The maximum linear-elastic determined stress in MPa
sig_e_max = 1200;
% As we calculate the load-notch-strain-curve for the hysteresis branch the
% stress range is the doubled maximum linear-elastic determined stress
AST(2:200,2) = (sig_e_max/((199)/2):sig_e_max/((199)/2):(2*sig_e_max))';

for k=2:1:200
    tic
    [AST(k,4), steps] = SeegerBeste_quotient_numdiff(E,...
                                                     K_,...
                                                     n_,...
                                                     K_p,...
                                                     AST(k,2),...
                                                     2,...        % termination if abs(f) < eps 
                                                     eps);       
    time = toc;
    
    AST(k,3) = ( AST(k,4) / E + 2 * (AST(k,4)/(2*K_))^(1/n_));
end

