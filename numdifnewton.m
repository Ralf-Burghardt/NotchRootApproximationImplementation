function [ x, i ] = numdifnewton( f, x, modus, eps )
% numdifnewton applies Newton's method to the function 
% using a numerical approximation of the derivative and a startvalue. 
% NOTE: This version of Newton's method is adapted to perform
% notch root approximations and may only be suitable for this application.
%
% [ x, i ] = numdifnewton( f, x, mode, eps )
%
%Parameter:
%   f       - functionhandler of the function 
%   df      - functionhandler of the derivative 
%   x       - startvalue
%   mode    - termination criterion 
%             1 --> 10 steps (Standard in FKM-Guideline nonlinear) 
%             2 --> abs(f(x)) < tol 
%   eps     - startvalue of the tolerance to be reached
%
%Returns:
%   x       - resulting value which meets the termination criterion
%   i       - number of performed steps
%
%Created by: Ralf Burghardt, Dr.-Ing. Michael Waechter

tol = eps;
switch modus
    case 1
        for i=1:1:10      
            x_values = [x - 0.0001, x, x + 0.0001];
            y_values = f(x_values);
            
            x = x - (y_values(2)/((y_values(3) - y_values(1)) / (0.0002)));
        end
    case 2
        startvalue = x;
        x_values = [x - 0.0001, x, x + 0.0001];
        y_values = f(x_values);
        i = 1; 
        while abs(y_values(2)) > tol %
            x_values = [x - 0.0001, x, x + 0.0001];
            y_values = f(x_values);
            x = x - (y_values(2)/((y_values(3) - y_values(1)) / (0.0002)));
            i = i + 1;
            if i > 5000
                % When using Newton's method for performing 
                % notch root approximations, the resulting stresses are always
                % smaller than the starting values, therefore, 
                % for a working convergence, the starting value must 
                % be reset when the accuracy barrier is raised.
                if x > startvalue
                    x = startvalue;
                end
                i = 0;
                tol = tol * 10;
            end
            if (tol> 1.0e-2)
                error('Tolerance limit exceeded');
            end
        end 
    otherwise 
        error("Non supported Mode for the application of Newton's method selected");
end
if(tol>eps)
    fprintf(sprintf('Newton applied with reduced tolerance %d  \n', tol));
end
end
