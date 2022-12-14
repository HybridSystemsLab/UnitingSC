function xdot = f(x)
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: f.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Strongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2015 1:34:00
   
% The global variables
global c

% state
z1 = x(1);
z2 = x(2);
tau = x(3);

% z1dot = z2;
% z2dot = -(3/tau)*z2 - 4*c*GradientL(z1);

z1dot = (2/tau)*(z2 - z1);
z2dot = -2*c*tau*GradientL(z1);

xdot = [z1dot;z2dot;1]; 
end