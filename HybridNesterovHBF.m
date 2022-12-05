%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HybridNesterovHBF.m
%--------------------------------------------------------------------------
% Project: Uniting Nesterov's accelerated gradient descent globally with
% heavy ball locally. Strongly convex version.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 09/29/2015 1:34:00

clear all

set(0,'defaultTextInterpreter','tex');

% global variables
global delta d beta M gamma lambda a c_0 c_10 cTilde_0 cTilde_10 d_0 d_10 tauMin tauMax c alpha V0

%%%%%%%%%%%%%%%%%%%%%
% setting the globals
%%%%%%%%%%%%%%%%%%%%%
setMinima();

% Nesterov constants
kappa = 1; % kappa >= 1, kappa = M/mu = 2/2
M = 2;
d = 1/(sqrt(kappa) + 1);
beta = (sqrt(kappa) - 1)/(sqrt(kappa) + 1);
a = d + beta/(2*kappa);

% Heavy Ball constants
gamma = 2/3; 
lambda = 40; 

% Uniting parameters for \mathcal{U}_0 and \mathcal{T}_{1,0}:
c_0 = 3000;  
c_10 = 402.83351565; 

% These are the same, since L = x^2
alpha = 1;

% eps_0 has to be bigger than eps_10
eps_0 = 20;
eps_10 = 15;

cTilde_0 = eps_0*alpha
cTilde_10 = eps_10*alpha
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)
d_10 = c_10 - (a^2)*(cTilde_10/alpha)^2 - (1/M)*((cTilde_10^2)/alpha)

% Poveda-Na constants
tauMin = 3; 
tauMax = 5.1; %5.63
c = 1/8; 

% Convergence within one percent from min
delta = 0.5;

%%%%%%%%%%%%%%%%%%%%%
% setting the locals
%%%%%%%%%%%%%%%%%%%%%

deltaVecNest = [0,0,0];
deltaVecHBF = [0,0,0];
deltaVecUniting = [0,0,0];
deltaVecPN = [0,0,0];

% initial conditions
z1_0 = 50; 
z2_0 = 0;
z2_00 = 50;
q_0 = 1;
tau_0 = tauMin; 

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];
x00 = [z1_0;z2_00;tau_0];

% simulation horizon
TSPAN_HBF=[0 400];
TSPAN=[0 10];
JSPAN = [0 5000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

[tNest,jNest,xNest] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options);

% Find the L values for Nesterov:
lNest = zeros(1,length(tNest));
[deltaVecNest lNest lDeltaNest] = timeToConv(xNest,tNest);

lNestAvg = lNest.';
% This is the dotted line indicating the average value:
PO = polyfit(tNest,log10(lNestAvg(:,1)),1);
yO = polyval(PO,tNest);

[tHBF,jHBF,xHBF] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x0,TSPAN_HBF,JSPAN,rule,options);

% Find the L values for HBF:
lHBF = zeros(1,length(tHBF));
[deltaVecHBF lHBF lDeltaHBF] = timeToConv(xHBF,tHBF);

[tPN,jPN,xPN] = HyEQsolver(@f,@g,@C,@D,...
    x00,TSPAN,JSPAN,rule,options);

% Find the L values for HAND-1:
lHAND = zeros(1,length(tPN));
[deltaVecPN lHAND lDelta] = timeToConv(xPN,tPN);

lHANDAvg = lHAND.';
% This is the dotted line indicating the average value:
PH = polyfit(tPN,log10(lHANDAvg(:,1)),1);
yH = polyval(PH,tPN);

% Finally, simulate the hybrid closed-loop heavy ball system
[tUniting,jUniting,xUniting] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options);

% Find the L values for Uniting:
lUniting = zeros(1,length(tUniting));
[deltaVecUniting lUniting lDeltaUniting] = timeToConv(xUniting,tUniting);

minarc = min([length(xUniting),length(xPN),length(xHBF),length(xNest)]);
ta = [tUniting(1:minarc),tPN(1:minarc),tHBF(1:minarc),tNest(1:minarc)];
ja = [jUniting(1:minarc),jPN(1:minarc),jHBF(1:minarc),jNest(1:minarc)];
xa = [xUniting(1:minarc,1),xPN(1:minarc,1),xHBF(1:minarc,1),xNest(1:minarc,1)];
xb = [xUniting(1:minarc,2),xPN(1:minarc,2),xHBF(1:minarc,2),xNest(1:minarc,2)];

figure(1)
clf
semilogy(tUniting,lUniting,'LineWidth',1.5);
hold on
semilogy(tHBF,lHBF,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
semilogy(tNest,lNest,'Color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'LineStyle','--');
semilogy(tNest,10.^(yO(:,1)),'k--','LineWidth',1.5);
semilogy(tPN,lHAND,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
semilogy(tPN,10.^(yH(:,1)),'k:','LineWidth',1.5);

% Plotting times to convergence:
semilogy(deltaVecPN(3),lDelta,'k.','MarkerSize', 20)
strDelta = [num2str(deltaVecPN(3)), 's'];
text(deltaVecPN(3),lDelta,strDelta,'HorizontalAlignment','left','VerticalAlignment','bottom');
semilogy(deltaVecUniting(3),lDeltaUniting,'k.','MarkerSize', 20)
strDeltaUniting = [num2str(deltaVecUniting(3)), 's'];
text(deltaVecUniting(3),lDeltaUniting,strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom');
semilogy(deltaVecNest(3),lDeltaNest,'k.','MarkerSize', 20)
strDeltaNest = [num2str(deltaVecNest(3)), 's'];
text(deltaVecNest(3),lDeltaNest,strDeltaNest,'HorizontalAlignment','left','VerticalAlignment','bottom');

hold off
axis([0 10 10^(-32) 10^(6)]);
ylabel('L(z_1)-L^*','FontSize',20)
xlabel('t','FontSize',20)
legend({'Hybrid','Heavy ball','Nesterov','Nesterov, average','HAND-2','HAND-2, average'},'Location','southwest')

axes('Position',[0.7 0.6 0.15 0.1])
box on
hold on
semilogy(tHBF,lHBF,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
semilogy(deltaVecHBF(3),lDeltaHBF,'k.','MarkerSize', 20)
strDeltaHBF = [num2str(deltaVecHBF(3)), 's'];
text(deltaVecHBF(3),lDeltaHBF,strDeltaHBF,'HorizontalAlignment','right','VerticalAlignment','bottom');
hold off
set(gca,'xtick',[0 100 200])
set(gca,'ytick',[10^(-2) 10^(2) 10^(6)])
axis([0 200 10^(-2) 10^(6)])
hold off

saveas(gcf,'Plots\Semilog','epsc')
saveas(gcf,'Plots\Semilog','png')