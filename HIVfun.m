function dy=HIVfun(t,y,q)
% -- model function for the HIV model --
% P. Bauer 2017/01/24
dy=y;

% INPUT
T1 = y(1);  % uninfected T-lymphocytes
T2 = y(2);  % uninfected macrophages
T1i= y(3);  % infected T-lymphocytes
T2i= y(4);  % ubfected macrophages
V  = y(5);  % free virus
E  = y(6);  % immune effector cells

% PARAMETERS (UQ book page 55)
lambda1 = 1e4;
lambda2 = 31.98;
c=13;
deltaE=0.1;
d2=0.01;
m1=1e-5;
p1=1;
Kb=100;
epsilon=0;
f=0.34;
m2=1e-5;
p2=1;
dE=0.25;
k1=8e-7;
Nt=100;
lambdaE=1;
Kd=500;

% ESTIMATED PARAMETERS 
d1= q(1);
k2= q(2);
delta = q(3);
bE = q(4);

% OUTPUT
dy(1) = lambda1 - d1*T1 - (1-epsilon)*k1*V*T1;
dy(2) = lambda2 - d2*T2 - (1-epsilon)*k2*V*T2;
dy(3) = (1-epsilon)*k1*V*T1 - delta*T1i - m1*E*T1i;
dy(4) = (1-f*epsilon)*k2*V*T2 - delta*T2i - m2*E*T2i;
dy(5) = Nt*delta*(T1i+T2i) - c*V - ((1-epsilon)*p1*k1*T1 + ...
 (1-f*epsilon)*p2*k2*T2)*V;
dy(6) = lambdaE + (bE*(T1i+T2i))/(T1i+T2i+Kb)*E - ...
(dE*(T1i+T2i))/(T1i+T2i+Kd)*E - deltaE*E;
