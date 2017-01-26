function y=HIVpredict(time,q)
% prediction function calling the ODE15s solver
% P. Bauer 2017/01/24

y0=[0.9e6, 4000, 0.1, 0.1, 1, 12]; 
[t,y] = ode15s(@HIVfun,time,y0,[],q);
