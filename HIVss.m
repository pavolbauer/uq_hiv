function y=HIVss(q,data)
% evaluate the sum of squares

tdata = data.tdata;
ydata = data.ydata;

% initial values for the ODE
y0=[0.9e6, 4000, 0.1, 0.1, 1, 12]; 

disp(q)

% evaluate model with proposed Q  vector
[~, ymodel] = ode15s(@HIVfun,tdata,y0,[],q);

% sum of root squares
try
    y = sum((ydata-ymodel).^2);
catch
    disp('error!')
    y=NaN;
end