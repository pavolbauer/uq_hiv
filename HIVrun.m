% Parameter estimation of the HIV model using 
% Delayed Rejection Adaptive MCMC (DRAM)
%
% P. Bauer, 2017/01/25
%
% Synthetic data for the HIV model is provided at
% http://www4.ncsu.edu/%7Ersmith/UQ_TIA/CHAPTER9/hiv_data.mat
%
% Uses the Delayed Rejection MCMC method from the MCMC toolbox,
% available at http://helios.fmi.fi/~lainema/mcmc/

% load synthetic data - contains measurement for all unknowns, taken every
% 5th day, for 200 days
load hiv_data

% seperate time and data columns into a struct
data.tdata=hiv_data(:,1);
data.ydata=hiv_data(:,2:7);
data.labels={'T1','T2','T1i','T2i','V','E'};

% initial guess - initial values of the ODE model
% Q=[d1,k2,delta,bE]
q0=[0.01, 1e-4, 0.7, 0.3];  
 
% optimize Q based on the data   
%options = optimset('Display','iter','MaxIter',400);
%[qopt,rss] = fminsearch(@HIVss,q0,options,data);

% optimal values found by fminsearch:
qopt=[0.0098,0.0001,0.6989,0.2941];
rss=6.8318e+09;

% now using DRAM
model.ssfun = @HIVss;

%parameters
params = {
    {'d1',      qopt(1), 0}
    {'k2',      qopt(2), 0}
    {'delta',   qopt(3), 0}
    {'bE',      qopt(4), 0}
    };

%prior information
model.S20 = [1];
model.N0  = [4];

% generate an initial "burn-in" chain
options.nsimu = 1000;
[results, chain, s2chain]= mcmcrun(model,data,params,options);

% observe properties and stats of the initial chain
figure
mcmcplot(chain,[],results,'pairs');
figure
mcmcplot(chain,[],results,'denspanel',2);
chainstats(chain,results)

% Now re-run starting from the results of the previous run,
% this will take about 10 minutes.
options.nsimu = 5000;
[results, chain, s2chain] = mcmcrun(model,data,params,options, results);

% observe properties and stats of the burned-in chain
figure
mcmcplot(chain,[],results,'pairs');
figure
mcmcplot(chain,[],results,'denspanel',2);
chainstats(chain,results)

% predictive model function
modelfun = @(t,q) HIVpredict(t,q);

% sample 500 parameters and calculate predictive plots
nsample = 500;
out = mcmcpred(results,chain,s2chain,data.tdata,modelfun,nsample);

%prediction plot with overlayed data
mcmcpredplot(out);
hold on
for i=1:6
  subplot(3,2,i)
  hold on
  plot(data.tdata,data.ydata(:,i),'.r'); 
  ylabel(''); title(data.labels{i});
  hold off
  plcnt=plcnt+1;
end
