function stoch_disc=computeStochDiscHW(today,simulations_grid, i, ACT_365,a,sigmaBach,xts,fwd_disc)
% Computes the stochastic discount factor in the Hull-White model (Bachelier version),
% between times simulations_grid(i-1) and simulations_grid(i), based on the value of the process xts

% Calculation of the year fractions from today to the grid dates i and i-1
yf_lunga=yearfrac(today,simulations_grid(i),ACT_365);          % Year fraction up to t_i
yf_corta=yearfrac(today,simulations_grid(i-1),ACT_365);        % year fraction uo to t_{i-1}

% Definition of the integrand function for computing the deterministic part of the volatility
% (difference between the two time-integrated terms used to obtain the variance)
integranda = @(u) sigmaBach^2.*((1-exp(-a*(yf_lunga-u)))./a).^2 ...
                   - sigmaBach^2.*((1-exp(-a*(yf_corta-u)))./a).^2;

% Numerical integration of the integrand between 0 and t_{i-1} to obtain the incremental variance
int=quadgk(integranda, 0, yearfrac(today ,simulations_grid(i-1),ACT_365));

% Compute the coefficient C that multiplies the process xts in the discount factor
C=(1-exp(-a*yearfrac(simulations_grid(i-1),simulations_grid(i),ACT_365)))/a;

% Deterministic part A of the discount factor, including the forward discount and the variance term
A=fwd_disc(i)*exp(-0.5*int);

% Final computation of the stochastic discount factor: deterministic part times stochastic part
stoch_disc = A.*exp(-xts(:,i).*C);
