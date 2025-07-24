function stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365,indice3m,indice1y)
% Computes the stepwise stochastic discount factor in the LGM (Bachelier) model,
% using a piecewise constant volatility defined between the dates in date_sigma.

% Determine the index of the volatility to use based on the current time horizon
if i <= indice3m
    indice = 1;                            % Up to 3 months, the first volatility is used
elseif i > indice3m && i <= indice1y
    indice = 2;                            % Between 3 months and 1 year, the second volatility is used
else
    indice = floor(i/indice1y) + 2;        % Beyond 1 year
end
if indice > 11
   indice = 11;                            % Cap at 11 to avoid going outside the sigma vector
end

int = 0;                                   % Initialization of the integral

% Calculation of the integral over multiple intervals (up to the current interval)
if indice > 1
    for ids = 1:indice-1
        % Calculation of year fractions relative to today
        yf_lunga = yearfrac(today, simulations_grid(i), ACT_365);
        yf_corta = yearfrac(today, simulations_grid(i-1), ACT_365);

        % Definition of the integrand corresponding to the interval \[date\_sigma(ids), date\_sigma(ids+1)]
        integranda = @(u) sigmaLGM_Bachelier(ids)^2 .* ((1 - exp(-a * (yf_lunga - u))) / a).^2 ...
                        - sigmaLGM_Bachelier(ids)^2 .* ((1 - exp(-a * (yf_corta - u))) / a).^2;

        % Integration over that interval
        int = int + quadgk(integranda, yearfrac(today, date_sigma(ids), ACT_365), yearfrac(today, date_sigma(ids+1), ACT_365));
    end

    % Last incomplete interval, from date\_sigma(index) to $t_{i-1}$
    integranda = @(u) sigmaLGM_Bachelier(indice)^2 .* ((1 - exp(-a * (yf_lunga - u))) / a).^2 ...
                    - sigmaLGM_Bachelier(indice)^2 .* ((1 - exp(-a * (yf_corta - u))) / a).^2;
    int = int + quadgk(integranda, yearfrac(today, date_sigma(indice), ACT_365), yearfrac(today, simulations_grid(i-1), ACT_365));
else
    % If still within the first volatility interval, integrate directly from 0 to $t_{i-1}$
    yf_lunga = yearfrac(today, simulations_grid(i), ACT_365);
    yf_corta = yearfrac(today, simulations_grid(i-1), ACT_365);
    integranda = @(u) sigmaLGM_Bachelier(1)^2 .* ((1 - exp(-a * (yf_lunga - u))) / a).^2 ...
                    - sigmaLGM_Bachelier(1)^2 .* ((1 - exp(-a * (yf_corta - u))) / a).^2;
    int = quadgk(integranda, 0, yearfrac(today, simulations_grid(i-1), ACT_365));
end

% Calculation of the coefficient C that multiplies the stochastic factor
C = (1 - exp(-a * yearfrac(simulations_grid(i-1), simulations_grid(i), ACT_365))) / a;

% Deterministic part A, includes the forward discount and the variance factor
A = fwd_disc(i) * exp(-0.5 * int);

% Calculation of the stochastic discount factor
stoch_disc = A .* exp(-xtsLGM(:,i) .* C);
end
