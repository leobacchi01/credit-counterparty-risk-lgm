function [collateral_Delta, collateral_Epsilon] = Weekly_LGM_Collateral(Portfolio, jamshResults, dates, discounts, simulations_num, aLGM, sigmaBach)
% Set the seed for random number generation to ensure reproducibility
rng(42);

% Take the starting date from the first date provided
today = dates(1);

% Constants for calculating year fractions
EU_30_360 = 6;
ACT_365 = 3;

% Mean reversion parameter of the LGM model
a = aLGM;

% Define an initial time grid with today and today + 3 months
time_grid = [today, today + calmonths(3)];

% Add to the time grid the business day offset dates from today + years 1 to 10
time_grid = businessdayoffset([time_grid, today + calyears(1):calyears(1):today + calyears(10)]);

% Extract the fixed rate of the swap from the portfolio
swapRate = Portfolio.FixedLegRate;

% Extract the swap maturity
maturity = Portfolio.Maturity;

% Extract the type of fixed leg (Receive or Pay)
fixedLeg = Portfolio.FixedLeg;

% Initialize a zero vector to classify the type of swap
swapType = zeros(size(fixedLeg));

% Set +1 for Receive
swapType(strcmpi(fixedLeg,'Receive')) = 1;

% Set -1 for Pay
swapType(strcmpi(fixedLeg,'Pay')) = -1;

% Assign volatility dates to the time_grid
date_sigma = time_grid;

% Calculate notional multiplied by 1 million (scale size)
notional = 10^6 * Portfolio.NotionalAmount;

% Take the Bachelier volatility calculated from jamshResults
sigmaLGM_Bachelier = jamshResults.BACHELIER;

% Generate fixed leg payment dates at weekly intervals up to 10 years
fixed_leg_payments_dates = businessdayoffset(today:calweeks(1):today + calyears(10));

% Calculate interpolated discount factors for payment dates
discounts_timegrid = interpolation(discounts, dates, fixed_leg_payments_dates);

% Calculate time deltas (year fractions) between payment dates according to EU_30_360 convention
deltas = yearfrac(fixed_leg_payments_dates(1:end-1), fixed_leg_payments_dates(2:end), EU_30_360);

% Define the simulation grid (excluding the first payment date)
simulations_grid = fixed_leg_payments_dates(2:end);

% Number of simulation steps
N = length(simulations_grid);

% Take the counterparty name from the portfolio
counterparty = Portfolio.CounterpartyName;

% Count how many counterparties are of type 'Epsilon'
nEpsilon = sum(strcmp(counterparty, 'Epsilon'));

% Count how many counterparties are of type 'Delta'
nDelta = sum(strcmp(counterparty, 'Delta'));

% Calculate fixed cash flow multiplying rate by deltas
fixed_cash_flow = swapRate .* deltas;

% Generate a matrix of standard normal random numbers for simulations
z = randn(simulations_num, N);

% Define dimensions for output matrices
r = simulations_num; c = N;    

% Initialize 3D matrices for LGM simulations separated by counterparties Epsilon and Delta
TLGM_Epsilon = zeros(r, c, nEpsilon);
TLGM_Delta = zeros(r, c, nDelta);

% Initialize the state matrix x_t for all simulations and all times
xtsLGM = zeros(simulations_num, N);

% Counters for managing Delta and Epsilon counterparties
countDelta = 1;
countEpsilon = 1;

% Loop over the dates of the simulation grid
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx); % current simulation date
    
    if idx == 1
        prev_sim_date = today; % at the first iteration the previous date is today
        x_t = zeros(1, simulations_num); % initial LGM state at zero
    else
        prev_sim_date = simulations_grid(idx-1); % previous simulation date
    end
    
    % Simulate the state x_t at the current time with volatility varying by index
    if idx < 12
        % first 11 steps with volatility sigmaLGM_Bachelier(1)
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(1), x_t, prev_sim_date, sim_date);
    elseif idx >= 12 && idx < 52
        % from step 12 to 51 with volatility sigmaLGM_Bachelier(2)
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(2), x_t, prev_sim_date, sim_date);
    elseif idx < N
        % after the first year, volatility is taken from the index position, max 11
        indice = floor(idx/52) + 2;
        if indice > 11
            indice = 11;
        end
        x_t = simulateLGM(z(:, idx), aLGM, sigmaLGM_Bachelier(indice), x_t, prev_sim_date, sim_date);
    elseif idx == N
        % at the last step, use sigmaBach passed as parameter
        x_t = simulateLGM(z(:, idx), aLGM, sigmaBach, x_t, prev_sim_date, sim_date);
    end
    
    % Save the simulated x_t value for all simulations at the date idx
    xtsLGM(:, idx) = x_t;
end

for j=1:length(maturity)
    
    % Initialize MtM (Mark-to-Market) matrices for Epsilon and Delta counterparties for each simulation and time step
    MtMs_Epsilon = zeros(simulations_num, N);
    MtMs_Delta = zeros(simulations_num, N);

    % Build the simulation time grid with weekly steps based on the current maturity
    if maturity >= 1
        % If maturity is at least 1 year, create grid with weekly steps up to 'maturity(j)' years
        simulations_grid = businessdayoffset(today:calweeks(1):today+calyears(maturity(j)));
    else
        % If maturity is less than 1 year, create grid with weekly steps up to 'maturity(j)' months
        simulations_grid = businessdayoffset(today:calweeks(1):today+calmonths(maturity(j)*12));
    end

    % Remove the first date (today) from the simulation grid
    simulations_grid = simulations_grid(2:end);
    Nstep = numel(simulations_grid);

    % First simulation date
    sim_date = simulations_grid(1);

    % Calculate coefficients A and C via the function affine_trick_LGM,
    % which depend on the simulation date, LGM parameters, discount curves, and time grid
    [A, C] = affine_trick_LGM(sim_date, simulations_grid, aLGM, sigmaLGM_Bachelier, discounts_timegrid, date_sigma, dates, discounts);

    for idx = 1:length(simulations_grid)
        
        % Current simulation date
        sim_date = simulations_grid(idx);
        % State of the LGM process (x_t) for all scenarios, transposed to row vector
        x_t = xtsLGM(:, idx)';
        
        % Case maturity >= 1 year and simulation date before maturity
        if (maturity(j) >= 1) && sim_date < businessdayoffset(today + calyears(maturity(j)))

            % If counterparty is Epsilon
            if counterparty(j) == "Epsilon"
                % Calculate discounted cash flows using the affine formula for all scenarios
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);  
                % Sum of discounted flows
                sum_term = sum(term1, 1);
                % Terminal term, representing the final value of the exponential component
                terminal = A(end) * exp(-C(end) .* x_t);
                % Adjustment for intermediate simulation steps (adjust=1 if idx < Nstep)
                adjust = (idx < Nstep);
                % Calculate MtM multiplying by notional and swap type (pay/receive)
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                % Assign MtM to Epsilon counterparties matrix
                MtMs_Epsilon(:, idx) = MtM;
                % Shift vectors A and C for next step
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];

            % If counterparty is Delta
            elseif counterparty(j) == "Delta"
                % Similar calculations for Delta counterparties
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                MtMs_Delta(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
            end

        % Case maturity less than 1 year and sim_date before maturity in months
        elseif sim_date < businessdayoffset(today + calmonths(maturity(j)*12))

            % Logic identical to the previous block, but with maturity in months
            if counterparty(j) == "Epsilon"
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                MtMs_Epsilon(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];

            elseif counterparty(j) == "Delta"
                term1 = exp(log(fixed_cash_flow(j,1:length(A))' .* A) - C * x_t);
                sum_term = sum(term1, 1);
                terminal = A(end) * exp(-C(end) .* x_t);
                adjust = (idx < Nstep);
                MtM = notional(j)*swapType(j) * (sum_term + terminal - adjust);
                MtMs_Delta(:, idx) = MtM;
                A = [0; A(1:end-1)];
                C = [0; C(1:end-1)];
            end
        end

    end
% Assign the calculated MtM results to the TLGM tensors based on the counterparty
if counterparty(j) == "Epsilon"
    TLGM_Epsilon(:,:,countEpsilon) = MtMs_Epsilon;
    countEpsilon = countEpsilon + 1;
else
    TLGM_Delta(:,:,countDelta) = MtMs_Delta;
    countDelta = countDelta + 1;
end

end

% Sum along the third dimension (number of contracts) to obtain aggregated MtM for each counterparty
MtM_Epsilon = sum(TLGM_Epsilon, 3);
MtM_Delta = sum(TLGM_Delta, 3);

% Initialize collateral matrices (collateralization)
collateral_Epsilon = zeros(simulations_num, N);
collateral_Delta = zeros(simulations_num, N);

% Copies of MtM for collateralized calculation
MtM_coll_Epsilon = MtM_Epsilon;
MtM_coll_Delta = MtM_Delta;

% Time grid of fixed payments (excluding the first date)
simulations_grid = fixed_leg_payments_dates(2:end);

% Calculate discount factors between steps in the simulation grid
disc_simgrid = interpolation(discounts, dates, simulations_grid);
fwd_disc = disc_simgrid(2:end) ./ disc_simgrid(1:end-1);

% Reference indices for quarters and years (12 weeks = 3 months, 52 weeks = 1 year)
indice3m = 12;
indice1y = 52;

% Loop to update collateralized exposures for counterparty Epsilon
for i=1:N-1
    if i == 1
        % First step: collateral is not updated, MtM remains unchanged
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i);
    else
        % Calculate stochastic discount factor for step i
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        % Add the effect of the previously collateralized discounted collateral
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i-1) ./ stoch_disc;
        % Update collateral as the negative of the collateralized MtM exposure
        collateral_Epsilon(:, i) = -MtM_coll_Epsilon(:, i);
        % Reset current exposure to zero to avoid double counting
        MtM_coll_Epsilon(:, i) = 0;
    end
end

% Same logic for counterparty Delta
for i=1:N-1   
    if i == 1
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i);
    else
        stoch_disc = computeStochDiscStep(i, fwd_disc, simulations_grid, sigmaLGM_Bachelier, date_sigma, a, xtsLGM, today, ACT_365, indice3m, indice1y);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i-1) ./ stoch_disc;
        collateral_Delta(:, i) = -MtM_coll_Delta(:, i);
        MtM_coll_Delta(:, i) = 0;
    end
end
