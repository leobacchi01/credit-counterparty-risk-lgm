function [MtM_coll_Delta, MtM_coll_Epsilon, collateral_Delta, collateral_Epsilon] = Weekly_HW_Collateral(Portfolio, dates, discounts, simulations_num, aBach, sigmaBach)
% Set the random number generator seed for reproducibility
rng(42);

% Defines today's date as the first element of the array 'dates'
today = dates(1);

% Definitions of day count conventions
EU_30_360 = 6;
ACT_365 = 3;

% Mean reversion parameter for the LGM model
a = aBach;

% data estrapolation
swapRate = Portfolio.FixedLegRate;      % swap fixed rate
maturity = Portfolio.Maturity;           % swap maturity
fixedLeg = Portfolio.FixedLeg;           % fixed rate type (Pay or Receive)

% Encoding of swap type: +1 for Receive, -1 for Pay, 0 otherwise
swapType = zeros(size(fixedLeg));             
swapType(strcmpi(fixedLeg,'Receive')) =  1;   
swapType(strcmpi(fixedLeg,'Pay'))     = -1;   

% Calculation of the total notional (scale 10^6)
notional = 10^6 * Portfolio.NotionalAmount; 

% Construction of the fixed payment schedule with weekly steps
fixed_leg_payments_dates = businessdayoffset(today:calweeks(1):today+calyears(10));

% Interpolation of discount factors for the fixed payment dates
discounts_timegrid = interpolation(discounts, dates, fixed_leg_payments_dates);

% Calculation of year fractions between payments, according to EU_30_360 convention
deltas = yearfrac(fixed_leg_payments_dates(1:end-1), fixed_leg_payments_dates(2:end), EU_30_360);

% Defines the simulation time grid (payment dates excluding the first)
simulations_grid = fixed_leg_payments_dates(2:end);
N = length(simulations_grid);   % Number of time steps

% Extraction of counterparty names and count of how many are Delta and Epsilon
counterparty = Portfolio.CounterpartyName;
nEpsilon = sum(strcmp(counterparty, 'Epsilon'));
nDelta = sum(strcmp(counterparty, 'Delta'));

% Calculation of fixed cash flows multiplying the rate by the time deltas
fixed_cash_flow = swapRate .* deltas;

% Generation of matrix of standard Gaussian variables for Monte Carlo simulation
z = randn(simulations_num, N);

% Initialization of matrices for LGM results for both counterparties
r = simulations_num; c = N;      % Dimensions of each simulation matrix
TLGM_Epsilon = zeros(r, c, nEpsilon);
TLGM_Delta = zeros(r, c, nDelta);

% Matrix to save simulated x_t values of the Hull-White (LGM) model
xtsHW = zeros(simulations_num, N);

% Counters for Delta and Epsilon counterparties (for possible data population)
countDelta = 1; 
countEpsilon = 1;


% Temporal loop over simulation dates
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx);
    if idx == 1
        % For the first simulation date, the previous time is "today"
        prev_sim_date = today;
        % Initial value of the LGM simulation is zero (initial condition)
        x_t = zeros(1, simulations_num);
    else
        % Otherwise the previous time is the previous simulation date
        prev_sim_date = simulations_grid(idx - 1);
    end
    % Simulation of the LGM value at time sim_date using the discretization step
    x_t = simulateLGM(z(:, idx), aBach, sigmaBach, x_t, prev_sim_date, sim_date);
    % Saving simulated values into the xtsHW matrix
    xtsHW(:, idx) = x_t;
end
for j = 1:length(maturity)
    % Initialize the MtM matrices for the two counterparties at each iteration
    MtMs_Epsilon = zeros(simulations_num, N);
    MtMs_Delta = zeros(simulations_num, N);
    
    % Construction of the simulation time grid based on maturity
    if maturity >= 1
        % If maturity is at least 1 year, weekly grid up to "maturity" years
        simulations_grid = businessdayoffset(today : calweeks(1) : today + calyears(maturity(j)));
    else
        % If maturity is less than 1 year, weekly grid up to "maturity" months
        simulations_grid = businessdayoffset(today : calweeks(1) : today + calmonths(maturity(j)*12));
    end

   % The first dates are excluded from the grid (starts from the second element)
simulations_grid = simulations_grid(2:end);
Nstep = numel(simulations_grid);  % Number of steps in the simulation

% Calculation of affine coefficients A and C for each simulation date
sim_date = simulations_grid(1);
[A, C] = affine_trick(sim_date, simulations_grid, aBach, sigmaBach, discounts_timegrid, dates, discounts);

% Loop to calculate MtM at each simulation date
for idx = 1:length(simulations_grid)
    sim_date = simulations_grid(idx);
    
    % Extraction of the simulated x_t value of the LGM model for date idx
    x_t = xtsHW(:, idx)';
    
    % Condition if maturity >= 1 year and sim_date before maturity
    if (maturity(j) >= 1) && sim_date < businessdayoffset(today + calyears(maturity(j)))
        
        if counterparty(j) == "Epsilon"
            % Calculation of the exponential term for fixed cash flows weighted by affine coefficients
            term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
            sum_term = sum(term1, 1);                  % Sum of terms (1×N scenarios)
            terminal = A(end) * exp(-C(end) * x_t);    % Terminal term (1×N scenarios)
            adjust = (idx < Nstep);                     % Boolean adjustment for last step
            MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
            MtMs_Epsilon(:, idx) = MtM;                 % Saving MtMs for Epsilon
            
            % Shift A and C for the next step (prepend zero)
            A = [0; A(1:end-1)];
            C = [0; C(1:end-1)];
            
        elseif counterparty(j) == "Delta"
            % Same operations for Delta counterparties
            term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
            sum_term = sum(term1, 1);
            terminal = A(end) * exp(-C(end) * x_t);
            adjust = (idx < Nstep);
            MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
            MtMs_Delta(:, idx) = MtM;
            A = [0; A(1:end-1)];
            C = [0; C(1:end-1)];
        end
    
    % Condition if maturity < 1 year and sim_date before maturity in months
    elseif sim_date < businessdayoffset(today + calmonths(maturity(j) * 12))
        
        if counterparty(j) == "Epsilon"
            % Calculation identical to the previous case
            term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
            sum_term = sum(term1, 1);
            terminal = A(end) * exp(-C(end) * x_t);
            adjust = (idx < Nstep);
            MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
            MtMs_Epsilon(:, idx) = MtM;
            A = [0; A(1:end-1)];
            C = [0; C(1:end-1)];
            
        elseif counterparty(j) == "Delta"
            % Same operations for Delta counterparties
            term1 = exp(log(fixed_cash_flow(j, 1:length(A))' .* A) - C * x_t);
            sum_term = sum(term1, 1);
            terminal = A(end) * exp(-C(end) * x_t);
            adjust = (idx < Nstep);
            MtM = notional(j) * swapType(j) * (sum_term + terminal - adjust);
            MtMs_Delta(:, idx) = MtM;
            A = [0; A(1:end-1)];
            C = [0; C(1:end-1)];
        end
    end
end

   % Save the results in the TLGM matrices depending on the counterparty
if counterparty(j) == "Epsilon"
    TLGM_Epsilon(:, :, countEpsilon) = MtMs_Epsilon;
    countEpsilon = countEpsilon + 1;
else
    TLGM_Delta(:, :, countDelta) = MtMs_Delta;
    countDelta = countDelta + 1;
end

end

% Sum along the third dimension to get the total MtM values
MtM_Epsilon = sum(TLGM_Epsilon, 3);
MtM_Delta = sum(TLGM_Delta, 3);

% Initialization of collateral matrices for both counterparties
collateral_Epsilon = zeros(simulations_num, N);
collateral_Delta = zeros(simulations_num, N);

% Copy of total MtM for collateral calculation
MtM_coll_Epsilon = MtM_Epsilon;
MtM_coll_Delta = MtM_Delta;

% Recalculation of the time grid and forward discount factors
simulations_grid = fixed_leg_payments_dates(2:end);
disc_simgrid = interpolation(discounts, dates, simulations_grid);
fwd_disc = disc_simgrid(2:end) ./ disc_simgrid(1:end-1);

% Collateral calculation and MtM adjustment for the Epsilon counterparty
for i = 1:N-1
    if i == 1
        % For the first step, keep the MtM value unchanged
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i);
    else
        % Calculation of the stochastic discount factor for step i
        stoch_disc = computeStochDiscHW(today, simulations_grid, i, ACT_365, a, sigmaBach, xtsHW, fwd_disc);
        % Update MtM including the previous collateral discounted
        MtM_coll_Epsilon(:, i) = MtM_coll_Epsilon(:, i) + collateral_Epsilon(:, i - 1) ./ stoch_disc;
        % Calculation of the current collateral as the negative of the adjusted MtM
        collateral_Epsilon(:, i) = -MtM_coll_Epsilon(:, i);
        % Reset MtM to zero after collateral adjustment
        MtM_coll_Epsilon(:, i) = 0;
    end
end

% Collateral calculation and MtM adjustment for the Delta counterparty (analogous to Epsilon)
for i = 1:N-1
    if i == 1
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i);
    else
        stoch_disc = computeStochDiscHW(today, simulations_grid, i, ACT_365, a, sigmaBach, xtsHW, fwd_disc);
        MtM_coll_Delta(:, i) = MtM_coll_Delta(:, i) + collateral_Delta(:, i - 1) ./ stoch_disc;
        collateral_Delta(:, i) = -MtM_coll_Delta(:, i);
        MtM_coll_Delta(:, i) = 0;
    end
end
