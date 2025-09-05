function j_vec = calc_j(trained_smpl_file, validation_data_file, lambda, Nskip)
% Usage: j_vec = calc_j(trained_smpl_file, validation_data_file, lambda, Nskip)
% This function returns a vector j_vec of values of j(x1, x1_bar), where
% the x's are or are found in validation_data_file, using the trained model at 
% trained_smpl_file.
% lambda is the rate constant of the exponential reference prior in days ^ (-1).

% Change Log:
%
%     1.1          12:sep:24    mx243      First version.
%     1.40         20:sep:24    rfs34      Edited with actions for mx243 after review by rfs34.
%     1.41         23:sep:24    mx243      Changes made following the
%                                          instructions. This is modified
%                                          from the first part of calc_ASI.
%     1.42         29:sep:24    mx243      Now can take in
%                                          validation_data_file as the data itself.
%     1.45         03:oct:24    rfs34      Layout improved and prior on lambda changed to match tow24's.
%     1.46         03:oct:24    rfs34      Calculation of lambda moved to run_calc_j in order to be 
%                                          sure that the same value of lambda is always used by 
%                                          both models and all splits; adapted to new indices coming from
%                                          primary sample file.
%     1.50         08:oct:24    mx243      Fixed error in the calculation
%                                          of j(x, x_\bar{1}).
%     1.62         28:oct:24    mx243      Added Nskip option; rfs34 reduced nalpha to 1e3.
%     1.63         29:oct:24    mx243      Bug fixed by including the
%                                          coefficients from chain rule for uncensored patients.
%     1.71         25:jan:25    rfs34      Switched to numerical integration over alpha for the ddf.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.calc_j.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 4 || isempty(Nskip),
    Nskip = 1;
end

% Load the trained model.
load(trained_smpl_file, 'Samples', 'LN0', 'LN1');
Phi = Samples.phi(LN0 : LN1);
C = Samples.C(LN0 : LN1);

if ischar(validation_data_file)

    % This is only used for testing. x_unseen is loaded from data file made by run_calc_logj.
    load(validation_data_file, 'x_unseen'); 
    x_ob = x_unseen;

else

    % If validation_data_file is not a char vector, the it's the data itself.
    x_ob = validation_data_file; 

end

K = size(x_ob, 2);
D = size(x_ob, 1);

is_censored = x_ob(1, :) < 0;
n_uncensored = K - sum(is_censored);

smpl_num = LN1 - LN0 + 1; % Number of samples we extracted at t_SA = 1.
p_dis = 0.15; % Discarding the first 15% of samples.

d = abs(x_ob(1, :));
x_ob(1, :) = log(abs(x_ob(1, :)));

log_P2_vec = NaN(1, K); % P2(x) ~ Exp(lambda).
for k = 1 : K
    if is_censored(k)
        log_P2_vec(k) = - lambda * d(k); % log(P2(x > d)).
    else
        log_P2_vec(k) = log(lambda) - lambda * d(k); % log(P2(x = d)).
    end
end

log_Q_c = -Inf(1, K);
norm_const = -Inf(1, K);
T = eye(D);
T = T(2 : D, :);
nalpha = -1e3; % Number of samples of alpha drawn for the integration;
               % when negative the number of log-spaced points for the integration.
cnt = 0;
for num = ceil(p_dis * smpl_num) : Nskip : smpl_num 
    for it_c = 1 : C(num)
    
        tmp_log_Q_c = NaN(1, K);
        mu = Phi(num).mu(:, it_c);
        S = Phi(num).S(:, :, it_c);
        m = Phi(num).m(it_c);
        nu = Phi(num).nu(:, it_c);
        mu_bar = mu(2 : D);
        S_bar = (T * (S \ T')) ^ -1;
        nu_bar = nu(2 : D);

        % Pdf of x1_bar evaluated at x_ob(2 : D, :). Note r = m in the model for primary data analysis.
        tmp_norm_const = log_density_skewStu(x_ob(2 : D, :), mu_bar, S_bar, m, m, nu_bar, D - 1); 
        tmp_norm_const = tmp_norm_const + log(Phi(num).p(it_c)); % Weighted average wrt p.

        norm_const = -Eadd(-norm_const, -tmp_norm_const);

        tmp_log_Q_c(~is_censored) = log_density_skewStu(x_ob(:, ~is_censored), mu, S, m, m, nu, D);
        tmp_log_Q_c(is_censored) = log_ddf_skewStu(x_ob(:, is_censored), mu, S, m, m, nu, D, nalpha);
        tmp_log_Q_c = tmp_log_Q_c + log(Phi(num).p(it_c)); % Weighted average wrt p.
        
        log_Q_c = -Eadd(-log_Q_c, -tmp_log_Q_c);

    end

    cnt = cnt + 1;
    fprintf('\rCalc_Q: ...done %d of %d...', cnt, (smpl_num - ceil(p_dis * smpl_num)) / Nskip + 1);
end

log_Q_c = log_Q_c - norm_const; % No need to subtract log(cnt) from the two terms, they cancel out.

% Coefficient from chain rule. x_ob(1, :) has already undergone log(abs()).
log_Q_c(~is_censored) = log_Q_c(~is_censored) - x_ob(1, ~is_censored); 

j_vec = log_Q_c - log_P2_vec; 

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
