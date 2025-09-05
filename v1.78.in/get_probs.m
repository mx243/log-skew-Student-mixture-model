function [log_Px, log_ddf, log_Px1bar] = get_probs(Phi, C, plotx, typ, usealpha, x1_bar)
% Usage: [log_Px, log_ddf, log_Px1bar] = get_probs(Phi, C, plotx, typ, usealpha, x1_bar)
% If x1_bar is not empty, then this function returns 
%  log(P(x1 @ log(plotx), x1_bar(:, it_Np) | Phi(it_Ns)))
% and
%  log(P(x1_bar(:, it_Np) | Phi(it_Ns)))
% at the values x1 @ log(plotx), for each sample Phi and each
% patient x1_bar. 
% plotx should be a row vector of times.
% The value returned is the density at the vector x (or x1bar), of which the first 
%  component is the log of the lifetime, not the lifetime itself.
% typ is the type of proGamma used for m.
%
% If x1_bar is empty, then it returns log(P(x1 @ log(plotx) | Phi(it_Ns)). 
%
% log_Px is an Nsmpl * (Npt + 1) * Npatient matrix.  
% log_Px1bar is an Nsmpl * 1 * Npatient matrix or empty.
%
% If usealpha ~= 0, return log_ddf an Nsmpl * (Npt + 1) * Npatient matrix;
%  if usealpha > 0 then usealpha samples are drawn from the relevant Gamma in calculating
%  the ddf of the skew Student;
%  if usealpha < 0 then abs(usealpha) regualrly log-spaced values of alpha are used.
% where log_ddf(i, j, k) = 
% log(P(x1 @ log(plotx(j)), x1bar @ x1bar of patient k | the ith sample of Phi)), if x1_bar is nonempty,
% or
% log(P(x1 @ log(plotx(j)) | the ith sample of Phi)), if x1_bar is empty (in this case k can only be 1),
% i.e. this is just the ddf of log_Px.
% If usealpha = 0, return log_ddf = [].

% Change Log:
%
%     1.1          09:oct:24    mx243      First version.
%     1.51         10:oct:24    mx243      Change return values to
%                                          log_Px and log_Px1bar
%     1.52         14:oct:24    rfs34      Now requires typ as a parameter rather than assuming it;
%                                          avoids squeeze(); comments corrected and extended;
%                                          repeated calculation of log(plotx) avoided.
%     1.53         16:oct:24    rfs34      Fixed bug in setting T in line 51.
%     1.56         17:oct:24    mx243      Added usealpha.
%     1.61         26:oct:24    rfs34      Generalised usealpha to be what used to be nalpha.
%     1.71         25:jan:25    rfs34      Comments only changed.
%     1.75         19:mar:25    mx243      Now allows for general plotx, comments changed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.get_probs.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 5 || isempty(usealpha),
   usealpha = 0;
end
if nargin < 6 || isempty(x1_bar),
    x1_bar = [];
end

D = size(Phi(1).mu(:, 1), 1);
Nsmpl = size(Phi, 2);
Npatient = max(size(x1_bar, 2), 1); % Set Npatient = 1 when x1_bar = [].

% plotx = 0 : T_max / Npt : T_max;
Npt = size(plotx, 2) - 1;
logplotx = log(plotx);

log_Px = -Inf(Nsmpl, Npt + 1, Npatient);
log_ddf = -Inf(Nsmpl, Npt + 1, Npatient);

if isempty(x1_bar) % Marginal distribution of x1.   
    T = [1, zeros(1, D - 1)];
    for it_Ns = 1 : Nsmpl 
        for it_c = 1 : C(it_Ns)
            mu_1 = Phi(it_Ns).mu(1, it_c);
            nu_1 = Phi(it_Ns).nu(1, it_c);
            S = Phi(it_Ns).S(:, :, it_c);
            S_1 = (T * (S \ T')) ^ -1;
            m = Phi(it_Ns).m(it_c);
            r = r_from_m(m, typ);

            tmp_log_Px = log_density_skewStu(logplotx, mu_1, S_1, m, r, nu_1, 1); % Marginal distribution of x1.
            tmp_log_Px = tmp_log_Px + log(Phi(it_Ns).p(it_c)); % Weighted average wrt p.

            log_Px(it_Ns, :, 1) = -Eadd(-log_Px(it_Ns, :, 1), -tmp_log_Px); % P(x1 @ log(plotx) | Phi(it_Ns))

            if usealpha
                tmp_log_ddf = log_ddf_skewStu(logplotx, mu_1, S_1, m, r, nu_1, 1, usealpha); % Ddf of the marginal distribution of x1.
                tmp_log_ddf = tmp_log_ddf + log(Phi(it_Ns).p(it_c)); % Weighted average wrt p.

                log_ddf(it_Ns, :, 1) = -Eadd(-log_ddf(it_Ns, :, 1), -tmp_log_ddf); % P(x1 > log(plotx) | Phi(it_Ns))
            end
        end

        % Factor from the chain rule going from the pdf of log(x) to the pdf of x
        % is *not* added here because it is added where this function is called in get_prediction.m .
        % log_Px(it_Ns, :, 1) = log_Px(it_Ns, :, 1) - logplotx; 
        % log_Px(it_Ns, 1, 1) = -Inf; % Otherwise would be NaN due to Inf - Inf.

        fprintf('\rget_Probs: ...done %d of %d...', it_Ns, Nsmpl);
    end
    fprintf('\n');

    log_Px1bar = zeros(Nsmpl, 1, Npatient); % No x1_bar.

else % Joint distribution of x1 and x1bar along with marginal distribution of x1bar.
    log_Px1bar = -Inf(Nsmpl, 1, Npatient); % Dimensions chosen for easier manipulation.

    tmpx = NaN(D, Npt + 1, Npatient); % Feed into log_density_skewStu.
    for it_Np = 1 : Npatient
        tmp_x1_bar = x1_bar(:, it_Np);
        tmpx(1, :, it_Np) = logplotx;
        tmpx(2 : D, :, it_Np) = tmp_x1_bar(:, ones(1, Npt + 1));
    end

    T = eye(D);
    T = T(2 : D, :);    
    for it_Ns = 1 : Nsmpl 
        for it_c = 1 : C(it_Ns)
            mu = Phi(it_Ns).mu(:, it_c);
            S = Phi(it_Ns).S(:, :, it_c);
            m = Phi(it_Ns).m(it_c);
            r = r_from_m(m, typ);
            nu = Phi(it_Ns).nu(:, it_c);
            mu_bar = mu(2 : D);
            S_bar = (T * (S \ T')) ^ -1;
            nu_bar = nu(2 : D);

            % Pdf of x1_bar evaluated at the given x1_bar.
            % Note r = m in the model for primary data analysis, but best not to assume that.
            tmp_log_Px1bar = log_density_skewStu(x1_bar, mu_bar, S_bar, m, r, nu_bar, D - 1);
            tmp_log_Px1bar = tmp_log_Px1bar + log(Phi(it_Ns).p(it_c)); % Weighted average wrt p.

            % P(x1_bar(:, it_Np) | Phi(it_Ns)) for each it_Np:
            log_Px1bar(it_Ns, 1, :) = -Eadd(-log_Px1bar(it_Ns, 1, :), -reshape(tmp_log_Px1bar, [1, 1, Npatient])); 

            for it_Np = 1 : Npatient
                tmp_log_Px = log_density_skewStu(tmpx(:, :, it_Np), mu, S, m, r, nu, D);
                tmp_log_Px = tmp_log_Px + log(Phi(it_Ns).p(it_c)); % Weighted average wrt p.

                log_Px(it_Ns, :, it_Np) = -Eadd(-log_Px(it_Ns, :, it_Np), -tmp_log_Px); % P(x1 @ log(plotx), x1_bar(:, it_Np) | Phi(it_Ns))

                if usealpha
                    tmp_log_ddf = log_ddf_skewStu(tmpx(:, :, it_Np), mu, S, m, r, nu, D, usealpha); % Ddf of the conditional distribution of x1.
                    tmp_log_ddf = tmp_log_ddf + log(Phi(it_Ns).p(it_c)); % Weighted average wrt p.

                    log_ddf(it_Ns, :, it_Np) = -Eadd(-log_ddf(it_Ns, :, it_Np), -tmp_log_ddf); % P(x1 > log(plotx), x1_bar(:, it_Np) | Phi(it_Ns))
                end
            end
        end
        fprintf('\rget_Probs: ...done %d of %d...', it_Ns, Nsmpl);
    end   
    fprintf('\n');
end

if ~usealpha
    log_ddf = [];
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
