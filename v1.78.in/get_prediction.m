function [lifetime, survival, hazardrt, MClifetime, MCsurvival, MChazardrt] = ...
         get_prediction(log_Px, log_Px1bar, plotx, mode, if_centiles, usealpha)
% Usage: [lifetime, survival, hazardrt, MClifetime, MCsurvival, MChazardrt] = ...
%        get_prediction(log_Px, log_Px1bar, plotx, mode, if_centiles, usealpha)
% This function returns different typs of the predictions of the three 
% values, and optionally their means and centiles.
% log_Px is an Nsmpl * (Npt + 1) * Npatient matrix.  
% log_Px1bar is an Nsmpl * 1 * Npatient matrix.  
%
% mode = A: return values are Npatient * (Npt + 1) matrices. Each row of
% lifetime is P(x1 | x_ob, x1_bar) = P(x | x_ob) / P(x1_bar | x_ob) 
% where the numerator and the denominator are obtained by averaging over 
% all samples of Phi.
%
% mode = B: return values are Nsmpl * (Npt + 1) matrices. Each row of
% lifetime is mean_{all patients} {P(x1 | x1_bar, Phi)} for the
% corresponding sample of Phi.
%
% mode = C: lifetime, survival, hazardrt same as when mode = B. If
% if_centiles = 1 (or 2), return the centiles of predictions with all pairs of
% (Phi, patient). If additionally the centiles of the Nsmpl curves in mode
% B are also wanted, they can be obtained by calling get_centiles
% separately.
%
% When usealpha > 0, the input log_Px is actually the decumulative
% distribution of x1. We do the averaging operations above to this 
% decumulative distribution, and get the survival data which is more
% analytic than what we'd otherwise get with numerical integration.
% In this case, the lifetime and MClifetime returned are actually this more 
% analytic version of survival data. The only change required is to suppress 
% the coefficient from the chain rule when calculating pdf.

% Change Log:
%
%     1.1          09:oct:24    mx243      First version.
%     1.52         14:oct:24    rfs34      Minor speedup in mode C, and avoided squeeze() throughout.
%     1.56         17:oct:24    mx243      Added usealpha.
%     1.71         25:jan:25    rfs34      Comments only changed.
%     1.75         19:mar:25    mx243      Now allows for general plotx, fixed matrix size on line 106.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.get_prediction.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if mode == 'C' && if_centiles == 0
    mode = 'B';
end

Nsmpl = size(log_Px, 1);
Npatient = size(log_Px, 3);

% plotx = 0 : T_max / Npt : T_max;
Npt = size(plotx, 2) - 1;

if mode == 'A'
    % Averaged across samples of Phi, giving an Npatient * (Npt + 1) matrix:
    log_Px = reshape(-Esum(-log_Px, 1), [(Npt + 1), Npatient]).'; 

    % Averaged across samples of Phi, giving an Npatient * 1 matrix.
    % Note: This was actually wrong in v1.51, and would have resulted in patient 1 being
    %       subtracted from all patients in the following paragraph.
    log_Px1bar = reshape(-Esum(-log_Px1bar, 1), [1, Npatient]).'; 

    % Make the conditional distribution. No need to subtract log(Nsmpl) from the two terms, they cancel out.
    log_lifetime = log_Px - log_Px1bar(:, ones(1, Npt + 1));
    
    if ~usealpha % In this case log_lifetime is the pdf of x1, so need this coefficient. In the other case, log_lifetime is the ddf of x1, so no need for this coefficient.
        log_lifetime = log_lifetime - log(plotx(ones(Npatient, 1), :)); % Factor from the chain rule going from the pdf of log(x) to the pdf of x.
        if plotx(1) == 0
            log_lifetime(:, 1) = -Inf; % Otherwise would be NaN due to Inf - Inf.
        end
    end

    [lifetime, survival, hazardrt] = get_LSH(log_lifetime, plotx, usealpha); % If usealpha ~= 0, survival and hazardrt here are empty, while lifetime is actually the survival data.
end

if mode == 'B'
    log_lifetime = log_Px - log_Px1bar(:, ones(1, Npt + 1, 1), :); % Conditional distribution. Nsmpl * (Npt + 1) * Npatient matrix.
    log_lifetime = -Esum(-log_lifetime, 3) - log(Npatient); % Averaged across all patients, Nsmpl * (Npt + 1) matrix.

    if ~usealpha % In this case log_lifetime is the pdf of x1, so need this coefficient. In the other case, log_lifetime is the ddf of x1, so no need for this coefficient.
        log_lifetime = log_lifetime - log(plotx(ones(Nsmpl, 1), :)); % Factor from the chain rule going from the pdf of log(x) to the pdf of x.
        if plotx(1) == 0
            log_lifetime(:, 1) = -Inf; % Otherwise would be NaN due to Inf - Inf.
        end
    end

    [lifetime, survival, hazardrt] = get_LSH(log_lifetime, plotx, usealpha);
end

if mode == 'C'
    logplotx = log(plotx);
    coeff1 = logplotx(ones(Nsmpl, 1), :);
    coeff2 = coeff1(:, :, ones(1, 1, Npatient)); % Factor from the chain rule going from the pdf of log(x) to the pdf of x.

    log_lifetime = log_Px - log_Px1bar(:, ones(1, Npt + 1, 1), :); % Conditional distribution. Nsmpl * (Npt + 1) * Npatient matrix.

    if ~usealpha % In this case log_lifetime is the pdf of x1, so need this coefficient. In the other case, log_lifetime is the ddf of x1, so no need for this coefficient.
        log_lifetime = log_lifetime - coeff2; % No averaging.
        if plotx(1) == 0
            log_lifetime(:, 1, :) = -Inf; % Otherwise would be NaN due to Inf - Inf.
        end
    end

    % perm_log_lifetime = reshape(permute(log_lifetime, [2, 1, 3]), Npt + 1, Npatient * Nsmpl)'; % Shape into a (Npatient * Nsmpl) * (Npt + 1) 2D matrix.

    % Shape into a (Npatient * Nsmpl) * (Npt + 1) 2D matrix. Works a bit
    % faster this way on my machine.
    perm_log_lifetime = NaN(Npatient * Nsmpl, Npt + 1);
    for it_Np = 1 : Npatient
        perm_log_lifetime(Nsmpl * (it_Np - 1) + 1 : Nsmpl * it_Np, :) = log_lifetime(:, :, it_Np);
    end

    [lifetime, survival, hazardrt] = get_LSH(perm_log_lifetime, plotx, usealpha);
end

if if_centiles
    [MClifetime, MCsurvival, MChazardrt] = get_mean_and_cent(lifetime, survival, hazardrt, usealpha);
    if mode == 'C'
        log_lifetime = -Esum(-log_lifetime, 3) - log(Npatient); % Averaged across all patients, Nsmpl * (Npt + 1) matrix.

        [lifetime, survival, hazardrt] = get_LSH(log_lifetime, plotx, usealpha); % Same as the 'B' case.
    end
else
    MClifetime = [];
    MCsurvival = [];
    MChazardrt = [];
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
