function [lifetime, survival, hazardrt, MClifetime, MCsurvival, MChazardrt] = ...
         get_prediction_use_alpha(log_Px, log_ddf, log_Px1bar, plotx, mode, if_centiles, usealpha)
% Usage: [lifetime, survival, hazardrt, MClifetime, MCsurvival, MChazardrt] = ...
%        get_prediction_use_alpha(log_Px, log_ddf, log_Px1bar, plotx, mode, if_centiles, usealpha)
% This function returns different typs of the predictions of the three 
% values, and optionally their means and centiles, when the survival data
% is obtained using Monte-Carlo integration over alpha.
% log_ddf is the decumulative distribution of the correspoding log_Px, see
% get_probs.
% This function calls get_prediction to avoid writing repeated codes, but 
% in a fairly ugly way. See comments in get_prediction for detail.


% Change Log:
%
%     1.1          17:oct:24    mx243      First version.
%     1.61         26:oct:24    rfs34      Generalised usealpha to be what used to be nalpha.
%     1.71         25:jan:25    rfs34      Comments only changed.
%     1.75         19:mar:25    mx243      Now allows for general plotx.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.get_prediction_use_alpha.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if ~usealpha
    error('get_prediction_use_alpha called with usealpha = 0');
end

% plotx = 0 : T_max / Npt : T_max;

% Get lifetime prediction by calling get_prediction normally, the trash
% bits will be survival data from numerical integration and the corresponding hazard rate.
[lifetime, trash, trash, MClifetime, trash, trash] = ...
get_prediction(log_Px, log_Px1bar, plotx, mode, if_centiles, 0);

% Get survival data by putting log_ddf as the first argument and usealpha > 0 or < 0. 
% The trash bits will be empty.
[survival, trash, trash, MCsurvival, trash, trash] = ...
get_prediction(log_ddf, log_Px1bar, plotx, mode, if_centiles, usealpha);

hazardrt = - gradient(survival, plotx, 1) ./ survival; % Numerical differentiation.

if if_centiles
    % Below copied from get_mean_and_cent, used to calculate the mean and
    % centiles of hazard rate.
    Ncurves = size(lifetime, 1);
    LB = ceil(0.025 * Ncurves);
    UB = floor(0.975 * Ncurves);

    hazardrt = sort(hazardrt, 1);
    MChazardrt.mean = mean(hazardrt, 1);
    MChazardrt.p025 = hazardrt(LB, :);
    MChazardrt.p975 = hazardrt(UB, :);
else
    MChazardrt = [];
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
