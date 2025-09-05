function [v, vlow, vhigh] = calc_ASI(ASI_data_file, ASI_prior_file, ASI_smpl_num, ASI_ctrl_file, ...
                                     save_dir, ASI_id, analysis_seed)
% Usage: [v, vlow, vhigh] = calc_ASI(ASI_data_file, ASI_prior_file, ASI_smpl_num, ASI_ctrl_file, ...
%                                    save_dir, ASI_id, analysis_seed)
% This function calculates the ASI of the Student mixture model. 
% We train another Stu model with the values of j(x1, x1_bar) to 
% infer the distribution of j(x1, x1_bar), and then return its posterior mean as ASI 
% and the posterior centiles as the posterior confidence limits of the ASI.

% Change Log:
%
%     1.1          12:sep:24    mx243      First version.
%     1.40         20:sep:24    rfs34      Edited with actions for mx243 after review by rfs34.
%     1.41         23:sep:24    mx243      Changes made following the
%                                          instructions. This is the last part of calc_ASI.
%     1.45         03:oct:24    rfs34      Merged with changes made by rfs34 in the meantime.
%     1.50         08:oct:24    mx243      Changed line 43 so that the
%                                          first 15% of the samples are discarded.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.calc_ASI.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 7
    analysis_seed = [];
end

smpl_file = sample_Stu_mix(ASI_data_file, ASI_prior_file, ASI_smpl_num, save_dir, ASI_id, ...
                           ASI_ctrl_file, [], analysis_seed, 0, 1);

load(smpl_file, 'Samples', 'LN0', 'LN1');
smpl_num = LN1 - LN0 + 1;
Phi = Samples.phi(LN0 : LN1);
C = Samples.C(LN0 : LN1);
p_dis = 0.15; % Discarding the first 15% of samples.

Num = ceil((1 - p_dis) * (LN1 - LN0));
vsamples = NaN(Num, 1);
cnt = 0;
for num = smpl_num - Num + 1 : smpl_num,
    vsample = 0;
    for it_c = 1 : C(num)
        vsample = vsample + Phi(num).p(it_c) * ...
                            mean_skewStu(Phi(num).mu(it_c), Phi(num).nu(it_c), ... % mu, nu are 1D.
                                         Phi(num).m(it_c), Phi(num).m(it_c) - 1); % r = m - 1 in the model for ASI.
    end
    cnt = cnt + 1;
    vsamples(cnt) = vsample;
end
v = mean(vsamples, 1);

vsamples = sort(vsamples);
nlow = ceil(0.025 * Num);
nhigh = floor(0.975 * Num);
vlow = vsamples(nlow);
vhigh = vsamples(nhigh);

M = matfile(smpl_file, 'Writable', true);
M.ASI = v; % Record ASI.
M.ASI025 = vlow;
M.ASI975 = vhigh;
M.ASImeansamples = vsamples; % For calculating other quantiles of the posterior distribution of the mean ASI.

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
