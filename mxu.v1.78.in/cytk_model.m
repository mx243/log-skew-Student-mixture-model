function [cytk_mean, cytk_mean_smpls] = cytk_model(log_cytk, y_bar, sigma)
% Usage: [cytk_mean, cytk_mean_smpls] = cytk_model(log_cytk, y_bar, sigma)
% This function returns 3400 posterior samples of the mean log cytokine
% value in a certain group of patients. 
% log_cytk is a row vector of the log cytokine values of the patients 
% we want to investigate. 
% y_bar, sigma are the mean and standard deviation of log cytokine across 
% all patients in the dataset.
% Note a .mat file containing the MCMC samples will be stored as
% ASI_smpls.cytk.mat.

% Change Log:
%
%     1.1          06:jun:25    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.cytk_model.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

z = (log_cytk - y_bar) / sigma;
smpl_num = 4e3;

calc_ASI(z, 'cytokine_Whitworth_prior.mat', smpl_num, 'sa_ctrl_Stu_mix.rv0.mat', '.', 'cytk', []);

load('ASI_smpls.cytk.mat', 'ASI', 'ASImeansamples');

cytk_mean = ASI * sigma + y_bar;
cytk_mean_smpls = ASImeansamples * sigma + y_bar;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
