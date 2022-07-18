function[rho_all, p_all, rho_intra, p_intra, rho_inter, p_inter] = ...
    fn_corrModelFC(FC, modelFC, cfg)

LH_parc = cfg.LH_parc; 
RH_parc = cfg.RH_parc; 
corrType = cfg.corrType;

% Extract and vectorize FC elements
[FC_vec,FC_intra_vec,FC_inter_vec] = fn_interIntraConnVec(FC, LH_parc, RH_parc);

% Extract and vectorize model FC elements
[modelFC_vec,modelFC_intra_vec,modelFC_inter_vec] = fn_interIntraConnVec(modelFC, LH_parc, RH_parc);

% Correlate between actual and model FC
[rho_all, p_all] = ...
    corr(modelFC_vec,FC_vec,'type',corrType,'rows','pairwise');

[rho_intra, p_intra] = ...
    corr(modelFC_intra_vec,FC_intra_vec,'type',corrType,'rows','pairwise');

[rho_inter, p_inter] = ...
    corr(modelFC_inter_vec,FC_inter_vec,'type',corrType,'rows','pairwise');


end