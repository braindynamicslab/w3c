function [x_attr_matrix, fixpts_int] = fn_calc_xAttr(fixpts_attr, cfg)
% fn_calc_xAttr: calculate xAttr matrix based on seleted attractor fixpts

nbins = 30;

fighist = figure;

histSe = histogram(fixpts_attr(:,1:cfg.num_parc),nbins);

min_bins = find(islocalmin(histSe.Values));

regionEdges = [0 histSe.BinEdges(min_bins) histSe.BinEdges(nbins+1)];

close(fighist);

num_states = size(fixpts_attr,1);

fixpts_int = nan(num_states, cfg.num_parc);

% assign inter bin number to each parcel value for every state
for k_state = 1:num_states

    fixpts_int(k_state,:) = discretize(fixpts_attr(k_state,1:cfg.num_parc),regionEdges);

end

x_attr_matrix = corr(fixpts_int,fixpts_int,'type',cfg.FC_corrType);


end