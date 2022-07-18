function [E_max, E_mean, E_median, E_std] = fn_eGap(fixpts_attr, cfg)
% Calculate Se difference or energy gap between adjacent attractor fixed points

Se = mean(fixpts_attr(:, 1:cfg.num_parc),2);

eGap = diff(sort(Se));

if isempty(eGap)
    E_max = nan;
    E_median = nan;
    E_mean = nan;
    E_std = nan;
else
    E_max = max(eGap,[],1);
    E_mean = mean(eGap,1);
    E_median = median(eGap,1);
    E_std = std(eGap,0,1);
end

if length(E_max)>1
    error('eGap dimensions incorrect');
end

end