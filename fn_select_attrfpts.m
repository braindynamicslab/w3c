function [attrIdx] = fn_select_attrfpts(idx_fixpt, att_type, Gidx)
%fn_select_attrfpts: find attractor indices of fixpts at specified G
% Input(s):
% att_type: 'all', 'stSpiral_limCycl', 'stSpiral', 'limCycl'
% kG: selected G index
% Output(s):
% attrIdx: indices of attractors


if strcmp(att_type,'all')
    % All attractors: stable nodes, stable spirals, limit cycles
    attrIdx = [idx_fixpt.stSpiral{Gidx}; idx_fixpt.stNode{Gidx}; idx_fixpt.limCycl{Gidx}];
        
elseif strcmp(att_type,'stSpiral_limCycl')
    % Use only stable spiral and limit cycles
    attrIdx = [idx_fixpt.stSpiral{Gidx}; idx_fixpt.limCycl{Gidx}];

elseif strcmp(att_type,'stSpiral')
    % Use only stable spiral
    attrIdx = idx_fixpt.stSpiral{Gidx};
    
elseif strcmp(att_type,'limCycl')
    % Use limit cycles
    attrIdx = idx_fixpt.limCycl{Gidx};
end

end