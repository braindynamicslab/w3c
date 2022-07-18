function [fixpts_discret,state_types,state_max] = discretefp(fixpts,nbins,varargin)
%DISCRETEFP discretize fixed points using a histogram. The coordinates will
%become integers corresponds to the smallest to the largest local peak in
%the histogram, separated by neighboring troughs. 
%   [fixpts_discret,state_type,state_max] = discretefp(fixpts,nbins,...)
% input:
%   fixpts: N-by-D matrix (or a cell array of such matrices). N is the
%   number of fixed points (in each cell), D is the dimension of the state
%   space.
%   nbins: number of bins in the histogram, an integer. 
% output:
%   fixpts_discret: N-by-d matrix (or a cell array of such matrices). N is
%   still the number of fixed points, d is dimension of the discretized
%   state (d<=D). Discretization can use a d-dimensional subspace of the
%   original state space, which is defined by parameter `dimIdx` below. 
%   state_type: Ns-by-1 vector of integer indices from 1 to Ns, where Ns is
%   the number of discrete states.
%   state_max: Ns.
% parameters:
%   dimIdx: d-by-1 vector, indices of the columns of matrix `fixpts` to be
%   used to compute the discretized states.
%   binlims: the [min, max] values defining the boundary of the histogram.
%   plot: whether or not to plot the histogram, and boundary of
%   discretization. 
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 07-08-2020 ~
%}

% -- aggregate all states
if iscell(fixpts)
    allfp=cell2mat(fixpts);
else
    allfp = fixpts;
end
ND = size(fixpts,2);

% -- options
p = inputParser;
p.addParameter('dimIdx',1:ND) % indices of the dimensions to include
p.addParameter('binlims',[min(allfp(:)), max(allfp(:))])% [min max] of bins
p.addParameter('plot',true) % plot histogram (for diagnostics)
p.addParameter('minprominence',0) % minimal prominence of the troughs.
p.parse(varargin{:})
par=p.Results;


allse=allfp(:,par.dimIdx);% state in selected dimensions

% -- histogram
binedges = linspace(par.binlims(1),par.binlims(2),nbins+1);
ctrs = mean([binedges(1:end-1);binedges(2:end)]);
stateCount = histcounts(allse(:),binedges,'Normalization','pdf');
figure
plot(ctrs,stateCount); hold on
xlabel('continuous state')
ylabel('probability density')

% -- find classification boundaries
[~, bdr_idx] = findpeaks(-stateCount,'MinPeakProminence',par.minprominence);
bdr = ctrs(bdr_idx);
arrayfun(@(x) plot([x x],ylim,'--k','linewidth',2),bdr);

% -- discretize state variables
if iscell(fixpts)
    fixpts_discret = cellfun(@(x) discretize(x(:,par.dimIdx),[par.binlims(1) bdr par.binlims(2)]),fixpts,'UniformOutput',0);
    state_types = unique(cell2mat(fixpts_discret));
else
    fixpts_discret = discretize(fixpts(:,par.dimIdx),[par.binlims(1) bdr par.binlims(2)]);
    state_types = unique(fixpts_discret);
end

state_max = max(state_types);
end

