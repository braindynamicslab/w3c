%classifying fixpts and concatentate across all attractors
function [idx_fixpt, idx_fixpt_overall] = fn_processFpts(sim_filepath)

load(sim_filepath);

tolr = 0.05; % tolerance for considering two points are the same
tolr_osc = 0.2; % use a different asymmetric tolerance for limit cycle
% -- find which fixed points are at steady state
fp = cell2mat(fixpts);
ev = cell2mat(Jeigs);
se = mean(fp(:,1:a.N),2);

% -- control parameter
Npts = cell2mat(cellfun(@(x) size(x,1),fixpts,'uniformoutput',0));
ies = cell2mat(cellfun(@(x,y) repmat(x,y,1), num2cell(Gopts'),num2cell(Npts),'UniformOutput',0));

% -- classification
tolrJeig = 1e-4; % set tolerance for eigenvalues
[stNodeIdx, unstNodeIdx, stSpiralIdx, sdlSpiralIdx, unstSpiralIdx] = classifyFpts(ev,tolrJeig);

% -- compute steady state
xss = cell2mat(Xss);
dfpXss = max(abs(fp-xss),[],2);
xssIdx = dfpXss < tolr;% stable node tolerance
xssIdx(unstSpiralIdx) = xssIdx(unstSpiralIdx) | dfpXss(unstSpiralIdx)<tolr_osc;% limit-cycle tolerance

%% Adjust classification test result based on additional steady state perturbation test
fixpt_test = [];
fixpt_test.stNode = xssIdx & stNodeIdx;
fixpt_test.stSpiral = (xssIdx & stSpiralIdx);
fixpt_test.limCycl = (xssIdx & unstSpiralIdx);
fixpt_test.other = ~(xssIdx | stNodeIdx | stSpiralIdx);

%% create stacked fix points with overall index
overall_count = 0;
idx_fixpt_overall = [];
for k1 = 1:size(fixpts,1)
    for k2 = 1:size(fixpts{k1},1)
        overall_count = overall_count + 1;
        idx_fixpt_overall{k1,1}(k2) = overall_count;
    end
end

%% Create separate indices for each G
idx_fixpt = [];
for kG = 1:size(fixpts,1)

    idx_fixpt.stNode{kG} = find(fixpt_test.stNode(idx_fixpt_overall{kG}));
    idx_fixpt.stSpiral{kG} = find(fixpt_test.stSpiral(idx_fixpt_overall{kG}));
    idx_fixpt.limCycl{kG} = find(fixpt_test.limCycl(idx_fixpt_overall{kG}));
    idx_fixpt.other{kG} = find(fixpt_test.other(idx_fixpt_overall{kG}));    
    
end
    
end