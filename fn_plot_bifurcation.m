function [] = fn_plot_bifurcation(sim_filename)

% load all workspace variables from simulation output file
load(sim_filename);

%% ===== plot steady states
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

% -- plotting
dotsize =40;
figure
hold on

% non-steady states
scatter(ies(~(xssIdx | stNodeIdx | stSpiralIdx)),se(~(xssIdx | stNodeIdx | stSpiralIdx)),dotsize,'o','filled','k')

% steady states
scatter(ies(xssIdx & stNodeIdx),se(xssIdx & stNodeIdx),dotsize,'o','filled','r')
scatter(ies(xssIdx & stSpiralIdx),se(xssIdx & stSpiralIdx),dotsize,'o','filled','MarkerFaceColor','#00A651')
scatter(ies(xssIdx & unstSpiralIdx),se(xssIdx & unstSpiralIdx),dotsize,'o','filled','b')


xlabel('G', 'fontsize',20)
ylabel('$\bar{S}_E$','Interpreter','latex', 'fontsize',20)
set(gca,'fontsize',20);
% legend('others', 'stable node', 'stable spiral','unstable spiral', ...
%     'location','best','fontsize',20)
legend('others', 'stable node', 'stable spiral','limit cycle', ...
    'location','best','fontsize',20)
set(gcf, 'Position',[423 325 1178 713]);
title(sprintf('w_{EE}=%g, w_{EI}=%g',a.w_EE,a.w_EI),'fontsize',20)


end