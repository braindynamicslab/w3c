% bifurcation_global_model.m
% ------------------------
% systematically study the bifurcation of the global WWWC model (Zhang, Sun &
% Saggar, 2022, NeuroImage). This program could take more than a week to run. If you
% simply want to get the gist of the computation of bifurcation diagrams.
% Please see the computation for the local model. 
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 01-13-2020 ~
%}
% Please cite: Zhang, M., Sun, Y., & Saggar, M. (2022). Cross-attractor repertoire provides new 
% perspective on structure-function relationship in the brain. NeuroImage, 119401.

%%
clear 
close all
clc

%% set paths based on configured default paths
% w3cpath = fn_w3c_setenvBox();
% scriptDir = w3cpath.scriptDir;
% BoxMainDir = w3cpath.BoxMainDir;
ProjMainDir = '/Users/saggar/Dropbox/matlab_code/w3c-main-ys-dec22_2021';%[BoxMainDir,'/YinmingSun/GitHub_data/w3c'];
addpath models/
addpath utils/

DWIinputDir = 'data/';
% DWIinputDir = [ProjMainDir, '/LA5C_subDWI_DKT/'];

simOutputDir = [ProjMainDir, '/RUNS/'];

%% Select SC matrix and set configurations

% Group SC
project_str = 'HCP11avg';
population_str = 'hlt';
C_name = 'connAvg66';
conn_filename = 'HCP11avg_66.mat';

% % Individual SC
% project_str = 'ucla_la5c';
% sub_str = 'sub-50004';
% population_str = 'scz';
% C_name = 'DWI_infNorm';
% conn_filename = [sub_str, '_DWI_infNorm.mat'];

load([DWIinputDir, conn_filename]);
eval(['C = ', C_name, ';']);

%% Set variable simulation parameters 
w_EE_select = 2;
w_EI_select = 1;
G1 = 1.7; %start G sweep
G2 = 3; %end G sweep
stepsize = 0.5; %G sweep stepsize

%% create output folder based on variable simulation parameters

config_str = ['wEE',num2str(w_EE_select),'wEI',num2str(w_EI_select),'step',num2str(stepsize)];
config_str = strrep(config_str,'.','_');

datetime_str = char(datetime('now','Format','yyyyMMDDHHmm'));
foldername = ['run_',project_str, '_', population_str,'_',config_str,'_',datetime_str];

folderpath = [simOutputDir, foldername];

mkdir(folderpath);

%% -- set up model

a = WWWC(C);

% -- setting fixed model parameters
a.w_EE = w_EE_select;
a.w_EI = w_EI_select;
a.w_IE = a.w_EE;

% -- simulation parameters
a.dt = 0.001;
a.T = 10;
t_trans = 5;

% -- tolerance for difference between fixed points
tolrfpts = 1e-4;

% -- initial conditions for simulations (for finding seed steady states)
gridreso=0.1;
[XX,YY]=meshgrid(0:gridreso:1,0:gridreso:1);
ICs=repelem([XX(:),YY(:)],1,a.N);
Nics = size(ICs,1);

% -- use fsolve
opt=optimset('fsolve');
opt.Display = 'none';

% -- define range of G
% Gopts_stepsize = stepsize;
% Gopts = 0:Gopts_stepsize:5; % !!step size must be small enough b/c using fixed points from previous iteration
% N_Gs = length(Gopts);

Gopts = G1:stepsize:G2;
N_Gs = length(Gopts);

% -- preallocate memory
[fixpts,Jeigs]=deal(cell(N_Gs,1));
fpts_past = nan(0,a.N*2);

%% -- plotting fixed points on the flight
figure('position',[10 10 1000 1000])
hold on
t0=tic;
for ng=1:N_Gs
    a.G = Gopts(ng);
    disp(['G=' num2str(Gopts(ng))])
    ttmp=tic;
    Xss = nan(a.N*2,Nics);
    % -- compute steady state as initial seed for finding fixed points   
    for nic=1:Nics
        a = a.HeunSolver(ICs(nic,:)');
        Xss(:,nic)=mean(a.X(a.t>t_trans,:))';
    end
    Xss = [Xss';fpts_past];
    [~,seOrder] = sort(mean(Xss(:,1:a.N),2));% sort solutions by average S_E
    Xss = uniqueTolr(Xss(seOrder,:),tolrfpts);
    
    % -- find fixed points recursively
    f=@(x) a.D([],x);
    [fixpts{ng,1},Jeigs{ng,1}] = recurfindFpts(f,Xss,'fsolveOpt',opt,'maxFpts',200);
    toc(ttmp)
    
    % -- remember past fixed points for reuse
    fpts_past = fixpts{ng,1};
    
    % -- plotting the solutions
    Nfpts = size(fixpts{ng,1},1);
    unstabIdx = any(real(Jeigs{ng,1})>0,2); % unstable and saddle nodes are plotted in blue
    cmap = repmat([1 0 0],Nfpts,1);
    cmap(unstabIdx,:)= repmat([0 0 1], sum(unstabIdx),1);
    scatter(repmat(a.G(0),Nfpts,1),mean(fixpts{ng,1}(:,1:a.N),2),[],cmap,'filled');
    drawnow
    save([folderpath sprintf('/fixed_points_cont_%02d', ng)]);
    
    % remove previous iteration of fixed points to save storage space
    file_prev = [folderpath sprintf('/fixed_points_cont_%02d.mat', ng-1)];
    if exist(file_prev, 'file')==2
        delete(file_prev);
    end
    
end
toc(t0)

% save fix points
save([folderpath, '/', project_str, '_', population_str, '_conn_w3c']);

%% ===== compute steady states (attractors) ===== %%
% -- define transient
t_transient = 5;
% -- define perturbation on initial condition
sigma0=0.5;

% -- compute steady state following small perturbation from fixed points
NGs = length(Gopts);
Xss = cell(size(fixpts));
figure;
for ng = 1:NGs
    disp(Gopts(ng))
    a.G = Gopts(ng);
    Nfpts = size(fixpts{ng},1);
    Xss{ng} = nan(size(fixpts{ng}));
    tic
    for nfp = 1:Nfpts
        % -- simulation
        a = a.HeunSolver(fixpts{ng}(nfp,:)'+sigma0*rand(size(a.X_0)));
        % -- compute steady state
        Xss{ng}(nfp,:) = mean(a.X(a.t>t_transient,:),1)';
        plot(a.X);
        pause()
        clf;
    end
    toc
end

% save fix points with steadystate info
save([folderpath, '/', project_str, '_', population_str, '_conn_w3c_wXss']);

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

%% Plot fixed points in separate subplots based on type 
figure, hold on;

% non-steady states
subplot(4,1,1), scatter(ies(~(xssIdx | stNodeIdx | stSpiralIdx)),se(~(xssIdx | stNodeIdx | stSpiralIdx)),dotsize,'o','filled','k');
title('non-steady states');

% stable nodes
subplot(4,1,2), scatter(ies(xssIdx & stNodeIdx),se(xssIdx & stNodeIdx),dotsize,'o','filled','r');
title('stable nodes');

% stable spirals
subplot(4,1,3), scatter(ies(xssIdx & stSpiralIdx),se(xssIdx & stSpiralIdx),dotsize,'o','filled','MarkerFaceColor','#00A651');
title('stable spirals');

% unstable spirals
subplot(4,1,4), scatter(ies(xssIdx & unstSpiralIdx),se(xssIdx & unstSpiralIdx),dotsize,'o','filled','b');
title('limit cycles');

hold off;
