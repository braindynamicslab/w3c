% Main script for exploring w3c results with HCP11avg or one subject
clear all;
clc;

%% Set paths and variables
w3cpath = fn_w3c_setenvBox();
scriptDir = w3cpath.scriptDir;
BoxMainDir = w3cpath.BoxMainDir;

ProjMainDir = [BoxMainDir,'/YinmingSun/GitHub_data/w3c'];

%% Set processing configuration parameters
cfg = [];

% DKT atlas
cfg.num_parc = 66;
cfg.RH_parc = 1:33;
cfg.LH_parc = 34:66;

cfg.att_type = 'all'; %types of attractors to include for xAttr calc
cfg.FC_corrType = 'Spearman'; %correlation to generate model FC
cfg.corrType = 'Spearman'; %correlation between model and real FC

%% add path for GIFTI toolbox
% can be downloaded locally from: https://www.artefact.tk/software/matlab/gifti/
addpath([ProjMainDir, '/Matlab/gifti-master']);

%% Plot selected configurations on brain map
load('fsaverage_mesh_wAtlases.mat', 'cortex_327684V_wSchaeferAtlas');
cortical_surface = cortex_327684V_wSchaeferAtlas;
clear cortex_327684V_wSchaeferAtlas;
gifti_surface = gifti();

gifti_surface.faces = cortical_surface.Faces;
gifti_surface.vertices = cortical_surface.Vertices;

% Find index of selected atlas
Atlas_name = 'Desikan-Killiany';
iAtlas = find(strcmp(Atlas_name, {cortical_surface.Atlas.Name}));

Selected_Atlas = cortical_surface.Atlas(iAtlas).Scouts;
% need to ignore 2 parcels (L and R insula)

%% Load ROI mapping from cortical surface order to w3c paper order
load('DK_cortex2w3c_sort.mat');
roinames_sorted = {Selected_Atlas(w3c_sort).Label}';

[~, unsort_order] = sort(w3c_sort); %revese sort results for plotting on cortex

%% Load and plot SC
load('HCP11avg_66.mat');
SC = connAvg66;

figure, imagesc(SC);
xticks(1:cfg.num_parc);
yticks(1:cfg.num_parc);
xticklabels(roinames_sorted);
xtickangle(45);
yticklabels(roinames_sorted);
set(gcf,'Position',[742 314 1188 998]);

%% Load and plot example rsfMRI parcel connectivity matrix
subID = '100307';
FC_filename = [subID, '_DK62parc_conn.mat'];
FC = fn_load_DK62parc_conn(FC_filename, cfg);

figure, imagesc(FC);

xticks(1:cfg.num_parc);
yticks(1:cfg.num_parc);
xticklabels(roinames_sorted);
xtickangle(90);
yticklabels(roinames_sorted);
set(gca, 'TickLabelInterpreter', 'none','fontsize',12);
set(gcf,'Position',[742 314 1188 998]);

title('RSFC','fontsize',20);

%% Load saved simulation results
sim_filename = [ProjMainDir, '/RUNS/run_HCP11avg_wEE2wEI1step0_1_2021082381325/', ...
    'HCP11avg_conn_w3c_wXss.mat'];

% plot the bifurcation diagram
fn_plot_bifurcation(sim_filename);

% load bifurcation parameters
simOutputs = load(sim_filename);

% classify fixed points based on jacobian and steady state perturbation test
[idx_fixpt, idx_fixpt_overall] = fn_processFpts(sim_filename);
% 'idx_fixpt': index of fixed points for each G
%              include subfields: stNode, stSpiral, limCycl, idx_fixpt
% 'idx_fixpt_overall': index of fixed points when concatenated across all G

%% Plot xAttr matrix for selected G index
kG = 23;

% find attractor indices
[kG_attrIdx] = fn_select_attrfpts(idx_fixpt, cfg.att_type, kG);

% define attractor fixed points
fixpts_attr = simOutputs.fixpts{kG}(kG_attrIdx,:);

% calculate xAttr coordination based on seleted fixed points
[x_attr_matrix, fixpts_int] = fn_calc_xAttr(fixpts_attr, cfg);

figure, imagesc(x_attr_matrix(1:cfg.num_parc,1:cfg.num_parc));

colorbar;
caxis([0 1]);

xticks(1:cfg.num_parc);
yticks(1:cfg.num_parc);
xticklabels(roinames_sorted);
xtickangle(90);
yticklabels(roinames_sorted);
set(gcf,'Position',[742 314 1188 998]);
set(gca,'fontsize',12);

title(['G = ', num2str(simOutputs.Gopts(kG))],'fontsize',20);

%% Simulate rfMRI timeseries for ONE selected G and fixpt
% see 'calc_wAttr_coord_allG.m' for calculating result for multiple fixpts

w3c_obj = simOutputs.a; % copy object from bifurcation simulation

kG = 23;

% find attractor indices
[kG_attrIdx] = fn_select_attrfpts(idx_fixpt, cfg.att_type, kG);

% define attractor fixed points
fixpts_attr = simOutputs.fixpts{kG}(kG_attrIdx,:);

nfp_select= 5; %one randomly selected fixed point
fixpt_select = fixpts_attr(nfp_select,:);

w3c_obj.G = simOutputs.Gopts(kG);

% run simulation for 864s to match data
w3c_obj.T = 864;

% noise level
w3c_obj.sigma = 0.01;

% can take a few minutes
tic
w3c_obj = w3c_obj.HeunSolver(fixpt_select);
toc

%% plot simulated timeseries
timerange = 500:10000; %.mlx may not plot correctly with full range
figure, plot(w3c_obj.t(timerange), w3c_obj.X(timerange, 1:w3c_obj.N));
xlabel('Time (s)');
xlim([min(w3c_obj.t(timerange)) max(w3c_obj.t(timerange))]);
ylabel('Se');
title('Simulated rfMRI timeseries');

%% Calculate and plot wAttr matrix for simulated rfMRI timeseries

% transient (s)
t_transient = 5;

% remove transient and z-score
wAttr_matrix = corr(zscore(w3c_obj.X(w3c_obj.t>t_transient,1:w3c_obj.N)),'type',cfg.FC_corrType);

[FCfit_rho,FCfit_p] = corr(belowdiag(wAttr_matrix),belowdiag(FC),'type',cfg.corrType,'rows','pairwise');

figure, imagesc(wAttr_matrix);

xticks(1:cfg.num_parc);
yticks(1:cfg.num_parc);
xticklabels(roinames_sorted);
xtickangle(90);
yticklabels(roinames_sorted);
set(gca, 'TickLabelInterpreter', 'none','fontsize',12);
set(gcf,'Position',[742 314 1188 998]);

title('','fontsize',20);
title(['wAttr Conn, G = ', num2str(simOutputs.Gopts(kG))]);
    
%% Examine FC fit and energy gap across all G

% initialize loop variables
E_max = []; 
E_median = []; 
E_mean = []; 
E_std = [];

x_attr_matrix =[];

rho_all = []; p_all = [];
rho_intra = []; p_intra = [];
rho_inter = []; p_inter = [];
rho_all_condSC = []; p_all_condSC = [];
rho_intra_condSC = []; p_intra_condSC = [];
rho_inter_condSC = []; p_inter_condSC = [];

%% Loop through all G
for kG = 1:size(simOutputs.fixpts,1)
    
    % select attractors fixed points
    [kG_attrIdx] = fn_select_attrfpts(idx_fixpt, cfg.att_type, kG);
    fixpts_attr = simOutputs.fixpts{kG}(kG_attrIdx,:);
            
    % calculate energy gap properties at each G value    
    [E_max(kG), E_mean(kG), E_median(kG), E_std(kG)] = fn_eGap(fixpts_attr, cfg);
        
    % calculate xAttr matrix based on seleted fixpts and G
    [x_attr_matrix{kG,1}, fixpts_int] = fn_calc_xAttr(fixpts_attr, cfg); 
    
    if ~isempty(x_attr_matrix{kG,1})
        
        x_attr_full = x_attr_matrix{kG,1};
        
        % replace nan values with zeroes
        x_attr_full_wZeros = x_attr_full;
        x_attr_full_wZeros(isnan(x_attr_full)) = 0;
        
        % compute correlation between actual and model FC
        [rho_all(kG), p_all(kG), rho_intra(kG), p_intra(kG), rho_inter(kG), p_inter(kG)] = ...
            fn_corrModelFC(FC, x_attr_full_wZeros, cfg);
        
        % compute partial correlation between actual and model FC while controlling for structural connectivity
        [rho_all_condSC(kG),p_all_condSC(kG),rho_intra_condSC(kG),p_intra_condSC(kG),rho_inter_condSC(kG),p_inter_condSC(kG)] = ...
            fn_parcorrModelFC(FC, x_attr_full_wZeros, SC, cfg);
        
    end
    
end

%% Create subject data structure
ksub = 1;

disp(subID);
subData_all(ksub,1).subID = subID;
subData_all(ksub,1).num_parc = cfg.num_parc;
subData_all(ksub,1).parc_names = roinames_sorted;

subData_all(ksub,1).SC = SC;
subData_all(ksub,1).FC = FC;

subData_all(ksub,1).simParams.G_range = simOutputs.Gopts;
subData_all(ksub,1).simParams.G_step = simOutputs.stepsize;
subData_all(ksub,1).simParams.w_EE = simOutputs.w_EE_select;
subData_all(ksub,1).simParams.w_EI = simOutputs.w_EI_select;
subData_all(ksub,1).simParams.sigma = simOutputs.sigma0;
subData_all(ksub,1).fixpts = simOutputs.fixpts;

subData_all(ksub,1).idx_fixpt = idx_fixpt;
%     subfields: stNode, stSpiral, limCycl, idx_fixpt

% index of fixed points when concatenated across all G
subData_all(ksub,1).idx_fixpt_overall = idx_fixpt_overall;

subData_all(ksub,1).E_max = E_max;
subData_all(ksub,1).E_mean = E_mean;
subData_all(ksub,1).E_median = E_median;
subData_all(ksub,1).E_std = E_std;

subData_all(ksub).xAttrFC = x_attr_matrix;
subData_all(ksub).fitFC_xAttr_r = rho_all;
subData_all(ksub).fitFC_xAttr_p = p_all;
subData_all(ksub).fitFC_xAttrIntra_r = rho_intra;
subData_all(ksub).fitFC_xAttrIntra_p = p_intra;
subData_all(ksub).fitFC_xAttrInter_r = rho_inter;
subData_all(ksub).fitFC_xAttrInter_p = p_inter;

subData_all(ksub).fitFC_xAttr_r_condSC = rho_all_condSC;
subData_all(ksub).fitFC_xAttr_p_condSC = p_all_condSC;
subData_all(ksub).fitFC_xAttrIntra_r_condSC = rho_intra_condSC;
subData_all(ksub).fitFC_xAttrIntra_p_condSC = p_intra_condSC;
subData_all(ksub).fitFC_xAttrInter_r_condSC = rho_inter_condSC;
subData_all(ksub).fitFC_xAttrInter_p_condSC = p_inter_condSC;

%% Plot energy gap profile across all G
fn_plot_eGap(subData_all(ksub));

%% Plot FC fitting results across all G
fn_plot_fitProfile(subData_all(ksub));

%% Examine spatial distribution of fixed points
Xss_all = [];
G_all = [];

for kG = 1:size(simOutputs.Xss,1) 
    G_all = [G_all; repmat(simOutputs.Gopts(kG),size(simOutputs.Xss{kG},1),1)];
    Xss_all = [Xss_all; simOutputs.Xss{kG}];
end

figure, hold on;
subplot(1,2,1),imagesc(Xss_all);
caxis([.5 1]);
title('Parcel Se (Excitatory & Inhibitory)');
  
% Cluster and plot fixed points by distribution across all G values
num_C = 5;
% Clabel = clusterdata(Xss_all,num_C);
Clabel = kmeans(Xss_all,num_C);
subplot(1,2,2), imagesc(Clabel);
title('cluster assignment');
hold off;

%% assign cluster colors
Clabel_RGB = [];
template_RGB = fn_get_templateRGB();

for k = 1:size(Clabel,1)     
   Clabel_RGB(k,:) = template_RGB(Clabel(k),:);   
end

%% Plot fixed points with cluster labels
Se_cluster = [];
Clabel_legend = [];
figure, hold on;
hold on;
for kC = 1:num_C
    Se_cluster(kC,:) = nanmean(Xss_all(find(Clabel==kC),1:cfg.num_parc),1);
    Clabel_legend{kC,1} = ['Cluster ', num2str(kC)];        
    
    scatter(G_all(Clabel==kC),nanmean(Xss_all(Clabel==kC,1:cfg.num_parc),2),[],Clabel_RGB(Clabel==kC,:),'filled');
end
legend(Clabel_legend, 'location', 'best');
set(gca, 'fontsize',16);
set(gcf, 'Position',[423 325 1178 713]);
hold off;

%% Plot brainplot of each cluster
k_atlas = iAtlas;
for kC = 1:num_C
    atlas_feature = Se_cluster(kC,unsort_order);
    
    feature_median = nanmedian(atlas_feature);
    feature_std = nanstd(atlas_feature);
    
    cdata = zeros(size(gifti_surface.vertices,1),1); %initialized to all zeros
    
    for k_parc = 1:cfg.num_parc
        cdata(cortical_surface.Atlas(k_atlas).Scouts(k_parc).Vertices) = atlas_feature(k_parc);
    end
    
    gifti_data = gifti();
    gifti_data.cdata = cdata;
    figure,plot(gifti_surface, gifti_data);
    
%     caxis([feature_median-1.5*feature_std feature_median+1.5*feature_std]);
    caxis([0 1]);
    
    title(['Cluster ', num2str(kC)]);
    colorbar;
    view(270,90);
end

%% Plot fixed points without label (with ginput callback)
k_atlas = iAtlas;
Se_all = nanmean(Xss_all(:,1:cfg.num_parc),2);
ButtonDown_config = ['g = ginput(1);D = pdist2([G_all Se_all],g); [~,ix] = min(D);'...
    'G_select = G_all(ix); Se_select = Se_all(ix); atlas_feature = Xss_all(ix,unsort_order);'... 
    'figure_text=[''G = '', num2str(G_select), '', Se = '', num2str(Se_select)],'...
    'PlotNewGIFTI(atlas_feature, cortical_surface, k_atlas, cfg.num_parc, figure_text)'];
figure, scatter(G_all,Se_all,[],'filled',...
        'ButtonDownFcn', ButtonDown_config);

%% Plot fixed points with cluster labels (with ginput callback)
k_atlas = iAtlas;
Se_all = nanmean(Xss_all(:,1:cfg.num_parc),2);

Se_cluster = [];
Clabel_legend = [];
figure, hold on;
hold on;

for kC = 1:num_C
    Se_cluster(kC,:) = nanmean(Xss_all(find(Clabel==kC),1:cfg.num_parc),1);
    Clabel_legend{kC,1} = ['Cluster ', num2str(kC)];        
        
    ButtonDown_config = ['g = ginput(1);D = pdist2([G_all Se_all],g); [~,ix] = min(D);'...
        'G_select = G_all(ix); Se_select = Se_all(ix); atlas_feature_temp = Xss_all(ix,unsort_order);'...
        'atlas_feature = nan(1,cfg.num_parc); atlas_feature(1:6) = atlas_feature_temp(1:6); atlas_feature(9:68) = atlas_feature_temp(7:66);'...
        'figure_text=[''G = '', num2str(G_select), '', Se = '', num2str(Se_select)],'...
        'PlotNewGIFTI(atlas_feature, cortical_surface, k_atlas, cfg.num_parc, figure_text)'];
    
    scatter(G_all(Clabel==kC),nanmean(Xss_all(Clabel==kC,1:cfg.num_parc),2),[],Clabel_RGB(Clabel==kC,:),'filled',...
        'ButtonDownFcn', ButtonDown_config);
end

legend(Clabel_legend);
set(gca, 'fontsize',16);
set(gcf, 'Position',[423 325 1178 713]);
hold off;