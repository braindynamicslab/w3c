% Calculate within attractor coordination for fixed points associated with all available G
% in a precomupted simulation file
% Note: simulation takes several hours even with parfor, and much long without parallel

clear all;
clc;

%% Set paths and variables
w3cpath = fn_w3c_setenvBox();
scriptDir = w3cpath.scriptDir;
BoxMainDir = w3cpath.BoxMainDir;

ProjMainDir = [BoxMainDir,'/YinmingSun/GitHub_data/w3c'];
outputDir = [ProjMainDir, '/sim_wAttr/'];

opt_save = 1; %Save flag

cfg.num_parc = 66;
cfg.FC_corrType = 'Spearman';

%% Load rsfMRI parcel connectivity matrix
subID = '100307';
FC_filename = [subID, '_DK62parc_conn.mat'];
FC = fn_load_DK62parc_conn(FC_filename, cfg);

%% Load saved simulation results
sim_filename = [ProjMainDir, '/RUNS/run_HCP11avg_wEE2wEI1step0_1_2021082381325/', ...
    'HCP11avg_conn_w3c_wXss.mat'];

% Note: loads ALL workspace variables from simulation output file
load(sim_filename);

%% Set additional parameters for within attractor calculation

% run simulation for 864s to match data
a.T = 864;

% noise level
a.sigma = 0.01;

% transient (s)
t_transient = 5;

%% compute within attractor coordination for fixed points of all G

within_attr_matrix_allG=[];
FCfit_rho_allG=[];
FCfit_p_allG=[];

% create cell array before parfor loop
for kG = 1:size(fixpts,1)
    within_attr_matrix_allG{kG}=[];
    FCfit_rho_allG{kG}=[];
    FCfit_p_allG{kG}=[];
end

tic
parfor kG = 1:size(fixpts,1)
    
    Nfpts = size(fixpts{kG},1);
    within_attr_matrix = nan(cfg.num_parc,cfg.num_parc,Nfpts);
    FCfit_rho = nan(Nfpts,1);
    FCfit_p = nan(Nfpts,1);

    % create a working copy of w3c object for parfor
    a_copy = a;
    
    % change the G parameter associated with working copy
    a_copy.G = Gopts(kG);
    
    for nfp = 1:Nfpts

        disp([kG, nfp]);

        a_copy = a_copy.HeunSolver(fixpts{kG}(nfp,:)');
        
        % remove transient and z-score
        within_attr_matrix(:,:,nfp) = corr(zscore(a_copy.X(a_copy.t>t_transient,1:a_copy.N)),'type',cfg.FC_corrType);
        
        [FCfit_rho(nfp),FCfit_p(nfp)] = corr(belowdiag(within_attr_matrix(:,:,nfp)), ...
            belowdiag(FC),'type',cfg.FC_corrType,'rows','pairwise');

    end
    
    within_attr_matrix_allG{kG}=within_attr_matrix;
    FCfit_rho_allG{kG}=FCfit_rho;
    FCfit_p_allG{kG}=FCfit_p;
    
end
toc

%% save results
if opt_save
    filename = [outputDir, subID, '_simFC_', cfg.FCcorr_type, '_allG.mat']
    save(filename,'subID', 'Gopts',...
    'within_attr_matrix_allG', 'FCfit_rho_allG', 'FCfit_p_allG');
end