function [FC, DK62labels_sorted] = fn_load_DK62parc_conn(filepath, cfg)

num_parc = cfg.num_parc;

% sort order for rsfMRI
load('rsfMRI_DK62parc_order.mat'); %DK62parc_order

load(filepath); %DK62parc_conn, DK62parc_labels

FC = nan(num_parc,num_parc);
DK62labels_sorted = cell(66,1); %for double check of sort order

for k1 = 1:num_parc
    if ~isnan(DK62parc_order(k1))
        DK62labels_sorted(k1) = DK62parc_labels(DK62parc_order(k1));
        for k2 = 1:num_parc
            if ~isnan(DK62parc_order(k2))
                FC(k1,k2)=DK62parc_conn(DK62parc_order(k1),DK62parc_order(k2));
            end
        end
    end
end

end