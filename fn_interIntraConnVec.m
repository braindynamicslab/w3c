function [full_vec,intra_vec,inter_vec] = fn_interIntraConnVec(conn_matrix, LH_parc, RH_parc)
    full_vec = belowdiag(conn_matrix);

    intra_vec = [belowdiag(conn_matrix(RH_parc,RH_parc));...
        belowdiag(conn_matrix(LH_parc,LH_parc))];

    inter_vec = reshape(conn_matrix(LH_parc,RH_parc),...
        length(LH_parc)*length(RH_parc),1);
end