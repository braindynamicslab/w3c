function [cmap] = rankcolor(cmap,val)
%RANKCOLOR re-order a color map based on the rank of the values in a vector
%   [cmap] = rankcolor(cmap,val)
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 08-09-2019 ~
%}
[~,rank]= sort(val); % rank by degree
cmap(rank,:) = cmap;
end

