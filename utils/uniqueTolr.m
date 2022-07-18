function [X,idx] = uniqueTolr(X,tolr)
%UNIQUETOLR find unique elements from an array given a certain level of
%tolerance, such that any elements are close to each other within that
%tolerance are considered identical. 
%   [X,idx] = uniqueTolr(X,tolr)

[X,idx]=unique(round(X,floor(-log10(tolr))),'rows');

end

