function centers = edge2ctr( edges )
%EDGE2CTR convert edges of histogram to centers
%   (2016-05-11, by MZ)

centers=conv(edges,[1/2 1/2],'valid');
end

