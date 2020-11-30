function [window] = window1(sigma,filtersize)
ind = -floor(filtersize/2) : floor(filtersize/2);
[X Y] = meshgrid(ind, ind);

%// Gaussian Mask
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma))./(2*pi);

%// Normalize 
h = h / sum(h(:));
window=h;
end

