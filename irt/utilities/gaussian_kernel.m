 function kern = gaussian_kernel(fwhm, nk_half)
%function kern = gaussian_kernel(fwhm, nk_half)
% samples of a gaussian kernel at [-nk_half:nk_half]
% with given FWHM in pixels
% uses integral over each sample bin so that sum is very close to unity
%
% Copyright 2001-9-18, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), return, end
if streq(fwhm, 'test'), gaussian_kernel_test, return, end
if nargin < 2, nk_half = 2 * ceil(fwhm); end

if fwhm == 0
	kern = zeros(nk_half*2+1, 1);
	kern(nk_half+1) = 1;
else
	sig = fwhm / sqrt(log(256));
	x = [-nk_half:nk_half]';
	kern = normcdf_jf(x+1/2, 0, sig) - normcdf_jf(x-1/2, 0, sig);
end

% avoid matlab's stat toolbox...
function p = normcdf_jf(x, mu, sig)
z = (x - mu) / sig;
p = 0.5 * erfc(-z ./ sqrt(2));


function gaussian_kernel_test
kern = gaussian_kernel(3);
plot(kern, '-o')
