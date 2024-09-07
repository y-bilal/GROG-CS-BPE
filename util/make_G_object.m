function G=make_G_object(kx,ky,N)

% For ESMRMB Parallel Imaging Workshop 2010 in Wuerzburg, Germany
% Written by Nicole Seiberlich on 06/23/2010
% Please contact nes30@case.edu with comments, questions, suggestions.
% Uses Fessler toolbox :  http://www.eecs.umich.edu/~fessler/code/


kspace=double([ky(:) kx(:)]);
size(kspace);
mask=true([N N]);
size(mask);
Nm=size(mask);
nufft_args = {Nm, [3 3], 2*Nm, Nm/2, 'table', 2^12, 'minmax:kb'};
G = Gmri(kspace, mask, 'fov', 1, 'basis', {'rect'}, 'nufft', nufft_args);
