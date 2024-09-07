  function img = feldkamp(cg, ig, proj, varargin)
%|function img = feldkamp(cg, ig, proj, varargin)
%|
%| FBP reconstruction of cone-beam tomography data collected with
%| a circular source trajectory.
%| See feldkamp_example.m for example.
%|
%| in
%|	cg			ct_geom()
%|	ig			image_geom()
%|	proj	[ns,nt,na]	cone-beam projection views (line integrals)
%|
%| options
%|	'window' [npad]		'ramp' (default), or 'hann', or array.
%|				if array, then use samples [-K/2, K/2).
%|	'offset_source'		distance from isocenter to perpendicular ray
%|				[the same units (e.g., mm) as pixel_size etc.]
%|				caution: probably should not be used
%|	'ia_skip' [int]		downsample in angle to save time for tests
%|	'use_mex' 0|1		backprojector: 0 for matlab, 1 for mex (default)
%|	'nthread' 		default: jf('ncore')
%|
%| out
%|	img	[nx,ny,nz]	reconstructed image
%|
%| References: Feldkamp, Davis, Kress, JOSA-A, 1(6):612-9, June 1984.
%| Notation here follows Fessler tomography book chapter (ask if interested).
%|
%| Copyright 2004-8-28 Nicole Caparanis, Patty Laskowsky, Taka Masuda,
%| and Jeff Fessler, University of Michigan
%| arc detector case contributed by Yingying Zhang 2005-6-13

if nargin == 1 && streq(cg, 'test')
	run_mfile_local('feldkamp_example')
return
end
if nargin < 3, help(mfilename), error(mfilename), end

% defaults
arg.use_mex = 1;
arg.window = 'ramp';
arg.offset_source = 0; % r_off, distance between rotation iso-center
			% and ray from source that is orthogonal to detector.
arg.ia_skip = 1;
arg.nthread = jf('ncore');
arg = vararg_pair(arg, varargin);

if cg.pitch ~= 0 || any(cg.zshifts ~= 0)
	fail('sorry, helical CT unsupported')
end

img = feldkamp_do(proj, ...
	cg, ig, ...
	cg.ds, cg.dt, cg.offset_s, cg.offset_t, arg.offset_source, ...
	cg.dsd, cg.dso, cg.dfs, cg.orbit, cg.orbit_start, ...
	ig.mask_or, ig.nz, ig.dx, ig.dy, ig.dz, ...
	[ig.offset_x ig.offset_y ig.offset_z], ...
	arg.window, arg.ia_skip, arg.use_mex, arg.nthread);
end % feldkamp()


%
% feldkamp_do()
%
function img = feldkamp_do(proj, ...
	cg, ig, ...
	ds, dt, offset_s, offset_t, offset_source, ...
	dsd, dso, dfs, orbit, orbit_start, ...
	mask, nz, dx, dy, dz, offset_xyz, ...
	window, ia_skip, use_mex, nthread)

% step 1: weight the projections as in fan-beam case
proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, dsd, dso, dfs);

% step 2: filter the (zero padded) projections
[ns nt na] = size(proj);
if use_mex
	pad1 = 0;
else
	pad1 = 2;
end
proj = feldkamp_filter(proj, window, dsd, dfs, ds, pad1);

if pad1 % trick: zero at end saves indexing in loop
%	proj = [proj; zeros(2, nt, na)];
	proj(:,end+1,:) = 0;
end

% step 3: cone-beam backprojection of the filtered views
cpu etic
img = cbct_back(proj, cg, ig, ...
	'offset_source', offset_source, ...
	'ia_skip', ia_skip, ...
	'use_mex', use_mex, ...
	'nthread', nthread);
cpu etoc 'fdk backprojection time:'

end % feldkamp_do()


%
% feldkamp_weight1()
% step 1: weight the projections as in fan-beam case
%
function proj = feldkamp_weight1(proj, ds, dt, offset_s, offset_t, ...
	dsd, dso, dfs);
[ns nt na] = size(proj);
ss = ([-(ns-1)/2:(ns-1)/2]' - offset_s) * ds;
tt = ([-(nt-1)/2:(nt-1)/2]' - offset_t) * dt;

[ss tt] = ndgrid(ss, tt);
if isinf(dfs) % flat
	ww1 = dso * sqrt(1 + (tt/dsd).^2) ./ sqrt(dsd^2 + ss.^2 + tt.^2);
elseif dfs == 0 % arc
	ww1 = (dso/dsd) * cos(ss ./ (dsd * sqrt(1 + (tt/dsd).^2)));
else
	error 'other configurations not implemented'
end

for ia=1:na % same weighting for each view angle
	proj(:,:,ia) = proj(:,:,ia) .* ww1;
end
end % feldkamp_weight1()


%
% feldkamp_filter()
% step 2: filter the (zero padded) projections
%
function proj = feldkamp_filter(proj, window, dsd, dfs, ds, pad1)
[ns nt na] = size(proj);
npadh = 2^ceil(log2(2*ns-1));
printf('ns=%d npadh=%d', ns, npadh)

if isinf(dfs)
	H = fan_filter('flat', npadh, ds, [], window);	% [nb,1]
elseif dfs == 0
	H = fan_filter('arc', npadh, ds, dsd, window);
end
H = ds * H; % differential for discrete-space convolution vs integral

proj = ifft_sym( fft(proj, npadh, 1) .* repmat(H, [1 nt na]) );
proj = proj(1:(ns+pad1),:,:); % trick: extra zero at end saves indexing in loop
end % feldkamp_filter()


%
% fan_filter()
% apodized filter frequency response
%
function H = fan_filter(type, n, ds, dsd, window)

if streq(type, 'flat')
	h = fbp_ramp('flat', n, ds);
else
	h = fbp_ramp('arc', n, ds, dsd);
end
H = reale(fft(fftshift(h)));

if ischar(window)
	if streq(window, 'ramp')
		window = ones(n,1);
	elseif streq(window, 'hann')
		window = hann(n, 'periodic');
	else
		error 'unknown window'
	end
elseif length(window) ~= n
	error 'bad window length'
end

H = H .* fftshift(window);
end % fan_filter()
