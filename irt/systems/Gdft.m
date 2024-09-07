  function ob = Gdft(varargin)
%|function ob = Gdft([args])
%| Construct Gdft object that computes samples of the DFT of a signal
%| with dimensions [(Nd)].  This is useful for "under-sampled" MRI.
%|
%| options (at least one must be provided):
%|	'mask'	logical [(Nd)]	image-domain mask, usually: true(nx,ny)
%|	'samp'	logical [(Nd)]	which frequency-domain samples to return
%|				(default: all samples)
%|
%| out
%|	ob	[M np]	Fatrix object, where M = sum(samp(:)), np = sum(mask(:))
%|
%| See Gdft_test.m for example usage.
%|
%| Basically, you create a system matrix object by calling:
%|	A = Gdft( ... )
%| and then you can use it thereafter by typing commands like
%|	y = A * x;
%| which will auto-magically evaluate the DFT samples.
%| This is useful for iterative image reconstruction in MRI.
%|
%| Besides simple utilities like display, there are the following
%| capabilities of this object:
%|	y = A * x		forward operation
%|	x = A' * y		adjoint operation
%|	
%| Copyright 2003-6-1, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), Gdft_test, return, end

arg.mask = [];
arg.samp = [];
arg = vararg_pair(arg, varargin);

if ~islogical(arg.mask), error 'mask must be logical', end
if ~islogical(arg.samp), error 'samp must be logical', end
if isempty(arg.mask) && isempty(arg.samp), error 'must give mask or samp', end

if isempty(arg.mask)
	arg.mask = true(size(arg.samp));
elseif isempty(arg.samp)
	arg.samp = true(size(arg.mask));
else
	if ~isequal(size(arg.mask), size(arg.samp)), fail 'mismatch', end
end

arg.Nd = size(arg.mask);
arg.ndim = length(arg.Nd);
if arg.Nd(end) == 1 % 1D case
	arg.ndim = arg.ndim - 1;
end
arg.M = sum(arg.samp(:));
arg.np = sum(arg.mask(:));
arg.dim = [arg.M arg.np];

%
% build Fatrix object
%
ob = Fatrix(arg.dim, arg, ...
	'forw', @Gdft_forw, 'back', @Gdft_back, ...
	'caller', mfilename);


%
% Gdft_forw(): y = G * x
% in
%	x	[np L] or [(Nd) L]
% out
%	y	[M L]
%
function y = Gdft_forw(arg, x)

if size(x,1) == arg.np		% [np (L)]
	x = embed(x, arg.mask);	% [(Nd) (L)]
end

if isequal(size(x), size(arg.mask)) % [(Nd)]
	y = fftn(x); % [(Nd)]
	y = y(arg.samp); % [M 1]
else
	xdim = size(x);
	NN = prod(arg.Nd);
	LL = xdim(end);

	if arg.ndim == 1
		y = fft(x, [], 1); % [Nd (L)]
		y = y(arg.samp,:); % [M (L)]

	else
		y = zeros([NN LL]);
		for ll=1:LL
			tmp = stackpick(x, ll);
			y(:,ll) = col(fftn(tmp));
		end
		y = y(arg.samp,:); % [M L]
	end
end


%
% Gdft_back(): x = G' * y
% in
%	y	[M L]
% out
%	x	[np L]
%
function x = Gdft_back(arg, y)

y = embed(y, arg.samp); % [(Nd) L]
ydim = size(y);
NN = prod(arg.Nd);

if isequal(size(y), size(arg.samp)) % [(Nd)]
	x = ifftn(y); % [(Nd)] adjoint
	x = col(x); % [*Nd]
else
	if arg.ndim == 1
		x = ifft(y, [], 1); % [Nd L]
	else
		LL = ydim(end);
		x = zeros([NN LL]);
		for ll=1:LL
			tmp = stackpick(y,ll);
			x(:,ll) = col(ifft2(tmp));
		end
	end
end

% note the "NN" factor that is needed to make it an adjoint operation:
x = NN * x(arg.mask(:),:); % [np *L]

%x = reshape(x, [Ns Ld]); % [np (L)] % not needed if L can be scalar only!


function Gdft_test

Nd = [10 6];
samp = rand(Nd) > 0.3;
mask = true(Nd); mask(1) = false; % stress

A = Gdft('mask', mask, 'samp', samp);

Fatrix_test_basic(A, mask, 'complex', 1)
test_adjoint(A, 'complex', 1);

x = rand(Nd);
y = fftn(x);
y = y(samp(:));
jf_equal(y, A * x)
