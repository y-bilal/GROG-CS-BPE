  function y = downsample3(x, m)
%|function y = downsample3(x, m)
%| downsample by averaging by integer factors
%| m can be a scalar (same factor for both dimensions)
%| or a 3-vector
%| in
%|	x	[nx ny nz]
%| out
%|	y	[nx/m ny/m nz/m]
%| function modified to 3D space by: Taka Masuda

if nargin == 1 & streq(x, 'test'), downsample3_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if length(m) == 1
	m = m * ones(ndims(x),1);
end
if length(m) ~= ndims(x), error 'bad m', end

% downsample along each dimension
y = downsample1(x, m(1));
y = downsample1(permute(y, [2 1 3]), m(2));
y = downsample1(permute(y, [3 2 1]), m(3)); % [3 1 2] order
y = permute(y, [2 3 1]);


%
% 1d down sampling of each column
%
function y = downsample1hide(x, m)
'ok'
if m == 1, y = x; return; end
n1 = floor(size(x,1) / m);
n2 = size(x,2);
n3 = size(x,3);
y = zeros(n1,n2,n3);
for ii=0:m-1
	y = y + x(ii+[1:m:m*n1],:,:);
	ticker(mfilename, ii+1, m)
end
y = y / m;


%
% downsample3_test
%
function downsample3_test
x = [6 5 2];
x = reshape(2*[1:prod(x)], x);
y = downsample3(x, 2);
jf_equal(y, [39 63; 43 67; 47 71])

if has_aspire
	filex = [test_dir filesep 'testx.fld'];
	filey = [test_dir filesep 'testy.fld'];
	fld_write(filex, x)
	delete(filey)
	com = ['op sample3 mean ' filey ' ' filex ' 2 2 2'];
	os_run(com)
	z = fld_read(filey);
	if ~isequal(y, z), error 'aspire/matlab mismatch', end
end

if 0 % big
	x = zeros(2.^[7 8 9]);
	cpu etic
	y = downsample3(x, 2);
	cpu etoc 'downsample3 time:'
end
