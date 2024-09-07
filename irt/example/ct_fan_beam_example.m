% ct_fan_beam_example.m
% compare FBP and iterative reconstruction for a 2D fan-beam CT problem
% Copyright 2009-10-05, Jeff Fessler, University of Michigan

if ~isvar('sino'), printm 'setup geometry, image, sinogram'
	down = 4;
	ig = image_geom('nx', 512, 'fov', 500, 'down', down);
	sg = sino_geom('ge1', 'down', down);

	% read image
	ddir = path_find_dir([filesep 'data']);
	xtrue256 = fld_read([ddir filesep 'ncat,256,slice,140,ct,x100.fld']);

	if 1 % more realistic sinogram from finer image
		ig_big = image_geom('nx', 512, 'fov', 500, 'down', 2);
		Abig = Gtomo2_dscmex(sg, ig_big);
		sino = Abig * xtrue256;
	end
	xtrue = downsample2(xtrue256, 2);

	% system object
	A = Gtomo2_dscmex(sg, ig);

	im clf, im pl 2 2
	clim = [0 200];
	im(1, xtrue, 'x', clim), cbar
	im(2, sino, 'sino'), cbar

	clear ddir ig_big Abig
prompt
end

if ~isvar('fbp'), printm 'fbp 2d fan-beam reconstruction'
	tmp = fbp2(sg, ig);
	fbp = fbp2(sino, tmp);
	im(3, fbp, 'FBP', clim), cbar
return
end

% todo: iterative reconstruction
