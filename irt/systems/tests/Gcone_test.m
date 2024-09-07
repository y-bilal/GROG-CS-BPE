% Gcone_test.m
% test the Gcone object
% Copyright 2008-1-1, Jeff Fessler, University of Michigan

%ptypes = {'nn1', 'pd1', 'sf1'};
ptypes = {'sf1'};
if 1 && exist('dd_ge1_mex') == 3
	ptypes{end+1} = 'dd1'; % UM only
	ptypes{end+1} = 'dd2'; % UM only
end
nn = length(ptypes);

% small systems for basic tests
if 0 || ~isvar('A1'), printm 'setup small'
	f.down = 16;
	igs = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 1, 'dz', 0.5, ...
		'offset_x', 2.9, 'offset_y', 3.5, 'offset_z', -3.4, ...
		'down', f.down);

	dfs_list = [0 inf inf]; % arc flat parallel
	dsd_list = [949.075 949.075 inf]; % arc flat parallel
	for kk = 1:length(dfs_list)
		cgs = ct_geom('ge1', 'nt', 320, 'dt', +1.0239, ...
       			'pitch', 0, ... % todo: test helix later
			'dfs', dfs_list(kk), ... % arc or flat
			'dsd', dsd_list(kk), ... % fan or parallel beam
			'down', f.down);
		if im, cgs.plot(igs); end

		clear A1 Ac
		for ii=1:nn
			ptype = ptypes{ii};
			if (streq(ptype, 'dd1') || streq(ptype, 'dd2')) ...
				&& isinf(cgs.dsd)
				A1{ii} = [];
				Ac{ii} = [];
				continue
			end
			A1{ii} = Gcone(cgs, igs, 'type', ptype, 'nthread', 1);
			Ac{ii} = Gcone(cgs, igs, 'type', ptype);
		end

		printm 'test small'
		for ii=1:nn
			printm('testing type %s dfs=%g dsd=%g', ...
				ptypes{ii}, cgs.dfs, cgs.dsd)
			if isempty(A1{ii}), continue, end
			tester_tomo2(A1{ii}, igs.mask, 'G2', Ac{ii}) % paces
			test_adjoint(A1{ii}, 'big', 1, 'tol', 5e-5)
			test_adjoint(Ac{ii}, 'big', 1, 'tol', 5e-5)
		end
	end
end

if ~isvar('x0'), printm 'x0 big'
	f.down = 4;
	igb = image_geom('nx', 512, 'ny', 480, 'nz', 416, ...
		'dx', 1, 'dz', 0.5, ...
		'offset_x', 12.9, 'offset_y', 3.5, 'offset_z', -3.4, ...
		'down', f.down);
	ell = [3*igb.dx 5*igb.dx -2*igb.dz ...
		igb.dx*igb.nx/3 igb.dy*igb.ny/4 igb.zfov/4 ...
		0 0 10];
	x0 = ellipsoid_im(igb, ell, 'oversample', 2);
end

% big systems for accuracy tests
if ~isvar('Ab'), printm 'setup big'
	cgb = ct_geom('ge1', 'nt', 320, 'dt', +1.0239, ...
       		'pitch', 0, ... % todo: test helix later
		'down', f.down);
%		'dfs', inf, ... % flat detector
%		'dsd', inf, 'dfs', inf, ... % parallel beam
%		'na', 1, 'orbit_start', 13, ...

	clear Ab Ah
	for ii=1:nn
		ptype = ptypes{ii};
		if streq(ptype, 'dd', 2) && isinf(cgb.dsd)
			Ab{ii} = [];
%			Ah{ii} = [];
			continue
		end
		Ab{ii} = Gcone(cgb, igb, 'type', ptype);
%		Ah{ii} = Gcone(cgb, igb, 'type', ptype, 'use_hct2', true);
	end
end

if ~isvar('ya'), printm 'analytical projections'
	ya = ellipsoid_proj(cgb, ell, 'oversample', 2);
	im clf, im(ya)
prompt
end

if ~isvar('yb'), printm 'discrete projections'
	nrmse = @(x,y) norm(y(:)-x(:)) / norm(x(:)) * 100;
	for ii=1:nn
		if ~isempty(Ab{ii})
			cpu etic
			yb{ii} = Ab{ii} * x0;
			f.time(ii) = cpu('etoc');
			printm('nrmse %s %g %%', ptypes{ii}, nrmse(ya, yb{ii}))
		else
			f.time(ii) = 0;
			yb{ii} = cgb.zeros;
		end
	end
end

if 0, printm 'look at error in worst views'
	im('pl', 2, nn)
	for ii=1:length(ptypes)
		err = yb{ii} - ya;
		tmp = reshape(err, [], cgb.na);
		tmp = sum(tmp.^2); % error in each view
		ia = imax(tmp); % worst view
		im(ii, err(:,:,ia)), cbar h
		titlef('%s ia=%d', ptypes{ii}, ia)
		im('subplot', ii+nn)
		plot(tmp), axis tight
	end
return
end

if 1, printm 'projection profiles'
	it = cgb.nt;
	it = round(cgb.nt/2); % middle
	it = it + [-2 0 2];
	ia = imin(abs(cgb.ad - 45));
	ia = ceil(ia/2);
%	ia = 1;
	pro = @(y) col(y(:,it,ia));
	arg = [];
	for ii=1:length(ptypes)
		arg = [arg pro(yb{ii})];
	end
	if im
		clf, plot([arg pro(ya)])
		text(10, 200, sprintf('ia=%d', ia))
		text(10, 400, sprintf('ang=%g', cgb.ad(ia)))
		legend(ptypes{:}, 'true')
		axisy(0, 1.2 * max(ya(:)))
		grid
	end
return
end


if 0 % explore scale factor
	pr [igb.dx igb.dz cgb.ds cgb.dt]
	pr rad2deg(atan(cgb.dt * cgb.nt/2 / cgb.dsd)) % cone angle
	for ii=1:length(ptypes)
		tmp = sum(col(yb{ii} .* ya)) / sum(ya(:).^2);
		printm('scale needed %s %g', ptypes{ii}, 1/tmp)
	end
return
end

if 0 % check shifts for yb from Gcone - ok
	ia = [1 17 100 cgb.na];
	ii = 4
	im_toggle(yb{ii}(:,:,ia), ya(:,:,ia))
return
end

if 0 && jf('isum')
	if ~isvar('yh'), printm 'hct2 projections'
		for ii=1:nn
			yh{ii} = Ah{ii} * x0;
			if streq(ptypes{ii}, 'nn1') || streq(ptypes{ii}, 'pd1')
				thresh = 2e-2; % todo: why big?
			else
				thresh = 7e-6;
			end
			equivs(yb{ii}, yh{ii}, 'thresh', thresh)
		end
	end

	if 1 % check shifts for hct2: all OK
		ii = 4
		if 0 % hct2 vs mat
			tmp = yh{ii} - yb{ii};
			ia = sum(reshapee(abs(tmp), [], cgb.na)) ~= 0;
			im(tmp(:,:,ia))
		end

		ia = [1 17 100 cgb.na];
		im_toggle(yh{ii}(:,:,ia), yb{ii}(:,:,ia))
%		im_toggle(yh{ii}(:,:,ia), ya(:,:,ia))
	end
end

if 0 % check shifts - ok
	ia = [1 17 100 cgb.na];
	im_toggle(ya(:,:,ia), yb{4}(:,:,ia))
return
end

if 0 % dd1 vs dd2 - they match well
	i_dd1 = 3; i_dd2 = 4;
	if streq(ptypes{i_dd1}, 'dd1') && streq(ptypes{i_dd2}, 'dd2')
		im clf, im(yb{i_dd1} - yb{i_dd2}), cbar
		equivs(yb{i_dd1}, yb{i_dd2}, 'thresh', 9.8e-6)
	else
		fail 'bug'
	end
return
end

if 1, printm 'show projections and differences'
	im clf, im('pl', 2, 1+nn)
	im(1, x0)
	ia = round([1 cgb.na/4 cgb.na/2 cgb.na]);
	im(nn+2, ya(:,:,ia))

	for ii=1:nn
		tmp = yb{ii};
		im(ii+1, tmp(:,:,ia))
		xlabel(ptypes{ii})
		im(ii+2+nn, tmp(:,:,ia)-ya(:,:,ia))
	end
prompt
end

if 1, printm 'show back-projections'
	im clf, im('pl', 1, nn)
	iz = round([1 igb.nz/4 igb.nz/2 igb.nz]);
	iz = 1:2:igb.nz;
	for ii=1:nn
%		tmp = cgb.ones;
		tmp = cgb.zeros; tmp(:,:,20) = 1;
		tmp = Ab{ii}' * tmp;
		im(ii, tmp(:,:,iz))
		xlabel(ptypes{ii})
	end
end

if 0
	im pl 1 3
	ia = round([1 cg.na/4 cg.na/2 cg.na]);
	im row 4
	im(1, ya(:,:,ia));
	im(2, yc(:,:,ia));
	im(3, yc(:,:,ia)-ya(:,:,ia));
	im reset
end
