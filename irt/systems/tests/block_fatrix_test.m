% block_fatrix_test.m
% Test the block_fatrix object

if 1 | ~isvar('A5'), printm 'A1'
	rand('state', 0)
%	A1 = Gsparse(sparse(rand(10,20)));
	A1 = rand(10,20);
	A2 = magic(20);
	A3 = rand(10,30);
	A4 = Gnufft({'epi', [5 4], [4 4], 2*[5 4]});
	A5 = rand(size(A4));
end

Ad = block_fatrix({A1, A4}); % diag
Ac = block_fatrix({A1, A2}, 'type', 'col');
Ak = block_fatrix({A1}, 'type', 'kron', 'Mkron', 2);
Ar = block_fatrix({A1, A3}, 'type', 'row');
As = block_fatrix({A4, A5}, 'type', 'sum');

if 1 % test col
	x = [1:ncol(A1)]';
	y1 = A1 * x;
	y2 = A2 * x;
	yy = Ac * x;
	printm('col forw error %g', max_percent_diff([y1; y2], yy))

	x1 = A1' * y1;
	x2 = A2' * y2;
	xx = Ac' * [y1; y2];
	printm('col back error %g', max_percent_diff(x1+x2, xx))
end

if 1 % test diag
	x1 = [1:ncol(A1)]';
	x2 = [1:ncol(A4)]';
	x = [x1; x2];
	y1 = A1 * x1;
	y2 = A4 * x2;
	yy = Ad * x;
	printm('diag forw error %g', max_percent_diff([y1; y2], yy))

	x1 = A1' * y1;
	x2 = A4' * y2;
	xx = Ad' * yy;
	printm('diag back error %g', max_percent_diff([x1; x2], xx))

	Td = build_gram(Ad, [], 0);
	y1 = Td * xx;
	y2 = [A1' * A1 * x1; A4' * (A4 * x2)];
	printm('diag gram error %g%%', max_percent_diff(y1, y2))
end

if 1 % test kron
	x1 = [1:ncol(A1)]';
	x2 = [1:ncol(A1)]';
	xx = [x1; x2];
	y1 = A1 * x1;
	y2 = A1 * x2;
	yy = Ak * xx;
	printm('kron forw error %g', max_percent_diff([y1; y2], yy))

	x1 = A1' * y1;
	x2 = A1' * y2;
	xx = Ak' * yy;
	printm('kron back error %g', max_percent_diff([x1; x2], xx))
end

if 1 % test row
	x1 = [1:ncol(A1)]';
	x2 = [1:ncol(A3)]';
	y1 = A1 * x1;
	y2 = A3 * x2;
	yy = Ar * [x1; x2];
	printm('row forw error %g', max_percent_diff(y1+y2, yy))

	x1 = A1' * yy;
	x2 = A3' * yy;
	xx = Ar' * yy;
	printm('row back error %g', max_percent_diff([x1; x2], xx))
end

if 1 % test sum
	x = [1:ncol(A4)]';
	y1 = A4 * x;
	y2 = A5 * x;
	yy = As * x;
	printm('sum forw error %g', max_percent_diff(y1+y2, yy))
	x1 = A4' * yy;
	x2 = A5' * yy;
	xx = As' * yy;
	printm('sum back error %g', max_percent_diff(x1+x2, xx))
end

tester = @(A) test_adjoint(Ac, 'complex', 1);
tester(Ac);
tester(Ad);
tester(Ak);
tester(Ar);
tester(As);
