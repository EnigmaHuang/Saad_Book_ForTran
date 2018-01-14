function [x, converged, iter_cnt, res_norm] = GS_Iter(A, b, res_tol, max_iter)
% Gauss-Seidel Iteration
	n = size(A, 1);
	
	if (nargin < 3)	res_tol  = 1e-9; end
	if (nargin < 4)	max_iter = 1000; end
	
	x = zeros(n, 1);
	r = b - A * x;
	rn_stop = norm(r, 2) * res_tol;
	iter_cnt = 1;
	res_norm(iter_cnt) = norm(r, 2);
	
	converged = 0;
	while ((iter_cnt < max_iter) && (res_norm(iter_cnt) > rn_stop))
		for i = 1 : n
			A_ii = A(i, i);
			x(i) = b(i) - A(i, :) * x + A_ii * x(i);
			x(i) = x(i) / A_ii;
		end
		r = b - A * x;
		iter_cnt = iter_cnt + 1;
		res_norm(iter_cnt) = norm(r, 2);
	end
	if (res_norm(iter_cnt) <= rn_stop) converged = 1; end
end