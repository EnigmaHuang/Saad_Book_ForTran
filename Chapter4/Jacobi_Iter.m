function [x, converged, iter_cnt, res_norm] = Jacobi_Iter(A, b, res_tol, max_iter)
% Jacobi Iteration
	n = size(A, 1);
	
	if (nargin < 3)	res_tol  = 1e-9; end
	if (nargin < 4)	max_iter = 1000; end
	
	L = -tril(A, -1);
	U = -triu(A, 1);
	D = diag(diag(A));
	for i = 1 : n
		D(i, i) = 1.0 / D(i, i);
	end
	
	B = D * (L + U);
	g = D * b;
	
	x = zeros(n, 1);
	r = b - A * x;
	rn_stop = norm(r, 2) * res_tol;
	iter_cnt = 1;
	res_norm(iter_cnt) = norm(r, 2);
	
	converged = 0;
	while ((iter_cnt < max_iter) && (res_norm(iter_cnt) > rn_stop))
		x = B * x + g;
		r = b - A * x;
		iter_cnt = iter_cnt + 1;
		res_norm(iter_cnt) = norm(r, 2);
	end
	if (res_norm(iter_cnt) <= rn_stop) converged = 1; end
end