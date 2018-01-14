function [x, converged, iter_cnt, res_norm] = Block_Jacobi_Iter(A, b, block_spos, res_tol, max_iter)
% Block Jacobi Iteration
% Diagonal block (i, i) row & col range: block_spos(i) : block_spos(i+1)-1
	n = size(A, 1);
	
	if (nargin < 4) res_tol  = 1e-9; end
	if (nargin < 5) max_iter = 1000; end
	
	nblocks = size(block_spos, 2) - 1;
	EF = -A;
	% Remove diagonal blocks
	for i = 1 : nblocks
		spos = block_spos(i);
		epos = block_spos(i + 1) - 1;
		EF(spos : epos, spos : epos) = 0;
	end
	
	x = zeros(n, 1);
	r = b - A * x;
	rn_stop = norm(r, 2) * res_tol;
	iter_cnt = 1;
	res_norm(iter_cnt) = norm(r, 2);
	
	converged = 0;
	while ((iter_cnt < max_iter) && (res_norm(iter_cnt) > rn_stop))
		x0 = x;
		
		% Notice: here "for" can be replaced by "parfor", which is the parallel block Jacobi preconditioner
		for i = 1 : nblocks
			spos = block_spos(i);
			epos = block_spos(i + 1) - 1;
			x1 = EF(spos : epos, :) * x0 + b(spos : epos);
			% Using direct solver here, in real-world applications we should use other solvers
			x(spos : epos) = A(spos : epos, spos : epos) \ x1;
		end
		
		r = b - A * x;
		iter_cnt = iter_cnt + 1;
		res_norm(iter_cnt) = norm(r, 2);
	end
	if (res_norm(iter_cnt) <= rn_stop) converged = 1; end
end