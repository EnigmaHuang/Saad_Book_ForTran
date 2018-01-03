function x = LUSolve_CSR(n, rhs, lu_val, row_ptr, col_idx, u_ptr)
% Perform a forward and backward solve for an ILU factorization, L * U * x = rhs
% n      : Dimension of problem
% rhs    : Right hand side
% lu_val : Values of the LU matrix. L and U are stored together in the CSR format,
%          the diagonal elements of U is NOT INVERTED. In each row, the L values 
%          are followed by the U values.
% u_ptr  : Pointer to the diagonal elements in lu_val and col_idx
% row_ptr, col_idx : Standard CSR arrays indicating the start position of each 
%                    row and the column indices of each rows' nnz.

	x = zeros(n, 1);
	
	% Forward solve L * y = rhs
	for i = 1 : n
		x(i) = rhs(i);
		for k = row_ptr(i) : u_ptr(i) - 1
			x(i) = x(i) - lu_val(k) * x(col_idx(k));
		end
	end
	
	% Backward solve U * x = y
	for i = n : -1 : 1
		for k = u_ptr(i) + 1 : row_ptr(i + 1) - 1
			x(i) = x(i) - lu_val(k) * x(col_idx(k));
		end
		x(i) = x(i) / lu_val(u_ptr(i));
	end
end