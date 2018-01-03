function [lu_val, u_ptr, icode] = MILU0(n, row_ptr, col_idx, a_val, sum_dir)
% Modified Incomplete LU Factorization with zero fill-in,
% Using CSR format as the input and output of the matrix.
% Modified from Yousef Saad's ``Iterative Methods for Sparse Linear Systems'', 
% FORTRAN code list in section 10.3.2, variable names are changed 
% according to ALgorithm 10.4 in the same section
% row_ptr, col_idx, a_val : Standard CSR arrays, column index start from 1
% n       : Dimension of matrix
% sum_dir : 'row' for row sum compensation, 'col' for column sum compensation
% lu_val  : The CSR value array to replace a_val after ILU(0). L and U are stored 
%          together. In each row, the L values are followed by the U values. The 
%          diagonal elements of U is NOT INVERTED. 
% u_ptr   : Pointer to the diagonal elements in the value array (iu_val)
% icode   : 0 for normal return, k for zero pivot at step k

	lu_val  = a_val;
	u_ptr   = zeros(n, 1);
	row_sum = zeros(n, 1);
	col_sum = zeros(n, 1);
	col_idx_mapper = zeros(n, 1);
	
	for i = 1 : n
		j1 = row_ptr(i);
		j2 = row_ptr(i + 1) - 1;
		
		% Set mapper
		for j = j1 : j2
			col_idx_mapper(col_idx(j)) = j;
		end
		
		% ILU factorization IKJ version k-loop
		for j = j1 : j2
			k = col_idx(j);
			if (k < i)
				% lu_val(j) == a_{ik} in Algorithm 10.4
				lu_val(j) = lu_val(j) / lu_val(u_ptr(k));
				
				for jj = u_ptr(k) + 1 : row_ptr(k + 1) - 1;
					% col_idx(jj) == j in Algorithm 10.4
					jw = col_idx_mapper(col_idx(jj));
					% lu_val(jj) == a_{kj} in Algorithm 10.4
					delta = lu_val(j) * lu_val(jj);
					if (jw ~= 0)
						lu_val(jw) = lu_val(jw) - delta;
					else
						row_sum(i) = row_sum(i) + delta;
						col_sum(col_idx(jj)) = col_sum(col_idx(jj)) + delta;
					end
				end
			else
				break;
			end
		end
		u_ptr(i) = j;
		
		% Error: zero pivot
		if ((k ~= i) || (lu_val(j) == 0.0))
			icode = i;
			return;
		end
		
		% Reset mapper
		for j = j1 : j2
			col_idx_mapper(col_idx(j)) = 0;
		end
	end
	icode = 0;
	
	if (strcmp(sum_dir, 'row') == 1)
		for i = 1 : n
			lu_val(u_ptr(i)) = lu_val(u_ptr(i)) - row_sum(i);
		end
	end
	
	if (strcmp(sum_dir, 'col') == 1)
		for i = 1 : n
			lu_val(u_ptr(i)) = lu_val(u_ptr(i)) - col_sum(i);
		end
	end
end

% The output of this function is not the same as MATLAB, and 
% I think MATLAB is wrong. Consider the following matrix:
% [4 -1 -1 0; -1 4 0 -1; -1 0 4 -1; 0 -1 -1 4];
% MATLAB full LU result is:
% [4 -1 -1 0; -0.25 3.75 -0.25 -1; -0.25 -0.0667 3.7333 -1.0667; 0 -0.2667 -0.2857 3.4286]
% which shows that no fill-in in A(1, 4) and A(4, 1)
% But the MATLAB ilu() using MILU with row compensation gives:
% [4 -1 -1 0; -0.25 3.5 0 -1; -0.25 0 3.5 -1; 0 -0.2875 -0.2875 3.4286]
% which means that there is one or more fill-in in the 4th row.
% Actually there is no fill-in in the 4th row, so I believe the result
% of the last row should be [0 -0.2667 -0.2667 3.4667]