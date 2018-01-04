function [lu_val, lu_index, lu_uptr, ierr] = ILUT(n, val, row_ptr, col_idx, max_row_nnz, drop_tol, lu_buff_size)
% Incomplete LU factorization with dual truncation mechanism
% Translating Prof. Saad's fortran 77 code in SPARSKITv2 to MATLAB,
% should be easy to translate into C (indices need to start from 0).
% ===== Input parameters =====
% n : Row dimension of A
% val, row_ptr, col_idx : Standard CSR arrays storing A
% max_row_nnz  : Fill-in paramater, the maximum number of nnz each row of L and U 
%               can have (excluding diagoal), must be >= 0
% drop_tol     : Dropping threshold
% lu_buff_size : The length of arrays lu_val and lu_index. If the arrays are not
%                big enough ILUT will stop with an error message
% =====  Output values   =====
% lu_val   : Matrix stored in Modified Sparse Row (MSR) format containing
%             the L and U values together. Diagonal elements are stored in 
%             lu_val(1:n). Each row of the lu_val contains the i-th row of L
%             followed by the i-th row of U.
% lu_index : (1 : n) is the counter part of row_ptr, 
%             (n + 2 : end) is the counter part of col_idx
% lu_uptr  : Pointer to the beginning of U part in each row
% ierr     : Error message: 
%             >0 --> Zero pivot encountered at step number ierr
%              0 --> Success
%             -1 --> Error, input matrix may be wrong
%             -2 --> The matrix L overflows the array lu_val
%             -3 --> The matrix U overflows the array lu_val
%             -4 --> Illegal value for max_row_nnz
%             -5 --> Zero row encountered
% Variables naming mapping
% SPARSKITv2       This function
% ==============|==============
% a, ia, ja    --> val, row_tr, col_idx
% alu, jlu, ju --> lu_val, lu_index, lu_rowptr
% lfil, iwk    --> max_row_nnz, lu_buff_size
% jw(1 : n)    --> w_idx
% jw(n+1 : 2n) --> nz_idx
% ju0          --> lu_val_ptr
% tnorm        --> row_1norm

	% Allocate space for output
	lu_val   = zeros(lu_buff_size, 1);
	lu_index = zeros(lu_buff_size, 1);
	lu_uptr  = zeros(n, 1);
	
	if (max_row_nnz < 0)
		ierr = -4;
		return;
	end
	
	% Working arrays
	w      = zeros(n+1, 1);  % Working array, 1:ii-1 = L-part, ii:n = U-part
	w_idx  = zeros(n,   1);  % Indices for working array w
	nz_idx = zeros(n,   1);  % Nonzero indicators
	
	lu_val_ptr  = n + 2;    % Points to next element to be added to lu_val, lu_index
	lu_index(1) = lu_val_ptr;
	
	nz_idx(1 : n) = 0;     % Explicitly initial nnz indicator array
	
	% Beginning of the main loop
	for ii = 1 : n
		j1 = row_ptr(ii);
		j2 = row_ptr(ii + 1) - 1;
		
		% Computing 1-norm of current row for truncation
		row_1norm = 0;
		for k = j1 : j2
			row_1norm = row_1norm + abs(val(k));
		end
		if (row_1norm == 0)
			ierr = -5;
			return;
		end
		row_1norm = row_1norm / (j2 - j1 + 1);
		
		% Unpack L-part and U-part of row in A in array w
		% Nonzeros in L-part and U-part are stored contiguously in w
		len_l = 0; 
		len_u = 1;
		w(ii) = 0;
		w_idx(ii)  = ii;
		nz_idx(ii) = ii;
		for j = j1 : j2
			k = col_idx(j);
			t = val(j);
			if (k < ii)
				len_l = len_l + 1;
				w_idx(len_l) = k;
				w(len_l)     = t;
				nz_idx(k)    = len_l;
			end
			if (k == ii)
				w(ii) = t;
			end
			if (k > ii)
				len_u = len_u + 1;
				jpos  = ii + len_u - 1;
				w_idx(jpos) = k;
				w(jpos)     = t;
				nz_idx(k)   = jpos;
			end
		end
		
		% Eliminate previous rows
		len = 0;
		for jj = 1 : n       % Label 150
			if (jj > len_l)  % len_l may be changed during the loop
				break;       % Goto label 160
			end    
			
			% To do the elimination in the correct order we must select
			% the smallest column index among w_idx(jj+1 : len_l)
			jrow = w_idx(jj);
			k = jj;
			for j = jj + 1 : len_l
				if (w_idx(j) < jrow)
					jrow = w_idx(j);
					k = j;
				end
			end
			if (k ~= jj) 
				% Exchange w_idx()
				j         = w_idx(jj);
				w_idx(jj) = w_idx(k);
				w_idx(k)  = j;
				% Exchange nz_idx()
				nz_idx(jrow) = jj;
				nz_idx(j)    = k;
				% Exchange in w
				s     = w(jj);
				w(jj) = w(k);
				w(k)  = s;
			end
			
			nz_idx(jrow) = 0; % Zero out element in row by setting nz_idx(jrow) = 0 ?
			
			% Get the multiplier for row to be eliminated
			factor = w(jj) * lu_val(jrow);
			if (abs(factor) <= drop_tol)    % First truncation, if cannot pass then check next row
				continue;                   % Goto label 150
			end
			
			% Combining current row and row jrow
			for k = lu_uptr(jrow) : lu_index(jrow + 1) - 1 
				j = lu_index(k);                % Update position
				s = factor * lu_val(k);
				jpos = nz_idx(j);
				if (j >= ii)       % Dealing with upper part
					if (jpos == 0) % This is a fill-in element
						len_u = len_u + 1;
						if (len_u > n)
							ierr = -1;
							return;
						end
						i = ii + len_u - 1;
						w_idx(i)  = j;
						nz_idx(j) = i;
						w(i)      = -s;
					else           % Not a fill-in
						w(jpos) = w(jpos) - s;
					end
				else               % Dealing with lower part
					if (jpos == 0) % This is a fill-in element
						len_l = len_l + 1;
						if (len_l > n)
							ierr = -1;
							return;
						end
						w_idx(len_l) = j;
						nz_idx(j)    = len_l;
						w(len_l)     = -s;
					else           % Not a fill-in
						w(jpos) = w(jpos) - s;
					end
				end
			end
			
			% Store this pivot element (from left to right, no danger of overlapping
			% with the working elements in L)
			len = len + 1;
			w(len) = factor;
			w_idx(len) = jrow;
		end                  % Label 160
		
		% Reset U-part double-pointer to 0
		for k = 1 : len_u
			nz_idx(w_idx(ii + k - 1)) = 0;
		end
		
		% Update L-matrix
		len_l = len;
		len = min(len_l, max_row_nnz);
		
		% "Sort" by quick split, second truncation
		[w(1 : len_l), w_idx(1 : len_l)] = QuickSplit(w(1 : len_l), w_idx(1 : len_l), len_l, len);
		
		% Copy L-part to output buffer
		if (len + lu_val_ptr > lu_buff_size)
			ierr = -2;
			return;
		end
		for k = 1 : len
			lu_val(lu_val_ptr)    = w(k);
			lu_index(lu_val_ptr) = w_idx(k);
			lu_val_ptr = lu_val_ptr + 1;
		end
		
		% Save pointer to beginning of row ii of U-part
		lu_uptr(ii) = lu_val_ptr;
		
		% Update U-part (applying truncations)
		len = 0;
		% First truncation, based on dropping threshold
		for k = 1 : len_u - 1
			if (abs(w(ii + k)) > drop_tol * row_1norm)
				len = len + 1;
				w(ii + len)     = w(ii + k);
				w_idx(ii + len) = w_idx(ii + k);
			end
		end
		% "Sort" by quick split, second truncation
		len_u = len + 1;
		len = min(len_u, max_row_nnz);
		[w(ii + 1 : end), w_idx(ii + 1 : end)] = QuickSplit(w(ii + 1 : end), w_idx(ii + 1 : end), len_u - 1, len);
		
		% Copy U-part to output buffer
		if (len + lu_val_ptr > lu_buff_size)
			ierr = -3;
			return;
		end
		% t = abs(w(ii));                     % Seems to be uesless
		for k = ii + 1 : ii + len - 1
			lu_val(lu_val_ptr)    = w(k);
			lu_index(lu_val_ptr) = w_idx(k);
			% t = t + abs(w(k));              % Seems to be uesless
			lu_val_ptr = lu_val_ptr + 1;
		end
		
		% Store inverse of diagonal element of U
		if (w(ii) == 0)
			w(ii) = (0.0001 + drop_tol) * row_1norm;
		end
		lu_val(ii) = 1.0 / w(ii);
		
		% Update pointer to beginning of next row of L
		lu_index(ii + 1) = lu_val_ptr;
	end
	lu_val    = lu_val(1 : lu_val_ptr - 1);
	lu_index = lu_index(1 : lu_val_ptr - 1);
	ierr = 0;
end