function A = MSR2COO(n, val, ind)
% Modified Sparse Row format to COO format,
% for showing the results of ILUT().
% n   : Row dimension of A (A should be square)
% val : Matrix stored in Modified Sparse Row (MSR), (1 : n) is INVERSED
%       diagonal elements, (n + 2 : end) are off-diagonal elements
% ind : (1 : n) is the beginning position of each row in val,
%       (n + 2 : end) is the column index for each off-diagonal elements
	nnz  = max(size(val)) - 1;
	vals = zeros(nnz, 1);
	rows = zeros(nnz, 1);
	cols = zeros(nnz, 1);
	cnt  = 0;
	
	% Copy the diagonal elements
	for i = 1 : n
		cnt = cnt + 1;
		rows(cnt) = i;
		cols(cnt) = i;
		vals(cnt) = 1.0 / val(i);
	end
	
	% Copy the values of each row
	for i = 1 : n
		for jcol = ind(i) : ind(i + 1) - 1
			cnt = cnt + 1;
			rows(cnt) = i;
			cols(cnt) = ind(jcol);
			vals(cnt) = val(jcol);
		end
	end
	
	A = sparse(rows, cols, vals, n, n, nnz);
end