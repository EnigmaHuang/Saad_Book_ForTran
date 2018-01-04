function [val, idx] = QuickSplit(val0, idx0, n, n_cut)
% Do a quick-sort split for a real array s.t. |val(i)| >= |val(n_cut)|
% and |val(j)| <= |val(n_cut)| for any 1 <= i <= n_cut and n_cut <= j.
% val   : Array to be split
% idx   : Array to be permuted in the same way as val
% n     : Length of val
% n_cut : Number of elements we need to split from val
	val = val0;
	idx = idx0;
	
	first = 1;
	last  = n;
	if ((n_cut < first) || (n_cut > last))
		return;
	end
	
	while (1) % Outer loop
		mid = first;
		abskey = abs(val(mid));
		for j = first + 1 : last
			if (abs(val(j)) > abskey)
				mid = mid + 1;
				% Swap j-th and mid-th elements
				val_tmp  = val(mid);
				val(mid) = val(j);
				val(j)   = val_tmp;
				idx_tmp  = idx(mid);
				idx(mid) = idx(j);
				idx(j)   = idx_tmp;
			end
		end
			
		% Swap first and mid-th elements
		val_tmp    = val(mid);
		val(mid)   = val(first);
		val(first) = val_tmp;
		idx_tmp    = idx(mid);
		idx(mid)   = idx(first);
		idx(first) = idx_tmp;
		
		% Test for while loop
		if (mid == n_cut)
			break;
		end
		if (mid > n_cut)
			last  = mid - 1;
		else
			first = mid + 1;
		end
	end
end