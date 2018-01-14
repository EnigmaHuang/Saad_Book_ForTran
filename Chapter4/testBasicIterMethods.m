function testBasicIterMethods(p)
% Test Jacobi, Gauss-Seidel, blocked Jacobi and blocked Gauss-Seidel methods
% Test matrix: 2D Possion equation with FDM on (p-2) * (p-2) square grid, random RHS
	A = delsq(numgrid('S', p));
	n = size(A, 1);
	b = rand(n, 1);
	eps = rand(n) * 0.001;
	A = A + eps + eps';
	
	% Orginal Jacobi and Gauss-Seidel
	[~, ~, ic_ja, rn_ja] = Jacobi_Iter(A, b);
	[~, ~, ic_gs, rn_gs] = GS_Iter(A, b);
	
	% Create the block partitioning for blocked Jacobi and blocked Gauss-Seidel
	block_spos = 0 : p-2 : n;
	block_spos = block_spos + 1;
	
	[~, ~, ic_bja, rn_bja] = Block_Jacobi_Iter(A, b, block_spos);
	[~, ~, ic_bgs, rn_bgs] = Block_GS_Iter(A, b, block_spos);
	
	% Plot results
	rn_ja  = rn_ja  ./ rn_ja(1);
	rn_gs  = rn_gs  ./ rn_gs(1);
	rn_bja = rn_bja ./ rn_bja(1);
	rn_bgs = rn_bgs ./ rn_bgs(1);
	semilogy(1:ic_ja, rn_ja, 'r-', 1:ic_gs, rn_gs, 'b-', 1:ic_bja, rn_bja, 'g-', 1:ic_bgs, rn_bgs, 'c-'), hold on
	grid on, xlabel('Iterations'), ylabel('Relative Residual Norm'), hold on
	legend('Jacobi', 'Gauss-Seidel', 'Block Jacobi', 'Block Gauss-Seidel'), hold on
	title_str1 = 'Basic Iteration Methods & Their Block Version for Solving Ax = b';
	title_str2 = ['Matrix: delsq(numgrid(''S'', ' int2str(p) ')) + 0.002 * rand(' int2str((p-2)*(p-2)) ')'];
	title({[title_str1]; [title_str2]}), hold off
end