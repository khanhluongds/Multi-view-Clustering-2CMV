function L = constructL(W, alpha, options)
	nSmp = size(W,1);
	if alpha > 0
		%the degree is sum of rows for W, and create diagonal matrix with degree values
		W = alpha*W;
		DCol = full(sum(W,2));
		D = spdiags(DCol,0,nSmp,nSmp);
		L = D - W;
		%normalize L, 
		%L = D^-1*(D-W);
		if isfield(options,'NormW') && options.NormW ==1
			D_mhalf = spdiags(DCol.^-.5,0,nSmp,nSmp) ;
			L = D_mhalf*L*D_mhalf;
		end
		if isfield(options,'NormW') && options.NormW == 2
			L = D^-1*(D-W);
		end
		if isfield(options,'NormW') && options.NormW == 3
			CKSym = abs(W+W')/2;
			N = size(CKSym,1);

			DKS = ( diag( sum(CKSym) ) )^(-1/2);
			L = speye(N) - DKS * CKSym * DKS;
		end
	else
		L = [];
	end