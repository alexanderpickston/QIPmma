function opts = sdp_StateTomo(v,S)
% Solves the SDP to fit a postive semi-definite matrix to measured
% tomographic data.
%
% Arguments:
%	- S:    The (n x d) matrix corresponding to quadratic form
%   - v:    vector measurements
%
% Outputs:
%	- opts:     The (d x d) matrix that is the solution to the primal problem
%	


dims = size(S);
d = sqrt(dims(2));
di = log2(d);
I = eye(di);


cvx_precision high
cvx_begin sdp
% Tomo SDP
	variable X(d,d) hermitian
	minimize norm(S*reshape(X,[d*d,1])-v)
	subject to
      X >= 0
	  trace(X)-1 == 0
cvx_end

opts = X;


end