function X=SINDy(A,B,N,lambda)
% Sparse linear least squares system solver based on 
% program sparsifyDynamics.m
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
	NA=diag(sqrt(diag(A'*A)));
	A=A/NA;
	X = NA\(A\B);  % initial guess: Least-squares
	m=size(B,2);
	n=size(A,2);

% lambda is our sparsification knob.
K=sqrt(m*n)*norm(A,'fro');
X0=X;
for k=1:N
    smallinds = (abs(X)<lambda);   % find small coefficients
    X(smallinds)=0;                % and threshold
    for ind = 1:m                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        X(biginds,ind) = A(:,biginds)\B(:,ind); 
    end
    X=NA\X;
    Error=norm(X-X0,'inf');
    X0=X;
    if Error<=lambda
	break;
	return;
    end
end
end
