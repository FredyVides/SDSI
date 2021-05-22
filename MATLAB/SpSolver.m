% function x=SpSolver(A,y,L,tol,delta)
%SPSOLVER  Sparse linear regression solver
%   Code by Fredy Vides
%   For Paper, "Sparse system identification by low-rank approximation"
%   by F. Vides
function X=SpSolver(A,Y,L,tol,delta)
if nargin<=4 
	delta=tol;
end
N=size(Y,2);
if length(L)==1, L=L*ones(1,N);end
[u,s,v]=svd(A,0);
s=diag(s);
rk=sum(s>tol);
u=u(:,1:rk);
s=diag(1./s(1:rk));
v=v(:,1:rk);
A=u'*A;
Y=u'*Y;
X0=v*(s*Y);
M=size(A,2);
w=sparse(M,1);
for k=1:N
	K=1;
	c=X0(:,k);
	x0=c;
	ac=abs(c);
	[~,f]=sort(-ac);
	N0=max(sum(ac(f)>delta),1);
	Error = 1+tol;
	while K<=L & Error>tol
		ff=sort(f(1:N0));
		X(:,k)=w;
		c=A(:,ff)\Y(:,k);
		X(ff,k)=c;
		Error=norm(x0-X(:,k),'inf');
		x0=X(:,k);
		ac=abs(x0);
		[~,f]=sort(-ac);
		N0=max(sum(ac(f)>delta),1);
		K=K+1;
	end
end
end
