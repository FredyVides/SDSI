#Created on Wed Mar 31 02:57:52 2021
#SPSOLVER  Sparse linear regression solver
#   Code by Fredy Vides
#   For Paper, "Sparse system identification by low-rank approximation"
#   by F. Vides
function spsolver(A,Y,L,tol,delta)
    M=size(A,2);
    N=size(Y,2);
    X=zeros(M,N);
    u,s,v=svd(A);
    rk=sum(s.>tol);
    u=u[:,1:rk];
    s=s[1:rk];
    s=s.^(-1);
    s=diagm(s);
    v=v[:,1:rk];
    A=u'*A;
    Y=u'*Y;
    X0=v*(s*Y);
    w=zeros(M,1);
    for k=1:N
        K=1;
	    c=X0[:,k];
	    x0=c;
	    ac=abs.(c);
	    f=sortperm(-ac[:,1]);
	    N0=max(sum(ac[f].>delta),1);
        Error = 1+tol;
        while (K<=L) && (Error>tol)
            ff=f[1:N0,1];
            X[:,k]=w;
            c=A[:,ff]\Y[:,k];
            X[ff,k]=c;
            Error=norm(X[:,k]-x0,Inf);
            x0=X[:,k];
            ac=abs.(x0);
            f=sortperm(-ac[:,1]);
            N0=max(sum(ac[f].>delta),1);
            K=K+1;
        end
    end
    return X
end
