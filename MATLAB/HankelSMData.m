function [H,H0,H1]=HankelSMData(x,n,L)
[M,N]=size(x);
H=[];
H0=H;
H1=H;
T=floor(N/L);
for j=1:L
x0=x(:,(1+(j-1)*T):j*T);
Ht=[];
for k=1:(T-n+1)
	Y=x0(:,k:(k+n-1));
	Ht=[Ht Y(:)];
end
H=[H Ht];
H0=[H0 Ht(:,1:(T-n))];
H1=[H1 Ht(:,2:(T-n+1))];
end
end
