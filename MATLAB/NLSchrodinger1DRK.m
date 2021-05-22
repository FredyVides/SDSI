function [t,x,w,wdata]=NLSchrodinger1DRK(A,c,p,V,L,N,m,t,ss,sn)
%%
% Example:
% [t,x,w,wdata]=NLSchrodinger1DRK(1,1,1,1,20,160,800,[0 8],400,1);
% [t,x,w,wdata1]=NLSchrodinger1DRK(wdata(:,end),1,1,1,20,160,800,[0 8],400,2);
% wdata=[wdata(:,1:(end-1)) wdata1];
% x=x(1,:)';
% Wdata=[x [zeros(1,801);wdata;zeros(1,801)]];
% csvwrite('NLSEqData.csv',Wdata)
%%
hx=2*L/N;
ht=diff(t)/m;
L1D=c*spdiags(ones(N-1,1)*[1 -2 1],-1:1,N-1,N-1)/hx^2;
x=-L:hx:L;
if sn==1
	w0L=A/sqrt(p)*sech(-A/sqrt(2)*L).*exp(-i/2*V*L);
        a=-i/L*sech(L/sqrt(2))*sin(L/2);
        b=-sech(L/sqrt(2))*cos(L/2);
        w0=A/sqrt(p)*sech(A/sqrt(2)*x).*exp(i/2*V*x)+a*x+b;
	w=w0';
	w0=w(2:N);
else
	w0=A;
	A=max(abs(A(1:(N-1))));
	w=[0;w0(1:(N-1),1);0];
end
wdata=w0;
ss=floor(m/ss);
for k=1:m
subplot(211),plot(x,[0;real(w0(1:(N-1),1));0]);
axis([-L,L,-1.5*A,1.5*A]);
subplot(212),plot(x,[0;imag(w0(1:(N-1),1));0]);
axis([-L,L,-1.5*A,1.5*A]);
pause(.2);
W1=SLRHS(L1D,w0,p);
W2=SLRHS(L1D,w0+W1*ht/2,p);
W3=SLRHS(L1D,w0+W2*ht/2,p);
W4=SLRHS(L1D,w0+W3*ht,p);
w0=w0+ht*(W1+2*(W2+W3)+W4)/6;
if mod(k,ss)==0
w=[w,[0;w0(1:(N-1),1);0]];
wdata=[wdata,w0];
end
end
Nw=size(w,2);
[x,t]=meshgrid(x,t(1):diff(t)/(Nw-1):t(2));
waterfall(x,t,abs(w').^2)
end
function y=SLRHS(L,x,p)
	y=i*L*x+i*p*(abs(x).^2).*x;
end
