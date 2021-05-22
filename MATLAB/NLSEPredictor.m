function [w,wdata]=NLSEPredictor(ht,hp,w0,x,L,p,Nt)
%%
% Example:
% Import and preprocess data (data generated with NLSchrodigner1DRK.m):
% wdata=csvread('NLSEqData.csv');
% tol=1e-4;Threshold= 6e-1;noise=1e-6;
% [M,N]=size(wdata);
% x=wdata(:,1)';
% Wdata=wdata(:,2:N)+noise*randn(161,4401);
% wdata=Wdata(2:(M-1),:);
% dt=.01;
% Compute approximate transition operators
% L=200;S=15;
% tic,[H,h0,c,c0]=NLSESpModelID(2*i*dt,wdata,4,4,S,tol,L,Threshold);toc
% subplot(211),spy(H),subplot(212),stem(abs(c))
% Visualize predictions
% [w,wdata1]=NLSEPredictor(dt,h0,x,H,c(4:203),800);
% 
%
N=length(w0);
A=w0;
A=max(abs(A(1:(N-1))));
w=[0;w0;0];
wdata=w0;
j=1;
mx=max(x);
for k=1:Nt
if mod(k-1,hp)==0    
    w=[w,[0;w0;0]];
    wdata=[wdata,w0];
    plot(x,[0;abs(w0);0],'r');
    axis([-mx,mx,-2.5*A,2.5*A]);
    pause(.2);
    j=j+1;
end
W1=SLRHS(L,w0,p);
W2=SLRHS(L,w0+W1*ht/2,p);
W3=SLRHS(L,w0+W2*ht/2,p);
W4=SLRHS(L,w0+W3*ht,p);
w0=w0+ht*(W1+2*(W2+W3)+W4)/6;
end
end
function y=SLRHS(L,x,p)
	y=i*L*x;
	for k=1:200
	y=y+i*p(k)*(abs(x).^k).*x;
	end
end
