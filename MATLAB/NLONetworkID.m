function [c,c0,t0,x0]=NLONetworkID(t,data,L,sp,tol,thr)
	% Example:
	% [t,x]=DuffingNetwork(1,-36,0,.2,[8 7 4 15 14 9],[0 20]);
	% [st,sx]=DataSpliner(t,x);
	% [c,c0,t0,x0]=NLONetworkID(st,sx,9,.15,1e-5,1e-5);
	% [rt0,rx0]=NLONetwork(c0,x0,[t0 20],9);
	% [rt,rx]=NLONetwork(c,x0,[t0 20],9);
	% [t,x]=DuffingNetwork(1,-36,0,.2,x0,[t0 20]);
	% E=eye(3);            
	% g1=kron(eye(2),E([2 3 1],:));
	% g2=kron(eye(2),E([1 3 2],:));
	% subplot(311),plot(t,x(:,1:3).'),axis tight;
	% subplot(312),plot(rt,rx(:,1:3).'),axis tight;
	% subplot(313),plot(rt0,rx0(:,1:3).'),axis tight;
	% Error=sqrt(sum((x-rx).^2,2));
	% Error0=sqrt(sum((x-rx0).^2,2));
	% SError=sqrt(sum((NLONetworkRhs(c,x*g1.',9)-g1*NLONetworkRhs(c,x,9)).^2));
	% SError0=sqrt(sum((NLONetworkRhs(c0,x*g1.',9)-g1*NLONetworkRhs(c0,x,9)).^2));
	% SError1=sqrt(sum((NLONetworkRhs(c,x*g2.',9)-g2*NLONetworkRhs(c,x,9)).^2));
	% SError01=sqrt(sum((NLONetworkRhs(c0,x*g2.',9)-g2*NLONetworkRhs(c0,x,9)).^2));
	% Lt=1:length(rt);
	% Lt1=1:length(t);
	% figure,
	% subplot(311),semilogy(Lt,Error,'k',Lt,Error0,'r-.'),axis tight,grid on;
	% subplot(312),semilogy(Lt1,SError,'k',Lt1,SError0,'r-.'),axis tight,grid on;
	% subplot(313),semilogy(Lt1,SError1,'k',Lt1,SError01,'r-.'),axis tight,grid on;
	%
	M=floor(size(data,2)/2);	
	sp=min(sp,1);
	SData=data;
	for k=1:L
		SData=[SData data(:,1:M).^(k+1)];
	end
	N=floor(sp*size(data,1));
	DS=SData(1:N,:);
	dt=t(2)-t(1);
	rhs=DtData(data.',dt,4).';
	rhs=rhs(1:N,:);
	disp('Running SDSI:');
	tic,c=SpSolver(DS,rhs,M*(L+2),tol,thr).';toc
	disp('Running SINDy:');
	tic,c0=SINDy(DS,rhs,500,tol).';toc
	x0=data(N,:);
	t0=t(N);
end
function y=Fy(x,N,L)
	y=x.';
	for k=1:L
		y=[y (x(1:N).').^(k+1)];
	end
end
