function [t,sdata]=DataSpliner(t,data)
	[M,N]=size(data);
	x=cos(pi*((N-1):-1:0)/(N-1));
	[X,T]=meshgrid(x,t);
	t=t(1):(t(M)-t(1))/(M-1):t(M);
	[X,S]=meshgrid(x,t);
	sdata=interp2(X,T,data,X,S,'spline');
end
