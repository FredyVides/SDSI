function [t,x]=NLONetwork(C,x0,T,L)
	% Example:
	%
	% C=spdiags(ones(4,1)*[.2 -36-.4 .2],-1:1,4,4);
	% C=[sparse(4,4) speye(4);C sparse(4,4)];      
	% C=[C [sparse(4,4);-speye(4)]];  
	% [t,x]=NLONetwork(C,[4 2 1 3 0 0 0 0],[0 5],1);
	%
	N=max(size(x0))/2;
	f=@(t,y)C*Fy(y,N,L).';
	opt=odeset ("RelTol", eps);
	[t,x] = ode45 (f, [T(1), T(2)], x0,opt);
end
function y=Fy(x,N,L)
	y=x.';
	for k=1:L
		y=[y (x(1:N).').^(k+1)];
	end
end
