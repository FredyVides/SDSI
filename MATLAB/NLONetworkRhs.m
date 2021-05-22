function y=NLONetworkRhs(c,x,L)
	N=floor(size(x,2)/2);
	y=x.';
	for k=1:L
		y=[y;(x(:,1:N).').^(k+1)];
	end
	y=c*y;
end
