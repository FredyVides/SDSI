function [q,r]=rk(A,delta)
	[q,r,~]=svd(A,0);
	r=diag(r);
	r=sum(r>delta);
	q=q(:,1:r);
end
