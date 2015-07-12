function [t]=odeexample(g)
    %input is intial condition
	x=lsode("f",g,(t=linspace(0,5,50)'));
    t(:,2)=x;
    %output is first column time, 2nd column x

	function xdot = f (x,t)
		xdot=-exp(-t)*x^2;
	endfunction




endfunction

