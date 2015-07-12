function odeexampleplot
	x=lsode("f",2,(t=linspace(0,5,50)'));
	plot(t,x)


	function xdot = f (x,t)
		xdot=-exp(-t)*x^2;
	endfunction




endfunction

