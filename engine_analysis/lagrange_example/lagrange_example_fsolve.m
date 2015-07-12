function lagrange_example_fsolve
    %my guesses
    var0(1) = 1;
    var0(2) = 1;
    var0(3) = 1;
    var = fsolve(@fun,var0);
    var
    function F=fun(var)
       x = var(1);
       y = var(2);
       lam = var(3);
       F1 = 1-lam*2*x;
       F2 = 1-lam*2*y;
       F3 = x^2 + y^2 -1;
       F = [F1 ; F2 ; F3];
    end
end



