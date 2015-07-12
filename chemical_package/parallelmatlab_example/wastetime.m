function [r]=wastetime
    x=1;
    for i=1:300000000
    x=x*99999999^-99999999;
    x=x*99999999^99999999;
    end
end