function [J]=getfacobian()


% assemble JC
JC = [DWDC DWDT; ...
      JCTC JCTT];
end