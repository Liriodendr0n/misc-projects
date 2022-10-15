function out = DSvar(Nvar, Mder, varnum, val)
%DSvar is a helper function that removes the redundant DS argument
    out = makeDSvar(DS, Nvar, Mder, varnum, val);
end

