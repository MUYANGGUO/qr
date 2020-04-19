% Return householder matrix 
function A = HH(v)

I = eye(length(v));
A = I - ( (2/(v' * v)) * (v * v') );

end