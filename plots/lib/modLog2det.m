% Purpose: get the log2 (Square of det(Matrix))
function logdet = log2det(A)
ei = diag(chol(A));
%disp([ei]);
en = length(ei);
expo = zeros(1,en);
for k=1:en
    while(ei(k)<1)
        ei(k) = ei(k)*2;
        expo(k) = expo(k)+1;                                               %Numerical Method
    end
end

logdet = 2*(log2(prod(ei))-sum(expo));

%disp ([logdet]);

end