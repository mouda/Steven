function A = covm(x,y,sigma,a)

size = max([length(x), length(y)]);
A = zeros(size,size);

for k = 1:size
    for l = 1:size
        A(k,l) = sigma^2*exp(-a*((x(k)-x(l))^2+(y(k)-y(l))^2));
    end
end

end
