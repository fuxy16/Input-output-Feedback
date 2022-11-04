function M = calculate_M(A,B,C,p)
O=C;
CC=B;
[m,n]=size(C*B);
T=zeros(m,n);
for i=1:(p-1)
    O=[C*(A^i);O];
    CC=[CC,(A^i)*B];
    T=[T,C*A^(i-1)*B];
end
for i=0:(p-2)
    T=[T;zeros(m,n),T((i*m+1):(i*m+m),1:(p*n-n))];
end

O_inv=inv(O'*O)*O';
M=[CC-(A^p)*O_inv*T,(A^p)*O_inv];
end

