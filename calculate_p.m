function p = calculate_p(A,B,C,state_dim)
Observability=C;
Controllability=B;
p=1;
while rank(Observability)<state_dim || rank(Controllability)<state_dim
    Observability=[C*(A^p);Observability];
    Controllability=[Controllability,(A^p)*B];
    p=p+1;
end
end

