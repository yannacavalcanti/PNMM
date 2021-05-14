function A = generate_abundances(A)

A(1,A(1,:)<0.1)=0;