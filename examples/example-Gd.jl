#
println("Gd) Symbolic simplification of Racah expressions including spherical harmonics and Wigner rotation matrices.")
#
using SymEngine
j1  = Basic(:j1);    j2  = Basic(:j2);    j3  = Basic(:j3);    j4  = Basic(:j4);    j5  = Basic(:j5);   j  = Basic(:j)
l1  = Basic(:l1);    l2  = Basic(:l2);    l3  = Basic(:l3);    l4  = Basic(:l4);    l5  = Basic(:l5)
m1  = Basic(:m1);    m2  = Basic(:m2);    m3  = Basic(:m3);    m4  = Basic(:m4);    m5  = Basic(:m5)
p1  = Basic(:p1);    p2  = Basic(:p2);    p3  = Basic(:p3);    p4  = Basic(:p4);    p5  = Basic(:p5);   p  = Basic(:p)
q1  = Basic(:q1);    q2  = Basic(:q2);    q3  = Basic(:q3);    q4  = Basic(:q4);    q5  = Basic(:q5);   q  = Basic(:q)

alpha  = Basic(:alpha);    beta  = Basic(:beta);    gamma  = Basic(:gamma)
alpha1 = Basic(:alpha1);   beta1 = Basic(:beta1);   gamma1 = Basic(:gamma1)
theta  = Basic(:theta);    phi   = Basic(:phi)
theta1 = Basic(:theta1);   phi1  = Basic(:phi1)
theta2 = Basic(:theta2);   phi2  = Basic(:phi2)
J   = Basic(:J)

if  true
    # Last successful:  unknown ...
    # Compute 
    w1 = RacahAlgebra.DFunction(j1, p1, q1, alpha, beta, gamma)
    w2 = RacahAlgebra.Djpq(j,p,q,beta)
    wa = RacahAlgebra.Ylm(l1,m1,theta1,phi1)
    RacahAlgebra.symmetricForms(w2)
    #
elseif false
    # Last successful:  unknown ...
    # Compute 
    #
end
