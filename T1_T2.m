function [ T1, T2 ] = T1_T2( t )
    %T1_T2 : finds T1 & T2
    m = cosd(t); m2 = m.^2;
    n = sind(t); n2 = n.^2;
    T1 = [m2 n2 2*m*n; n2 m2 -2*m*n; -m*n m*n m2-n2];
    T2 = [m2 n2 m*n; n2 m2 -m*n; -2*m*n 2*m*n m2-n2];
end