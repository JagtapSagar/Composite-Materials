function [ Epsilon_Th ] = Epsilon_Thermal( t, dT, a1, a2 )
    %Epsilon_Thermal : Calculates Thermal coefficients
    m = cosd(t); m2 = m.^2;
    n = sind(t); n2 = n.^2;
    Epsilon_Th = dT*[m2*a1+n2*a2; n2*a1+m2*a2; 2*m*n*(a1-a2)];
end