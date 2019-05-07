function [ Qbar,Sbar ] = Qbar_Sbar( E1,E2,v12,G12,theta )
    % QBAR_SBAR : Calculates Qbar & Sbar matrices
    S11 = 1/E1;
    S22 = 1/E2;
    S66 = 1/G12;
    S12 = -v12/E1;
    S21=S12;
    S = [S11 S12 0; S21 S22 0; 0 0 S66];
    Q = inv(S);
    [T1,T2] = T1_T2(theta);
    Sbar = inv(T2)*S*T1;
    Qbar = inv(T1)*Q*T2;
end
