function [ABD, Force] = ABDmatrix_Forcematrix(thetas, H, layer_thickness, E1, E2, v12, G12, a1, a2, dT)
    %ABDmatrix_Forcematrix : Calculates ABD & Force Matrices for all values of theta
    
    N = length(thetas);  % Unit: degrees
    Nf = [5000 5000 0]'; % Unit: N/m
    M = [-3000 0 0]';    % Unit: N
    Nt = zeros(3,1);
    Mt = zeros(3,1);
    A = zeros(3); B = zeros(3); D = zeros(3);
    for k=1:N
        [Qbk,~] = Qbar_Sbar(E1,E2,v12,G12,thetas(k));
        zk = -H + k*layer_thickness;
        zk1 = -H + (k-1)*layer_thickness;
        Ak = Qbk*(zk-zk1);         % Unit: Pa.m
        % Bk = Qbk*(zk^2-zk1^2)/2;   % Unit: Pa.m^2
        Dk = Qbk*(zk^3-zk1^3)/3;   % Unit: Pa.m^3
        A = A+Ak;
        % B = B+Bk;
        D = D+Dk;
        ABD = [A,B;B,D];
        [Epsilon_Th] = Epsilon_Thermal(thetas(k),dT,a1,a2);
        Nt = Nt + Qbk*Epsilon_Th*(zk-zk1); % N/m
        Mt = Mt + 0.5*Qbk*Epsilon_Th*(zk^2-zk1^2); % N
        Force = [Nf+Nt; M+Mt];
    end
end
