% Charakterystyka liniowa,
% obiekt w postaci ss(A,B,C,D)
% 1 wejscie
% 2 wyjscia

function plant = servo()
    Kt = 1280.2;  % Torsional rigidity
    Km = 10;      % Motor constant
    Jm = 0.5;     % Motor inertia
    Jl = 50*Jm;   % Load inertia
    N = 20;       % Gear ratio
    Bm = 0.1;     % Rotor viscous friction
    Bl = 25;      % Load viscous friction
    R = 20;       % Armature resistance

    A = [        0       1             0                0;
            -Kt/Jl  -Bl/Jl     Kt/(N*Jl)                0;
                 0       0             0                1;
         Kt/(Jm*N)       0  -Kt/(Jm*N^2)  -(Bm+Km^2/R)/Jm];
    B = [0; 0; 0; Km/(R*Jm)];
    C = [  1  0       0  0; 
          Kt  0   -Kt/N  0];
    D = [0; 0];
    
    plant = ss(A,B,C,D);
end