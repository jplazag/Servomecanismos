%% MODELO CONTROL
% Constante proporcional motor1
P1 = 1;
% Constante proporcional motor2
P2 = 1;

%% MODELO MOTORES

% Modelo motor DC054B–2%
% VA = 24V
Va = 24; % V
Ra = 1.73; % Ohms
La = 2.5; % mH
Kt = 0.0551; % Constante de torque (N m)/A
b = 1.1E-05; % (N m) / (rad/s)
J = 1.6E-05; % kg*m^2
Kb = 0.0551; % V/(rad/s)
TL = 0;


%% Caracterización
A = [2.03034366646370,338.129178581801];
B = [1.20000573358804,170.780488526573];

x = linspace(0,10);
m = (B(2)-B(1))/(A(2)-A(1));
n = B(2)*m - A(2);
y = m*x + n;

line(x,y)
hold on