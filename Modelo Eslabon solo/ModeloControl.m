%% MODELO CONTROL
% Constante proporcional motor1
P1 = tf(50);
% Constante proporcional motor2
P2 = tf(50);

%% MODELO MOTORES

% Modelo motor 1 - REF: EC042B-3
% VA = 24V
Va1 = 24; % V
Ra1 = 1.09; % Ohms
La1 = 1.7E-03; % H
Kt1 = 0.0539; % Constante de torque (N m)/A
b1 = 6.7E-06; % (N m) / (rad/s)
J1 =  2.1E-05; % kg*m^2
Kb1 = 0.0539; % V/(rad/s)
TL1 = 0;

% Modelo motor 1 - REF: EC042B-3
% VA = 24V
Va2 = Va1; % V
Ra2 = Ra1; % Ohms
La2 = La1; % mH
Kt2 = Kt1; % Constante de torque (N m)/A
b2 = b1; % (N m) / (rad/s)
J2 = J1; % kg*m^2
Kb2 = Kb1; % V/(rad/s)
TL2 = 0;

input_a1 = timeseries(angulo_e1+pi,tiempo);
input_a2 = timeseries(angulo_e2+pi,tiempo);

%% Caracterización 

% data1 = load('CaracterizacionPlanta1.mat').data ;
planta1 = extractTimetable(load('CaracterizacionPlanta1.mat').data );
planta2 = extractTimetable(load('CaracterizacionPlanta2.mat').data );

[p1, tau_m1, km1] = identificarPlanta(seconds(planta1.Time), planta1.theta1, 20);
[p2, tau_m2, km2] = identificarPlanta(seconds(planta2.Time), planta2.theta2, 20);

s = tf('s');

sys1 = km1/(s*(tau_m1*s +1));
sys2 = km2/(s*(tau_m2*s +1));

%% Functions
function [p, tau_m, Km] = identificarPlanta(x, y, A)
    figure

    % zona de rampa para tiempo de caract 2.5
    init = 0.12; % 0.3 s inicio
    endr = 0.80; % 2 s fin
    samples = size(x, 1);
    x_ramp = x(round(init*samples):round(endr*samples));
    y_ramp = y(round(init*samples):round(endr*samples));
    p = polyfit(x_ramp,y_ramp,1); 
    Km = p(1)/A;
    tau_m = roots(p);
    f = polyval(p,x_ramp); 
    plot(x,y,x_ramp,f,'--') 
    title('Respuesta al paso')
    ylabel('Posición (rad)')
    xlabel('Tiempo (s)')
    grid on
end