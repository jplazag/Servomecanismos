ent_car=20;

%Caracterizacion en tiempo real
%car1=sim('PlantaAmpliada1.slx',0.2);
%car2=sim('PlantaAmpliada2.slx',0.2);

%Caracterizacion del informe
planta1 = extractTimetable(load('CaracterizacionPlanta1.mat').data );
planta2 = extractTimetable(load('CaracterizacionPlanta2.mat').data );

%% Caracterización por posición

%Caracterizacion del informe
[p1, tau_m1, km1] = identificarPlantaPosicion(seconds(planta1.Time), planta1.theta1, ent_car);
[p2, tau_m2, km2] = identificarPlantaPosicion(seconds(planta2.Time), planta2.theta2, ent_car);

%Caracterizacion en tiempo real
%[p1, tau_m1, km1] = identificarPlantaPosicion(car1.tout, car1.simth1.data, ent_car);
%[p2, tau_m2, km2] = identificarPlantaPosicion(car2.tout, car2.simth2.data, ent_car);

s = tf('s');

sys1 = km1/(s*(tau_m1*s +1));
sys2 = km2/(s*(tau_m2*s +1));

%% Caracterización por velocidad

s = tf('s');

%Caracterizacion del informe
[tau_m1, km1] = identificarPlantaVelocidad(seconds(planta1.Time), planta1.Wm1, ent_car);
[tau_m2, km2] = identificarPlantaVelocidad(seconds(planta2.Time), planta2.Wm2, ent_car);

%Caracterizacion en tiempo real
%[tau_m1, km1] = identificarPlantaVelocidad(car1.tout, car1.simwm1.data, ent_car);
%[tau_m2, km2] = identificarPlantaVelocidad(car2.tout, car2.simwm2.data, ent_car);

tf_sys1 = km1/(tau_m1*s +1);
tf_sys2 = km2/(tau_m2*s +1);

tf_pos1 = tf_sys1 * (1/s);
tf_pos2 = tf_sys2 * (1/s);
zeta = 1; % Criticamente amortiguado
% Constante proporcional motor1
P1 = 1 / (4*km1*tau_m1*zeta^2);
% Constante proporcional motor2
P2 = 1 / (4*km2*tau_m2*zeta^2);

%% Variables de Simulacion
Den_Acrilico = 1.18 / 1000; % kg/cm^3

run('Robot_eslabones_DataFile.m');

%% MODELO MOTORES PARA LA MODELACION

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
N_catalogo_1 = 200;

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
N_catalogo_2 =  16;

input_a1 = timeseries(angulo_e1+pi/2,tiempo);
input_a2 = timeseries(angulo_e2-pi,tiempo);

%% Functions
function [tau, Km] =  identificarPlantaVelocidad(x,y, Vin)
    graficar = false;
    
    maxVel = max(y);
    bin = y >= 0.9502129 * maxVel;
    pos_tau = find(bin, 1, "first");
    tau = x(pos_tau - 1)/3;
    Km = y(end)/Vin;
    if graficar
        figure
        plot(x,y);
        title('Respuesta al paso')
        ylabel('Velocidad (rad/s)')
        xlabel('Tiempo (s)')
    
        dim = [0.6 0.2 0.3 0.3];
        str = { strcat('W_{max} = ', num2str(maxVel)), strcat('\tau_m = ', num2str(tau)), strcat('K_m = ',num2str(Km))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    %     dt = datatip(s,pos_tau, y(pos_tau));
    %     s.DataTipTemplate.DataTipRows(1).Label = 'Drivers (in thousands)';
    %     s.DataTipTemplate.DataTipRows(2).Label = 'Fatalities'; 
        grid on
    end
end

function [p, tau_m, Km] = identificarPlantaPosicion(x, y, A)
    graficar = false;

    % zona de rampa para tiempo de caract 2.5
    init = 0.03; % 0.1 s inicio
    endr = 0.80; % 2 s fin
    samples = size(x, 1);
    x_ramp = x(round(init*samples):round(endr*samples));
    y_ramp = y(round(init*samples):round(endr*samples));
    p = polyfit(x_ramp,y_ramp,1); 
    Km = p(1)/A;
    tau_m = roots(p);
    if graficar
        figure
        f = polyval(p,x_ramp); 
        plot(x,y,x_ramp,f,'--') 
        title('Respuesta al paso')
        ylabel('Posición (rad)')
        xlabel('Tiempo (s)')
        grid on
    end
end
