close all
clear

dt = 0.05;

r0 = 10;
ra = r0*3/10;
n = 3;
x0 = 25;
y0 = 15;
theta0_f = 45*pi/180;
theta_f = flip(0:dt:2*pi); %linspace(0,2*pi,nt);
r = r0-ra*sin(n*(theta_f + theta0_f));
scale = 1;

% R = sqrt(x.^2+y.^2);
% disp(max(R));
% [dmax,in]=max(R);

L1 =21.5; %1.0*ceil(max(R))/2;    cm
L2 =21.5; %1.0*ceil(max(R))/2;    cm


% Selección de la pendiente de la trayectoria recta de aproximación -------

x_inicio = min(1.3*0.8*(r0-ra*sin(n*(theta_f + theta0_f))) .* sin(theta_f)+ x0) - 5;

theta_pre_dev = (2*pi:-dt:0) - theta0_f;

x_pre_dev = scale*0.8*(r0-ra*sin(n*(theta_pre_dev + theta0_f))) .* sin(theta_pre_dev)+ x0 ;
y_pre_dev = scale*0.8*(r0-ra*sin(n*(theta_pre_dev + theta0_f))) .* cos(theta_pre_dev)+ y0 ;

derivadas_trebol = diff( [y_pre_dev y_pre_dev(1)] )./diff([x_pre_dev x_pre_dev(1)]);

devs_dist = abs( derivadas_trebol  - (y_pre_dev - 0)./(x_pre_dev - x_inicio) ) .* (sqrt((y_pre_dev - 0).^2 + (x_pre_dev - x_inicio).^2));

indice_r = find(devs_dist == min(devs_dist));

x_final_r = x_pre_dev(indice_r);
y_final_r = y_pre_dev(indice_r);


pendiente_r = (y_final_r - 0) / (x_final_r - x_inicio);


% Selección de velocidad de empalme y la aceleración necesaria ------------

v_r_x = diff(x_pre_dev)./abs(diff(theta_pre_dev));

v_r_x_f = v_r_x(indice_r);

a_r = v_r_x_f^2 / (2 * (x_final_r - x_inicio) );


% Construcción de la trayectoria recta de aproximación --------------------


t_r = 0:dt: sqrt(2*(x_final_r - x_inicio)/a_r);

if t_r(end) ~= sqrt(2*(x_final_r - x_inicio)/a_r)
    t_r = [t_r sqrt(2*(x_final_r - x_inicio)/a_r)];
end

x_r = a_r*t_r.^2/2 + x_inicio;
y_r = pendiente_r * (x_r - x_inicio);


% Trayectoria de interés --------------------------------------------------

theta_f = (2*pi:-dt:0) + theta_pre_dev(indice_r) - dt ;

x = scale*0.8*(r0-ra*sin(n*(theta_f + theta0_f))) .* sin(theta_f)+ x0;
y = scale*0.8*(r0-ra*sin(n*(theta_f + theta0_f))) .* cos(theta_f)+ y0;

% figure
% plot([x_inicio x_final_r],[0 y_final_r], 'r'); hold on;
% plot(x,y);

% Preparación de variables para almacenar los datos a graficar ------------

tiempo = [t_r t_r(end)+dt:dt: t_r(end) + (size(theta_f,2))*dt];



%% Pruebas de velocidad lineal 

rutas_x = [x_r x];
rutas_y = [y_r y]; 

limit = 0.8;
s = 0;

while true
    velocidad_x = gradient(rutas_x)./gradient(tiempo);
    velocidad_y = gradient(rutas_y)./gradient(tiempo);
    
    velocidad = sqrt(velocidad_x.^2 + velocidad_y.^2);
    
    velocidad_m = rms(velocidad( size(x_r,2):end ));

    if (round(velocidad_m,3) <= 10.0)     
       break 
    end
    
    tiempo = [t_r t_r(end)+dt:dt: t_r(end) + (size(theta_f,2))*dt]*velocidad_m/10;

end

disp(velocidad_m);

figure
plot(tiempo, velocidad);

%% Simulación y cálculo de ángulos

% v = VideoWriter('Movimiento2.avi');
% open(v);



figure('Name','Movimiento','NumberTitle','off');

for j = 1:60
    plot(x,y,'b');hold on;

    plot([-10 max(x)+15], [y0 y0], 'k--'); hold on;
    plot([min(x) min(x)], [-25 -25+(max(x)+10)+15], 'k--'); hold on;


    theta1_i = -1.4122; %theta1(1)
    theta2_i = 2.8244; %theta2(1)

    plot([0 L1 * cos(theta1_i)], [0 L1 * sin(theta1_i)], 'ro-');hold on;
    plot([L1 * cos(theta1_i)  L1 * cos(theta1_i) + L2 * cos(theta2_i + theta1_i)], [L1 * sin(theta1_i) L1 * sin(theta1_i) + L2 * sin(theta2_i + theta1_i)], 'bo-');hold on;

    title('Recorrido')
    ylabel('Distancia (cm)')
    xlabel('Distancia (cm)')

    axis([-10 max(x)+15 -25 -25+(max(x)+10)+15]);

    % frame = getframe(gcf);
    % writeVideo(v,frame);

pause(0.01);

end

hold off;

for ii = 1:2
    
    if ii == 1
        rutas = [x_r; y_r];
    elseif ii == 2
        rutas = [x; y];
    end
    
    for i = 1:size(rutas,2)
        [theta1, theta2] = inversa(rutas(1,i), rutas(2,i), L1, L2);
        angulo_e1(i + size(x_r,2)*(ii>1)) = theta1;
        angulo_e2(i + size(x_r,2)*(ii>1)) = theta2;
        if ii == 2
            plot(x_r, y_r,'r--'); hold on;
        end
        
        plot(rutas(1,1:i), rutas(2,1:i),'g-');hold on;
        plot([0 L1 * cos(theta1)], [0 L1 * sin(theta1)], 'ro-');hold on;
        plot([L1 * cos(theta1)  L1 * cos(theta1) + L2 * cos(theta2 + theta1)], [L1 * sin(theta1) L1 * sin(theta1) + L2 * sin(theta2 + theta1)], 'bo-');hold off;
        
        axis([-10 max(x)+15 -25 -25+(max(x)+10)+15]);
        
        title('Recorrido')
        ylabel('Distancia (cm)')
        xlabel('Distancia (cm)')
        
%         for l = 1:2
%            frame = getframe(gcf);
%             writeVideo(v,frame); 
%         end
        

        pause(0.001)
        
        %         if i==nt
        %             plot([x(in)],[y(in)],'x','MarkerSize',10,'MarkerEdgeColor','k');hold on
        %         end
    end
    
end

% close(v);

%% Análisis cinemático


% Velocidad angular--------------------------------------------------------

omega_1 = gradient(angulo_e1)./gradient(tiempo);
omega_2 = gradient(angulo_e2)./gradient(tiempo) + omega_1;


% Velocidad lineal---------------------------------------------------------

angulo_1_v = angulo_e1(1:end);
angulo_2_v = angulo_e2(1:end) + angulo_1_v;

vxA = L1*omega_1.*sin(angulo_1_v);
vyA = L1*omega_1.*cos(angulo_1_v);

vxBA = L2*omega_2.*sin(angulo_2_v);
vyBA = L2*omega_2.*cos(angulo_2_v);

vxB = vxA + vxBA;
vyB = vyA + vyBA;



% Aceleración angular------------------------------------------------------


alpha_1 = gradient(omega_1)./gradient(tiempo(1:end));
alpha_2 = gradient(omega_2)./gradient(tiempo(1:end)) + alpha_1;




% Aceleración lineal ------------------------------------------------------

angulo_1_a = angulo_e1(1:end);
angulo_2_a = angulo_e2(1:end) + angulo_1_a;

axA = L1*alpha_1.*sin(angulo_1_a) - L1*omega_1.^2.*cos(angulo_1_a);
ayA = L1*alpha_1.*cos(angulo_1_a) - L1*omega_1.^2.*sin(angulo_1_a);

axBA = L2*alpha_2.*sin(angulo_2_a) - L2*omega_2.^2.*cos(angulo_2_a);
ayBA = L2*alpha_2.*cos(angulo_2_a) - L2*omega_2.^2.*sin(angulo_2_a);

axB = axA + axBA;
ayB = ayA + ayBA;

axG1 = axA/2;
ayG1 = ayA/2;

axG2 = axA + axBA/2;
ayG2 = ayA + ayBA/2;




% Gráficas ---------------------------------------------------------

figure('Name','Cinemática','NumberTitle','off');
tiledlayout(3,2)

ax1 = nexttile;
plot(ax1,tiempo, angulo_e1*180/pi);
title(ax1,'Theta_1')
ylabel(ax1,'Ángulo (°)')
xlabel(ax1,'Tiempo (s)')

ax2 = nexttile;
plot(ax2,tiempo, angulo_e2*180/pi + angulo_e1*180/pi);
title(ax2,'Theta_2')
ylabel(ax2,'Ángulo (°)')
xlabel(ax2,'Tiempo (s)')

ax3 = nexttile;
plot(ax3,tiempo(1:end),omega_1);
title(ax3,'Omega_1')
ylabel(ax3,'Velocidad angular (rad/s)')
xlabel(ax3,'Tiempo (s)')

ax4 = nexttile;
plot(ax4,tiempo(1:end),omega_2);
title(ax4,'Omega_2')
ylabel(ax4,'Velocidad angular (rad/s)')
xlabel(ax4,'Tiempo (s)')

ax5 = nexttile;
plot(ax5,tiempo(1:end),alpha_1);
title(ax5,'Alpha_1')
ylabel(ax5,'Aceleración angular (rad/s^2)')
xlabel(ax5,'Tiempo (s)')

ax6 = nexttile;
plot(ax6,tiempo(1:end),alpha_2);
title(ax6,'Alpha_2')
ylabel(ax6,'Aceleración angular (rad/s^2)')
xlabel(ax6,'Tiempo (s)')


%% Análisis de fuerzas

% Momentos de inercia -----------------------------------------------------
Den_Acrilico = 1.18 / 1000; %g/cm^3

m1 = Den_Acrilico * L1 * 4 * 1 ; % kg

m2 = Den_Acrilico * L2 * 4 * 1 ; % kg

IG1 = 1/12 * m1 * L1^2 / 100^2; % kg*m^2

IG2 = 1/12 * m2 * L2^2 / 100^2; % kg*m^2

g = 9.81; % cm/s^2

% Ecuaciones --------------------------------------------------------------

% Se considera como eslabón 0 a la bancada

% Eslabon #1
% F01x + F21x = m1 * axG1
% F01y + F21y = m1 * ayG1 + m1 * g
% 
% T01 + (R01x * F01y - R01y * F01x) + (R21x * F21y - R21y * F21x) = IG1 * alpha_1

% Eslabón #2
% F12x = m2 * axG2
% F12y = m2 * ayG2 + m2 * g
% 
% T12 + (R12x * F12y - R12y * F12x) = IG2 * alpha_2

%Matricez para hallar incógnitas

% X = [F01x;
%     F01y;
%     T01;
%     F21x;
%     F21y;
%     T12];
for iii = 1:size(angulo_e1,2)
    
    R01x = -L1/2*cos(angulo_e1(iii)) / 100; % m
    R01y = -L1/2*sin(angulo_e1(iii)) / 100; % m

    R21x = L1/2*cos(angulo_e1(iii)) / 100; % m
    R21y = L1/2*sin(angulo_e1(iii)) / 100; % m

    R12x = -L2/2*cos(angulo_e1(iii) + angulo_e2(iii)) / 100; % m
    R12y = -L2/2*sin(angulo_e1(iii) + angulo_e2(iii)) / 100; % m


    A = [1      0       0       1       0       0;
        0       1   	0       0       1       0;
        -R01y   R01x    1       -R21y   R21x    0;
        0       0       0       -1      0       0;
        0       0       0       0       -1      0;
        0       0       0       R12y    R12x    1];


    B = [m1*axG1(iii);
        m1 * ayG1(iii) + m1 * g;
        IG1 * alpha_1(iii);
        m2 * axG2(iii);
        m2 * ayG2(iii) + m2 * g;
        IG2 * alpha_2(iii)];
    
%     plot([0 -R01x], [0 -R01y], 'ro-');hold on;
%     plot([-R01x  -R01x+R21x], [-R01y  -R01y+R21y], 'bo-');hold on;
%     plot([-R01x+R21x -R01x+R21x+(-R12x)], [-R01y+R21y -R01y+R21y+(-R12y)], 'bo-');hold off;
    
%     axis([-10 max(x)+15 -25 -25+(max(x)+10)+15]);
    
%     pause(0.0001)
    
    X = A\B;
    
    F01x(iii) = X(1);
    F01y(iii) = X(2);
    T01(iii) = X(3);
    F21x(iii) = X(4);
    F21y(iii) = X(5);
    T12(iii) = X(6);
     
end

figure;
tiledlayout(2,1)

ax1 = nexttile;
plot(ax1,tiempo,T01);
title(ax1,'Torque_1')
ylabel(ax1,'Torque (N-m)')
xlabel(ax1,'Tiempo (s)')

ax2 = nexttile;
plot(ax2,tiempo,T12);
title(ax2,'Torque_2')
ylabel(ax2,'Torque (N-m)')
xlabel(ax2,'Tiempo (s)')

% figure;
% plot(X(3,:));

%% Código para selección de servmotores
rad2rpm = 9.5493;

% Velocidad angular maxima
omega_1_max = max(omega_1); % rad/s
omega_2_max = max(omega_2); % rad/s

omega_1_rms = rms(omega_1); % rad/s
omega_2_rms = rms(omega_2); % rad/s

% Se realiza un suavizado del vector para eliminar picos
T01_peak = max(T01);
T01_rms = rms(T01);

T12_peak = max(T12);
T12_rms = rms(T12);

% Potencia
power_1 = omega_1.*T01; % w
power_1_peak = max(power_1); % w
power_1_rms = rms(power_1); % w

power_2 = omega_2.*T12; % w
power_2_peak = max(power_2); % w
power_2_rms = rms(power_2); % w

% Intertia ratio

motor_names = ["DC022C-1" "DC022C-2" "DC022C-3" "DC026C-1" "DC026C-2" "DC026C-3" "DC030B-1" "DC030B-2" "DC030B-3" "DC030C-1" "DC030C-2" "DC030C-3" "DC040B-1" "DC040B-2" "DC040B-3" "DC040B-4" "DC040B-5" "DC040B-6" "DC054B-1"	"DC054B-2" "DC054B-3" "DC054B-4" "DC054B-5" "DC054B-6" "DC054B-7" "DC083A-1" "DC083A-2" "DC083A-3" "DC083A-4"];
motor_inertias = [0.00000052 0.00000068	0.00000081 0.00000099	0.0000012	0.0000016 0.00000099	0.0000012	0.0000016 0.000002	0.0000037	0.0000058 0.0000019	0.0000032	0.0000042	0.0000056	0.0000071	0.0000085 0.000011	0.000016	0.000021	0.000026	0.000031	0.000037	0.000047 1.31E-04 	2.27E-04 	3.30E-04 	4.33E-04];
motor_continuos_torques = [0.0057	0.0093	0.014 0.014	0.017	0.022 0.011	0.014	0.018 0.019	0.041	0.060 0.017	0.033	0.043	0.049	0.067	0.081 0.071	0.099	0.15	0.18	0.22	0.26	0.35 0.53 	0.85  	1.20 	1.59];
motor_peak_torques = [0.018	0.037	0.066 0.059	0.084	0.13 0.045	0.065	0.10 0.068	0.22	0.36 0.086	0.20	0.26	0.32	0.40	0.50 0.39	0.67	1.0	1.3	1.4	1.8	2.6 2.65 	4.24  	6.00 	7.94];
IG1_O1 = IG1 + m1*L1/2*1e-2;
IG2_O2 = IG2 + m2*L2/2*1e-2;
IG2_O1 = IG2 + m2*(L1+L2/2);
intertia_ratios_motor_1 = (IG1_O1+IG2_O1)./motor_inertias;
intertia_ratios_motor_2 = (IG2_O2)./motor_inertias;

intertia_ratio_aparente = 5;

N_motor_1 = sqrt(intertia_ratios_motor_1/intertia_ratio_aparente);
N_motor_2 = sqrt(intertia_ratios_motor_2/intertia_ratio_aparente);

T_aparante_motor_1 = T01_rms./N_motor_1;
T_aparante_motor_2 = T12_rms./N_motor_2;

factores_seguridad_torque_motor_1 = motor_continuos_torques./T_aparante_motor_1;
factores_seguridad_torque_motor_2 = motor_continuos_torques./T_aparante_motor_2;

seleccionado_1 = 20;
N_catalogo_1 = 218.4;
omega_max_motor_1 = 6000;

disp("Motor_1 seleccionado: " + motor_names(seleccionado_1))
if omega_1_max*N_catalogo_1 < omega_max_motor_1/rad2rpm
    disp("----El motor cumple el requisito de velocidad angular");
    disp("--------Velocidad angular maxima alcanzada: " + omega_1_max*N_catalogo_1*rad2rpm + " RPM")
    disp("--------Velocidad angular limite: " + omega_max_motor_1 + " RPM")
end

if T01_rms/N_catalogo_1 < motor_continuos_torques(seleccionado_1)
    disp("----El motor cumple el requisito de torque");
    disp("--------Torque rms aplicado: " + T01_rms/N_catalogo_1 + " N m")
    disp("--------Torque rms limite: " + motor_continuos_torques(seleccionado_1) + " RPM")
end

disp("La relación de inercias es " + (intertia_ratios_motor_1(seleccionado_1)/N_catalogo_1^2))


%% Código para selección de rodamientos

%Con rodamiento SKF 618/7   % Rodamiento de dimetro interno 7mm
%C0=0.26*1000; %N
%C= 0.78*1000; %N

%Con rodamiento SKF 617/7   % Rodamiento de dimetro interno 7mm
C0=0.04*1000; %N
C= 0.26*1000; %N

P01 = sqrt(F01x.^2+F01y.^2);
P21 = sqrt(F21x.^2+F21y.^2);

P01f=rms(P01)*1.3*5; %1.3 de correa dentada, 5 como factor de seguridad propio
P21f = rms(P21)*1.3*5;

s01=C0/P01f;
s21=C0/P21f;

L01= (C/P01f)^3; %millones de revoluciones
L21= (C/P21f)^3; %millones de revoluciones

function[theta1, theta2] = inversa(Px,Py,L1,L2)
    r1 = sqrt(Px^2 + Py^2);
    M = (r1^2 - L1^2 - L2^2)/(2 * L1 * L2);
    B1 = atan2(Py,Px);
    M2 = (r1^2 + L1^2 - L2^2)/(2 * r1 * L1);

    theta1 = B1 - acos(M2);
    theta2 = acos(M);
end
