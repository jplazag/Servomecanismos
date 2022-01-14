close all
clear all

% syms t
% a = scale*(r0-ra*sin(n*(t+theta0))).* sin(t)+ x0;
% b = scale*(r0-ra*sin(n*(t+theta0))).* cos(t)+ y0;
% 
% t = [pi/2 7*pi/6 11*pi/6];
% picos_a = subs(a);
% picos_b = subs(b);
% 
% da = vpa(subs(a,t,2));
% db = vpa(subs(b,t,2));

rad2rpm = 9.5493;
dt = 0.01;
tf = 5; %Maximo 35 segundos y minimo 3 segundos
nt = tf/dt;

r0 = 10;
ra = r0*3/10;
n = 3;
x0 = 25;
y0 = 15;
theta0_f = 30*pi/180;
theta_f = flip(0:dt:2*pi); %linspace(0,2*pi,nt);
r = r0-ra*sin(n*(theta_f + theta0_f));
scale = 1;

L1 =21.5; %1.0*ceil(max(R))/2;    cm
L2 =21.5; %1.0*ceil(max(R))/2;    cm


t = [pi/2 7*pi/6 11*pi/6];

picos_a = scale*(r0-ra*sin(n*(t + theta0_f))).* sin(t) + x0;
picos_b = scale*(r0-ra*sin(n*(t + theta0_f))).* cos(t) + y0;

vectores_a_picos = sqrt(picos_a.^2 + picos_b.^2);

punto_i = dot(t, floor(vectores_a_picos/max(vectores_a_picos)));

theta_f = flip(0:dt:2*pi) + punto_i - theta0_f;



% Trayectoria de interés --------------------------------------------------

x = scale*0.8*(r0-ra*sin(n*(theta_f + theta0_f))) .* sin(theta_f)+ x0;
y = scale*0.8*(r0-ra*sin(n*(theta_f + theta0_f))) .* cos(theta_f)+ y0;

fun_interes = sqrt( (diff(x)./diff(theta_f)).^2 + (diff(y)./diff(theta_f)).^2 );

arco_interes = sum(fun_interes)*abs(mean(diff(theta_f)));

 
% Variables de interés ----------------------------------------------------
vx = diff(x)/dt;
vxm = mean(abs(vx));

vy = diff(y)/dt;
vym = mean(abs(vy));
vt = sqrt(vx.^2+vy.^2);

vt_media = mean(vt);
vt_maxima = max(vt);

p = punto_i - theta0_f;

% Punto de inicio ---------------------------------------------------------
x_p = scale * 0.8 * (r0 - ra * sin(n*(p + theta0_f)))* sin(p) + x0;
y_p = scale * 0.8 * (r0 - ra * sin(n*(p + theta0_f)))* cos(p) + y0;

% Pendiente de conexión ---------------------------------------------------
syms time

a = scale * 0.8 * (r0 - ra * sin(n*(time + theta0_f)))* sin(time) + x0;
b = scale * 0.8 * (r0 - ra * sin(n*(time + theta0_f)))* cos(time) + y0;

derivada = double(vpa(subs(diff(b),time,p))/vpa(subs(diff(a),time,p)));

% Trayectoria circular ----------------------------------------------------
a_c = 1;
b_c = -2*y_p * (1 + 1/(derivada^2));
c_c = y_p^2 * (1 + 1/(derivada^2));

if derivada > 0
    r_c = (-b_c + sqrt(b_c^2 - 4*a_c*c_c))/(2*a_c);
else
    r_c = (-b_c - sqrt(b_c^2 - 4*a_c*c_c))/(2*a_c);
end

x0_c = x_p - sqrt(r_c^2 - (y_p - r_c)^2);
% y_c = linspace(0, y_p, nt);
% x_c = sqrt(r_c^2 - (y_c - r_c).^2) + x0_c;

p_final_c = asin((y_p - r_c)/r_c);

theta_c = -pi/2:dt*1.1:p_final_c;

x_c = r_c*cos(theta_c) + x0_c;
y_c = r_c*sin(theta_c) + r_c;

x_inicio = min(1.3*0.8*(r0-ra*sin(n*(theta_f + theta0_f))) .* sin(theta_f)+ x0) - 5;

vt_c = ((y_p)/nt)./diff(x_c);

% Trayectoria recta -------------------------------------------------------
x_r = linspace(x_inicio-dt,x_c(1)-dt,nt);


% Arreglo de las rutas que componen la trayectoria a recorrer -------------
% rutas = [size(x_r); size(x_c); size(x)];

% Preparación de variables para almacenar los datos a graficar ------------

tiempo = 0:dt:(size(x_r,2) + size(x_c,2) + size(x,2)-dt)*dt ;

v_p = [0.0001 0.0001 0.0001];

%% Simulación y cálculo de ángulos

for ii = 1:3
    
    if ii == 1
        rutas = [x_r; zeros(size(x_r))];
    elseif ii == 2
        rutas = [x_c; y_c];
    elseif ii == 3
        rutas = [x; y];
    end

    for i = 1:size(rutas,2)
        [theta1, theta2] = inversa(rutas(1,i), rutas(2,i), L1, L2);
        angulo_e1(i + size(x_r,2)*(ii>1) + size(x_c,2)*(ii>2) ) = theta1;
        angulo_e2(i + size(x_r,2)*(ii>1) + size(x_c,2)*(ii>2) ) = theta2;   
    end
    
end

% close(v);


%% Análisis cinemático

% Velocidad angular--------------------------------------------------------

omega_1 = [0 diff(angulo_e1)/dt];
omega_2 = [0 diff(angulo_e2)/dt] + omega_1;

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


alpha_1 = [diff(omega_1)/dt 0];
alpha_2 = [diff(omega_2)/dt 0] + alpha_1;

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
    
    X = A\B;
    
    F01x(iii) = X(1);
    F01y(iii) = X(2);
    T01(iii) = X(3);
    F21x(iii) = X(4);
    F21y(iii) = X(5);
    T12(iii) = X(6); 
end

%% Código para selección de servmotores

% Velocidad angular maxima
omega_1_max = max(omega_1); % rad/s
omega_2_max = max(omega_2); % rad/s

omega_1_rms = rms(omega_1); % rad/s
omega_2_rms = rms(omega_2); % rad/s

% Se realiza un suavizado del vector para eliminar picos
T01_s = medfilt1(T01,3);
T01_peak = max(T01_s);
T01_rms = rms(T01_s);

T12_s = medfilt1(T12,3);
T12_peak = max(T12_s);
T12_rms = rms(T12_s);

% Potencia
power_1 = omega_1.*T01_s; % w
power_1_peak = max(power_1); % w
power_1_rms = rms(power_1); % w

power_2 = omega_2.*T12_s; % w
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


