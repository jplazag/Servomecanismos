close all
clear

%Variables para modelado del trebol
scale = 1;
rota=90; 
vueltas = 4;
v_max = 10

dt = 0.05;

scale_limit = 0.8
r0 = 10;
ra = r0*3/10;
n = 3;
x0 = 25;
y0 = 15;
theta0_f = rota*pi/180;
theta_f = flip(0:dt:2*pi); %linspace(0,2*pi,nt);
r = r0-ra*sin(n*(theta_f + theta0_f));


% R = sqrt(x.^2+y.^2);
% disp(max(R));
% [dmax,in]=max(R);

L1 =21.5; %1.0*ceil(max(R))/2;    cm
L2 =21.5; %1.0*ceil(max(R))/2;    cm


% Selección de la pendiente de la trayectoria recta de aproximación -------

x_inicio = min(1.3*scale_limit*(r0-ra*sin(n*(theta_f + theta0_f))) .* sin(theta_f)+ x0) - 5;

theta_pre_dev = (2*pi:-dt:0) - theta0_f;

x_pre_dev = scale*scale_limit*(r0-ra*sin(n*(theta_pre_dev + theta0_f))) .* sin(theta_pre_dev)+ x0 ;
y_pre_dev = scale*scale_limit*(r0-ra*sin(n*(theta_pre_dev + theta0_f))) .* cos(theta_pre_dev)+ y0 ;

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

theta_f = (2*vueltas*pi:-dt:0) + theta_pre_dev(indice_r) - dt ;

x = scale*scale_limit*(r0-ra*sin(n*(theta_f + theta0_f))) .* sin(theta_f)+ x0;
y = scale*scale_limit*(r0-ra*sin(n*(theta_f + theta0_f))) .* cos(theta_f)+ y0;

% figure
% plot([x_inicio x_final_r],[0 y_final_r], 'r'); hold on;
% plot(x,y);

% Preparación de variables para almacenar los datos a graficar ------------

tiempo = [t_r t_r(end)+dt:dt: t_r(end) + (size(theta_f,2))*dt];

%% Pruebas de velocidad lineal 

rutas_x = [x_r x];
rutas_y = [y_r y]; 

s = 0;

while true
    velocidad_x = gradient(rutas_x)./gradient(tiempo);
    velocidad_y = gradient(rutas_y)./gradient(tiempo);
    
    velocidad = sqrt(velocidad_x.^2 + velocidad_y.^2);
    
    velocidad_m = rms(velocidad( size(x_r,2):end ));

    if (round(velocidad_m,3) <= v_max)     
       break 
    end
    
    tiempo = [t_r t_r(end)+dt:dt: t_r(end) + (size(theta_f,2))*dt]*velocidad_m/v_max;

end

%% Cálculo de ángulos

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
    end
    
end

%% Funciones

function[theta1, theta2] = inversa(Px,Py,L1,L2)
    r1 = sqrt(Px^2 + Py^2);
    M = (r1^2 - L1^2 - L2^2)/(2 * L1 * L2);
    B1 = atan2(Py,Px);
    M2 = (r1^2 + L1^2 - L2^2)/(2 * r1 * L1);

    theta1 = B1 - acos(M2);
    theta2 = acos(M);
end
