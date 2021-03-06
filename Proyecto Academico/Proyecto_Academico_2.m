close all

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



dt = 0.01;
tf = 5; %Maximo 35 segundos y minimo 3 segundos
nt = tf/dt;

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


% fplot(diff(a),[0-dt+ punto_i - theta0_f, 2*pi-dt+ punto_i - theta0_f]);
% hold on;
% plot(linspace(0-dt+ punto_i - theta0_f,2*pi -dt+ punto_i - theta0_f,nt-1),flip(diff(x)./diff(theta_f)));

derivada = double(vpa(subs(diff(b),time,p))/vpa(subs(diff(a),time,p)));

prueba = round(int(sqrt(diff(a)^2 + diff(b)^2),0+ punto_i - theta0_f-dt,2*pi+ punto_i - theta0_f-dt));
% fplot(diff(b)/diff(a),[1 size(vt,2)]);hold on


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

% figure;
% plot(vt_c);

% Trayectoria recta -------------------------------------------------------
x_r = linspace(x_inicio-dt,x_c(1)-dt,nt);


% Arreglo de las rutas que componen la trayectoria a recorrer -------------
% rutas = [size(x_r); size(x_c); size(x)];

% Preparación de variables para almacenar los datos a graficar ------------

tiempo = 0:dt:(size(x_r,2) + size(x_c,2) + size(x,2)-dt)*dt ;

v_p = [0.0001 0.0001 0.0001];

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

for ii = 1:3
    
    if ii == 1
        rutas = [x_r; zeros(size(x_r))];
    elseif ii == 2
        rutas = [x_c; y_c];
    elseif ii == 3
        rutas = [x; y];
    end
    disp(size(rutas,2))
    for i = 1:size(rutas,2)
        [theta1, theta2] = inversa(rutas(1,i), rutas(2,i), L1, L2);
        angulo_e1(i + size(x_r,2)*(ii>1) + size(x_c,2)*(ii>2) ) = theta1;
        angulo_e2(i + size(x_r,2)*(ii>1) + size(x_c,2)*(ii>2) ) = theta2;
        disp(i + size(x_r,2)*(ii>1) + size(x_c,2)*(ii>2) )
        if ii == 2
            plot(x_r, zeros(size(x_r)),'r--'); hold on;
        elseif ii == 3
            plot(x_r, zeros(size(x_r)),'r--'); hold on;
            plot(x_c, y_c,'r--'); hold on;
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
        

        pause(v_p(ii))
        
        %         if i==nt
        %             plot([x(in)],[y(in)],'x','MarkerSize',10,'MarkerEdgeColor','k');hold on
        %         end
    end
    
end

% close(v);


%% Análisis cinemático


% Posición angular---------------------------------------------------------

figure('Name','Ángulos','NumberTitle','off');
tiledlayout(2,1)

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


% Velocidad angular--------------------------------------------------------

omega_1 = [0 diff(angulo_e1)/dt];
omega_2 = [0 diff(angulo_e2)/dt] + omega_1;

figure;
tiledlayout(2,1)

ax1 = nexttile;
plot(ax1,tiempo(1:end),omega_1);
title(ax1,'Omega_1')
ylabel(ax1,'Velocidad angular (rad/s)')
xlabel(ax1,'Tiempo (s)')

ax2 = nexttile;
plot(ax2,tiempo(1:end),omega_2);
title(ax2,'Omega_2')
ylabel(ax2,'Velocidad angular (rad/s)')
xlabel(ax2,'Tiempo (s)')

% Velocidad lineal---------------------------------------------------------

angulo_1_v = angulo_e1(1:end);
angulo_2_v = angulo_e2(1:end) + angulo_1_v;

vxA = L1*omega_1.*sin(angulo_1_v);
vyA = L1*omega_1.*cos(angulo_1_v);

vxBA = L2*omega_2.*sin(angulo_2_v);
vyBA = L2*omega_2.*cos(angulo_2_v);

vxB = vxA + vxBA;
vyB = vyA + vyBA;

figure();
tiledlayout(2,1)

ax1 = nexttile;

plot(ax1,-vxB(101:149));
hold on;
plot(ax1,vx);
hold on;
legend(ax1,"vxB","vx");


ax2 = nexttile;

plot(ax2,vy);
hold on;
plot(ax2,vyB(101:149));
legend(ax2,"vy","vyB");



% Aceleración angular------------------------------------------------------


alpha_1 = [diff(omega_1)/dt 0];
alpha_2 = [diff(omega_2)/dt 0] + alpha_1;


figure;
tiledlayout(2,1)

ax1 = nexttile;
plot(ax1,tiempo(1:end),alpha_1);
title(ax1,'Alpha_1')
ylabel(ax1,'Aceleración angular (rad/s^2)')
xlabel(ax1,'Tiempo (s)')

ax2 = nexttile;
plot(ax2,tiempo(1:end),alpha_2);
title(ax2,'Alpha_2')
ylabel(ax2,'Aceleración angular (rad/s^2)')
xlabel(ax2,'Tiempo (s)')

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
figure;
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
    
    plot(tiempo(1:iii),F21x(1:iii),'r-'); hold on;
    
    
end



% figure;
% plot(X(3,:));
%% Codigo para seleccion de rodamientos

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


