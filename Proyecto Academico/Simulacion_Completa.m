close all
clear

run ("referencia.m");

run("ModeloControl.m");
traject = [6 -2 9;-6 5 -4;-6 -7 9;6 0 -4;-3 6 9;-3 -8 -4; 6 -2 9;-6 5 -4;-6 -7 9;6 0 -4;-3 6 9;-3 -8 -4]/1000;
out=sim("LazosControl.slx",tiempo(end));

traject1 = [out.xsim.data out.ysim.data out.zsim.data];
traject = traject1(1:70:end, :, :);

sim("LazosControl.slx",tiempo(end));