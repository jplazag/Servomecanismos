close all
clear

run ("referencia.m");

run("ModeloControl.m");

out=sim("LazosControl.slx",tiempo(end));