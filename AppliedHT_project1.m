%{ 
Genevieve Gaudin and Austin Warren 
Applied Heat transfer Project 1 
%} 

clc 
clear all 
close all 


%% Material properties 

% Aluminum 
rho_Al= 2702; %kg/m^3 
Cp_Al= 903; % j/kgK 
K_Al=237; %W/mK 
alpha_Al= 97.1 *10^-6; %m^2/s


%Cast Iron 
rho_CI= 7854; %kg.m^3 
Cp_CI=434; %J/kgK 
K_CI= 60.5; %W/mK 
alpha_CI=17.7; %m^2/s

%Brick <- this one's properties are at 478K, I don't know why theya were
%liek that in the book 

rho_Br= 2645; %kg/m^3
K_Br= 1.0; %W/mK 
Cp_Br= 960; %J/KgK 

%% 

