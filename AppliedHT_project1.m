%{ 
Genevieve Gaudin and Austin Warren 
Applied Heat transfer Project 1 
%} 
%Fernando is the best goddam  person ever
%so strong
%so masculine
%so incredible

%fernando has a super bald head 
%genevieve is a short meanue poopie head
%fernandos bald head is actually relaly great and hes not self conscious
%about it
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

%% I need to add shit 

%edit 2 that fernando needed to teach me how to do this 


%% Initial Conditions and Givens
h = 15; % convection heat transfer coefficient in W/m^2/K
Ti = 305; % initial temp of plancha in K
T_goal = 190+273; % goal temp of planch top in K
T_amb = 27+273; % ambient air temp in K
q_flux = 4500; % heat flux applied to bottom of planch in W/m^2
width = 0.45; % width given by problem in meters
thickness = 0.01; % thickness given by problem in meters


%% Domain stuff (Genevieve already did this, but not in this edit?) -- whack


% x direction
dx = 0.045; % dx for 10 steps in x (0.45m/10) in meters
x = 0:dx:width; % total x array

% y direction
dy = 0.001; % dy for 10 steps in y (0.01/10) in meters

y = 0:dy:thickness; % total y array


% this part gets fun
% find the most restrictive time step from the stability analysis
rho = [rho_Al; rho_CI; rho_Br];
c = [Cp_Al; Cp_CI; Cp_Br];
k = [K_Al; K_CI; K_Br];


[dt] = find_dt_please(rho, c, k, h, dx, dy); % function to find the most restrictive dt


%% LOOPS!! -- for loop for materials, while loop for time/converging on answer, for loop for x, for y

% first for loop -- materials (3)
% set up output arrays -- these may not all be used for finite diff, but
% all are on the spreadsheet
time = zeros(length(rho),1); % time to get to 190 deg C on cook surface (s)
index = zeros(length(rho),1); % index for time to 190 deg C
tol = [0.001, 0.0000001, 0.0000001];
Q_absorbed = zeros(length(rho),1); % total heat absorbed in this time (kJ/m)
Q_lost = zeros(length(rho),1); % total heat lost to convection (W/m)
Q_lost_ss = zeros(length(rho),1); % heast lost to convection at steady state temp (W/m)
T_ss = zeros(length(rho),1); % steady state temp at stated flux (deg C)
q_flux_new = zeros(length(rho),1); % new flux to required to maintain at 190 deg C (W/m^2)
T_half_time = zeros(length(rho),1); % temp at top center at half the time to 190 deg C
time_half_h = zeros(length(rho),1); % time to 190 deg C with half h value (s)

% plot title cell array to generate plot titles in the loop
plot_title = {'Top Center Temperature vs Time for Aluminum','Top Center Temperature vs Time for Cast Iron','Top Center Temperature vs Time for Ceramic Brick'};

max_count = 100000;
% actually start the loop now
for i = 1:length(rho)
    
    alpha = k(i) / rho(i) / c(i); % alpha for the material
    beta = h / rho(i) / c(i); % made up constant that shows up alot, prob something to do with Fo number lol
    
    
    %count = 1; % reset counter for while loop so we dont explode
    
    % need to create T array
    T = zeros(length(x), length(y)); % create array of zeros
    T = T + Ti; % add the initial condition
    T_top_mid = zeros(max_count+1, 1); 
    T_top_mid(1) = T(ceil(length(x)/2), ceil(length(y))); % find the initial top mid value
    %err = abs(T_top_mid - T_goal); % stop condition
    % start while loop ( actually using for loop because idexing is used
    % for keeping track of the temp)
    for count = 1:max_count
    
        % for loop for x direction
        for m = 1:length(x)
            
            
            % for loop for y direction
            for n = 1:length(y)
                
                % where are we? -- find the right finite diff based on
                % position
                % bottom left corner
                if m==1 && n==1 % bottom left corner
                    T(m,n) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) - (2*beta*dt(i)/dx) )*T(m,n) + (2*dt(i)/rho(i)/c(i)/dy*q_flux) + (2*alpha*dt(i)/dy^2*T(m,n+1)) + (2*alpha*dt(i)/dx^2*T(m+1,n)) + (2*beta*dt(i)/dx*T_amb);
                % bottom between corners
                elseif m>1 && m<length(x) && n==1 % bottom between corners
                    T(m,n) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(m,n) + (2*dt(i)/rho(i)/c(i)/dy*q_flux) + (2*alpha*dt(i)/dy^2*T(m,n+1)) + (alpha*dt(i)/dx^2*(T(m+1,n) + T(m-1,n)));
                % bottom right corner
                elseif m==length(x) && n==1 % bottom right corner
                    T(m,n) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) - (2*beta*dt(i)/dx) )*T(m,n) + (2*dt(i)/rho(i)/c(i)/dy*q_flux) + (2*alpha*dt(i)/dy^2*T(m,n+1)) + (2*alpha*dt(i)/dx^2*T(m-1,n)) + (2*beta*dt(i)/dx*T_amb);
                % left side between corners
                elseif m == 1 && n>1 && n<length(y) % left side between corners
                    T(m,n) = (1 - (2*beta*dt(i)/dx) - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(m,n) + (2*beta*dt(i)/dx*T_amb) + (alpha*dt(i)/dy^2*(T(m,n+1) + T(m,n-1))) + (2*alpha*dt(i)/dx^2*T(m+1,n));
                % right side between corners
                elseif m==length(x) && n>1 && n<length(y) % right side between corners
                    T(m,n) = (1 - (2*beta*dt(i)/dx) - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(m,n) + (2*beta*dt(i)/dx*T_amb) + (alpha*dt(i)/dy^2*(T(m,n+1) + T(m,n-1))) + (2*alpha*dt(i)/dx^2*T(m-1,n));
                % top left corner
                elseif m == 1 && n==length(y) % top left corner
                    T(m,n) = (1 - (2*beta*dt(i)) - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2))*T(m,n) + (2*alpha*dt(i)/dy^2*T(m,n-1)) + (2*alpha*dt(i)/dx^2*T(m+1,n)) + (2*beta*dt(i)/dx*T_amb) + (2*beta*dt(i)/dy*T_amb);
                % top between corners
                elseif m>1 && m<length(x) && n==length(y) % top between corners
                    T(m,n) = (1 - (2*beta*dt(i)/dy) - (2*alpha*dt(i)/dx^2) -(2*alpha*dt(i)/dy^2))*T(m,n) + (2*beta*dt(i)/dy*T_amb) + (alpha*dt(i)/dy^2*(T(m+1,n) + T(m-1,n))) + (2*alpha*dt(i)/dx^2*T(m,n-1));
                % top right corner
                elseif m == length(x) && n==length(y) % top right corner
                    T(m,n) = (1 - (2*beta*dt(i))- (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2))*T(m,n) + (2*alpha*dt(i)/dy^2*T(m,n-1)) + (2*alpha*dt(i)/dx^2*T(m-1,n)) + (2*beta*dt(i)/dx*T_amb) + (2*beta*dt(i)/dy*T_amb);
                % interior points
                elseif m>1 && m<length(x) && n>1 && n<length(y) % interior points
                    T(m,n) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(m,n) + (alpha*dt(i)/dy^2*(T(m,n+1) + T(m,n-1))) + (alpha*dt(i)/dx^2*(T(m+1,n) + T(m-1,n)));
                else
                    fprintf('You messed up your indexing, dummy\n');
                end
                    
                    
            end
            
            
            
            
            
        end
        
        
        
        % update while loop conditions
        T_top_mid(count+1) = T(ceil(length(x)/2), ceil(length(y)));
        %err = abs(T_top_mid(count+1) - T_goal);
        %if T_top_mid(count+1) > T_goal
        %    break
        %end
    end
    
    % find the time
    tot_time = (length(T_top_mid)-1) * dt(i);
    t = 0:dt(i):tot_time;
    F = find(T_top_mid < T_goal);
    index(i) = F(end) + 1;
    time(i) = (index(i)-1)*dt(i);
    %[time(i),index(i)] = max(T_top_mid(find(T_top_mid < T_goal))); % find the time to reach desired temp -- may need to use find function
    
    % plot some stuff
    figure(i)
    plot(t,T_top_mid, '-k', 'linewidth',2)
    hold on
    plot(time(i), T_top_mid(index(i)), '*r', 'linewidth',2)
    title(plot_title{i})
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    
    
    
    
end

% print the times
fprintf('Time to 190 deg C for Al: %f seconds\n',time(1));
fprintf('Time to 190 deg C for Cast Iron: %f seconds\n',time(2));
fprintf('Time to 190 deg C for Ceramic: %f seconds\n',time(3));

%% TO DO:

% lumped capacitance lol
% all of the other answers
% plot for temp at centerline at one-half time
% h is halfed
% steady state stuff
% heat lost to convection -- lumped capacitance?





