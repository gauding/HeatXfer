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
dx = 0.001; % dx for 10 steps in x (0.45m/10) in meters

x = 0:dx:width; % total x array

% y direction
dy = 0.001; % dy for 10 steps in y (0.01/10) in meters

y = 0:dy:thickness; % total y array


% this part gets fun
% find the most restrictive time step from the stability analysis
rho = [rho_Al; rho_CI; rho_Br];
c = [Cp_Al; Cp_CI; Cp_Br];
k = [K_Al; K_CI; K_Br];

% lumped capacitance things
a = h.*(2./width + 1./thickness)./(rho.*c);
b = q_flux./(rho.*c.*thickness);
% Bi calcs
Bi = h * thickness ./ k;


[dt] = find_dt_please(rho, c, k, h, dx, dy); % function to find the most restrictive dt

% Fo
Fo = (k.*dt)./(rho.*c.*dx.^2); %k./rho./c.*dt./dx^2; % fourier number in the x direction -- may need to change to y


%% LOOPS!! -- for loop for materials, while loop for time/converging on answer, for loop for x, for y

% first for loop -- materials (3)
% set up output arrays -- these may not all be used for finite diff, but
% all are on the spreadsheet
time = zeros(length(rho),1); % time to get to 190 deg C on cook surface (s)
index = zeros(length(rho),1); % index for time to 190 deg C
time_LC = zeros(length(rho),1); % time to get to 190 deg C on cook surface (s) (for lumped capacitance)
T_ss_LC = zeros(length(rho),1); % steady state temp for lumped capacitance
T_ss = zeros(length(rho),1); % steady state temp for finite diff
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


max_count = [250000000, 100000000, 20000000];

%max_count = 1e5; % commenting out to reduce conflict -- AMW

% actually start the loop now
for i = 1:length(rho)
    
    alpha = k(i) / rho(i) / c(i); % alpha for the material
    beta = h / rho(i) / c(i); % made up constant that shows up alot, prob something to do with Fo number lol
    
    
    %count = 1; % reset counter for while loop so we dont explode
    
    
    % need to create T array
    T = zeros(length(y), length(x)); % create array of zeros
    T = T + Ti; % add the initial condition
    
    %T_top_mid = zeros(max_count(i)+1, 1); 
    T_top_mid(1) = T(ceil(length(y)),ceil(length(x)/2)); % find the initial top mid value
    %err = abs(T_top_mid - T_goal); % stop condition
    % start while loop ( actually using for loop because idexing is used
    % for keeping track of the temp)
    %for count = 1:max_count(i)
    
     
    count=1; 
    error=1; 
    while count<max_count(i) && error>1e-7
        % for loop for x direction
        for m = 1:length(x)
            
            
            % for loop for y direction
            for n = 1:length(y)
                
                % where are we? -- find the right finite diff based on
                % position
                % bottom left corner
                
                %GG: matlab indicies are weird! The first index is the row
                %and the second is the column, so these need to be T(n, m) 
                if m==1 && n==1 % bottom left corner
                    T(n,m) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) - (2*beta*dt(i)/dx) )*T(n,m) + (2*dt(i)/rho(i)/c(i)/dy*q_flux) + (2*alpha*dt(i)/dy^2*T(n+1,m)) + (2*alpha*dt(i)/dx^2*T(n,m+1)) + (2*beta*dt(i)/dx*T_amb);
                % bottom between corners
                elseif m>1 && m<length(x) && n==1 % bottom between corners
                    T(n,m) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(n,m) + (2*dt(i)/rho(i)/c(i)/dy*q_flux) + (2*alpha*dt(i)/dy^2*T(n+1,m)) + (alpha*dt(i)/dx^2*(T(n,m+1) + T(n,m-1)));
                % bottom right corner
                elseif m==length(x) && n==1 % bottom right corner
                    T(n,m) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) - (2*beta*dt(i)/dx) )*T(n,m) + (2*dt(i)/rho(i)/c(i)/dy*q_flux) + (2*alpha*dt(i)/dy^2*T(n+1,m)) + (2*alpha*dt(i)/dx^2*T(n,m-1)) + (2*beta*dt(i)/dx*T_amb);
                % left side between corners
                elseif m == 1 && n>1 && n<length(y) % left side between corners
                    T(n,m) = (1 - (2*beta*dt(i)/dx) - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(n,m) + (2*beta*dt(i)/dx*T_amb) + (alpha*dt(i)/dy^2*(T(n+1,m) + T(n-1,m))) + (2*alpha*dt(i)/dx^2*T(n,m+1));
                % right side between corners
                elseif m==length(x) && n>1 && n<length(y) % right side between corners
                    T(n,m) = (1 - (2*beta*dt(i)/dx) - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(n,m) + (2*beta*dt(i)/dx*T_amb) + (alpha*dt(i)/dy^2*(T(n+1,m) + T(n-1,m))) + (2*alpha*dt(i)/dx^2*T(n,m-1));
                % top left corner
                elseif m == 1 && n==length(y) % top left corner
                    T(n,m) = (1 - (2*beta*dt(i)) - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2))*T(n,m) + (2*alpha*dt(i)/dy^2*T(n-1,m)) + (2*alpha*dt(i)/dx^2*T(n,m+1)) + (2*beta*dt(i)/dx*T_amb) + (2*beta*dt(i)/dy*T_amb);
                % top between corners
                elseif m>1 && m<length(x) && n==length(y) % top between corners
                    T(n,m) = (1 - (2*beta*dt(i)/dy) - (2*alpha*dt(i)/dx^2) -(2*alpha*dt(i)/dy^2))*T(n,m) + (2*beta*dt(i)/dy*T_amb) + (alpha*dt(i)/dy^2*(T(n,m+1) + T(n,m-1))) + (2*alpha*dt(i)/dx^2*T(n-1,m));
                % top right corner
                elseif m == length(x) && n==length(y) % top right corner
                    T(n,m) = (1 - (2*beta*dt(i))- (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2))*T(n,m) + (2*alpha*dt(i)/dy^2*T(n-1,m)) + (2*alpha*dt(i)/dx^2*T(n,m-1)) + (2*beta*dt(i)/dx*T_amb) + (2*beta*dt(i)/dy*T_amb);
                % interior points
                elseif m>1 && m<length(x) && n>1 && n<length(y) % interior points
                    T(n,m) = (1 - (2*alpha*dt(i)/dy^2) - (2*alpha*dt(i)/dx^2) )*T(n,m) + (alpha*dt(i)/dy^2*(T(n+1,m) + T(n-1,m))) + (alpha*dt(i)/dx^2*(T(n,m+1) + T(n,m-1)));
                else
                    fprintf('You messed up your indexing, dummy\n');
                end
                    
                    
            end
        end
        
        
        
        % update while loop conditions
        T_top_mid(count) = T(ceil(length(y)), ceil(length(x)/2));
        %err = abs(T_top_mid(count+1) - T_goal);
        %if T_top_mid(count+1) > T_goal
        %    break
        %end
        if count>1 && T_top_mid(count) > (Ti + 10)
            
            error= abs(T_top_mid(count)- T_top_mid(count-1)); 
        end 
        count=count+1; 
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
    %plot(t,T_top_mid(index(i)), '-k', 'linewidth',2)
    hold on
    plot(t(1:count-1) , T_top_mid(1:count-1), '-k', 'linewidth',2)
    grid on 
    title(plot_title{i})
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    
    % Lumped Capacitance
    time_LC(i) = 1/a(i) * log((Ti-T_amb-(b(i)/a(i)))/(T_goal-T_amb-(b(i)/a(i)))); % time to desired temp using lumped capacitance method
    

    % Lumped Capacitance steady state temp
    T_ss_LC(i) = T_amb + (b(i)/a(i));
    % finite diff steady state temp
    T_ss(i) = T_top_mid(end);

    %% Lumped Capacitance calculations 
    
    %Find the Biot number 
    
    %Bi(i)= h*thickness./(k(i)); % commenting out bc it is calc'ed further up and to reduce conflict --AMW
    
    %lumped Capacitance doesn't actually work for the brick because of the
    %low k 
    
    %variables used in the lumped capacitance solution 
    Asc= width+2*thickness ;
    Asq= width ; 
    
    
    Vol= width*thickness; 
    
    %a= (h*Asc)/ (rho(i)* Vol* c(i)) ; % commenting out bc it is calc'ed further up and to reduce conflict --AMW
    
    %b= (q_flux*Asq) / (rho(i)*Vol*c(i)); % commenting out bc it is calc'ed further up and to reduce conflict --AMW
    
    
    %fcn=@(time) exp(-a*time)+((b/a)/(Ti-T_amb))*(1-exp(-a*time)) -((T_goal-T_amb)/ (Ti-T_amb)) ; 
    %time(i)= fzero(fcn, 1) % commenting out for now to reduce conflict, this did not want to run when I tried, even with changing variable names -- AMW
    
    %{
    figure(4) 
    plot ( T(:,ceil(length(y)/2) ), y) 
    hold on 
    xlabel("Temperature (K) " ) 
    ylabel ("Thickness of Plancha") 
    %}

end

% print the times
fprintf('Time to 190 deg C for Al: %f seconds\n',time(1));
fprintf('Time to 190 deg C for Cast Iron: %f seconds\n',time(2));
fprintf('Time to 190 deg C for Ceramic: %f seconds\n',time(3));
fprintf('LC Time to 190 deg C for Al: %f seconds\n',time_LC(1));
fprintf('LC Time to 190 deg C for Cast Iron: %f seconds\n',time_LC(2));
fprintf('LC Time to 190 deg C for Ceramic: %f seconds\n',time_LC(3));

%% More Lumped Capacitance Calcs
Q_conv_loss_goal = h*(2*thickness + width)*(T_goal - T_amb); % heat lost to convection at operating temp (W/m)
Q_conv_loss_ss = h*(2*thickness + width)*(T_ss_LC(1) - T_amb); % heat lost to convection at steady state temp (W/m)
Q_absorbed_goal = rho.*thickness.*width.*c.*(T_goal-Ti)/1000; % heat absorbed by plancha at operating temp (kJ/m)
q_flux_goal = (T_goal- T_amb) * (h*thickness*(2/width + 1/thickness)); % heat flux required to maintain goal temp
T_half_t = (Ti - T_amb - b./a).*exp(-a.*time_LC./2) + T_amb + b./a; % temp at half time to 190 deg C
time_LC_half_h = 1./(a./2) .* log((Ti-T_amb-(b./(a./2)))./(T_goal-T_amb-(b./(a./2)))); % time to 190 deg C at half h


%% Write the data into the spreadsheet
data_out = zeros(17,3);
data_out(1,:) = thickness; % thickness of plancha
data_out(2,:) = k; % conductivity
data_out(3,:) = rho; % density
data_out(4,:) = c; % specific heat
data_out(5,:) = 0;%'Textbook'; % source for properties -- need to figure out how to do text - maybe cell array?
data_out(6,:) = Bi; % Bi
data_out(7,:) = dx; % delta x
data_out(8,:) = Fo; % Fo
data_out(9,:) = dt; % delta t
data_out(10,:) = time_LC; % time to 190 deg C
data_out(11,:) = Q_absorbed_goal; % absorbed heat at operating temp
data_out(12,:) = T_ss_LC-273; % steady state temp at given flux (deg C)
data_out(13,:) = Q_conv_loss_goal; % heat lost to convection at operating temp (W/m)
data_out(14,:) = Q_conv_loss_ss; % heat lost to convection at steady state temp (W/m)
data_out(15,:) = q_flux_goal; % heat flux required to maintain goal temp as ss (W/m^2)
data_out(16,:) = T_half_t-273; % temp at top mid at half time (deg C)
data_out(17,:) = time_LC_half_h; % time to 190 deg C for half h value (s)

% this might not work on earlier versions of MATLAB -- let me know and we
% can figure it out
%
filename = 'Project 1 Results Template.xlsx'; % name of spreadsheet file
writematrix(data_out,filename,'Sheet','Sheet1','Range','B4:D20') % writes data to spreadsheet
writecell({'Textbook','Textbook','Textbook'},filename,'Sheet','Sheet1','Range','B8:D8') % writes source of properties to spreadsheet
%
%% TO DO:

% plot for temp at centerline at one-half time -- Genevieve claimed
% plot for temp over time with LC? -- finite diff and lumped capacitance
% are different



%% Genevieve's stuff that got put at the bottom during the merge
q_loss= h*Asc*(T_goal-T_amb) 

%% making the fun plot 


