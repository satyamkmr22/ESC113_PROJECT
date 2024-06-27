clear();

% Tank parameters
n = input("Enter number of tanks:"); %Number of tanks
initial_temp = input("Enter the initial temperature of the oil in all tanks:"); % Initial temperature of the oil in all tanks (in °C)
oil_mass = 1000; % Mass of oil in each tank (in kg)
flow_rate = 100; % Oil flow rate (in kg/min)
heat_capacity = 2.0; % Heat capacity of the oil (in KJ/kg·°C)
UA = 10; % Heat transfer coefficient times area of the coil (in kJ/min·°C)

% Steam parameters
steam_temp = input("Enter the temperature of the steam:"); % Temperature of the saturated steam (in °C)


h=0.5; %time_step


tolerance = 0.1; % Tolerance for convergence
max_iterations = 1000; % Maximum number of iterations


% Explicit Euler

T=ones(n,1)*initial_temp;
T_values= zeros(max_iterations,n); %to store values of Temp at each time_step

for i = 1:max_iterations

    T_values(i,:)=T; %storing value of T

    steam_temp=ones(n,1).*steam_temp;
    Q=UA*(steam_temp-T); % Calculate heat transfer rate for each tank
    
    delta_T=zeros(n,1);
    delta_T(1)=(flow_rate * (initial_temp - T(1)) + Q(1)) / (oil_mass * heat_capacity); % Calculate temperature change in each tank

    %to show convergence
    figure(1);
    title("Showing Temperatures reaching steady state by convergence");
    xlabel("Tank Number");
    ylabel("Change in Temperature");
    plot(1,delta_T(1),'+');
    hold on;
    for j=2:n
        delta_T(j)=(flow_rate * (T(j-1) - T(j)) + Q(j)) / (oil_mass * heat_capacity);
        plot(j,delta_T(j),'+');
        hold on;
    end
    
    % Update temperatures
    T=T+h*delta_T;

    %plotting values of T at each time_step
    figure(2);
    title("Temperature in each tank");
    xlabel("Tank Number");
    ylabel("Temperature");
    for j=1:n
        plot(j,T(j),'+');
        hold on;
    end

    
    % Check convergence
    if vecnorm(delta_T) < tolerance
        break;
    end
end

%Deleting extra rows
for j=i+1:max_iterations
    T_values(i+1,:) = [];
end

%Calculation of Steady State Time
SteadyStateTime = i*h; 
fprintf("Steady State Time: %d seconds\n", SteadyStateTime);

% Display the results
fprintf('Steady State Temperatures:\n');
for j=1:n
    fprintf('Tank %d: %.2f°C\n',j,T(j));
end