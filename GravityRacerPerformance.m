%Ryan Fermoyle
%Velocity vs Time
%MATLAB
%25-05-2023
clc

%Inputs
m=300; %Mass
g=9.81; %Gravity
Theta=(37.5*pi)/180; %Slope      
Crr=0.007845; %Rolling Friction Coefficient
N=g*cos(Theta); %Normal Force without mass (Mass is divided)
BearingBore=2.7; %Bore in mm 
Cd=0.15; %Drag Coefficiemt
Fdownhill=g*sin(Theta); %Sloped Gravitational Force without mass (Mass is divided)
p=1.293; %Air Density
Area=1; %Cross sectional area
a=0; %Velocity is maximum when acceleration is 0

%Calculations
Vmax=sqrt(-1*(m*(a-Fdownhill+Crr*N+0.5*0.0015*g*2.7))/(0.5*p*Cd*Area)); %m/s
steps=linspace(0,Vmax,1000000);
Veloc= Vmax*3.6; %m/s to km/h
Velocity=linspace(0,Veloc,1000000);
Acceleration=[];
for i=1:length(steps)
Acceleration(end+1)=Fdownhill - Crr*N - 0.5*0.0015*g*BearingBore-(0.5*p*Cd*Area*steps(i)^2)/m;
end

plot(Velocity,Acceleration,'K--'); xlabel('Velocity (km/h)'); ylabel('Acceleration'); title('Velocity vs Acceleration')
uiwait

%Gianfranco Noriega Del Valle
%Maximum Velocity, Final Time, Final Acceleration, Maximum Displacement
%MATLAB
%29-05-2023

%Inputs

MaximumDisplacement = 466.38;

Time = 0;

Displacement = 0;

Counter = 0;

VelocityMetresPerSecond = Velocity * 3.6;

StoppingDistance = 75;

%Calculations

while Displacement < MaximumDisplacement

Counter = Counter + 1;

a = (Acceleration(Counter) + Acceleration(Counter + 1)) / 2; 
% ^^^ The change in acceleration from one step to the next is unbeliveably
% small (you can open up the 'Acceleration' variable in the workspace to
% see for yourself), so taking the average of the two is an accurate method
% of determing the change in acceleration during one step to the next

T = (VelocityMetresPerSecond(Counter + 1) - VelocityMetresPerSecond(Counter)) / a;

Time = Time + T;

X = (VelocityMetresPerSecond(Counter) * T) + (0.5 * a * (T^2));

Displacement = Displacement + X;

end

Time %#ok<NOPTS> 

FinalDisplacement = Displacement %#ok<NOPTS> 

FinalVelocity = VelocityMetresPerSecond(Counter + 1) %#ok<NOPTS> 

MinimumRequiredDeceleration = -FinalVelocity^2 / (2 * StoppingDistance) %#ok<NOPTS> 

if MinimumRequiredDeceleration > -50

    msgbox('BADA BING BADA BOOM BABY')

end