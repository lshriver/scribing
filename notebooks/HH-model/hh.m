% Integration of Hodgkin-Huxley equations with Euler method
% Source: Fundamentals in Computatinal Neuroscience by Thomas Trappenberg
clear; clc; close all

% maximal conductances (units of mS/cm^2); 1=K, 2=Na, 3=R
g(1)=36; g(2)=120; g(3)=0.3;

% battery voltage (mV); 1=n, 2=m, 3=h
E(1)=-12; E(2)=115; E(3)=10.613;

% initialize some of the variables
I_ext=0; V=-10; x=zeros(1,3); x(3)=1; t_rec=0;

% time steps for integration
dt = 0.01;

% integration with Euler method
for t=-30:dt:50
    if t==10; I_ext=10; end % turns external current on at t=10
    if t==40; I_ext=0; end  % turns external curren off at t=0

    % alpha functions used by Hodgkin-and-Huxely
    alpha(1)=(10-V)/(100*(exp((10-V)/10)-1));
    alpha(2)=(25-V)/(10*(exp((25-V)/10)-1));
    alpha(3)=0.07*exp(-V/20);

    % beta functiond used by Hodgkin and Huxely
    beta(1)=0.125*exp(-V/80);
    beta(2)=4*exp(-V/18);
    beta(3)=1/(exp((30-V)/10)+1);

    % Tau_x and x_0 are defined with respect to alpha and beta
    tau=1./(alpha+beta);
    x_0=alpha.*tau;

    % leaky integration with Euler's method
    x=(1-dt./tau).*x+dt./tau.*x_0;

    % calculate actual conductances g with g iven n, m, h
    gnmh(1)=g(1)*x(1)^4;
    gnmh(2)=g(2)*x(2)^3*x(3);
    gnmh(3)=g(3);
    
    % Ohm's law
    I=gnmh.*(V-E);

    % update voltage of membrane
    V=V+dt*(I_ext-sum(I));

    % record some variables for plotting after equilibrium
    if t>=0;
        t_rec=t_rec+1;
        x_plot(t_rec)=t;
        y_plot(t_rec)=V;
    end
end
plot(x_plot,y_plot); xlabel('Time'); ylabel('Voltage')
title('Hodgkin-Huxely Model')
grid on
