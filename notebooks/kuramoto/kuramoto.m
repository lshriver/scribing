close all
clc
clear

tic

set(0,'DefaultFigurePosition', [0 0 800 800])
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultTextFontSize',14);

colormap hsv
%whitebg([1 1 1])
%set(gca,'Color','k')
reset(gcf)
DOT_SIZE = 250;

%% Parameters
N = 500; %number of oscillators

%Create a beta distribution of intrinsic frequrencies for each oscillator 
    a1 = 2;  
    b1 = 2; 
    W = (betarnd(a1,b1,1,N));
%Create beta distribution of initial conditions for oscillator position
    a2 = 1;
    b2 = 1;
    P0 = (2*pi*betarnd(a2,b2,1,N)) - pi; 
    
K = 50; %coupling constant
Kn = K/N; %coupling strength
tmax = 50; 

%% Plot the PDFs for the intrinsic frequencies and the initial conditions
figure
reset(gcf)
whitebg([1 1 1])
subplot(2,1,1)
h_W = histogram(W);
grid on
h_W.FaceColor = [0, 0.4470, 0.7410];
xlabel('Intrinsic Frequency Distribution (g(\omega))');

subplot(2,1,2)
h_P0 = histogram(P0,[-pi:2*pi/20:pi]);
grid on
h_P0.FaceColor = [0.4940, 0.1840, 0.5560];
xlabel('Initial Conditions for \theta_i')

%% Define the phase frustration matrix
p_alpha = 0;
alpha_0 = 0;
alpha = (alpha_0)*pi;
alpha_matrix = rand(N,N); 
for i = 1:N
    alpha_matrix(i,i) = 0;
    for s = i + 1:N
        alpha_matrix(i,s) = alpha_matrix(i,s) < p_alpha;
        alpha_matrix(i,s) = alpha_matrix(s,i);
    end
end

alpha_matrix = alpha.*alpha_matrix;    

%Print out phase frustration parameters on workspace
T_alpha = table([alpha],[p_alpha],'VariableNames',{'alpha';'p_alpha'})

%% Create a random modular matrix
C = 10; %number of clusters/modules
p_edge = 0.05; %overall probability of attatchment
l = 0.75; %proportion of links within modules
R = random_modular_graph(N,C,p_edge,l);

A = full(R);
A_K = A.*K;
R = graph(R); %creates graph G(V,E) in MATLAB 

%Create tables for the random modular graph
    deg = degree(R);
    E = numedges(R);
    hub = (2*E)+(2*E*(0.1));
    d = [mean(mean(distances(R))), nnz(degree(R)>=hub),...
        E,isconnected(A),selfloops(A),rich_club_metric(A,(E+(E/20)))];
    T1 = table(d(:,1),d(:,3),d(:,4),d(:,5),d(:,6),...
        'VariableNames',{'AveragePathLength' ...
        'NumberofEdges' 'Connected' 'NumberofSelfLoops' 'RichClubMetric'})

%Plot the random modular graph
    figure
    colormap hsv
    whitebg([1 1 1])
    set(gcf,'inverthardcopy','off')
    deg = degree(R);
    nSizes = 2*sqrt(deg-min(deg)+10);
    nColors = deg;
    plot(R,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.4)
    title(['Random Modular Graph:',' ',...
        'N = ' num2str(N),', ',...
        'c = ' num2str(C),', ',...
        'p = ' num2str(p_edge),', ',...
        'r = ' num2str(l)]) 
    colorbar

%% Solve the system
%dPdt = @(P,W,A_K) W'+(1./N).*sum((A_K).*sin(meshgrid(P)-meshgrid(P)'),2);
dPdt = @(P,W,A_K) W'+(1./N).*sum((A_K).* sin(meshgrid(P)-(meshgrid(P)')-alpha_matrix),2);

[T,P] = ode45(@(t,p)dPdt(p,W,A_K),[0 tmax],P0);


%% Plot the results 
P2 = mod(P,2*pi); %rewrite P matrix entries in mod P
q = exp(1i*P); 
q_sum = sum(q,2)/N;
r_sum = abs(sum(q,2)/N); 

figure
whitebg([0 0 0])
rows = size(P,1);
for idx = 1:rows %length(P)
    
    subplot(2,2,1)
    colormap hsv 
    scatter(cos(P(idx,:)), sin(P(idx,:)), DOT_SIZE, 1:N,'.');
    grid on
    xlabel('cos(\theta_i)');
    ylabel('sin(\theta_i)');
    xlim([-1.5 1.5]); 
    ylim([-1.5 1.5]);
    
    subplot(2,2,2)
    scatter(T(idx)*ones(size(P(idx,:))), P(idx, :), DOT_SIZE, 1:N, '.'); 
    hold on
    xlabel('time(t)');
    ylabel('\theta_i(t)');
    
    subplot(2,2,3)
    t = 0:pi/50:2*pi;
    R = 1;
    x = R.*cos(t);
    y = R.*sin(t);
    plot(x,y,'color',[0.8500, 0.3250, 0.0980])
    grid on
    hold on
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
    plot(q_sum(idx),'o','color',[0, 0.4470, 0.7410]); 
    xlabel('$e^{i\theta_j}$', 'interpreter', 'latex')
    ylabel('$\sum{e^{i\theta_j}}/N$', 'interpreter', 'latex')
    
    subplot(2,2,4)
    plot(T(idx), r_sum(idx),'o','color',[0, 0.4470, 0.7410]);
    grid on
    hold on
    xlabel('time (t)');
    ylabel('$R = |\sum{e^{i\theta_j(t)}}/N|$', 'interpreter', 'latex')
    ylim([0 1.5]);
    
    drawnow; 
    pause(0.1); 
end
hold off
toc