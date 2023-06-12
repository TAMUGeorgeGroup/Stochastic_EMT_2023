%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@author annice
%description: In the following program a stochastic model is used to
%simulate the EMT transition, good match is observed between simulations
%and theory
%Date: 02/22
%Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 100;
M_E_init = M;
M_H_init = 0;
M_M_init = 0;

mu_vect = [0.01, 0.01, 0.1, 0.1, 0.5, 0.5];
%mu_M_vect = [0.01, 0.01, 0.1, 0.1, 0.5, 0.5];
lambda_M_vect = [0.01, 1, 0.01, 1, 0.01, 1];
lambda_E_vect = [1, 0.01, 1, 0.01, 1, 0.01];

for i=1:6
    mu = mu_vect(i);
    lambda_E = lambda_E_vect(i);
    lambda_M = lambda_M_vect(i);
    [t,M_E, M_H, M_M] =  gillespie_simulate(M_E_init, M_M_init, M_H_init,M, lambda_M, lambda_E, mu);

    figure(9)
    subplot(3, 3, i)
    marker='.';
    plot(t, M_E, 'Marker', marker, 'MarkerFaceColor',[0.10198 0.4435 0.429752], 'MarkerEdgeColor',[0.10198 0.4435 0.429752],LineWidth=0.5)
    hold on
    plot(t, M_H, 'Marker', marker, 'MarkerFaceColor',[0.18656 0.30597 0.57462], 'MarkerEdgeColor',[0.18656 0.30597 0.57462],LineWidth=0.5)
    plot(t, M_M, 'Marker', marker, 'MarkerFaceColor',[0.4666 0.20645 0.32688], 'MarkerEdgeColor',[0.4666 0.20645 0.32688],LineWidth=0.5)
    xlabel("Time")
    ylabel("Number of Cells")


    p0 = [1 0 0];
    E=0;
    H=0;
    M=0;
    E(1) = 1;
    H(1) = 0;
    M(1) = 0;
    j= 2;

    %Please feel free to change these parameters
    if i==1 || i==2
        numnum = 15;
    elseif i==3 || i==4
        numnum = 1.5;
    else
        numnum = 0.45;
    end

    for k = [2:1:100]
    p = StochasticEMTMatrix_func(lambda_E, lambda_M, mu, numnum*(k-1));
    vect = p0 * p;
    E(j) = vect(1);
    H(j) = vect(2);
    M(j) = vect(3);
    j=j+1;
    p0 = vect;
    end
    marker = "*";
    time_line = linspace(0, t(length(t)), 100);

    time_line_x_E = [time_line(1:2:100)];
    E_y = [E(1:2:100)]*100;
    time_line_x_H = [time_line(1:2:100)];
    H_y = [H(1:2:100)]*100;
    time_line_x_M = [time_line(1:2:100)];
    M_y = [M(1:2:100)]*100;
    

    
    plot(time_line_x_E, E_y, 'Marker', marker, 'MarkerFaceColor',[0.51 0.32 0.71], 'MarkerEdgeColor',[0.51 0.32 0.71], LineWidth=3)
    hold on;
    plot(time_line_x_H, H_y, 'Marker', marker, 'MarkerFaceColor',[1 0.48 0.85], 'MarkerEdgeColor',[1 0.48 0.85],LineWidth=3)
    plot(time_line_x_M, M_y, 'Marker', marker, 'MarkerFaceColor',[0.28 0.62 0.56], 'MarkerEdgeColor',[0.28 0.62 0.56],LineWidth=3)
    legend("Epithelial Simulations", "Hybrid Simulations", "Mesenchymal Simulations", "Epithelial Theory", "Hybrid Theory", "Mesenchymal Theory")



end

function [t, M_E, M_H, M_M] = gillespie_simulate(M_E_init, M_M_init, M_H_init,M, lambda_M, lambda_E, mu)
%{
params:
M_E
M_M
M_H
lambda_M
lambda_E
mu
----
returns [t, M_E, M_H, M_M] : time and population size at three EMT related
phenotypic states
%}
    %if (M_E_init + M_H_init + M_M_init) ~=M
    %    disp("Warning, does not equate to M")
    %end
    

    t(1)=0; %starting at time 0
    i = 1; %while loop counter

    M_E(1) = M_E_init;
    M_H(1) = M_H_init;
    M_M(1) = M_M_init;

    while i<10000
            
        rate_E_H = mu;
        rate_H_E = lambda_E;
        rate_H_M = lambda_M;
        rate_M_H = mu;

        rand_num_1 = rand([M_E(i) 1]); %random number from a uniform dist between 0 and 1
        rand_num_2 = rand([M_H(i) 1]); %random number from a uniform dist between 0 and 1
        rand_num_3 = rand([M_H(i) 1]); %random number from a uniform dist between 0 and 1
        rand_num_4 = rand([M_M(i) 1]); %random number from a uniform dist between 0 and 1

        tau1 = (1/rate_E_H).*log(1./rand_num_1);
        tau2 = (1/rate_H_M).*log(1./rand_num_2);
        tau3 = (1/rate_H_E).*log(1./rand_num_3);
        tau4 = (1/rate_M_H).*log(1./rand_num_4);
        tau = min([tau1; tau2; tau3; tau4]);

        t(i+1) = t(i) + tau;

        if sum(tau == tau1)==1 && M_E(i)~=0
           M_E(i+1) = M_E(i) - 1;
           M_H(i+1) = M_H(i) + 1;
           M_M(i+1) = M_M(i);
        elseif sum(tau == tau2)==1 && M_H(i)~=0
           M_E(i+1) = M_E(i);
           M_H(i+1) = M_H(i) - 1;
           M_M(i+1) = M_M(i) + 1;
        elseif sum(tau == tau3)==1 && M_H(i)~=0
           M_E(i+1) = M_E(i) + 1;
           M_H(i+1) = M_H(i) - 1;
           M_M(i+1) = M_M(i);
        elseif sum(tau == tau4)==1 && M_M(i)~=0
           M_E(i+1) = M_E(i);
           M_H(i+1) = M_H(i) + 1;
           M_M(i+1) = M_M(i) - 1;
        else
           M_E(i+1) = M_E(i);
           M_H(i+1) = M_H(i);
           M_M(i+1) = M_M(i);
        end
       

        i = i + 1; %add to counter
   
    %track whether it acquires phenotype prior to extinction
    end

end
function P = StochasticEMTMatrix_func(lambdaE,lambdaM,mu,t)
    %Description: Returns the continuous-time stochastic matrix for
    %dynamics of a symmetric transitions between E,H,M states.
    %INPUT:
    %   lambdaE : rate of H->E
    %   lambdaM : rate of H->M
    %   mu      : rate of E->H and M->H
    %   t       : time.

    D = lambdaE+lambdaM+mu;
    piE = lambdaE/D;
    piM = lambdaM/D;
    piH = mu/D;
    
    k1=mu;
    k2=lambdaE+lambdaM+mu;
    k3=lambdaE+lambdaM;
    
    P =[piE+exp(-k1*t).*(1-lambdaE/k3*(1-piH*exp(-k3*t))) ...
        piH*(1-exp(-k2*t)) ...
        piM-exp(-k1*t)*(piM+lambdaM/k3*piH*(1-exp(-k3*t))); ...
        ...
        piE*(1-exp(-k2*t)) ...
        piH+exp(-k2*t)*(1-piH) ...
        piM*(1-exp(-k2*t)); ...
        ...
        piE-exp(-k1*t)*(piE+lambdaE/k3*piH*(1-exp(-k3*t))) ...
        piH*(1-exp(-k2*t)) ...
        piM+exp(-k1*t)*(1-lambdaM/k3*(1-piH*exp(-k3*t)))];
end
