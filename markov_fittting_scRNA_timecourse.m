%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@author annice
%description: In the following program we extract the
%percentage of E, H and M phenotypes of time-dependent data and fit the
%CTMC model
%Date: 07/04/22
%Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Constants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_path = '/Users/annicenajafi/Desktop/Markov_fitting/'; %Path to context directory
transcription_factor = ["TGFB1", "TNF", "EGF"];
cell_line = ["A549", "DU145","MCF7","OVCA420"];
cutoffs = ["45", "40", "50", "55"];
treatment_withdrawal_day = 5; %Day 7, the fifth timepoint
treatment_withdrawal_day_actual_day = 7;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for u=1:length(cell_line)

%Read the file
file_path = '/Users/annicenajafi/Desktop/DTW_based_align/timecourse_data_' + cell_line(u)+ "_"+ transcription_factor(1)+"_"+cutoffs(u)+ ".csv";    
data = csvread(file_path, 1, 0);

%initial value for alpha
alph=1;

%fun = @(alpha_hold)find_min_alph_exp(data, alpha_hold);
%[alph, res] = fminsearch(fun,alpha_hold);

%Extract fractions from file
data_time = data(:, 1);
data_H = data(:, 3); %Extract hybrid
E_cad = data(:, 2);
hybrid = data_H; 
ZEB = data(:, 4);

[max_H, max_ind] = max(data_H(1:treatment_withdrawal_day)); %Find where the the maximum happens pre withdrawal


pi_E = E_cad(max_ind);
pi_H = hybrid(max_ind);
pi_M = ZEB(max_ind);

E(1) = E_cad(1);
H(1) = hybrid(1);
M(1) = ZEB(1);
numnum=0.04;
p0 = [E(1) H(1) M(1)];
t = 150;
min_E = inf;
min_H = inf;
min_M = inf;
min_total = inf;
hold_val_E = 0; %alpha factor stored here
hold_val_H = 0;
hold_val_M = 0;
hold_val = 0;
mse_hold(1)=0;
count = 1;
arbit_scale = 20/3; %arbitrarily chosen, feel free to change

%Define lambdaE, lambdaM, and Mu for first Markov chain
lambda_E = pi_E/pi_H;
mu = 1;
lambda_M = pi_M/pi_H;
if(E_cad(treatment_withdrawal_day)==0)
%Define scaling factors for rates for the second Markov chain
M_sc = (E_cad(treatment_withdrawal_day)+1e-2)/ZEB(treatment_withdrawal_day);
Mu_sc = (E_cad(treatment_withdrawal_day)+1e-2)/hybrid(treatment_withdrawal_day);
else
M_sc = (E_cad(treatment_withdrawal_day))/ZEB(treatment_withdrawal_day);
Mu_sc = (E_cad(treatment_withdrawal_day))/hybrid(treatment_withdrawal_day);
end


timepoints = data_time(max_ind:treatment_withdrawal_day)';
fun = @(alph)find_min_alph_exp(alph, E_cad, hybrid, ZEB, M_sc, Mu_sc, max_ind-1, treatment_withdrawal_day, timepoints);
[alph, res] = fminsearch(fun,alph);




for alpha = alph
   
E_st(1)=E(1);
H_st(1)=H(1);
M_st(1)=M(1);
    %Before reaching hybrid peak
    if(u==1)
        E_st=NaN;
        H_st=NaN;
        M_st=NaN;
        loop_to = arbit_scale*(treatment_withdrawal_day);
    else
    for jj = [1:1:arbit_scale*data_time(max_ind)] %arbitrarily normalized, feel free to change *20/3
            p = StochasticEMTMatrix_func(lambda_E, lambda_M, mu, numnum*(jj-1));
            vect = p0 * p;
            E_st(jj) = vect(1);
            H_st(jj) = vect(2);
            M_st(jj) = vect(3);
            p0 = vect;
            loop_to = arbit_scale *(treatment_withdrawal_day_actual_day-data_time(max_ind))/(data_time(max_ind));
    end
    end
    %After passing hybrid peak
    for k = [1:1:loop_to]
            lambda_E = alph;
            lambda_M = alph/M_sc;
            mu = alph/Mu_sc;
            p = StochasticEMTMatrix_func(lambda_E, lambda_M, mu, numnum*(k-1));
            vect = p0 * p;
            E(k) = vect(1);
            H(k) = vect(2);
            M(k) = vect(3);
            p0 = vect;
    end

    alpha_rev=1;
    fun = @(alph)find_min_alph_exp(alph, E_cad, hybrid, ZEB, (ZEB(1)/E_cad(1)), (hybrid(1)/E_cad(1)), 4, 8, [7 7.33 8 10]);
    [alph_rev, res_rev] = fminsearch(fun,alpha_rev);
    %After treatment withdrawal
    for j = [1:1:20]
            lambda_E = alph_rev;
            lambda_M = alph_rev*(ZEB(1)/E_cad(1));
            mu = alph_rev*(hybrid(1)/E_cad(1));
            p = StochasticEMTMatrix_func(lambda_E, lambda_M, mu, numnum*(j-1));
            vect = p0 * p;
            E_end(j) = vect(1);
            H_end(j) = vect(2);
            M_end(j) = vect(3);
            p0 = vect;
    end
    marker = "_";
    time_line = linspace(0, t(length(t)), 100);
    time_line_x_E = [linspace(0, data_time(max_ind), arbit_scale*data_time(max_ind)) linspace(data_time(max_ind), data_time(treatment_withdrawal_day), loop_to) linspace(data_time(treatment_withdrawal_day), data_time(end), 20)];
    E_y = [E_st E E_end];
    time_line_x_H = [linspace(0, data_time(max_ind), arbit_scale*data_time(max_ind)) linspace(data_time(max_ind), data_time(treatment_withdrawal_day), loop_to) linspace(data_time(treatment_withdrawal_day), data_time(end), 20)];
    H_y = [H_st H H_end];
    time_line_x_M = [linspace(0, data_time(max_ind), arbit_scale*data_time(max_ind)) linspace(data_time(max_ind), data_time(treatment_withdrawal_day), loop_to) linspace(data_time(treatment_withdrawal_day), data_time(end), 20)];
    M_y = [M_st M M_end];

    if(u==1)
        E_y = E_y(2:end);
        H_y = H_y(2:end);
        M_y = M_y(2:end);
        max_ind = 5;
    end
    mse_E = mean((E_y(1:length(E_y)/treatment_withdrawal_day:length(E_y))'-E_cad(1:treatment_withdrawal_day)).^2)/treatment_withdrawal_day;
    mse_H = mean((H_y(1:length(H_y)/treatment_withdrawal_day:length(H_y))'-hybrid(1:treatment_withdrawal_day)).^2)/treatment_withdrawal_day;
    mse_M = mean((M_y(1:length(M_y)/treatment_withdrawal_day:length(M_y))'-ZEB(1:treatment_withdrawal_day)).^2)/treatment_withdrawal_day;
    mse_total = mse_E + mse_M + mse_H;
    mse_hold_E(count) = mse_E;
    mse_hold_M(count) = mse_M;
    mse_hold_H(count) = mse_H;
    mse_hold(count) = mse_total;
    count = count+1;

    disp(mse_total)
    if mse_E<min_E
        min_E = mse_E;
        hold_val_E = alpha;
    end
    if mse_H<min_H
        min_H = mse_H;
        hold_val_H = alpha;
    end
    if mse_M<min_M
        min_M = mse_M;
        hold_val_M = alpha;
    end
    if mse_total<min_total
        min_total = mse_total;
        hold_val = alpha;
    end
end
   
    figure(31); 
    subplot(2, 4, u)
    xmin=0;
    xmax=7;
    ymin=0;
    ymax=1;
    patch([xmin xmax xmax xmin],[ymin ymin ymax ymax], [0.84,0.91,0.94])
    xline(data_time(max_ind),'k--', LineWidth=2)
     hold on;
    plot(time_line_x_E , E_y, ':', 'Color', [0.58,0.36,0.63], LineWidth=3)
    plot(time_line_x_H , H_y, ':', 'Color', [0.93,0.74,0.96], LineWidth=3)
    plot(time_line_x_M , M_y, ':', 'Color', [0.43,0.73,0.67], LineWidth=3)
    plot(data_time, E_cad, '-*', 'Color', [0.58,0.36,0.63],  LineWidth=3)
    plot(data_time, hybrid, '-*', 'Color', [0.93,0.74,0.96],  LineWidth=3)
    plot(data_time, ZEB, '-*', 'Color', [0.43,0.73,0.67],  LineWidth=3)

    xlim([0 10])
    title("Experimental EMT Phenotypic Ratio")
    ylabel("Fraction")
    xlabel("Time (Days)")
    legend("","Regime Change", "Data Epithelial", "Data Hybrid", "Data Mesenchymal", "Epithelial Theoretic Fraction", "Hybrid Theoretic Fraction", "Mesenchymal Theoretic Fraction")
 
    clearvars -except dir_path transcription_factor cell_line cutoffs treatment_withdrawal_day treatment_withdrawal_day_actual_day
end
