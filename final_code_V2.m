clear;
clc;
tb = 1/360;
fs = 360;
choice = input("For different files of noise and signal press '1' otherwise press '2' : ");
if choice == 1
    %% Input taking ECG 
    path = uigetdir('','select the folder containing csv files:');
    if path == 0
        error('No folder selected !!');
    end
    files = dir(fullfile(path,'*.csv'));
    if isempty(files)
        error('No csv files found !');
    end
    file_names = {files.name};
    [idx, tf] = listdlg("PromptString",'Select a file:','SelectionMode','single','ListString',file_names);
    if tf ==0
        error('No file selected!');
    end
    filepath = fullfile(path,files(idx).name);
    data = readtable(filepath);
    headers = data.Properties.VariableNames;
    data_mat_col = table2array(data);
    
    %% ECG Signal calculations
    data_mat = data_mat_col'; 
    E1 = data_mat(1,:)'; % Transposed to column vector
    E2 = data_mat(2,:)'; % Transposed to column vector
    title1 = strrep(headers{1}, '_', ' '); 
    title2 = strrep(headers{2}, '_', ' ');
    %% Input taking noise
    path2 = uigetdir('','select the folder containing csv files:');
    if path2 == 0
        error('No folder selected !!');
    end
    files2 = dir(fullfile(path2,'*.csv'));
    if isempty(files2)
        error('No csv files found !');
    end
    file_names2 = {files2.name};
    [idx2, tf2] = listdlg("PromptString",'Select a file:','SelectionMode','single','ListString',file_names2);
    if tf2 ==0
        error('No file selected!');
    end
    filepath2 = fullfile(path2,files2(idx2).name);
    data2 = readtable(filepath2);
    headers2 = data2.Properties.VariableNames;
    data_mat_col2 = table2array(data2);
    
    %% noise Signal calculations
    data_mat2 = data_mat_col2'; 
    N1 =data_mat2(1,:)';
    N2= data_mat2(2,:)';
    
    %% Signal Corruption
    target_SNR_dB = 6; 
    N = length(E1);
    t = (0:N-1)'/fs;
    [x1,alpha_x1] = corrupt_channel(E1, N1, target_SNR_dB);
    [x2,alpha_x2] = corrupt_channel(E2, N1, target_SNR_dB);
     
    corrupted_ecg = [x1, x2];

elseif choice == 2
    %% Input taking ECG 
    path = uigetdir('','select the folder containing csv files:');
    if path == 0
        error('No folder selected !!');
    end
    files = dir(fullfile(path,'*.csv'));
    if isempty(files)
        error('No csv files found !');
    end
    file_names = {files.name};
    [idx, tf] = listdlg("PromptString",'Select a file:','SelectionMode','single','ListString',file_names);
    if tf ==0
        error('No file selected!');
    end
    filepath = fullfile(path,files(idx).name);
    data = readtable(filepath);
    headers = data.Properties.VariableNames;
    data_mat_col = table2array(data);
    
    %% ECG Signal calculations
    data_mat = data_mat_col'; 
    x1 = data_mat(1,:)'; % Transposed to column vector
    x2 = data_mat(2,:)'; % Transposed to column vector
    title1 = strrep(headers{1}, '_', ' '); 
    title2 = strrep(headers{2}, '_', ' ');

else
    disp("Invalid Choice !!!");

end

N = length(x1);
t = (0:N-1)'/fs;
disp_samples = min(600 * fs, N);

function [y, alpha] = corrupt_channel(s, n, target_snr)
    sac = s - mean(s);
    nac = n - mean(n);
    Ps = var(sac);
    Pn = var(nac);
    alpha = sqrt(Ps/(Pn*10^(target_snr/10)));
    y = s+alpha *nac;
end

%% 2nd Derivative using Central Finite Difference Method
dx1_2nd =zeros(size(x1));
dx2_2nd =zeros(size(x2));

dx1_2nd(2:end-1) = (x1(3:end) -2*x1(2:end-1)+x1(1:end-2))*(fs^2);
dx2_2nd(2:end-1) = (x2(3:end) - 2*x2(2:end-1)+x2(1:end-2))*(fs^2);

dx1_2nd(1) = (2*x1(1) -5*x1(2) + 4*x1(3) - x1(4))*(fs^2);
dx1_2nd(end) =(2*x1(end)- 5*x1(end-1) + 4*x1(end-2) -x1(end-3))*(fs^2);

dx2_2nd(1) = (2*x2(1)-5*x2(2) +4*x2(3) -x2(4))*(fs^2);
dx2_2nd(end) = (2*x2(end)-5*x2(end-1)+4*x2(end-2)-x2(end-3))*(fs^2);

d2_corrupted_ecg = [dx1_2nd,dx2_2nd];

%% Filter Initialization
fc =0.5;
wn = 2*pi*fc;
zeta =0.707;
dt =tb;

z1_x1 = 0; 
z2_x1 = 0; 
z1_x2 = 0; 
z2_x2 = 0;

y1_rk4 =zeros(N, 1);
y2_rk4 =zeros(N, 1);

f1 = @(z2) z2;
f2 = @(z1, z2, d2x) -2*zeta*wn*z2 - wn^2*z1 + d2x;

%% 4th order RK4 Implementation
for i = 1:N-1
    d2x1_i = d2_corrupted_ecg(i, 1);
    d2x2_i = d2_corrupted_ecg(i, 2);
    
    if i <N-1
        d2x1_mid = (d2_corrupted_ecg(i,1) +d2_corrupted_ecg(i+1,1))/2;
        d2x2_mid = (d2_corrupted_ecg(i,2) + d2_corrupted_ecg(i+1,2))/2;
        d2x1_next = d2_corrupted_ecg(i+1,1);
        d2x2_next = d2_corrupted_ecg(i+1,2);
    else
        d2x1_mid =d2x1_i; 
        d2x2_mid = d2x2_i;
        d2x1_next = d2x1_i; 
        d2x2_next = d2x2_i;
    end

    k1_1_1 = f1(z2_x1);
    k2_1_1 =f2(z1_x1, z2_x1, d2x1_i);
    k1_2_1 = f1(z2_x1 + (dt/2)*k2_1_1);
    k2_2_1 = f2(z1_x1 + (dt/2)*k1_1_1, z2_x1+(dt/2)*k2_1_1,d2x1_mid);   
    k1_3_1 = f1(z2_x1 + (dt/2)*k2_2_1);
    k2_3_1 = f2(z1_x1 + (dt/2)*k1_2_1, z2_x1+(dt/2)*k2_2_1,d2x1_mid);    
    k1_4_1 = f1(z2_x1 + dt*k2_3_1);
    k2_4_1 = f2(z1_x1 + dt*k1_3_1, z2_x1+dt*k2_3_1,d2x1_next);

    z1_x1 = z1_x1 + (dt/6)*(k1_1_1 + 2*k1_2_1 + 2*k1_3_1 + k1_4_1);
    z2_x1 = z2_x1 + (dt/6)*(k2_1_1 +2*k2_2_1 + 2*k2_3_1 + k2_4_1);
    y1_rk4(i+1) = z1_x1;

    k1_1_2 = f1(z2_x2);
    k2_1_2 = f2(z1_x2, z2_x2, d2x2_i);    
    k1_2_2 = f1(z2_x2 +(dt/2)*k2_1_2);
    k2_2_2 = f2(z1_x2 + (dt/2)*k1_1_2, z2_x2 +(dt/2)*k2_1_2, d2x2_mid);    
    k1_3_2 = f1(z2_x2 +(dt/2)*k2_2_2);
    k2_3_2 = f2(z1_x2 + (dt/2)*k1_2_2, z2_x2 +(dt/2)*k2_2_2, d2x2_mid);    
    k1_4_2 = f1(z2_x2 +dt*k2_3_2);
    k2_4_2 = f2(z1_x2 + dt*k1_3_2, z2_x2 + dt*k2_3_2, d2x2_next);
    z1_x2 = z1_x2 + (dt/6)*(k1_1_2 + 2*k1_2_2 + 2*k1_3_2 + k1_4_2);
    z2_x2 = z2_x2 + (dt/6)*(k2_1_2 + 2*k2_2_2 + 2*k2_3_2 +k2_4_2);
    y2_rk4(i+1) =z1_x2;
end

%% Parametric Analysis: Euler's Method for comparison on Signal 1
y1_euler = zeros(N, 1);
z1_e = 0; z2_e = 0;
for i = 1:N-1
    z1_next = z1_e + dt * f1(z2_e);
    z2_next = z2_e + dt * f2(z1_e, z2_e, d2_corrupted_ecg(i, 1));
    z1_e = z1_next;
    z2_e = z2_next;
    y1_euler(i+1) = z1_e;
end

%% Expected Outputs and Plots
% 1. Time-Domain Comparison Plot
figure('Name', 'Time-Domain Comparison', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t(1:disp_samples), x1(1:disp_samples), 'r', 'LineWidth', 0.8); hold on;
plot(t(1:disp_samples), y1_rk4(1:disp_samples), 'b', 'LineWidth', 1);
title([title1, ': Corrupted Input vs Cleaned Output (RK4)']);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend('Corrupted x(t)', 'Cleaned y(t)'); grid on;
xlim([0,5]);

subplot(2,1,2);
plot(t(1:disp_samples), x2(1:disp_samples), 'k', 'LineWidth', 0.8); hold on;
plot(t(1:disp_samples), y2_rk4(1:disp_samples), 'b', 'LineWidth', 1);
title([title2, ': Corrupted Input vs Cleaned Output (RK4)']);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend('Corrupted x(t)', 'Cleaned y(t)'); grid on;
xlim([0,5]);

% 2. Error and Residue Plot
noise_extracted1 = x1 - y1_rk4;
noise_extracted2 = x2 - y2_rk4;

figure('Name', 'Error and Residue Plot', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t(1:disp_samples), noise_extracted1(1:disp_samples), 'k--', 'LineWidth', 1);
title([title1, 'Extracted Baseline Wander']);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
grid on;
xlim([0,5]);

subplot(2,1,2);
plot(t(1:disp_samples), noise_extracted2(1:disp_samples), 'k--', 'LineWidth', 1);
title([title2, 'Extracted Baseline Wander']);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
grid on;
xlim([0,5]);

% 3. Parametric Analysis Plot
figure('Name', 'Parametric Analysis: RK4 vs Euler', 'NumberTitle', 'off');
plot(t(1:disp_samples), y1_rk4(1:disp_samples), 'b', 'LineWidth', 1.2); hold on;
plot(t(1:disp_samples), y1_euler(1:disp_samples), 'g--', 'LineWidth', 1);
title([title1, ': Output Stability (RK4 vs Euler Method)']);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend('RK4 Output', 'Euler Output'); grid on;
xlim([0,5]);