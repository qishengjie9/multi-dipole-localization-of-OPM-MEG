close all
clear 
%load FIXED_MEG_DATA
G1 = load('./template_array/SPMgainmatrix_Temp_dualsim_opm_30mm_1.mat');
G1 = G1.G;
load ./template_array/Temp_dualsim_opm_30mm.mat
Cortex.Vertices = D.other.inv{1, 1}.mesh.tess_mni.vert;
Cortex.Faces = D.other.inv{1, 1}.mesh.tess_mni.face;
N_sensors = size(G1,1);
ndipoles  = size(Cortex.Vertices,1);
signalPeak = 500;
sensors = D.sensors.meg.chanpos;  
%size(sensors,1)


Gain = G1;
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
MC_repetitions = 500; % Monte-Carlo repetitions
nsources = 2; % number of sources
corrs = 1; % inter-sources correlation
T = 500; % time samples
dist_AP = zeros(MC_repetitions,1);
RMSE_AP = zeros(MC_repetitions,1);
dist_dISW = zeros(MC_repetitions,1);
RMSE_dISW = zeros(MC_repetitions,1);
dist_APMU = zeros(MC_repetitions,1);
RMSE_APMU = zeros(MC_repetitions,1);
M=cell(MC_repetitions,1);
% choose sources combinations
si = randi(ndipoles,nsources,MC_repetitions);
AP_max_iters=6; % AP maximal number of iterations
SNR_range = [-10 0 10 20];
%SNR_range = -10:10; % dB
crosstalk_range = [0 0.01 0.02 0.03 0.04 0.05];
Nerror = length(crosstalk_range);
crosstalkerr_dual_2S_corr1 = zeros(2,3,length(SNR_range),MC_repetitions,Nerror);
% 2 stands for two evaluation metrics; 3 for the three inverse solutions
%ap_gain_result{1,1} =  aTdist_AP;

for i=1:MC_repetitions
    i
    for SNRindex = 1:length(SNR_range) % SNR loop
        SNR_range(SNRindex);
        for JJ = 1:Nerror          
            d_AP=zeros(nsources,1);
            d_APM=zeros(nsources,1);
            d_ISW=zeros(nsources,1);
            dist_ap=zeros(nsources,1);
            dist_isw=zeros(nsources,1);
            dist_APM = zeros(nsources,1);
            S = gen_correlated_sources(corrs,T,nsources);
            M = Gain(:,si(:,i)') * S;
            scale = signalPeak/max(abs(M(:)));
            Ms = M*scale;
            MEG_energy = trace(Ms*Ms')/(N_sensors*(T));
            noise_var = MEG_energy/(10^(SNR_range(SNRindex)/10));
            Noise = randn(N_sensors,T).*sqrt(noise_var);
            Msnr = Ms + Noise;
            Msnr = crosstalk_dual_error(Msnr,sensors,crosstalk_range(JJ));
            
            % AP 
            [~,S_AP] = alternating_projections(Msnr, ndipoles, Gain, nsources,AP_max_iters,'AP');
            % AP-MUSIC 
            [~,S_AP_MUSIC] = alternating_projections(Msnr, ndipoles, Gain, nsources,AP_max_iters,'AP-MUSIC');
%           % AP-wMUSIC 
            [~,S_AP_W_MUSIC] = alternating_projections(Msnr, ndipoles, Gain, nsources,AP_max_iters,'AP-w-MUSIC');
            curr_pos=si(:,i);
            for s=1:nsources
                for n=1:nsources
                    dist_ap(n)       = norm(diff(Cortex.Vertices([curr_pos(s) S_AP(n)],:))); %distance in mm
                    dist_isw(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_W_MUSIC(n)],:))); %distance in mm
                    dist_APM(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_MUSIC(n)],:))); %distance in mm
                end
                [d_AP(s), index_AP(s)] = min(dist_ap);             
                [d_ISW(s), index_ISW(s)] = min(dist_isw);
                [d_APM(s), index_APM(s)] = min(dist_APM);
            end

            crosstalkerr_dual_2S_corr1(1,1,SNRindex,i,JJ) = mean(d_AP);
            crosstalkerr_dual_2S_corr1(2,1,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP,index_AP);
            crosstalkerr_dual_2S_corr1(1,2,SNRindex,i,JJ) = mean(d_ISW);
            crosstalkerr_dual_2S_corr1(2,2,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP_W_MUSIC,index_ISW);
            crosstalkerr_dual_2S_corr1(1,3,SNRindex,i,JJ) = mean(d_APM);
            crosstalkerr_dual_2S_corr1(2,3,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP_MUSIC,index_APM);
            
        end
    end
end
save('crosstalkerr_dual_2S_corr1.mat','crosstalkerr_dual_2S_corr1');

close all
clear 
%load FIXED_MEG_DATA
G1 = load('./template_array/SPMgainmatrix_Temp_dualsim_opm_30mm_1.mat');
G1 = G1.G;
load ./template_array/Temp_dualsim_opm_30mm.mat
Cortex.Vertices = D.other.inv{1, 1}.mesh.tess_mni.vert;
Cortex.Faces = D.other.inv{1, 1}.mesh.tess_mni.face;
N_sensors = size(G1,1);
ndipoles  = size(Cortex.Vertices,1);
signalPeak = 500;
sensors = D.sensors.meg.chanpos;  
%size(sensors,1)


Gain = G1;
%-----------------------------------------------------------------------
%-----------------------------------------------------------------------
MC_repetitions = 500; % Monte-Carlo repetitions
nsources = 5; % number of sources
corrs = 1; % inter-sources correlation
T = 500; % time samples
dist_AP = zeros(MC_repetitions,1);
RMSE_AP = zeros(MC_repetitions,1);
dist_dISW = zeros(MC_repetitions,1);
RMSE_dISW = zeros(MC_repetitions,1);
dist_APMU = zeros(MC_repetitions,1);
RMSE_APMU = zeros(MC_repetitions,1);
M=cell(MC_repetitions,1);
% choose sources combinations
si = randi(ndipoles,nsources,MC_repetitions);
AP_max_iters=6; % AP maximal number of iterations
SNR_range = [-10 0 10 20];
%SNR_range = -10:10; % dB
crosstalk_range = [0 0.01 0.02 0.03 0.04 0.05];
Nerror = length(crosstalk_range);
crosstalkerr_dual_5S_corr1 = zeros(2,3,length(SNR_range),MC_repetitions,Nerror);
% 2 stands for two evaluation metrics; 3 for the three inverse solutions
%ap_gain_result{1,1} =  aTdist_AP;

for i=1:MC_repetitions
    i
    for SNRindex = 1:length(SNR_range) % SNR loop
        SNR_range(SNRindex);
        for JJ = 1:Nerror          
            d_AP=zeros(nsources,1);
            d_APM=zeros(nsources,1);
            d_ISW=zeros(nsources,1);
            dist_ap=zeros(nsources,1);
            dist_isw=zeros(nsources,1);
            dist_APM = zeros(nsources,1);
            S = gen_correlated_sources(corrs,T,nsources);
            M = Gain(:,si(:,i)') * S;
            scale = signalPeak/max(abs(M(:)));
            Ms = M*scale;
            MEG_energy = trace(Ms*Ms')/(N_sensors*(T));
            noise_var = MEG_energy/(10^(SNR_range(SNRindex)/10));
            Noise = randn(N_sensors,T).*sqrt(noise_var);
            Msnr = Ms + Noise;
            Msnr = crosstalk_dual_error(Msnr,sensors,crosstalk_range(JJ));
            
            % AP 
            [~,S_AP] = alternating_projections(Msnr, ndipoles, Gain, nsources,AP_max_iters,'AP');
            % AP-MUSIC 
            [~,S_AP_MUSIC] = alternating_projections(Msnr, ndipoles, Gain, nsources,AP_max_iters,'AP-MUSIC');
%           % AP-wMUSIC 
            [~,S_AP_W_MUSIC] = alternating_projections(Msnr, ndipoles, Gain, nsources,AP_max_iters,'AP-w-MUSIC');
            curr_pos=si(:,i);
            for s=1:nsources
                for n=1:nsources
                    dist_ap(n)       = norm(diff(Cortex.Vertices([curr_pos(s) S_AP(n)],:))); %distance in mm
                    dist_isw(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_W_MUSIC(n)],:))); %distance in mm
                    dist_APM(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_MUSIC(n)],:))); %distance in mm
                end
                [d_AP(s), index_AP(s)] = min(dist_ap);             
                [d_ISW(s), index_ISW(s)] = min(dist_isw);
                [d_APM(s), index_APM(s)] = min(dist_APM);
            end

            crosstalkerr_dual_5S_corr1(1,1,SNRindex,i,JJ) = mean(d_AP);
            crosstalkerr_dual_5S_corr1(2,1,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP,index_AP);
            crosstalkerr_dual_5S_corr1(1,2,SNRindex,i,JJ) = mean(d_ISW);
            crosstalkerr_dual_5S_corr1(2,2,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP_W_MUSIC,index_ISW);
            crosstalkerr_dual_5S_corr1(1,3,SNRindex,i,JJ) = mean(d_APM);
            crosstalkerr_dual_5S_corr1(2,3,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP_MUSIC,index_APM);
            
        end
    end
end
save('crosstalkerr_dual_5S_corr1.mat','crosstalkerr_dual_5S_corr1');


