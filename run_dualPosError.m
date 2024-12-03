clc
clear all
%load FIXED_MEG_DATA
G1 = load('template_array/SPMgainmatrix_Temp_dualsim_opm_30mm_1.mat');
G1 = G1.G;
G2 = load('template_array/SPMgainmatrix_Temp_pos1Temp_dualsim_opm_30mm_1.mat');
G2 = G2.G;
G3 = load('template_array/SPMgainmatrix_Temp_pos2Temp_dualsim_opm_30mm_1.mat');
G3 = G3.G;
G4 = load('template_array/SPMgainmatrix_Temp_pos3Temp_dualsim_opm_30mm_1.mat');
G4 = G4.G;
G5 = load('template_array/SPMgainmatrix_Temp_pos4Temp_dualsim_opm_30mm_1.mat');
G5 = G5.G;
G6 = load('template_array/SPMgainmatrix_Temp_pos5Temp_dualsim_opm_30mm_1.mat');
G6 = G6.G;
load template_array/Temp_dualsim_opm_30mm.mat
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

M=cell(MC_repetitions,1);
% choose sources combinations
si = randi(ndipoles,nsources,MC_repetitions);
AP_max_iters=6; % AP maximal number of iterations
SNR_range = [-10 0 10 20];
%SNR_range = -10:10; % dB
pos_range = [0 1 2 3 4 5];
Nerror_pos = length(pos_range);
poserr_dual_2S_corr1 = zeros(2,length(SNR_range),MC_repetitions,Nerror_pos);
% 2 stands for two evaluation metrics; 3 for the three inverse solutions
%ap_gain_result{1,1} =  aTdist_AP;

for i=1:MC_repetitions
    i
    for SNRindex = 1:length(SNR_range) % SNR loop
        SNR_range(SNRindex);
        for JJ = 1:Nerror_pos
            Gain_err = sprintf('G%d',JJ);
            Gain_err = eval(Gain_err,'=Gain_err');
            d_AP=zeros(nsources,1);
            dist_ap=zeros(nsources,1);
            S = gen_correlated_sources(corrs,T,nsources);
            M = Gain(:,si(:,i)') * S;
            scale = signalPeak/max(abs(M(:)));
            Ms = M*scale;
            MEG_energy = trace(Ms*Ms')/(N_sensors*(T));
            noise_var = MEG_energy/(10^(SNR_range(SNRindex)/10));
            Noise = randn(N_sensors,T).*sqrt(noise_var);
            Msnr = Ms + Noise;
            %Msnr_crosstalkerr = crosstalk_err(Msnr,sensors,crosstalk_range(JJ));
            
            % AP 
            [~,S_AP] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP');
            % AP-MUSIC 
%             [~,S_AP_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-MUSIC');
% %           % AP-wMUSIC 
%             [~,S_AP_W_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-w-MUSIC');
            curr_pos=si(:,i);
            for s=1:nsources
                for n=1:nsources
                    dist_ap(n)       = norm(diff(Cortex.Vertices([curr_pos(s) S_AP(n)],:))); %distance in mm
%                     dist_isw(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_W_MUSIC(n)],:))); %distance in mm
%                     dist_APM(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_MUSIC(n)],:))); %distance in mm
                end
                [d_AP(s),  index_AP(s)] = min(dist_ap);             
%                 [d_ISW(s), index_ISW(s)] = min(dist_isw);
%                 [d_APM(s), index_APM(s)] = min(dist_APM);
            end

            poserr_dual_2S_corr1(1,SNRindex,i,JJ) = mean(d_AP);
            poserr_dual_2S_corr1(2,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP,index_AP);
%             poserr_500result_5S(1,2,SNRindex,i,JJ) = mean(d_ISW);
%             poserr_500result_5S(2,2,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_W_MUSIC,index_ISW);
%             poserr_500result_5S(1,3,SNRindex,i,JJ) = mean(d_APM);
%             poserr_500result_5S(2,3,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_MUSIC,index_APM);
            
        end

    end
end
save('poserr_dual_2S_corr1.mat','poserr_dual_2S_corr1');

clc
close all
clear all
%load FIXED_MEG_DATA
G1 = load('template_array/SPMgainmatrix_Temp_dualsim_opm_30mm_1.mat');
G1 = G1.G;
G2 = load('template_array/SPMgainmatrix_Temp_pos1Temp_dualsim_opm_30mm_1.mat');
G2 = G2.G;
G3 = load('template_array/SPMgainmatrix_Temp_pos2Temp_dualsim_opm_30mm_1.mat');
G3 = G3.G;
G4 = load('template_array/SPMgainmatrix_Temp_pos3Temp_dualsim_opm_30mm_1.mat');
G4 = G4.G;
G5 = load('template_array/SPMgainmatrix_Temp_pos4Temp_dualsim_opm_30mm_1.mat');
G5 = G5.G;
G6 = load('template_array/SPMgainmatrix_Temp_pos5Temp_dualsim_opm_30mm_1.mat');
G6 = G6.G;
load template_array/Temp_dualsim_opm_30mm.mat
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

M=cell(MC_repetitions,1);
% choose sources combinations
si = randi(ndipoles,nsources,MC_repetitions);
AP_max_iters=6; % AP maximal number of iterations
SNR_range = [-10 0 10 20];
%SNR_range = -10:10; % dB
pos_range = [0 1 2 3 4 5];
Nerror_pos = length(pos_range);
poserr_dual_5S_corr1 = zeros(2,length(SNR_range),MC_repetitions,Nerror_pos);
% 2 stands for two evaluation metrics; 3 for the three inverse solutions
%ap_gain_result{1,1} =  aTdist_AP;

for i=1:MC_repetitions
    i
    for SNRindex = 1:length(SNR_range) % SNR loop
        SNR_range(SNRindex);
        for JJ = 1:Nerror_pos
            Gain_err = sprintf('G%d',JJ);
            Gain_err = eval(Gain_err,'=Gain_err');
            d_AP=zeros(nsources,1);
            dist_ap=zeros(nsources,1);
            S = gen_correlated_sources(corrs,T,nsources);
            M = Gain(:,si(:,i)') * S;
            scale = signalPeak/max(abs(M(:)));
            Ms = M*scale;
            MEG_energy = trace(Ms*Ms')/(N_sensors*(T));
            noise_var = MEG_energy/(10^(SNR_range(SNRindex)/10));
            Noise = randn(N_sensors,T).*sqrt(noise_var);
            Msnr = Ms + Noise;
            %Msnr_crosstalkerr = crosstalk_err(Msnr,sensors,crosstalk_range(JJ));
            
            % AP 
            [~,S_AP] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP');
            % AP-MUSIC 
%             [~,S_AP_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-MUSIC');
% %           % AP-wMUSIC 
%             [~,S_AP_W_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-w-MUSIC');
            curr_pos=si(:,i);
            for s=1:nsources
                for n=1:nsources
                    dist_ap(n)       = norm(diff(Cortex.Vertices([curr_pos(s) S_AP(n)],:))); %distance in mm
%                     dist_isw(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_W_MUSIC(n)],:))); %distance in mm
%                     dist_APM(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_MUSIC(n)],:))); %distance in mm
                end
                [d_AP(s),  index_AP(s)] = min(dist_ap);             
%                 [d_ISW(s), index_ISW(s)] = min(dist_isw);
%                 [d_APM(s), index_APM(s)] = min(dist_APM);
            end

            poserr_dual_5S_corr1(1,SNRindex,i,JJ) = mean(d_AP);
            poserr_dual_5S_corr1(2,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP,index_AP);
%             poserr_500result_5S(1,2,SNRindex,i,JJ) = mean(d_ISW);
%             poserr_500result_5S(2,2,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_W_MUSIC,index_ISW);
%             poserr_500result_5S(1,3,SNRindex,i,JJ) = mean(d_APM);
%             poserr_500result_5S(2,3,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_MUSIC,index_APM);
            
        end

    end
end
save('poserr_dual_5S_corr1.mat','poserr_dual_5S_corr1');

