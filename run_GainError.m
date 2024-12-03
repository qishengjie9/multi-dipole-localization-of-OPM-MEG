close all
clear 
%load FIXED_MEG_DATA
G1 = load('./template_array/SPMgainmatrix_Temp_dualsim_opm_30mm_1.mat');
G1 = G1.G;
load pos_40_snr_0sim_opm_30mm.mat
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
corrs = 0.5; % inter-sources correlation
T = 500; % time samples

M=cell(MC_repetitions,1);
% choose sources combinations
si = randi(ndipoles,nsources,MC_repetitions);
AP_max_iters=6; % AP maximal number of iterations
SNR_range = [-10 0 10 20];
%SNR_range = -10:10; % dB
gain_range = [0 0.02 0.04 0.06 0.08 0.1];
Nerror = length(gain_range);
gainerr_500result_2S_corr1 = zeros(2,3,length(SNR_range),MC_repetitions,Nerror);
% 2 stands for two evaluation metrics; 3 for the three inverse solutions
%ap_gain_result{1,1} =  aTdist_AP;

for i=1:Nerror
    for SNRindex = 1:length(SNR_range) % SNR loop
        SNR_range(SNRindex);
        for JJ = 1:MC_repetitions
            i
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
            Msnr = gain_error(Msnr,gain_range(JJ));
            
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

            gainerr_500result_2S_corr1(1,1,SNRindex,i,JJ) = mean(d_AP);
            gainerr_500result_2S_corr1(2,1,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP,index_AP);
            gainerr_500result_2S_corr1(1,2,SNRindex,i,JJ) = mean(d_ISW);
            gainerr_500result_2S_corr1(2,2,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP_W_MUSIC,index_ISW);
            gainerr_500result_2S_corr1(1,3,SNRindex,i,JJ) = mean(d_APM);
            gainerr_500result_2S_corr1(2,3,SNRindex,i,JJ) = calculate_RMSE(S,Gain,Msnr,S_AP_MUSIC,index_APM);
            
        end
    end
end
save('gainerr_500result_5S.mat','gainerr_500result_2S_corr1');

