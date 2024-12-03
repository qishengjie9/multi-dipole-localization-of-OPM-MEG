close all
clear 
%load FIXED_MEG_DATA
G1 = load('SPMgainmatrix_0pos_1_snr_10sim_opm_30mm_pos_0_1.mat');
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
nsources = 2; % number of sources
corrs = 0; % inter-sources correlation
T = 10000; % time samples

M=cell(MC_repetitions,1);
% choose sources combinations
si = randi(ndipoles,nsources,MC_repetitions);
AP_max_iters=6; % AP maximal number of iterations
SNR_range = [20];
%SNR_range = -10:10; % dB

for i=1:6
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
        end
    end
end
save('gainerr_500result_5S.mat','gainerr_500result_2S_corr1');

