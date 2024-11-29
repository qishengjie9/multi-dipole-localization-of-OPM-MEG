close all
clear 
%load FIXED_MEG_DATA
G1 = load('SPMgainmatrix_pos_1_snr_10sim_opm_10mm_angular_1_1.mat');
G1 = G1.G;
G2 = load('SPMgainmatrix_2pos_1_snr_0sim_opm_10mm_angular_2_1.mat');
G2 = G2.G;
G3 = load('SPMgainmatrix_3pos_1_snr_0sim_opm_10mm_angular_3_1.mat');
G3 = G3.G;
G4 = load('SPMgainmatrix_4pos_1_snr_0sim_opm_10mm_angular_4_1.mat');
G4 = G4.G;
G5 = load('SPMgainmatrix_5pos_1_snr_0sim_opm_10mm_angular_5_1.mat');
G5 = G5.G;
G6 = load('SPMgainmatrix_6pos_1_snr_0sim_opm_10mm_angular_6_1.mat');
G6 = G6.G;
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
angular_range = [0 1 2 3 4 5];
Nerror_angular = length(angular_range);
angularerr_500result_5S = zeros(2,3,length(SNR_range),MC_repetitions,Nerror_angular);
% 2 stands for two evaluation metrics; 3 for the three inverse solutions
%ap_gain_result{1,1} =  aTdist_AP;

for i=1:MC_repetitions
    i
    for SNRindex = 1:length(SNR_range) % SNR loop
        SNR_range(SNRindex);
        for JJ = 1:Nerror_angular
            Gain_err = sprintf('G%d',JJ);
            Gain_err = eval(Gain_err,'=Gain_err');
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
            %Msnr_crosstalkerr = crosstalk_err(Msnr,sensors,crosstalk_range(JJ));
            
            % AP 
            [~,S_AP] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP');
            % AP-MUSIC 
            [~,S_AP_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-MUSIC');
%           % AP-wMUSIC 
            [~,S_AP_W_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-w-MUSIC');
            curr_pos=si(:,i);
            for s=1:nsources
                for n=1:nsources
                    dist_ap(n)       = norm(diff(Cortex.Vertices([curr_pos(s) S_AP(n)],:))); %distance in mm
                    dist_isw(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_W_MUSIC(n)],:))); %distance in mm
                    dist_APM(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_MUSIC(n)],:))); %distance in mm
                end
                [d_AP(s),  index_AP(s)] = min(dist_ap);             
                [d_ISW(s), index_ISW(s)] = min(dist_isw);
                [d_APM(s), index_APM(s)] = min(dist_APM);
            end

            angularerr_500result_5S(1,1,SNRindex,i,JJ) = mean(d_AP);
            angularerr_500result_5S(2,1,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP,index_AP);
            angularerr_500result_5S(1,2,SNRindex,i,JJ) = mean(d_ISW);
            angularerr_500result_5S(2,2,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_W_MUSIC,index_ISW);
            angularerr_500result_5S(1,3,SNRindex,i,JJ) = mean(d_APM);
            angularerr_500result_5S(2,3,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_MUSIC,index_APM);
            
        end

    end
end
save('angularerr_500result_5S.mat','angularerr_500result_5S');

close all
clear 
%load FIXED_MEG_DATA
G1 = load('SPMgainmatrix_pos_1_snr_10sim_opm_10mm_angular_1_1.mat');
G1 = G1.G;
G2 = load('SPMgainmatrix_2pos_1_snr_0sim_opm_10mm_angular_2_1.mat');
G2 = G2.G;
G3 = load('SPMgainmatrix_3pos_1_snr_0sim_opm_10mm_angular_3_1.mat');
G3 = G3.G;
G4 = load('SPMgainmatrix_4pos_1_snr_0sim_opm_10mm_angular_4_1.mat');
G4 = G4.G;
G5 = load('SPMgainmatrix_5pos_1_snr_0sim_opm_10mm_angular_5_1.mat');
G5 = G5.G;
G6 = load('SPMgainmatrix_6pos_1_snr_0sim_opm_10mm_angular_6_1.mat');
G6 = G6.G;
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
corrs = 1; % inter-sources correlation
T = 500; % time samples

M=cell(MC_repetitions,1);
% choose sources combinations
si = randi(ndipoles,nsources,MC_repetitions);
AP_max_iters=6; % AP maximal number of iterations
SNR_range = [-10 0 10 20];
%SNR_range = -10:10; % dB
angular_range = [0 1 2 3 4 5];
Nerror_angular = length(angular_range);
angularerr_500result_5S_corr1 = zeros(2,3,length(SNR_range),MC_repetitions,Nerror_angular);
% 2 stands for two evaluation metrics; 3 for the three inverse solutions
%ap_gain_result{1,1} =  aTdist_AP;

for i=1:MC_repetitions
    i
    for SNRindex = 1:length(SNR_range) % SNR loop
        SNR_range(SNRindex);
        for JJ = 1:Nerror_angular
            Gain_err = sprintf('G%d',JJ);
            Gain_err = eval(Gain_err,'=Gain_err');
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
            %Msnr_crosstalkerr = crosstalk_err(Msnr,sensors,crosstalk_range(JJ));
            
            % AP 
            [~,S_AP] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP');
            % AP-MUSIC 
            [~,S_AP_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-MUSIC');
%           % AP-wMUSIC 
            [~,S_AP_W_MUSIC] = alternating_projections(Msnr, ndipoles, Gain_err, nsources,AP_max_iters,'AP-w-MUSIC');
            curr_pos=si(:,i);
            for s=1:nsources
                for n=1:nsources
                    dist_ap(n)       = norm(diff(Cortex.Vertices([curr_pos(s) S_AP(n)],:))); %distance in mm
                    dist_isw(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_W_MUSIC(n)],:))); %distance in mm
                    dist_APM(n)      = norm(diff(Cortex.Vertices([curr_pos(s) S_AP_MUSIC(n)],:))); %distance in mm
                end
                [d_AP(s),  index_AP(s)] = min(dist_ap);             
                [d_ISW(s), index_ISW(s)] = min(dist_isw);
                [d_APM(s), index_APM(s)] = min(dist_APM);
            end

            angularerr_500result_5S_corr1(1,1,SNRindex,i,JJ) = mean(d_AP);
            angularerr_500result_5S_corr1(2,1,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP,index_AP);
            angularerr_500result_5S_corr1(1,2,SNRindex,i,JJ) = mean(d_ISW);
            angularerr_500result_5S_corr1(2,2,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_W_MUSIC,index_ISW);
            angularerr_500result_5S_corr1(1,3,SNRindex,i,JJ) = mean(d_APM);
            angularerr_500result_5S_corr1(2,3,SNRindex,i,JJ) = calculate_RMSE(S,Gain_err,Msnr,S_AP_MUSIC,index_APM);
            
        end

    end
end
save('angularerr_500result_5S_corr1.mat','angularerr_500result_5S_corr1');