function RMSE_wave = calculate_RMSE(wave_true,Gain,MEG_data,dipole_estimate,index)
% calculate_RMSE calculate the relative mean squared error between the simulated wave and the
% estimated wave
%% get the normalized simulated wave in the source space
rowNorm_true = sqrt(sum(wave_true.^2,2));
rowNorm_true(rowNorm_true==0)=1;
wave_true_norm = bsxfun(@rdivide,wave_true,rowNorm_true);
%% get the normalized estimated wave in the source space
Ap = Gain(:,dipole_estimate);
wave_recon = inv(Ap'*Ap)*Ap'*MEG_data;
wave_recon = wave_recon(index,:);
rowNorm_estimate = sqrt(sum(wave_recon.^2,2));
rowNorm_estimate(rowNorm_estimate==0)=1;
wave_recon_norm = bsxfun(@rdivide,wave_recon,rowNorm_estimate);
%% calculate the mean RMSE between the simulated wave and the estimated wave
for i = 1:size(dipole_estimate,2)
    RMSE(i) = (norm((wave_true_norm(i,:)-wave_recon_norm(i,:)),2)/norm(wave_true_norm(i,:),2))^2;
end
RMSE_wave = mean(RMSE);