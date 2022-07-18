function data_collect_Noise_Mag(T,noise_magnitude,ID)
N = 500;
ZZ = {};
ZX = {};
ZT_R = {};

for i=1:N-1
    [ZZ{i},ZX{i},ZT_R{i},~,~,~,~,~,~] = func_data_collect(T,noise_magnitude);
end
[ZZ{i+1},ZX{i+1},ZT_R{i+1},A,B,C,Q,R,K] = func_data_collect(T,noise_magnitude);
% save(['/scratch/bp1471/Ergodic/Recht_Average_1.mat'],'ZZ','ZX','K',...
%     'A','B','C','Q','R','ZT_R','T','noise_magnitude');
save(['Recht_Average_Noise_Mag_T1e',num2str(log10(T)),'-',num2str(ID),'.mat'],...
    'ZZ','ZX','K','A','B','C','Q','R','ZT_R','T','noise_magnitude');
end