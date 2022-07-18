clc;
clear;

T = 1e3;
noise_magnitude = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,...
    1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8];
ID = length(noise_magnitude);

for i = 1:ID
   data_collect_Noise_Mag(T,sqrt(noise_magnitude(i)),i); 
end