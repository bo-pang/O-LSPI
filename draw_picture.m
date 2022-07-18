clc;
clear;
N = 17;
J_rel_avg = zeros(N,1);
J_rel_var = zeros(N,1);
stab_frac = zeros(N,1);
noise_magnitude = zeros(N,1);
for i=1:N
    load(['Recht_Average_Noise_Mag_result_T1e2-',num2str(i),'.mat']);
    sum = 0;
    noise_magnitude(i) = rs.noise_magnitude^2;
    for j=1:100
        if rs.tr{j}.stable
            sum = sum + 1;
        end
    end
    stab_frac(i) = sum/100;
    J_rel = [];
    j=1;
    jhis = [];
    while length(J_rel)<100
        if rs.tr{j}.stable
            J_rel = [J_rel,rs.tr{j}.J_rel(end)];
            jhis = [jhis j];
        end
        j=j+1;
    end
    J_rel_avg(i) = mean(J_rel);
    J_rel_var(i) = var(J_rel);
end

% Delta_G = zeros(rs.I-1,100);
% for i=1:rs.I-1
%     for j=1:100
%         Delta_G(i,j) = rs.tr{j}.Delta_G{i};
%     end
% end
% Delta_G_avg = mean(Delta_G,2);
% Delta_G_var = var(Delta_G,0,2);

N = 17;
J_rel_second_avg = zeros(N,1);
J_rel_second_var = zeros(N,1);
stab_second_frac = zeros(N,1);
for i=1:N
    load(['Recht_Average_Noise_Mag_result_T1e4-',num2str(i),'.mat']);
    sum = 0;
    for j=1:100
        if rs.tr{j}.stable
            sum = sum + 1;
        end
    end
    stab_second_frac(i) = sum/100;
    J_rel = [];
    j=1;
    jhis = [];
    while length(J_rel)<100
        if rs.tr{j}.stable
            J_rel = [J_rel,rs.tr{j}.J_rel(end)];
            jhis = [jhis j];
        end
        j=j+1;
    end
    J_rel_second_avg(i) = mean(J_rel);
    J_rel_second_var(i) = var(J_rel);
end

N1 = 17;
J_rel_LSTD_avg = zeros(N1,1);
J_rel_LSTD_var = zeros(N1,1);
stab_frac_LSTD = zeros(N1,1);
clear sum;
for i=1:N1
    load(['LSTD_Q_',num2str(i),'.mat']);
    stab_frac_LSTD(i) = (100-sum(isinf(fraction(51:150))))/100;
    fraction = fraction(~isinf(fraction));
    J_rel = fraction(51:150);
    J_rel_LSTD_avg(i) = mean(J_rel);
    J_rel_LSTD_var(i) = var(J_rel);
end

N1 = 17;
J_rel_LSTD_second_avg = zeros(N1,1);
J_rel_LSTD_second_var = zeros(N1,1);
stab_frac_LSTD_second = zeros(N1,1);
clear sum;
for i=1:N1
    load(['LSTD_Q_second_',num2str(i),'.mat']);
    stab_frac_LSTD_second(i) = (100-sum(isinf(fraction(51:150))))/100;
    fraction = fraction(~isinf(fraction));
    J_rel = fraction(51:150);
    J_rel_LSTD_second_avg(i) = mean(J_rel);
    J_rel_LSTD_second_var(i) = var(J_rel);
end

figure(1)
subplot(2,3,1);
plot(log10(noise_magnitude),stab_frac,'-bx');
hold on;
plot(log10(noise_magnitude),stab_frac_LSTD,'-ro');
ylim([0 1.2]);
title('Fraction Stable');
ax = gca;
ax.TitleFontSizeMultiplier = 0.95;
xlabel('$\log_{10}(\sigma^2_u)$','interpreter','latex');
yticks([0,0.5,1,1.2]);
xlim([-8,8]);
xticks([-8,-4,0,4,8]);
% legend('M = 10^2');
legend('O-LSPI','LSPIv1');

subplot(2,3,4);
plot(log10(noise_magnitude),stab_second_frac,'-bx');
hold on;
plot(log10(noise_magnitude),stab_frac_LSTD_second,'-ro');
xlim([-8,8]);
ylim([0 1.2]);
title('Fraction Stable');
ax = gca;
ax.TitleFontSizeMultiplier = 0.95;
xlabel('$\log_{10}(\sigma^2_u)$','interpreter','latex');
yticks([0,0.5,1,1.2]);
xticks([-8,-4,0,4,8]);
% legend('M = 10^2');
legend('O-LSPI','LSPIv1');

subplot(2,3,2);
plot(log10(noise_magnitude),J_rel_avg,'-bx');
hold on;
plot(log10(noise_magnitude),J_rel_LSTD_avg,'-ro');
title('Relative Error (Avg.)');
ax = gca;
ax.TitleFontSizeMultiplier = 0.95;
xlabel('$\log_{10}(\sigma^2_u)$','interpreter','latex');
xlim([-8,8]);
xticks([-8,-4,0,4,8]);
% legend('M = 10^2');
legend('O-LSPI','LSPIv1');

subplot(2,3,5);
plot(log10(noise_magnitude),J_rel_second_avg,'-bx');
hold on;
plot(log10(noise_magnitude),J_rel_LSTD_second_avg,'-ro');
title('Relative Error (Avg.)');
ax = gca;
ax.TitleFontSizeMultiplier = 0.95;
xlabel('$\log_{10}(\sigma^2_u)$','interpreter','latex');
xlim([-8,8]);
xticks([-8,-4,0,4,8]);
% legend('M = 10^4');
legend('O-LSPI','LSPIv1');

subplot(2,3,3);
plot(log10(noise_magnitude),J_rel_var,'-bx');
hold on;
plot(log10(noise_magnitude),J_rel_LSTD_var,'-ro');
title('Relative Error (Var.)');
ax = gca;
ax.TitleFontSizeMultiplier = 0.95;
xlabel('$\log_{10}(\sigma^2_u)$','interpreter','latex');
xlim([-8,8]);
xticks([-8,-4,0,4,8]);
% legend('M = 10^4');
legend('O-LSPI','LSPIv1');

subplot(2,3,6);
plot(log10(noise_magnitude),J_rel_second_var,'-bx');
hold on;
plot(log10(noise_magnitude),J_rel_LSTD_second_var,'-ro');
title('Relative Error (Var.)');
ax = gca;
ax.TitleFontSizeMultiplier = 0.95;
xlabel('$\log_{10}(\sigma^2_u)$','interpreter','latex');
xlim([-8,8]);
xticks([-8,-4,0,4,8]);
% legend('M = 10^4');
legend('O-LSPI','LSPIv1');
