%{
***************************************************************************
Code for Paper:
Bo Pango, and Zhong-Ping Jiang. "Robust reinforcement learning: A case 
study in linear quadratic regulation." Proceedings of the AAAI Conference
 on Artificial Intelligence. Vol. 35. No. 10. 2021.
***************************************************************************
Please cite this paper if you find this code helpful!
%}

clc;
clear;

for l=1:17
    load(['Recht_Average_Noise_Mag_T1e4-',num2str(l),'.mat']);
%     load('Power_Lewis_v6_1e5.699.mat');
    rs.I = 5;     % number of policy iteration
    rs.Ipi = 45;   % number of iteration for policy evaluation
    rs.A = A;
    rs.B = B;
    rs.C = C;
    rs.Q = Q;
    rs.R = R;
    rs.rollout = T;
    rs.Kini = K;
    rs.tr = {};
    rs.noise_magnitude = noise_magnitude;
    rs.data.zx = ZX;
    rs.data.zz = ZZ;
    rs.data.zt_r = ZT_R;
    rs.Ntrial = length(ZX);


    [n,m] = size(B);
    n1 = m+n+1;
    n2 = n1*(n1+1)/2;
    n3 = n*(n+1)/2;

    [rs.Pstar,~,rs.Kstar] = dare(A,B,Q,R);
    rs.Jstar = trace(C'*rs.Pstar*C);


    for j=1:rs.Ntrial

        tr.zz = rs.data.zz{j};
        tr.zx = rs.data.zx{j};
        tr.zt_r = rs.data.zt_r{j};
        tr.stable = 1;
        data_mat_1 = (tr.zz)\(tr.zx);
        data_mat_2 = (tr.zz)\(tr.zt_r);

        Khat = {};
        Khat{1} = rs.Kini;
        Qhat = {};
        P = {};
        Phat = {};
        F = {};
        Fhat = {};
        for i=1:rs.I-1
            if any(any(isnan(Khat{i}))) || ~all(eig(A-B*Khat{i})<1)
                disp('Instability occur');
                i
                tr.stable = 0;
                break;
        %         return;
            end
            Ptmp = zeros(n);
            for k=1:rs.Ipi-1
                ftmp = data_mat_1*sm2vec(Ptmp)+data_mat_2;
                Ftmp = vec2sm(ftmp,n1);
                Qtmp = Ftmp(1:end-1,1:end-1);%+[Q zeros(n,m);zeros(m,n) R];
                Ptmp = [eye(n),-Khat{i}']*Qtmp*[eye(n);-Khat{i}];
            end
            ftmp = data_mat_1*sm2vec(Ptmp)+data_mat_2;
            Ftmp = vec2sm(ftmp,n1);
            Qtmp = Ftmp(1:end-1,1:end-1);%+[Q zeros(n,m);zeros(m,n) R];
            Qhat{i} = Qtmp;
            Fhat{i} = Ftmp;
            Phat{i} = Ptmp;
            if any(any(isnan(Qtmp(n+1:end,n+1:end)))) || rank(Qtmp(n+1:end,n+1:end)) ~= m
                disp('Singularity occur');
                i
                tr.stable = 0;
            end
            Khat{i+1} = Qtmp(n+1:end,n+1:end)\Qtmp(n+1:end,1:n);

            if any(any(isnan(Khat{end})))
                disp('Instability occur');
                i
                tr.stable = 0;
                break;
            end

            P{i} = dlyap((A-B*Khat{i})',Q+Khat{i}'*R*Khat{i});
            F{i} = [Q+A'*P{i}*A, A'*P{i}*B, zeros(n,1);
                B'*P{i}*A, R+B'*P{i}*B, zeros(m,1);
                zeros(1,n), zeros(1,m), trace(C'*P{i}*C);];
        end
        Delta_G = {};
        K_rel = zeros(rs.I,1);
        J_rel = zeros(rs.I,1);
        if tr.stable==1
            P{i+1} = dlyap((A-B*Khat{i+1})',Q+Khat{i+1}'*R*Khat{i+1});
            F{i+1} = [Q+A'*P{i+1}*A, A'*P{i+1}*B, zeros(n,1);
                B'*P{i+1}*A, R+B'*P{i+1}*B, zeros(m,1);
                zeros(1,n), zeros(1,m), trace(C'*P{i+1}*C);];
            % Plot figures
            for i=1:rs.I-1
                Delta_G{i} = Qhat{i}-F{i}(1:end-1,1:end-1)+[-Phat{i}+P{i},zeros(n,m);
                    zeros(m,n),zeros(m)];
                Delta_G{i} = norm(Delta_G{i},'fro');
                J_rel(i) = norm(trace(C'*P{i}*C)-rs.Jstar,'fro')/rs.Jstar;
                K_rel(i) = norm(Khat{i}-rs.Kstar,'fro')/norm(rs.Kstar,'fro');
            end
            J_rel(i+1) = norm(trace(C'*P{i+1}*C)-rs.Jstar,'fro')/rs.Jstar;
            K_rel(i+1) = norm(Khat{i+1}-rs.Kstar,'fro')/norm(rs.Kstar,'fro');
        end

        tr.Qhat = Qhat;
        tr.Fhat = Fhat;
        tr.Phat = Phat;
        tr.Khat = Khat;
        tr.P = P;
        tr.F = F;
        tr.Delta_G = Delta_G;
        tr.K_rel = K_rel;
        tr.J_rel = J_rel;

        rs.tr{j} = tr;
    end
    save(['Recht_Average_Noise_Mag_Result_T1e4-',num2str(l),'.mat'],'rs');
%     save('Power_Lewis_result.mat','rs');
end