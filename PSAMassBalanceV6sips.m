clc
clear all
close all

% e = linspace(0.2, 0.5, 4);
% for i=1:length(e)
%     getBrC(e(i), i)
% end
% 
% function none=getBrC(e, itera)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Variable Definition
    L = 10; % bed length m
    e = 0.36; % voidage
    Q = 3.626621; % gas flowrate m3/s
    D = 0.65; % bed internal diameter, m
    A = 0.25*pi*D^2; % bed cross section area, m2
    u = Q/A; % superficial velocity m/

    DR2 = [6.2e-4, 338e-4, 4.2e-4]; % rate cst, linear driving force [N2, O2, Ar]

    KL = [3.83e2, 1.46e2, 3.09e2]; % Langmuir constant /bar
    qs = [2.573, 4.066, 2.088]; % Saturation uptake mol/kg
    n = [0.806, 0.894, 0.889];

    P = 8; % pressure bar
    y = [75.56, 23.15, 1.29]/100; %mol fraction 
    po = y.*P; %partial pressure initial
    %po = [2.0114, 1.4274e-12, 0.0203];
    
    R = 8.3145; % J/mol/K
    T = 293; % K
    
    N = 101;
    z = linspace(0, L, N);
    dz = L/N;
    tmax = 300;
    t = linspace(0, tmax, tmax+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solving
    %initial guess
    IC = zeros(1, 6*N);
    
    [t, y] = ode15s(@f,t,IC, [], N, po, KL, qs, DR2, dz, u, e, R, T, P, Q, A,n);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data Extraction and Setting inital boundary conditions
    pN2 = real(y(:,1:N));
    pN2(:,1)=po(1);
    pN2(:,end) = 2*pN2(:,end-1)-pN2(:,end-2);
    
    
    pO2 = real(y(:,N+1:2*N));
    pO2(:,1)=po(2);
    pO2(:,end) = 2*pO2(:,end-1)-pO2(:,end-2);
    
    pAr = real(y(:,2*N+1:3*N));
    pAr(:,1)=po(3);
    pAr(:,end) = 2*pAr(:,end-1)-pAr(:,end-2);
    
    qN2 = y(:,3*N+1:4*N);
    qO2 = y(:,4*N+1:5*N);
    qAr = y(:,5*N+1:6*N);
    
    
    dimCell = num2cell(size(pN2));
    [r, c] = dimCell{:};
    N2pur=zeros(size(pN2));
    for i=1:r
        for j=1:c
            if pN2(i,j) + pO2(i,j) + pAr(i,j) < 1
                N2pur(i,j) = (pN2(i,j) + pO2(i,j) + pAr(i,j))*0.25+0.75;
            else
                % N2pur(i,j) = pN2(i,j)/(pN2(i,j) + pO2(i,j)+ pAr(i,j));
                N2pur(i,j) = pN2(i,j)/(pN2(i,j) + pO2(i,j));
            end
        end
    end
    psum = pN2+pO2+pAr;
    psumArFree = pN2+pO2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data Analysis
    % get height and cycle time
    targ_pur = 0.99999;
    for i=1:c-1
        if N2pur(2,i) > targ_pur && N2pur(2,i) > N2pur(2,i+1)
            H_idx = i;
            H = z(H_idx);
        end
    end

    cum_pur=0;
    N2_cum = 0;
    ArFree_cum = 0;
    for i=1:r
        N2_cum = N2_cum + pN2(i,H_idx);
        ArFree_cum = ArFree_cum + psumArFree(i,H_idx);
        cum_pur = N2_cum/ArFree_cum;
        if cum_pur > targ_pur
            cycle_t = t(i);
            t_idx = i;
        end
    end
    
    nN2=pN2(1:t_idx,H_idx)*1e5*Q/P/R/T;
    nO2=pO2(1:t_idx,H_idx)*1e5*Q/P/R/T;
    nAr=pAr(1:t_idx,H_idx)*1e5*Q/P/R/T;
    n = nN2+nO2+nAr;
    n_avg = sum(n) / cycle_t;
    nN2_avg = sum(nN2) / cycle_t;
    nO2_avg = sum(nO2) / cycle_t;
    nAr_avg = sum(nAr) / cycle_t;
    purity = nN2_avg/(nO2_avg+nN2_avg);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plotting
    %figure(itera)
    tiledlayout(3,1)
    nexttile
    imagesc(z,t,N2pur)
    colormap summer
    hold on
    plot([0 H], [cycle_t, cycle_t], "r-")
    plot([H H], [0, cycle_t], "r-")
    legend('Operating Region')
    ylabel('Time (s)')
    xlabel('Length (m)')
    titletext = sprintf('Bed voidage %.1f, Optimal Height: %.1f m, Half Cycle time %.1f s, mol flow %.1f mol/s, Mean Effluent Purity (AR free): %.7f',e, H, cycle_t, n_avg, purity);
    title(titletext)
    colorbar
    xlim([0 L])
    hold off
    % O2 breakthrough curve
    nexttile
    plot(t,pO2(:,H_idx)/po(2))
    ylabel('O2 p/po')
    xlabel('Time (s)')
    hold on 
    plot(cycle_t, pO2(t_idx,H_idx)/po(2), 'r*')
    plot([cycle_t cycle_t], [1.1 0], 'b-')
    legend('Breakthrough Curve', 'Breakthrough Point', 'Cycle Time')
    ylim([0 1.1])
    hold off
    % N2 breakthrough curve
    nexttile
    plot(t,N2pur(:,H_idx))
    ylabel('Effluent N2 purity (Argon Free)')
    xlabel('Time (s)')
    hold on 
    plot(cycle_t, N2pur(t_idx,H_idx), 'r*')
    plot([cycle_t cycle_t], [1.1 0], 'b-')
    legend('Breakthrough Curve', 'Breakthrough Point', 'Cycle Time')
    ylim([0 1.1])
    hold off
    

    function dydt=f(t, y, N, po, KL, qs, DR2, dz, u, e, R, T, P, Q, A,n)
    dydt = zeros(length(y), 2*length(po));
    dpdt = zeros(N, length(po));
    dqdt = zeros(N, length(po));
    
    
    % 1 = N2, 2 = O2, 3 = Ar
    for i=1:3
        p(:,i) = y((i-1)*N+1:i*N);
        q(:,i) = y((i+2)*N+1:(i+3)*N);
    end
    p(1, :) = po;
    
    
    
    for i = 2:N-1
        
        denom = 1;
        for j = 1:length(po)
            denom = denom + p(i,j)^n(j)* KL(j);
        end
        
    
        for j = 1:length(po)
            q_eq = qs(j)*KL(j)*p(i,j)^n(j)/denom;
            dqdt(i,j) = 15*DR2(j)*(q_eq-q(i,j));
            dpdz(i,j) = (p(i,j)-p(i-1,j)) / (dz);
            dpdt(i,j) = -u*dpdz(i,j)-R*T*((1-e)/e)*dqdt(i,j);
        end
    
    end
    
    dydt=[dpdt(:,1);dpdt(:,2);dpdt(:,3);dqdt(:,1);dqdt(:,2);dqdt(:,3)];
    end
%     none = 0;
% end