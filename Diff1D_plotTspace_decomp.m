%%% Calculate temperature binned heat fluxes and plot:

% Restrict times:
% $$$ d1 = 105;
% $$$ d2 = 165;
d1 = 0;
d2 = 45;
[tmp tII] = min(abs(t/86400-d1));
[tmp tFF] = min(abs(t/86400-d2));
Nt = tFF-tII+1;

T = T(:,tII:tFF);
kt = kt(:,tII:tFF);

% Calculate J and dJdz for boundary layer:
DTDZ = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
Jbl = zeros(Nz+1,Nt);
Jbl(2:end-1,:) = -kt_bl(2:(end-1),:).*DTDZ*Cp*rho0;
dJdzbl = (Jbl(2:end,:)-Jbl(1:(end-1),:)); % Wm-2 within each grid cell

% Calculate J and dJdz for interior:
DTDZ = (T(2:end,:)-T(1:(end-1),:))./repmat(z_rho(2:end)-z_rho(1:(end-1)),[1 ...
                    Nt]);
Jint = zeros(Nz+1,Nt);
Jint(2:end-1,:) = -kt_int(2:(end-1),:).*DTDZ*Cp*rho0;
dJdzint = (Jint(2:end,:)-Jint(1:(end-1),:)); % Wm-2 within each grid cell


% Calculate J non-local:
Jnl = zeros(Nz+1,Nt);
Jnl = -gamt.*kt*Cp*rho0;
dJnldz = (Jnl(2:end,:)-Jnl(1:(end-1),:)); % Wm-2 within each grid cell

% Define temperature grid:
dTi = 0.5;
Tie = 10:dTi:30;
Ti = avg(Tie);
NTi = length(Ti);

% Do temperature binning:
JblT = zeros(NTi,Nt);
JintT = zeros(NTi,Nt);
JnlT = zeros(NTi,Nt);

for ii = NTi:-1:1
    ii    
    for zi = 1:Nz
        inds = T(zi,:)>Tie(ii) & T(zi,:)<=Tie(ii+1);
        JblT(ii,:) = JblT(ii,:)+dJdzbl(zi,:).*(inds*1);
        JintT(ii,:) = JintT(ii,:)+dJdzint(zi,:).*(inds*1);
        JnlT(ii,:) = JnlT(ii,:)+dJnldz(zi,:).*(inds*1);
    end
end
JblT = cat(1,zeros(1,Nt),cumsum(JblT));
JintT = cat(1,zeros(1,Nt),cumsum(JintT));
JnlT = cat(1,zeros(1,Nt),cumsum(JnlT));

% Do temperature binning of mean:
Tmean = mean(T,2);
dJdzblmean = mean(dJdzbl,2);
JTblmean = zeros(NTi,1);
dJdzintmean = mean(dJdzint,2);
JTintmean = zeros(NTi,1);
dJnldzmean = mean(dJnldz,2);
JnlTmean = zeros(NTi,1);
for ii = NTi:-1:1
    for zi = 1:Nz
        inds = Tmean(zi)>Tie(ii) & Tmean(zi)<=Tie(ii+1);
        JTblmean(ii) = JTblmean(ii)+dJdzblmean(zi).*(inds*1);
        JTintmean(ii) = JTintmean(ii)+dJdzintmean(zi).*(inds*1);
        JnlTmean(ii) = JnlTmean(ii)+dJnldzmean(zi).*(inds*1);
    end
end
JTblmean = cat(1,zeros(1,1),cumsum(JTblmean));
JTintmean = cat(1,zeros(1,1),cumsum(JTintmean));
JnlTmean = cat(1,zeros(1,1),cumsum(JnlTmean));

% Do some plotting in T and z space:
figure;
tind = 50;
subplot(1,2,1);
plot(Jbl(:,tind),z_w,'--k', 'DisplayName','BL Diffusive Snapshot');
hold on;
plot(mean(Jbl,2),z_w,'-k','linewidth',2,'DisplayName','BL Diffusive Mean');

plot(Jint(:,tind),z_w,'--b', 'DisplayName','Int Diffusive Snapshot');
plot(mean(Jint,2),z_w,'-b','linewidth',2,'DisplayName','Int Diffusive Mean');

plot(Jnl(:,tind),z_w,'--r','DisplayName','Nonlocal Snapshot');
plot(mean(Jnl,2),z_w,'-r','linewidth',2, 'DisplayName','Nonlocal Mean');
legend('Location','southeast','fontsize',9);

subplot(1,2,2);
plot(JblT(:,tind),Tie,'--k', 'DisplayName','BL Diffusive Snapshot');
hold on;
plot(mean(JblT,2),Tie,'-k','linewidth',2,'DisplayName','BL Diffusive Online');
plot(JTblmean,Tie,':k','linewidth',2,'DisplayName','BL Diffusive Mean');
plot(JintT(:,tind),Tie,'--b', 'DisplayName','Int Diffusive Snapshot');
plot(mean(JintT,2),Tie,'-b','linewidth',2,'DisplayName','Int Diffusive Online');
plot(JTintmean,Tie,':b','linewidth',2,'DisplayName','Int Diffusive Mean');
plot(JnlT(:,tind),Tie,'--r','DisplayName','Nonlocal Snapshot');
plot(mean(JnlT,2),Tie,'-r','linewidth',2,'DisplayName','Nonlocal Online');
plot(JnlTmean,Tie,':r','linewidth',2,'DisplayName','Nonlocal Mean');

legend('Location','southeast','fontsize',9);
