% This code extract lifetimes from each monomer's population evolution of RF results

%% Get residue information
C = calc_coupling(atom);

%% Get notable connection
don = {'N','Y'}; acc = {'S'}; % Donor and acceptor chain
kcut = 1/50;
exsite = (RF.U.^2)==max(RF.U.^2); %Exciton-site-basis convertor
for i = 1:length(don)
    for j = 1:length(acc)
        fprintf('\nFinding %s-%s connection faster than %.1f ps:\n',don{i},acc{j},1/kcut); 
        idon = RF.chain==don{i};
        iacc = RF.chain==acc{j};
        Kforw = RF.K.*(iacc.*idon');
        Kback = RF.K.*(idon.*iacc');
        Ktot = Kforw + Kback';
        [row_ex,col_ex] = find(Ktot>kcut);
        [row_site,~] = find(exsite(:,row_ex)==1);
        [col_site,~] = find(exsite(:,col_ex)==1);
        for n = 1:length(row_ex)
            fprintf('%s and %s - %.1f ps - %.1f Angtrom\n',[C.chain(col_site(n)) '.' C.resname{col_site(n)} num2str(C.molid(col_site(n)))],[C.chain(row_site(n)) '.' C.resname{row_site(n)} num2str(C.molid(row_site(n)))],1/Ktot(row_ex(n),col_ex(n)),10*C.R(row_site(n),col_site(n)));
        end
    end
end
   
%% Get average distance between complexes
don = {'Y'}; acc = {'C'}; % Donor and acceptor chain
for i = 1:length(don)
    for j = 1:length(acc)
        if length(don{i})>1
            idon = zeros(length(C.chain),1);
            for k = 1:length(don{i})
                idon = idon + (C.chain==don{i}(k));
            end
            idon = logical(idon);
        else
            idon = C.chain==don{i};
        end
        if length(acc{j})>1
            iacc = zeros(length(C.chain),1);
            for k = 1:length(acc{j})
                iacc = iacc + (C.chain==acc{j}(k));
                sum(iacc)
            end
            iacc = logical(iacc);
        else
            iacc = C.chain==acc{j};
        end
        R = C.R(idon,iacc)*10;
        Ndon = sum(idon); Nacc = sum(iacc); 
        avgR = sum(sum(R))/2/Ndon/Nacc;
        fprintf('\nAverage distance between %s and %s : %.1f Angtrom\n',don{i},acc{j},avgR);
    end
end

%% Avg EET
% Based on averaging the don-acc connection in K matrix
don = {'S'}; acc = {'Y'}; % Donor and acceptor chain
for i = 1:length(don)
    for j = 1:length(acc)
        forw = 1./mean(sum(RF.K(RF.chain==acc{j},RF.chain==don{i})));
        back = 1./mean(sum(RF.K(RF.chain==don{j},RF.chain==acc{i})));
        fprintf('Chain %s to %s: %.2f ps\n',don{i},acc{j},forw);
        fprintf('Chain %s to %s: %.2f ps\n',acc{i},don{j},back);
        fprintf('Total rate: %.2f ps\n',1/(1/forw+1/back));
    end
end

% Based on simulated kinetics
tt = Par.t;
for i = 1:length(don)
    for j = 1:length(acc)
        idon = RF.chain==don{i};
        iacc = RF.chain==acc{j};
        id = logical(idon+iacc);
        p0 = id;
        for nid = 1:length(id)
            if strcmp(C.resname{nid},'CHL')
                p0(nid) = p0(nid)*15/21;
            end
        end
        p = zeros(length(id),length(tt));
        for nK = 1:size(RF.Kfull,3)
            Ktrim = squeeze(RF.Kfull(:,:,nK)).*(id.*id');
            Ktrim = Ktrim - diag(sum(Ktrim)) - diag(diag(Ktrim));
            p = p + popkin(Ktrim,tt,p0);
        end
        p = p/nK;
        pdon = sum(p(idon,:));
        pacc = sum(p(iacc,:));
        fitresdon = fitexp(tt,pdon);
        fitresacc = fitexp(tt,pacc);
        figure; plot(tt,pdon,'o',tt,pacc,'o');
        hold on;
        set(gca,'colororderindex',1,'xscale','log');
        plot(tt,feval(fitresdon,tt),tt,feval(fitresacc,tt));
        title(sprintf('Chain %s and %s: %.2f ps',don{i},acc{j},fitresdon.b));
        legend(don{i},acc{j},'location','best')
    end
end

%% Population from RF.P
ysel_ls = [650];
fitresult = cell(length(unique(RF.chain)),length(ysel_ls));
gof = fitresult;
for nysel = 1:length(ysel_ls)
    ysel = ysel_ls(nysel);
    if ysel<=650
        dy = 50;
    else
        dy = 3;
    end
    exc_range = unit_conv([ysel_ls(nysel)-dy ysel_ls(nysel)+dy],'nm to cm-1');
    [fitresult(:,nysel),gof(:,nysel),pop] = RFPpop(exc_range,Par,RF,1);
end

% Plot
figure; formatfig('default');
tt = Par.t;
ch = unique(RF.chain);
sparse = 1:2:length(tt);
for i = 1:length(ch)
    h = plot(tt(sparse),pop(i,sparse),'o','displayname',['Chain ' ch(i) ', ' sprintf('%.0f ps',fitresult{i}.b)]);
    hold on;
    plot(tt,feval(fitresult{i},tt),'color',h.Color,'linewidth',1.2,'displayname','none');
end
set(gca,'xscale','log','fontsize',14);
xlabel('Time (ps)'); ylabel('Population');
xlim([10^-2,10^3]);
xticks([0.01, 0.1, 1, 10, 100, 1000]);
xticklabels(string([0.01, 0.1, 1, 10, 100, 1000]));
ylim([0 .6]);
box on;
legend('location','best');

%% Population averaged from RF.Kfull (gives same results as RF.P)
tt = Par.t;
ch = unique(RF.chain);
pt0 = ones(length(RF.K),1);
for nid = 1:length(pt0)
    if strcmp(C.resname{nid},'CHL')
        pt0(nid) = pt0(nid)*15/21;
    end
end
pt = zeros(length(pt0),length(tt));
for nK = 1:size(RF.Kfull,3)
    K = RF.Kfull(:,:,nK);
    % Block connection
    %------
    idon = RF.chain=='Y';
    iacc = RF.chain=='C';
    K(idon,iacc) = 0; K(iacc,idon) = 0;
    %------
    K = K - diag(diag(K)) - diag(sum(K));
pt = pt + popkin(K,tt,pt0);
end
pop = zeros(length(ch),length(tt));
for i = 1:length(ch)
    pop(i,:) = sum(pt(RF.chain==ch(i),:));
end
pop = pop/sum(pop(:,1));

% Fit exp
fitresult = cell(size(ch)); gof = cell(size(ch));
for i = 1:length(ch)
    [fitresult{i}, gof{i}] = fitexp(tt, pop(i,:));
end

% Plot
figure; formatfig('default');
sparse = 1:2:length(tt);
for i = 1:length(ch)
    h = plot(tt(sparse),pop(i,sparse),'o','displayname',ch(i));
    hold on;
    plot(tt,feval(fitresult{i},tt),'color',h.Color,'linewidth',1.2,'displayname',sprintf('%.1f ps',fitresult{i}.b));
end
set(gca,'xscale','log','fontsize',14);
xlabel('Time (ps)'); ylabel('Population');
xlim([10^-2,10^3]);
box on;
legend('location','best');

%% Distribution of pairwise Chl distance
chainID = 'C'; 
include_b = false;
idon = RF.chain==chainID;
if ~include_b
    idon = idon.*strcmp(C.resname,'CLA');
    idon = logical(idon);
end
Rmatrix = triu(C.R(idon,idon));
Rmatrix = Rmatrix(:);
Rmatrix(Rmatrix==0) = [];
figure;
histogram(Rmatrix,7);
xlabel('Distance (nm)');
title(sprintf('Distribution of pairwise Chl distance in chain %c',chainID));



%% Script function
function [fitresult,gof,pop] = RFPpop(exc_range,Par,RF,tstart)
% Fit population kinetic from RF.P
exci = ind(Par.exc,exc_range(2)):ind(Par.exc,exc_range(1));
ch = unique(RF.chain);
tt = Par.t;
pop = zeros(length(ch),length(tt));
for i = 1:length(ch)
    pop(i,:) = sum(sum(RF.P(RF.chain==ch(i),:,exci),3));
end
pop = pop/sum(pop(:,1));
% pop = pop + randn(size(pop))*0.005; % Add noise

% Fit exp
fitresult = cell(size(ch)); gof = cell(size(ch));
for i = 1:length(ch)
    [fitresult{i}, gof{i}] = fitexp(tt(tt>=tstart), pop(i,tt>=tstart));
end
end