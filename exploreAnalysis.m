function  exploreAnalysis(Mech,Geo,prox,neuralWord,cNum)
% function [R2,vars,meanTh,meanRho] = exploreAnalysis(Mech,Geo,prox,neuralWord,cNum)
X = Mech.all;
X2 = [Mech.Fx sqrt(Mech.Fy.^2+Mech.Fz.^2) unwrap(atan2(Mech.Fz,Mech.Fy)) Mech.Mx sqrt(Mech.My.^2+Mech.Mz.^2) unwrap(atan2(Mech.Mz,Mech.My))];
X3 = Geo.all;

C = all(~isnan(X),2);
C = C(~prox);
X = X(~prox,:);
X2 = X2(~prox,:);
X3 = X3(~prox,:);
neuralWord = neuralWord(~prox);


cstart = find(diff(C)==1)+1;
cend = find(diff(C)==-1);
cLength = [];
cMech = {};
cSpikes = {};
cFR = [];
toRM  =[];
for ii = 1:length(cstart)
    clength(ii) = cend(ii)-cstart(ii);
    if clength(ii)<30 
        toRM = [toRM ii];
        continue
    end
    cSpikes{ii} = neuralWord(cstart(ii):cend(ii));
    
%     cFR(ii) = sum(cSpikes{ii})/length(cSpikes{ii})*1000;
%     if cFR(ii) ==0
%         toRM = [toRM ii];
%     end
end
cFR = [];cSpikes = {};
clength(toRM) = [];cstart(toRM) = []; cend(toRM) = [];

for ii = 1:length(cstart)
    clength(ii) = cend(ii)-cstart(ii);
    
    cMech{ii} = Mech.all(cstart(ii):cend(ii),:);
    cSpikes{ii} = neuralWord(cstart(ii):cend(ii));
    cFR(ii) = sum(cSpikes{ii})/length(cSpikes{ii})*1000;
    meanISI(ii) = mean(diff(find(cSpikes{ii})));
end

maxMech = nan(length(cstart),size(Mech.all,2));
minMech = maxMech;
rangeMech = nan(length(cstart),size(Mech.all,2));
for ii = 1:length(cMech)
    
    maxMech(ii,:) = nanmax(cMech{ii});
    minMech(ii,:) = nanmin(cMech{ii});
    
    rangeMech(ii,:) = nanmax(cMech{ii})-cMech{ii}(1,:);
end
absMaxMech = maxMech;
absMaxMech(abs(minMech)>maxMech) = minMech(abs(minMech)>maxMech);
Fdir = [];
Mdir = [];
FMag = [];
MMag = [];
RR = [];
MX = [];
FX = [];
R = [];
TH = [];
PHI = [];
dum = [];

for ii = 1:length(cstart)
    clip3 = X3(cstart(ii):cend(ii),:)-repmat(X3(cstart(ii),:),cend(ii)-cstart(ii)+1,1);
    
    R(ii) = nanmean(X3(cstart(ii):cend(ii),1));
    
    [~,dum(ii)] = nanmax(abs(clip3(:,2)));
    TH(ii) =  X3(cstart(ii)+dum(ii)-1,2)-X3(cstart(ii),2);
    
    [~,dum] = nanmax(abs(clip3(:,4)));
    PHI(ii) = X3(cstart(ii)+dum-1,4)-X3(cstart(ii),4);
    
    
    
    
    [~,dum] = nanmax(abs(X(cstart(ii):cend(ii),6)));
    FZ(ii) =  X(cstart(ii)+dum-1,6)-X(cstart(ii),6);
    
    [~,dum] = nanmax(abs(X(cstart(ii):cend(ii),3)));
    MZ(ii) =  X(cstart(ii)+dum-1,3)-X(cstart(ii),3);
    
    [~,dum] = nanmax(abs(X(cstart(ii):cend(ii),5)));
    FY(ii) =  X(cstart(ii)+dum-152)-X(cstart(ii),5);
    
    [~,dum] = nanmax(abs(X(cstart(ii):cend(ii),2)));
    MY(ii) = X(cstart(ii)+dum-1,2)-X(cstart(ii),2);
    
    FMag(ii) = max(X2(cstart(ii):cend(ii),2));
    Fdir(ii) = mean(X2(cstart(ii):cend(ii),3));
     MMag(ii) = max(X2(cstart(ii):cend(ii),5));
    Mdir(ii) = mean(X2(cstart(ii):cend(ii),6));
    
    clip2 = X2(cstart(ii):cend(ii),:)-repmat(X2(cstart(ii),:),cend(ii)-cstart(ii)+1,1);

    [~,fxdum(ii)] = nanmax(abs(clip2(:,1)));
    FX(ii) = X2(cstart(ii)+fxdum(ii)-1,1);
    [~,mxdum] = nanmax(abs(clip2(:,4)));
    MX(ii) = X2(cstart(ii)+mxdum-1,4)-X2(cstart(ii),4);
    
    
end
% 
% mdlFX = fitlm(FX,cFR);
% mdlMX = fitlm(MX,cFR);
% mdlMMag = fitlm(MMag,cFR);
% mdlMdir = fitlm(Mdir,cFR);
% mdlFMag = fitlm(FMag,cFR);
% mdlFdir = fitlm(Fdir,cFR);
% % 
% fullMDL = stepwiselm([MMag' Mdir' FMag' Fdir' MX' FX'],cFR,'Criterion','aic','Upper','linear');
% 
% R2.FX = mdlFX.Rsquared.Adjusted;
% R2.MX = mdlMX.Rsquared.Adjusted;
% R2.MMag = mdlMMag.Rsquared.Adjusted;
% R2.Mdir = mdlMdir.Rsquared.Adjusted;
% R2.FMag = mdlFMag.Rsquared.Adjusted;
% R2.Fdir = mdlFdir.Rsquared.Adjusted;
% R2.full = fullMDL.Rsquared.Adjusted;
% vars = fullMDL.CoefficientNames
% noSpikeMY = MY(cFR ==0);
% noSpikeMZ = MZ(cFR ==0);
noSpikeMX = MX(cFR ==0);
% 
% noSpikeFY = FY(cFR ==0);
% noSpikeFZ = FZ(cFR ==0);
noSpikeFX = FX(cFR ==0);

noSpikeMMag = MMag(cFR==0);
noSpikeFMag = FMag(cFR==0);

noSpikeMdir = MMag(cFR==0);
noSpikeFdir = FMag(cFR==0);

% figure
% scatter(noSpikeFX,noSpikeMMag,ones(size(noSpikeFX))*20,'k','filled')
% ho
% scatter(FX(cFR~=0),MMag(cFR~=0),ones(size(FX(cFR~=0)))*20,cFR(cFR~=0),'filled')
% axis square
% colormap('cool')
% colorbar
% title(['Cell ' num2str(cNum)])
% xlabel('F_x')
% ylabel('Magnitude of Bending Moment')
%% Outlier removal
toRM = [];
if cNum == 8 
    toRM = find(FX<-80/100000);
elseif cNum == 10;
    toRM = find(FX<-30/100000);
end
FX(toRM) = [];
cFR(toRM) = [];
TH(toRM) = [];
R(toRM) = [];
PHI(toRM) =[];
MMag(toRM) = [];
Mdir(toRM) = [];

out.noSpikeFX = noSpikeFX;
out.noSpikeMdir = noSpikeMdir;
out.noSpikeMMag = noSpikeMMag;
out.cFR = cFR;
out.spikeMdir = Mdir(cFR~=0);
out.MMag = MMag(cFR~=0);
out.FX = FX(cFR~=0);
[aa,bb] = pol2cart(Mdir,MMag);
%% Theta v R 
% TH = deg2rad(TH);
[THTH,RR] = pol2cart(TH,R);

h = figure
% polar(TH,R)
set(gca,'fontsize',24)
ho
% rm(1)
scatter(TH(cFR==0),R(cFR==0),ones(size(cFR(cFR==0)))*5,'k','filled')
scatter(TH(cFR~=0),R(cFR~=0),ones(size(cFR(cFR~=0)))*20,cFR(cFR~=0),'filled')
colormap('cool')
t = findall(gcf,'type','text');
% delete(t)
% colorbar
% title(['Cell ' num2str(cNum)])
xx = xlim;
yy = ylim;
set(gca,'ylim',yy)

print(['TH_R_Cell_' num2str(cNum) '.tif'],'-dtiff','-r600')

%% Phi v R
% PHI = deg2rad(PHI);
[PHPH,RR] = pol2cart(PHI,R);

h2= figure
% polar(PHI,R)

ho
% rm(1)
scatter(PHI(cFR==0),R(cFR==0),ones(size(cFR(cFR==0)))*5,'k','filled')
scatter(PHI(cFR~=0),R(cFR~=0),ones(size(cFR(cFR~=0)))*20,cFR(cFR~=0),'filled')

t = findall(gcf,'type','text');
% delete(t)
colormap('cool')
colorbar
set(gca,'ylim',yy)

% title(['Cell ' num2str(cNum)])
print(['PHI_R_Cell_' num2str(cNum) '.tif'],'-dtiff','-r600')

%% 2D plots 

figure
scatter(FX(cFR==0)*100000,MMag(cFR==0)*100000,ones(size(FX(cFR==0)))*5,'k','filled')
ho
scatter(FX(cFR~=0)*100000,MMag(cFR~=0)*100000,ones(size(FX(cFR~=0)))*20,cFR(cFR~=0),'filled')
colormap('cool')
title(['Cell ' num2str(cNum)])
print(['FXvMMag_Cell_' num2str(cNum) '.tif'],'-dtiff','-r600')
figure
polar(aa*100000,bb*100000);rm(1)

t = findall(gcf,'type','text');
% delete(t)
ho
scatter(aa(cFR==0)*100000,bb(cFR==0)*100000,ones(size(FX(cFR==0)))*5,'k','filled');
ho
scatter(aa(cFR~=0)*100000,bb(cFR~=0)*100000,ones(size(FX(cFR~=0)))*20,out.cFR(cFR~=0),'filled')
colormap('cool')
title(['Cell ' num2str(cNum)])
print(['Mdir_v_Mmag_Cell_' num2str(cNum) '.tif'],'-dtiff','-r600')


% %% cylinder 
% figure
% scatter3(aa(cFR==0),bb(cFR==0),FX(cFR==0),ones(size(FX(cFR==0)))*20,'k','filled');
% ho
% scatter3(aa(cFR~=0),bb(cFR~=0),FX(cFR~=0),ones(size(FX(cFR~=0)))*20,out.cFR(cFR~=0),'filled')
% colormap('cool')
% title(['Cell ' num2str(cNum)])
%%
% % 
% 
% 
% %% Fy v FZ
% figure
% subplot(2,2,1)
% scatter(FY(cFR==0),FZ(cFR==0),ones(size(FY(cFR==0)))*20,'k','filled')
% ho
% scatter(FY(cFR~=0),FZ(cFR~=0),ones(size(FX(cFR~=0)))*20,cFR(cFR~=0),'filled')
% axis square
% colormap('cool')
% colorbar
% title(['Cell ' num2str(cNum)])
% xlabel('FY')
% ylabel('FZ')
% 
% %% My v Mz
% subplot(2,2,2)
% scatter(MY(cFR==0),MZ(cFR==0),ones(size(MY(cFR==0)))*20,'k','filled')
% ho
% scatter(MY(cFR~=0),MZ(cFR~=0),ones(size(MY(cFR~=0)))*20,cFR(cFR~=0),'filled')
% axis square
% colormap('cool')
% colorbar
% title(['Cell ' num2str(cNum)])
% xlabel('MY')
% ylabel('MZ')
% %% FY v My
% 
% subplot(2,2,3)
% scatter(FY(cFR==0),MY(cFR==0),ones(size(MY(cFR==0)))*20,'k','filled')
% ho
% scatter(FY(cFR~=0),MY(cFR~=0),ones(size(MY(cFR~=0)))*20,cFR(cFR~=0),'filled')
% axis square
% colormap('cool')
% colorbar
% title(['Cell ' num2str(cNum)])
% xlabel('FY')
% ylabel('MY')
% %% Fz Mz
% 
% subplot(2,2,4)
% scatter(FY(cFR==0),MZ(cFR==0),ones(size(MZ(cFR==0)))*20,'k','filled')
% ho
% scatter(FY(cFR~=0),MZ(cFR~=0),ones(size(MZ(cFR~=0)))*20,cFR(cFR~=0),'filled')
% axis square
% colormap('cool')
% colorbar
% title(['Cell ' num2str(cNum)])
% xlabel('FZ')
% % ylabel('MZ')
% 
% figure
% scatter3(FX(cFR==0),MY(cFR==0),MZ(cFR==0),ones(size(FX(cFR==0)))*5,'k','filled')
% ho
% scatter3(FX(cFR~=0),MY(cFR~=0),MZ(cFR~=0),ones(size(FX(cFR~=0)))*20,cFR(cFR~=0),'filled')
% colormap('cool')
 % 
% print(['ComponentCompare_cell_' num2str(cNum) '.tif'],'-dtiff','-r600') 
% % % figure
% % h = polar(Fdir,cFR,'k.')
% [aa,bb] = pol2cart(Fdir,cFR);
% n = norm([aa' bb']);
% aa = sum(aa);bb = sum(bb);
% aaa = [aa,bb]/n;
% [meanTh,meanRho]=cart2pol(aaa(1),aaa(2));
% ho
% polar([0 meanTh],[0 meanRho],'r-*')
% 
% set(h,'markersize',15)
% title(['Cell ' num2str(cNum)])
% % print(['Fdir_cell_' num2str(cNum) '.tif'],'-dtiff','-r600') 
