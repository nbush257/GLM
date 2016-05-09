
function concatVars(ratNum,trialNum,whisker,RCCR)
%% function concatVars(ratNum,trialNum,whisker,RCCR)
% %Concatenate all mechanical variables
% clear
% %choose the data to concatenate
% For example
% ratNum = '15'
% trialNum = 't01'
% whisker = 'B1'
% % choose the frames of RCCR:
% RCCR = [1 60000+11800];
%%


% find the relevant files
dE2D = dir(['*E2D*' ratNum '*' whisker '*' trialNum '*E2D*.mat']);
dProc = dir(['*PROC*' ratNum '*' whisker '*' trialNum '*.mat']);
if length(dE2D)~=length(dProc)
    error('Mismatched files')
end

lastFileName = dE2D(end).name;
%pull last mechanical frame from the filenames
lastFrame = regexp(lastFileName,'F\d{6}F\d{6}');lastFrame = str2num(lastFileName(lastFrame+8:lastFrame+13));
TAGstart = regexp(lastFileName,'rat');TAGend = regexp(lastFileName,'_t\d{2}');
TAG = lastFileName(TAGstart:TAGend+3)


%%

E2D.FX = nan(lastFrame,1);
E2D.FY = nan(lastFrame,1);
E2D.M = nan(lastFrame,1);
E2D.BP = nan(lastFrame,2);
E2D.C = nan(lastFrame,1);
E2D.TH = nan(lastFrame,1);
E2D.xs{lastFrame} = [];
E2D.ys{lastFrame} = [];
E2D.CP = nan(lastFrame,2);

for ii = 1:length(dE2D)
    load(dE2D(ii).name);
    load(dProc(ii).name);
    
    clipStart = regexp(dE2D(ii).name,'F\d{6}');clipStart = str2num(dE2D(ii).name(clipStart(1)+1:clipStart(1)+6));
    clipEnd = regexp(dE2D(ii).name,'F\d{6}');clipEnd = str2num(dE2D(ii).name(clipEnd(2)+1:clipEnd(2)+6));
%     if ii ==1
%         FX =FX(1:20000);
%         FY =FY(1:20000);
%         M =M(1:20000);
%         BP =BP(1:20000,:);
%         C =C(1:20000,:);
%         TH =TH(1:20000,:);
%         CP =CP(1:20000,:);
%         xs = xs(1:20000);
%         ys = ys(1:20000);
%     end
    
    E2D.FX(clipStart:clipEnd) = FX;
    E2D.FY(clipStart:clipEnd) = FY;
    E2D.M(clipStart:clipEnd) = M;
    E2D.BP(clipStart:clipEnd,:) = BP;
    E2D.C(clipStart:clipEnd) = C;
    E2D.TH(clipStart:clipEnd) = TH;
    E2D.xs(clipStart:clipEnd) = xs;
    E2D.ys(clipStart:clipEnd) = ys;
    E2D.CP(clipStart:clipEnd,:) = CP;
end

fNames = fieldnames(E2D);
for ii = 1:length(fNames)
    if ismember(fNames{ii},{'C','xs','ys'})
        E2D_medfilt.(fNames{ii}) = E2D.(fNames{ii});
        continue
    end
    
    E2D_medfilt.(fNames{ii}) = medfilt1(E2D.(fNames{ii}));
end
newRCCR = zeros(length(E2D.C),1);
for ii = 1:size(RCCR,1)
    newRCCR(RCCR(ii,1):RCCR(ii,2)) = 1;
end
RCCR = newRCCR;
ow = 'y';
if exist([TAG '_varConcat.mat'],'file')
    ow = input('Overwrite the existing file? (y/n)','s')
end

if strcmp(ow,'y')
    save([TAG '_varConcat.mat'],'E2D*','RCCR')
end


