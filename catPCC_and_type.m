allprefDir = [];type = [];    rG_dis_all = [];rG_med_all = [];rM_dis_all = [];rM_med_all = []; rB_dis_all = []; rB_med_all = [];
rG_all = [];rM_all = [];
rB_all = [];
for dd = 1:40

    try
        dN = dir([cellType(dd).name(1:end-10) '*simGLM.mat'])
        
        load(dN.name,'rB*','rG*','rM*')
        
    catch
        rG = NaN;
        rM = NaN;
        rB = NaN;
        rG_dis = NaN;
        rG_med = NaN;
        rM_dis = NaN;
        rM_med = NaN;
        rB_dis = NaN;
        rB_med = NaN;
    end
    
    rG_all(dd) = rG;
    rB_all(dd) = rB;
    rM_all(dd) = rM;
    
        if exist('rG_dis','var')
            rG_dis_all = [rG_dis_all rG_dis];
            rG_med_all = [rG_med_all rG_med];
            rM_dis_all = [rM_dis_all rM_dis];
            rB_dis_all = [rB_dis_all rB_dis];
            rM_med_all = [rM_med_all rM_med];
            rB_med_all = [rB_med_all rM_med];
        else
            rG_dis_all = [rG_dis_all NaN];
            rG_med_all = [rG_med_all NaN];
            rM_dis_all = [rM_dis_all NaN];
            rB_dis_all = [rB_dis_all NaN];
            rM_med_all = [rM_med_all NaN];
            rB_med_all = [rB_med_all NaN];

        end
%     
%         if strcmp(cellType(dd).type,'R')
%             type = [type 1];
%         elseif strcmp(cellType(dd).type,'S')
%             type = [ type 0];
%         elseif strcmp(cellType(dd).type,'M')
%             type = [type 2];
%         else
%             type = [type NaN];
%         end
    clear rG rM rB rG_dis rB_dis rM_dis rG_med rM_med rB_med
    
end
type = [0;2;0;2;0;1;2;0;1;1;0;0;1;2;0;2;0;0;0;0;0;1;1;1;0;0;1;1;2;0;0;0;0;0;0;1;1;1;0;0]

allPrefDir =[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,90,225,135,270,45,135,NaN,NaN,NaN,135,180,NaN,NaN,NaN,315,180,135,270,225,270,90,NaN,NaN,180,135]