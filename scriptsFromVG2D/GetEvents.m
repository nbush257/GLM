function [Mevents,s,seh,d,deh,t,teh,q,qeh,quints,quintseh] = GetEvents(spikevec); 
% [Mevents,s,seh,d,deh,t,teh,q,qeh,quints,quintseh] = GetEvents(spikevec); 

s = strfind(spikevec,1);
seh = strfind(spikevec,[0 1 0]);seh = seh +1;
d =strfind(spikevec,[1 1]);
deh = strfind(spikevec,[0 1 1 0]);deh = deh+1;
t = strfind(spikevec,[1 1 1]);
teh = strfind(spikevec,[0 1 1 1 0]);teh = teh+1;
q = strfind(spikevec,[1 1 1 1]);
qeh = strfind(spikevec,[0 1 1 1 1 0]);qeh = qeh+1; 

quints = strfind(spikevec,[1 1 1 1 1 ]);quints   = quints +1; 
quintseh = strfind(spikevec,[0 1 1 1 1 1 0]);quintseh = quintseh+1; 
if ~isempty(quints)
  disp('     Heads up, there''s a quint');
end;

sextuplets = strfind(spikevec,[1 1 1 1 1 1 ]);quints   = quints +1; 
sextupleteh = strfind(spikevec,[0 1 1 1 1 1 1 0]);quintseh = quintseh+1; 
if ~isempty(sextuplets)
  disp('     Heads up, there''s a sextuplet');
end;

septuplets = strfind(spikevec,[1 1 1 1 1 1 1 ]);quints   = quints +1; 
septupleteh = strfind(spikevec,[0 1 1 1 1 1 1 1 0]);quintseh = quintseh+1; 
if ~isempty(septuplets)
  disp('     Heads up, there''s a septuplet');
end;

Mevents.s = s;
Mevents.seh = seh;
Mevents.d = d;  
Mevents.deh = deh;
Mevents.t = t;
Mevents.teh = teh;
Mevents.q = q;
Mevents.qeh = qeh;

% % Mevents{1} = s;
% % Mevents{2} = seh;
% % Mevents{3} = d;
% % Mevents{4} = deh;
% % Mevents{5} = t;
% % Mevents{6} = teh;
% % Mevents{7} = q;
% % Mevents{8} = qeh;
