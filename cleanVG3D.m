Mf = {'Mx','My','Mz','Fx','Fy','Fz'};
Gf = {'R','TH','Zeta','Phi'};
thresh = 6;

for ii = 1:length(Mf)
    rem = Mech.(Mf{ii})>nanstd(Mech.(Mf{ii}))*thresh | Mech.(Mf{ii})<nanstd(Mech.(Mf{ii}))*-thresh;
    
    Mech.(Mf{ii})(rem) = nan;
    Mech.(Mf{ii}) = InterpolateOverNans(Mech.(Mf{ii}),20);
end


for ii = 1:length(Gf)
    rem = Geo.(Gf{ii})>nanstd(Geo.(Gf{ii}))*thresh | Geo.(Gf{ii})<nanstd(Geo.(Gf{ii}))*-thresh;
    Geo.(Gf{ii})(rem) = nan;
end
   
Mech.all = [Mech.Mx Mech.My Mech.Mz Mech.Fx Mech.Fy Mech.Fz];
Geo.all = [Geo.R Geo.TH Geo.Zeta Geo.Phi];

    