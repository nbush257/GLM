function x_idx = getIndexByTime(x,time)
a = x;clear x
b = time;clear time
if ~isrow(a);a = a';end
if ~isrow(b);b = b';end



 m = size(a,2); n = size(b,2);
 [c,p] = sort([a,b]);
 q = 1:m+n; q(p) = q;
 t = cumsum(p>m);
 r = 1:n; r(t(q(m+1:m+n))) = r;
 s = t(q(1:m));
 id = r(max(s,1));
 iu = r(min(s+1,n));
 [d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
 ib = id+(it-1).*(iu-id);
 x_idx = ib;
end
