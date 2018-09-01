K = 3 ;
data = uint8(rand(2,1000) * 255) ;
[C,A] = vl_ikmeans(data,K) ;
datat = uint8(rand(2,10000) * 255) ;
AT = vl_ikmeanspush(datat,C) ;
cl = get(gca,'ColorOrder') ;
ncl = size(cl,1) ;
for k=1:1
  sel  = find(A  == k) ;
  selt = find(AT == k) ;
  plot(data(1,sel),  data(2,sel),  '.',...
       'Color',cl(mod(k,ncl)+1,:)) ;
  hold on ;
  plot(datat(1,selt),datat(2,selt),'+',...
       'Color',cl(mod(k,ncl)+1,:)) ;
end