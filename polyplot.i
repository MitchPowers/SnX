#include "~/Tinkel/sn4.i"

fma;

yarr = span(0.0,1.0,41);

xcoords = yarr*0.0;
xcoords(1) = 0.0;

dx = span(-2.5,2.5,9);
pdist = 2^-(dx/3)^2;
pdist/= sum(pdist);
pcum = pdist(cum);

for(i=2; i<=numberof(xcoords); i++){
  ww = where(random() >= pcum)(0);
  delx = dx(ww);
  xcoords(i) = xcoords(i-1) + delx;
 }

xcoords -= xcoords(9);
ogcoord = xcoords;
xpath = savgol(xcoords,7,4);
ypath = yarr;
plotl,xpath(9:-8),ypath(9:-8),color="black",width=4;  


for(l=9; l<=numberof(yarr)-8; l+=4){
  gencol = reddish();
  xcoords = ogcoord;
  for(i=l; i<=numberof(xcoords); i++){
    ww = where(random() >= pcum)(0);
    delx = dx(ww);
    xcoords(i) = xcoords(i-1) + delx;
  }
  xpath(l:) = savgol(xcoords,7,4)(l:);
  ypath = yarr;
  plotl,xpath(9:-8),ypath(9:-8),color=gencol,width=4;  
 }
    
  
  
