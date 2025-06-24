/******************************************************************

                      Coords      Data
                     At All        Block 
                    Times           Analysis

******************************************************************/

#include "~/TinkerExp/sn3.i"
extern NCOL
NCOL = 30;

func unwrap_traj(cat,dimarr){

  color=reddish();
  for(n=1; n<=dimsof(cat)(3); n++){//for all molecs
    trajx = cat(1,n,);
    trajy = cat(2,n,);
    trajz = cat(3,n,);
    if(n%20 == 1){
      if( (n/20)%3 == 1 ){
        color=blueish();
      }else if( (n/20)%3 == 2 ){
        color=reddish();
      }else{
        color=greenish();
      }
    }
    
    while(anyof(abs(trajx(dif))>20)){
      ww=where(abs(trajx(dif))>20);
      trajx(ww+1) -= DIM_ARR(1,ww)*sign(trajx(dif)(ww));
    }
    while(anyof(abs(trajy(dif))>20)){
      ww=where(abs(trajy(dif))>20);
      trajy(ww+1) -= DIM_ARR(2,ww)*sign(trajy(dif)(ww));
    }
    while(anyof(abs(trajz(dif))>20)){
      ww=where(abs(trajz(dif))>20);
      trajz(ww+1) -= DIM_ARR(3,ww)*sign(trajz(dif)(ww));
    }
    plotl,,trajz,color=color;
    pause,20;
    cat(3,n,) = trajz;
    cat(2,n,) = trajy;
    cat(1,n,) = trajx;
  }
  return cat;
  //cat analysis is performed by whisker.i
}

func pack_traj(cat,dimarr){
  moln = NATO/NMOL;
  
  if(anyof( abs(cat(1,,)) > DIM_ARR(1,-:1:moln,)/2.)){
    ww = where( abs(cat(1,,)) > DIM_ARR(1,-:1:moln,)/2.);
    tind = ww/moln;
    cat(1,ww) -= DIM_ARR(1,tind)*sign(cat(1,ww));
  }
  if(anyof( abs(cat(2,,)) > DIM_ARR(2,-:1:moln,)/2.)){
    ww = where( abs(cat(2,,)) > DIM_ARR(2,-:1:moln,)/2.);
    tind = ww/moln;
    cat(2,ww) -= DIM_ARR(2,tind)*sign(cat(2,ww));
  }
  return cat;
}
  

func ni_t(void,cat=){
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "CAT?";

  niz = cos(cat(4,,));
  nix = sin(cat(4,,))*cos(cat(5,,));
  niy = sin(cat(4,,))*sin(cat(5,,));

  nit = cat(:3,,);
  nit(1,,) = nix;
  nit(2,,) = niy;
  nit(3,,) = niz;
  
  return nit;
}

func mszd(void,cat=){
  /*DOCUMENT calcs msd along z-axis */
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "CAT?";

  zat = cat(3,,);

  TTT = numberof(zat(1,)); //length of time dimension
  dt_arr = array(0.,TTT);
  dt_arrs= array(0.,[2,NATO/NMOL,TTT]);
  
  dt_arr(1) = avg(zat(,dif)^2);
  dt_arrs(,1) = (zat(,dif)^2)(,avg);
  
  for(dt=2; dt<TTT; dt++){
    write,format="\b\b\b\b\b%.4g",int(10000.*dt/TTT)/100.;
    delt = zat(,::dt)(,dif);
    bign = numberof(delt);
    liln = numberof(delt(1,));
    delts= (delt^2)(,sum);
    delt = sum(delt^2);
    
    for(nt=2; nt<=min(dt,TTT-dt); nt++){
      dtmp = (zat(,nt::dt)(,dif));
      delt +=sum(dtmp^2);
      delts+=(dtmp^2)(,sum);
      bign += numberof(dtmp);
      liln += numberof(dtmp(1,));
    }
    dt_arr(dt) = delt/(bign-1.);
    dt_arrs(,dt)=delts/(liln - (liln>1));
  }
  write,"\b\b\b\b\b\b";
  dt_arr(0) = avg((zat(,0)-zat(,1))^2);
  dt_arrs(,0) = (zat(,0)-zat(,1))^2;
  
  for(i=1; i<=NATO/NMOL; i++){
    plotl,,dt_arrs(i,:-1),color=reddish();
  }
  
  q1 = int(NATO/NMOL * .25);
  q3 = int(NATO/NMOL * .75);
  out = array(0.,[2,3,TTT]);
  out(1,) = dt_arrs( sort(dt_arrs,1)(q1,) );
  out(2,) = dt_arr;
  out(3,) = dt_arrs( sort(dt_arrs,1)(q3,) );
  plotl,,out(1,:-1),color="grayc";
  plotl,,out(3,:-1),color="grayc";
  plotl,,out(2,:-1),color="black",width=4;

  //out is [Q1,avg,Q3] such that 50% of traces are between Q1 and Q3
  return out;
}


func mzd(void,cat=,plot_out=){
  /*DOCUMENT calcs msd along z-axis */
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "CAT?";

  zat = cat(3,,);

  TTT = numberof(zat(1,)); //length of time dimension
  dt_arrs= array(0.,[2,NATO/NMOL,TTT]);
  
  dt_arrs(,1) = zat(,dif)(,avg);

  for(dt=2; dt<TTT; dt++){
    write,format="\b\b\b\b\b%.4g",int(10000.*dt/TTT)/100.;
    delt = zat(,::dt)(,dif);
    liln = numberof(delt(1,));
    delts= (delt)(,sum);
    
    for(nt=2; nt<=min(dt,TTT-dt); nt++){
      dtmp  = (zat(,nt::dt)(,dif));
      delts+= (dtmp)(,sum);
      liln += numberof(dtmp(1,));
      
    }
    dt_arrs(,dt)=delts/(liln*1.);
  }
  write,"\b\b\b\b\b\b";
  
  for(i=1; i<=NATO/NMOL; i++){
    if(i%20==1) color=int(random(3)*255);
    plotl,,dt_arrs(i,),color=color;
    if(i%20==1) pause,100;
  }
  plotl,,abs(dt_arrs)(avg,),color="magenta",width=6;

  if(plot_out==1) speed_plot,dt_arrs,cat=cat;
  return dt_arrs;
}
func speed_plot(zat,cat=){
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "CAT?";

  tt = dimsof(zat)(0);
  slope_arr = avgx = avgy =array(0.,numberof(tt));
  for(i=1; i<=NATO/NMOL; i++){
    slope_arr(i) = regress(zat(i,),[tt^0,tt])(2);
    avgx(i) = cat(1,i,avg);
    avgy(i) = cat(2,i,avg);
  }
  slope_arr /= max(abs(slope_arr));
  for(i=1; i<=NATO/NMOL; i++){
    clr = slope_arr(i);
    if(clr>0){
      plot,avgx(i),avgy(i),color=int(clr*[255,0,0]);
    }else{
      plot,avgx(i),avgy(i),color=int(abs(clr*[0,0,255]));
    };
  }
}


func oriS_cat(void,cat=){

  nit = ni_t(cat=cat);
  AA = array(0.,[2,3,3]);
  num = NATO/NMOL;
  tim = dimsof(nit)(0);

  oriSt = array(0.,tim);

  n_z = array(0.,tim);
  
  for(t=1; t<=tim; t++){
    AA *= 0.;
    write,format="\b\b\b\b\b\b %g",int(10000*t/tim)/100.;
    for(n=1; n<=num; n++){
      ui=nit(,n,t);
      AA+=3*ui*ui(-:1:3,)-[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    }
    AA/=2.*num;
    eval=SVdec(AA,evec);
    oriSt(t) = max(eval);

    n_z(t) = abs(evec(3,where(eval==max(eval))));
  }
  write,"\b\b\b\b";
  write,format="n.z ~ \{%g,%g,%g\}\n",min(n_z),median(n_z),max(n_z);
  plotl,,n_z;
  return oriSt;
}



func gzt(void,cat=){
  /*DOCUMENT calcs G(z,t) for motion along z-axis  s-VdH??*/
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "CAT?";

  zat = cat(3,,);

  TTT = numberof(zat(1,)); //length of time dimension

  extern z_arr;
  extern t_arr;
  extern f_arr;
  z_arr = span(0.,10,251);
  t_arr = indgen(TTT);
  f_arr = z_arr(,-:1:TTT)*0.;
  f_tmp = f_arr(,1);

  fast_bin,abs(zat(,dif)),z_arr,f_tmp;
  f_arr(,1) = f_tmp/sum(f_tmp*1.);
  
  for(dt=2; dt<TTT; dt++){
    write,format="\b\b\b\b\b%.4g",int(10000.*dt/TTT)/100.;
    delt = zat(,::dt)(,dif);

    f_tmp*=0.;
    fast_bin,abs(delt),z_arr,f_tmp;
    
    for(nt=2; nt<=min(dt,TTT-dt); nt++){
      dtmp  = (zat(,nt::dt)(,dif));      
      fast_bin,abs(dtmp),z_arr,f_tmp;
    }
    f_arr(,dt) = f_tmp/sum(f_tmp*1.);
  }
  write,"\b\b\b\b\b\b";
  f_tmp *= 0.;
  fast_bin,abs(zat(,0)-zat(,1)),z_arr,f_tmp;
  f_arr(,0) = f_tmp/sum(f_tmp*1.);

  write,"Externals: z_arr, t_arr and f_arr";
  plfc,f_arr,z_arr(,-:1:TTT),t_arr(-:1:251,);
}



func cat_md(void,cat=,plot_out=){
  /*DOCUMENT calcs mean displacement of position cat params */
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "CAT?";

  TTT = numberof(cat(1,1,)); //length of time dimension
  ttt = TTT + TTT%2;
  l1 = indgen(ttt)(::-1);
  l2 = indgen(ttt/2); grow,l2,l2(::-1);
  est = (l1*l2)(cum)/sum(l1*l2*1.);
  dt_arrs= array(0.,[3,3,NATO/NMOL,TTT]);
  
  dt_arrs(,,1) = cat(:3,,dif)(,,avg);

  for(dt=2; dt<TTT; dt++){
    write,format="\b\b\b\b\b%.4g",int(100000.*est(dt))/1000.;
    delt = cat(:3,,::dt)(,,dif);
    liln = numberof(delt(1,1,));
    delts= (delt)(,,sum);
    
    for(nt=2; nt<=min(dt,TTT-dt); nt++){
      dtmp  = (cat(:3,,nt::dt)(,,dif));
      delts+= (dtmp)(,,sum);
      liln += numberof(dtmp(1,1,));
      
    }
    dt_arrs(,,dt)=delts/(liln*1.);
  }
  dt_arrs(,,0) = cat(:3,,0)-cat(:3,,1);
  
  write,"\b\b\b\b\b\b";
  
  for(i=1; i<=NATO/NMOL; i++){
    if(i%20==1) color=reddish();
    plotl,dt_arrs(1,i,),dt_arrs(2,i,),color=color;
    if(i%20==1) pause,100;
  }
  plotl,dt_arrs(1,avg,),dt_arrs(2,avg,),color="magenta",width=6;

  if(plot_out==1) speed_plot,dt_arrs,cat=cat;
  return dt_arrs;
}



func column_tracker(void,cat=,dim=){
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "You're gonna need a CAT";

  dim = is_void(dim)?DIM_ARR:dim;
  if(is_void(dim)) return "You're gonna need a DIM_ARR";

  if(dimsof(cat)(0)!=dimsof(dim)(0)) return "check your arrays";
     
  NUM = dimsof(cat)(0);
  MOL = int(NATO/NMOL);


  col_lists = array(0,[2,MOL,NUM]);
  colsat = array(0.,[2,2,NCOL]); //coords of the columns
  molsat = indgen(NCOL)(-:1:MPC,)(:MOL);  //indexes of columns
  //note that this assumes the first fram has naive ordering
  
  for(n=1; n<=NUM; n++){
    write,format="\b\b\b\b\b\b%g%%",int(10000.*n/NUM)/100.;
    //write,format="\b\b\b\b\b\b%d",n;
    cat_frame = cat(,,n);
    DIMS = DIM_ARR(:2,n);
    
    for(ci=1; ci<=NCOL; ci++){
      colvecs = cat_frame(,where(molsat==ci))(:2,);
      wx = where( abs(colvecs(1,) -colsat(1,ci))>DIM_ARR(1,n)/2.);
      wy = where( abs(colvecs(2,) -colsat(2,ci))>DIM_ARR(2,n)/2.);
      if(numberof(wx) ){
        colvecs(1,wx) -= DIM_ARR(1,n)*sign(colvecs(1,wx));
      }
      if(numberof(wy) ){
        colvecs(2,wy) -= DIM_ARR(2,n)*sign(colvecs(2,wy));
      }
      colsat(,ci) = colvecs(,avg);
    }
    for(i=1; i<=NATO/NMOL; i++){
      dist = norm(pbc(cat_frame(:2,i) - colsat));
      if(molsat(i)==0){//check if orphan has joined column
        molsat(i) = where(dist == min(dist))(1);
      }else if( min(dist)+1.25 < dist(molsat(i)) ){ 
        molsat(i) = where(dist == min(dist))(1); //closest
        molsat(i)*= min(dist)<4.5; //close enough
      }
    }
    orphans = where(abs(cat_frame(4,)-pi/2.)<=pi/6);
    if(numberof(orphans)) molsat(orphans) = 0;
    col_lists(,n) = molsat;
  }

  jumps = numberof( where( col_lists(,dif)));
  write,format="\n<Jump frequency> = %g nHz\n",1000.*jumps/(MOL*NUM*1.);
  return col_lists;
}

func col_sort(col_list,cat=){
  num = dimsof(col_list)(0); //number of times
  mol = dimsof(col_list)(2); //number of molecules
  cat = is_void(cat)?CAT:cat;
  oat = array(0.,[3,7,mol,num]);
  if(is_void(cat)) return "cat?";

  for(n=1; n<=num; n++){
    write,format="\b\b\b\b\b\b%g%%",int(1000.*n/num)/10.; 
    temp = cat(,,n)*0.;
    for(ci=1; ci<=mol/MPC; ci++){
      col = where(col_list(,n)==ci); //catable index of cith column
      srt = sort(cat(3,col,n));
      ind = where(temp(3,)==0)(1);
      temp(,ind:ind+numberof(col)-1) = cat(,col(srt),n);
    }
    orphans = where( col_list(,n)==0);
    ind = where(temp(3,)==0);
    if( numberof(orphans)>=1) temp(,ind(1):) = cat(,orphans,n);
    oat(:6,,n) = temp;
    if(numberof(orphans)>=1) oat(7,orphans,n) = 1; 
  }
  write,format="%s\n","\b\b\b\b\bZorted cat."
  return oat;
}
  


func zr_corr(col_list,cat=,dim=,plot_em=){
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "You're gonna need a CAT";

  dim = is_void(dim)?DIM_ARR:dim;
  if(is_void(dim)) return "You're gonna need a DIM_ARR";

  if(dimsof(cat)(0)!=dimsof(dim)(0)) return "check your arrays";

  tim = numberof(cat(3,1,));
  col = 30;
  num = dimsof(cat)(3);

  liln = int(sqrt(tim*(MPC-1)*MPC))+1;
  medn = int(sqrt(tim*col*MPC*col))+1;
  extern dz_arr;
  extern fz_arr;
  dz_arr = span(0.1,min(dim(3,))/2.,liln);
  fz_arr = dz_arr*0;
  extern dr_arr;
  extern fr_arr;
  dr_arr = span(0.1,min(dim(:2,))*sqrt(2)/2.,medn);
  fr_arr = dr_arr*0;

  if(!is_void(plot_em)){
    animate,1;
    fma;
  }

  zdim = dim(3,)(-:1:600,);
  ydim = dim(2,)(-:1:600,);
  xdim = dim(1,)(-:1:600,);
  
  for(i=1; i<=num; i++){
    write,format="\b\b\b\b%d",i;

    dzi = abs( cat(3,,) - cat(3,i,)(-:1:num,) );
    dri = abs( cat(:2,,) - cat(:2,i,)(,-:1:num,) );
    
    coli = col_list(i,)(-:1:num,);
    colat = where( (col_list == coli)*(coli!=0) );

    wz = where( dzi >= zdim/2.);
    if(numberof(wz)){
      dzi(wz) = zdim(wz) - dzi(wz);
    }
    wx = where( dri(1,) >= xdim/2.);
    if(numberof(wx)){
      dri(1,wx) = xdim(wx) - dri(1,wx);
    }
    wy = where( dri(2,) >= ydim/2.);
    if(numberof(wy)){
      dri(2,wy) = ydim(wy) - dri(2,wy);
    }
    dri = abs( dri(1,,) , dri(2,,) );
    taloc = where( (dri < dr_arr(0))*(dzi < 2.) );
    
    fast_bin,dzi(colat),dz_arr,fz_arr;
    fast_bin,dri(taloc),dr_arr,fr_arr;
    if(!is_void(plot_em)){
      plh,fz_arr/max(fz_arr*1.),dz_arr,color=greenish();
      plh,fr_arr/max(fr_arr*1.),dr_arr,color=blueish();
      fma;
    }
  }
  if(!is_void(plot_em)){
    fma;
    animate,0;
    fma;
  }
  
  gpara = (fz_arr(:-1)/sum(fz_arr*1.))*(2*MPC/dz_arr(dif)(1));
  zpara = dz_arr(zcen);

  //gperp = (fr_arr(:-1)/sum(fr_arr*1.))*(num*4/(pi*(dr_arr^2)(dif)));
  gperp = fr_arr(:-1)/sum(fr_arr*1.);
  zperp = dr_arr(zcen);
  /* cylinder sweep in a square...*/
  a = min(dim(2:,))/2.;
  cx = where(dr_arr>=a);
  ph = acos(a/dr_arr(cx));
  dang = 1 - ph(zcen)/pi;
  darr = (dr_arr(cx)^2)(dif);
  dA = (dr_arr^2)(dif);
  dA(cx(2:)-1) = dang*darr;
  gperp *= ((col-1)*dr_arr(0))/dA;
  
  plh,gpara,zpara;
  plh,gperp,zperp;
  fn1 = create("gpara.asc");
  write,fn1,"BinCenter normedP";
  write,fn1,format="%g %g\n",zpara,gpara;
  close,fn1;
  fn2 = create("gperp.asc");
  write,fn2,"BinCenter normedP";
  write,fn2,format="%g %g\n",zperp,gperp;
  close,fn2;
  
}
      
      
        
/*

  Recipes!

*/

func z_spacing(eqtime,col_list=,file_out=){
  //eqtime is the index to start collecting data from
  //assumes parse_traj and no more

  if(is_void(col_list)){
    write,"Tracking columns...";
    col_list = column_tracker();
  }
  zat = col_sort(col_list,cat=CAT);
  zat = zat(,,eqtime:);

  orphat = where(zat(7,,)==0);
  dz = zat(3,,)(orphat)(dif);
  dz = abs(dz);
  dz = dz(where(dz<10.));

  num = int(sqrt(numberof(dz)))+1;
  xz_arr = span(2.,6.,num);
  fz_arr = xz_arr*0;
  fast_bin,dz,xz_arr,fz_arr;

  zat = dz = orphat = col_list = [];

  extern prob;
  prob = fz_arr(:-1)/sum(1.*fz_arr);
  extern zarr;
  zarr = xz_arr(zcen);

  plh,prob,zarr,color="black";
  if(is_void(file_out)==0){
    file = create(file_out);
    write,file,"BinCenter Prob";
    write,file,format="%g %g\n",zarr,prob;
  }

  zexp = sum(zarr * prob);
  write,format="<dz> = %g\n",zexp;
  
  zsig = sqrt( sum(zarr^2 * prob) - zexp^2 );
  write,format="zsig = %g\n",zsig;

  HM = max(prob)/2.;
  FW = zarr(where( prob>HM ));
  FWHM = FW(0)-FW(1);
  write,format="FWHM = %g\n",FWHM;
  plotl,[FW(1),FW(0)],HM*[1,1];

  write,"Externs: zarr prob";

  return [zexp,zsig];
}
  
func cat_dph(eqtime,col_list=,file_out=){
  //eqtime is the index to start collecting data from
  //assumes parse_traj and no more

  if(is_void(col_list)){
    write,"Tracking columns...";
    col_list = column_tracker();
  }
  zat = col_sort(col_list,cat=CAT); //z sorted within cols
  zat = zat(,,eqtime:); //trimmed for time
  incol = where(zat(7,,)==0);

  dph = zat(6,incol)(dif);
  jmp = where(abs(zat(3,incol)(dif))<10.);
  dph = dph(jmp);

  signs = abs((zat(4,incol)>pi/2)(dif));
  signs = signs(jmp);
  dph(where(signs)) *= -1;

  if(anyof(abs(dph)>pi)){
    ww=where(abs(dph)>pi);
    dph(ww) -= 2*pi*sign(dph(ww));
  }

  extern x_dph;
  x_dph = span(-pi,pi,int(sqrt(numberof(dph)))+1);
  extern f_dph;
  f_dph = x_dph*0;
  fast_bin,dph,x_dph,f_dph;

  pdens = f_dph(:-1)/sum(f_dph*x_dph(dif)(1));
  pharr = x_dph(zcen);

  if(is_void(file_out)==0){
    file = create(file_out);
    write,file,"BinCenter Prob";
    write,file,format="%g %g\n",pdens,pharr;
  }

  return [];
}

  
  
