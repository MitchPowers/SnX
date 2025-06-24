/******************************************************************

                      Coords      Data
                     At All        Block 
                    Times           Analysis

******************************************************************/

#include "~/TinkerExp/sn3.i"
extern NCOL
NCOL = 30;
MPC = 20;
extern NATO
NATO = NCOL*MPC*NMOL;

if(is_void(MOLS)) MOLS=MPC*NCOL;


/*
CAT Processing
 */
func unwrap_traj(Cat,dim_arr){
  /* DOCUMENT unwrap_traj(cat,dimarr)
     de-periodicizes CAT
  */
  cat = Cat;
  color=reddish();
  for(n=1; n<=dimsof(cat)(3); n++){//for all molecs
    write,format="\r%d / %d",n,dimsof(cat)(3);
    trajx = cat(1,n,);
    trajy = cat(2,n,);
    trajz = cat(3,n,);
    
    traji = cat(7,n,);
    trajj = cat(8,n,);
    trajk = cat(9,n,);

    while(anyof(abs(trajx(dif))>20)){
      ww=where(abs(trajx(dif))>20)(1);
      write,format="\r%3d / %3d X %6d",n,dimsof(cat)(3),ww;
      trajx(ww+1:) -= dim_arr(1,ww)*sign(trajx(dif)(ww));
    }
    while(anyof(abs(trajy(dif))>20)){
      ww=where(abs(trajy(dif))>20)(1);
      write,format="\r%3d / %3d Y %6d",n,dimsof(cat)(3),ww;
      trajy(ww+1:) -= dim_arr(2,ww)*sign(trajy(dif)(ww));
    }
    while(anyof(abs(trajz(dif))>20)){
      ww=where(abs(trajz(dif))>20)(1);
      write,format="\r%3d / %3d Z %6d",n,dimsof(cat)(3),ww;
      trajz(ww+1:) -= dim_arr(3,ww)*sign(trajz(dif)(ww));
    }
    cat(3,n,) = trajz;
    cat(2,n,) = trajy;
    cat(1,n,) = trajx;

    while(anyof(abs(traji(dif))>20)){
      ww=where(abs(traji(dif))>20)(1);
      write,format="\r%3d / %3d I %6d",n,dimsof(cat)(3),ww;
      traji(ww+1:) -= dim_arr(1,ww)*sign(traji(dif)(ww));
    }
    while(anyof(abs(trajj(dif))>20)){
      ww=where(abs(trajj(dif))>20)(1);
      write,format="\r%3d / %3d Y %6d",n,dimsof(cat)(3),ww;
      trajj(ww+1:) -= dim_arr(2,ww)*sign(trajj(dif)(ww));
    }
    while(anyof(abs(trajk(dif))>20)){
      ww=where(abs(trajk(dif))>20)(1);
      write,format="\r%3d / %3d Z %6d",n,dimsof(cat)(3),ww;
      trajk(ww+1:) -= dim_arr(3,ww)*sign(trajk(dif)(ww));
    }
    cat(9,n,) = trajk;
    cat(8,n,) = trajj;
    cat(7,n,) = traji;
  }
  return cat;
}

func column_tracker(void,rcut=,cat=,dim=,ok=,nit=){
  /* Tracks column membership.
     NEEDS TO START WITH ORDERED COLUMNS
  */
  rcut = is_void(rcut)?6.0:rcut;
  
  cat = is_void(cat)?CAT:cat;
  if(is_void(cat)) return "You're gonna need a CAT";

  dim = is_void(dim)?DIM_ARR:dim;
  if(is_void(dim)) return "You're gonna need a DIM_ARR";

  if(dimsof(cat)(0)!=dimsof(dim)(0)) return "check your arrays";
     
  NUM = dimsof(cat)(0);
  MOL = dimsof(cat)(-1);

  col_lists = array(0,[2,MOL,NUM]);
  colsat = array(0.,[2,2,NCOL]); //coords of the columns
  molsat = indgen(NCOL)(-:1:MPC,)(:MOL);  //indexes of columns
  //note that this assumes the first frame has naive ordering
  if(ok==1){
    animate,1;
    fma;
    pause,500;
  }

  rmscol = array(0.,[2,NCOL,NUM]);
  popcol = array(0 ,[2,NCOL,NUM]);
  
  for(n=1; n<=NUM; n++){
    write,format="\b\b\b\b\b\b%g%%",int(10000.*n/NUM)/100.;
    //write,format="\b\b\b\b\b\b%d",n;
    cat_frame = cat(,,n);
    DIMS = dim(:2,n);
    
    for(ci=1; (ci<=NCOL)*(anyof(molsat)); ci++){
      //Loop over columns to find their average x-y coord
      mat = where(molsat==ci);
      if(numberof(mat)<=2){
        colsat(,ci) = [-999.,-999.];
        continue;
      }
      colvecs = cat_frame(,mat)(:2,);
      wx = where(abs(colvecs(1,)-colsat(1,ci))>dim(1,n)/2.);
      wy = where(abs(colvecs(2,)-colsat(2,ci))>dim(2,n)/2.);
      if(numberof(wx) ){
        colvecs(1,wx) -= dim(1,n)*sign(colvecs(1,wx));
      }
      if(numberof(wy) ){
        colvecs(2,wy) -= dim(2,n)*sign(colvecs(2,wy));
      }
      colsat(,ci) = colvecs(,avg);
      if(ok==1){
        plot,cat_frame(1,mat),cat_frame(2,mat),color="grayd";
        plot,colsat(1,ci),colsat(2,ci),color="red";
      } 
    }
    if( (ok==1)*(!allof(molsat)) ) plot,cat_frame(1,where(molsat==0)),cat_frame(2,where(molsat==0)),color="black";
    if(ok==1){
      plot,colsat(1,),colsat(2,),color="red";
      plotlc,colsat(1,),colsat(2,),rcut,color="grayc";
      //plotlc,colsat(1,),colsat(2,),3.25,color="black";
      limits,-dim(1,n)/2.,dim(1,n)/2.,-dim(2,n)/2.,dim(2,n)/2.;
      fma;
    }
    orphani = where(molsat==0);
    for(i=1; i<=MOLS; i++){
      //Loop over molecules to find the distance to the columns
      dist = pbc(cat_frame(:2,i) - colsat);
      dist = abs(dist(1,),dist(2,));

      //check if orphan has joined column
      if( (molsat(i)==0)*(min(dist)<rcut) ){
        molsat(i) = where(dist == min(dist))(1);
      }

      //check for adoptions
      if( where(dist==min(dist))(1) != molsat(i) ){ //i.e. near col != curr col
        if(min(dist)+1.0 < dist(molsat(i)) ){
          molsat(i) = where(dist == min(dist))(1);
        }
        else molsat(i) = 0;
        //so, when a mol is closer to a col than it's current col
        //if it's closer by an angstrom, it's adopted
        //otherwise it's orphaned
      }

      //check for wayward mols
      molsat(i) *= dist(molsat(i))<rcut;

      //check for verts
      molsat(i) *= abs(cat_frame(4,i)-pi/2) > pi/6;

      if(molsat(i)){
        rmscol(molsat(i),n) += dist(molsat(i))^2;
        popcol(molsat(i),n) += 1;
      }      
    }
        
    if(is_void(nit)!=1){
      for(ci=1; ci<=NCOL; ci++){
        colat = where(molsat==ci);
        if(numberof(colat) <= 1) continue;
        AA = array(0.,[2,3,3]);
        for(m=1; m<=numberof(colat); m++){
          ui = nit(,colat(m),n);
          AA+=3*ui*ui(-:1:3,)-[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
        }
        AA/=2.*numberof(colat);
        eval = SVdec(AA);
        if(!is_void(cS_arr)*(dimsof(cS_arr)(0)>=n)){
          cS_arr(ci,n) = max(eval);
        }
        if(max(eval)<=0.125){
          molsat(colat) = 0;
        }
      }
    }
    
    //if(numberof(orphans)>=1) molsat(orphans) = 0;
    col_lists(,n) = molsat;
    if( noneof(molsat) ){
      col_lists(,n:) = 0;
      break;
    }

  }
  rmscol(where(popcol))/= popcol(where(popcol));  
  rmscol = sqrt(rmscol(avg,avg) );
  pli,popcol;

  write,format="\navg col rms disp: %g",rmscol;
  if(rcut<rmscol*1.5){
    return column_tracker(rcut=rmscol*1.5,cat=cat,dim=dim,ok=ok,nit=nit);
  }

  jumps = numberof( where( col_lists(,dif)));
  write,format="\n<Jump frequency> = %g nHz\n",1000.*jumps/(MOL*NUM*1.);

  if(ok==1){
    animate,0;
    limits,e,e,e,e;
  }
  
  return col_lists;
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


/*
  Dynamics
  arguably the only things that need CATs
*/


func dynarun(dirs,lab=,maxT=,rat=){

  maxT = is_void(maxT)?3200:maxT;
  maxT = int(maxT);
  
  for(d=1; d<=numberof(dirs); d++){
    write,format="%s\n",dirs(d);
    qq = openb(dirs(d)+"cat.b");
    maxt = dimsof(qq.DOG)(0);
    eqt = int(qq.eqt);
    dog = qq.DOG(,,eqt:);
    //dog = qq.DOG(,,);
    close,qq;

    //dynacat,dog,lab=lab;
    dynacatALT,dog,maxT=maxT;

    file = is_void(lab)?"cat_dyn":"cat_dynLAB";
    qq=createb(dirs(d)+file+".b");
    save,qq,z2avg,r2avg,z2var,r2var,z4avg,r4avg,z4var,r4var;
    close,qq;
 
    dynout,dir=dirs(d),outdir=dirs(d),lab=lab;
  }

}

    
    
func dynacat(dog,lab=){
  //Direct calculation of
  // <r(t)>, <r(t)^2>, <r(t)^4>
  //for separations perp and para

  Tmax = dimsof(dog)(0);
  Nmol = MPC*NCOL;

  nit = is_void(lab)?ni_t(cat=dog):[0,0,1.](,-:1:Nmol,-:1:Tmax);
  
  Nwrt = array(0,Tmax-1);

  extern z2avg;  extern r2avg;  extern z2var;  extern r2var;
  extern z4avg;  extern r4avg;  extern z4var;  extern r4var;
  r2avg = r2var = z2avg = z2var = array(0.,[2,Nmol,Tmax-1]);
  r4avg = r4var = z4avg = z4var = array(0.,[2,Nmol,Tmax-1]);

  winkill,9;
  window,9,style="boxed.gs",wait=1;
  logxy,1,1;
  animate,1;
  
  for(t1=1; t1<Tmax; t1++){
    ti = 1./t1;
    write,format="\r%d/\%d ",t1,Tmax-1;
    dri = dog(:3,,t1) - dog(:3,,t1+1:);
    dzi = para2(dri, nit(,t1));
    dri-= dzi;

    Nwrt(:Tmax-t1)++;

    dr2 = (dri^2)(sum,,);
    dz2 = (dzi^2)(sum,,);
    dr4=dr2^2; dz4=dz2^2;
    
    // Welford rolling avg and var
    welf = dr2 - r2avg(,:Tmax-t1);
    r2avg(,:Tmax-t1) += welf * ti;
    welf*= dr2 - r2avg(,:Tmax-t1);
    r2var(,:Tmax-t1) += welf;//Ns2

    welf = dr4 - r4avg(,:Tmax-t1);
    r4avg(,:Tmax-t1) += welf * ti;
    welf*= dr4 - r4avg(,:Tmax-t1);
    r4var(,:Tmax-t1) += welf;//Ns2

    welf = dz2 - z2avg(,:Tmax-t1);
    z2avg(,:Tmax-t1) += welf * ti;
    welf*= dz2 - z2avg(,:Tmax-t1);
    z2var(,:Tmax-t1) += welf;//Ns2

    welf = dz4 - z4avg(,:Tmax-t1);
    z4avg(,:Tmax-t1) += welf * ti;
    welf*= dz4 - z4avg(,:Tmax-t1);
    z4var(,:Tmax-t1) += welf;//Ns2

    time = indgen(Tmax-1);

    plotl,time(:Tmax-t1),sqrt(r2avg(avg,:Tmax-t1)),color="grayc";
    plotl,time(Tmax-t1:),sqrt(r2avg(avg,Tmax-t1:)),color="blue",width=4;
    plotl,time(:Tmax-t1),sqrt(z2avg(avg,:Tmax-t1)),color="grayc";
    plotl,time(Tmax-t1:),sqrt(z2avg(avg,Tmax-t1:)),color="red",width=4;
    plxtitle,"Time displacement";
    plytitle,"rmsd";
    fma;   
    
  }
  
  animate,0;

  Nwrt = (1./Nwrt)(-:1:Nmol,);
  r2var *=  Nwrt;
  r4var *=  Nwrt;
  z2var *=  Nwrt;
  z4var *=  Nwrt;
}


func dynacatALT(dog,maxT=){
  //Direct calculation of
  // <r(t)>, <r(t)^2>, <r(t)^4>
  //for separations perp and para

  Tmax = dimsof(dog)(0);
  maxT = is_void(maxT)?3200:maxT;
  maxT = min(maxT,Tmax-1);
  Nmol = MPC*NCOL;

  Nwrt = array(0,maxT);

  extern z2avg;  extern r2avg;  extern z2var;  extern r2var;
  extern z4avg;  extern r4avg;  extern z4var;  extern r4var;
  r2avg = r2var = z2avg = z2var = array(0.,[2,Nmol,maxT]);
  r4avg = r4var = z4avg = z4var = array(0.,[2,Nmol,maxT]);

  time = indgen(maxT);
      
  winkill,9;
  window,9,style="boxed.gs",wait=1;
  logxy,1,1;
  animate,1;
  
  for(t1=1; t1<Tmax; t1++){
    ti = 1./t1;
    write,format="\r%d/\%d ",t1,Tmax-1;
    dri = dog(:3,,t1) - dog(:3,,t1+1:);
    if(dimsof(dri)(0)>maxT) dri = dri(,,:maxT);
    lent = dimsof(dri)(0);
    dzi = dri(3,,);
    dri = dri(:2,,);

    Nwrt(:lent)++;

    dr2 = (dri^2)(sum,,);
    dz2 = (dzi^2);
    dr4=dr2^2; dz4=dz2^2;
    
    // Welford rolling avg and var
    welf = dr2 - r2avg(,:lent);
    r2avg(,:lent) += welf * ti;
    welf*= dr2 - r2avg(,:lent);
    r2var(,:lent) += welf;//Ns2

    welf = dr4 - r4avg(,:lent);
    r4avg(,:lent) += welf * ti;
    welf*= dr4 - r4avg(,:lent);
    r4var(,:lent) += welf;//Ns2

    welf = dz2 - z2avg(,:lent);
    z2avg(,:lent) += welf * ti;
    welf*= dz2 - z2avg(,:lent);
    z2var(,:lent) += welf;//Ns2

    welf = dz4 - z4avg(,:lent);
    z4avg(,:lent) += welf * ti;
    welf*= dz4 - z4avg(,:lent);
    z4var(,:lent) += welf;//Ns2


    plotl,,sqrt(r2avg(avg,)),color="grayc";
    plotl,time(lent:),sqrt(r2avg(avg,lent:)),color="blue",width=4;
    plotl,,sqrt(z2avg(avg,)),color="grayc";
    plotl,time(lent:),sqrt(z2avg(avg,lent:)),color="red",width=4;
    plotl,[1,3200],sqrt(5.6)*[1,1],color="black";
    plxtitle,"Time displacement";
    plytitle,"rmsd";
    fma;
    
  }
  
  animate,0;

  Nwrt = (1./Nwrt)(-:1:Nmol,);
  r2var *=  Nwrt;
  r4var *=  Nwrt;
  z2var *=  Nwrt;
  z4var *=  Nwrt;
}




func dynout(void,dir=,outdir=,lab=){
  dir = is_void(dir)?"./":dir;
  outdir = is_void(outdir)?dir:outdir;

  file = is_void(lab)?"cat_dyn":"cat_dynLAB";
  file+= ".b";
  dat = openb(dir+file);
  restore,dat;

  z2 = z2avg(avg,:-100);
  z4 = z4avg(avg,:-100);
  za = z4/(3*z2^2) - 1.;
  r2 = r2avg(avg,:-100);
  r4 = r4avg(avg,:-100);
  ra = r4/(2*r2^2) - 1.;

  if(!is_void(lab)) dat = create(outdir+"msdoutLAB.asc");
  else              dat = create(outdir+"msdout.asc");
  write,dat,format="%s %s %s %s %s\n","Time","z2","z_ngp","r2","r_ngp";
  write,dat,format="%d %g %g %g %g\n",indgen(numberof(z2)),z2,za,r2,ra;
  close,dat;

  window,0,wait=1;
  logxy,1,1;
  plotl,,z2,color="red";
  plotl,,r2,color="black";

}
    
  
  
    
    
  







func neighbors(ref,tind,z_cut=,r_cut=){
  /*DOCUMENT returns index of molecules near the ref molec
    uses a cylinder of height z_cut= and radius r_cut=
    defaults to 6x16
  */

  z_cut = is_void(z_cut)?6.0:z_cut;
  r_cut = is_void(r_cut)?16.0:r_cut;
  com = cat(1:3,,tind);
  dir = cat(4:6,,tind);
  drr = pbc(com - com(,ref)); //separation vectors w/ enforced pbc
  inz = where(norm(para2(drr,dIR(,ref)))<=z_cut/2.);
  inr = where(norm(perp2(drr,dIR(,ref)))<=r_cut);
  uni = union(inz,inr);
  uni = uni(where(uni!=ref));
  return int(uni);
}



func oriS_cat(void,cat=){

  nit = ni_t(cat=cat);
  AA = array(0.,[2,3,3]);
  num = dimsof(cat)(2);
  tim = dimsof(nit)(0);

  oriSt = array(0.,tim);

  n_z = array(0.,tim);
  
  for(t=1; t<=tim; t++){
    AA *= 0.;
    write,format="\r %g",int(10000*t/tim)/100.;
    for(n=1; n<=num; n++){
      ui=nit(,n,t);
      AA+=3*ui*ui(-:1:3,)-[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    }
    AA/=2.*num;
    eval=SVdec(AA,evec);
    oriSt(t) = max(eval);

    n_z(t) = abs(evec(3,where(eval==max(eval))));
  }
  write,format="%s","\r";
  //write,format="n.z ~ \{%g,%g,%g\}\n",min(n_z),median(n_z),max(n_z);
  return oriSt;
}

func herd_md(dirs){
  if( anyof(strpart(dirs,0:)!="/") ){
    ww= where(strpart(dirs,0:)!="/");
    dirs(ww) += "/";
  }
  for(d=1; d<=numberof(dirs); d++){
    write,format=" Opening %s...\n",dirs(d);
    qq=openb(dirs(d)+"cat.b");
    dog = unwrap_traj(qq.CAT(,,qq.eqt:),qq.DIM_ARR(,qq.eqt:));
    close,qq;
    mds = cat_md(cat=dog);
    qq=updateb(dirs(d)+"cat_md_redog.b");
    save,qq,mds;
    close,qq;
  }
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
  dt_arrs= array(0.,[3,6,MPC*NCOL,TTT]);
  
  dt_arrs(,,1) = cat(:6,,dif)(,,avg);

  for(dt=2; dt<TTT; dt++){
    write,format="\r%.4g",int(100000.*est(dt))/1000.;
    delt = cat(:6,,::dt)(,,dif);
    liln = numberof(delt(1,1,));
    delts= (delt)(,,sum);
    
    for(nt=2; nt<=min(dt,TTT-dt); nt++){
      dtmp  = (cat(:6,,nt::dt)(,,dif));
      delts+= (dtmp)(,,sum);
      liln += numberof(dtmp(1,1,));      
    }
    dt_arrs(,,dt)=delts/(liln*1.);
  }
  dt_arrs(,,0) = cat(:6,,0)-cat(:6,,1);
  
  write,"\b\b\b\b\b\b";
  
  for(i=1; i<=NCOL*MPC; i++){
    if(i%MPC==1) color=reddish();
    plotl,dt_arrs(1,i,),dt_arrs(2,i,),color=color;
    if(i%MPC==1) pause,100;
  }
  plotl,dt_arrs(1,avg,),dt_arrs(2,avg,),color="magenta",width=6;

  if(plot_out==1) speed_plot,dt_arrs,cat=cat;
  return dt_arrs;
}

func speed_plot(dt_arr,cat){
  trimat = int(sqrt(10)*1000);
  trimat = min(trimat,dimsof(cat)(0));
  cat = cat(,,:trimat);
  mds = dt_arr(,,:trimat);
  tim = indgen(trimat)*1.;

  zlmfs = array(0.,[2,2,30]);
  colat = array(0.,[2,2,30]);

  for(ci=1; ci<=30; ci++){
    mems = where( (cat(10,,)==ci)(,avg)>.8 );
    mdsc = mds(3,mems,)(avg,);
    zlmfs(,ci) = regress( mdsc , [tim^0,tim] );
    colat(,ci) = cat(:2,mems,avg)(,avg);

    zcliq = round( (zlmfs(2,ci)*3160)/1.9 ); //approx number of half steps traversed by the column
    label = "";

    if(zcliq==0){
      label+="~";
      plotlc,colat(1,ci),colat(2,ci),2,width=4,color="grayc";
      plt,label,colat(1,ci),colat(2,ci),tosys=1,justify="CH";
    }
    if(zcliq >0){
      label+="+"; label+=totxt(int(zcliq));
      plotlc,colat(1,ci),colat(2,ci),2,width=4,color="red";
      plt,label,colat(1,ci),colat(2,ci),justify="CH",tosys=1;
    }
    if(zcliq< 0){
      label+=totxt(int(zcliq));
      plotlc,colat(1,ci),colat(2,ci),2,width=4,color="blue";
      plt,label,colat(1,ci),colat(2,ci),justify="CH",tosys=1;
    }
  }
}
  

  

func ori_stuff(dirs){
  for(i=1; i<=numberof(dirs); i++){
    fn = dirs(i)+"cat.b";
    qq=openb(fn);
    restore,qq;
    write,format="%d / %d\n",eqtime,dimsof(CAT)(0);
    Sarr = oriS_cat(cat=CAT(,,eqtime:));
    xarr = span(max(0.,min(Sarr)-.01),min(1.,max(Sarr)+.01), int(sqrt(numberof(Sarr))+1) );
    farr = xarr*0;
    bb = fast_bin(Sarr,xarr,farr);
    window,1,wait=1;
    plh,bb(,2)/sum(bb(,2)*bb(dif,1)(1)),bb(,1),color=["black","green","blue","yellow","red"](i);
    distro_dat,farr,xarr;
    window,0,wait=1;
    out = create("oriS_tmp"+totxt(i)+".asc");
    write,out,"Sbin Pdens";
    write,out,format="%g %g\n",bb(zcen,1),bb(:-1,2)/sum(bb(,2)*bb(dif,1)(1));
    close,out;
  }
}
    

func col_oriS(cat){
  nit = ni_t(cat=cat);
  ncol = 30;
  Sout = array(0.,[2,ncol,dimsof(cat)(0)]);
  for(t=1; t<=dimsof(cat)(0); t++){
    for(ci=1; ci<=ncol; ci++){
      colat = where(cat(10,,t)==ci);
      if(numberof(colat) <= 1) continue;
      AA = array(0.,[2,3,3]);
      for(n=1; n<=numberof(colat); n++){
        ui = nit(,colat(n),t);
        AA+=3*ui*ui(-:1:3,)-[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
      }
      AA/=2.*numberof(colat);
      eval = SVdec(AA);
      Sout(ci,t) = max(eval);
    }
  }
  return Sout;
}



func pair_or(ind,jnd,tind){
  //USES CAT//
  //cm's are geometric centers based on the given KERN
  //use check_kern to verify
  jnd = jnd(where(jnd!=ind)); //can't orient to itself
  if(!numberof(jnd)) return "no pair!";
  if(is_void(CAT) ) return "no cat!";
  rc = CAT(7:9,,tind);
  or = CAT(6,,tind);
  nz = sign(CAT(4,,tind)-pi/2.); //+:"up" -:"down"

  
  //angle between vec R_12 and labframe x
  Rij= pbc( rc(,jnd)-rc(,ind));
  dr = abs( Rij(1,), Rij(2,) );
  orij = acos(Rij(1,)/dr) * sign(Rij(2,));

  //thi is the orientation of ind facing jnd
  thi = (or(ind) - orij)*sign(nz(ind));
  if(anyof(abs(thi)>pi)){
    ww=where(abs(thi)>pi);
    thi(ww) -= 2*pi*sign(thi(ww));
  }
  
  //thj is the orientation of jnd facing ind
  thj = (pi + or(jnd) - orij)*sign(nz(jnd));
  if(anyof(abs(thj)>pi)){
    ww=where(abs(thj)>pi);
    thj(ww) -= 2*pi*sign(thj(ww));
  }

  return [thi,thj,dr];
}

func autotranVH(dir){
  write,format="\r%s","Opening CAT and getting DOG";
  qq=updateb(dir+"cat.b");
  /*
    DOG = qq.CAT;
  DOG = unwrap_traj(DOG,qq.DIM_ARR);
  save,qq,DOG;
  */
  dog = qq.DOG(,,int(qq.eqt):);
  close,qq;
  write,format="\r%s\n","VHing "+dir;
  
  //sVDH,dog,cutat=17.5;
  sVDHalt,dog;
  //qq = createb(dir+"cat_vdhLAB.b");
  qq = createb(dir+"cat_vdhLABalt.b");
  save,qq,dr_arr,dz_arr;
  save,qq,fr_arr,fz_arr;
  save,qq,pzden,prden;
  close,qq;
}
func sVDH(dog,cutat=){
  if(is_void(MOLS)) MOLS=MPC*NCOL;

  dog = dog(7:9,,); //COR!

  dTmax = 3201;
  dTmax = min(dTmax,dimsof(dog)(0));
  Tmax = dimsof(dog)(0);
  numB = int( (Tmax*MOLS)^(1/2.) ) + 1;
  numB = min(numB,int(300*8.5/2.54));

  
  zmax = rmax = 0.0;
  cutat = is_void(cutat)?17.5:cutat;
  write,dimsof(dog);
  zmax = max( dog(3,,max) - dog(3,,min) ) + 0.01;
  zmax = min(cutat,zmax)+0.01;
  rmax = cutat+0.01;
    //for(t1=1; t1<Tmax; t1++){
      //write,format="\r%d/\%d ",t1,Tmax;
      //dri = dog(:3,,t1) - dog(:3,,t1+1:);
      //dzi = para2(dri, nit(,t1));
      //dzi = abs(dri(3,,));
      //dri = dri - dzi;
      //dzi = abs(dzi(1,,),dzi(2,,),dzi(3,,));
      //dri = abs(dri(1,,),dri(2,,));
      //zmax = max(zmax,max(dzi));
      //rmax = max(rmax,max(dri));
      //t1 += int(random()*10);
    //}
  extern dz_arr;
  extern dr_arr;
  dz_arr = span(0.,(cutat)^(.8),numB)^(1.25);  //NONLINEAR SPACING
  dr_arr = span(0.,(cutat)^(.8),numB)^(1.25);  //NONLINEAR SPACING
  grow,dz_arr,zmax;
  grow,dr_arr,rmax;
  extern fz_arr;
  extern fr_arr;
  fz_arr = dz_arr(,-:1:dTmax) * 0;
  fr_arr = dr_arr(,-:1:dTmax) * 0;
  extern pzden;
  extern prden;
  pzden = fz_arr*0.0;
  prden = fr_arr*0.0;

  shuf = indgen(dTmax);

  numz = numberof(dz_arr);
  numr = numberof(dr_arr);
  animate,1;
  for(dt0=1; dt0<dTmax; dt0++){
    dt = shuf(dt0);
    write,format="\r%d/\%d ",dt0,dTmax;

    fzt = fz_arr(,dt);
    frt = fr_arr(,dt);
    cap = min(dt,Tmax-dt);
    for(t0=1; t0<=cap; t0++){
      dri = dog(,,t0::dt)(,,dif);
      dzi = abs( dri(3,,) );
      dri = abs( dri(1,,), dri(2,,) );
      fast_bin,dzi,dz_arr,fzt;
      fast_bin,dri,dr_arr,frt;
    }
    fr_arr(,dt)+= frt;
    fz_arr(,dt)+= fzt;
    frt*=0;
    fzt*=0;

    Plotl,dr_arr(zcen),fr_arr(:-1,dt)/(dr_arr^2)(dif),color="black";
    
  }
  animate,0;

  dz = dz_arr(dif);
  zsums = fz_arr(sum,)(-:1:numz-1,);
  ww = where(zsums);
  pzden = fz_arr(:-1,)/dz;
  pzden(ww)/= zsums(ww);

  dr = (dr_arr^2)(dif);
  rsums = fr_arr(sum,)(-:1:numr-1,);
  ww = where(rsums);
  prden = fr_arr(:-1,)/(dr*pi);
  prden(ww)/= rsums(ww);
}
  
func sVDHagate(dirs,outdir){
  outl=strpart(dirs,9:15);
  for(d=1; d<=numberof(dirs); d++){
    write,dirs(d);
    qq=updateb(dirs(d)+"cat_vdhLAB.b");
    restore,qq;
    sVDH_gauss,outdir=outdir;//,labl=outl(d);
    save,qq,pzden,prden;
  }
  write,"\n";
}


func sVDH_gauss(void,outdir=,labl=,no_plot=){
  //uses externals currently in memory
  labl = is_void(labl)?totxt(int(random()*10^3)):labl;

  extern pzden;
  extern prden;
  write,format=" %s\n","External Var: pzden","External Var: prden";
  ztots = fz_arr(sum,)*1.;
  trimt = (ztots==0);
  trimt = anyof(trimt)?where(trimt)(1)-1:dimsof(fr_arr)(0)-1;
  write,trimt;
  pzden = fz_arr(,:trimt);
  pzden/= ztots(-:1:dimsof(pzden)(2),:trimt);
  pzden(2:)/= dz_arr(dif);
  pzden(1,)/= dz_arr(1);
  sig_z = (fz_arr(sum,:trimt));

  rtots = fr_arr(sum,)*1.;
  trimt = (rtots==0);
  trimt = anyof(trimt)?where(trimt)(1)-1:dimsof(fr_arr)(0)-1;
  write,trimt;
  prden = fr_arr(,:trimt);
  prden/= rtots(-:1:dimsof(prden)(2),:trimt);
  dA    = pi*dr_arr(1)^2;
  grow,dA,pi*(dr_arr^2)(dif);
  prden/= dA;
  sig_r = (fr_arr(sum,:trimt));

  if(is_void(no_plot)){
    window,0,wait=1;
    plotl,,sig_z;
    plotl,,sig_r;
  }
  
  zcens = dz_arr(1)/2.; grow,zcens,dz_arr(zcen);
  rcens = dr_arr(1)/2.; grow,rcens,dr_arr(zcen);
  zcuts = min(pzden(where(fz_arr(,:trimt)))/fz_arr(where(fz_arr(,:trimt))));
  rcuts = min(prden(where(fr_arr(,:trimt)))/fr_arr(where(fr_arr(,:trimt))));

  outfit = array(0.,[2,4,trimt]);

  //fit to gaussian
  //Y   = A exp(-B*x^2 )
  //logY = log(A) - B x^2
  //logy = log(A) - B x^2
  oom = 0;
  
  for(i=1; i<=dimsof(outfit)(0); i++){
    ww = where(pzden(,i)>zcuts);
    zlmf = regress(log(pzden(ww,i)),[zcens(ww)^0.,zcens(ww)^2],sigy=(pzden(ww,i)/max(pzden(ww,i))/sig_z(ww))^-.3,ax,va,ch);

    ww = where(prden(,i)>rcuts);
    rlmf = regress(log(prden(ww,i)),[rcens(ww)^0.,rcens(ww)^2],sigy=(prden(ww,i)/max(prden(ww,i))/sig_r(ww))^-.3,ax,va,ch);

    outfit(2,i) =zb= -1.*zlmf(2); //b
    outfit(1,i) =za= exp( zlmf(1) ); //a

    outfit(4,i) =rb= -1.*rlmf(2); //b
    outfit(3,i) =ra= exp( rlmf(1) ); //a

    if(i%int(10^oom)==0){
      if(!is_void(outdir)){
        out = create(outdir+labl+totxt(i)+"_Pdens.asc");
        write,out,format="%s %s %s %s %s %s\n","Z(A)","Pzdens","Fit","R(A)","PrdensX2pi","Fit";
        write,out,format="%g %g %g %g %g %g\n",zcens,pzden(,i),za*exp(-zb*(zcens)^2),rcens,prden(,i)*2*pi,ra*exp(-rb*(rcens)^2)*2*pi;
        close,out;
      }
      oom+=0.5; 
      if(is_void(no_plot)){
        window,7,wait=1;
        ww = where(pzden(,i)>zcuts);
        maxr = max(pzden(ww,i));
        plotl,-zcens(ww),pzden(ww,i)/maxr,color="black";
        plotl,-zcens(ww),za*exp(-zb*(zcens(ww))^2)/maxr,color="blue";

        ww = where(prden(,i)>rcuts);
        maxr = max(prden(ww,i));
        plotl,rcens(ww),prden(ww,i)/maxr,color="grayd";
        plotl,rcens(ww),ra*exp(-rb*(rcens(ww))^2)/maxr,color="red";
      }
    }
  }
  return 0; 
  zW = 1./sqrt(outfit(2,));
  rW = 1./sqrt(outfit(4,));
  tt = indgen(numberof(zW));

  lzf=regress(zW,[tt^0.,sqrt(tt)],sigy=fz_arr(sum,:trimt)^-.75); zD=(lzf(0)*0.5)^2;
  lrf=regress(rW,[tt^0.,sqrt(tt)],sigy=fr_arr(sum,:trimt)^-.75); rD=(lrf(0)*0.5)^2;
  write,format="D_z = %g AA/s\n",zD;
  write,format="LMF Z: (%g,%g) to (%g,%g)\n",1.,lzf(1)+lzf(2),60.,lzf(1)+60*lzf(2);
  write,format="D_r = %g AA/s\n",rD;
  write,format="LMF R: (%g,%g) to (%g,%g)\n",1.,lrf(1)+lrf(2),60.,lrf(1)+60*lrf(2);  

  if(is_void(no_plot)){
    limits,e,e,cutat,e;
    window,0,wait=1;
  
    window,window()+1,wait=1;
    plotl,sqrt(tt),zW,color="cyan",width=4;
    plotl,sqrt(tt),lzf(1)+2*sqrt(zD*tt),color="blue";
    plotl,sqrt(tt),rW,color="magenta",width=4;
    plotl,sqrt(tt),lrf(1)+2*sqrt(rD*tt),color="red";
  }

  if(!is_void(outdir)){
    out = create(outdir+labl+"GaussWidths.asc");
    write,out,format="%s %s %s\n","Time(ps)","Wz(Ang)","Wr(Ang)";
    write,out,format="%g %g %g\n",tt*1.,zW,rW;
  }
  
  return outfit;
}

func sVDH_uhlen(void,outdir=,trim=){
  trim = is_void(trim)?0:trim;
  //uses externals currently in memory
  extern pzden;
  extern prden;
  write,format=" %s\n","External Var: pzden","External Var: prden";
  fz_arr= fz_arr(,where(fz_arr(sum,)));
  pzden = (fz_arr/(fz_arr(sum,)(-:1:dimsof(fz_arr)(2),)*1.));
  pzden*= 1./dz_arr(dif)(1);
  fr_arr= fr_arr(,where(fr_arr(sum,)));
  prden = (fr_arr/(fr_arr(sum,)(-:1:dimsof(fr_arr)(2),)*1.));
  dA    = dr_arr(1)^2; grow,dA,(dr_arr^2)(dif);
  prden*= 1./dA;
  rcuts = 1./sqrt(sum(fr_arr));

  zcens = dz_arr(1)/2.; grow,zcens,dz_arr(zcen);
  rcens = dr_arr(1)/2.; grow,rcens,dr_arr(zcen);
  zcuts = 1./sqrt(sum(fz_arr));

  if(trim){
    pzden = pzden(,:trim);
    prden = prden(,:trim);
  }

  outfit = array(0.,[2,6,dimsof(pzden)(0)]);

  //fit to gaussian
  //  Y  = A exp(-B*(x-C)^2 )
  //logY = log(A) - B(C^2 + x^2) + 2BC x 
  //logy = (log(A)-BC^2) + 2BC x - B x^2
  //  B == -lmf(3)
  //  C == lmf(2)/(2*B)
  //  A = exp(lmf(1)+C*B*C)
  //sig = sqrt(1./2*B)
  oom = 0;
  
  for(i=1; i<=dimsof(outfit)(0); i++){
    ww = where(pzden(,i)>zcuts);
    zlmf = regress(log(pzden(ww,i)),[zcens(ww)^0.,zcens(ww),zcens(ww)^2],sigy=1./pzden(ww,i)^1.75,ax,va,ch);

  
    ww = where(prden(,i)>rcuts);
    rlmf = regress(log(prden(ww,i)),[rcens(ww)^0.,rcens(ww),rcens(ww)^2],sigy=1./prden(ww,i)^1.75,ax,va,ch);

    outfit(2,i) =zb= -1.*zlmf(3); 
    outfit(3,i) =zc= zlmf(2)/(2*zb);
    outfit(1,i) =za= exp( zlmf(1) - zb*zc*zc );

    
    outfit(5,i) =rb= -1.*rlmf(3);
    outfit(6,i) =rc= rlmf(2)/(2*rb);
    outfit(4,i) =ra= exp( rlmf(1) - rb*rc*rc );

    if((outfit(2,i)<=0)+(outfit(5,i)<=0) ){
      write,"trim before "+totxt(i);
    }
    outfit(2,i) = 1./sqrt(outfit(2,i)); //W
    outfit(5,i) = 1./sqrt(outfit(5,i)); //W
  
    if(i%10^oom==0){
      if(i==1) fma;
      window,7,wait=1;
      oom++;
      //plotl,-zcens,pzden(,i),color="black";
      plotl,-zcens,pzden(,i) - za*exp(-zb*(zcens-zc)^2),color="red";

      //plotl,rcens,prden(,i),color="grayd";
      plotl,rcens,prden(,i) - ra*exp(-rb*(rcens-rc)^2),color="blue";
    }
  }
  window,0,wait=1;
  plotl,,outfit(3,),color=blueish();
  plotl,,outfit(6,),color=reddish();
  pltitle,"M - z is blue";
  window,1,wait=1;
  plotl,,outfit(2,),color=blueish();
  plotl,,outfit(5,),color=reddish();
  pltitle,"W - z is blue";


  tt = indgen(numberof(outfit(1,)))*1.;

   if(!is_void(outdir)){
    out = create(outdir);
    write,out,format="%s %s %s %s %s\n","Time(ps)","Wz(Ang)","Wr(Ang)","Mz(ang)","Mr(ang)";
    write,out,format="%g %g %g %g %g\n",tt*1.,outfit(2,),outfit(5,),outfit(3,),outfit(6,);
  }
   
  return outfit;
}




func wedgie(xyz_mol_tim){
  animate,1;
  t1=t2=[0.,0.,0.];
  for(t=1; t<=dimsof(xyz_mol_tim)(0); t++){
    x_tim = xyz_mol_tim(1,,t);
    y_tim = xyz_mol_tim(2,,t);
    z_tim = xyz_mol_tim(3,,t);
    ro_tim= abs(x_tim,y_tim);
    xp_tim= ro_tim*cos(acos(x_tim/ro_tim)%(pi/3));
    yp_tim= ro_tim*sin(acos(x_tim/ro_tim)%(pi/3));

    plot,xp_tim,yp_tim,color="red";
    while( (t2-t1)(0)<0.01){
      timer,t2;
    }
    fma;
    if(t==dimsof(xyz_mol_tim)(0)) t=1;
  }
  animate,0;
}


func orphan_census(cat,color=){
  write,format="Sim Length: %d",dimsof(cat)(0);
  melt = cat(10,sum,)==0;
  if(anyof(melt)) melt = where(melt)(1);
  else melt = numberof(melt);
  cat = cat(,,:melt);
  write,format="\nNot Melted: %d",dimsof(cat)(0);
  
  oavg = cat(10,1,);
  time = indgen(numberof(oavg));
  osum=(cat(10,,)==0)(sum,);

  Nlib = 0;
  trap = indgen(200);
  frap = trap*0;
  lrap = log(trap);
  for(i=1; i<=600; i++){
    qorf = anyof(cat(10,i,)==0);
    Nlib += qorf;
    out_in = where((cat(10,i,)==0)(dif));
    qorf = numberof(out_in)>=2;
    if(qorf){
      fast_bin, out_in(dif)(::2),trap,frap;
    } 
  }
  write,format="\n %d of %d were orphaned.",Nlib,dimsof(cat)(3);
  if(sum(frap)>2){
    yy = frap/sum(frap*1.);
    ww = where(frap);
    lmf = regress(log(yy(ww)),[lrap(ww)^0.,lrap(ww)],sigy=trap(ww)^.5);
    nexp = lmf(2);
    texp = exp(lmf(1)/(-1.*nexp)); 
    write,format="\n Trap time exponent, n = %.3g",nexp;
    write,format="\n                     To= %.3g",texp
    
    window,2,wait=1;
    plot,trap,yy,color=color;
    plotl,trap,exp(lmf(1))*trap^nexp;
  }

    window,1,wait=1;
    plotl,time*.001,osum,color=color;
    
 
  return osum;
}



func orphan_plotter(eqt){

  fma; 

  zoff = max(DIM_ARR(3,eqt:))/2. + max(DIM_ARR(2,eqt:))/2. + 5;

  reds = int(random([2,3,600])*256);
  reds(1,)*=0;
  
  for(t=eqt+2; t<dimsof(CAT)(0); t++){
    DIMS=DIM_ARR(,t);
    cat=pbc(CAT(,,t));
    ww=where(CAT(10,,t));//not orphans
    if(numberof(ww)){
      plot,cat(1,ww),cat(2,ww),color=(int((t-eqt)*200./(dimsof(CAT)(0)-eqt))+25)(-:1:3),sym='\1';
      plot,cat(1,ww),cat(3,ww)-zoff,color=(int((t-eqt)*200./(dimsof(CAT)(0)-eqt))+25)(-:1:3),sym='\1';
      pause,1;
    }
  }
  lt=4; rt=6;
  //lt=rt=dt=3;
  
  for(t=eqt+lt+2; t<dimsof(CAT)(0)-rt; t++){
    ww = where(CAT(10,,t)==0);//orphans
    if(numberof(ww)){
      for(n=1; n<=numberof(ww); n++){
        for(tt=lt*-1; tt<=rt; tt++){
          tail=pbc( CAT(:3,ww(n),t+tt-1:t+tt) );
          if( abs(tail(1,dif)(1))>DIMS(1)/2. ){
            tail(1,1) += DIMS(1)*sign(tail(1,dif));
          }
          if( abs(tail(2,dif)(1))>DIMS(2)/2. ){
            tail(2,1) += DIMS(2)*sign(tail(2,dif));
          }
          if( abs(tail(3,dif)(1))>DIMS(3)/2. ){
            tail(3,1) += DIMS(3)*sign(tail(3,dif));
          }
          plotl,tail(1,),tail(2,),color=reds(,ww(n)),width=1;
          plotl,tail(1,),tail(3,)-zoff,color=reds(,ww(n)),width=1;
        }
        plot, tail(1,0),tail(2,0),color="red",sym='\1';
        plot, tail(1,0),tail(3,0)-zoff,color="red",sym='\1';
        pause,1;
      }
    }
  }
  
}

/*
  
func klass_one(kits,cmap=){
  //kits(1,) == colID
  //kits(2,) == cRank
  //kits(3,) == up/dn

  //Looking for durable intercolumn motion
  thresh = 50;

  tmax = dimsof(kits)(0);
  mols = dimsof(kits)(-1);

  ii = 23;
  window,1,wait=1;
  fma;

  ok=-1;

  //shuffle columns
  colc = is_void(cmap)?indgen(30):cmap(:30);
  kitc = kits(1,,);

  extern SIEVE;
  SIEVE = CAT(10,,)*0;
  extern f_dura;
  extern t_dura;
  t_dura = indgen(200)(::4)-1;
  f_dura = t_dura*0;
  extern f_zdst;
  extern x_zdst;
  extern f_rdst;
  extern x_rdst;
  x_zdst=x_rdst=span(0,20,81);
  f_zdst=f_rdst=0*x_zdst;
  extern v_velor;
  extern f_velor;
  extern v_veloz;
  extern f_veloz;
  v_veloz = v_velor = span(0,2,101);  //ang/ps 
  f_veloz = f_velor = v_velor*0;
  extern x_jmp;
  extern f_jmp;
  x_jmp = [0,1,2];
  f_jmp = [0,0,0];
  
  for(cc=1; cc<=30; cc++){
    kitc(where(kits(1,,)==cc)) = colc(cc);
  }
  kits(1,,) = kitc;
  
  for(ii=1; ii<=600; ii++){
    inin = indgen(tmax)(where(kits(1,ii,)));
    sgkk = savgol(kits(1,ii,where(kits(1,ii,))),2*thresh+1,1);
    rsgkk = round(sgkk);
    foxx = kits(1,ii,inin);
        
    if(abs(ok)==1){
      Plot ,inin,kits(1,ii,inin),color="grayc";
      down = where( kits(3,ii,inin)<0);
      if(numberof(down)){
        plot,inin(down),kits(1,ii,inin)(down),color="black";
      }
      plotl,inin,sgkk,color=greenish(),width=4;
      plotl,inin,rsgkk,color="blue",width=2;
    }
    
    for(t=1; t<=numberof(inin)-thresh; t++){
      if(
         allof( sgkk(t+1:t+thresh) != sgkk(t) )*
         allof( foxx(max(1,t-thresh):t) == rsgkk(t) )
         ){

        mark = kits(2,ii,inin(t));
        dura = (where( foxx(t+1:)==rsgkk(t+1:) )(1) );
        home = int( foxx(t) );
        dest = int( foxx(t+dura) );
        if(home==dest) continue;

        if(abs(ok)==1){
          plot,inin(t),kits(1,ii,inin(t));
          plotl,inin([t,t+dura]),[home,dest],width=4,color="red";
        }

        
        write,format="\nMol %3d at %5d from col %2d --> %2d",ii,t,home,dest;

        //Jump Out
        coli = where( kits(1,,inin(t))==home );
        ncol = numberof(coli);
        topq =  kits(2,ii,inin(t))==ncol;
        botq = (kits(2,ii,inin(t))==1)*numberof(coli);
        wnor = where( kits(2,coli,inin(t))==!topq*(kits(2,ii,inin(t))+1)+topq)(1);
        wsou = where( kits(2,coli,inin(t))==!botq*(kits(2,ii,inin(t))-1)+botq)(1);
        nort = kits(3,coli,inin(t))(wnor);
        sout = kits(3,coli,inin(t))(wsou);
        cent = kits(3,ii,inin(t));

        
        SIEVE(ii,inin(t):inin(t+dura)) = 2;

        
        cgap = kits(2,coli,inin(t):inin(t+dura))([wnor,wsou],);
        cgap = abs(cgap(dif,));
        if(noneof(cgap==1)){
          write,format="  %s\n","Gavroche!";
          window,window()+1,wait=1;
        } else {
          tout = where(cgap==1)(1);
          if(abs(ok)==1) plot,inin(t:t+tout),kits(1,ii,inin(t:t+tout)),color="cyan";
          SIEVE(wnor,inin(t):inin(t+tout)) = 3;
          SIEVE(wsou,inin(t):inin(t+tout)) = 1;

          write,format="  ( %+2g : %+2g : %+2g )",sout,kits(3,ii,inin(t)),nort;
          f_jmp( int(.5*sum(abs([sout,cent,nort](dif))) )+1 )++;

          //Jump in
          coli = where( kits(1,,inin(t+dura))==dest );
          ncol = numberof(coli);
          topq =  kits(2,ii,inin(t+dura))==ncol;
          botq = (kits(2,ii,inin(t+dura))==1)*numberof(coli);
          wnor = where( kits(2,coli,inin(t+dura))==!topq*(kits(2,ii,inin(t+dura))+1)+topq)(1);
          wsou = where( kits(2,coli,inin(t+dura))==!botq*(kits(2,ii,inin(t+dura))-1)+botq)(1);
          nort = kits(3,coli,inin(t+dura))(wnor);
          sout = kits(3,coli,inin(t+dura))(wsou);
          cent = kits(3,ii,inin(t+dura)); 
          
          cgap = abs(kits(2,coli,inin(t:t+dura))([wnor,wsou],)(dif,));
          if(anyof(cgap==1)){
            t_in = where(cgap==1)(0);
            if(abs(ok)==1) plot,inin(t+t_in:t+dura),kits(1,ii,inin(t+t_in:t+dura)),color="magenta";
            SIEVE(wnor,inin(t+t_in):inin(t+dura)) = -3;
            SIEVE(wsou,inin(t+t_in):inin(t+dura)) = -1;
            write,format="  ~~  ( %+2g : %+2g : %+2g )",sout,kits(3,ii,inin(t+dura)),nort;
            f_jmp( int(.5*sum(abs([sout,cent,nort](dif))) )+1 )--;
          }
        }
        
        if(ok==2){
          xx=CAT(1,ii,inin(t):inin(t+dura)); 
          yy=CAT(2,ii,inin(t):inin(t+dura)); 
          zz=CAT(3,ii,inin(t):inin(t+dura)); 
          plotl,xx+yy*sqrt(3)/2,zz+yy*1/2.,color=colr;
        }

        fast_bin,dura,t_dura,f_dura;
        DIMS=DIM_ARR(,inin(t));
        dr = pbc( CAT(:3,ii,inin(t+dura))-CAT(:3,ii,inin(t)) );
        fast_bin,sqrt((dr(:2)^2)(sum)),x_rdst,f_rdst;
        fast_bin,dr(3),x_zdst,f_zdst;
        fast_bin,sqrt((dr(:2)^2)(sum))/dura,v_velor,f_velor;
        fast_bin,dr(3)/dura,v_veloz,f_veloz;
          
        flag=1;
        if(numberof(dura)==1) t+=dura;
                      
      }
    }
    if(flag){
      if(abs(ok)==1){
        plxtitle,"Molecule "+totxt(ii);
        pause,1000;
        if(ok>0) window,window()+1,wait=1;
        fma;
        write,"";
      }
    }
    flag = 0;
  }
  //winkill,window();
  write,"";
  lim = limits();
  plotl,[-10,10],[0,0],color="red",width=4;
  plotl,[0,0],[-10,10],color="red",width=4;
  plotl,[-10,10]*sqrt(3)/2,[-5,5],color="red",width=4;

  fma;
  plh,f_rdst,x_rdst;
  plh,f_zdst,x_zdst;
}

*/



func colog(t,tmax){
  leg = int(log(t)/log(10));
  t1 = 10^leg - 1;
  t2 = 10^(leg+1);

  /*
  step= 1500;
  leg = t/step;
  t1 = step*leg;
  t2 = step*(leg+1);
  */

  if(leg==0){
    //black to red
    return int([255.*(t-t1)/(t2-t1),0,0]);
  }else if(leg==1){
    //red to red/green
    return int([255.,255.*(t-t1)/(t2-t1),0]);
  }else if(leg==2){
    //red/green to green
    return int([255.*(t2-t)/(t2-t1),255.,0]);
  }else if(leg==3){
    //green to green/blue
    return int([0,255.,255.*(t-t1)/(t2-t1)]);
  }else if(leg==4){
    //green/blue to blue
    return int([0,255.*(t2-t)/(t2-t1),255.]);
  }else if(leg==5){
    //blue to blue/red
    return int([255.*(t-t1)/(t2-t1),0,255.]);
  }else if(leg==6){
    //blue/red to white
    return int([255.,255.*(t-t1)/(t2-t1),255.]);
  }
}

func colspin(t,tmax){
  return int([ 120 + 120*cos(2*pi*t/tmax), 120 + 120*sin(2*pi*t/tmax), 255.*t/tmax]);
}



func ffuunncc(fn,color=,nofma=){
  color = is_void(color)?"red":color;
  
  qq=openb(fn);
  restore,qq;

  window,8,wait=1;
  if(is_void(nofma)) fma;
  if(dimsof(prden)(0)<1000){
    plotl,dr_arr,prden(,0),color=color,width=2;
    write,format=" WARNING: %s is only %d long\n",fn,dimsof(prden)(0);
    xmax = dr_arr(where(prden(,0)>0)(0) );
  }
  else{
    plotl,dr_arr,prden(,1000),color=color,width=4;
    xmax = dr_arr(where(prden(,1000)>0)(0) );
  }
  plotl,dr_arr,prden(,100),color=color;
  plotl,dr_arr,prden(,10),color=color;
  plotl,dr_arr,prden(,1),color=color;
  plxtitle,"r_perp";
  plytitle,"2 pi G(r)";
  logxy,0,1;
  limits,0,xmax;

  window,9,wait=1;
  if(is_void(nofma)) fma;
  if(dimsof(pzden)(0)<1000){
    plotl,dz_arr,pzden(,0),color=color,width=2;
    xmax = dz_arr(where(pzden(,0)>0)(0) );
  }
  else{
    plotl,dz_arr,pzden(,1000),color=color,width=4;
    xmax = dz_arr(where(pzden(,1000)>0)(0) );
  }
  plotl,dz_arr,pzden(,100),color=color;
  plotl,dz_arr,pzden(,10),color=color;
  plotl,dz_arr,pzden(,1),color=color;
  plxtitle,"z_para";
  plytitle,"G(z)";
  logxy,0,1;
  limits,0,xmax;
  close,qq;
}




func gausstail(xx,yy){

  plotl,xx,yy,color="grayd",width=4;
  
  zero1 = (anyof(yy==0))?where(!yy)(1)-1:0;
  xx = xx(:zero1);
  yy = yy(:zero1);
  plotl,xx,yy,color="black",width=4;
  
  lmf = regress( log(yy), [xx^0.,xx^2.], sigy=1./xx );
  ao = exp(lmf(1));
  bb = -lmf(2);
  plotl,xx,ao*exp(-bb*xx^2),color="red";

  lmf = regress( log(yy), [xx^0.,xx^2.], sigy=xx );
  ao = exp(lmf(1));
  bb = -lmf(2);
  plotl,xx,ao*exp(-bb*xx^2),color="blue";
 

  //plotl,xx,yy - ao*exp(-bb*xx^2),color="blue";

  return sqrt(avg( (yy - ao*exp(-bb*xx^2))^2 ));

}

func plotl_color(xarr,yarr,c1,c2,width=){
  //plotl's xarr v yarr with color linearly shifting from c1 to c2
  dimx = dimsof(xarr);
  dimy = dimsof(yarr);
  numx = dimx(0);

  if(allof(dimx!=dimy)) return -999;
  /*
  if(allof(c1==c2)){
    if(dimx(1)==1) plotl,xarr,yarr,color=c1,width=width;
    else{
      for(s=1; s<=dimx(2); s++){
        plotl,xarr(s,),yarr(s,),color=c1,width=width;
      }
    }
  }*/
  carr = span(c1,c2,numx-1);
  dith = (random([2,dimx(2),3])*20)-10;
  for(l=numx-1; l>=1; l--){
    for(s=1; s<=dimx(2); s++){
      coll = carr(l,)+dith(s,);
      if(anyof(coll<0)) coll(where(coll<0)) = 0;
      if(anyof(coll>255)) coll(where(coll>255)) = 255;
      coll = int(coll);
      plotl,xarr(s,[l,l+1]),yarr(s,[l,l+1]),color=coll,width=width;
    }
  }
}

func dyn_suite(dirs,rat=,skip1=){
  rat = is_void(rat)?1:rat;
  rat = int(rat);
  if(skip1!=1){
    for(d=rat; d<=numberof(dirs); d++){
      write,format="\n%s\n",dirs(d);
      dir = dirs(d);
      
      qq=openb(dir+"cat.b");
      pp=createb(dir+"rotocat.b");
      
      tla = qq.CAT(4:6,,);
      dim_arr = qq.DIM_ARR;
      maxt = dimsof(qq.CAT)(0);
      eqt = int(qq.eqt);
      //eqt=1;
      
      tla = tla(,,eqt:);
      dim_arr = dim_arr(,eqt:);
      
      tla = unroll(tla);
      save,pp,tla;
      
      rotocat,tla;
      save,pp,f2avg,f2var,f4avg,f4var;
      close,pp;
    }
    rat = 1;
  }
  rat = abs(rat);
  for(d=rat; d<=numberof(dirs); d++){
    write,format="\n%s\n",dirs(d);
    dir = dirs(d);
    autorotVH3,dir;
  }
}

  
  
  
  
  
func unroll(tla){
  fhi = tla(3,,);
  nover = (abs(fhi(,dif))>pi)(,sum);
  mm = where(nover==max(nover))(1);
  nover = sum(nover);
  
  if(!nover) return 0;
  for(t=2; t<=dimsof(fhi)(0); t++){
    ww = abs(fhi(,t)-fhi(,t-1))>pi;
    mover = sum(ww);
    ww = where(ww);
    if(mover){
      fhi(ww,t:) -= 2*pi*sign( fhi(ww,t)-fhi(ww,t-1));
    }
    nover = sum( abs(fhi(,dif))>pi );
    write,format="\r%6d",nover;
    fma;
    plotl,,fhi(mm,dif),color="black";
    plotl,,fhi(mm,dif)(:t-1),color="red";
    if(!nover) break;
  }
  tla(3,,) = fhi;
  return tla;
}

      

func rotocat(tla,maxT=){
  //Direct calculation of

  Tmax = dimsof(tla)(0);
  maxT = is_void(maxT)?3200:maxT;
  maxT = min(maxT,Tmax-1);
  Nmol = MPC*NCOL;
 
  Nwrt = array(0,maxT);

  extern f2avg;  extern f2var;  
  extern f4avg;  extern f4var;  
  f2avg = f2var = array(0.,[2,Nmol,maxT]);
  f4avg = f4var = array(0.,[2,Nmol,maxT]);

  time = indgen(maxT);

  winkill,9;
  window,9,style="boxed.gs",wait=1;
  animate,1;
  
  for(t1=1; t1<Tmax; t1++){
    ti = 1./t1;
    write,format="\r%d/\%d ",t1,Tmax-1;
    dfi = tla(3,,t1) - tla(3,,t1+1:);
    if(dimsof(dfi)(0)>maxT) dfi = dfi(,:maxT);
    lent = dimsof(dfi)(0);
    
    Nwrt(:lent)++;

    df2 = (dfi^2);
    df4=df2^2;
    
    // Welford rolling avg and var
    welf = df2 - f2avg(,:lent);
    f2avg(,:lent) += welf * ti;
    welf*= df2 - f2avg(,:lent);
    f2var(,:lent) += welf;//Ns2

    welf = df4 - f4avg(,:lent);
    f4avg(,:lent) += welf * ti;
    welf*= df4 - f4avg(,:lent);
    f4var(,:lent) += welf;//Ns2

    plotl,time,sqrt(f2avg(avg,)),color="grayc";
    plotl,time(lent:),sqrt(f2avg(avg,lent:)),color="blue",width=4;
    plxtitle,"Time displacement";
    plytitle,"rmsd";
    fma;   
    
  }
  
  animate,0;

  Nwrt = (1./Nwrt)(-:1:Nmol,);
  f2var *=  Nwrt;
  f4var *=  Nwrt;
}

func autorotVH(dir){
  qq=updateb(dir+"rotocat.b");
  tla = qq.tla;
  if( anyof(abs(tla(3,,dif)>pi)) ){
    tla = unroll(tla);
    save,qq,tla;
  }
  Tmax = 1 + int(sqrt(10)*1000);
  
  Tmax = min(Tmax,dimsof(tla)(0));
  Tmax = int(Tmax);
  Nbin = int(1+(600*Tmax)^(1./3));
  Gmax = pi*4;
  dg_arr = (span(0,sqrt(Gmax),Nbin+1)^2)(2:); //NON CONSTANT DX
  fg_arr = array(0,[2,Nbin,Tmax]);

  gam = tla(3,,);
  
  for(dt=1; dt<=Tmax; dt++){
    write,format="\r %5d / %5d",dt,Tmax;
    for(t=1; t<=dt; t++){
      ftmp = array(0,Nbin);
      dgam = gam(,t::dt);
      if(dimsof(dgam)(0)==1){
        fast_bin,abs(dgam),dg_arr,ftmp;
      } else{
        fast_bin,abs(dgam(,dif)),dg_arr,ftmp;
      }
      fg_arr(,dt) += ftmp;
    }
    if(random()>0.9){
      plotl,dg_arr(zcen)*180/pi,dt+fg_arr(:-1,dt)*dg_arr(0)/dg_arr(dif)/fg_arr(sum,dt),color=reddish();
      pause,0;
    }
  }
  save,qq,fg_arr,dg_arr;
  close,qq;
}

func autorotVH2(dir){
  qq=updateb(dir+"rotocat.b");
  tla = qq.tla;
  if( anyof(abs(tla(3,,dif)>pi)) ){
    tla = unroll(tla);
    save,qq,tla;
  }
  gam = tla(3,,);
  Tmax = 1 + int(sqrt(10)*1000);
  Tmax = min(Tmax,dimsof(tla)(0));
  Tmax = int(Tmax);
  Nbin = int(1+(600*Tmax)^(1./3));

  Gmax = max(gam(,max)-gam(,min));
  gmax = min(Gmax,2*pi);
  dg_arr = (span(0,gmax^(.8),Nbin+1)^(1.2))(2:); //NON CONSTANT DX
  //dg_arr = span(0,gmax,Nbin+1)(2:); //CONSTANT DX
  grow,dg_arr,Gmax+0.01;
  Nbin+=1;
  fg_arr = array(0,[2,Nbin,Tmax]);

  shuf = Tmax-1;
  grow,shuf,sort(random(Tmax-2));

  ndg = numberof(dg_arr);
  for(t1=1; t1<Tmax; t1++){
    t0 = shuf(t1);
    write,format="\r%d/\%d ",t1,Tmax;
    dg = abs( gam(,t0) - gam(,t0+1:) );
    digg = digitize(dg,dg_arr);

    for(i=1; i<=600; i++){
      for(dt=1; dt<=Tmax-t0; dt++){
        if(digg(i,dt)<=ndg) fg_arr(digg(i,dt),dt)++;
      }
    }
    if(random()>0.9){
      fma;
      limits,1,dg_arr(-1)*180/pi,1,max(shuf(:t1));
      plf,log0(Gmax*fg_arr(:-1,)/dg_arr(dif)(,-:1:Tmax)/(fg_arr(:-1,)(sum,)(-:1:Nbin-1,) +1) ),indgen(Tmax)(-:1:Nbin-1,),dg_arr(zcen)(,-:1:Tmax)*180/pi;
      //plf,log0(Gmax*fg_arr(:-1,)/dg_arr(dif,-:1:Tmax)),indgen(Tmax)(-:1:Nbin-1,),dg_arr(zcen)(,-:1:Tmax);
      pause,2;
    }
  }

  dg = dg_arr(1);
  grow,dg,dg_arr(dif);
  dg = dg(:-1);
  pg_arr = Gmax*fg_arr(:-1,)/dg(,-:1:Tmax);
  den = fg_arr(sum,)(-:1:Nbin-1,);
  pg_arr(where(den))/= den(where(den));
  save,qq,pg_arr,dg_arr,fg_arr;
  close,qq;
}



func autorotVH3(dir){
  qq=updateb(dir+"rotocat.b");
  tla = qq.tla;
  if( anyof(abs(tla(3,,dif)>pi)) ){
    close,qq;
    qq=createb(dir+"rotocat.b");
    tla = unroll(tla);
    save,qq,tla;
  }
  gam = tla(3,,);
  dTmax = 3201;
  dTmax = min(dTmax,dimsof(tla)(0));
  Tmax = dimsof(tla)(0);
  Tmax = int(Tmax);
  Nbin = int(1+(600*Tmax)^(1./3));
  Nbin = min(300*8.5/2.54,Nbin);
  Nbin = int(Nbin);

  Gmax = max(abs(gam(,max)-gam(,min)));
  gmax = min(Gmax,4*pi);
  dg_arr = (span(0,gmax^.8,Nbin+1)^(1.25)); //NON CONSTANT DX
  grow,dg_arr,Gmax+0.01;
  Nbin+=1;
  fg_arr = array(0,[2,Nbin+1,dTmax]);

  shuf = indgen(dTmax);

  ndg = numberof(dg_arr);
  animate,1;
  for(dt0=1; dt0<=dTmax; dt0++){
    dt = shuf(dt0);
    write,format="\r%d/\%d ",dt0,dTmax;
    
    fgt = fg_arr(,dt);
    cap = min(dt,Tmax-dt);
    for(t0=1; t0<=cap; t0++){
      dg = abs( gam(,t0::dt)(,dif) );
      fast_bin,dg,dg_arr,fgt;
    }
    fg_arr(,dt)+= fgt;
    fgt*= 0;

    Plotl,dg_arr(zcen)*180/pi,fg_arr(:-1,dt)/dg_arr(dif),color="black";
    pltitle,dir+" ROTO";

  }
  animate,0;
  dg = dg_arr(dif);
  dg = dg(,-:1:dTmax);
  pg_arr = fg_arr(:-1,)/dg;
  den = fg_arr(sum,)(-:1:Nbin-1,);
  pg_arr(where(den))/= den(where(den));
  f2avg = qq.f2avg;
  f4avg = qq.f4avg;
  f2var = qq.f2var;
  f4var = qq.f4var;
  close,qq;
  qq = createb(dir+"rotocat.b");
  save,qq,dg_arr,fg_arr,pg_arr;
  save,qq,tla;
  save,qq,f2avg,f2var;
  save,qq,f4avg,f4var;
  close,qq;
}


func plot_alot_colors(xxx,yyy,col1=,col2=){
  tmax = dimsof(xxx)(0);
  num = dimsof(xxx)(-1);
  col1 = is_void(col1)?[0,57,118]:col1;
  col2 = is_void(col2)?[239,171,0]:col2;


  colt = int(span(col1,col2,tmax) );

  ot = int( 10^(indgen( 0:int( log(4000)/log(10)  *2) )*.5) );
  for(n=2; n<=numberof(ot); n++){
    col = colt(ot(n-1),);
    //colt(ot(n-1):ot(n)-1, ) =  col(-:1:ot(n)-ot(n-1),);
  }
  //colt(ot(0):,) = colt(ot(0) ,)(-:1:tmax-ot(n-1),);
      

  for(t=tmax-1; t>=1; t--){
    for(i=1; i<=num; i++){
      colr = colt(t+1,);
      //colr = colt(t+1,avg)(-:1:3);
      plotl,xxx(i,[t,t+1]),yyy(i,[t,t+1]),color=colr,width=4;
    }
  }

}
                               


func otrarotVH2(dir){
  qq=openb(dir+"rotocat.b");
  tla = qq.tla;
  if( anyof(abs(tla(3,,dif)>pi)) ){
    tla = unroll(tla);
    save,qq,tla;
  }
  if(is_void(qq.dg_arr)) return -999;
  close,qq;

  return "CURRENTLY UNDER CONSTRUCTION";
  
  gam = tla(3,,);
  //dg1dg2dt = span(-pi,pi,121)(,-:1:122,-:1:nT+1);
  dg1dg2dt = array(0,[3,121,122,nT+1]);
  dg1dg2dtx = array(0,[3,121,122,nT+1]);

  animate,1;
  ravg = 0.;
  navg = 0;
  for(t=1; t<=nT; t++){
    if(navg) write,format="\r %g",ravg/navg;
    for(i=1; i<=600; i++){

      for(k=1; k<=2; k++){

        j = nnz(k,i,t);
        if(j==0) continue;
        tchk = where( (nnz(,i,t::5)==j)(sum,));
        if(!numberof(tchk)) continue;
        treal = (tchk-1)*5+1;

        dgi = gam(i,t) - gam(i,treal);
        dgj = gam(j,t) - gam(j,treal);

        ravg+= dgi(sum);
        navg+= numberof(dgi);
        
        dgi = int( (dgi+pi)*60/pi);
        dgi = min(dgi,122);
        if(anyof(dgi<0) ) dgi(where(dgi<=0)) = 122;
        dgj = int( (dgj+pi)*60/pi);
        dgj = min(dgj,122);
        if(anyof(dgj<0) ) dgj(where(dgj<=0)) = 122;

        if(allof(dgi==122)) continue;

        for(tn=1; tn<=numberof(tchk); tn++){
          if(dgi(tn)!=122){
            dg1dg2dt( dgi(tn), dgj(tn), tchk(tn) )++;
          }
        }
        

      }

    }

    if(random()>.95){
      fma;
      pli,dg1dg2dt(,:121,200),-180,-180,180,180;
    }
    
  }

    ravg = 0.;
  navg = 0;
  for(t=1; t<=nT; t++){
    if(navg) write,format="\r %g",ravg/navg;
    for(i=1; i<=600; i++){

      for(k=1; k<=18; k++){

        j = nnr(k,i,t);
        if(j==0) continue;
        tchk = where( (nnr(,i,t::5)==j)(sum,));
        if(!numberof(tchk)) continue;
        treal = (tchk-1)*5+1;

        dgi = gam(i,t) - gam(i,treal);
        dgj = gam(j,t) - gam(j,treal);

        ravg+= dgi(sum);
        navg+= numberof(dgi);
        
        dgi = int( (dgi+pi)*60/pi);
        dgi = min(dgi,122);
        if(anyof(dgi<0) ) dgi(where(dgi<=0)) = 122;
        dgj = int( (dgj+pi)*60/pi);
        dgj = min(dgj,122);
        if(anyof(dgj<0) ) dgj(where(dgj<=0)) = 122;

        if(allof(dgi==122)) continue;

        for(tn=1; tn<=numberof(tchk); tn++){
          if(dgi(tn)!=122){
            dg1dg2dtx( dgi(tn), dgj(tn), tchk(tn) )++;
          }
        }
        

      }

    }

    if(random()>.95){
      fma;
      pli,log0(dg1dg2dtx(,:121,200)),-180,-180,180,180;
    }
    
  }


  

  animate,0;

  qq=updateb(dir+"rotocat.b");
  save,qq,dg1dg2dt,nnz,dg1dg2dtx,nnr;
  close,qq;
  
}



func meet_neighbors(dir){
  
  qq = openb(dir+"cat.b");
  CAT = qq.CAT;
  ti = int(qq.eqt);
  tf = dimsof(CAT)(0);
  comt = CAT(1:3,,);
  cent = CAT(7:9,,);
  dimt = qq.DIM_ARR;
  dirt = ni_t(cat=CAT);
  close,qq;
  
  nnz = array(0,[3,2,600,tf-ti+1]);
  nnr = array(0,[3,18,600,tf-ti+1]);
  
  for(t=ti; t<=tf; t++){
    write,format="\r%d / %d",t,tf;
    coms = comt(,,t);
    cens = cent(,,t);
    voro = voroneigh(coms,dimt(,t));
    dirs = dirt(,,t);
    DIMS = dimt(,t);
    
    for(i=1; i<=MOLS; i++){//assign mols to cols
      neigh = voro(,i,1);
      neigh = neigh(where(neigh));
      neigh = neigh(where(neigh!=i));
      neigh = int(neigh);  //voro neighbor absolute indices

      dndot = abs(dirs(+,i)*dirs(+,neigh));
      neigh = neigh(where(dndot>.25));
      if(!numberof(neigh)) continue;
      //cofacial within 75degrees
      
      dr = pbc( cens(,neigh) - cens(,i) );
      dpara = para2(dr,dirs(,i));
      dperp = perp2(dr,dirs(,i));
      
      sgn = sign( dirs(+,i)*dpara(+,) );
      //above or below REL TO director
      
      dpara = norm(dpara);
      dperp = norm(dperp);
      
      incol = ( (dpara<5.)*(dperp<5.25) );
      if(anyof(incol)){
        ineigh = neigh(where(incol));
        inpara = (dpara*sgn)(where(incol));
        if(anyof(inpara>0)){
          lo = min(inpara(where(inpara>0)) );
          top = ineigh(where(inpara==lo)(1));
          nnz(1,i,t-ti+1) = top;
        }
        if(anyof(inpara<0)){
          hi = max(inpara(where(inpara<0)) );
          bot = ineigh(where(inpara==hi)(1));
          nnz(2,i,t-ti+1) = bot;
        }
      }

      outneigh = ( (neigh!=top)*(neigh!=bot) );
      outpara = dpara(where(outneigh));
      outperp = dperp(where(outneigh));
      outneigh = neigh(where(outneigh)); //absolute indices

      numneigh = numberof(outneigh);
      if(numberof(numneigh<=18)){
        nnr(:numneigh,i,t-ti+1) = numneigh;
      }else{
        aa = outpara <= outpara((sort(outpara)(:18))(0));
        bb = outperp <= outperp((sort(outperp)(:18))(0));
        outneigh = outneigh(where(aa*bb));
        numneigh = numberof(outneigh);
        nnr(:numneigh,i,t-ti+1) = numneigh;
      }

    }
  }

  qq = createb(dir+"pound.b");
  save,qq,nnz,nnr;
  close,qq;

  write,format="Average in neighbors: %g\n",avg((nnz!=0)(sum,,));
  write,format="Average ex neighbors: %g\n",avg((nnr!=0)(sum,,));
}




func orphan_mobility(dir){

  qq = openb(dir+"cat.b");

  tilt = qq.CAT(4,,qq.eqt:);
  dogz = qq.DOG(3,,qq.eqt:);

  sideQ = abs(tilt - pi/2)<pi/5;
  sideQ+= abs(tilt - pi/2)<pi/4;

  Tlen = dimsof(dogz)(0);

  trapt = array(0,Tlen);
  orpht = array(0,600);

  //orphan data
  dzzdt = array(0.,Tlen);
  numdt = array(0,Tlen);
  //reg data
  dzzdtr = array(0.,Tlen);
  numdtr = array(0,Tlen);

  animate,1;
  //Orphans
  for(m=1; m<=600; m++){
    write,format="\r%d",m;
    if(noneof(sideQ(m,))) continue;
    tilted = min(1,2*ceil(savgol(sideQ(m,),51,6)));
    if(noneof(tilted)) continue;

    startat = endedat = [];
    if(anyof(tilted(dif)>0)) startat = where(tilted(dif)>0);
    if(anyof(tilted(dif)<0)) endedat = where(tilted(dif)<0)+1;
    if(tilted(1)==1){
      tmp = 1;
      grow,tmp,startat;
      startat = tmp;
      tmp = [];
    }
    if(tilted(0)==1){
      grow,endedat,numberof(tilted);
    }

   
    write,format="\r %d ! (%d)",m,numberof(startat);
    for(n=1; n<=numberof(startat); n++){
      twindow = [startat(n),endedat(n)];
      if(twindow(0)==twindow(1)) continue;
      if(twindow(0)!=numberof(tilted)){
        duration = int(twindow(dif)(1));
        trapt(int(duration))++;
        if(duration>2) orpht(m)++;
      }
      deltaz = dogz(m,twindow(1):twindow(0));
      for(t=1; t<numberof(deltaz); t++){
        dzz = (deltaz(t) - deltaz(t+1:))^2;
        dzzdt(:numberof(dzz)) += dzz;
        numdt(:numberof(dzz)) += 1;
      }
    }

    Plotl,indgen(Tlen)(where(numdt)),dzzdt(where(numdt))/numdt(where(numdt)),color="red";
   
  }

  write,format="\n %d Orphans",sum(!(!orpht));
  if(anyof(numdt)) dzzdt(where(numdt))/=numdt(where(numdt));

  animate,0;
  return [dzzdt,numdt,dzzdtr,numdtr,trapt*1.];

    //Non-orphans
  for(m=1000; m<=600; m++){
    write,format="\r%d",m;
    if(allof(sideQ(m,))) continue;
    tilted = min(1,2*ceil(savgol(sideQ(m,),51,6)));
    if(allof(tilted)) continue;

    startat = endedat = [];
    if(anyof(tilted(dif)<0)) startat = where(tilted(dif)<0);
    if(anyof(tilted(dif)>0)) endedat = where(tilted(dif)>0)+1;
    if(tilted(1)==0){
      tmp = 1;
      grow,tmp,startat;
      startat = tmp;
      tmp = [];
    }
    if(tilted(0)==0){
      grow,endedat,numberof(tilted);
    }
    
    //write,format="\r %d !\n",m;
    for(n=1; n<=numberof(startat); n++){
      twindow = [startat(n),endedat(n)];
      deltaz = dogz(m,twindow(1):twindow(0));
      for(t=1; t<numberof(deltaz); t++){
        dzz = (deltaz(t) - deltaz(t+1:))^2;
        dzzdtr(:numberof(dzz)) += dzz;
        numdtr(:numberof(dzz)) += 1;
      }
    }
    

    
    //Plotl,indgen(Tlen)(where(numdtr)),dzzdtr(where(numdtr))/numdtr(where(numdtr)),color="blue";
    plotl,indgen(Tlen)(where(numdt)),dzzdt(where(numdt))/numdt(where(numdt)),color="red";
   

  }

  dzzdtr(where(numdtr))/=numdtr(where(numdtr));
  animate,0;
  return [dzzdt,numdt,dzzdtr,numdtr,trapt*1.];
}
        
      
    
