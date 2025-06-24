#include "~/TinkerExp/sn3.i"
#include "~/TinkerExp/whisker.i"
/*

  Le Petite Donéeoduc

  built to interact with the output of the sn3 function parse_traj

  Goals

  Parse trajectories into CATs
  Estimate Equillibrium Time

  Principle Funcs
  avant_leduc,prefix,dirs=,outdirs=
  ->  prepares CAT etc
  merge_cat,dirs,outdir
  ->  merges CATS etc from dirs to outdir
  ->  also parses outs and runs energy width
  ->  can be used with just one cat
  add_eqs,dirs,eqs
  ->  adds eqt to cat.b
  turbine,dirs
  ->  build msd in both lab and mol frame
  
*/

 write,format="\n %s","\bRun something like:","\n avant_leduc,prefix,dirs=trajdirs,outdirs=destdirs,temps=temperatures","\n merge_cats,dirs,outdir","\nThen get eqts and","\n add_eqs,dirs,eqts","\n turbine,dirs";

func avant_leduc(prefix,dirs=,outdirs=,temps=){
  /* DOCUMENT avant_leduc(dirs=,outdirs=)
     prefixes is the trajectory root
     dirs is directory with trajectory
     outdirs is where to put the cat
   */
  dirs = is_void(dirs)?"./":dirs;
  if(allof(strpart(dirs,0:)!="/")) dirs+="/";
  outdirs = is_void(outdirs)?dirs:outdirs;
  if(allof(strpart(outdirs,0:)!="/")) outdirs+="/";
  for(n=1; n<=numberof(dirs); n++){
    faire_lechat(prefix,dirs(n),outdir=outdirs(n));
    outs = ls(dirs(n)+prefix+".[0-9].out");
    qq = createb(outdirs(n)+"out_list.b");
    save,qq,outs;
    close,qq;
    if(!is_void(temps)){
      window,8,wait=1;
      alph = int(min(1000,dimsof(DIM_ARR)(0)*0.5))
      plot,temps(n),exp(log(DIM_ARR(,:-alph))(sum)),color=["red","blue"]( int(strchar(prefix)(sum)%2));
      window,0,wait=1;
    }
  }
  write,"Finished making cats.";
  
    
}
func faire_lechat(prefix,dir,outdir=){
  //some numbers, in case they haven't been said
  extern MPC;
  MPC = 20;
  extern NCOL;
  NCOL = 30;
  extern MOLS;
  MOLS = MPC*NCOL;

  outdir = is_void(outdir)?dir:outdir;

  NCF = count_files(prefix,dir=dir);
  //set up arrays for order parameters
  extern oS_arr; //theta order param xarr                    //coded
  oS_arr = array(0.,NCF);
  extern cS_arr;
  cS_arr = array(0.,[2,NCOL,NCF]);

  write,format="Parsing %s in %s\n",prefix,dir;
  oS_arr = parse_traj(prefix,dir=dir,funct=oriS);
  write,format="\rParsed and saved to %s\n",outdir;
  
  winkill,0;
  window,0,style="work.gs",wait=1;
  fma;
  write,format="%s\n","   Avg OP's:";
  plxtitle,"OP at various times (500ps each)";
  for(n=1; n<=numberof(oS_arr); n+=500){
    alph = n; omeg = min(n+500,numberof(oS_arr));
    tmp_min = max(0., min(oS_arr(alph:omeg) * 0.9) );
    tmp_max = min(1., max(oS_arr(alph:omeg) * 1.1) );
    tmp_clr = [int(240*n/(numberof(oS_arr)*1.)),0,0];
    bb = fast_bin(oS_arr(alph:omeg),span(tmp_min,tmp_max,100) );
    plh,bb(,2)+(n/50.),bb(,1),color=tmp_clr;
    limits,0,1,e,e;
    Savg = oS_arr(alph:omeg)(avg);
    Slong= oS_arr(alph:)(avg);
    Sstd = sqrt(avg((oS_arr(alph:omeg) - Slong)^2));
    write,format="    %d : %d | %.3g  ( %.3g )\n",alph,omeg,Savg,Sstd;
  }

  write,"~~~Saving CAT/DOG/DIMS~~~";
  qq=createb(outdir+"cat.b");
  save,qq,CAT,DIM_ARR;
  DOG = CAT;//unwrap_traj(CAT,DIM_ARR);
  save,qq,DOG;
  close,qq;
  write,format="\rSaved CAT to %s\n",outdir+"cat.b";


  write,"~~~Saving Arrays~~~";
  qq=createb(outdir+"cat_data.b");
  save,qq,oS_arr,cS_arr;
  write,format="\rSaved oS_arr array to %s\n",outdir+"cat_data.b";

  nit = ni_t();

  CAT(10,,) = column_tracker(nit=nit);
  qq = updateb(outdir+"cat.b");
  save,qq,CAT; close,qq;
  write,format="\n%s","Tracked columns.";
}



func merge_cat(dirs,outdir){
  write,format="\n%s\n\n","DIRS MUST BE IN ORDER";
  if(anyof(strpart(dirs,0:)!="/")) dirs(where(strpart(dirs,0:)!="/"))+="/";
  if(strpart(outdir,0:)!="/") outdir+="/";
  /*while(catch(-1)){
    write,format="ERR: No cat in %s. Moving to next\n",dirs(1);
    dirs = drop(dirs,1);
    }*/
  qq=openb(dirs(1)+"cat.b");
  cat = qq.CAT;
  dim_arr = qq.DIM_ARR;
  close,qq;
  pp = openb(dirs(1)+"out_list.b");
  outs = pp.outs;
  close,pp;
  write,format="\r CAT is %d long.",dimsof(cat)(0);
  for(d=2; d<=numberof(dirs); d++){
    /*while(catch(-1)){
      write,format="ERR: No cat in %s. Moving to next\n",dirs(d);
      if(d==numberof(dirs)) exit;
      dirs = drop(dirs,d);
    }*/
    qq = openb(dirs(d)+"cat.b");
    grow,cat,qq.CAT;
    grow,dim_arr,qq.DIM_ARR;
    close,qq;
    pp = openb(dirs(d)+"out_list.b");
    grow,outs,pp.outs;
    close,pp;
    write,format="\r CAT is %d long.",dimsof(cat)(0);
  }
  CAT=cat; cat=[];
  DIM_ARR=dim_arr; dim_arr=[];
  OUTS=outs; outs=[];
  qq = createb(outdir+"cat.b");
  save,qq,CAT,DIM_ARR,OUTS;
  if( catch(-1)){
    write,"Dog Problems";
    return 0;
  }
  DOG = unwrap_traj(CAT,DIM_ARR);
  save,qq,DOG;
  close,qq;
  if( catch(-1)){
    write,"Out Problems";
    return 0;
  }
  winkill,0;
  window,0,wait=1,style="boxed.gs";
  parse_out,OUTS,concat=1;
  energy_width;
  freeat = where(!GEOE)(1)*0.5;
  plotl,freeat*[1,1],limits()(3:4),color="magenta";
  png,outdir+"ewidth";
  fma;
  plotl,,DIM_ARR(1,),color="red";
  plotl,,DIM_ARR(2,),color="blue";
  plotl,,DIM_ARR(3,),color="green";
  plotl,,(DIM_ARR(1,)*DIM_ARR(2,)*DIM_ARR(3,))^(1./3),color="black",width=4;
  plotl,freeat*[1,1],limits()(3:4),color="magenta";
  plxtitle,"Time (ps)";
  plytitle,"Dims (A)";
  png,outdir+"dimst";
}

func add_eqs(dirs,eqts){
  for(d=1; d<=numberof(dirs); d++){
    qq = updateb(dirs(d)+"cat.b");
    len = dimsof(qq.CAT)(0);
    eqt = len;
    if(len>eqts(d)) eqt = eqts(d);
    else write,format="ERR: %s has short CAT!!\n",dirs(d);
    save,qq,eqt;
    close,qq;
  }
}

func turbine(dirs){
  write,"Dynamics in lab frame";
  dynarun,dirs,lab=1;
  write,"\rDynamics in mol frame";
  dynarun,dirs;
}

























func le_duc(prefix,dir,eqtime=){

  extern CAT;
  extern DIM_ARR;
  CAT=DIM_ARR=DOG=[]; //kill any current CAT
  write,format="Opening %s for CAT\n",dir+"cat.b";
  dat = openb(dir+"cat.b");
  if(dat.CAT==[]){
    write,"Did not find CAT.";
    return -999;
  }
  if(is_void(eqtime)){
    if(dat.eqt==[]){
      write,"ERR: Ain't got no time for this.";
      return -998;
    }
    eqtime = dat.eqt;
  }
  if(eqtime<0) eqtime += dimsof(dat.CAT)(0);
  if(abs(eqtime) >= dimsof(dat.CAT)(0)){
    write,"ERR: Bad Time.";
    return -997;
  }
  
  
  //trim relevant arrs
  restore,dat;
  CAT = CAT(,,eqtime:);
  DIM_ARR = DIM_ARR(,eqtime:);
  close,dat;
  write,format="This CAT is %d ps long.\n",dimsof(CAT)(0);

  qq = openb(dir+"cat_data.b");
  restore,qq;
  close,qq;

  extern MOLS;
  MOLS = dimsof(CAT)(3);
  times = dimsof(DIM_ARR)(0);


  oS_arr = oriS_cat(cat=CAT);
  cS_arr = array(0.,[2,NCOL,times]);
  oH_arr = oS_arr;
  cH_arr = cS_arr;
  HAT = CAT(10,,)*0.0;
 

  //Estimate number of data points
  minD = int(min(DIM_ARR)*5.)/10.;
  vfrac= 1.333*pi*minD^3.;
  vfrac/= DIM_ARR(1,)*DIM_ARR(2,)*DIM_ARR(3,);
  vfrac = vfrac(avg);
  //pairN is number of bins to use for pairwise interactions
  //based on the Rice Rule
  //alt mins to approximately 0.025A bins or 2^10 bins 
  pairN = int( 2. * (vfrac*MOLS*(MOLS-1)*times)^(1/3.) )+1;
  pairN = min(pairN , int(minD*40) , 1024);
  
  //cageN is number of bins to use for intra cage interactions
  //based on the Rice Rule
  //alt mins to 2^10 bins
  cageN = int( 2 * (6*MOLS*times)^(1/3.) ) + 1;
  cageN = min(1024 , cageN);

  //bulkN is number of bins for per frame quantities
  //strictly for 1D distros
  //based on sqrt(N)
  //alt mins to 2^10 bins
  bulkN = int( sqrt( MOLS*times)) + 1;
  bulkN = min(1024 , bulkN);

  //Set up arrays for 2D Pair Correlation Function
  extern gx_arr; //perp array for g pcf                      //coded
  extern gy_arr; //para array for g pcf                      //coded
  extern gf_arr; //freq array for g pcf                      //coded
  gy_arr = span(0.,minD,pairN+1)(2:);  
  gx_arr = gy_arr(,-:1:pairN);
  gy_arr = gy_arr(-:1:pairN,);
  gf_arr = gx_arr*0;

  //set up arrays for 2D orrientation correlations
  extern lx_arr; //ith molecule lambda array                 //coded
  extern ly_arr; //jth molecule lambda array                 //coded
  extern lf_arr; //lambda corr heatmap freqs                 //coded
  extern lr_arr; //lat distance w.r.t lambda
  extern lf2_arr;//lat distance freqcy array 
  lx_arr = ly_arr = span(-pi,pi,cageN);
  lx_arr = lx_arr(,-:1:cageN);
  ly_arr = ly_arr(-:1:cageN,);
  lf_arr = lx_arr*0;
  lr_arr = span(5.,15.,cageN)(-:1:cageN,);
  lf2_arr= lf_arr*0;

  //set up arrays for "cage" relative orientations
  extern na_arr; //neighbor angles theta jik                 //coded
  extern nf_arr; //'cause hexatic op is hard                 //coded
  na_arr = span(0.,pi,cageN);
  nf_arr = na_arr*0.;
  extern ma_arr; //intracol nn dot products                  //coded
  extern mf_arr; //named out of convenience                  //coded
  ma_arr = span(-1.,1.,cageN);
  mf_arr = ma_arr*0;
  extern zx_arr; //intercol lab z offsets                    //coded
  extern zf_arr; //the x suggests crosses                    //coded
  zx_arr = span(0.,6.,cageN);
  zf_arr = zx_arr*0;
  extern xd_arr; //intercol dot product                      //coded
  extern xf_arr; //do names even matter                      //coded
  xd_arr = span(-1.,1.,cageN);
  xf_arr = xd_arr*0.;

  //set up arrays for col vibes and tilt
  extern cr_arr; //bins for distance from colaxis            //coded
  extern cf_arr; //freqs to go with the bin array            //coded
  cr_arr = span(0.,4.0,bulkN)(2:);
  cf_arr = cr_arr*0.;
  extern th_arr; //bins for angle between mols and cols      //coded
  extern tf_arr; //frequency bins for intracolumn tilts      //coded
  th_arr = span(0.,1*pi/3.,bulkN);
  tf_arr = th_arr*0;
  extern ph_arr; //bins for intracolumn delta phi angle
  extern pf_arr; //frequency distributions for the same
  extern pf2_arr; //freq dist for signed, in column phi
  ph_arr = span(-pi,pi,bulkN);                               //coded
  pf_arr = ph_arr*0;                                         //coded
  pf2_arr= ph_arr*0;
  extern dz_arr; //bins for intracol nn sep para to col 
  extern df_arr; //frequency distributions for the same
  dz_arr = span(2.,6.,bulkN);                                //coded
  df_arr = dz_arr*0;                                         //coded

  write,format="%s\n","Numerous externals have been setup.\n   (Voir apres le duc si tu-plait.)";

  write,"ni_ting..."
  dirt = ni_t(cat=CAT);
  
  write,"Allons a le loop";
  for(t=1; t<=times; t++){
    write,format="\r%d/%d",t,times;
    DIMS = DIM_ARR(,t);
    coms = CAT(1:3,,t);
    angs = CAT(4:6,,t);
    dirs = dirt(,,t);

    for(i=1; i<=MOLS; i++){
      dr = pbc( coms - coms(,i) );

      paradr = para2(dr,dirs(,i));
      perpdr = dr - paradr;
      
      paradr = abs(paradr(1,),paradr(2,),paradr(3,));
      perpdr = abs(perpdr(1,),perpdr(2,),perpdr(3,));
      
      separr = [perpdr,paradr];
      fast_bin2,separr,gy_arr(1,),gx_arr(,1),gf_arr;
      separr = [];

      //build neighbor list and do intracol calcs
      //neighbors are from cylinder: 2.0 < |z| < 6.0
      //                                    r <= 6.0
      uni = where((perpdr<=6.0)*(paradr>2.0)*(paradr<=6.0));
      if(numberof(uni)>=1){
        indot = dirs(+,uni)*dirs(+,i);
        fast_bin,indot,ma_arr,mf_arr;
        
        nflip = sign(indot);
        d_phi = CAT(6,uni,t)-CAT(6,i,t); //fhj - fhi
        ww = where(abs(d_phi)>pi);   //positive for ccw rots of mol
        if(numberof(ww)){            //when dir is +z in lab frame.
          d_phi(ww) -= 2*pi*sign(d_phi(ww));
        }
        d_phi *= nflip; //flip when orientations are reversed
        
        //POSITIVE D_PHI IS FROM 4F TOWARD 1F
        fast_bin,d_phi,ph_arr,pf_arr;
        
      }
    
      //build neighbor list and do intercol calcs
      //neighbors are from an annulus: 7.0 <  r <= 15.0
      //                                     |z|<= 2.5
      //(zcut is used for the r ID because it is arbitrarily
      //smaller than that dimension of the zone of exclusion
      uni = where((perpdr<=15.0)*(perpdr>7.0)*(paradr<=2.5));
      hex_run = 0.;
      if(numberof(uni)>=2){
        gal = uni( where( abs(dirs(+,i)*dirs(+,uni))>=0.5 ));
        if(numberof(gal)>=1){
          lamij = pair_or(i,gal,t);
          fast_bin2,lamij(,[1,2]),ly_arr(1,),lx_arr(,1),lf_arr;
          fast_bin2,lamij(,[1,3]),lr_arr(1,),lx_arr(,1),lf2_arr;
          //intercol lab z offest
          dzx = abs(dr(3,gal));
          fast_bin,dzx,zx_arr,zf_arr;
          //intercol dot products
          xdt = dirs(+,gal)*dirs(+,i);
          fast_bin,xdt,xd_arr,xf_arr;
        }
        //use neighbor list to build hexatic OP
        Rij = perp2( dr(,uni) , dirs(,i) ); //in ith molec plane!
        Rij/= norm(Rij)(-:1:3,); //unit vecs rij
        for(j=1; j<numberof(uni); j++){
          d_jk = Rij(+,j)*Rij(+,j+1:);
          d_jk = min(d_jk,1.0000);
          d_jk = max(d_jk,-1.0000);
          thjk = acos(d_jk);
          fast_bin,thjk,na_arr,nf_arr;
          hex_run += sum( cos(6.*thjk) );
        }
        hex_run /= (numberof(uni)^2 - numberof(uni) )/2;
        
      }
      HAT(i,t) = hex_run;
    }
    oH_arr(t) = HAT(avg,t);

    for(ci=1; ci<=NCOL; ci++){
      coli= where(CAT(10,,t)==ci);
      if(numberof(coli)>2){
        xoz = regress(coms(1,coli),[coms(3,coli)^0,coms(3,coli)]);//x(z)
        yoz = regress(coms(2,coli),[coms(3,coli)^0,coms(3,coli)]);//y(z)
        cvec= [xoz(2),yoz(2),1.]/sqrt(1.+xoz(2)^2+yoz(2)^2);
        cavg= [xoz(1),yoz(1),0.];
        rcol= pbc( cavg - coms(:3,coli) );
        rcol= perp2( rcol, cvec );
        rcol= abs(rcol(1,),rcol(2,),rcol(3,));
        dots= dirt(+,coli)*cvec(+,);
        dots= acos(abs(dots));
        fast_bin,rcol,cr_arr,cf_arr;
        fast_bin,dots,th_arr,tf_arr;

        //build columnar order parameter for each column
        AA = array(0.,[2,3,3]);
        for(cj=1; cj<=numberof(coli); cj++){
          ui = dirs(,coli(cj));
          AA+= 3*ui*ui(-:1:3,)-[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
        }
        AA /= 2.*numberof(coli);
        eval = max(SVdec(AA));
        cS_arr(ci,t) = eval;
        cH_arr(ci,t) = HAT(coli,t)(avg);


        //nn separation parallel to colaxis
        zcoli = sort(CAT(3,coli,t));
        grow,zcoli,zcoli(1); //looped
        dz_col = CAT(:3,zcoli,t)(,dif);
        dz_col = para2( dz_col, cvec);
        dz_col = abs(dz_col(1,),dz_col(2,),dz_col(3,));
        fast_bin,dz_col,dz_arr,df_arr;


        //signed nearest neighbor phi difference
        nflip = sign(cos(CAT(4,zcoli,t)));
        nflip = nflip*roll(nflip,-1);
        d_phi = CAT(6,zcoli,t)(dif); //angle between forients def'd
        ww = where(abs(d_phi)>pi);   //positive for ccw rots of mol
        if(numberof(ww)){            //when dir is +z in lab frame.
          d_phi(ww) -= 2*pi*sign(d_phi(ww));
        }
        d_phi *= nflip(:-1);
        //POSITIVE D_PHI IS FROM 4F TOWARD 1F
        fast_bin,d_phi,ph_arr,pf2_arr;
        
      }
    }
  }

  write,format="\r Looped over %d ps",times;

  write,"~~~Saving Arrays~~~";
  qq=createb(dir+"cat_data.b");
  save,qq,gx_arr,gy_arr,gf_arr;
  save,qq,lx_arr,ly_arr,lf_arr;
  save,qq,oS_arr,cS_arr,oH_arr,cH_arr;
  save,qq,cr_arr,cf_arr,th_arr,tf_arr;
  save,qq,ph_arr,pf_arr,dz_arr,df_arr;
  save,qq,na_arr,nf_arr,lr_arr,lf2_arr;
  save,qq,ma_arr,mf_arr,pf2_arr;
  save,qq,zx_arr,zf_arr,xd_arr,xf_arr;
  close,qq;
  write,format="\rSaved arrays to %s\n",dir+"cat_data.b";
  
}






func le_cud(exec,dirs,lab=){
  for(i=1; i<=numberof(dirs); i++){
    if(exec(i)){
      dat = openb(dirs(i)+"cat.b");
      dog = unwrap_traj(dat.CAT(,,dat.eqt:),dat.DIM_ARR(,dat.eqt:));
      nit = ni_t(cat=dog);
      nit = is_void(lab)?!(!(nit*[0,0,1])):nit;
      sVDH,dog,nit*1.;
      label = "_vdh";
      if(!is_void(lab)) label+="LAB";
      label+=".b";
      dat = createb(dirs(i)+"cat"+label);
      save,dat,dz_arr,fz_arr,dr_arr,fr_arr;
      close,dat;
    }
  }
  //alle_cat(exec,dirs,lab=lab);
}
func alle_cat(exec,dirs,lab=){
  for(i=1; i<=numberof(dirs); i++){
    if(exec(i)){
      dat = openb(dirs(i)+"cat.b");
      dog = unwrap_traj(dat.CAT(,,EQTIMS(i):), dat.DIM_ARR(,EQTIMS(i):) );
      close,dat;
      label = "_dyn";
      if(!is_void(lab)) label+="LAB";
      label+=".b";
      dat = createb(dirs(i)+"cat"+label);
      dynacat,dog,lab=lab;
      save,dat,r2avg,r4avg,z2avg,z4avg;
      save,dat,r2var,r4var,z2var,z4var;
      write,"saved square displacements";
      close,dat;   
    }
  }
}




func detailizer(dirs,eqts){
  //accesses cat.b and uses dim_arr and eqts to determine the average number density over the run
  if(numberof(dirs)!=numberof(eqts)){
    write,format=" Err: %s\n","dirs don't match eqts";
    return "no op";
  }
  if(allof(strpart(dirs,0:)!="/")){
    ww=where(strpart(dirs,0:)!="/");
    dirs(ww) += "/";
  }
  for(d=1; d<=numberof(dirs); d++){
    eqt = eqts(d);
    qq = updateb(dirs(d)+"cat.b");
    dims = qq.DIM_ARR(,eqt:);
    vols =dims(1,)*dims(2,)*dims(3,);
    rhos = 600./vols;
    rhoo = rhos(avg);
    pstd = sqrt( avg( (rhos-rhoo)^2 ) );
    write,format=" %s : rho_o = %3g (%3g)\n",dirs(d),rhoo,pstd;
    save,qq,rhoo,eqt;
    close,qq;
  }
} 
    

  








func apres_ledeluge(void,dir=,outdir=,nosav=,reset=){
  out_arr = 0.;
  dir = is_void(dir)?"":dir;
  outdir = is_void(outdir)?dir:outdir;
  //outdir+= totxt(int(random(3)*99))(sum)+"_"
  dat = openb(dir+"cat_data.b");
  restore,dat;
  close,dat;
  if(!is_void(reset)) window,0,wait=1;

  
  //prepare g(perp,para) plot
  gy_arr-= gy_arr(1,1); //shifts bins so that zcen
  gx_arr-= gx_arr(1,1); //will give the bin center
  y_step = gy_arr(2,2);
  x_step = gx_arr(2,2); 

  write,format="%s","~~~Preparing g(perp,para)~~~";
  keep = gx_arr^2 + gy_arr^2 <= max(gx_arr)^2 ;
  dz = 2*gy_arr(:-1,dif);
  dA = pi*((gx_arr+x_step)^2)(dif);
  dV = dA*dz;
  VV = 4 * pi/3. * max(gy_arr+y_step)^3.;

  cw = window();
  window,cw+1,style="heatmap.gs",wait=1;
  fma; limits,e,e,e,e;
  palette,"heat.gp";
  prel = (gf_arr*keep)(:-1,:-1) / sum((gf_arr*keep)(:-1,:-1)*1.);
  pdens= prel/dV;
  prel = pdens*VV;
  pp   = where(pdens>0);
  Sgpp = -sum( (pdens*dV)(pp)*log((pdens*dV)(pp)) );
  out_arr = Sgpp;
  write,format="\nSgpp = %6.3g nats\n",Sgpp;

  
  levs = 0.; grow,levs,2.^(indgen(11)-6);
  plfc,prel,gy_arr(zcen,zcen),gx_arr(zcen,zcen),levs=levs;
  plc, prel,gy_arr(zcen,zcen),gx_arr(zcen,zcen),levs=[1.],marks=0,color="black";
  plxtitle,"r_perp (A)";
  plytitle,"r_para (A)";
  //build color bar - contours are spaced like 2^n and labelled as n
  plfc_levs(2:) = log(plfc_levs(2:))/log(2);
  plfc_levs(1) = -999;
  color_bar,vert=1,labs=[3,4];
  if(is_void(nosav)){
    write,format="\r~~~Saving as %s.png~~~",outdir+"gparaperp.png";
    png,outdir+"gparaperp";
    write,format="\rSaved to %s\n",outdir+"gparaperp.png";

    write,format="\rSaving 2D g(para,perp) data to %s...",outdir+"g_paraperp.asc";
    datfn = create(outdir+"g_paraperp.asc");
    write,datfn,format="%s %s %s\n","rperp(zcen)","rpara(zcen)","Prel";
    len = numberof(prel);
    write,datfn,format="%g %g %g\n",gx_arr(zcen,zcen),gy_arr(zcen,zcen),prel;
    close,datfn;
    write,format="\rSaved 2D g(para,perp) data to %s. \n",outdir+"g_paraperp.asc";
  }

  write,format="%s","\r~*~ Nearest Neighbor Values ~*~";

  //First peak parallel to director
  //checks cylinder: r<=5 and 3 < z <= 5
  xx = where(gx_arr(,1)<=5.);
  yy = where( (gy_arr(1,)>3.)*(gy_arr(1,)<=5.) );
  denom = sum((pdens*dV)(xx,yy));
  numer = (pdens*dV)(xx,yy)(sum,);
  probd = numer/denom;
  z_avg = sum(probd*gy_arr(1,yy));
  z_var = sum(probd*(gy_arr(1,yy)-z_avg)^2);
  z_std = sqrt(z_var);
  z_skw = sum(probd*(gy_arr(1,yy)-z_avg)^3)/z_std^3;
  write,format="\n  1st Intra Peak (std , skew) | %6.4g (%6.4g , %6.4g)",z_avg,z_std,z_skw;
  grow,out_arr,z_avf,z_std,z_skw;
  y_exc = z_avg;

  //First peak perp to director
  //checks cylinder: z<=5 and 5 < r <= 15
  xx = where( (gx_arr(,1)>5.)*(gy_arr(1,)<=15.));
  yy = where( (gy_arr(1,)<=5.) );
  denom = sum((pdens*dV)(xx,yy));
  numer = (pdens*dV)(xx,yy)(,sum);
  probd = numer/denom;
  z_avg = sum(probd*gx_arr(xx,1));
  z_var = sum(probd*(gx_arr(xx,1)-z_avg)^2);
  z_std = sqrt(z_var);
  z_skw = sum(probd*(gx_arr(xx,1)-z_avg)^3)/z_std^3;
  write,format="\n  1st Inter Peak (std , skew) | %6.4g (%6.4g , %6.4g)\n",z_avg,z_std,z_skw;
  grow,out_arr,z_avg,z_std,z_skw;
  x_exc = z_avg;

  exclude = where( (!gf_arr)*(gx_arr<x_exc)*(gy_arr<y_exc) );
  VolExcl = sum(dV(exclude));
  hc_para = gy_arr(1,where(gf_arr(1,)>=1)(1)); // || distance
  hc_perp = gx_arr(where(gf_arr(,1)>=1),1)(1); // -| distance
  ee_dist = (hc_perp-hc_para/2.)*2;            // -- distance
  write,format="  Hard Core ff dist ~ %6.4g\n  Hard Core fe dist ~ %6.4g\n  Hard Core ee dist ~ %6.4g\n  Excluded Core Vol ~ %6.4g A^3\n",hc_para,hc_perp,ee_dist,VolExcl;
  grow,out_arr,hc_para,hc_perp,ee_dist,VolExcl;
  yy = span(0.,hc_para,100);
  plotl,sqrt(1-(yy/hc_para)^2)*hc_perp,yy,color="magenta";
  plotl,sqrt(1-(yy/hc_para)^2)*ee_dist,yy,color="cyan";

  rcut = hc_perp/2.;
  zcut = hc_para/2.;

  write,format="%s","\r~~~Preparing g(perp)~~~";
  lims = limits();
  limits,lims(1),lims(2),lims(3),lims(4);
  plotl,lims(1:2),zcut*[1.,1.],color="blue";
  kp_cy = where(gy_arr(1,)+y_step < zcut);
  gperp = (gf_arr*keep)(:-1,kp_cy)(,sum);
  dz_cy = (gy_arr*keep)(:-1,kp_cy)(,max)*2;
  dA_cy = pi*((gx_arr+x_step)^2)(dif,1);
  dV_cy = dz_cy*dA_cy;
  gxarr = gx_arr(zcen,1);
  VV_cy = sum(dV_cy);
  VV_cy = 2*pi*zcut*(gx_arr(max,1)^2 - zcut*zcut/3.);
  prel = (gperp / (sum(gperp)*1.)) * (VV_cy/dV_cy);

  
  window,cw+2,style="boxed.gs",wait=1;
  plotl,lims(:2),[1.,1.],color="graya";
  plh,prel,gxarr,color="blue";
  plxtitle,"r (A)";
  plytitle,"Relative Density";
  window,cw+1,wait=1;

  if(is_void(nosav)){
    write,format="\rSaving 1D g perp data to %s...",outdir+"g_perp.asc";
    datfn = create(outdir+"g_perp.asc");
    write,datfn,format="%s, %s\n","rperp(zcen)","Prel";
    len = numberof(prel);
    write,datfn,format="%g, %g\n",gxarr,prel+1e-20;
    close,datfn;
    write,format="\rSaved 1D g perp data to %s. \n",outdir+"g_perp.asc";
  }
  
  write,format="%s","\r~~~Preparing g(para)~~~";
  plotl,rcut*[1.,1.],lims(3:4),color="blue";
  kp_cy = where(gx_arr(,1) < rcut);
  gpara = (gf_arr*keep)(kp_cy,:-1)(sum,);
  dz_cy = 2*gy_arr(1,dif);
  dr_cy = (gx_arr*keep)(kp_cy,:-1)(max,) + gx_arr(2,2);
  dA_cy = pi*dr_cy*dr_cy;
  dV_cy = dz_cy*dA_cy;
  gyarr = gy_arr(1,zcen);
  VV_cy = sum(dV_cy);
  rrad = gx_arr(max,1);
  zpiv = sqrt(rrad^2 - rcut^2 );
  VV_cy = 4*pi/3 *(rrad^3 - zpiv*rrad^2 + zpiv*rcut^2)
  prel = (gpara/(1.*sum(gpara))) * (VV_cy/dV_cy);

  window,cw+2,wait=1;
  plh,prel,gyarr,color="red";
  plxtitle,"r (A)";
  plytitle,"Relative Density";

  if(is_void(nosav)){
    write,format="\rSaving 1D g para data to %s...",outdir+"g_para.asc";
    datfn = create(outdir+"g_para.asc");
    write,datfn,format="%s, %s\n","rpara(zcen)","Prel";
    len = numberof(prel);
    write,datfn,format="%g, %g\n",gyarr,prel+1e-20;
    close,datfn;
    write,format="\rSaved 1D g para data to %s. \n",outdir+"g_para.asc";
  }

  write,"\r**For some reason the g curves aren't normalizing** \n";

  if(anyof(lf_arr)){
    write,format="%s","~~~Preparing orientation heatmap~~~";
    window,cw+3,style="heatmap.gs",wait=1;
    palette,"earth.gp";
    dlx = lx_arr(dif,1)(1);
    dly = ly_arr(1,dif)(1);
    prel = (lf_arr/(sum(lf_arr(:-1,:-1))*1.)) * (4*pi*pi / (dlx*dly) );
    plfc,prel(:-1,:-1),180/pi*ly_arr(zcen,zcen),180/pi*lx_arr(zcen,zcen),levs=(2.^(indgen(8)-4));
    plc, prel(:-1,:-1),180/pi*ly_arr(zcen,zcen),180/pi*lx_arr(zcen,zcen),levs=[1.],color="black",marks=0;
  
    plfc_levs = log(plfc_levs)/log(2);
    color_bar,vert=1,labs=[1,3];
 

    //Plotting guide lines. Monofluoro to port!
    plotl,[-180, 180],[ 120, 120],width=5,color="black";  //neighbor no-F
    plotl,[ 120, 120],[-180, 180],width=5,color="black";  //reference no-F
    plotl,[-180, 180],[ 120, 120],width=2,color="green";
    plotl,[ 120, 120],[-180, 180],width=2,color="green";
    plotl,[-180, 180],[-120,-120],width=5,color="black";  //neighbor monoF
    plotl,[-120,-120],[-180, 180],width=5,color="black";  //reference monoF
    plotl,[-180, 180],[-120,-120],width=2,color="red";
    plotl,[-120,-120],[-180, 180],width=2,color="red";
    plotl,[-180, 180],[   0,   0],width=5,color="black";  //neighbor perF
    plotl,[   0,   0],[-180, 180],width=5,color="black";  //reference perF
    plotl,[-180, 180],[   0,   0],width=2,color="grayc";
    plotl,[   0,   0],[-180, 180],width=2,color="grayc";

    if(is_void(nosav)){
      
      write,format="\r~~~Saving as %s.png~~~",outdir+"lambdalambda.png";
      png,outdir+"lambdalambda.png";
      write,format="\rSaved to %s\n",outdir+"lambdalambda.png";

      write,format="\rSaving 2D lambda lamda data to %s...",outdir+"lamlam.asc";
      datfn = create(outdir+"lambdalambda.asc");
      write,datfn,format="%s %s %s\n","lambda1(rad)","lambda2(rad)","Prel";
      len = numberof(prel);
      write,datfn,format="%g %g %g\n",ly_arr,lx_arr,prel+1e-20;
      close,datfn;
      write,format="\rSaved 2D lambda lambda data to %s. \n",outdir+"g_para.asc";
    }
    write,format="\r%s\n","Orientational Data:";
    pdens = prel/(4*pi*pi);
    write,format=" Integrates to %g\n",sum(pdens*dlx*dly);
    pdens = pdens(where(pdens));
    write,format=" S = %6.4g nats\n",-sum(pdens*dlx*dly*log(pdens*dlx*dly));
    pair_integrate; //see below 
  }

  write,format="\n%s","~~~Preparing order param analysis~~~";

  if(is_void(nosav)*noneof(oS_arr)){
    cc = openb(dir+"cat.b");
    oS_arr = oriS_cat(cat=cc.CAT(,,-numberof(oS_arrr):));
    close,cc;
  }
  avgS = oS_arr(avg); avgH = oH_arr(avg);
  stdS = sqrt(avg( (oS_arr - avgS)^2 ) );
  stdH = sqrt(avg( (oH_arr - avgH)^2 ) );
  write,format="\r%s\n"," Col ##           avgS (stdS) | avgH (stdH)";
  write,format=  "  Total  S= %-6.3g ( %6.3g ) | %-6.3g ( %6.3g )\n",avgS,stdS,avgH,stdH;
  
  avgSc = cS_arr(,avg);
  stdSc = sqrt(avg( (cS_arr - avgSc)^2));
  avgHc = cH_arr(,avg);
  stdHc = sqrt(avg( (cH_arr - avgHc)^2));
  write,format="Col Avg  S= %-6.3g ( %-6.3g ) | %-6.3g ( %-6.3g )\n",avgSc(avg),stdSc,avgHc(avg),stdHc;
  

  //orientation analysis
  window,cw+4,wait=1,style="bare.gs";
  palette,"yarg.gp";
  limits,-15,15,-15,15;
  //xxx
  //parallel neighbor phorientation
  dphi = ph_arr(dif)(1);
  pden = pf_arr(:-1)/sum(pf_arr*dphi);
  prad = sqrt(2*pden);
  cphi = cos(ph_arr(zcen));
  sphi = sin(ph_arr(zcen));

  pden = pden(where(pden));
  sent = sum( pden*dphi*log(pden*dphi) );
  write,format="%s\n","Nearest parallel neighbor values"; 
  write,format="   PhiS: %6.3g nats\n",-1.*sent;
  write,format="   Peak dphi: %6.3g\n",ph_arr(zcen)(where(prad==max(prad)));

  write,format="%s\n","Nearest perpular neighbor values";
  //perp neighbor distance vs orientation
  prel = lf2_arr(,:-1); prel(0,)=lf2_arr(1,:-1);
  dphi = lx_arr(dif,)(1);
  prel/= sum(lf2_arr*1.);
  xcir = (lx_arr - dphi/2.)(,:-1); //zcen'd and looped
  drr  = (lr_arr^2)(,dif);
  dA   = dphi*drr;
  prel = prel * (pi*(15^2 - 5^2))/dA;
  plfc,prel,sin(xcir)*lr_arr(,zcen),cos(xcir)*lr_arr(,zcen),levs=2^(indgen(9)-4.);
//plc ,prel,sin(xcir)*lr_arr(,zcen),cos(xcir)*lr_arr(,zcen),levs=[1.0],color="black",marks=0;
  pp = createb("test.b");
  save,pp,prel,xcir,lr_arr;
  close,pp;

  arb = 5.0; //scale value to multiply pdens by for plotting
  plotl,arb*prad*cphi,arb*prad*sphi,width=3,color="blue";
  plotl,arb*cphi/sqrt(pi),arb*sphi/sqrt(pi),width=2,color="grayd";



  pdens = prel / (pi*(15^2 - 5^2));
  rmsex = sqrt( (pdens*dA*lr_arr(,zcen)^2)(,sum) / (pdens*dA)(,sum));
  plotl,rmsex*cos(xcir(,1)),rmsex*sin(xcir(,1)),color="black",width=6;
  plotl,rmsex*cos(xcir(,1)),rmsex*sin(xcir(,1)),color="magenta",width=3;
  plt,"{r_perp_: 5 - 15 A}",0,-13.75,tosys=1,justify=9;

  rmin = min(rmsex);
  rmat = xcir(where(rmsex==rmin),1)*180/pi;
  write,format="   min <r^2>: %6.4g at %6.4g deg\n",rmin,rmat; 
  rmax = max(rmsex);
  rmat = xcir(where(rmsex==rmax),1)*180/pi;
  write,format="   max <r^2>: %6.4g at %6.4g deg\n",rmax,rmat;
  
  if(is_void(nosav)){
    write,format="\r~~~Saving as %s.png~~~",outdir+"lambdaphi.png";
    png,outdir+"lambdaphi.png";
    write,format="\rSaved to %s\n",outdir+"lambdaphi.png";

    write,format="\r~~~Saving as %s~~~",outdir+"rnnexp.asc";
    datfn = create(outdir+"rnnexp.asc");
    write,datfn,format="%s, %s\n","acos(Rij.f)","r.exp";
    write,datfn,format="%g, %g \n",xcir(,1)*180/pi,rmsex;
    close,datfn;
    write,format="Saved as %s\n",outdir+"rnnexp.asc";
  }

  //plot grid lines
  for(n=1; n<12; n++){
    plotl,[5,15]*cos(n*pi/6.),[5,15]*sin(n*pi/6.),color="black";
  }
  for(r=7.5; r<15; r+=2.5){
    plotlc,0,0,r,color="black";
  }
  plotl,15*[0,cos(4*pi/3)],15*[0,sin(4*pi/3)],width=6,color="black";
  plotl,15*[0,cos(2*pi/3)],15*[0,sin(2*pi/3)],width=6,color="black";
  plotl,15*[0,cos(0*pi/3)],15*[0,sin(0*pi/3)],width=6,color="black";
  plotl,15*[0,cos(4*pi/3)],15*[0,sin(4*pi/3)],width=3,color="green";
  plotl,15*[0,cos(2*pi/3)],15*[0,sin(2*pi/3)],width=3,color="red";
  plotl,15*[0,cos(0*pi/3)],15*[0,sin(0*pi/3)],width=3,color="grayd";

  //skip column analysis when over 80% of CAT are orphans
  if( numberof(where(cS_arr>0))*5>numberof(cS_arr) ){

    //Column analysis
    write,format="\n%s\n","~~~Column Analysis~~~";
    window,cw+5,wait=1;
    
    //Analyzing perp displacement from colaxis
    prel = cf_arr(:-1)/sum(cf_arr(:-1)*1.);
    dr = 0.; grow,dr,cr_arr;
    dA = pi*(dr^2)(dif)(:-1); 
    AA = pi*dr(0)^2;
    prel *= (AA/dA);
    pdens = prel/AA;

    rpeak = where(pdens==max(pdens))(1);
    r_exp = sum(pdens*dA*cr_arr(zcen));
    r_var = sum(pdens*dA*cr_arr(zcen)^2);
    write,format="%s\n","Displacements perp to colaxis";
    write,format="  Peak @ bin # %d\n  <r>   = %6.3g A\n  <r^2> =(%6.3g A)^2\n  sig   = %6.3g A\n",rpeak,r_exp,sqrt(r_var),sqrt(abs(r_exp^2 - r_var ));
    return sqrt(r_var); //deletethis
    
    if(is_void(nosav)){
      write,format="\rSaving column axis displacement data to %s...",outdir+"column_perpD.asc";
      datfn = create(outdir+"column_perpD.asc");
      write,datfn,format="%s, %s, %s\n","rperp(zcen)","Prel","Count";
      write,datfn,format="%g, %g, %g\n",cr_arr(zcen),prel+1e-20,cf_arr(:-1);
      close,datfn;
      write,format="\rSaved column axis displacement data to %s. \n",outdir+"column_perpD.asc";
    }

    //tilt analysis
    dA = 2*2*pi*sin(th_arr(zcen))*th_arr(dif)(1);
    degs = th_arr(zcen)*180/pi;
    pdens = (tf_arr(:-1)/sum(tf_arr(:-1)*1.));
    pdens*= (1./dA);
    pp = where(pdens>0);
    Sent = -sum((pdens*dA)(pp)*log((pdens*dA)(pp)));
    th_pk = where(pdens==max(pdens))(1);
    thexp = sum(pdens*dA*degs);
    thvar = sum(pdens*dA*degs^2);

    window,0,wait=1;
    prad = sqrt(pdens*2);
    colr = int(random(3)*255)+1;
    plotl,prad*sin(th_arr(zcen)),prad*cos(th_arr(zcen)),color=colr;
    plotl,-prad*sin(th_arr(zcen)),prad*cos(th_arr(zcen)),color=colr;
    window,cw+5,wait=1;
    
    write,format="%s\n  S_th  = %6.3g nats\n","Angle off colaxis:",Sent;
    write,format="  Peak @ bin # %d\n  <th>  = %6.3g deg\n  <th^2>=(%6.3g deg)^2\n  sig   = %6.3g deg\n",th_pk,thexp,sqrt(thvar),sqrt(abs(thexp^2 - thvar ));

    smove = filter_bessel(th_arr(zcen),pdens,2e2);
    peaks = where( (sign(smove(dif))(dif)<0)*(tf_arr(:-3)>10) );
    write,format="  Detected %d potential peaks\n",numberof(peaks);
    write,format="   %6.3g deg\n",th_arr(zcen)(peaks)*180/pi;
    
    plh,pdens*sum(dA),degs;
    plxtitle,"Intracolumn Tilt (degs)";
    plytitle,"P_rel";


    if(is_void(nosav)){
      write,format="\rSaving column axis tilt data to %s...",outdir+"column_tilt.asc";
      datfn = create(outdir+"column_tilt.asc");
      write,datfn,format="%s, %s, %s\n","angle(deg) Pdens Count";
      write,datfn,format="%g, %g, %g\n",degs,pdens+1e-20,tf_arr(:-1);
      close,datfn;
      write,format="\rSaved column axis tilt data to %s. \n",outdir+"column_tilt.asc";
    }    
    
  } else write,format="\n%s\n","Skipped column specific analysis.";


  
  
  window,cw,wait=1;


}
