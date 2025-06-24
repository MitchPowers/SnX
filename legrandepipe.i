#include "~/TinkerExp/sn3.i"
#include "~/TinkerExp/whisker.i"
/*

  Le Grande Donéeoduc

  built to interact with the output of the sn3 function parse_traj

  Calculates static quantities as the traj is parsed

  Takes prefix as arg,
  all else is external
  Only return 0 or err
  (err is an int flag)

  Assumes:
  CoM is molecular com arr
  TlA is theta phi fhi arr
  CoR is naive molecular cor arr
  DIMS holds fresh DIMS

*/
extern EQTIMS;
QQ = 0;
read,prompt="\nF@ ",format="%d",QQ;
if(QQ==8)/*F8*/EQTIMS = [5000,5000,6000,4500,2000];
if(QQ==7)/*F7*/EQTIMS = [15000,4000,8000,5000,6000,2000];
if(QQ==6)/*F6*/EQTIMS = [10000,2500,5000,5000,2500];
if(QQ==5)/*F5*/EQTIMS = [4000,5000,3000];


func avant_donneeoduc(void,eqts=,skip=,lebin=,outdirs=){
  lebin = is_void(lebin)?faire_lebin():lebin;
  pres = lebin(,1);
  dirs = lebin(,2);
  outdirs = is_void(outdirs)?dirs:outdirs;
  if(is_void(skip)){
    for(n=1; n<=numberof(dirs); n++){
      faire_lechat(pres(n),dirs(n),outdir=outdirs(n));
    }
    write,"Finished herding cats.";
  }
  if(is_void(eqts)){
    write," ~>Determine eq times, then run donneeduc";
    return [pres,dirs];
  }
  if(numberof(eqts)==numberof(pres)){
    donneeoduc(pres,dirs,eqts);
  }
    
}





func faire_lebin(void){
  write,format="USER SHOULD EDIT legrandepipe.i\n Currently set dir ~ %s\n                fn ~ %s\n","F*res*/","F*res*_free";
  dirs = "./"+ls(" F*res* -d")+"/";
  pres = strpart(ls(" F*res* -d"),:6)+"_free";
  return [pres,dirs];
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
  extern cS_arr; //order param of single columns             //coded
  cS_arr = array(0.,[2,NCOL,NCF]);
  extern oH_arr; //hexatic order x arrays                    //coded
  oH_arr = array(0.,NCF);
  extern cH_arr; //hexatic param of single columns           //coded
  cH_arr = array(0.,[2,NCOL,NCF]);
  extern HAT;    //hex OP of each molecule at all times      //coded
  HAT = array(0.,[2,MOLS,NCF]);
 

  write,format="Parsing traj: %s\n",prefix;
  oS_arr = parse_traj(prefix,dir=dir,funct=oriS);
  write,format="\rParsed traj: %s\n",prefix;

  
  winkill,0;
  window,0,style="work.gs",wait=1; fma;
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
    write,format="    %d : %d | %g  ( %g )\n",alph,omeg,Savg,Sstd;
  }

  
  DOG = unwrap_traj(CAT,DIM_ARR);


  write,"~~~Saving CAT/DOG/DIMS~~~";
  qq=createb(outdir+"cat.b");
  save,qq,CAT,DIM_ARR,DOG;
  close,qq;
  write,format="\rSaved CAT to %s\n",outdir+"cat.b";

  
  write,"~~~Saving Arrays~~~";
  qq=createb(outdir+"cat_data.b");
  save,qq,oS_arr,cS_arr,oH_arr,cH_arr,HAT;
  write,format="\rSaved arrays to %s\n",outdir+"cat_data.b";

  nit = ni_t();

  CAT(10,,) = column_tracker(nit=nit);
  NCF = dimsof(CAT)(0);
  ii = indgen(MOLS)(,-:1:NCF);
  tt = indgen(NCF)(-:1:MOLS,);
  qq = updateb(outdir+"cat.b");
  save,qq,CAT; close,qq;
  write,format="\rTracked columns. %s","Preparing graphic...\n";


  
  winkill,1;
  window,1,height=MOLS,wait=1,style="nobox.gs";//"bare.gs";
  fma;
  palette,"stern.gp";
  pli,transpose(CAT(10,,));
  plotl,span(0.,NCF,NCF),oS_arr*MOLS,color="red";
  plxtitle,"Free Time (ps)";
  limits,e,e,e,e;
  png,outdir+prefix+"colcolor";
  for(i=1; i<=NCF; i++) plot,i,oS_arr(i:)(avg)*MOLS,sym='\1',color="grayd";
  
  eq = where( CAT(10,sum,)==0 );
  if(numberof(eq)>=1) write,format="\rLe temp isotrope est plus que %d\n",eq(1);
  else write,"\rBon chance avec cette.";


  
  write,"~~~Saving Column Plot values~~~";
  qq=create(outdir+"Column_plot.asc");
  write,qq,"t i col";
  len = numberof(tt);
  write,qq,format="%d %d %g\n",tt(:len),ii(:len),CAT(10,,)(:len);
  close,qq;
  write,format="\rSaved column data to %s\n",outdir+"Column_plot.asc";
}







func chatoduc(exec,pres,dirs,eqts,stops=){
  //exec is an array of 0 or 1, for wether to run that chat in leduc
  ww=where(exec==1);
  stops = is_void(stops)?array(0,dimsof(eqts)):stops;
  write,format="Nous vollons les CAT en la duc de %s\n",pres(ww);
  donneeoduc,pres(ww),dirs(ww),eqts(ww),stops=stops(ww);
  return "Ferme.";
}
func donneeoduc(pres,dirs,eqts,stops=){
  /*DOCUMENT donneeoduc(pres,dirs,eqts) */
  stops = is_void(stops)?eqts*0:stops;
  for(n=1; n<=numberof(pres); n++){
    write,format="\n Piping %s at %s starting from %d\n",pres(n),dirs(n),eqts(n);
    //faire_lechat(pres(n),dirs(n));
    le_duc,pres(n),dirs(n),eqts(n),stops=stops(n);
    apres_ledeluge,dir=dirs(n),nosav=1;
  }
}










func le_duc(prefix,dir,eqtime,stops=){

  extern CAT;
  extern DIM_ARR;
  CAT=DIM_ARR=[]; //kill any current CAT
  write,format="Opening %s for CAT\n",dir+"cat.b";
  dat = openb(dir+"cat.b");
  if(dat.CAT==[]){
    write,"Did not find CAT.";
    return -999;
  }
  if(eqtime > dimsof(dat.CAT)(0)){
    write,"ERR: Bad Time."
      return -998;
  }
  stops = is_void(stops)?0:stops;
  if((stops==0)+(stops>dimsof(dat.CAT)(0))) stops = dimsof(dat.CAT)(0);
  
  //trim relevant arrs
  restore,dat;
  CAT = CAT(,,eqtime:stops);
  DIM_ARR = DIM_ARR(,eqtime:stops);
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
  extern HAT
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
      dog = unwrap_traj(dat.CAT(,,EQTIMS(i):),dat.DIM_ARR(,EQTIMS(i):));
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
    write,datfn,format="%s %s\n","rperp(zcen)","Prel";
    len = numberof(prel);
    write,datfn,format="%g %g\n",gxarr,prel+1e-20;
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
    write,datfn,format="%s %s\n","rpara(zcen)","Prel";
    len = numberof(prel);
    write,datfn,format="%g %g\n",gyarr,prel+1e-20;
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
    write,datfn,format="%s\n","acos(Rij.f) r.exp";
    write,datfn,format="%g %g \n",xcir(,1)*180/pi,rmsex;
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
      write,datfn,format="%s %s %s\n","rperp(zcen)","Prel","Count";
      write,datfn,format="%g %g %g\n",cr_arr(zcen),prel+1e-20,cf_arr(:-1);
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
      write,datfn,format="%s %s %s\n","angle(deg) Pdens Count";
      write,datfn,format="%g %g %g\n",degs,pdens+1e-20,tf_arr(:-1);
      close,datfn;
      write,format="\rSaved column axis tilt data to %s. \n",outdir+"column_tilt.asc";
    }    
    
  } else write,format="\n%s\n","Skipped column specific analysis.";


  
  
  window,cw,wait=1;


}

  
  


func pair_integrate(void){
  //works with lambda plot externs presently in memory
  tmpx=lx_arr(zcen,zcen);
  tmpy=ly_arr(zcen,zcen);
  tmpf=lf_arr(:-1,:-1);
  dxdy=tmpx(dif,)(1)^2;

  write,format="\n%s\n","60deg squares";
  //vert-vert
  sieve = (cos(3*tmpx)>0)*(cos(3*tmpy)>0);
  AAwide = sum(sieve*tmpf)/sum(tmpf*1.);
  write,format=" >< : %.3g  |",AAwide*100;
 
  //face-face
  sieve = (cos(3*tmpx)<0)*(cos(3*tmpy)<0);
  BBwide = sum(sieve*tmpf)/sum(tmpf*1.);
  write,format=" <> : %.3g\n",BBwide*100;

  //face-vert
  sieve = (cos(3*tmpx)>0)*(cos(3*tmpy)<0);
  sieve+= (cos(3*tmpx)<0)*(cos(3*tmpy)>0);
  ABwide = sum(sieve*tmpf)/sum(tmpf*1.);

  write,format="%s\n",array("-",32)(sum);
  write,format="Homoiety  : %.3g%%\n",(AAwide+BBwide)*100;
  write,format="Heteroiety: %.3g%%\n",ABwide*100;

  
  write,format="\n%s\n","30deg squares";
  //sharp ><
  sieve = (cos(3*tmpx)>0)*(cos(6*tmpx)>0)*(cos(3*tmpy)>0)*(cos(6*tmpy)>0);
  AApri = sum(sieve*tmpf)/sum(tmpf*1.);
  write,format=" >< : %g  |",AApri;

  //sharp <>
  eveis = sieve;
  sieve = (cos(3*tmpx)<0)*(cos(6*tmpx)>0)*(cos(3*tmpy)<0)*(cos(6*tmpy)>0);
  eveis+= sieve;
  BBpri = sum(sieve*tmpf)/sum(tmpf*1.);
  write,format=" <> : %g\n",BBpri;

  //sharp >>/<<
  sieve = (cos(3*tmpx)<0)*(cos(6*tmpx)>0)*(cos(3*tmpy)>0)*(cos(6*tmpy)>0);
  sieve+= (cos(3*tmpx)>0)*(cos(6*tmpx)>0)*(cos(3*tmpy)<0)*(cos(6*tmpy)>0);
  eveis+= sieve;
  ABpri = sum(sieve*tmpf)/sum(tmpf*1.);
  write,format=" >> : %g\n",ABpri;
       

  
  //offpeak
  sieve = (cos(6*tmpx)<=0)*(cos(6*tmpy)>=0);
  sieve+= (cos(6*tmpx)>=0)*(cos(6*tmpy)<=0);
  offpk = sum(sieve*tmpf)/sum(tmpf*1.);

  offoff = 1. - sum((sieve+eveis)*tmpf)/sum(tmpf*1.);
  
  write,format="%s\n",array("-",32)(sum);
  write,format="Homoiety  : %.3g%%\n",(AApri+BBpri)*100;
  write,format="Heteroiety: %.3g%%\n",ABpri*100;
  write,format="PeakPeak  : %.3g%%\n",(AApri+ABpri+BBpri)*100;
  write,format="Peak-off  : %.3g%%\n",offpk*100;
  write,format=" Off-Off  : %.3g%%\n",offoff*100;

  
  //plf,(sieve+eveis)*tmpf/sum(tmpf*dxdy),tmpy,tmpx;
   
}

func double_lambda(cat_dat1,cat_dat2){  
    palette,"earth.gp";
    fn = openb(cat_dat1);
    lx_arr = fn.lx_arr;
    ly_arr = fn.ly_arr;
    lf_arr = fn.lf_arr;
    close,fn;
    dlx = lx_arr(dif,1)(1);
    dly = ly_arr(1,dif)(1);
    prel_1 = (lf_arr(:-1,:-1)/(sum(lf_arr(:-1,:-1))*1.)) * (4*pi*pi / (dlx*dly) );
    lxarr1 = 180/pi * lx_arr(zcen,zcen);
    lyarr1 = 180/pi * ly_arr(zcen,zcen);
    keep_1 = where(lxarr1>lyarr1);
    
    fn = openb(cat_dat2);
    lx_arr = fn.lx_arr;
    ly_arr = fn.ly_arr;
    lf_arr = fn.lf_arr;
    close,fn;
    dlx = lx_arr(dif,1)(1);
    dly = ly_arr(1,dif)(1);
    prel_2 = (lf_arr(:-1,:-1)/(sum(lf_arr(:-1,:-1))*1.)) * (4*pi*pi / (dlx*dly) );
    lxarr2 = 180/pi * lx_arr(zcen,zcen);
    lyarr2 = 180/pi * ly_arr(zcen,zcen);
    keep_2 = where(lxarr1<lyarr1);

    for(n=1; n<max(dimsof(lxarr1)(0),dimsof(lxarr2)(0)); n++){
      if(n<dimsof(lxarr1)(2)){
        plfc,prel_1(n:,n:n+1),lyarr1(n:,n:n+1),lxarr1(n:,n:n+1),levs=2.^(indgen(7)-3);
        plc,prel_1(n:,n:n+1),lyarr1(n:,n:n+1),lxarr1(n:,n:n+1),levs=[1.],marks=0,color="black";
      }
      if(n<dimsof(lxarr2)(2)){
        plfc,prel_2(n:n+1,n:),lyarr2(n:n+1,n:),lxarr2(n:n+1,n:),levs=2.^(indgen(7)-3);
        plc,prel_2(n:n+1,n:),lyarr2(n:n+1,n:),lxarr2(n:n+1,n:),levs=[1.],marks=0,color="black";
      }
    }
    plotl,[-180,180],[-180,180],color="black";
    
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

}

func triple_gpp(cat_dat1,cat_dat2,cat_dat3,nocon=){
  lims = [0.,0.,0.,0.];
  fma;
  
  bx = span(-50,50,129)(,-:1:129);
  by = span(-50,50,129)(-:1:129,);
  ff = ((bx>=0)*(by<=0));
  plf,ff,by,bx;

  if(cat_dat1!="SKIP"){
    fn = openb(cat_dat1);
    restore,fn;
    lims(2) = max(gx_arr);
    lims(4) = max(gy_arr);
    plot_gpp,nocon=nocon;
    plt,"r_perp_(A)",15.0,-3.75,tosys=1,justify="CH";
  }

  if(cat_dat2!="SKIP"){
    fn = openb(cat_dat2);
    restore,fn;
    gx_arr*=-1;
    lims(1) = min(gx_arr);
    lims(4) = max(max(gy_arr),lims(4));
    plot_gpp,nocon=nocon;
  }

  if(cat_dat3!="SKIP"){
    fn = openb(cat_dat3);
    restore,fn;
    gx_arr*=-1;
    gy_arr*=-1;
    lims(1) = min(min(gx_arr),lims(1));
    lims(3) = min(gy_arr);
    plot_gpp,nocon=nocon;
    plt,"r_para_ (A)",3.75,-15.0,tosys=1,justify="CH",orient=1;
  }
  
  plfc_levs(2:) = log(plfc_levs(2:))/log(2);
  plfc_levs(1) = -999;
  color_bar,vert=1,labs=[3,4];
  
  limits,lims(1),lims(2),lims(3),lims(4);
  plotl,[0.,0.],lims(3:4),color="grayd",width=2;
  plotl,lims(1:2),[0.,0.],color="grayd",width=2;
  for(i=1; i<=50; i++){
    len = 0.25;
    len += (i%5==0)*len;
    len += (i%10==0)*len;
    plotl,i*[1,1],len*[-1,1],color="grayd",width=2;
    plotl,i*[-1,-1],len*[-1,1],color="grayd",width=2;
    plotl,len*[-1,1],i*[1,1],color="grayd",width=2;
    plotl,len*[-1,1],i*[-1,-1],color="grayd",width=2;
  }
  plt,"10",10,-1.25,justify="CT",tosys=1;
  plt,"10",1.25,-10,justify="LH",tosys=1;
  plt,"20",20,-1.25,justify="CT",tosys=1;
  plt,"20",1.25,-20,justify="LH",tosys=1;
}
func plot_gpp(void,nocon=){
  //uses current external vars

  gy_arr -= gy_arr(1,1);//*1.5; //shifts bins so zcen
  gx_arr -= gx_arr(1,1);//*1.5; //will give bin center

  keep = gx_arr^2 + gy_arr^2 < max(gx_arr^2) ;
  dz = abs(2*gy_arr(:-1,dif));
  dA = pi*(gx_arr^2)(dif);
  dV = dA*dz;
  VV = 4 * pi/3. * max(abs(gy_arr))^3.;

  prel = (gf_arr*keep)(:-1,:-1) / sum((gf_arr*keep)(:-1,:-1)*1.);
  pdens= prel/dV;
  prel = pdens*VV;
  
  levs = 0.; grow,levs,2.^(indgen(11)-6);
  plfc,prel,gy_arr(zcen,zcen),gx_arr(zcen,zcen),levs=levs;
  if(is_void(nocon)) plc, prel,gy_arr(zcen,zcen),gx_arr(zcen,zcen),levs=[1.],marks=0,color="black";
}
 
  
func output_svdh(indirs,outdirs=,outlabs=){
  outdirs = is_void(outdirs)?indirs:outdirs;
  if( numberof(outdirs)==1 ) outdirs=outdirs(-:1:numberof(indirs));
  dumbout = totxt(int(random(3,numberof(indirs))*99))(sum,)+"_";
  outlabs = is_void(outlabs)?dumbout:outlabs; //avoid clobber
  for(i=1; i<=numberof(indirs); i++){
    window,20+i,wait=1; fma;
    inf = openb(indirs(i)+"cat_vdh.b");
    restore,inf;
    tmax= dimsof(fr_arr)(0);
    noom= 10^(indgen(int(log(tmax)/log(10))+1)-1);
    half= 0;
    if(tmax>=noom(0)*5 ){
      grow,noom,int(noom(0)*sqrt(10));
      half = 1;
    }
    write,format="Plotting %d ps\n",noom;
    drcen = dr_arr - dr_arr(1)/2.;
    pdens = fr_arr(,noom)/(1.*fr_arr(sum,noom));
    dA = pi*(drcen(1)*2)^2; grow,dA,pi*(dr_arr^2)(dif);
    AAperp = sum(dA);
    pdens /= dA(,-:1:numberof(noom));
    prel_perp = pdens * AAperp;
    prel_perp += 1e-20;
    //plotl,drcen,prel_perp;

    dzcen = dz_arr - dz_arr(1)/2.;
    pdens = fz_arr(,noom)/(1.*fz_arr(sum,noom));
    dA = (2*dzcen(1))(-:1:numberof(dzcen));
    AApara = sum(dA);
    pdens /= dA(1);
    prel_para = pdens * AApara;
    prel_para += 1e-20;
    //plotl,dzcen,prel_para;

    write,format=" Writing to %s\n",outdirs(i)+outlabs(i)+"sVDH.asc";
    ouf = create(outdirs(i)+outlabs(i)+"sVDH.asc");
    write,ouf,format="Relative densities - prel factors : %g , %g\n",AAperp,AApara;

    head = "perp_disp 1ps 10ps";
    line = "para_disp 1ps 10ps";
    frmt = "%g %g %g";
    if(max(noom)>=100){
      head += " 100ps"; line += " 100ps"; frmt += " %g";
    }
    if(max(noom)>=1000){
      head += " 1ns";   line += " 1ns";    frmt += " %g";
    }
    if(max(noom)>=10000){
      head += " 10ns";  line += " 10ns";   frmt += " %g";
    }
    if(half){
      head += " +.5";  line += " +.5";   frmt += " %g";
    }
    headline = head+" "+line+"\n";
    formatte = frmt+" "+frmt;

    outs = dimsof(prel_perp)(0);
    write,ouf,format="%s",headline;
    if(outs==2)  write,ouf,format=formatte,drcen,prel_perp(,1),prel_perp(,2),dzcen,prel_para(,1),prel_para(,2);
    if(outs==3)  write,ouf,format=formatte,drcen,prel_perp(,1),prel_perp(,2),prel_perp(,3),dzcen,prel_para(,1),prel_para(,2),prel_para(,3);
    if(outs==4)  write,ouf,format=formatte,drcen,prel_perp(,1),prel_perp(,2),prel_perp(,3),prel_perp(,4),dzcen,prel_para(,1),prel_para(,2),prel_para(,3),prel_para(,4);
    if(outs==5)  write,ouf,format=formatte,drcen,prel_perp(,1),prel_perp(,2),prel_perp(,3),prel_perp(,4),prel_perp(,5),dzcen,prel_para(,1),prel_para(,2),prel_para(,3),prel_para(,4),prel_para(,5);
    if(outs==6)  write,ouf,format=formatte,drcen,prel_perp(,1),prel_perp(,2),prel_perp(,3),prel_perp(,4),prel_perp(,5),prel_perp(6),dzcen,prel_para(,1),prel_para(,2),prel_para(,3),prel_para(,4),prel_para(,5),prel_para(,6);
    for(p=1; p<=outs; p++) plotl,drcen,prel_perp(,p);
    logxy,0,1;
    limits,0,e,1e-3,e;
    close,ouf;

    sVDH_gauss,outdir=outdirs(i)+"GaussFits_"+outlabs(i);
  }
}
        
    


func asc_vdh(indirs,outdirs){
  for(i=1; i<=numberof(indirs); i++){
    fin = openb(indirs(i));
    restore,fin;
    dimr = dimsof(fr_arr);
    dimz = dimsof(fz_arr);
    rt_arr = dr_arr(,-:1:dimr(0));
    tr_arr = indgen(dimr(0))(-:1:dimr(2),);
    zt_arr = dz_arr(,-:1:dimz(0));
    tz_arr = indgen(dimz(0))(-:1:dimz(2),);

    fot = create(outdirs(i)+"_sVDH_perp.asc");
    write,fot,format="%s\n","r(t) (A), t (ps), counts";
    len = numberof(fr_arr);
    write,fot,format="%g, %d, %g\n",rt_arr(:len),tr_arr(:len),fr_arr(:len);
    close,fot;
    fot = create(outdirs(i)+"_sVDH_para.asc");
    write,fot,format="%s\n","z(t) (A), t (ps), counts";
    len = numberof(fz_arr);
    write,fot,format="%g, %d, %g\n",zt_arr(:len),tz_arr(:len),fz_arr(:len);
    close,fot;
  }
}

    
      
