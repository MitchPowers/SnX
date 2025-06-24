/*****************************************************************
 *
 *  MITCH'S TINKER ANALYSIS SUITE - FUNDAMENTALS                 *
                                                                 *
 *****************************************************************/

#include "/ycode/lab.i"


extern NMOL   //number of atoms per molecule
extern KERN   //array setting atom ID's for finding director
extern BEAC   //short for beacon, use to track a substituent.

NMOL = 30;
KERN = [30,20,11,21];
if(is_void(NMS)){
  NMS=array("F",NMOL);
  NMS(20:)="x";
 }
BEAC = [8,9];

/**************************************************************
 ******                                                  ******
 ****           File I/O and limited processing            ****
 ******                                                  ******
 *************************************************************/

func count_files(fn,dir=){
  /*DOCUMENT count_files counts the files*/
  ll=ls(dir,hush=1);
  sift=where(strpart(ll,:strlen(fn)+1)==fn+".");
  ll=ll(sift);
  n_frames=numberof(where(tonum(strpart(ll,-3:))>=0.));
  return n_frames;
}


func read_frame(fn,arcfile=){
  /*DOCUMENT read_frame(file,arcfile=)
    file is .xyz or general traj file to be read
    optional arcfile is a flag to read an archive
  */
  
  extern NATO   //number of atoms
  extern IDS    //atom ID
  extern IDK    //type, I think
  extern NMS    //names of atoms
  extern COORDS //[x,y,z] coordinates of all atoms
  extern BONDS  //list of interatomic bonds
  extern DIMS   //dimensions of bounding box
  extern MASS   //mass holds the mass of each atom type
  extern ENRG   //array of [vdw,coulumb,restraint] energy

    if(is_void(arcfile)){
      file=open(fn);
      file=rdfile(file);
    }
    else file=arcfile;
  NATO=0;
  ev=ec=eg=0.;
  sread,file(1),format="%d %g %g %g",NATO,ev,ec,eg;
  ENRG=[ev,ec,eg];
  DIMS=[0.,0.,0.];
  sread,file(2),format="%g %g %g",DIMS(1),DIMS(2),DIMS(3);
  if(numberof(file)<=NATO) return -789;
  file=file(-NATO+1:);
  IDS=array(6,NATO);
  NMS=array("merde",NATO);
  COORDS=array([9.,9.,9.],NATO);
  IDK=IDS*3/6;
  BONDS=array([0,0,0,0],NATO);
  MASS=IDK*1/3.;
  for(i=1; i<=NATO; i++){
    sread,file(i),format="%d %s %g %g %g %d %d %d %d %d",IDS(i),NMS(i),COORDS(1,i),COORDS(2,i),COORDS(3,i),IDK(i),BONDS(1,i),BONDS(2,i),BONDS(3,i),BONDS(4,i);
  }
  MASS(where(strpart(NMS,1:1)=="C")) = 12.011;
  MASS(where(strpart(NMS,1:1)=="H")) = 1.008;
  MASS(where(strpart(NMS,1:1)=="F")) = 18.998;
return COORDS;
}

func write_xyz(output){
  /*DOCUMENT write_xyz(output_file)
    writes .xyz file based on current external variables
    check COORDS for latest coordinates
  */
  ff=create(output);

  write,ff,format="%d %s\n",NATO,"User Made";
  write,ff,format="%g %g %g %g %g %g\n",DIMS(1),DIMS(2),DIMS(3),90.,90.,90.;
  for(i=1; i<=NATO; i++){
    write,ff,format=" %d %s %g %g %g %d %s\n",IDS(i),NMS(i),COORDS(1,i),COORDS(2,i),COORDS(3,i),IDK(i),sum(totxt(BONDS(where(BONDS(,i)),i))+" "(..));
  }
}

func check_kern(void,kern=){
  /*DOCUMENT check_kern displays current molecular kernel*/
  cm=aa_cm()(,1);
  offset=cm;
  COORDS -= offset;
  window,9;
  fma;
  molec_plot,1;
  if(is_void(kern)) kern = KERN;
  ni=orients(kern=kern);
  print,kern;
  plotl,COORDS(1,kern(1:2)),COORDS(2,kern(1:2)),color="red",width=4;
  plotl,COORDS(1,kern(3:4)),COORDS(2,kern(3:4)),color="red",width=4;
  circ=span(0.,2*pi,25);
  for(i=1; i<=4; i++) plotl,COORDS(1,kern(i))+cos(circ),COORDS(2,kern(i))+sin(circ),color="gray"+["a","b","c","d"](i);
  for(i=1; i<=NMOL; i++) plot,COORDS(1,i),COORDS(2,i),sym=strchar(totxt(i))(-1),color="blue";
  plot,0,0,sym='\5',color="green";
  plot,avg(COORDS(1,kern)),avg(COORDS(2,kern)),color="cyan";
  plotl,[0.,ni(1)],[0.,ni(2)],color="green",width=4;
  plot,0.,0.,color="black";
  COORDS += offset;
  write,format="%s\n%s\n%s\n%s\n","Red lines indicate kernel.","Green vector is director.","Cyan is avg coord.","Black is center of mass.";
}


func recut_key(fnkey,fnout,ncol,new_coords=,unmoor=){
  file = rdfile(fnkey);
  pllr = where( strpart(file, :14)=="RESTRAIN-PLLR ");
  if( numberof(pllr) != ncol) return -999;
  cms = split(aa_cm(),ncol)(:2,avg,);
  new_coords = is_void(new_coords)?cms:new_coords;
  for(i=1; i<=ncol; i++){
    kword = "bird";
    coli = i;
    colx=coly= 0.0;
    kk = 2.85;  // ~ <r^2> = 1 Angstrom^2 at 500K
    dr = 1.0;
    sread,file(pllr(i)),format="%s %d %g %g %g %g\n",kword,coli,colx,coly,kk,dr;
    newx = new_coords(1,i);
    newy = new_coords(2,i);
    kk = is_void(unmoor)?kk:-kk;
    file(pllr(i)) = swrite(format="%s   %2d   %8g   %8g   %3g   %3g","RESTRAIN-PLLR  ",coli,newx,newy,kk,dr);
  }
  new_key = create(fnout);
  for(i=1; i<=numberof(file); i++){
    write,new_key,format="%s\n",file(i);
  }
}
                 

func hex_em(prefix,fnout,qarr,num,rand=){
  /* DOCUMENT builds a grid of molecules
    xarr is an integer array of the number of molecs in each direction
    qarr is a double array of the molecular spacing in each direction
  */
  read_frame,prefix+".xyz";
  file=create(fnout);
  write,file,format=" %d Generated from %s",NATO*7*num,prefix;
  write,file,format="\n %g  %g  %g  90.000 90.000 90.000",3*qarr(1),3*qarr(1),num*qarr(2);
  offset=[0.,0.,0.];
  coords=COORDS;
  for(i=0; i<=6; i++){
    for(j=1; j<=num; j++){
        coords=COORDS+offset;
        th=2.*pi*random();
        if(rand==0) th=rand;
        cm=aa_cm();
        coords-=cm+offset;
        coords=[[cos(th),-sin(th),0.],[sin(th),cos(th),0.],[0.,0.,1.]](+,)*coords(+,);
        if(random()>.5) coords=[[1,0,0],[0,-1,0],[0,0,-1]](+,)*coords(+,);
        coords+=.1*(random(3)>.5)*(random(3)-.5);

        coords+=cm+offset;
        
        for(n=1; n<=NATO; n++){
          write,file,format="\n %-d   %s   %g   %g   %g   %d",IDS(n),NMS(n),coords(1,n),coords(2,n),coords(3,n),IDK(n);
          for(m=1; m<=numberof(where(BONDS(,n))); m++){
            write,file,format="   %d",BONDS(m,n);
            BONDS(m,n)+=NATO;
          }
          IDS(n)+=numberof(NMS);
        }
        offset+=[0.,0.,qarr(2)];
      }
    offset=qarr(1)*[cos(i*pi/3),sin(i*pi/3),0.];
  }
  
  close,file;
  read_frame,fnout;
  atom_plot;
  bond_plot;
}


/********************************************************
 ********************************************************
 ****            Plotting and Visualization          ****
 ********************************************************
 *******************************************************/

func bond_plot(void,color=,width=){
  if(noneof(BONDS)) return -1;
  if(color==[]) color="black";
  width = is_void(width)?2:width;
  for(i=1; i<=dimsof(COORDS)(0); i++){
    bonds=BONDS(where(BONDS(,i)),i);
    dummy=array(i,dimsof(bonds));
    bonds=transpose([dummy,bonds]);
    bonds=bonds(where(bonds));
    plotl,COORDS(1,bonds),COORDS(2,bonds),color=color;
  }
}

func atom_plot(void){
  /*DOCUMENT no longer supported
    load sn.1 for atom_plot
  */
  return "ERR: ATOM PLOT NO LONGER SUPPORTED AS OF sn2.i";
}

func molec_plot(mol_num,color=,width=){
  /*DOCUMENT molec_plot(mol_num,color=,width=)
    plots a single molecule, or an array of molecules
  */
  coords=COORDS;
  bonds=BONDS; 
  for(i=1; i<=numberof(mol_num); i++){
    COORDS = coords(,NMOL*(mol_num(i)-1)+1:NMOL*mol_num(i));
    BONDS = bonds(,NMOL*(mol_num(i)-1)+1:NMOL*mol_num(i));
    BONDS(where(BONDS))-=min(BONDS(where(BONDS)))-1;
    bond_plot,color=color,width=width;
    COORDS = coords;
    BONDS = bonds;
  }
}

func polar_plot_prob(prob,xprb,color=,mag=,offset=){
  // A = int( r(phi) dr dphi
  mag = is_void(mag)?3:mag;
  offset = is_void(offset)?[0.,0.]:offset;
  color = is_void(color)?reddish():color;
  //plot,offset(1),offset(0),color=color;
  dith = random()*pi;
  plotl,offset(1) + mag*cos(xprb + dith)/sqrt(pi),offset(0) + mag*sin(xprb + dith)/sqrt(pi),width=2,color="grayc";
  plotl,offset(1) + mag*sqrt(2*prob)*cos(-xprb),offset(0) + mag*sqrt(2*prob)*sin(-xprb),width=2,color=color;
}



/********************************************************
 ********************************************************
 ****         Fundamental Analysis Tools             ****
 ********************************************************
 *******************************************************/

func aa_cm(void){
  /*DOCUMENT aa_cm() returns center of mass array of all molecules*/

  cm = split(COORDS*MASS(-:1:3,),NATO/NMOL)(,sum,);
  cm /= sum(MASS(:NMOL));
  return cm;
}
func aa_rc(void){
  /*DOCUMENT aa_rc() returns the naive centers of rotations
   defined as the average of the carbons*/
  cat = where(strpart(NMS(:NMOL),1:1)=="C");
  cor = split(COORDS,NATO/NMOL)(,cat,)(,avg,);
  return cor;
}
  
  

func orients(void,kern=){
  /*DOCUMENT orients returns array of molecular orientation vectors
    determined by KERN
    NOT DEFINED AS NORMAL TO MOLECULAR PLANE, but close to it
  */
  n_molecs=numberof(NMS)/NMOL;
  if(kern==[]) kern=KERN;
  v12=COORDS(,kern(1)::NMOL)-COORDS(,kern(2)::NMOL);
  v13=COORDS(,kern(3)::NMOL)-COORDS(,kern(4)::NMOL);
  d12=norm(v12);
  d13=norm(v13);
  orarr=array([0.,0.,0.],n_molecs);
  for(i=1; i<=n_molecs; i++){
    dir=cross(v12(,i),v13(,i));
    dir/=norm(dir);
    orarr(,i)=dir;
  }
  return orarr;
}

func forients(void){
  /*DOCUMENT forient returns array of vectors from central ring,
    through the perfluoro ring, and projected onto a plane
    perpendiculiar to orients() vecs
  */
  coords = split(COORDS,NATO/NMOL);
  center = coords(,KERN,)(,avg,);
  fluoro = coords(,BEAC,)(,avg,);
  fector = fluoro - center;
  
  ni = orients();
  fi = fector - (fector*ni)(sum,)(-:1:3,) * ni; //fi - (fi.ni)ni
  fo = fi/norm(fi)(-:1:3,);
  return fo;
}

func para2(vecarr,ref){
  /*DOCUMENT para2(vecarr,ref)
    returns vector array parallel to ref vector
  */
  ref /= norm(ref);
  if(dimsof(vecarr)(2)!=dimsof(ref)(2)) return "Nein!";
  return (ref*vecarr)(sum,)(-:1:3,)*ref;
}
func perp2(vecarr,ref){
  /*DOCUMENT perp2(vecarr,ref)
    returns vector array perpendiculiar to ref vector
  */
  if(dimsof(vecarr)(2)!=dimsof(ref)(2)) return "Nein!";
  return vecarr - para2(vecarr,ref);
}
  

func compass(void){
  /*DOCUMENT returns 'phi' angle in lab frame
    defined as the angle between the fluorinated ring
    and the X axis
  */
  dir = forients();
  theta = asin(dir(3,));
  plane = dir(:2,)/(cos(theta)(-,));
  phi = acos(plane(1,)) * sign(plane(2,));
  return phi;
}

func dir_angles(void){
  /* DOCUMENT returns theta and phi angles, in the lab frame,
     describing the molecular director orientation.
  */
  ni=orients();
  thi=acos(ni(3,));
  cphi=ni(1,)/sin(thi);
  sphi=ni(2,)/sin(thi);
  phi=acos(cphi)*sign(sphi);
  return transpose([thi,phi]);
}

func neighbors(ref){
  z_cut = 3.;
  r_cut = 16.0;
  com = is_void(com)?aa_cm():com;
  dIR = is_void(dIR)?orients():dIR;
  drr = pbc(com - com(,ref));
  inz = where(norm(para2(drr,dIR(,ref)))<=z_cut);
  inr = where(norm(perp2(drr,dIR(,ref)))<=r_cut);
  uni = union(inz,inr);
  uni = uni(where(uni!=ref));
  return int(uni);
}
  
  

/********************************************************
 ********************************************************
 ****                  Order Parameters              ****
 ********************************************************
 *******************************************************/
  
func oriS(void,vec_out=,ind=){
  //per Soft Matter, 2016, 12, 1295-1312
  ni=orients();
  AA=array(0.,[2,3,3]);
  num=int(numberof(ni)/3);
  for(i=1; i<=num; i++){
    ui=ni(,i);
    AA+=3*ui*ui(-:1:3,)-[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
  }
  AA/=2.*num;
  eval=SVdec(AA,evec);
  if(vec_out==1) return eval;
  return max(eval);
}


/********************************************************
 ********************************************************
 ****             Parsers and Trackers               ****
 ********************************************************
 *******************************************************/


func parse_out(fn,fs,ps,dumps){
  /*DOCUMENT parse_out(fn,fs,ps,dumps) parses fn tinker output
    for energies, temperature, pressure, and density,
    all stored as external arrays
    fs,ps and dumps are the respective time steps and
    the number of of dumps
  */
  pre=strpart(fn,:strfind(".out",fn)(1));
  file=open(fn);
  lil=int(max(1,10*ps/fs)); //outputs per dump cycle w/ 100 steps
  nnn=lil*dumps;
  extern TOTALE;
  extern POTE;
  extern KINE;
  extern GEOE;
  extern VDWE;
  extern TEMP;
  extern PRESS;
  extern DENS;
  TOTALE=POTE=KINE=GEOE=VDWE=TEMP=PRESS=DENS=array(-999.,nnn);

  for(i=1; i<=nnn; i++){
    line = rdline(file);
    while(strpart(line,2:29)!="Average Values for the Last "){
      line = rdline(file);        
    } //positions rdline before data output
    rdline,file,2;
    TOTALE(i) = tonum(strpart(rdline(file),25:37));
    POTE(i) = tonum(strpart(rdline(file),25:37));
    KINE(i) = tonum(strpart(rdline(file),25:37));
    GEOE(i) = tonum(strpart(rdline(file),25:37));
    VDWE(i) = tonum(strpart(rdline(file),25:37));
    TEMP(i) = tonum(strpart(rdline(file),25:37));
    PRESS(i) = tonum(strpart(rdline(file),25:37));
    DENS(i) = tonum(strpart(rdline(file),25:37));
  }
}

func parse_traj(prefix,function,start=,dir=,prefunc=,postfunc=,color=,stop=,pfunc_opt=){
  /* DOCUMENT parse_traj(prefix,scalr_function,start=,dir=)
     prefix should be the trajectory file prefix
     function is a function that takes no arguments and
          returns a scalar. All calculations in function
          should assume the current frame has been read in
          but any non-external variables will need to be calculated
          If a more complex function is desired, consider making
          function interact with it's own external vairables.
     start= is an optional starting frame for use when burning diseq frames
  */
  if(!is_func(function)){
    write,"what's your function?\n";
    return "duh";
  }
  dir = is_void(dir)?"./":dir;
  ogcf = count_files(prefix,dir=dir);
  if(is_void(start)) start = 1;
  if(is_void(stop)) stop = ogcf;
  offset = start - 1;
  extern NCF;
  NCF = 2*ogcf - offset - stop;

  suffix = ".";
  if(start<10) suffix+="0";
  if(start<100) suffix+="0";
  read_frame,dir+prefix+suffix+totxt(start);
  
  dat_arr = array(0.,count_files(prefix,dir=dir)-offset);
  extern COM; //center of mass columnar array
  extern com; //aa_cm() output
  extern DIR; //molecular director columnar array;
  extern dIR; 

  if(is_func(prefunc)) prefunc;
  
  for( i=1+offset; i<=stop; i++){
    fmt="\b\b\b\b";
    suffix=".";
    if(i<1000) fmt="\b\b\b";
    if(i<100){ suffix+="0"; fmt="\b\b"; }
    if(i<10) { suffix+="0"; fmt="\b"; }
    write,format=fmt+"%d",i;
    suffix+=totxt(i);
    read_frame,dir+prefix+suffix;

    
    com = aa_cm();
    COM = split(com,NCOL);
    dIR = orients()
    DIR = split(dIR,NCOL);
    FIR = split(forients(),NCOL);

    srt = sort(COM(3,,)); //z-sort within columns
    COM = COM(,srt);
    DIR = DIR(,srt);
    FIR = FIR(,srt);
     
    frame_dat = function();
    dat_arr(i-offset) = oriS(); //frame_dat(1);

  }
  if(is_func(postfunc)) return postfunc(color=color,offset=pfunc_opt);

  return dat_arr;
}

 
func track_dyn(fn,start=,dir=){
  /*vars to track: dimensions*/
  start=(is_void(start))?1:start;
  dx=DX=dy=DY=dz=DZ=999.0;
  dout=xout=yout=array(0.0,2*count_files(fn,dir=dir)-start+1);
  for(i=1; i<=count_files(fn,dir=dir)-start+1; i++){
    suf=".";
    if(i+start-1<100) suf+="0";
    if(i+start-1<10)  suf+="0";
    suf+=totxt(i+start-1);
    qq=open(fn+suf);
    header=rdline(qq);
    sread,rdline(qq),format="%g %g %g",dx,dy,dz;
    dout(i) = dz; DZ = min(dz,DZ);
    xout(i) = dx; DX = min(dx,DX);
    yout(i) = dy; DY = min(dy,DY);
    if(i>1){
      plotl,[i-1,i],dout([i-1,i])/DZ,color=reddish();
      plotl,[i-1,i],xout([i-1,i])/DX,color=blueish();
      plotl,[i-1,i],yout([i-1,i])/DY,color=blueish();
    }
    pause,1;
   }
  dout=dout(:i-1);
  yout=yout(:i-1);
  xout=xout(:i-1);
  plotl,,dout,color="red";
  plotl,,yout,color="blue";
  plotl,,xout,color="green";
  for(j=1; j<=i/100; j++){
    plotl,j*100*[1,1],dout(j*100)*[1,1]+[1,-1],color="grayd";
  }
  unzoom;
  write,format="\nMin x: %g  Min y: %g  Min z:  %g\n",min(xout),min(yout),min(dout);
  return transpose([xout,yout,dout]);
}   


/********************************************************
 ********************************************************
 ****                   Utilities                    ****
 ********************************************************
 *******************************************************/

func pbc(scriptr,wrap=){
  /*DOCUMENT pbc(dr,wrap=) applies pbc's to the array r_i - r_j
    returns shifted array
  */
  for(i=1; i<=dimsof(DIMS)(2); i++){
    fi = anyof(abs(scriptr(i,))>DIMS(i)/2.);
    if(fi){
      wi = where(abs(scriptr(i,))>DIMS(i)/2.);
      scriptr(i,wi) -= DIMS(i)*sign(scriptr(i,wi));
    }
  }
  
  return scriptr;
}


func glob_spin(th,raxis=){
  /*DOCUMENT glob_spin(th,raxis=) rotates COORDS by th about
    optional raxis= rotation axis.
    raxis = 1,2,3 for x y and z raxes respectively
    raxis defaults to z axis
  */    
  raxis = is_void(raxis)?3:raxis; //assumes rotion around z-axis
  if(raxis==1){
    rotm = [[1,0,0],[0,cos(th),sin(th)],[0,-sin(th),cos(th)] ];
  }
  else if(raxis==2){
    rotm = [[cos(th),0,-sin(th)],[0,1,0],[sin(th),0,cos(th)] ];
  }
  else if(raxis==3){
    rotm = [[cos(th),sin(th),0],[-sin(th),cos(th),0],[0,0,1] ];
  }
  else return "-999";
  COORDS=rotm(+,)*COORDS(+,);
}

func utili_plot(void){
  temp = [300,350,400,425,450,475,500,525,550];
  datarr = "K"+totxt(temp)+"_intraph.dat";

  window,9,dpi=100,width=1000,height=2000,style="nobox.gs";
  plt,"F5",-.5,-2,tosys=1;
  plt,"F8",1.5,-2,tosys=1;
  plt,"F6",3.5,-2,tosys=1;
  plt,"F7",5.5,-2,tosys=1;
  pltitle,"Intracolumnar nearest neighbor orientation";
  limits,-5,8,-3,25;
  
  for(i=1; i<=numberof(temp); i++){
    plt,"Temp = "+totxt(temp(i))+"K",-4,3*(i-1),tosys=1;
    dat = openb(datarr(i));
    restore,dat;
    
    prob=F5_fdph/sum(F5_fdph*F5_xdph(dif)(1));
    prob(0)=prob(1);
    polar_plot_prob,prob,F5_xdph,mag=1,offset=[0,3*(i-1)],color="blue";
    
    prob=F8_fdph/sum(F8_fdph*F8_xdph(dif)(1));
    prob(0)=prob(1);
    polar_plot_prob,prob,F8_xdph,mag=1,offset=[2,3*(i-1)],color="red";
    
    prob=F6_fdph/sum(F6_fdph*F6_xdph(dif)(1));
    prob(0)=prob(1);
    polar_plot_prob,prob,F6_xdph,mag=1,offset=[4,3*(i-1)],color="blue";
    
    prob=F7_fdph/sum(F7_fdph*F7_xdph(dif)(1));
    prob(0)=prob(1);
    polar_plot_prob,prob,F7_xdph,mag=1,offset=[6,3*(i-1)],color="red";

  }
  close,dat;
  return 0;
}

func ks_plots(prob1,prob2,x,y){
  //prob(x,y)
  //both probs must be on same x y mesh!
  //first, KS y slices wrt x
  ksx = prob1(,1)*0.;
  for(nx=1; nx<=numberof(x); nx++){
    bk = sum(prob1(nx,));
    cd1 = prob1(nx,)(cum);
    if(bk) cd1/=bk;

    bk = sum(prob2(nx,));
    cd2 = prob2(nx,)(cum);
    if(bk) cd2/=bk;

    ksx(nx) = max(abs(cd1 - cd2));
  }
  plotl,x,ksx,color=blueish(),width=2;


  
  ksy = prob1(1,)*0.;
  for(ny=1; ny<=numberof(y); ny++){
    bk = sum(prob1(,ny));
    cd1 = prob1(,ny)(cum);
    if(bk) cd1/=bk;

    bk = sum(prob2(,ny));
    cd2 = prob2(,ny)(cum);
    if(bk) cd2/=bk;

    ksy(ny) = max(abs(cd1 - cd2));
  }
  plotl,y,ksy,color=reddish(),width=2;

}


func avg_sd(dat_arr){
  mean = avg(dat_arr);
  sd = sqrt(avg((dat_arr - mean)^2));
  if(anyof(abs(dat_arr - mean)/sd > 5. )){
    //write,"dropping outliers > 5sd";
    w = where(abs(dat_arr-mean)/sd <= 5.);
    new_arr = dat_arr(w);
    return avg_sd(new_arr);
  }
  return [mean,sd];
}
