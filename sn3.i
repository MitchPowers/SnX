/*****************************************************************
 *
 *  MITCH'S TINKER ANALYSIS SUITE - FUNDAMENTALS
    houses basic functions used by all other tin*.i files        *
                                                                 *
 *****************************************************************/

#include "/ycode/lab.i"
//be aware, Mitch's lab.i has forked

extern NMOL   //number of atoms per molecule
extern KERN   //array setting atom ID's for finding director
extern BEAC   //short for beacon, use to track a fluorine ring.
extern MPC    //molecules per column
extern NCOL   //number of columns

NMOL = 30;
KERN = [30,20,11,21];
BEAC = [8,9]; //F2 and F3
MPC = 20;
NCOL = 7;

/**************************************************************
 ******                                                  ******
 ****           File I/O and limited processing            ****
 ******                                                  ******
 *************************************************************/

func count_files(fn,dir=){
  /*DOCUMENT count_files counts the files*/
  if(numberof(dir)>1){
    out = array(0,numberof(dir));
    for(d=1; d<=numberof(dir); d++){
      out(d) = count_files(fn,dir=dir(d));
    }
    return out;
  }
  ll=ls(dir,hush=1);
  sift=where(strpart(ll,:strlen(fn)+1)==fn+".");
  if(numberof(sift)==0) return 0;
  ll=ll(sift);
  n_frames=numberof(where(tonum(strpart(ll,-3:))>=0.));
  //Only counts files with a number as extension
  return n_frames;
}


func read_frame(fn){
  /*DOCUMENT read_frame(file)
    file is .xyz or general traj file to be read
  */
 
  extern NATO   //number of atoms total
  extern MOLS  
  extern IDS    //atom ID
  extern IDK    //type, I think
  extern NMS    //names of atoms
  extern COORDS //[x,y,z] coordinates of all atoms
  extern BONDS  //list of interatomic bonds
  extern DIMS   //dimensions of bounding box
  extern MASS   //mass holds the mass of each atom type
  extern ENRG   //array of [vdw,coulumb,restraint] energy  

  file=open(fn);
  file=rdfile(file);

  NATO=0;
  ev=ec=eg=0.;
  sread,file(1),format="%d %g %g %g",NATO,ev,ec,eg;
  MOLS = NATO/NMOL;
  ENRG=[ev,ec,eg];
  DIMS=[0.,0.,0.];
  sread,file(2),format="%g %g %g",DIMS(1),DIMS(2),DIMS(3);
  if(numberof(file)<=NATO) return -789;
  file=file(-NATO+1:); //burns header
  IDS=array(6,NATO);
  NMS=array("merde",NATO);
  COORDS=array([9.,9.,9.],NATO);
  IDK=IDS*3/6;
  BONDS=array([0,0,0,0],NATO);
  MASS=IDK*1/3.;
  for(i=1; i<=NATO; i++){
    sread,file(i),format="%d %s %g %g %g %d %d %d %d %d",IDS(i),NMS(i),COORDS(1,i),COORDS(2,i),COORDS(3,i),IDK(i),BONDS(1,i),BONDS(2,i),BONDS(3,i),BONDS(4,i);
  }
  MASS(where(strpart(NMS,1:1)=="C")) = 12.011; //<-Tinker | NIST-> 12.0106;
  MASS(where(strpart(NMS,1:1)=="H")) = 1.008; //<-Tinker | NIST -> 1.00798;
  MASS(where(strpart(NMS,1:1)=="F")) = 18.998; //<-Tinker | NIST ->18.9984;
  /* Tinker Masses from poledit.f
    NIST mass values taken from NIST table
     'Atomic Weights and Isotopic Compositions for All Elements'
     additional masses should be added as needed
  */
  NCOL = NATO/NMOL/MPC;
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
  /*DOCUMENT check_kern displays current molecular kernel
    the kernel is used to define the molecular normal by way of
    a cross product. This approach relies on the central ring
    being flat, and representative of the molecular plane.
    This should be practically equivalent between species.
    
    use optional kern= array to test new kernels
  */
  cm=aa_cm()(,1);
  offset=cm;
  COORDS -= offset;

  molec_plot,1,width=4;
  if(is_void(kern)) kern = KERN;
  ni=orients(kern=kern);
  print,kern;
  //plotl,COORDS(1,kern(1:2)),COORDS(2,kern(1:2)),color="red",width=4;
  //plotl,COORDS(1,kern(3:4)),COORDS(2,kern(3:4)),color="red",width=4;
  circ=span(0.,2*pi,25);
  for(i=1; i<=NMOL; i++) plt,totxt(i),COORDS(1,i),COORDS(2,i),tosys=1,justify="CC",color="blue";
  wf=where(NMS=="F");
  plotlc,COORDS(1,wf),COORDS(2,wf),2.85/2,color="blue";
  wh=where(NMS=="HA");
  plotlc,COORDS(1,wh),COORDS(2,wh),2.42/2,color="grayb";
  wc=where(NMS=="CA");
  plotlc,COORDS(1,wc),COORDS(2,wc),3.55/2,color="graya";
  plot,0,0,sym='\5',color="green";
  plot,COORDS(1,avg),COORDS(2,avg),color="cyan";
  plotl,[0.,ni(1)],[0.,ni(2)],color="green",width=4;
  COORDS += offset;
  write,format="%s\n%s\n%s\n%s\n","Red lines indicate kernel.","Green vector is director.","Cyan is avg coord.";
}
func molec_illus(void,shift=){
  /*DOCUMENT molec_illus orients a molecule and plots it in detail

    shift is the offset applied to the centroid
   */
  centroid = COORDS(,:NMOL)(,avg);
  shift = is_void(shift)?0.(-:1:3):shift;
  while( numberof(shift)<3 ) grow,shift,0.;
  offset = centroid;
  COORDS -= offset;
  gmole_spin,1;
  glob_spin,pi/2,raxis=1;
  glob_spin,pi/2,raxis=2;
  COORDS -= shift;
  offset += shift;

  cm=aa_cm()(,1);
  centroid = COORDS(,:NMOL)(,avg);

  molec_plot,1,width=4;
  bow = (COORDS(,[2, 3])(,avg) - centroid)*1.3 + centroid;
  port= (COORDS(,[5,16])(,avg) - centroid)*1.3 + centroid;
  sbrd= (COORDS(,[6,26])(,avg) - centroid)*1.3 + centroid;
  plotl,[centroid(1),bow(1)],[centroid(2),0*bow(2)],color="grayd",width=4;
  plotl,[centroid(1),port(1)],[centroid(2),0*port(2)],color="red",width=4;
  plotl,[centroid(1),sbrd(1)],[centroid(2),0*sbrd(2)],color="green",width=4;
  
  ni=orients();
  /*
  for(i=1; i<=4; i++){
    plt,["A","C","B","D"](i),COORDS(1,KERN(i)),COORDS(2,KERN(i))+[-.1,.4,.4,-.1](i),tosys=1,justify="CC",color="blue";
    }*/
  wf=where(NMS=="F");
  plotlc,COORDS(1,wf),COORDS(2,wf),2.85/2,color="blue";
  wh=where(NMS=="HA");
  plotlc,COORDS(1,wh),COORDS(2,wh),2.42/2,color="grayb";
  wc=where(NMS=="CA");
  plotlc,COORDS(1,wc),COORDS(2,wc),3.55/2,color="graya";
  centroid = COORDS(,avg);
  plot,centroid(1),centroid(2),color="grayd";
  plot,cm(1),cm(2),sym='\5',color=[240,0,240];
  plotl,cm(1)+[0.,ni(1)],cm(2)+[0.,ni(2)],color=[240,0,240],width=4;

  COORDS += offset;
  write,format="%s\n%s\n%s\n%s\n","Red lines indicate kernel.","Green vector is director.","Cyan is avg coord.";
}



func recut_key(fnkey,fnout,ncol,new_coords=,unmoor=,quiet=){
  /*DOCUMENT recut_key updates a .key file with PLLR coords
    YOU MUST READ IN THE CURRENT FRAME BEFORE USING
    
    recut_key(fnkey,fnout,ncol,new_coords=,unmoor=)
    first two args are source and dest. If they are the same
    the original key will get clobbered without warning

    ncol is the number of columns, and therefore the number
    of constraints. While this can be determined internally,
    keeping it in the input args helps prevent accidental use
    of this function.

    new_coords= is an array of [x,y] coords for the constraints
    unmoor flag toggles the potential's unmoored flag. It is only
    needed when changing the current state.
   */
  file = rdfile(fnkey);
  pllr = where( strpart(file, :14)=="RESTRAIN-PLLR ");
  if( numberof(pllr) != ncol) return -999;
  cms = split(aa_cm(),ncol)(:2,avg,);
  new_coords = is_void(new_coords)?cms:new_coords;
  for(i=1; i<=ncol; i++){
    kword = "bird";
    coli = i;
    colx = coly= 0.0;
    kk = 20.0;
    dr = 1.5;
    sread,file(pllr(i)),format="%s %d %g %g %g %g\n",kword,coli,colx,coly,kk,dr;
    /*
    dr = 0.0;
    kk = 1.0;
    */
    newx = new_coords(1,i);
    newy = new_coords(2,i);
    kk = is_void(unmoor)?kk:-kk; //note, this can 'remoor' as well
    file(pllr(i)) = swrite(format="%s   %2d   %6g   %6g   %4g   %4g","RESTRAIN-PLLR  ",coli,newx,newy,kk,dr);
  }
  new_key = create(fnout);
  for(i=1; i<=numberof(file); i++){
    write,new_key,format="%s\n",file(i);
    if(is_void(quiet)) write,format="%s\n",file(i);
  }
}
     
func hex_em(void){
  /*DOCUMENT hex_em is supported in sn2.i*/
  write,format="%s\n","hex_em is not currently supported beyon sn2.i";
}


/********************************************************
 ********************************************************
 ****            Plotting and Visualization          ****
 ********************************************************
 *******************************************************/

func bond_plot(void,color=,width=){
  /*DOCUMENT bond_plot plots molecules, showing connectivity*/
  if(noneof(BONDS)) return -1;
  if(color==[]) color="black";
  width = is_void(width)?2:width;
  for(i=1; i<=dimsof(COORDS)(0); i++){
    bonds=BONDS(where(BONDS(,i)),i);
    dummy=array(i,dimsof(bonds));
    bonds=transpose([dummy,bonds]);
    bonds=bonds(where(bonds));
    plotl,COORDS(1,bonds),COORDS(2,bonds),color=color,width=width;
  }
}
func molec_plot(mol_num,color=,width=){
  /*DOCUMENT molec_plot(mol_num,color=,width=)
    plots a single molecule, or an array of molecules
  */
  coords=COORDS;
  bonds=BONDS;
  if(catch(99)){
    COORDS = coords;
    BONDS = bonds;
    write,format="ERR: %s\n","Something amiss in molec_plot";
  }
  for(i=1; i<=numberof(mol_num); i++){
    COORDS = coords(,NMOL*(mol_num(i)-1)+1:NMOL*mol_num(i));
    center = COORDS(,int(random()*NMOL)+1);
    disp = abs(COORDS-center)
    for(xi=1; xi<=3; xi++){ //don't split mols
      if(anyof(disp(xi,)>DIMS(xi)*.5)){
        ww = where( disp(xi,)>DIMS(xi)*.5);
        COORDS(xi,ww) += DIMS(xi)*sign( COORDS(xi,ww)-center(xi) );
      }
    }
    BONDS = bonds(,NMOL*(mol_num(i)-1)+1:NMOL*mol_num(i));
    BONDS(where(BONDS))-=min(BONDS(where(BONDS)))-1;
    bond_plot,color=color,width=width;
    COORDS = coords;
    BONDS = bonds;
  }
}
func plot_molec(mol_num,color=,width=){
  write,format="%s\n","it's called 'molec_plot' you putz";
  molec_plot,mol_num,color=color,width=width;
}




/********************************************************
 ********************************************************
 ****                                                ****
 ****         Fundamental Analysis Tools             ****
 ****    ~~ These are worth rewriting in C ~~        ****
 ****                                                ****
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
    defined as normal to the molecular plane, determined by orthogonal regression
    determined by KERN
    Follows Eberly "Least Squares Fitting of Data by Linear or Quadratic Structures"
    
  */
  n_molecs=numberof(NMS)/NMOL;
  if(kern==[]) kern=KERN;
  v12=COORDS(,kern(1)::NMOL)-COORDS(,kern(2)::NMOL);
  v13=COORDS(,kern(3)::NMOL)-COORDS(,kern(4)::NMOL);
  d12=norm(v12);
  d13=norm(v13);
  orarr=array([0.,0.,0.],n_molecs);
  for(i=1; i<=n_molecs; i++){
    inds = indgen(1+NMOL*(i-1):NMOL*i);
    og = COORDS(,inds)(,avg);
    cog = COORDS(,inds) - og;
    c2c = cog(,+)*cog(,+);
    svd = SVdec(c2c,evecs);
    dir = evecs(,where(svd==min(svd)));
    //dir is via orthogonal regression
    orarr(,i)=dir;

    vdir=cross(v12(,i),v13(,i));
    vdir/=norm(dir);
    //vdir is via fixed cross product
    
    orarr(,i)*=sign( vdir(+)*dir(+) );
    //sign enforces up/down convention
  }
  if(anyof(abs(orarr)<1e-9)) orarr(where(abs(orarr)<1e-9)) = 0.;
  return orarr;
}
func forients(void,ni=){
  /*DOCUMENT forient returns array of vectors from central ring,
    through the perfluoro ring, and projected onto a plane
    perpendiculiar to orients() vecs
    Use ni= orients vectors to avoid recalculating them.
  */
  coords = split(COORDS,NATO/NMOL);
  center = coords(,avg,); //centroid
  fluoro = coords(,BEAC,)(,avg,); //between F2 and F3
  fector = fluoro - center;

  ni = is_void(ni)?orients():ni;
  fi = fector - (fector*ni)(sum,)(-:1:3,) * ni; //fi - (fi.ni)ni
  fo = fi/norm(fi)(-:1:3,);
  return fo;
}

func compass(void){
  /*DOCUMENT returns 'phi' angle in lab frame
    defined as the angle between the fluorinated ring
    and the +X axis
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
  flag=0;
  ni=orients();
  if(anyof(ni(3,)==1)){          //pointed straight up makes an
    ww = where( ni(3,)==1 );     //azimuthal bearing ill defined
    ni(,ww) = [1.,1.,1.]/sqrt(3);//and risks div by zero errors.
    flag=1; 
  }
  thi=acos(ni(3,));
  cphi=ni(1,)/sin(thi);
  sphi=ni(2,)/sin(thi);
  phi=acos(cphi)*sign(sphi);
  if(flag){
    thi(ww) = 0.;
    phi(ww) = 0.;
  }
  return transpose([thi,phi]);
}
func tout_langles(void){
  /*DOCUMENT returns array of angles defining the orientations
    specifically, [theta, phi, gamma] in the lab frame, as defined:
    theta = angle between z-axis and molec director
    phi = angle between x-axis and molec director
    gamma = angle between fluoraxis and x-axis
    This function is redundant, but convenient
  */
  ni = orients();
  thi = acos(ni(3,));
  thic= thi + !thi*pi/2; //thic is 'cleaned' to avoid /0 errors
  cphi = ni(1,)/sin(thic);
  sphi = ni(2,)/sin(thic);
  if(anyof(abs(cphi)>1.0000) ){
    write,format="Err in TlA cphi: %g\n",cphi(where(abs(cphi)>1.0000));
    cphi = min(cphi, 1.0000);
    cphi = max(cphi,-1.0000);
  }
  phi=acos(cphi);
  phi*=sign(sphi);

  dir = forients(ni=ni);
  theta = asin(dir(3,));
  plane = dir(:2,)/(cos(theta)(-,));
  if(anyof(abs(plane)>1.0000) ){
    plane = min(plane, 1.0000);
    plane = max(plane,-1.0000);
  }
  fhi = acos(plane(1,)) * sign(plane(2,));

  return transpose([thi,phi,fhi]);
}



/////
// //THESE ARE A HIGH PRIORITY FOR COMPILING
/////
func para2(vecarr,ref){
  /*DOCUMENT para2(vecarr,ref)
    returns vector array parallel to ref vector
  */
  if(dimsof(vecarr)(2)!=dimsof(ref)(2)) return "Nein!";
  ref /= norm(ref)(-,);
  return (ref*vecarr)(sum,)(-:1:3,)*ref;
}
func perp2(vecarr,ref){
  /*DOCUMENT perp2(vecarr,ref)
    returns vector array perpendiculiar to ref vector
  */
  if(dimsof(vecarr)(2)!=dimsof(ref)(2)) return "Nein!";
  ref /= norm(ref);
  return vecarr - (ref*vecarr)(sum,)(-:1:3,)*ref;
}
/////
// //
/////



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
  if(vec_out==1) return evec;
  return max(eval);
}


/********************************************************
 ********************************************************
 ****                   Utilities                    ****
 ********************************************************
 *******************************************************/

func pbc(scriptr){
  /*DOCUMENT pbc(dr) applies pbc's to the array r_i - r_j
    returns shifted array
  */
  if( dimsof(scriptr)(2) > numberof(DIMS) ) return "Too Few DIMS";
  
  for(i=1; i<=dimsof(scriptr)(2); i++){
    fi = anyof(abs(scriptr(i,))>DIMS(i)/2.);
    if(fi){
      wi = where(abs(scriptr(i,))>DIMS(i)/2.);
      scriptr(i,wi) -= DIMS(i)*sign(scriptr(i,wi));
    }
  }
  return scriptr;
}


func glob_spin(th,raxis=,matout=){
  /*DOCUMENT glob_spin(th,raxis=) rotates COORDS
    by th about optional raxis= rotation axis.
    raxis = 1,2,3 for x y and z raxes respectively
    raxis defaults to z axis
    matout flag returns the rotation matrix
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
  if(!is_void(matout)) return rotm;
  COORDS=rotm(+,)*COORDS(+,);
}
func gmole_spin(molec){
  /*DOCUMENT gmole_spin(molec)
    rotates molec so n==z and f==x
    rotations are global
  */
  glob_spin,tout_langles()(2,molec);
  glob_spin,tout_langles()(1,molec),raxis=2;
  glob_spin,tout_langles()(3,molec);
}

func distro_dat(f_arr,x_arr){
  /*DOCUMENT distro_dat(f,x) takes binned freqs and returns
    expectation value, variance, and FWHM
    ASSUMES SINGLE PEAK DISTROS
    ASSUMES idiosynchratic sets
    ASSUMES constant bin widths
  */
  dx = x_arr(dif)(1);
  pdens = f_arr(:-1)/sum(f_arr*dx);
  zcens = x_arr(zcen)-dx;
  
  expectation = sum( pdens * zcens * dx);
  sqrt_sigma2 = sqrt( sum( pdens * zcens^2 * dx) - expectation^2 );

  HM = max(pdens)/2.;
  FW = zcens(where(pdens>HM));
  FWHM = FW(0)-FW(1);

  write,format="< X > = %g\n sig = %g\n FWHM = %g\n",expectation,sqrt_sigma2,FWHM;
  return [expectation,sqrt_sigma2,FWHM];
}


func gauss(x,params){
  A_o = params(1,);
  sig = params(2,);
  x_o = params(3,);

  yar = -0.5*(x-x_o)^2/sig^2;
  yar = A_o*exp(yar)/(sig*sqrt(2*pi));
  return yar;
}
func gfit(dat,xarr,params){
  bg  = min(dat);
  dat-= bg;
  nlm = levmar(dat,xarr,gauss,params,AVAR,ACOV,wgt=dat);
  AVAR*=levmar_chi2;
  AVAR = sqrt(AVAR);
  write,format="A_o = %g (%g)\nsig = %g (%g)\nx_o = %g (%g)\n",nlm(1),AVAR(1),nlm(2),AVAR(2),nlm(3),AVAR(3);
  return gauss(xarr,nlm)+bg;  
}
/*****************************************************************
 *                                                         *******
 **                                                         ******
 ***                                                         *****
 ****        TRAJECTORY PARSER AND RELATED FUNCTIONS          ****
 *****                                                         ***
 ******                                                         **
 *******                                                         *
 *****************************************************************/

func parse_traj(prefix,start=,stop=,dir=,prefunct=,funct=,postfunct=,color=){
  /* DOCUMENT parse_traj(prefix,start=,stop=,dir=,prefunct=,funct=,postfunct=)
     prefix should be the trajectory file prefix
     start= is an optional starting frame for use when burning diseq frames
     stop= last frame to be read
     dir= is directory to use, if not the cwd

     prefunct is a function called before parsing the files
     funct is called with each frame
     postfunct is called after reading all frames

     the functs are all optional and take no arguments.
     prefunct's return is ignored
     funct may return a scalar
     postfunct can return anything

     These functs should interact with the data via externals,
      a few useful variables are calculated regardless of funcs
      
  */
  dir = is_void(dir)?"./":dir;
  ogcf = count_files(prefix,dir=dir);
  if(is_void(start)) start = 1;
  if(is_void(stop)) stop = ogcf;
  offset = start - 1;
  extern NCF;
  NCF = 2*ogcf - offset - stop; //Total number of files to be read

  suffix = ".";
  if(start<10) suffix+="0";
  if(start<100) suffix+="0";
  read_frame,dir+prefix+suffix+totxt(start);

  NCOL = is_void(NCOL)?NCOL:NATO/NMOL/MPC;
  
  dat_arr = array(0.,NCF);
  extern CAT; //Coords at all time
  CAT = array(0.,[3,10,NATO/NMOL,NCF]);
  extern CoM; //aa_cm() output -> CAT(1:3,,)
  extern TlA; //tout_langles() output -> CAT(4:6,,)
  extern CoR; //rotation center -> CAT(7:9,,);
              //CAT(10,,) is for flags

  extern DIM_ARR;
  DIM_ARR = array([0.,0.,0.],NCF);

  frame_dat = array(0.,NCF);

  if(is_func(prefunct)) prefunct;
 
  for( i=1+offset; i<=stop; i++){
    fmt="\b\b\b\b\b";
    suffix=".";
    if(i<1000) fmt="\b\b\b";
    if(i<100){ suffix+="0"; fmt="\b\b"; }
    if(i<10) { suffix+="0"; fmt="\b"; }
    write,format=fmt+"%d",i;
    suffix+=totxt(i);
  
    read_frame,dir+prefix+suffix;
    
    CoM = aa_cm();
    TlA = tout_langles();
    CoR = split(COORDS,NCOL*MPC)(,KERN,)(,avg,);
    CAT(1:3,,i-offset) = CoM;
    CAT(4:6,,i-offset) = TlA;
    CAT(7:9,,i-offset) = CoR;
    CAT( 10,,i-offset) = -999;

    DIM_ARR(,i-offset) = DIMS;
    
    frame_dat(i-offset) = funct();

    //delete this
    //plot,i,CoM(3,224);
    //pause,0;
    
    
  }
  if(is_func(postfunct)) postfunct;
  
  write,format="Externals:\n %s\n","DIM_ARR","CAT";
  
  return frame_dat;
}


func parse_out(fn,concat=,quiet=){
  /*DOCUMENT parse_out(fn) parses fn tinker output
    for energies, temperature, pressure, and density,
    all stored as external arrays

    concat will concatenate values on current externs
  */
  quiet = is_void(quiet)?0:quiet;
  if( (numberof(fn)==allof(dimsof(fn)==1) ) ){
    fn=fn(1);
    concat=[];
  }
  if(numberof(fn)>1){
    lens = array(0,numberof(fn));
    for(n=1; n<=numberof(fn); n++){
      if(n==1) len = parse_out(fn(n),quiet=quiet);
      else len = parse_out(fn(n),concat=1,quiet=quiet);
      if(!quiet) write,len;
      lens(n) = len;
    }
    return lens;
  }else{
    
    pre=strpart(fn,:strfind(".out",fn)(1));
    file=open(fn);
    if(!quiet) write,format="reading %s...\n",fn;
    extern TOTALE;
    extern POTE;
    extern KINE;
    extern GEOE;
    extern VDWE;
    extern TEMP;
    extern PRESS;
    extern DENS;
    extern ODIMS;
    if(is_void(concat)){
      TOTALE=POTE=KINE=GEOE=VDWE=TEMP=PRESS=DENS=ODIMS=array(-999.,1000);
      ODIMS = array(0.,[2,3,500]);
      i0=1;
    } else {
      i0=numberof(VDWE)+1;
    }
    extern line;
    line = "lorem ipsum";
    for(i=i0; line!=string(0); i++){
      line = rdline(file);
      while((line!=string(0))*(strpart(line,2:29)!="Average Values for the Last ")){
        if(strpart(line,2:17)=="Lattice Lengths "){
          if(ODIMS(sum,0)!=0.) grow,ODIMS,ODIMS*0.;
          mk = where(ODIMS(sum,)==0)(1);
          ODIMS(1,mk) = tonum(strpart(line,25:37));
          ODIMS(2,mk) = tonum(strpart(line,40:51));
          ODIMS(3,mk) = tonum(strpart(line,54:65));
        }
        line = rdline(file);
      } //positions rdline before data output
      rdline,file,2;
      if(i>=numberof(TOTALE)-1){
        if(!quiet) write,"embiggening externs";
        grow,TOTALE,TOTALE*0.;
        grow,POTE,POTE*0.;
        grow,KINE,KINE*0.;
        grow,GEOE,GEOE*0.;
        grow,VDWE,VDWE*0.;
        grow,TEMP,TEMP*0.;
        grow,PRESS,PRESS*0.;
        grow,DENS,DENS*0.;
      }
      if(line!=string(0)){
        TOTALE(i) = tonum(strpart(rdline(file),25:37));
        POTE(i) = tonum(strpart(rdline(file),25:37));
        KINE(i) = tonum(strpart(rdline(file),25:37));
        GEOE(i) = tonum(strpart(rdline(file),25:37));
        VDWE(i) = tonum(strpart(rdline(file),25:37));
        TEMP(i) = tonum(strpart(rdline(file),25:37));
        PRESS(i) = tonum(strpart(rdline(file),25:37));
        DENS(i) = tonum(strpart(rdline(file),25:37));
      } else i--;
    }
  if(DENS(i)<=0) i-=1;
  if(i==i0){
    TOTALE=POTE=KINE=GEOE=VDWE=TEMP=PRESS=DENS=ODIMS=[];
    if(!quiet) write,format="%s\n","Warning: Empty File";
    return 0;
  }
  if(!quiet) write,format="%s\n","Externs:";
  if(!quiet) write,"TOTALE","POTE","KINE","GEOE","VDWE","TEMP","PRESS","DENS","ODIMS";
  i = i-i%2; //Take out up to last picosecond to match dyn and data
  TOTALE = TOTALE(:i);
  POTE = POTE(:i);
  KINE = KINE(:i);
  GEOE = GEOE(:i);
  VDWE = VDWE(:i);
  TEMP = TEMP(:i);
  PRESS = PRESS(:i);
  DENS = DENS(:i);
  ODIMS = ODIMS(,:mk);
  }
  return numberof(TOTALE)/2;
}

func parse_dims(prefix,dir=){
  dir = is_void(dir)?"./":dir;
  ogcf = count_files(prefix,dir=dir);
  extern NCF;
  NCF = ogcf

  extern DIM_ARR;
  extern DIMS;
  DIM_ARR = array([0.,0.,0.],NCF);
  DIMS = [0.,0.,0.];
  
  for( i=1; i<=NCF; i++){
    write,format="\r%g",i*100./NCF;
    suffix=".";
    if(i<100) suffix+="0";
    if(i<10)  suffix+="0";
    suffix+=totxt(i);
    frame = open(dir+prefix+suffix);
    rdline,frame;
    sread,rdline(frame),format="%g %g %g",DIMS(1),DIMS(2),DIMS(3);
    DIM_ARR(,i) = DIMS;
  }
  write,format="\rExternals: %s\n","DIM_ARR"; 
  return min(DIM_ARR);
}

   
func energy_width(void,dir=,spec=){
  //winkill,7;
  //window,7,style="boxed.gs",wait=1;
  //fma;
  
  wind = 500; //window width in ps
  wind = int(wind/2.);
  tmax = numberof(POTE);

  navg = 0.;
  nsig = 0.;

  numm = NCOL*MPC;

  avgarr = array(0.0,[2,tmax,3]);
  fwdarr = array(0.0,[2,tmax,3]);
  avgarr(,3)=fwdarr(,3)=0.5*indgen(tmax);

  minmax = [999999.,-999999.];
  kb = 0.831446262/418.40; //kcal/mol/kelvin
  for(t=wind; t<=tmax-wind; t++){
    awin = t-max(1,t-wind+1); omeg = min(tmax,t+wind)-t;
    twin = min(awin,omeg); //makes window symmetric
    ewin = TOTALE( max(1,t-wind+1):min(tmax,t+wind) );
    /*
      pwin = POTE( t-twin:t+twin );
    kwin = KINE( t-twin:t+twin );
    ratw = pwin/kwin;
    navg = ratw(avg);
    nsig = sqrt( avg( (ratw - navg)^2 ) );
    avgarr(t,) = [navg,nsig,t*0.5];

    rfwd = (POTE(t:)/KINE(t:));
    nfwd = rfwd(avg);
    nwsg = sqrt( avg( (rfwd - nfwd)^2) );
    fwdarr(t,) = [nfwd,nwsg,t*0.5];
    */
    
    ewin = TOTALE( t-twin:t+twin )/NMOL;
    navg = ewin(avg);
    nsig = sqrt( avg( (ewin - navg)^2 ) );

    //kb_T = KINE( t-twin:t+twin )(avg)*(2./6); //equipartition ftw
    kb_T = kb*TEMP(t-twin:t+twin )(avg);
    nsig/=kb_T;
    navg/=kb_T;
    
    avgarr(t,) = [navg,nsig,t*0.5];
    
    nfwd = TOTALE(t:)(avg) / NMOL;
    nwsg = sqrt( avg( (TOTALE(t:)/NMOL - nfwd)^2) )/kb_T;
    nfwd/= kb_T;
    fwdarr(t,) = [nfwd,nwsg,t*0.5];
    

    //plotl,t*.5*[1,1]*.001,nfwd+nwsg*[1,-1],color="graya";
    plotl,t*.5*[1,1]*.001,navg+nsig*[1,-1],color="grayc";
    if( nfwd+nwsg>minmax(0) ) minmax(0) = nfwd+nwsg;
    if( navg+nsig>minmax(0) ) minmax(0) = navg+nsig;
    if( nfwd-nwsg<minmax(1) ) minmax(1) = nfwd-nwsg;
    if( navg-nsig<minmax(1) ) minmax(1) = navg-nsig;

    /*
      windo=window();
    window,16,wait=1;
    plot,t*.5,kb_T/TEMP( t-twin:t+twin )(avg),color="red",sym='\1';
    window,windo,wait=1;
    */
  }

  avgarr(:wind-1,:2) = avgarr(wind,:2)(-:1:wind-1,);
  fwdarr(:wind-1,:2) = fwdarr(wind,:2)(-:1:wind-1,);

  avgarr(tmax-wind+1:,:2) = avgarr(tmax-wind,:2)(-:1:wind,);
  fwdarr(tmax-wind+1:,:2) = fwdarr(tmax-wind,:2)(-:1:wind,);

  plotl,avgarr(,3)*.001,avgarr(,1),color="black",width=4;
  plotl,fwdarr(,3)*.001,fwdarr(,1),color="blue",width=4;

  free = where(!GEOE);
  if(numberof(free)){
    free = free(1);
    //plotl,[1,1]*avgarr(free,3)*.001,minmax,color="magenta";
  }
    
  plxtitle,"Time (ns)";
  plytitle,"Energy / <k_B_T> per molecule";
  
  plotl,avgarr(,3)*.001,avgarr(,1),color="black",width=4;
  plotl,fwdarr(,3)*.001,fwdarr(,1),color="blue",width=4;

  extern FINDUMONDE;
  FINDUMONDE = [avgarr(0,3)*.001,avgarr(0,1)];

  diffs = avgarr(,1)-fwdarr(,1);
  diffs/= avgarr(,2)+1e-9;
  mark = where( abs(diffs)>1.5 );
  if( (numberof(mark)==0)+(mark(0)<free) ){
      mark = where( abs(diffs)>1.5 );
  }
  if( (numberof(mark)==0)+(mark(0)<free) ){
      mark = where( abs(diffs)>.5 );
  }

  nomeq = avgarr(0,3);
  
  if(numberof(mark)){
    mark = mark(0);

    nomeq = avgarr(mark,3);
    nomeq = min(tmax,100.*ceil(.01*nomeq));
    //plotl,[1,1]*avgarr(mark,3)*.001,avgarr(mark,1)+avgarr(mark,2)*[-3.5,3.5],color="green",width=4;
  }
  if(numberof(spec)){
    plot,.001*avgarr(int(spec*2),3),avgarr(int(spec*2),1),color="red";
  }
    
  return [nomeq,avgarr(0,3)];
}

func outcat(void){
  //This is goings to be idiosynchratic to the point of being useless

  topdirs = ls("-d F*");
  for(n=1; n<=numberof(topdirs); n++){
    tempdirs = ls(topdirs(n)+"/");
    for(m=1; m<=numberof(tempdirs); m++){
      outfns = ls("./"+topdirs(n)+"/"+tempdirs(m)+"/*/F*out*");
      parse_out,outfns;
      qq=createb("./"+topdirs(n)+"/"+tempdirs(m)+"/long.out");
      save,qq,TOTALE,POTE,KINE,GEOE,VDWE,TEMP,PRESS,DENS;
      fma;
      energy_width;
      limits,e,e,e,e;
      png,"./"+topdirs(n)+"/"+tempdirs(m)+"/ewidth";
      fma;

      catfns = ls("./"+topdirs(n)+"/"+tempdirs(m)+"/*/cat.b");
      cdd = openb(catfns(1));
      CAT=cdd.CAT;
      DOG=cdd.DOG;
      DIM_ARR=cdd.DIM_ARR;
      for(c=2; c<=numberof(catfns); c++){
        cdd = openb(catfns(c));
        grow,CAT,cdd.CAT;
        grow,DOG,cdd.DOG;
        grow,DIM_ARR,cdd.DIM_ARR;
      }
      qq=createb("./"+topdirs(n)+"/"+tempdirs(m)+"/cat.b");
      save,qq,CAT,DOG,DIM_ARR;
      close,qq;
      CAT=DOG=DIM_ARR=[];

      catfns = ls("./"+topdirs(n)+"/"+tempdirs(m)+"/*/cat_data.b");
      cdd = openb(catfns(1));
      HAT = cdd.HAT;
      cH_arr = cdd.cH_arr;
      cS_arr = cdd.cS_arr;
      oH_arr = cdd.oH_arr;
      oS_arr = cdd.oS_arr;
      for(c=2; c<=numberof(catfns); c++){
        cdd = openb(catfns(c));
        grow,HAT,cdd.HAT;
        grow,cH_arr,cdd.cH_arr;
        grow,cS_arr,cdd.cS_arr;
        grow,oH_arr,cdd.oH_arr;
        grow,oS_arr,cdd.oS_arr;
      }
      qq=createb("./"+topdirs(n)+"/"+tempdirs(m)+"/cat_data.b");
      save,qq,HAT,cH_arr,cS_arr,oH_arr,oS_arr;
      close,qq;
      HAT=cH_arr=cS_arr=oH_arr=oS_arr=[];
    }
  }
}
      


func randosphere(dirs){
  eqts = array(0,numberof(dirs));
  cflag = eqts*0; //complete
  fflag = !cflag; //first
  max_t = int(eqts*0);
  iso_Q = max_t;
  for(n=1; n<=numberof(dirs); n++){
    qq=openb(dirs(n)+"cat.b");
    max_t(n) = dimsof(qq.CAT)(0);
    eqts(n) = qq.eqt;
    iso_Q(n) = oriS_cat(cat=qq.CAT(,,int(qq.eqt):))(avg)<.5;
  }
  //write,max_t;

  timeblock = indgen(max(max_t))(-:1:numberof(dirs),);
  timeblock*= timeblock>=eqts(,-:1:dimsof(timeblock)(0));
  timeblock*= timeblock<=max_t(,-:1:dimsof(timeblock)(0));

  wind = 500;

  extern CAT;
  extern DIM_ARR;

  nbx=nby=342;
  nbz=450;

  animate,1;
  while(anyof(cflag==0)){
    
    for(n=1; n<=numberof(dirs); n++){
      if(cflag(n)) continue;

      ntime = timeblock(n,);
      valid = where(ntime>0);
      picks = valid(sort(random(numberof(valid)))(:min(wind,numberof(valid))));
      qq=openb(dirs(n)+"cat.b");
      CAT=qq.CAT;
      DIM_ARR=qq.DIM_ARR;
      XX=YY=floor(min(DIM_ARR(:2,eqts(n):))*.5);
      XX = min(25,XX);
      YY = min(25,YY);
      ZZ=floor(min(DIM_ARR(3,eqts(n):))*.5);
      if(iso_Q(n)==1){
        nbx=nby=nbz=342;
        XX = floor(min(DIM_ARR(,eqts(n):))*.5);
        XX = min(25,XX);
        ZZ=YY=XX;
      }
        

      
      close,qq;
      if(fflag(n)){
        qq=createb(dirs(n)+"catball.b");
        f_array = g_array = 0(-:1:nbx,-:1:nby,-:1:nbz);
        f_array*=0;
        g_array*=0;
      }
      else{
        qq=updateb(dirs(n)+"catball.b");
        f_array = qq.f_array;
        g_array = qq.g_array;
      }
      write,format="\r Visiting %s (%g%% full)",dirs(n),sum(f_array!=0)*100./(nbx-1)/(nby-1)/(nbz-1);
      window,0,wait=1;
      limits,-XX,XX,-YY,YY;
      tmp = sphere_gr(ntime(picks),XX=XX,YY=YY,ZZ=ZZ,isoflag=iso_Q(n));
      f_array += tmp(,,,1);
      g_array += tmp(,,,2);
      save,qq,f_array,g_array,XX,YY,ZZ;
      close,qq;
      timeblock(n,picks) = 0;
      cflag(n) = noneof(timeblock(n,)>0);
      fflag(n) = 0;
      window,2,wait=1;
      limits,e,e,e,e;
      pli,timeblock;
      fma;
    }
  }
  animate,0;
}
func sphere_gr(times,XX=,YY=,ZZ=,isoflag=){

  nbin = 342; //lies!
  XX = is_void(XX)?20.:XX;
  YY = is_void(YY)?20.:YY;
  ZZ = is_void(ZZ)?20.:ZZ;

  if(isoflag==1){
    ZZ = XX;
    nbx=nby=nbz=342;
  }else{
    nbx=nby=342;
    nbz=450;
  }
  //write,XX,YY,ZZ;
  //write,format="[%g, %g, %g] with (%d, %d, %d) bins\n",XX,YY,ZZ,nbx,nby,nbz;
  x_grid = span(-XX,XX,nbx)(,-:1:nby,-:1:nbz);
  y_grid = span(-YY,YY,nbx)(-:1:nbx,,-:1:nbz);
  z_grid = span(-ZZ,ZZ,nbz)(-:1:nbx,-:1:nby,);
  f_grid = 0(-:1:nbx,-:1:nby,-:1:nbz);
  g_grid = 0(-:1:nbx,-:1:nby,-:1:nbz);
  dx = abs(x_grid(dif,1,1)(1));
  dy = abs(y_grid(1,dif,1)(1));
  dz = abs(z_grid(1,1,dif)(1));
  
  extern DIMS;
  extern COORDS;

  for(t=1; t<=numberof(times); t++){
    time = times(t);
    
    DIMS = DIM_ARR(,time);
    COORDS=CAT(,,time);
  
    for(i=1; i<=NCOL*MPC; i++){
      COORDS = pbc( CAT(:3,,time) - CAT(:3,i,time) );

      uni1= where( (abs(COORDS)<[XX,YY,ZZ])(sum,)==3 );
      uni1 = uni1(where(uni1!=i));
      inds1 = int((COORDS(,uni1)+[XX,YY,ZZ])/[dx,dy,dz])+1;
            
      glob_spin,CAT(5,i,time);
      glob_spin,CAT(4,i,time),raxis=2;
      glob_spin,CAT(6,i,time);
      
      uni2 = where( (abs(COORDS)<[XX,YY,ZZ])(sum,)==3 );
      uni2 = uni2(where(uni2!=i));
      inds2 = int((COORDS(,uni2)+[XX,YY,ZZ])/[dx,dy,dz])+1;
      for(j=1; j<=numberof(uni1); j++){
        f_grid(inds1(1,j),inds1(2,j),inds1(3,j))+=1;
      }
      for(j=1; j<=numberof(uni2); j++){
        g_grid(inds2(1,j),inds2(2,j),inds2(3,j))+=1;
      }
    }
    fma;
    tmp = f_grid;
    tmp(:nbx/2+1,,) = g_grid(:nbx/2+1,,);
    pli,log0(tmp(,,int(7*nbz/16):int(9*nbz/16))(,,sum)),XX,YY,-XX,-YY;
    plt,totxt(t),20,20,tosys=1,color="red";
  }
  return [f_grid,g_grid];
}



func brett2template(xyz,coords){
  xyzp = xyz - xyz(,avg);
  coords = coords - coords(,avg);

  moves=indgen(30);
  diffs=norm(coords-xyzp);
  while(anyof(diffs>1.1)){
    ii=where(diffs==max(diffs) );
    jj=sort(norm(coords - xyzp(,ii)))(1);
    xyzpp=xyzp;
    xyzpp(,jj)=xyzp(,ii);
    xyzpp(,ii)=xyzp(,jj);
    diffp = norm(coords-xyzpp);
    if(1){
      fma;
      pli,[moves],-5,-5,5,-4;
      xyzp=xyzpp;
      tmpi=moves(ii);
      moves(ii)=moves(jj);
      moves(jj)=tmpi;
      diffs=diffp;
      plotl,coords(1,),coords(2,),width=4,color="red";
      plotl,xyzp(1,),xyzp(2,),color="blue";
      pause,42;
    }
  }
  return moves;
}


func outster(temps,fnpre,col=){
  col = is_void(col)?"red":col;
  //assumes a specific directory strucutre
  //may be worth generalizing, but have not yet
  window,1,wait=1;
  fma;
  termarr = array(0.,[2,2,numberof(temps)]);
  for(tt=1; tt<=numberof(temps); tt++){
    outs=ls("./"+temps(tt)+"/[1-9]/"+fnpre+"*[1-9].out");
    if( catch(-1) ){
      write,"BOINK!";
      window,1,wait=1;
      continue;
    }
    //write,outs;
    parse_out,outs,concat=1,quiet=1;
    nomeq=int(energy_width());
    pause,1;
    termarr(,tt) = FINDUMONDE;

    window,9,wait=1;
    adens = DENS(2*nomeq:)(avg);
    sdens = sqrt(avg((DENS(2*nomeq:)-adens)^2));
    plot,tonum(temps(tt)),adens,color=col;
    plotl,tonum(temps(tt))*[1,1],adens+sdens*[1,-1],color=col;
    write,format="Density at %g: %g\n",tonum(temps(tt)),DENS(2*nomeq:)(avg);
    window,1,wait=1;
  }
  limits,e,e,e,e;
  lims = limits();
  xat = lims(2) + lims(dif)(1)*.1;
  doplot = where(termarr(1,)!=0.);
  ltemps = temps(doplot);
  termarr = termarr(,doplot);
  yat = span(lims(3),lims(4),numberof(ltemps));
  for(ll=1; ll<=numberof(ltemps); ll++){
    if(allof(termarr(,ll)==[0.,0.])) continue;
    plt,ltemps(ll),xat,yat(ll),tosys=1,justify="LH";
    plotl,[termarr(1,ll),xat],[termarr(2,ll),yat(ll)],color="black";
  }
  unzoom;
  limits,lims(1),lims(2)+lims(dif)(1)*.2,lims(3)-lims(dif)(3)*.05,lims(4)+lims(dif)(3)*.05;
}


func voroneigh(cms,dims){
  //assumes there's an external cat
  inds = indgen(600);
  vneigh = array(0.,[2,32,600]);
  vweigh = array(0.,[2,32,600]);
  mneigh = array(0,32);

  hdim = totxt(dims*.5);

  tag = totxt(int(random()*100000)+1);

  tmpf = create("tmp.voro"+tag);
  write,tmpf,format="%d %g %g %g\n",inds,cms(1,),cms(2,),cms(3,);
  close,tmpf;
  cstr = "voro++ -g -p -c \"%i     %n %f\"  -";
  cstr+=hdim(1)+" "+hdim(1)+" -"+hdim(2)+" "+hdim(2)+" -"+hdim(3)+" "+hdim(3)+" tmp.voro"+tag;
  system,cstr;
  pp = rdfile("tmp.voro"+tag+".vol");
  if(numberof(pp)!=600){
    write,format="ERR: %s","loose mol";
    return "";
  }
  inds = int(tonum(strpart(pp,:4)));
  pp = strpart(pp,4:);
  tmpstr = array(0.,64);
  for(mol=1; mol<=600; mol++){
    sread,pp(mol),format="%g ",tmpstr;
    nnnum = numberof(where(tmpstr))/2;
    vneigh(:nnnum,inds(mol)) = tmpstr(:nnnum);
    vweigh(:nnnum,inds(mol)) = tmpstr(nnnum+1:2*nnnum);
    vweigh(,inds(mol))/=sum(vweigh(,inds(mol))); // a/A
  }
  system,"rm tmp.voro"+tag+".vol";
  system,"rm tmp.voro"+tag+".gnu";
  system,"rm tmp.voro"+tag;


  return [vneigh,vweigh];
}

func hexOS(void){
  coms = aa_cm();
  cens = aa_rc();
  
  voro = voroneigh(cens,DIMS); //neighbors via voro++, padded w zeds
  out = array(0.,MOLS);
  for(i=1; i<=MOLS; i++){
    vni = voro(,i,1);
    wgt = voro(,i,2);
    wgt = wgt(where(wgt)>.001); //trim lil faces
    vni = int(vni(where(vni)));
    //write,"";
    //write,vni;
      
    dr = pbc( coms(,vni) - coms(,i) );

    paradr = para2(dr,[0,0,1.]);
    perpdr = dr - paradr;
    paradr = abs(paradr(1,),paradr(2,),paradr(3,));
    perpdr = abs(perpdr(1,),perpdr(2,),perpdr(3,));

    Rij = pbc(cens(,vni) - cens(,i) );
    Rij = perp2( Rij , [0,0,1.] );
    Rij/= norm(Rij)(-:1:3,); //unit vectors between centroids

    if(anyof(perpdr<4.)) wgt(where(perpdr<4.)) = 0.;
    wgt = wgt>0.;
    hex_run = 0.;
      
    for(j=1; j<numberof(vni); j++){
      d_jk = Rij(+,j)*Rij(+,j+1:);
      d_jk = min(d_jk, 1.0000);
      d_jk = max(d_jk,-1.0000);
      thjk = acos(d_jk);
      hex_run += sum( cos(6.*thjk) * wgt(j)*wgt(j+1:) );
    }
    hex_run /= (sum(wgt)^2 - sum(wgt) )/2;
    out(i) = hex_run;
    }

  return out(avg);
}

    
func splice_frames(fn1,fn2,out){
  read_frame,fn1;
  c1=COORDS;
  d1=DIMS;
  cm1 = aa_cm();
  read_frame,fn2;
  c2=COORDS;
  d2=DIMS;
  cm2 = aa_cm();

  DIMS = .5*(d1+d2); //easy

  //Distance from cm to aa for all atoms
  c2a1 = c1*0.;
  for(a=1; a<=30; a++) c2a1(,a::30) = cm1;
  c2a2 = c2*0.;
  for(a=1; a<=30; a++) c2a2(,a::30) = cm2;

  c2a = .5*(c2a1+c2a2);

  dm = .5*(c2a1-c1)+(c2a2-c2);

  for(i=1; i<=3; i++){//check for border crossers
    ww = where( (abs(c2a(i,)-c2a1(i,))>DIMS(i)*.25)*(abs(c2a2(i,)-c2a(i,))>DIMS(i)*.25) );
    if(numberof(ww)){
      c2a(i,ww) -= DIMS(i)*.5*sign( c2a(i,ww)-c2a1(i,ww));
    }
  }

  COORDS = c2a + dm;

  write_xyz,out;
}
  
  

func many2one(void){
  COORDS = COORDS(,:30);
  BONDS = BONDS(,:30);
  NMS = NMS(:30);
  IDS = IDS(:30);
  IDK = IDK(:30);
  MASS = MASS(:30);
  COORDS-=COORDS(,avg);
  gmole_spin,1;
  fma;
  bond_plot;
}



func read_charges(key){
  qq=rdfile(key);
  cc = where(strpart(qq,:7)=="CHARGE ");
  chg=array(0.,30);
  atm=indgen(30);
  sread(strpart(qq(cc),8:),format="%d %g",atm,chg);
        //chg = tonum(strpart(qq(cc),11:));
  extern CHG;
  CHG = chg;
  return chg;
}

func chg_pot(void){
  if(numberof(NMS)>30) return "no thanks";
  if(is_void(CHG)) return "what charges?";
  xx = span(-10,10,501)(,-:1:501);
  yy = span(-10,10,501)(-:1:501,);
  pot = xx*0.;
  for(aa=1; aa<=30; aa++){
    x = COORDS(1,aa);
    y = COORDS(2,aa);
    r = sqrt( (xx-x)^2 + (yy-y)^2 );
    r = max(1e-6,r);
    pot += CHG(aa)/r;
  }
  write,format="%d\r",aa;
  plfc,pot,yy,xx,levs=span(-1,1,21);
  return [pot,yy,xx];
}

func vdw_pot(void){
  if(numberof(NMS)>30) return "no thanks";
  va=vb=array(0.,30);
  va(where(IDK==659)) = 3.5500;
  va(where(IDK==661)) = 3.5500;
  va(where(IDK==90))  = 3.5500;
  va(where(IDK==660)) = 2.8500;
  va(where(IDK==662)) = 2.8500;
  va(where(IDK==91))  = 2.4200;

  vb(where(IDK==659)) = 0.0700;
  vb(where(IDK==661)) = 0.0700;
  vb(where(IDK==90))  = 0.0700;
  vb(where(IDK==660)) = 0.0610;
  vb(where(IDK==662)) = 0.0610;
  vb(where(IDK==91))  = 0.0300;

  xx = span(-10,10,501)(,-:1:501);
  yy = span(-10,10,501)(-:1:501,);
  pot = xx*0.;
  for(aa=1; aa<=30; aa++){
    x = COORDS(1,aa);
    y = COORDS(2,aa);
    r = sqrt( (xx-x)^2 + (yy-y)^2 );
    r = max(1e-6,r);
      //4*vb*( (va/r)^12 - (va/r)^6 )
    ir = va(aa)/r;
    r6 = ir^6;
    r12= r6*r6;
    pot += 4*vb(aa)*(r12-r6);
  }
  write,format="%d\r",aa;
  plfc,pot,yy,xx,levs=span(-1,0,10);
  return [pot,yy,xx];
}

func double_pot(void){
  pot = chg_pot();
  fma;
  vdw_pot;
  plc,pot(,,1),pot(,,2),pot(,,3),levs=[-1,-.1,0,.1,1],marks=0;
}
  
