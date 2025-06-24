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
  MASS=IDK*1/3;
  for(i=1; i<=NATO; i++){
    sread,file(i),format="%d %s %g %g %g %d %d %d %d %d",IDS(i),NMS(i),COORDS(1,i),COORDS(2,i),COORDS(3,i),IDK(i),BONDS(1,i),BONDS(2,i),BONDS(3,i),BONDS(4,i);
  }
  MASS(where(NMS=="CA")) = 12.011;
  MASS(where(NMS=="HA")) = 1.008;
  MASS(where(NMS=="F")) = 18.998;
return COORDS;
}

func write_xyz(output){
  ff=create(output);

  write,ff,format="%d %s\n",NATO,"User Made";
  write,ff,format="%g %g %g %g %g %g\n",DIMS(1),DIMS(2),DIMS(3),90.,90.,90.;
  for(i=1; i<=NATO; i++){
    write,ff,format=" %d %s %g %g %g %d %s\n",IDS(i),NMS(i),COORDS(1,i),COORDS(2,i),COORDS(3,i),IDK(i),sum(totxt(BONDS(where(BONDS(,i)),i))+" "(..));
  }
}


func bond_plot(void,color=,mode=){
  mode = is_void(mode)?[1,2]:mode;
  if(anyof(mode>3)+anyof(mode<1)) mode = [1,2];
  if(noneof(BONDS)) return -1;
  if(color==[]) color="black";
  for(i=1; i<=dimsof(COORDS)(0); i++){
    bonds=BONDS(where(BONDS(,i)),i);
    dummy=array(i,dimsof(bonds));
    bonds=transpose([dummy,bonds]);
    bonds=bonds(where(bonds));
    plotl,COORDS(mode(1),bonds),COORDS(mode(2),bonds),color=color;
  }
}

func atom_plot(void,colors=){
  atoms=swatch=NMS;
  if(colors==[]) colors=["black","grayc","green","blue","magenta"];
  rarr=(MASS*.75/pi)^(1./3);
  for(i=1; numberof(swatch)>=1; i++){
    atom=swatch(1);
    inds=where(atoms==atom);
    plot,COORDS(1,inds),COORDS(2,inds),sym='\1',color=colors(i);
    for(j=1; j<=numberof(inds); j++){
      ll=span(0.,2*pi,12);
      cc=[COORDS(1,inds(j))+rarr(i)*cos(ll),COORDS(2,inds(j))+rarr(i)*sin(ll),COORDS(2,inds(j))];
      plotl,cc(,1),cc(,2),color=colors(i);
    }
    swatch=swatch(where(swatch!=atom));
  }
}

func molec_plot(mol_num,color=){
  coords=COORDS;
  bonds=BONDS; 
  for(i=1; i<=numberof(mol_num); i++){
    COORDS = coords(,NMOL*(mol_num(i)-1)+1:NMOL*mol_num(i));
    BONDS = bonds(,NMOL*(mol_num(i)-1)+1:NMOL*mol_num(i));
    BONDS(where(BONDS))-=min(BONDS(where(BONDS)))-1;
    bond_plot,color=color;
    COORDS = coords;
    BONDS = bonds;
  }
}

func swap_moles(m1,m2){
  /* DOCUMENT swaps position of molecules m1 and m2
     Changes coordinates, not indexes or bonds.
  */
  coords=COORDS;
  ind1=NMOL*(m1-1);
  ind2=NMOL*(m2-1);
  COORDS(,ind1+1:ind1+NMOL)=coords(,ind2+1:ind2+NMOL);
  COORDS(,ind2+1:ind2+NMOL)=coords(,ind1+1:ind1+NMOL);
}

func trim_coord(inds){
  inds=inds(sort(inds));
  cutoff=numberof(inds)*NMOL;
  if(cutoff>NATO) return -999;
  for(i=1; i<=numberof(inds); i++) swap_moles,inds(i),i;
  IDS=IDS(:cutoff);
  NMS=NMS(:cutoff);
  COORDS=COORDS(,:cutoff);
  BONDS=BONDS(,:cutoff);
  IDK=IDK(:cutoff);
}

func check_kern(void,kern=){
  cm=aa_cm()(,1);
  offset=cm;
  COORDS -= offset;
  window,9;
  fma;
  molec_plot(1);
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
  COORDS += offset;
  write,format="%s\n%s\n%s\n","Red lines indicate kernel.","Green vector is director.","Cyan is avg coord.";
}

/********************************************************
 ********************************************************
 ****                                                ****
 ********************************************************
 *******************************************************/

func aa_cm(void){
  n_molec=NATO/NMOL;
  out=array([0.,0.,0.],n_molec);
  for(i=1; i<=n_molec; i++){
    molec=indgen(1+(i-1)*NMOL:i*NMOL);
    mwg = sum(MASS(molec));
    xcm = sum(COORDS(1,molec)*MASS(molec))/mwg;
    ycm = sum(COORDS(2,molec)*MASS(molec))/mwg;
    zcm = sum(COORDS(3,molec)*MASS(molec))/mwg;
    out(,i)=[xcm,ycm,zcm];
  }
  return out;
}

func orients(void,kern=){
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

func thiphi(void){
  ni=orients();
  thi=acos(ni(3,));
  cphi=ni(1,)/sin(thi);
  sphi=ni(2,)/sin(thi);
  phi=acos(cphi)*sign(sphi);
  return transpose([thi,phi]);
}

func forient(void){
  //returns normed vector from central ring through the fluorines
  coords = split(COORDS,NATO/NMOL);
  center = coords(,KERN,)(,avg,);
  fluoro = coords(,BEAC,)(,avg,);
  vector = fluoro - center;
  normed = vector/(norm(vector)(-,,));
  return normed;
}

func compass(void){
  //returns 'phi' angle in lab frame defined as the angle between the
  //fluorinated ring and the X axis
  dir = forient();
  theta = asin(dir(3,));
  plane = dir(:2,)/cos(theta(-,));
  phi = acos(dir(1,)) * sign(dir(2,));
  return phi;
}
  
func rap_sheet(chgfn,bohr=){
  qq=open(chgfn);
  chg=array(0.0,NMOL);
  xx=yy=zz=array(0.0,NMOL);
  inds=array(0,NMOL);
  do{
    line=rdline(qq);
  }while(!is_num(strpart(line,1:1)));
  rd=0;
  for(i=1; i<=NMOL; i++){
    rd+=sread(line,format="%d %g %g %g %g",inds(i),chg(i),xx(i),yy(i),zz(i));
    line=rdline(qq);
  }
  write,format="read %2d charges.\n",rd/5;
  if(rd/5<NMOL) return -999;
  rr=transpose([xx,yy,zz]);
  if(bohr==1){
    rr/=1.8897161646321;
    xx/=1.8897161646321;
    yy/=1.8897161646321;
    zz/=1.8897161646321;
    write,"Converted into Angstroms";
  }
  write,format="Net charge: %g\n",sum(chg);
  neg=where(chg<0.);
  extern ff;
  ff=where(inds==9);
  plot,xx,yy,color="blue";
  plot,xx(neg),yy(neg),color="red";
  plot,xx(ff),yy(ff),color="magenta",sym='\5';
  
  plot,10+zz,yy,color="blue";
  plot,10+zz(neg),yy(neg),color="red";
  plot,10+zz(ff),yy(ff),color="magenta",sym='\5';
  
  plot,xx,10+zz,color="blue";
  plot,xx(neg),10+zz(neg),color="red";
  plot,xx(ff),10+zz(ff),color="magenta",sym='\5';
  write,format="new F's at %d\n",ff;
  limits,e,e,e,e;
  return [xx,yy,zz,chg];
}


func check_key(keyfn){
  key=open(keyfn);
  key=rdfile(key);
  ward=where(strpart(key,:7)=="CHARGE ");

  
  xyzfn=strpart(keyfn,:-4)+".xyz";
  read_frame,xyzfn;
  cm=aa_cm()(,1);
  COORDS=COORDS(,:numberof(ward)) - cm;

  inds=array(0,numberof(ward));
  chgs=array(0.0,numberof(ward));
  key=strpart(key(ward),8:);
  sread,key,format="%d %g",inds,chgs;
  inds=abs(inds);
  neg=where(chgs<0.);
  extern ff;
  ff=where(NMS=="F")(:5);
  xx=COORDS(1,); yy=COORDS(2,); zz=COORDS(3,);
  
  plot,xx,yy,color="blue";
  plot,xx(neg),yy(neg),color="red";
  plot,xx(ff),yy(ff),color="magenta",sym='\5';
  
  plot,10+zz,yy,color="blue";
  plot,10+zz(neg),yy(neg),color="red";
  plot,10+zz(ff),yy(ff),color="magenta",sym='\5';
  
  plot,xx,10+zz,color="blue";
  plot,xx(neg),10+zz(neg),color="red";
  plot,xx(ff),10+zz(ff),color="magenta",sym='\5';
  write,format="Old F at %d\n",ff
  return [xx,yy,zz,chgs];
}


func spin(mat,phi){
  rot = [[cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0,0,1]];
  mat = rot(,+)*mat(+,);
  return mat;
}

func light_brigade(brettfn,mitchfn){
  newchg=transpose(rap_sheet(brettfn,bohr=1));
  f1=ff;
  oldchg=transpose(check_key(mitchfn));
  f0=ff;
  extern coor1; extern coor0;
  coor1= newchg(:3,) - newchg(:3,avg);
  coor0= oldchg(:3,) - oldchg(:3,avg);
  coor0(1,)*=-1; coor0(3,)*=-1;
  dist=array(0.0,30);
  window,7,wait=1;
  animate,1;
  mindist=999.; ang=0.;
  write,numberof(coor1(1,));
  write,numberof(coor0(1,));
  do{
    fma;
    rnd = exp(-random()/.25);
    ang += rnd;
    coor1 = spin(coor1,rnd*pi/180.);
    plot,coor1(1,),coor1(2,),color="red",sym='\2';
    plot,coor1(1,f1),coor1(2,f1),color="red";
    plot,coor0(1,),coor0(2,),color="blue",sym='\3';
    plot,coor0(1,f0),coor0(2,f0),color="blue";
    limits,e,e,e,e;
    for(i=1; i<=numberof(f1); i++){
      dist(i)=min(norm(coor1(:2,f1(i))-coor0(:2,f0)));
    }
    mindist=min(mindist,sum(dist));
    pause,1;
  }while( (ang<360)+( anyof(dist>.05)*sum(dist)>mindist*(1+(ang/360)*.001)) );
  write,ang%360;
  animate,0;

  window,8,wait=1;
  fma;
  dist=array(1,30);
  for(i=1; i<=30; i++){
    dist(i) = sort(norm(coor1(:2,i)-coor0(:2,)))(1);
    plot,i,-newchg(4,i),color="red",sym='\2';
    plot,i,oldchg(4,dist(i)),color="blue",sym='\3';
    plotl,[i,i],[-newchg(4,i),oldchg(4,dist(i))],color="black";
  }
  plot,f1,-newchg(4,f1),color="red";
  plot,f1,oldchg(4,dist(f1)),color="blue";
  limits,e,e,e,e;

  window,7,wait=1;
  animate,1;
  for(n=1; n<=100000; n++){
    roll=random();
    phi=(random()-.5)*(pi/90);
    if(roll<.333){
      rot = [[cos(phi),sin(phi),0],[-sin(phi),cos(phi),0],[0,0,1]];
    }else if(roll<.666){
      rot = [[cos(phi),0,sin(phi)],[0,1,0],[-sin(phi),0,cos(phi)]];
    }else{
      rot = [[1,0,0],[0,cos(phi),-sin(phi)],[0,sin(phi),cos(phi)]];
    }
      coor1q = rot(,+)*coor1(+,);
      if(avg(norm(coor1q-coor0(,dist)))<avg(norm(coor1-coor0(,dist)))){
        coor1=coor1q;
        plot,coor1(1,),coor1(2,),color="red",sym='\2';
        plot,coor1(1,f1),coor1(2,f1),color="red";
        plot,coor0(1,),coor0(2,),color="blue",sym='\3';
        plot,coor0(1,f0),coor0(2,f0),color="blue";
        fma;
        pause,10;
      }
  }
  animate,0;
  fma;
  coor0=coor0(,dist);
  plot,coor1(1,),coor1(2,),color="red",sym='\2';
  plot,coor1(1,f1),coor1(2,f1),color="red";
  plot,coor0(1,),coor0(2,),color="blue",sym='\3';
  plot,coor0(1,f1),coor0(2,f1),color="blue";
   
}
  
func track_dyn(fn,start=,dir=){
  /*vars to track: dimensions*/
  start=(is_void(start))?1:start;
  dx=dy=dz=0.0;
  dout=xout=yout=array(0.0,2*count_files(fn,dir=dir));
  for(i=start; i<=count_files(fn,dir=dir); i++){
    suf=".";
    if(i<100) suf+="0";
    if(i<10)  suf+="0";
    suf+=totxt(i);
    qq=open(fn+suf);
    header=rdline(qq);
    sread,rdline(qq),format="%g %g %g",dx,dy,dz;
    dout(i-start+1) = dz;
    xout(i-start+1) = dx;
    yout(i-start+1) = dy;
    if(i>1) plotl,[i-1,i],dout([i-1,i]),color=reddish();
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
}   


func parse_out(fn,fs,ps,dumps,arc=){
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
  if(is_string(arc)){
    afile=rdfile(arc);
    nato = 1;
    sread,afile(1),format="%d ",nato;
  }
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

    if(is_string(arc)){
      acoords=afile( 1+(i-1)*(nato+1):(i)*(nato+1) );
      read_frame,arcfile=acoords;
      minr=999;
      for(j=1; j<=NATO/2; j++){
        for(k=1; k<=NATO/2; k++){
          dr = COORDS(,j)-COORDS(,30+k);
          qpbc= (abs(dr)>[1,1,1]*30.2/2.);
          if(anyof(qpbc)){
            dr(where(qpbc)) -= 30.2*[1,1,1](where(qpbc))*sign(dr(where(qpbc)));
          }
          dr = norm(dr);
          minr=min(minr,min(dr));
        }
      }
      dr = aa_cm()(,2)-aa_cm()(,1);
      qpbc= (abs(dr)>[30.2,30.2,30.2]/2.);
      if(anyof(qpbc)){
        dr(where(qpbc)) -= [30.2,30.2,30.2](where(qpbc))*sign(dr(where(qpbc)));
      }
      dr = norm(dr)(1);

      if(i>1){
        window,1,wait=1;
        plotl,[ominr,minr],VDWE([i-1,i]),color=reddish();
        //plotl,[odr,dr],VDWE([i-1,i]),color="black";
        plotl,[ominr,minr],POTE([i-1,i])-GEOE([i-1,i])-VDWE([i-1,i]),color=blueish();
        window,9,wait=1;
        limits,-25,25,-25,25;
        COORDS -= aa_cm()(,1);
        fma;
        bond_plot;
      }
      ominr=minr;
      odr=dr;
    }
  }
}


func col_coords(void,ref=){
  cms = split(aa_cm(),NATO/(20*NMOL));
  col = cms(,avg,);
  for(i=1; i<=dimsof(col)(0); i++){
    write,format="Col %d: %#8g  %#8g\n",i,col(1,i),col(2,i);
    if(!is_void(ref)){
      dr = norm( col - col(,ref));
      write,format="      Col %d is %.8g Angstroms from Col %d\n",i,dr(i),ref;
    }
  }
  if(is_void(ref))  return col;
  else return dr;
}

func pbc(scriptr){
  /*DOCUMENT pbc(dr) applies pbc's to the array r_i - r_j
    returns shifted array */
  for(i=1; i<=min(dimsof(DIMS)(2),dimsof(scriptr)(2)); i++){
    fi = anyof(abs(scriptr(i,))>DIMS(i)/2.);
    if(fi){
      wi = where(abs(scriptr(i,))>DIMS(i)/2.);
      scriptr(i,wi) -= DIMS(i)*sign(scriptr(i,wi));
    }
  }
  return scriptr;
}

func stack_em(prefix,fnout,xarr,qarr){
  read_frame,prefix+".xyz";
  
  file=create(fnout);
  write,file,format=" %d Generated from %s",NATO*xarr(1)*xarr(2)*xarr(3),prefix;
  write,file,format="\n %g  %g  %g  90.000 90.000 90.000",(xarr(1))*qarr(2),(xarr(2))*qarr(1),(xarr(3))*qarr(3);
  offset=[0.,0.,0.];
  coords=COORDS;
  for(i=1; i<=xarr(1); i++){
    
    for(j=1; j<=xarr(2); j++){
      
      for(k=1; k<=xarr(3); k++){

        coords=COORDS+offset;
        th=2.*pi*random();   //th is angle about the zaxis
        if(rand==0) th=rand;
        cm=aa_cm();
        coords-=cm+offset;
        coords=[[cos(th),-sin(th),0.],[sin(th),cos(th),0.],[0.,0.,1.]](+,)*coords(+,);
        if(random()>.5) coords=[[1,0,0],[0,-1,0],[0,0,-1]](+,)*coords(+,);
        coords+=2.*(random(3)>.5)*(random(3)-.5);

        coords+=cm+offset;
        
        for(n=1; n<=NATO; n++){ //print molecule at i,j,k
          write,file,format="\n %-d   %s   %06.4g   %06.4g   %06.4g   %d",IDS(n),NMS(n),coords(1,n),coords(2,n),coords(3,n),IDK(n);
          for(m=1; m<=numberof(where(BONDS(,n))); m++){
            write,file,format="   %d",BONDS(m,n);
            BONDS(m,n)+=NATO;
          }
          IDS(n)+=numberof(NMS);
        }
        
        offset+=[0.,0.,qarr(3)];
      }
      offset(3)=0.;
      offset+=[0.,qarr(2),0.];
    }
    offset(2)=0.;
    offset+=[qarr(2)*sqrt(3)/2.,-(i%2)*qarr(2)*.5,0.];
  }
  close,file;
  read_frame,fnout;
  atom_plot;
  bond_plot;
}

func pllr_out(ncol,k=,d=,unmoor=){
  k = is_void(k)?20.0:k;
  d = is_void(d)?1.5:d;
  unmoor = is_void(unmoor)?0:unmoor;

  ncol = int(ncol);
  csom = split(aa_cm(),ncol)(:2,avg,);
  fmts = "RESTRAIN-PLLR    %d    %-6.6g   %-6.6g    ";
  if(unmoor) fmts+="-";
  fmts += "%-4.6g   %.3g\n";

  for(i=1; i<=ncol; i++){
    write,format=fmts,i,csom(1,i),csom(2,i),k,d;
  }
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
    kk = 20.0;
    dr = 1.5;
    sread,file(pllr(i)),format="%s %d %g %g %g %g\n",kword,coli,colx,coly,kk,dr;
    newx = new_coords(1,i);
    newy = new_coords(2,i);
    kk = is_void(unmoor)?kk:-kk;
    file(pllr(i)) = swrite(format="%s   %d   %g   %g   %g   %g","RESTRAIN-PLLR  ",coli,newx,newy,kk,dr);
  }
  new_key = create(fnout);
  for(i=1; i<=numberof(file); i++){
    write,new_key,format="%s\n",file(i);
  }
}
                 

func bin2d(xdat,ydat,&fmesh,&xmesh,&ymesh){
  /*DOCUMENT binner(indep var, dep var, &rhist, &ghist, &fhist)
   */
  if(xmesh==[]) xmesh=1.01*span(floor(min(0.,min(xdat))),ceil(max(abs(xdat))),101)(,-:1:101);
  if(ymesh==[]) ymesh=1.01*span(floor(min(0.,min(ydat))),ceil(max(abs(ydat))),101)(-:1:101,);
  if(fmesh==[]) fmesh=xmesh*0.;
  //plot,xdat,ydat,color="blue",sym='\1';
  for(i=1; i<dimsof(xmesh)(2); i++){
    for(j=1; j<dimsof(ymesh)(3); j++){
      freq=0;
      hits = (xdat>=xmesh(i,1))*(xdat<xmesh(i+1,1))*(ydat>=ymesh(1,j))*(ydat<ymesh(1,j+1));
      freq=sum(hits);
      fmesh(i,j) += freq;
    }
  }
  //  dxdy=(xmesh(2,1)-xmesh(1,1))*(ymesh(1,2)-ymesh(1,1));
  //  plfc,fmesh/sum(fmesh*dxdy),ymesh,xmesh;
  return fmesh;
}


func glob_spin(th){
    rotm = [[cos(th),sin(th),0.],[-sin(th),cos(th),0.],[0,0,1] ];
    COORDS=rotm(+,)*COORDS(+,);
}

func trigrid(void,rad=){
  d0 = is_void(rad)?180:pi;
  d1 = is_void(rad)?120:2*pi/3;
  d2 = is_void(rad)?30:pi/6;
  
  //Major lines
  plotl, d1*[1,1],d0*[-1,1],color="red",width=4;
  plotl,-d1*[1,1],d0*[-1,1],color="red",width=4;
  plotl,  0*[1,1],d0*[-1,1],color="red",width=4;

  plotl,d0*[-1,1], d1*[1,1],color="red",width=4;
  plotl,d0*[-1,1],-d1*[1,1],color="red",width=4;
  plotl,d0*[-1,1],  0*[1,1],color="red",width=4;

  //Minor lines
  plotl, (d1+d2)*[1,1],d0*[-1,1],color="green",width=2;
  plotl, (d1-d2)*[1,1],d0*[-1,1],color="green",width=2;
  plotl,(-d1+d2)*[1,1],d0*[-1,1],color="green",width=2;
  plotl,-(d1+d2)*[1,1],d0*[-1,1],color="green",width=2;
  plotl,      d2*[1,1],d0*[-1,1],color="green",width=2;
  plotl,     -d2*[1,1],d0*[-1,1],color="green",width=2;


  plotl,d0*[-1,1],  (d1+d2)*[1,1],color="green",width=2;
  plotl,d0*[-1,1],  (d1-d2)*[1,1],color="green",width=2;
  plotl,d0*[-1,1], (-d1+d2)*[1,1],color="green",width=2;
  plotl,d0*[-1,1], -(d1+d2)*[1,1],color="green",width=2;
  plotl,d0*[-1,1],       d2*[1,1],color="green",width=2;
  plotl,d0*[-1,1],      -d2*[1,1],color="green",width=2;

}

func check_dims(void){
  dd = DIMS/2.;
  bond_plot,color="grayc";
  plotl,dd(1)*[-1,1,1,-1,-1],dd(2)*[-1,-1,1,1,-1],color="black",width=2;

  COORDS(1,)+=DIMS(1);
  bond_plot,color="grayb";
  plotl,dd(1)*[1,3,3,1,1],dd(2)*[-1,-1,1,1,-1],color="black",width=2;

  COORDS(2,)+=DIMS(2);
  bond_plot,color="grayb";
  plotl,dd(1)*[1,3,3,1,1],dd(2)*[1,1,3,3,1],color="black",width=2;

  COORDS(1,)-=DIMS(1);
  bond_plot,color="grayb";
  plotl,dd(1)*[-1,1,1,-1,-1],dd(2)*[1,1,3,3,1],color="black",width=2;

  COORDS(2,)-=DIMS(2);
}
  
  
func oriS(void,vec_out=){
  //per Soft Matter, 2016, 12, 1295-1312
  ni=orients();
  AA=array(0.,[2,3,3]);
  num=int(numberof(ni)/3);
  for(i=1; i<=num; i++){
    ui=ni(,i);
    AA+=3*ui*ui(-,)-[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
  }
  AA/=2.*num;
  eval=SVdec(AA,evec);
  if(is_void(vec_out)) return max(eval);
  return eval;
}

func parse_traj(prefix,function,start=,stop=,dir=){
  /* DOCUMENT parse_traj(prefix,scalr_function,start=)
     prefix should be the trajectory file prefix
     scalr_function is a function that takes no arguments and
          returns a scalar. All calculations in scalar_function
          should assume the current frame has been read in
          but any non-external variables will need to be calculated
          If a more complex function is desired, consider making function
          interact with it's own external vairables.
     start= is an optional starting frame for use when burning diseq frames
  */
  if(!is_func(function)){
    write,"what's your function?\n";
    return "duh";
  }
  ogcf = count_files(prefix,dir=dir);
  if(is_void(start)) start = 1;
  if(is_void(stop)) stop = ogcf*10;
  offset = start - 1;

  suffix = ".";
  if(start<100) suffix+="0";
  if(start<1000) suffix+="0";
  read_frame,prefix+suffix+totxt(start);
  
  dat_arr = array(function(),count_files(prefix,dir=dir)-offset);
  
  for( i=1+offset; i<=min(stop,count_files(prefix,dir=dir)); i++){
    fmt="\b\b\b\b";
    suffix=".";
    if(i<1000) fmt="\b\b\b";
    if(i<100){ suffix+="0"; fmt="\b\b"; }
    if(i<10) { suffix+="0"; fmt="\b"; }
    write,format=fmt+"%d",i;
    suffix+=totxt(i);
    read_frame,prefix+suffix;

    frame_dat = function();
    dat_arr(i-offset) = frame_dat;

  }

  write,format="\nAvg S: %g",avg(dat_arr);
  write,format="\nVar S: %g\n",sqrt( avg((dat_arr - avg(dat_arr))^2));
  return dat_arr;
}


func mean_var(arr){
  mean = avg(arr);
  var = sqrt( sum((arr-mean)^2)/numberof(arr));
  return [mean,var];
}

