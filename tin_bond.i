#include "~/Tinkel/sn4.i"

if(is_void(COORDS)){
  read_frame,"~/TinkerData/525/data/9/F7_ramp.1234";
 }

func a_to_m(atomic){
  return (atomic-1)/NMOL + 1;
}

func peribonds(void){
  bf = (NMS(:NMOL)=="F");
  wf = where(bf);

  bh = (strpart(NMS(:NMOL),:1)=="H");
  wh = where(bh);

  peris = where(bh+bf);
  roots = BONDS(,:NMOL)(1,peris);

  return [roots,peris];
}

func bond_vec(mol_num, hat=){
  pbonds = peribonds() + (mol_num-1)*NMOL;

  b_vecs = COORDS(,pbonds)(,,dif)(,,1);

  if(is_void(hat)) return b_vecs;
  b_vecs/= norm(b_vecs)(-,);
  return b_vecs;
}






///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////

molecs = NATO/NMOL;
offs = (indgen(molecs)-1)*NMOL;

perii = peribonds()(,2);
all_peri = array(perii,molecs);
all_peri += offs(-,);
all_peri = all_peri(:numberof(all_peri));

rooti = peribonds()(,1);
all_root = array(rooti,molecs);
all_root += offs(-,);
all_root = all_root(:numberof(all_root));




///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////


func close_atoms(mol_num, rcut=){
  rcut = is_void(rcut)?2.95:rcut;

  mmin = (mol_num - 1)*NMOL + 1;
  mmax = mol_num * NMOL;
  
  perii = peribonds()(,2) + mmin - 1;


  
  peri2peri = norm(pbc(COORDS(,-,all_peri) - COORDS(,perii,-)));
  peri2peri(where(peri2peri==0)) = 999; //avoid self ID;

  isneigh = (peri2peri < rcut)*all_peri(-,);
  isneigh*= (isneigh<mmin)+(isneigh>mmax);
  ncount = (!(!isneigh))(,sum);

  if(max(ncount) < 1){
    return array(0, [2,1,numberof(perii)]);
  }

  neigh = array(0, [2,max(ncount),numberof(perii)]);
  for(m=1; m<=numberof(perii); m++){
    tmp = isneigh(m,where(isneigh(m,)));
    if(numberof(tmp)==0) continue;
    //if(anyof(tmp<min(minds))) tmp = tmp(where(tmp<min(minds)));
    //if(anyof(tmp>max(minds))) tmp = tmp(where(tmp>max(minds)));
    neigh(:numberof(tmp),m) = tmp;
  }
  return neigh;
  
}

func close_checks(mol_num,ok=){
  rcut = 2.95;
  dcut = 0.97; //~15 degrees off parallel
  
  neighbors = close_atoms(mol_num,rcut=rcut);

  output = neighbors * 0.0;

  bvi = bond_vec(mol_num,hat=1);

  if(ok==1) molec_plot,mol_num,color="grayd";

  for(n=1; n<=dimsof(neighbors)(0); n++){
    if( noneof(neighbors(,n) ) ){
      continue;
    }
    
    contacts = neighbors( where(neighbors(,n)) ,n);
    for(j=1; j<=numberof(contacts); j++){
      contact = contacts(j);
      rootj = BONDS(1,contact);
      bvj = COORDS(,rootj) - COORDS(,contact);
      bvj/= norm(bvj);
      dotq = abs(bvi(+)*bvj(+));

      if(ok==1){
        if(dotq>=dcut){
          plotl,COORDS(1,[rootj,contact]),COORDS(2,[rootj,contact]),color="green",width=4;
          plotlc,COORDS(1,rootj)+bvj(1),COORDS(2,rootj)+bvj(2),1.,color="green",width=2;
        } else {
          plotl,COORDS(1,[rootj,contact]),COORDS(2,[rootj,contact]),color="black",width=4;
          plotlc,COORDS(1,rootj)+bvj(1),COORDS(2,rootj)+bvj(2),1.,color="black",width=2;
        }
        plotl,COORDS(1,[(mol_num-1)*NMOL + perii(n),contact]),COORDS(2,[(mol_num-1)*NMOL + perii(n),contact]),color="blue",width=1;
        
      }

      output(j,n) = dotq;
    }
  }

  return output;
   
}


func clever_census(mol_num){
  dcut1 = 0.966; //~15 degrees off parallel
  dcut2 = 0.866; //~30 degrees off parallel
  
  neighbors = close_atoms(mol_num);
  /*
  contact = neighbors>0;

  gridij = digitize(neighbors%NMOL,perii)-1  + 12*indgen(0:11)(-,) * contact;
  if(noneof(gridij)) return 0;
  gridij = gridij(where(gridij));

  fgrid1(:numberof(fgrid1))(gridij)++;
  */
    ///
    
  neighbors = close_atoms(mol_num);
  C_i = indgen(12)(-:1:dimsof(neighbors)(2),);
  C_i = C_i(where(neighbors));
  neighbors = neighbors(where(neighbors));
  if(numberof(neighbors)<1) return [0];
  C_j = digitize(neighbors%NMOL,perii)-1;
  fgrid1(:numberof(fgrid1))(C_j + (C_i-1)*12)++;

  
  bvi = bond_vec(mol_num,hat=1)(,C_i);

  atomj = neighbors;
  rootj = BONDS(1,neighbors);
  bvj = COORDS(,rootj) - COORDS(,atomj);
  bvj/= norm(bvj)(-,);

  dotij = abs((bvi*bvj)(sum,));

  if(anyof(dotij>dcut2)){
    ww = where(dotij>dcut2);
    C_i2 = C_i(ww);
    neighbors2 = neighbors(ww);
    C_j2 = digitize(neighbors2%NMOL,perii)-1;
    fgrid2(:numberof(fgrid2))(C_j2 + (C_i2-1)*12)++;

    if(anyof(dotij>dcut1)){
      ww = where(dotij>dcut1);
      C_i3 = C_i(ww);
      neighbors3 = neighbors(ww);
      C_j3 = digitize(neighbors3%NMOL,perii)-1;
      fgrid3(:numberof(fgrid3))(C_j3 + (C_i3-1)*12)++;
    }
  }
  
  return dotij;
  
}

/*
  > pli,(fgrid2)([1,2,3,4,7,8,5,9,10,6,11,12],)(,[1,2,3,4,7,8,5,9,10,6,11,12])
> for(i=1; i<=12; i++) plt,strpart(NMS(perii([1,2,3,4,7,8,5,9,10,6,11,12])),:1)(i),i-.7,.5,tosys=1,color="red"
> for(i=2; i<=12; i++) plt,strpart(NMS(perii([1,2,3,4,7,8,5,9,10,6,11,12])),:1)(i),.3,i-.7,tosys=1,color="red"
*/

  



///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////

func setup_hist(void){
  cpm = 10; //estimated contacts per molecule
  liln = int(sqrt(cpm*NATO/NMOL) + 1);
  bign = int(liln*sqrt(1000)) + 1;
  
  extern d_arr;
  extern f_arr;

  d_arr = cos(span(0,60,liln)*pi/180)(::-1); //parallel to 60deg out
  f_arr = d_arr*0;

  extern fgrid1;
  extern fgrid2;
  extern fgrid3;

  fgrid1 = fgrid2 = fgrid3 = array(0,[2,12,12]);
  

}

func close_hist(void){
  tmp = f_arr;
  for(m=1; m<=NATO/NMOL; m++){
    dat = close_checks(m);
    dat = dat(where(dat));
    fast_bin,dat,d_arr,f_arr;
  }
  out = f_arr - tmp;
  Plotl,d_arr,f_arr;
  limits,e,e,0,e;
  return sum(out)/600.;
}
  

func clever_hist(void){
  tmp = f_arr;
  for(m=1; m<=NATO/NMOL; m++){
    dat = clever_census(m);
    dat = dat(where(dat));
    if(numberof(dat)>1) fast_bin,dat,d_arr,f_arr;
  }
  out = f_arr - tmp;
  Plotl,d_arr,f_arr;
  limits,e,e,0,e;
  return sum(out)/600.;
}
  

///////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////

func setup_census(void){
  cpm = 10; //estimated contacts per molecule
  liln = int(sqrt(cpm*NATO/NMOL) + 1);
  bign = int(liln*sqrt(1000)) + 1;
  
  extern d_arr;
  extern f_arr;

  fij_arr = array(0, [2,12,12])
}

func close_hist(void){
  tmp = f_arr;
  for(m=1; m<=NATO/NMOL; m++){
    dat = close_checks(m);
    dat = dat(where(dat));
    fast_bin,dat,d_arr,f_arr;
  }
  out = f_arr - tmp;
  Plotl,d_arr,f_arr;
  limits,e,e,0,e;
  return sum(out)/600.;
}
  


///////////////////////////////////////////////////////////////////////////////////

func make_plots(void){
  dcut1 = 0.966; //~15 degrees off parallel
  dcut2 = 0.866; //~30 degrees off parallel
 
  window,4,wait=1;
  fma;

  domega = sin(d_arr(zcen))*d_arr(dif);
  prob = (f_arr(:-1)/sum(1.*f_arr)) * (sum(domega) / domega);

  plotl,acos(d_arr(zcen))*180/pi,prob,color="black";
  limits,0,e,0,e;
  limits,0,e,0,ceil(limits()(-1));

  window,5,wait=1;
  fma;

  seq = [1,2,3,4,7,8,5,9,10,6,11,12]; //maps perii to radial order

  f1 = fgrid1(seq,)(,seq);
  f2 = fgrid2(seq,)(,seq);
  f3 = fgrid3(seq,)(,seq);

  coef = 12. / fgrid2(,sum);
  coef2 = domega(where(d_arr>=dcut2)-1)(sum);
  coef3 = domega(where(d_arr>=dcut1)-1)(sum);
  
  rad1 = 1.0; //sin(acos(0.5)); 
  rad2 = sin(acos(dcut2));
  rad3 = sin(acos(dcut1));

  //rcut2 = floor((rad2/rad1)*5*2)/2;
  //rcut3 = floor((rad3/rad1)*5*2)/2;
  rcut2 = 5;
  rcut3 = floor((rad3/rad2)*5*2)/2;

  extern fgrid;
  fgrid = array(0.,[2,120,120]);
  for(i=1; i<=120; i++){
    lili = (i-1)/10 + 1;
    for(j=1; j<=120; j++){
      lilj = (j-1)/10 + 1;
      if( (abs((i-1)%10-4.5)<=rcut3)*(abs((j-1)%10-4.5)<=rcut3) ){
        //fgrid(i,j) = fgrid3(lili,lilj) * coef(lili) * coef2 / coef3;
        fgrid(i,j) = (fgrid3(lili,lilj)*coef2) / (fgrid2(lili,lilj) * coef3);
        continue;
      }
      if( (abs((i-1)%10-4.5)<=rcut2)*(abs((j-1)%10-4.5)<=rcut2) ){
        fgrid(i,j) = fgrid2(lili,lilj) * coef(lili);
        continue;
      }
      //fgrid(i,j) = fgrid1(lili,lilj) * coef(lili);
      
    }
  }
  plfc,(fgrid),span(0,12,120)(,-:1:120),span(0,12,120)(-:1:120,),levs=floor(10.^span(-.5,.5,21)*100)*.01;
  for(i=1; i<=12; i++) plt,strpart(NMS(perii([1,2,3,4,7,8,5,9,10,6,11,12])),:1)(i),i-.5,-.5,tosys=1,color="red",justify="CH";
  for(i=1; i<=12; i++) plt,strpart(NMS(perii([1,2,3,4,7,8,5,9,10,6,11,12])),:1)(i),-.5,i-.5,tosys=1,color="red",justify="CH";
  limits,-1,12.5,-1,12.5;

  
  plfc_levs(11) = 1.0;
  color_bar,plfc_levs,vert=1,adjust=-.0125,labs=[1,5]

  

}
