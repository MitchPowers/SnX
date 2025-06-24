#include "~/TinkerExp/whisker.i"

func z2_UhlDiff(void,color=,temp=,ou=){
  //uses dyn info in memory, call sVDH_gauss

  ou = is_void(ou)?1:ou;
  temp = is_void(temp)?1.:temp;
  color = is_void(color)?blueish():color;
  nmol = dimsof(z2avg)(2);
  fitout = array(0.,[2,4,nmol]); //[uh_th,uh_zo,uh_dc,dif_dc];

  tmax = int(0.8*dimsof(z2avg)(0)); //paretto Cut
  xamt = int(0.2*dimsof(z2avg)(0));
  tmax = dimsof(z2avg)(0)-100;
  time = 1.*indgen(tmax);


  window,9,wait=1;
  limits,e,e,e,e;

  for(i=1; i<=nmol; i++){
    //plot current msd
    Plotl,,z2avg(i,),color="grayb";
    plotl,,z2avg(i,:tmax),color=reddish();
    plotl,,z2avg(i,)+.5*sqrt(z2var(i,)),color="graya";
    plotl,,z2avg(i,)-.5*sqrt(z2var(i,)),color="graya";

    nlm = levmar(z2avg(i,:tmax), time , z2uhl, [.5,0.,.1], amin=[1./(2*10000),1e-9,1e-9], wgt=1./z2var(i,:tmax));
    plotl,,z2uhl(time,nlm),color="blue",width=2;
    ch_ub = (z2avg(i,:tmax)-z2uhl(time,nlm))^2/z2var(i,:tmax);

    lmf = regress(z2avg(i,:tmax), [time,time^0.], sigy=1./sqrt(z2var(i,:tmax)), ax,va,ch2); //ax+b
    plotl,,time*lmf(1)+lmf(2),color="green",width=2;
    ch_df = (z2avg(i,:tmax)-lmf(1)*time-lmf(2))^2/z2var(i,:tmax);

    ch_ub = avg(ch_ub);
    ch_df = avg(ch_df(xmat:));
     
    if(lmf(1)<0.) ch_df = 999; //triggers if D<0
    if(nlm(1)==1./(2*1000)) ch_ub = 999;

    window,ou,wait=1;
    if(ch_ub < ch_df){
      fitout(:3,i) = nlm;
      plot,nlm(0)/temp,nlm(1)*temp/nlm(0),color=color;
      window,9,wait=1;
      plotl,,z2uhl(time,nlm),color="blue",width=8;
    }
    if(ch_df < ch_ub){
      fitout(4,i) = lmf(1);
      plot,lmf(1)/temp,1e-10+random()*1e-10,color=color;
      window,9,wait=1;
      plotl,,lmf(1)*time+lmf(2),color="green",width=8;
    }
    //pause,10;
  }


  return fitout;
}
func z2uhl(time,params){
  th = params(1);
  zo = params(2); //actually zo^2
  dc = params(3);
  return exp(-2*time*th)*(zo - dc/th) + dc/th;
}
func z2uhl_set(ou){
  ff = faire_lebin();
  colors=["black","red","blue","cyan","magenta","yellow"];
  window,ou,wait=1;
  logxy,1,1;

  levmar_itmax = 1024;

  for(i=numberof(ff(,2)); i>=1; i--){
    temp = tonum(strpart(ff(i,1),4:6));
    qq = openb(ff(i,2)+"cat_dynLAB.b");
    restore,qq;
    dat=z2_UhlDiff(color=colors(i),temp=temp,ou=ou);
    orphs = where(dat(4,));
    fosts = where(dat(3,));
    /*write,format="***%s***\n*  Orph: %g (%d)\n* !Orph: %g (%d)\n*\n*\n",strpart(ff(i,1),:6),dat(4,orphs)(avg),numberof(orphs),dat(3,fosts)(avg),numberof(fosts);

    window,ou+10,wait=1;
    plot,temp,dat(4,orphs)(avg),color="blue";
    plot,temp,dat(3,fosts)(avg),color="red";
    s1 = dat(3,fosts)(sum);
    s1+= dat(4,orphs)(sum);
    s1/=dimsof(z2avg)(2);
    plot,temp,s1,color="grayb";
    */
  
  }
  window,ou,wait=1;
  pltitle,strpart(ff(1,1),:2);
  plxtitle,"D/T (A^2^/psK)";
  plytitle,"D k/kb (K/pd)"; 
}
func msdfit2svdh(dat,rawz){
  //assumes pzden is loaded as external
  idiff = where( dat(4,));
  iorns = where(!dat(4,));

  zarr = rawz - rawz(1)*.5;

  Gans = zarr*0.;
  for(t=1; t<=3; t++){
    time = 10^t;
    for(i=1; i<=600; i++){
      if(dat(4,i)){
        DD = dat(4,i);
        Gans += (1./sqrt(4*DD*pi*time))*exp(-(zarr^2)/(4*DD*time));
      }
      else{
        DD = dat(3,i);
        zo = 0.;//sqrt(dat(2,i));
        th = dat(1,i);
        Gans += sqrt(th/(dd*2*pi*(1-exp(-2*th*time)))) * exp(-(zarr - zo*exp(-th*time))^2/(1-exp(-2*th*time)));
      }
    }
    Gans /= sum(Gans*rawz(1));
    plotl,zarr,pzden(,time),color="black",width=4;
    plotl,zarr,Gans,color="red";
  }
}


func msd_stats(msd_arr,subset=,no_plot=,color=,out=){
  //subset is for examining a, well, subset of the population
  subset = is_void(subset)?indgen(dimsof(msd_arr)(-1)):subset;
  msd_arr = msd_arr(subset,);
  nmol = numberof(subset);
  
  //ensemble averaging
  m1 = msd_arr(avg,);

  //ensemble variance
  m2 = ((msd_arr - m1(-:1:nmol,) )^2 )(avg,);
  s1 = sqrt(m2);
  if(anyof(m1-s1<0.)){
    sug = where(m1-s1<0)(1);
    write,format=" Suggest limiting time at %d\n",sug;
  }

  //ensemble skewness
  m3 = ((msd_arr - m1(-:1:nmol,))^3)(avg,);
  skew = m3/(m2^(3/2.));

  if(is_void(no_plot)){
    color = is_void(color)?"red":color;
    plotl,,m1,width=6,color=color;
    plotl,,m1 + s1,color="black",width=4;
    plotl,,m1 - s1,color="black",width=4;

  }
  if(!is_void(out)){
    fn = create(out);
    time = indgen(numberof(m1));
    write,out,format="Time(ps) msd unc%s","\n";
    write,out,format="%d %g %g\n",time,m1,s1;
  }
  
  return [m1,s1];
  
}
  
