functions {
    
  real ro_rng(real mu, real sigma) {
    real p = normal_cdf(0, mu, sigma);   // cdf for lb
    real u = uniform_rng(p, 1);
    real z = inv_Phi(u);
    real y = mu + sigma * z;
    return y;
  }
  
  real[] simCHF(real[] update, real intercept, real ro, real age_diab, 
  real atfib, real bmi, real egfr60l, real ld, real mmalb, real pvd, 
  real amp_hist, real ulcer, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+bmi*temp[13]+10*ld*temp[11]+
    amp_hist*temp[25]+ulcer*temp[24]+atfib*temp[20]+
    (egfr60l*(temp[15]<60?temp[15]:0))/10+mmalb*temp[19]+pvd*temp[16];
    
    H1=exp(predictor)*temp[7]^(ro);
    H2=exp(predictor)*(temp[7]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[29]=1;
      temp[33]=1;
    }
    return temp;
  }
  
  real[] simIHD(real[] update,real intercept,real ro,real age_diab,real female,
  real egfr,real hdl,real ldl,real pvd,real sbp,real amp_hist,real chf,
  real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+10*ldl*temp[11]+amp_hist*temp[25]+
    chf*temp[29]+(sbp*temp[9])/10+10*hdl*temp[12]+ female*temp[3]+
    (egfr*temp[15])/10+pvd*temp[16];
    H1=exp(predictor)*temp[7]^(ro);
    H2=exp(predictor)*(temp[7]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[28]=1;
      temp[34]=1;
    }
    return temp;
  }
  
  real[] simMIm(real[] update,real intercept,real afro,real age_diab,real asian,
  real hb1ac,real hdl,real ldl,real mmalb,real pvd, real sbp,real smoker,
  real wbc,real amp,real chf,real ihd,real stroke,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+10*ldl*temp[11]+amp*temp[25]+
    chf*temp[29]+hb1ac*temp[10]+(sbp*temp[9])/10+10*hdl*temp[12]+smoker*temp[8]+
    ihd*temp[28]+stroke*temp[27]+afro*(temp[4]==2?1:0)+asian*(temp[4]==1?1:0)+
    mmalb*temp[19]+pvd*temp[16]+wbc*temp[14];
    
    H1=exp(predictor)*temp[7];
    H2=exp(predictor)*(temp[7]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[21]=1;
      temp[26]=1;
      temp[35]=1;
    }
    return temp;
  }
  
  real[] simMIf(real[] update,real intercept,real ro,real afro,real age_diab,
  real egfr60l,real hb1ac,real ldl35,real mmalb,real pvd,real sbp,real smoker,
  real wbc,real chf,real ihd, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+10*ldl35*(temp[11]>3.5?temp[11]:0)+
    chf*temp[29]+(sbp*temp[9])/10+hb1ac*temp[10]+smoker*temp[8]+ihd*temp[28]+
    afro*(temp[4]==2?1:0)+(egfr60l*(temp[15]<60?temp[15]:0))/10+mmalb*temp[19]+
    pvd*temp[16]+wbc*temp[14];
    
    H1=exp(predictor)*temp[7]^(ro);
    H2=exp(predictor)*(temp[7]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[21]=1;
      temp[26]=1;
      temp[35]=1;
    }
    return temp;
  }
  
  real[] simMI2(real[] update,real intercept, real ldl,real mmalb,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+10*ldl*temp[11]+mmalb*temp[19];
    H1=exp(predictor)*temp[7];
    H2=exp(predictor)*(temp[7]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[21]=2;
      temp[26]=1;
      temp[35]=1;
    }
    return temp;
  }
  
  real[] simStroke(real[] update,real intercept,real ro,real age_diab,
  real female,real atfib,real egfr60l,real hb1ac,real ldl,real mmalb,real sbp,
  real smoker,real wbc,real amp,real ihd,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+10*ldl*temp[11]+
    (sbp*temp[9])/10+female*temp[3]+hb1ac*temp[10]+smoker*temp[8]+ihd*temp[28]+
    atfib*temp[20]+amp*temp[25]+(egfr60l*(temp[15]<60?temp[15]:0))/10+
    mmalb*temp[19]+wbc*temp[14];
    
    H1=exp(predictor)*temp[7]^(ro);
    H2=exp(predictor)*(temp[7]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[22]=1;
      temp[27]=1;
      temp[36]=1;
    }
    return temp;
  }
  
  real[] simStroke2(real[] update,real intercept,real ro,real age_diab,
  real mmalb,real smoker,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+smoker*temp[8]+mmalb*temp[19];
    H1=exp(predictor)*temp[7]^(ro);
    H2=exp(predictor)*(temp[7]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[22]=2;
      temp[27]=1;
      temp[36]=1;
    }
    return temp;
  }
  
  real[] simBlind(real[] update,real intercept, real age_diab,real hb1ac,
  real hr, real sbp, real wbc,real chf, real ihd, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+hb1ac*temp[10]+(hr*temp[18])/10+
    wbc*temp[14]+(sbp*temp[9])/10+chf*temp[29]+ihd*temp[28];
    
    H1=exp(predictor)*temp[7];
    H2=exp(predictor)*(temp[7]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[30]=1;
    }
    return temp;
  }
  
  real[] simUlcer(real[] update,real intercept, real age_diab,real female,
  real bmi, real hb1ac,real pvd, real rand){
    real predictor;
    real H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+hb1ac*temp[10]+pvd*temp[16]+
    bmi*temp[13]+female*temp[3];
    H1=exp(-predictor)/(1+exp(-predictor));
    prob=1 - H1;
    if (prob>rand){
      temp[24]=1;
    }
    return temp;
  }
  
  real[] simAmp_noUlcer(real[] update,real intercept,real ro, real age_diab, 
  real female, real atfib,real hb1ac, real hdl, real hr, real mmalb, real pvd, 
  real sbp, real wbc, real stroke, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]++(sbp*temp[9])/10+female*temp[3]+
    hb1ac*temp[10]+wbc*temp[14]+atfib*temp[20]+10*hdl*temp[12]+(hr*temp[18])/10+
    mmalb*temp[19]+pvd*temp[16]+stroke*temp[27];
    
    H1=exp(predictor)*temp[7]^(ro);
    H2=exp(predictor)*(temp[7]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[23]=1;
      temp[25]=1;
      temp[37]=1;
    }
    return temp;
  }

  real[] simAmp_Ulcer(real[] update,real intercept, real age_diab,real pvd,
  real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+pvd*temp[16];
    H1=exp(predictor)*temp[7];
    H2=exp(predictor)*(temp[7]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[23]=1;
      temp[25]=1;
      temp[37]=1;
    }
    return temp;
  }
  
  real[] simAmp2(real[] update,real intercept, real hb1ac,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+hb1ac*temp[10];
    H1=exp(predictor)*temp[7];
    H2=exp(predictor)*(temp[7]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[23]=2;
      temp[25]=1;
      temp[37]=2;
    }
    return temp;
  }
  
  real[] simRenal(real[] update,real intercept, real afro,real age_diab,
  real female,real bmi, real egfr60l,real egfr60m, real haem, real ldl,
  real mmalb, real sbp,real wbc,real amp, real blind, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+bmi*temp[13]+10*ldl*temp[11]+
    amp*temp[25]+blind*temp[30]+(egfr60l*(temp[15]<60?temp[15]:0))/10+
    (egfr60m*(temp[15]>=60?temp[15]:0))/10+female*temp[3]+afro*(temp[4]==2?1:0)+
    (sbp*temp[9])/10+haem*temp[17]+mmalb*temp[19]+wbc*temp[14];
    H1=exp(predictor)*temp[7];
    H2=exp(predictor)*(temp[7]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[31]=1;
      temp[38]=1;
    }
    return temp;
  }
  
  real[] simDeath_nohist(real[] update,real intercept, real phi,real female,
  real smoker, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+female*temp[3]+smoker*temp[8];
    H1=exp(predictor)*(exp(phi*(temp[1]))-1)/phi;
    H2=exp(predictor)*(exp(phi*(temp[1]+1))-1)/phi;
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[32]=0;
    }
    return temp;
  }
  
  real[] simDeath1st(real[] update,real intercept,real asian,real year,real age,
  real hr,real pvd,real smoker,real amp_event,real ihd_event,real mi_event,
  real renal_event,real stroke_event,real rand){
    real predictor;
    real H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+year*temp[7]+age*temp[1]+smoker*temp[8]+
    amp_event*(temp[37]==1?1:0)+ihd_event*temp[34]+mi_event*temp[35]+
    asian*(temp[4]==1?1:0)+pvd*temp[16]+renal_event*temp[38]+
    stroke_event*temp[36]+(hr*temp[18])/10;
    
    H1=exp(-predictor)/(1+exp(-predictor));
    prob=1 - H1;
    if (prob>rand){
      temp[32]=0;
    }
    return temp;
  }
  
  real[] simDeath_hist(real[] update,real intercept,real phi,real bm1,real bm3,
  real mmalb,real smoker,real wbc, real amp,real chf,real renal,real stroke,
  real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+bm1*(temp[13]<18.5?temp[13]:0)+
    bm3*(temp[13]>25?temp[13]:0)+smoker*temp[8]+amp*temp[25]+chf*temp[29]+
    renal*temp[31]+stroke*temp[27]+mmalb*temp[19]+wbc*temp[14];
    
    H1=exp(predictor)*(exp(phi*(temp[1]))-1)/phi;
    H2=exp(predictor)*(exp(phi*(temp[1]+1))-1)/phi;
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[32]=0;
    }
    return temp;
  }
  
  real[] simDeath_sub(real[] update,real intercept,real atfib,real age,real hdl,
  real pvd,real wbc,real amp_event,real amp,real amp_event2,real ihd_event,
  real ihd,real mi_event,real mi,real renal,real stroke_event,real rand){
    real predictor;
    real H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age*temp[1]+10*hdl*temp[12]+amp_event*(temp[37]==1?1:0)+
    amp*temp[25]+ihd_event*temp[34]+amp_event2*(temp[37]==2?1:0)+ihd*temp[28]+
    mi_event*temp[35]+mi*temp[26]+renal*temp[31]+stroke_event*temp[36]+
    atfib*temp[20]+pvd*temp[16]+wbc*temp[14];
    
    H1=exp(-predictor)/(1+exp(-predictor));
    prob=1 - H1;
    if (prob>rand){
      temp[32]=0;
    }
    return temp;
  }
}

data {
  //int<lower=0> N;
  int<lower=1> complicaciones;
  int<lower=1> time;
  int<lower=1> len_his;
  // real history[len_his];
  // real costs_input[11];
  real ag;       //current age of subject
  real age_dia;  //age at diagnostic
  real woma;      //if its a woman 1 for yes
  real eth;       //ethinicity 1 asian, 2 african, 0 white or other
  real smok;      //if smoker 1 for yes
  real sb;       //sbp
  real hba1;     //Hba1C
  real ld;       //LDL
  real hd;       //HDL
  real wbc;      //white blood cell count
  real haem;     //haemoglobine
  real hr;       //heart rate
  real weigh;    //weight
  real tall;     //height
  real atfib;     //if atrial fibrilation
  real egfr;      //egfr
  real micalb;    //albuminuria
  real pvd;       //if peripherial vascular disease
  real mii;      //mi history
  real strok;    //stroke history
  real ih;       //ihd history
  real ch;       //chf history
  real blin;     //blindness history
  real ampu;     //amount of amputation
  real rena;     //renal failure
  real ulce;     //ulcer hhistory
  
  real hba1c_red1_comp1_mean;
  real<lower=0> hba1c_red1_comp1_se;
  real weight_red1_comp1_mean;
  real<lower=0> weight_red1_comp1_se;
  real sbp_red1_comp1_mean;
  real<lower=0> sbp_red1_comp1_se;
  real ldl_red1_comp1_mean;
  real<lower=0> ldl_red1_comp1_se;
  real hdl_red1_comp1_mean;
  real<lower=0> hdl_red1_comp1_se;
  // real hba1c_red1_mean;
  // real<lower=0> hba1c_red1_se;
  // real weight_red1_mean;
  // real<lower=0> weight_red1_se;
  // real sbp_red1_mean;
  // real<lower=0> sbp_red1_se;
  // real ldl_red1_mean;
  // real<lower=0> ldl_red1_se;
  // real hdl_red1_mean;
  // real<lower=0> hdl_red1_se;
}

transformed data {
  real egfr60m;
  real egfr60l;
  real afro_carib;
  real indian_asian;
  real bm;
  real bmi1;
  real bmi3;
  real diab_years;
  real ldl35;
  real alive;
  real chf_even;
  real ihd_even;
  real mi_even;
  real stro_even;
  real amp_even;
  real renal_even;
  real mii_hist;
  real strok_hist;
  real amp_hist;
  real t_sbp_red1_comp1_mean;
  real t_sbp_red1_comp1_se;
  
  // real t_sbp_red1_mean;
  // real t_sbp_red1_se;
  t_sbp_red1_comp1_mean=sbp_red1_comp1_mean/10;
  t_sbp_red1_comp1_se=sbp_red1_comp1_se/100;
  // t_sbp_red1_mean=sbp_red1_mean/10;
  // t_sbp_red1_se=sbp_red1_se/100;
  
  alive = 1;
  chf_even = 0;
  ihd_even = 0;
  mi_even = 0;
  stro_even = 0;
  amp_even = 0;
  renal_even = 0;

  bm = weigh/(tall^2);
  diab_years = ag-age_dia;
  
  if (strok != 0){
    strok_hist = 1;
  }
  else {
    strok_hist = 0;
  }
  
  if (mii != 0){
    mii_hist = 1;
  }
  else {
    mii_hist = 0;
  }
  
  if (ampu != 0){
    amp_hist = 1;
  }
  else {
    amp_hist = 0;
  }
}

generated quantities {
  real mat_temp[time,len_his];
  int order[complicaciones]={1,2,3,4,5,6,7,8};
  simplex[complicaciones] theta=[0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]';
  real history[len_his];
  real update[len_his];
  real primera=0;
  history[1]=ag;                //current age of subject
  history[2]=age_dia;           //age diagnostic
  history[3]=woma;              //if its a woman 1 for yes
  history[4]=eth;               // ethinicity
  history[5]=weigh;
  history[6]=tall;
  history[7]=diab_years;        //time with diab
  history[8]=smok;               // if smoker 1 for yes
  history[9]=sb;                 // sbp
  history[10]=hba1;               //Hba1C
  history[11]=ld;                 //LDL 
  history[12]=hd;                 //HDL
  history[13]=bm;               //bmi
  history[14]=wbc;               //white cell
  history[15]=egfr;              //efgr
  history[16]=pvd;               //pvd
  history[17]=haem;              //haemoglobine
  history[18]=hr;                //heart rate
  history[19]=micalb;
  history[20]=atfib;             //if atrial fibrilation
  history[21]=mii;               //mii amount
  history[22]=strok;             //stroke amount
  history[23]=ampu;              //amputation amount
  history[24]=ulce;              //ulcer outcome
  history[25]=amp_hist;          //hitories
  history[26]=mii_hist;           
  history[27]=strok_hist;
  history[28]=ih;                //ihd outcome
  history[29]=ch;                // chf outcome  
  history[30]=blin;              // blindeness
  history[31]=rena;              //renal failure
  history[32]=alive;             // if alive 1 if yes
  history[33]=chf_even;         //CHF event in current cycle
  history[34]=ihd_even;         //IHD " "
  history[35]=mi_even;          // MI event
  history[36]=stro_even;        // stroke event
  history[37]=amp_even;         // amputation even
  history[38]=renal_even;       //renal failure event
  history[39]=smok;               // if smoker 1 for yes
  history[40]=sb;                 // sbp initial
  history[41]=hba1;               //Hba1C initial
  history[42]=ld;                 //LDL initial
  history[43]=hd;                 //HDL initial
  history[44]=bm;               //bmi initial
  history[45]=wbc;               //white cell initial
  history[46]=egfr;              //efgr initial
  history[47]=pvd;               //pvd initial
  history[48]=haem;              //haemoglobine initial
  history[49]=hr;                //heart rate initial
  history[50]=micalb;            //albumia initial
  history[51]=atfib;             //if atrial fibrilation initial

  update=history;
  for (i in 1:time){
    if (update[32]==0) {
      break; // en algun mumento tengo que cambiar esta condicion
    }
    for (j in 1:complicaciones){
      int temp;
      int k;
      k=categorical_rng(theta);
      temp=order[j];
      order[j]=order[k];
      order[k]=temp;
    }
    for (s in 1:complicaciones){

      if (order[s]==1){
        if (update[29]!= 1){
          real intercept;
          real ro;
          real age_diab;
          real bmi;
          real ldl;
          real amp;
          real ulcer;
          real atfi;
          real egf60l;
          real mmalb;
          real pv;
          real rand;
          intercept=normal_rng(-12.332,0.859);
          ro=ro_rng(1.514,0.096);
          age_diab=normal_rng(0.068,0.008);
          bmi=normal_rng(0.072,0.008);
          ldl=normal_rng(0.012,0.005);
          amp=normal_rng(0.658,0.334);
          ulcer=normal_rng(0.654,0.291);
          atfi=normal_rng(1.562,0.245);
          egf60l=normal_rng(-0.22,0.065);
          mmalb=normal_rng(0.771,0.116);
          pv=normal_rng(0.479,0.136);
          rand=uniform_rng(0,1);
          update=simCHF(update,intercept,ro,age_diab,atfi,bmi,egf60l,ldl,mmalb,
          pv,amp,ulcer,rand);
        }
    }
    else if (order[s]==2){
      if (update[28]!=1){
          real intercept;
          real ro;
          real age_diab;
          real chf;
          real ldl;
          real amp;
          real sbp;
          real hdl;
          real female;
          real efgr;
          real pv;
          real rand;
          intercept=normal_rng(-6.709,0.503);
          ro=ro_rng(1.276,0.059);
          age_diab=normal_rng(0.016,0.005);
          ldl=normal_rng(0.023,0.003);
          amp=normal_rng(0.526,0.266);
          chf=normal_rng(0.824,0.168);
          sbp=normal_rng(0.058,0.019);
          hdl=normal_rng(-0.065,0.014);
          female=normal_rng(-0.532,0.085);
          efgr=normal_rng(-0.053,0.023);
          pv=normal_rng(0.486,0.181);
          rand=uniform_rng(0,1);
          update=simIHD(update,intercept,ro,age_diab,female,efgr,hdl,ldl,pv,sbp,
          amp,chf,rand);
          }
    }
    else if (order[s]==3){
        if (update[26]==0){
          if (update[3]==0){
            real intercept;
            real age_diab;
            real chf;
            real ldl;
            real amp;
            real hba1c;
            real sbp;
            real ihd;
            real hdl;
            real smoker;
            real stroke;
            real afro;
            real asian;
            real mmalb;
            real pv;
            real wb;
            real rand;
            intercept=normal_rng(-8.791,0.486);
            afro=normal_rng(-0.830,0.237);
            asian=normal_rng(0.279,0.126);
            hba1c=normal_rng(0.108,0.023);
            age_diab=normal_rng(0.045,0.006);
            ldl=normal_rng(0.023,0.004);
            amp=normal_rng(0.743,0.241);
            chf=normal_rng(0.814,0.195);
            sbp=normal_rng(0.046,0.022);
            hdl=normal_rng(-0.049,0.016);
            ihd=normal_rng(0.846,0.117);
            smoker=normal_rng(0.277,0.091);
            stroke=normal_rng(0.448,0.173);
            mmalb=normal_rng(0.023,0.094);
            pv=normal_rng(0.34,0.111);
            wb=normal_rng(0.026,0.013);
            rand=uniform_rng(0,1);
            update=simMIm(update,intercept,afro,age_diab,asian,hba1c,hdl,ldl,
            mmalb,pv,sbp,smoker,wb,amp,chf,ihd,stroke,rand);
          }
          else {
            real intercept;
            real age_diab;
            real chf;
            real ldl;
            real hba1c;
            real sbp;
            real ihd;
            real smoker;
            real ro;
            real afro;
            real egf60l;
            real mmalb;
            real pv;
            real wb;
            real ld35;
            real rand;
            intercept=normal_rng(-8.708,0.844);
            ro=ro_rng(1.376,0.097);
            hba1c=normal_rng(0.078,0.030);
            age_diab=normal_rng(0.041,0.008);
            ldl=normal_rng(0.035,0.007);
            chf=normal_rng(0.853,0.200);
            sbp=normal_rng(0.056,0.027);
            ihd=normal_rng(0.876,0.163);
            smoker=normal_rng(0.344,0.138);
            afro=normal_rng(-1.684,0.506);
            egf60l=normal_rng(-0.28,0.062);
            ld35=normal_rng(0.035,0.007);
            mmalb=normal_rng(0.277,0.129);
            pv=normal_rng(0.469,0.132);
            wb=normal_rng(0.07,0.023);
            rand=uniform_rng(0,1);
            update=simMIf(update,intercept,ro,afro,age_diab,egf60l,hba1c,
            ld35,mmalb,pv,sbp,smoker,wb,chf,ihd,rand);
          }
        }
       else if (update[26]==1){
         real intercept;
         real ldl;
         real mmalb;
         real rand;
         intercept=normal_rng(-4.179,0.262);
         ldl=normal_rng(0.021,0.007);
         mmalb=normal_rng(0.344,0.162);
         rand=uniform_rng(0,1);
         update=simMI2(update,intercept,ldl,mmalb,rand);
       }
    }
      else if (order[s]==4){
        if (update[27]==0){
          real intercept;
          real age_diab;
          real ldl;
          real hba1c;
          real sbp;
          real ihd;
          real smoker;
          real ro;
          real female;
          real amp;
          real mmalb;
          real egf60l;
          real wb;
          real atfi;
          real rand;
          intercept=normal_rng(-13.053,0.722);
          ro=ro_rng(1.466,0.081);
          hba1c=normal_rng(0.092,0.026);
          age_diab=normal_rng(0.066,0.007);
          female=normal_rng(-0.420,0.098);
          ldl=normal_rng(0.016,0.004);
          sbp=normal_rng(0.170,0.022);
          ihd=normal_rng(0.876,0.163);
          smoker=normal_rng(0.331,0.111);
          amp=normal_rng(1.090,0.240);
          mmalb=normal_rng(0.42,0.101);
          egf60l=normal_rng(-0.19,0.056);
          wb=normal_rng(0.04,0.012);
          atfi=normal_rng(1.476,0.201);
          rand=uniform_rng(0,1);
          update=simStroke(update,intercept,ro,age_diab,female,atfi,egf60l,
          hba1c,ldl,mmalb,sbp,smoker,wb,amp,ihd,rand);
        }
        else if (update[27]==1){
          real intercept;
          real age_diab;
          real smoker;
          real ro;
          real mmalb;
          real rand;
          intercept=normal_rng(-9.431,1.569);
          ro=ro_rng(1.956,0.291);
          age_diab=normal_rng(0.046,0.018);
          smoker=normal_rng(0.656,0.263);
          mmalb=normal_rng(0.537,0.228);
          rand=uniform_rng(0,1);
          update=simStroke2(update,intercept,ro,age_diab,mmalb,smoker,rand);
        }
      }
      else if (order[s]==5){
        if (update[30]!=1){
          real intercept;
          real age_diab;
          real hba1c;
          real sbp;
          real ihd;
          real chf;
          real heart;
          real wb;
          real rand;
          intercept=normal_rng(-11.607,0.759);
          hba1c=normal_rng(0.171,0.032);
          age_diab=normal_rng(0.047,0.009);
          sbp=normal_rng(0.068,0.032);
          ihd=normal_rng(0.610,0.208);
          chf=normal_rng(0.841,0.287);
          heart=normal_rng(0.08,0.039);
          wb=normal_rng(0.052,0.019);
          rand=uniform_rng(0,1);
          update=simBlind(update,intercept,age_diab,hba1c,heart,sbp,wb,chf,ihd,
          rand);
        }
      }
      else if (order[s]==6){
        real intercept;
        real age_diab;
        real female;
        real bmi;
        real hba1c;
        real rand;
        real pv;
        intercept=normal_rng(-4.179,0.262);
        age_diab=normal_rng(0.021,0.007);
        female=normal_rng(0.021,0.007);
        bmi=normal_rng(0.021,0.007);
        hba1c=normal_rng(0.021,0.007);
        pv=normal_rng(0.968,0.258);
        rand=uniform_rng(0,1);
        update=simUlcer(update,intercept,age_diab,female,bmi,hba1c,pv,rand);
      }
      else if (order[s]==7){
        if(update[25]==0){
          if (update[24]==0) {
            real intercept;
            real age_diab;
            real hdl;
            real hba1c;
            real sbp;
            real stroke;
            real ro;
            real female;
            real atfi;
            real heart;
            real mmalb;
            real pv;
            real wb;
            real rand;
            intercept=normal_rng(-14.844,1.205);
            ro=ro_rng(2.067,0.193);
            hba1c=normal_rng(0.248,0.042);
            age_diab=normal_rng(0.023,0.011);
            female=normal_rng(-0.445,0.189);
            hdl=normal_rng(-0.059,0.032);
            sbp=normal_rng(0.086,0.043);
            stroke=normal_rng(1.299,0.245);
            atfi=normal_rng(1.088,0.398);
            heart=normal_rng(0.098,0.05);
            mmalb=normal_rng(0.602,0.18);
            pv=normal_rng(1.01,0.189);
            wb=normal_rng(0.04,0.017);
            rand=uniform_rng(0,1);
            update=simAmp_noUlcer(update,intercept,ro,age_diab,female,atfi,
            hba1c,hdl,heart,mmalb,pv,sbp,wb,stroke,rand);
          }
          else {
            real intercept;
            real age_diab;
            real rand;
            real pv;
            intercept=normal_rng(-0.881,1.390);
            age_diab=normal_rng(-0.065,0.027);
            pv=normal_rng(1.769,0.449);
            rand=uniform_rng(0,1);
            update=simAmp_Ulcer(update,intercept,age_diab,pv,rand);
          }
        }
        else if (update[25]==1){
          real intercept;
          real hb1ac;
          real rand;
          intercept=normal_rng(-3.455,0.565);
          hb1ac=normal_rng(0.127,0.060);
          rand=uniform_rng(0,1);
          update=simAmp2(update,intercept,hb1ac,rand);
        }
      }
      else{
        if (update[31]!= 1){
          real intercept;
          real age_diab;
          real female;
          real bmi;
          real ldl;
          real sbp;
          real amp;
          real blind;
          real efgr;
          real afro;
          real efr60l;
          real efr60m;
          real haemo;
          real mmalb;
          real wb;
          real rand;
          intercept=normal_rng(-13.549,1.480);
          age_diab=normal_rng(-0.029,0.013);
          female=normal_rng(-0.869,0.224);
          bmi=normal_rng(-0.054,0.020);
          ldl=normal_rng(0.027,0.007);
          sbp=normal_rng(0.085,0.047);
          amp=normal_rng(1.108,0.337);
          blind=normal_rng(0.732,0.290);
          afro=normal_rng(0.686,0.324);
          rand=uniform_rng(0,1);
          efr60l=normal_rng(-1.031,0.085);
          efr60m=normal_rng(-0.487,0.136);
          haemo=normal_rng(-0.268,0.061);
          mmalb=normal_rng(1.373,0.236);
          wb=normal_rng(0.029,0.013);
          update=simRenal(update,intercept,afro,age_diab,female,bmi,efr60l,
          efr60m,haemo,ldl,mmalb,sbp,wb,amp,blind,rand);
        }
      }
    }

    if(sum(update[33:38])>=1){
      primera+=1;
    }
    if(sum(update[24:31])==0 && sum(update[33:38])==0){
      real intercept;
      real phi;
      real female;
      real smoker;
      real rand;
      intercept=normal_rng(-10.908,0.304);
      phi=ro_rng(0.098,0.004);
      female=normal_rng(-0.229,0.077);
      smoker=normal_rng(0.379,0.089);
      update=simDeath_nohist(update,intercept,phi,female,smoker,rand);
    }
    else if (primera==1){
      real intercept;
      real year;
      real age;
      real smoker;
      real asian;
      real heart;
      real pv;
      real amp_event;
      real ihd_event;
      real mi_event;
      real renal_event;
      real stroke_event;
      real rand;
      intercept=normal_rng(-6.916,0.591);
      asian=normal_rng(-0.54,0.205);
      year=normal_rng(0.042,0.010);
      age=normal_rng(0.058,0.007);
      smoker=normal_rng(0.444,0.117);
      heart=normal_rng(0.124,0.032);
      pv=normal_rng(0.367,0.13);
      amp_event=normal_rng(-0.734,0.321);
      ihd_event=normal_rng(0.423,0.168);
      mi_event=normal_rng(1.309,0.158);
      renal_event=normal_rng(0.584,0.305);
      stroke_event=normal_rng(0.547,0.176);
      rand=uniform_rng(0,1);
      update=simDeath1st(update,intercept,asian,year,age,heart,pv,smoker,
      amp_event,ihd_event,mi_event,renal_event,stroke_event,rand);
    }
    else if (sum(update[24:31])>=1){
      real intercept;
      real phi;
      real bm1;
      real bm3;
      real smoker;
      real amp;
      real chf;
      real renal;
      real stroke;
      real mmalb;
      real wb;
      real rand;
      intercept=normal_rng(-9.207,0.591);
      phi=ro_rng(0.073,0.007);
      bm1=normal_rng(1.083,0.511);
      bm3=normal_rng(-0.293,0.114);
      smoker=normal_rng(0.374,0.133);
      amp=normal_rng(0.539,0.198);
      chf=normal_rng(0.632,0.130);
      renal=normal_rng(1.150,0.197);
      stroke=normal_rng(0.473,0.122);
      mmalb=normal_rng(0.348,0.107);
      wb=normal_rng(0.048,0.011);
      rand=uniform_rng(0,1);
      update=simDeath_hist(update,intercept,phi,bm1,bm3,mmalb,smoker,wb,amp,
      chf,renal,stroke,rand);
    }
    else {
      real intercept;
      real age;
      real hdl;
      real atfi;
      real pv;
      real wb;
      real amp_event;
      real amp;
      real amp_event2;
      real ihd_event;
      real ihd;
      real mi_event;
      real mi;
      real renal;
      real stroke_event;
      real rand;
      intercept=normal_rng(-4.868,0.828);
      age=normal_rng(0.050,0.010);
      atfi=normal_rng(1.081,0.396);
      pv=normal_rng(0.352,0.178);
      wb=normal_rng(0.089,0.038);
      hdl=normal_rng(0.068,0.030);
      amp_event=normal_rng(-1.267,0.344);
      amp=normal_rng(0.753,0.300);
      amp_event2=normal_rng(-1.727,0.467);
      ihd_event=normal_rng(0.583,0.243);
      ihd=normal_rng(-0.507,0.191);
      mi_event=normal_rng(0.982,0.230);
      mi=normal_rng(0.440,0.186);
      renal=normal_rng(0.961,0.396);
      stroke_event=normal_rng(-0.619,0.246);
      rand=uniform_rng(0,1);
      update=simDeath_sub(update,intercept,atfi,age,hdl,pv,wb,amp_event,amp,
      amp_event2,ihd_event,ihd,mi_event,mi,renal,stroke_event,rand);
    }
    if (update[32]!=0){
      real hba1c_red1_comp1;
      real weight_red1_comp1;
      real ldl_red1_comp1;
      real hdl_red1_comp1;
      real hba1c_red1;
      real weight_red1;
      real ldl_red1;
      real hdl_red1;
      real sbp_red1_comp1;
      real sbp_red1;

      hba1c_red1_comp1=normal_rng(hba1c_red1_comp1_mean,hba1c_red1_comp1_se);
      weight_red1_comp1=normal_rng(weight_red1_comp1_mean,weight_red1_comp1_se);
      ldl_red1_comp1=normal_rng(ldl_red1_comp1_mean,ldl_red1_comp1_se);
      hdl_red1_comp1=normal_rng(hdl_red1_comp1_mean,hdl_red1_comp1_se);
      sbp_red1_comp1=normal_rng(t_sbp_red1_comp1_mean,t_sbp_red1_comp1_se);

      // hba1c_red1=normal_rng(hba1c_red1_mean,hba1c_red1_se)
      // weight_red1=normal_rng(weight_red1_mean,weigth_red1_se)
      // ldl_red1=normal_rng(ldl_red1_mean,ldl_red1_se)
      // hdl_red1=normal_rng(hdl_red1_mean,hdl_red1_se)
      // sbp_red1=normal_rng(t_sbp_red1_mean,t_sbp_red1_se)
      //
      update[10]=history[10]+hba1c_red1_comp1;
      update[5]=history[5]+weight_red1_comp1;
      update[11]=history[11]+ldl_red1_comp1;
      update[12]=history[12]+hdl_red1_comp1;
      update[9]=history[9]+sbp_red1_comp1;
      update[1]+=1;
      update[7]+=1;
      update[13]=(history[5]+weight_red1_comp1)/(history[6])^2;
    }
    mat_temp[i]=update;
    update[33]=0;
    update[34]=0;
    update[35]=0;
    update[36]=0;
    update[37]=0;
    update[38]=0;
  }
}
