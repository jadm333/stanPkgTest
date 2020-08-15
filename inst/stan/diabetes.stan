functions {
  real ro_rng(real mu, real sigma) {
    real p = normal_cdf(0, mu, sigma);   // cdf for lb
    real u = uniform_rng(p, 1);
    real z = inv_Phi(u);
    real y = mu + sigma * z;
    return y;
}
  real[] simCHF(real[] update,real intercept,real ro, real age_diab,real bmi, real ldl, real amp,
  real ulcer, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+bmi*temp[13]+ldl*temp[8]+amp*temp[21]+ulcer*temp[15];
    H1=exp(predictor)*temp[12]^(ro);
    H2=exp(predictor)*(temp[12]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[19]=1;
      temp[24]=1;
    }
    return temp;
  }
  real[] simIHD(real[] update,real intercept,real ro, real age_diab, real ldl, real amp,
  real chf, real sbp,real hdl,real female, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+ldl*temp[8]+amp*temp[21]+chf*temp[19]+sbp*temp[14]+hdl*temp[9]+female*temp[3];
    H1=exp(predictor)*temp[12]^(ro);
    H2=exp(predictor)*(temp[12]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[18]=1;
      temp[25]=1;
    }
    return temp;
  }
  real[] simMIm(real[] update,real intercept,real hb1ac, real age_diab, real ldl, real amp,
  real chf, real sbp,real hdl, real smoker,real ihd,real stroke,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+ldl*temp[8]+amp*temp[21]+chf*temp[19]+sbp*temp[14]+hdl*temp[9]+hb1ac*temp[7]+
    smoker*temp[5]+ihd*temp[18]+stroke*temp[17];
    H1=exp(predictor)*temp[12];
    H2=exp(predictor)*(temp[12]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[16]=1;
      temp[26]=1;
    }
    return temp;
  }
  real[] simMIf(real[] update,real intercept,real ro,real hb1ac, real age_diab, real ldl,
  real chf, real sbp, real smoker,real ihd,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+ldl*(temp[8]>35 ? 1 : 0)+chf*temp[19]+sbp*temp[14]+hb1ac*temp[7]+
    smoker*temp[5]+ihd*temp[18];
    H1=exp(predictor)*temp[12]^(ro);
    H2=exp(predictor)*(temp[12]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[16]=1;
      temp[26]=1;
    }
    return temp;
  }
  real[] simMI2(real[] update,real intercept, real ldl,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+ldl*temp[8];
    H1=exp(predictor)*temp[12];
    H2=exp(predictor)*(temp[12]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[16]=2;
      temp[26]=1;
    }
    return temp;
  }
  real[] simStroke(real[] update,real intercept,real ro,real hb1ac, real age_diab, real ldl,
  real sbp, real smoker,real ihd, real female,real amp,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+ldl*temp[8]+sbp*temp[14]+female*temp[3]+hb1ac*temp[7]+
    smoker*temp[5]+ihd*temp[18]+amp*temp[21];
    H1=exp(predictor)*temp[12]^(ro);
    H2=exp(predictor)*(temp[12]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[17]=1;
      temp[27]=1;
    }
    return temp;
  }
  real[] simStroke2(real[] update,real intercept,real ro,real age_diab, real smoker,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+smoker*temp[5];
    H1=exp(predictor)*temp[12]^(ro);
    H2=exp(predictor)*(temp[12]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[17]=2;
      temp[27]=1;
    }
    return temp;
  }
  real[] simBlind(real[] update,real intercept, real age_diab,real hb1ac, real sbp, real chf,
  real ihd, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+hb1ac*temp[7]+sbp*temp[14]+
    chf*temp[19]+ihd*temp[18];
    H1=exp(predictor)*temp[12];
    H2=exp(predictor)*(temp[12]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[20]=1;
      temp[28]=1;
    }
    return temp;
  }
  real[] simUlcer(real[] update,real intercept, real age_diab,real female, real bmi, real hb1ac,
  real rand){
    real predictor;
    real H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+hb1ac*temp[7]+bmi*temp[13]+
    female*temp[3];
    H1=exp(-predictor)/(1+exp(-predictor));
    prob=1 - H1;
    if (prob>rand){
      temp[15]+=1;
      temp[29]=1;
    }
    return temp;
  }
  real[] simAmp_noUlcer(real[] update,real intercept,real ro,real hb1ac, real age_diab, real hdl,
  real sbp, real stroke, real female,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+hdl*temp[9]+sbp*temp[14]+female*temp[3]+hb1ac*temp[7]+
    stroke*temp[17];
    H1=exp(predictor)*temp[12]^(ro);
    H2=exp(predictor)*(temp[12]+1)^(ro);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[21]=1;
      temp[30]=1;
    }
    return temp;
  }
  real[] simAmp_Ulcer(real[] update,real intercept, real age_diab,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2];
    H1=exp(predictor)*temp[12];
    H2=exp(predictor)*(temp[12]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[21]=1;
      temp[30]=1;
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
    predictor=intercept+hb1ac*temp[7];
    H1=exp(predictor)*temp[12];
    H2=exp(predictor)*(temp[12]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[21]=2;
      temp[30]=1;
    }
    return temp;
  }
  real[] simRenal(real[] update,real intercept, real age_diab,real female,real bmi, real ldl, real sbp,real amp,
  real blind, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age_diab*temp[2]+bmi*temp[13]+ldl*temp[8]+amp*temp[21]+blind*temp[20]+
    female*temp[3]+sbp*temp[6];
    H1=exp(predictor)*temp[12];
    H2=exp(predictor)*(temp[12]+1);
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[22]=1;
      temp[31]=1;
    }
    return temp;
  }
  real[] simDeath_nohist(real[] update,real intercept, real phi,real female,real smoker, real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+female*temp[3]+smoker*temp[5];
    H1=exp(predictor)*(exp(phi*(temp[1]))-1)/phi;
    H2=exp(predictor)*(exp(phi*(temp[1]+1))-1)/phi;
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[23]=0;
    }
    return temp;
  }
  real[] simDeath1st(real[] update,real intercept,real year,real age,real smoker,real amp_event,real ihd_event,real mi_event,
            real renal_event,real stroke_event,real rand){
    real predictor;
    real H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+year*temp[12]+age*temp[1]+smoker*temp[5]+amp_event*temp[30]+ihd_event*temp[25]+mi_event*temp[26]+
    renal_event*temp[31]+stroke_event*temp[27];
    H1=exp(-predictor)/(1+exp(-predictor));
    prob=1 - H1;
    if (prob>rand){
      temp[23]=0;
    }
    return temp;
  }

  real[] simDeath_hist(real[] update,real intercept,real phi,real bm1,real bm2,real smoker,real amp,real chf,real renal,
            real stroke,real rand){
    real predictor;
    real H1;
    real H2=H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+bm1*(temp[13]<18.5?1:0)+bm2*(temp[13]>25?1:0)+smoker*temp[5]+amp*temp[21]+chf*temp[19]+
    renal*temp[22]+stroke*temp[17];
    H1=exp(predictor)*(exp(phi*(temp[1]))-1)/phi;
    H2=exp(predictor)*(exp(phi*(temp[1]+1))-1)/phi;
    prob=1 - exp(H1 - H2);
    if (prob>rand){
      temp[23]=0;
    }
    return temp;
  }
  real[] simDeath_sub(real[] update,real intercept,real age,real hdl,real amp_event,real amp,real amp_event2,real ihd_event,real ihd,
            real mi_event,real mi,real renal,real stroke_event,real rand){
    real predictor;
    real H1;
    real prob;
    real temp[size(update)];

    temp=update;
    predictor=intercept+age*temp[1]+hdl*temp[9]+amp_event*(temp[30]==1?1:0)+amp_event2*(temp[21]==2?1:0)+amp*(temp[21]==1?1:0)+
    ihd*temp[18]+mi_event*temp[26]+mi*temp[16]+renal*temp[22]+stroke_event*temp[27];
    H1=exp(-predictor)/(1+exp(-predictor));
    prob=1 - H1;
    if (prob>rand){
      temp[23]=0;
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
  real ag;
  real age_dia;
  real woma;
  real eth;
  real smok;
  real sb;
  real hba1;
  real ld;
  real hd;
  real weigh;
  real tall;
  real mii;
  real strok;
  real ih;
  real ch;
  real blin;
  real ampu;
  real rena;
  real ulce;


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
  real t_sbp_red1_comp1_mean;
  real t_sbp_red1_comp1_se;
  // real t_sbp_red1_mean;
  // real t_sbp_red1_se;
  t_sbp_red1_comp1_mean=sbp_red1_comp1_mean/10;
  t_sbp_red1_comp1_se=sbp_red1_comp1_se/100;
  // t_sbp_red1_mean=sbp_red1_mean/10;
  // t_sbp_red1_se=sbp_red1_se/100;
}



generated quantities {
  real mat_temp[time,len_his];
  int order[complicaciones]={1,2,3,4,5,6,7,8};
  simplex[complicaciones] theta = [0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125]';
  real history[len_his];
  real update[len_his];
  history[1]=ag;history[2]=age_dia;history[3]=woma;history[4]=eth;history[5]=smok;history[6]=sb;history[7]=hba1;history[8]=ld;
  history[9]=hd;history[10]=weigh;history[11]=tall;history[12]=ag-age_dia;history[13]=weigh/(tall^2);history[14]=sb/10;history[15]=ulce;history[16]=mii;
  history[17]=strok;history[18]=ih;history[19]=ch;history[20]=blin;history[21]=ampu;history[22]=rena;history[23]=1;history[24]=0;
  history[25]=0;history[26]=0;history[27]=0;history[28]=0;history[29]=0;history[30]=0;history[31]=0;
  update=history;
  //mat_temp[1]=history;
  for (i in 1:time){
    real primera=0;
    if (update[23]==0) break; // en algun mumento tengo que cambiar esta condicion
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
        if (update[19]!= 1){
          real intercept;
          real ro;
          real age_diab;
          real bmi;
          real ldl;
          real amp;
          real ulcer;
          real rand;
          intercept=normal_rng(-12.332,0.859);
          ro=ro_rng(1.514,0.096);
          age_diab=normal_rng(0.068,0.008);
          bmi=normal_rng(0.072,0.008);
          ldl=normal_rng(0.012,0.005);
          amp=normal_rng(0.658,0.334);
          ulcer=normal_rng(0.654,0.291);
          rand=uniform_rng(0,1);
          update=simCHF(update,intercept,ro,age_diab,bmi,ldl,amp,
          ulcer,rand);
        }
    }
    else if (order[s]==2){
      if (update[18]!=1){
          real intercept;
          real ro;
          real age_diab;
          real chf;
          real ldl;
          real amp;
          real sbp;
          real hdl;
          real female;
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
          rand=uniform_rng(0,1);
          update=simIHD(update,intercept,ro,age_diab,chf,ldl,amp,
          sbp,hdl,female,rand);
          }
    }
    else if (order[s]==3){
        if (update[16]==0){
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
            real rand;
            intercept=normal_rng(-8.791,0.486);
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
            rand=uniform_rng(0,1);
            update=simMIm(update,intercept,hba1c,age_diab,ldl,amp,chf,
            sbp,hdl,smoker,ihd,stroke,rand);
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
            rand=uniform_rng(0,1);
            update=simMIf(update,intercept,ro,hba1c,age_diab,ldl,chf,
            sbp,smoker,ihd,rand);
          }
        }
       else if (update[16]==1){
         real intercept;
         real ldl;
         real rand;
         intercept=normal_rng(-4.179,0.262);
         ldl=normal_rng(0.021,0.007);
         rand=uniform_rng(0,1);
         update=simMI2(update,intercept,ldl,rand);
       }
    }
      else if (order[s]==4){
        if (update[17]==0){
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
          rand=uniform_rng(0,1);
          update=simStroke(update,intercept,ro,hba1c,age_diab,ldl,
          sbp,smoker,ihd,female,amp,rand);
        }
        else if (update[17]==1){
          real intercept;
          real age_diab;
          real smoker;
          real ro;
          real rand;
          intercept=normal_rng(-9.431,1.569);
          ro=ro_rng(1.956,0.291);
          age_diab=normal_rng(0.046,0.018);
          smoker=normal_rng(0.656,0.263);
          rand=uniform_rng(0,1);
          update=simStroke2(update,intercept,ro,age_diab,smoker,rand);
        }
      }
      else if (order[s]==5){
        if (update[20]!=1){
          real intercept;
          real age_diab;
          real hba1c;
          real sbp;
          real ihd;
          real chf;
          real rand;
          intercept=normal_rng(-11.607,0.759);
          hba1c=normal_rng(0.171,0.032);
          age_diab=normal_rng(0.047,0.009);
          sbp=normal_rng(0.068,0.032);
          ihd=normal_rng(0.610,0.208);
          chf=normal_rng(0.841,0.287);
          rand=uniform_rng(0,1);
          update=simBlind(update,intercept,age_diab,hba1c,sbp,
          chf,ihd,rand);
        }
      }
      else if (order[s]==6){
        real intercept;
        real age_diab;
        real female;
        real bmi;
        real hba1c;
        real rand;
        intercept=normal_rng(-4.179,0.262);
        age_diab=normal_rng(0.021,0.007);
        female=normal_rng(0.021,0.007);
        bmi=normal_rng(0.021,0.007);
        hba1c=normal_rng(0.021,0.007);
        rand=uniform_rng(0,1);
        update=simUlcer(update,intercept,age_diab,female,bmi,hba1c,rand);
      }
      else if (order[s]==7){
        if(update[21]==0){
          if (update[15]==0) {
            real intercept;
            real age_diab;
            real hdl;
            real hba1c;
            real sbp;
            real stroke;
            real ro;
            real female;
            real rand;
            intercept=normal_rng(-14.844,1.205);
            ro=ro_rng(2.067,0.193);
            hba1c=normal_rng(0.248,0.042);
            age_diab=normal_rng(0.023,0.011);
            female=normal_rng(-0.445,0.189);
            hdl=normal_rng(-0.059,0.032);
            sbp=normal_rng(0.086,0.043);
            stroke=normal_rng(1.299,0.245);
            rand=uniform_rng(0,1);
            update=simAmp_noUlcer(update,intercept,ro,hba1c,age_diab,hdl,
            sbp,stroke,female,rand);
          }
          else {
            real intercept;
            real age_diab;
            real rand;
            intercept=normal_rng(-0.881,1.390);
            age_diab=normal_rng(-0.065,0.027);
            rand=uniform_rng(0,1);
            update=simAmp_Ulcer(update,intercept,age_diab,rand);
          }
        }
        else if (update[15]==1){
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
        if (update[22]!= 1){
          real intercept;
          real age_diab;
          real female;
          real bmi;
          real ldl;
          real sbp;
          real amp;
          real blind;
          real rand;
          intercept=normal_rng(3.549,1.480);
          age_diab=normal_rng(-0.029,0.013);
          female=normal_rng(-0.869,0.224);
          bmi=normal_rng(-0.054,0.020);
          ldl=normal_rng(0.027,0.007);
          sbp=normal_rng(0.085,0.047);
          amp=normal_rng(1.108,0.337);
          blind=normal_rng(0.732,0.290);
          rand=uniform_rng(0,1);
          update=simRenal(update,intercept,age_diab,female,bmi,ldl,sbp,amp,
          blind,rand);
        }
      }
    }
    if(sum(update[24:31])>=1)
    primera+=1;
    if(sum(update[15:22])==0 && sum(update[24:31])==0){
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
    else if (primera==1){ // quitar 5 y 6
      real intercept;
      real year;
      real age;
      real smoker;
      real amp_event;
      real ihd_event;
      real mi_event;
      real renal_event;
      real stroke_event;
      real rand;
      intercept=normal_rng(-6.916,0.591);
      year=normal_rng(0.042,0.010);
      age=normal_rng(0.058,0.007);
      smoker=normal_rng(0.444,0.117);
      amp_event=normal_rng(-0.734,0.321);
      ihd_event=normal_rng(0.423,0.168);
      mi_event=normal_rng(1.309,0.158);
      renal_event=normal_rng(0.584,0.305);
      stroke_event=normal_rng(0.547,0.176);
      rand=uniform_rng(0,1);
      update=simDeath1st(update,intercept,year,age,smoker,amp_event,ihd_event,mi_event,
      renal_event,stroke_event,rand);
    }
    else if (sum(update[15:22])>=1){
      real intercept;
      real phi;
      real bm1;
      real bm2;
      real smoker;
      real amp;
      real chf;
      real renal;
      real stroke;
      real rand;
      intercept=normal_rng(-9.207,0.591);
      phi=ro_rng(0.073,0.007);
      bm1=normal_rng(1.083,0.511);
      bm2=normal_rng(-0.293,0.114);
      smoker=normal_rng(0.374,0.133);
      amp=normal_rng(0.539,0.198);
      chf=normal_rng(0.632,0.130);
      renal=normal_rng(1.150,0.197);
      stroke=normal_rng(0.473,0.122);
      rand=uniform_rng(0,1);
      update=simDeath_hist(update,intercept,phi,bm1,bm2,smoker,amp,chf,renal,
      stroke,rand);
    }
    else {
      real intercept;
      real age;
      real hdl;
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
      update=simDeath_sub(update,intercept,age,hdl,amp_event,amp,amp_event2,ihd_event,ihd,mi_event,
      mi,renal,stroke_event,rand);
    }
    if (update[23]!=0){
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
      update[7]=history[7]+hba1c_red1_comp1;
      update[10]=history[10]+weight_red1_comp1;
      update[8]=history[8]+ldl_red1_comp1;
      update[9]=history[9]+hdl_red1_comp1;
      update[14]=history[14]+sbp_red1_comp1;
      update[1]+=1;
    }
    mat_temp[i]=update;
    update[24]=0;
    update[25]=0;
    update[26]=0;
    update[27]=0;
    update[28]=0;
    update[29]=0;
    update[30]=0;
    update[31]=0;

    }
  }

