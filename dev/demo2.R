ag=40
age_dia=19
woma=1
eth=0
smok=0
sb=127
hba1=11
ld=143
hd=47.73
weigh=70
tall=1.60
wbc=8
egfr=70
pvd=0
haem=15
hr=70



time=20
complicaciones=8


#### flag

micalb=0
atfib=0
mii=0
strok=0
ih=0
ch=0
blin=0
ampu=0
rena=0
ulce=0


#history=c(age,age_diab,woman,ethn,smoke,sbp,hba1c,ldl,hdl,weight,talla,age-age_diab,weight/(talla^2),sbp/10,ulcer,mi,stroke,ihd,
#chf,blind,amput,renal,1,0,0,0,0,0,0,0,0)

len_his=51


hba1c_red1_comp1_mean=-1.6
hba1c_red1_comp1_se=0.2

weight_red1_comp1_mean=2.30
weight_red1_comp1_se=0.2

sbp_red1_comp1_mean=2.9
sbp_red1_comp1_se=0.2

ldl_red1_comp1_mean=-6.4
ldl_red1_comp1_se=0.2

hdl_red1_comp1_mean=0
hdl_red1_comp1_se=0.2

sbp_red1_comp1_mean=sbp_red1_comp1_mean/10
sbp_red1_comp1_se=sbp_red1_comp1_se/100


stan_data=list(ag=ag,age_dia=age_dia,woma=woma,eth=eth,smok=smok,sb=sb,hba1=hba1,ld=ld,hd=hd,wbc=wbc,egfr=egfr,pvd=pvd,haem=haem,hr=hr,micalb=micalb,atfib=atfib,
               weigh=weigh,tall=tall,ulce=ulce,mii=mii,strok=strok,ih=ih,ch=ch,blin=blin,
               ampu=ampu,rena=rena,time=time,complicaciones=complicaciones,hba1c_red1_comp1_mean=hba1c_red1_comp1_mean,
               hba1c_red1_comp1_se=hba1c_red1_comp1_se,weight_red1_comp1_mean=weight_red1_comp1_mean,weight_red1_comp1_se=weight_red1_comp1_se,
               sbp_red1_comp1_mean=sbp_red1_comp1_mean,sbp_red1_comp1_se=sbp_red1_comp1_se,ldl_red1_comp1_mean=ldl_red1_comp1_mean,ldl_red1_comp1_se=ldl_red1_comp1_se,
               hdl_red1_comp1_mean=hdl_red1_comp1_mean,hdl_red1_comp1_se=hdl_red1_comp1_se,len_his=len_his)

library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)


#m=stan_model(file = '/home/antonio/diabetes2.stan')
m=stan_model(file = '/home/antonio/stanPkgTest/inst/stan/diabetes_update.stan')


fit <- sampling(m, data=stan_data,
                iter=10000, warmup=0, chains=1,
                algorithm="Fixed_param")



stan_fit <- extract(fit)
summary(stan_fit$update)
