libname glm "C:\Users\juach\Documents\UHasselt\Generalized Linear Models\Project";

/*Permanent Datasets*/
data glm.mice;
infile "C:\Users\juach\Documents\UHasselt\Generalized Linear Models\Project\EG.dat" firstobs=2;
input id dose response;
z2=0; z3=0; z4=0; 
if dose=0.75 then z2=1; if dose=1.5 then z3=1; if dose=3 then z4=1; 
resp_nonlive=(response=3);
resp_normal=(response=1);
resp_malform=(response ne 1);
dose2=dose*dose;
dosesqrt=sqrt(dose);
dosenew = dose*100;
dosenew2=dosenew*dosenew;

label id="Litter/Cluster"
        dose="Dosage Level g/kg"
        response="Status"
        resp_nonlive="Dead or Not"
        resp_normal="Normal or Not"
        resp_malform="Malformed/Dead or Not";
run;


/*For independent aggregated data*/
proc sql;
create table glm.mice_indep as
	select  dose
			,count(response) as total
			,sum(case when response=3 then 0 else 1 end) as alive
			,sum(resp_normal) as y
			,sum(case when response=1 then 0 else 1 end) as y_notnormal
			,sum(case when response=1 then 1 else 0 end) as y_normal
			,sum(case when response=2 then 1 else 0 end) as y_malform
			,sum(case when response=3 then 1 else 0 end) as y_dead
			,case when dose=0.75 then 1 else 0 end as z2
			,case when dose=1.5 then 1 else 0 end as z3
			,case when dose=3 then 1 else 0 end as z4
			,dose*dose as dose2
			,sqrt(dose) as dosesqrt
			,dose*100 as dosenew
			,(dose*100)*(dose*100) as dosenew2
	from glm.mice
	group by dose;
quit;


/*For clustered data*/
proc sql;
create table glm.mice_cluster as
	select ID
			,dose
			,count(response) as total
			,sum(case when response=3 then 0 else 1 end) as alive
			,sum(resp_normal) as y
			,sum(case when response=1 then 0 else 1 end) as y_notnormal
			,sum(case when response=1 then 1 else 0 end) as y_normal
			,sum(case when response=2 then 1 else 0 end) as y_malform
			,sum(case when response=3 then 1 else 0 end) as y_dead
			,case when dose=0.75 then 1 else 0 end as z2
			,case when dose=1.5 then 1 else 0 end as z3
			,case when dose=3 then 1 else 0 end as z4
			,dose*dose as dose2
			,sqrt(dose) as dosesqrt
	from glm.mice
	group by ID, dose;
quit;

data mice; set glm.mice; run;
data mice_indep; set glm.mice_indep; run;
data mice_cluster; set glm.mice_cluster; run;
proc export data=mice_cluster outfile="C:\Users\juach\Documents\UHasselt\Generalized Linear Models\Project\r_mice.csv" dbms=csv replace; run;
proc export data=mice outfile="C:\Users\juach\Documents\UHasselt\Generalized Linear Models\Project\r_mice_orig.csv" dbms=csv replace; run;


/*****************************Independent Binary**************************************************/

/*if you put descending option: normal response!*/
proc logistic data=mice desc;
model resp_normal = dose/lackfit aggregate scale=none;
run;
/*There is lack of fit as seen in deviance and pearson*/

/*correct for lack of fit using pearson*/
proc logistic data=mice ;
model resp_normal = dose/lackfit aggregate scale=pearson;
run;

/*****************************Clustered Binary**************************************************/

/*Binomial ML*/
proc logistic data=mice_cluster; 
model y_normal/total = dose/aggregate scale=none; run;
/*or:*/
proc logistic data=mice desc; 
model resp_normal= dose/aggregate scale=none; run;
/*There is lack of fit or overdispersion based on pearson or deviance*/

/*Quasi Likelihood*/
proc logistic data=mice_cluster; 
model y_normal/total = dose/aggregate scale=williams; run;
/*weight variable: 0.187005*/
/*The odds of fetus being malformed/dead is almost three(2.903)) times higher per unit increase of dose*/


/************GEE MODELS**************/
/* GEE with indep working correlation*/ 
proc genmod data=mice descending; 
class id; 
model resp_normal = dose / dist=bin link=logit lrci type3; 
repeated subject=id /type=indep modelse; 
run;
/* GEE with exchangeable working correlation*/
proc genmod data=mice descending; 
class id; 
model resp_normal = dose/ dist=bin link=logit lrci type3; 
repeated subject=id /type=exch modelse corrw; 
run;

proc genmod data=mice descending; 
class id; 
model resp_malform = dose/ dist=bin link=logit lrci type3; 
repeated subject=id /type=exch modelse corrw; 
run;
/*Exchangeable Working Correlation Correlation 0.2040437062 */

/************FULL LIKELIHOOD**************/

/* Full likelihood Beta-binomial approach */   /*EDIT!!!*/
/*REFER TO R CODES*/

/************RANDOM EFFECTS MODEL**************/

proc nlmixed data=mice qpoints=40; 
eta1 = i1 + beta1*dose + u; 
p1 = 1/(1 + exp(-eta1)); 
p2 = 1 - 1/(1 + exp(-eta1)); 
select; 
	when (response in (2,3)) y1=1; y2=0; 
	when (response=1) y1=0; y2=1; 
end; 
ll = y1*log(p1) + y2*log(p2); 
model response ~ general(ll); 
random u ~ normal(0,su*su) subject=id; 
run;

proc nlmixed data=mice qpoints=200; 
eta = alpha + beta2*dose + u ; 
p = exp(eta)/(1 + exp(eta)); 
model resp_malform ~ binary(p) ; 
random u ~ normal(0, sigma*sigma) subject=id; 
run;
/***************************************************************************************************
/********************************Trinomial Logit Model**********************************************
/**************************************************************************************************/

/*-------Baseline Logit----------*/
/* Fitting baseline category logit model, with baseline = Normal (has highest ni) */ 
proc logistic data=mice;
model response(ref='1') = dose/ link=glogit aggregate scale=none; 
run;

/*To correct lack of fit/overdispersion*/
proc logistic data=mice;
model response(ref='1') = dose/ link=glogit aggregate scale=pearson; 
run;

/*-------Adjacent Logit----------*/
/*A. Adjacent logit model with equal slopes */ 
proc logistic data=mice;
model response= dose/ link=alogit aggregate scale=none; 
run;
/*A. same model, different sign of parameters*/
proc logistic data=mice descending;
model response= dose/ link=alogit aggregate scale=none; 
run;
/*For Likelihood Ratio Test, check R Codes*/

/*B. Adjacent logit model with (un)equal slopes */ 
proc logistic data=mice descending;
model response= dose/ link=alogit aggregate unequalslopes scale=none; 
run;
/*B. Adjacent logit model with (un)equal slopes with scaling factor=pearson */ 
proc logistic data=mice descending;
model response= dose/ link=alogit aggregate unequalslopes scale=pearson; 
run;

/*-------Cumulative Link----------*/
/*equal slopes*/
proc logistic data=mice;
model response = dose/link=logit aggregate scale=none;
run;
/*unequal slopes*/
proc logistic data=mice descending;
model response = dose/link=logit aggregate unequalslopes scale=none;
run;

/*unequal slopes with scaling pearson*/
proc logistic data=mice descending;
model response = dose/link=logit aggregate unequalslopes scale=pearson;
run;/*REFER TO R CODES for Likelihood Ratio Test*/


/*A. Continuation logit model Equal Slopes*/ 
/*FINAL OUTPUT IN R CODES*/
	/*NONLIVE VS LIVE
proc logistic data=mice descending;
model resp_nonlive = dose/aggregate  scale=none;
run;
	/*MALFORMED VS NORMAL
proc logistic data=mice(where=(response ne 3)) descending;
model resp_malform = dose/aggregate scale=none;
run;
/*B. Continuation logit model Unequal Slopes*/

/*REFER TO R CODES for Adjacent Logit*/
/*REFER TO R CODES for Likelihood Ratio Test*/



/***************************************************************************************************
/*************************Clustered Trinomial Data Model***********************************************
/**************************************************************************************************/

/*GEE*/
proc genmod data=mice ; 
class id; 
model response= dose/ dist=multinomial link=clogit lrci type3 ; 
repeated subject=id/type=indep modelse corrw; 
run;/*FINAL*/

proc gee data=mice;
class id; 
model response= dose/ dist=multinomial link=clogit; 
repeated subject=id/logor=exch; 
run;

proc gee data=mice;
class id; 
model response= dose/ dist=multinomial link=clogit; 
repeated subject=id/logor=indep; 
run;

/*GLMM*/
proc nlmixed data=mice qpoints=40; 
bounds i2>0; 
eta1 = i1 + beta1*dose + u; 
eta2 = i1 + i2 + beta1*dose + u; 
p1 = 1/(1 + exp(-eta1)); 
p2 = 1/(1 + exp(-eta2)) - 1/(1 + exp(-eta1)); 
p3 = 1 - 1/(1 + exp(-eta2)); 
select; 
	when (response=1) y1=1; y2=0; y3=0; 
	when (response=2) y1=0; y2=1; y3=0; 
	when (response=3) y1=0; y2=0; y3=1;  
end; 
ll = y1*log(p1) + y2*log(p2) + y3*log(p3); 
model response ~ general(ll); 
estimate 'thresh2' i1+i2; 
random u ~ normal(0,su*su) subject=id; 
run;





