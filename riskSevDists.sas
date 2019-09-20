%macro riskSevDists;
   %put +-----------------------------------------------------------+;   
   %put | Severity distribution models supported:;                  |
   %put |                                       LOGN : Lognormal    |;
   %put |                                       WBL  : Weibull      |;
   %put +-----------------------------------------------------------+;   
%mend riskSevDists;

/* To extend to other to other distributions */
/* macros Params_XXX and Sim_XXX are need to */
/* be written for each model XXX             */ 

/******************************************************************/
/* Lognormal severity model                                       */
/******************************************************************/
  %macro Params_LOGN(v1, v2, p);
      param_1 = log(&v1.);
      param_2 = log(&v2. / &v1.)/probit(&p.);
  %mend Params_LOGN;

  %macro Sim_LOGN(loss);
      &loss = exp(rand("NORMAL") * Param_2 + Param_1);
  %mend SimLOGN;
  
/*  For lognormal, there may be a condition that needs to */
/*  be met, for the input to be meaningful.               */
/*                                                        */
/*  The condition is that:                                */
/*  Probit(MaxLevel)^2 > 2 x ln(SevMax/SevMean)           */
/*                                                        */
/*  In that case:                                         */ 
/*  sigma = Probit(MaxLevel) +                            */
/*          + (-)Sqrt[                                    */
/*                    Probit(MaxLevel)^2                  */
/*                    - 2 ln(SevMax/SevMean)              */
/*                      ]                                 */
/*  mu = ln(SevMax) - sigma Probit(MaxLevel)              */
/* I need to double check that restriction                */  
/******************************************************************/

/******************************************************************/
/* Weibull severity model                                         */
/******************************************************************/
  %macro Params_WBL(v1, v2, p);
      param_1 = log(log(.5) / log(1 - &p.)) / log(&v1. / &v2.);    
      param_2 = &v1. / (-log(.5))**(1/param_1);
  %mend Params_WBL;

  %macro Sim_WBL(loss);
      &loss = rand("WEIBULL", Param_2, Param_1);
  %mend SimWBL;
/*******************************************************************/
