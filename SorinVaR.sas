option mprint;

/************************************************/
/* Sample input data                            */
/* FreqMean     : loss frequency mean           */
/* SevMedian    : median for severity of loss   */
/* SevLikeMax   : severity likely maximal loss  */
/*.               i.e. a fixed high quantile.   */
/* LikeMaxLevel : level used to determine       */
/*                the SevLikeMax quantile.      */
/************************************************/
data models;
   FreqMean       = 2;
   SevMedian      = 100001;
   SevLikeMax     = 750002;
   SevModel       = "LOGN";
   LikeMaxLevel   = .9;
   output;

   FreqMean      = 1;
   SevMedian     = 200003;
   SevLikeMax    = 210004;
   SevModel      = "LOGN";
   LikeMaxLevel  = .95;
   output;

   FreqMean      = 5;
   SevMedian     = 100005;
   SevLikeMax    = 1750006;
   SevModel      = "WBL";
   LikeMaxLevel  = .90;
   output;
run;

/* Distribution model supported  */
/* LOGN : Lognormal              */
/* WBL  : Weibull                */

/* Macros Params_XXX and Sim_XXX */
/* are needed to extend to other */
/* distribution XXX.             */

/* LOGN severity model */
%macro Params_LOGN(v1, v2, p);
    param_1 = log(&v1.);
    param_2 = log(&v2. / &v1.)/probit(&p.);
%mend Params_LOGN;

%macro Sim_LOGN(loss);
    &loss = exp(rand("NORMAL") * Param_2 + Param_1);
%mend SimLOGN;


/* WBL severity model */
%macro Params_WBL(v1, v2, p);
    param_1 = log(log(.5) / log(1 - &p.)) / log(&v1. / &v2.);    
    param_2 = &v1. / (-log(.5))**(1/param_1);
%mend Params_WBL;

%macro Sim_WBL(loss);
    &loss = rand("WEIBULL", Param_2, Param_1);
%mend SimWBL;

/* Some useful FCMP functions */

proc fcmp outlib = work.utils.var;

   /*- Shell sort ------------------*/
   subroutine orva_sort_mat( mat[*,*], col);
      outargs mat;
   
      inc = 1;
      n = dim1(mat);
                /*Determine the starting increment */
      do while ( inc <= n);
         inc = inc*3;
         inc +=1;
      end;
      
      do while (inc > 1);      /* begin the main loop */
   
         inc = floor(inc/ 3);
         do i = inc+1 to n;
            v = mat[i, col];
            j = i;
          t =  mat[j-inc , col];
            do while (t > v);
               mat[j, col] = mat[j-inc, col];
               j = j - inc;
               if j <= inc then goto next;
               t =  mat[j-inc , col];
            end;
      next:      mat[j, col]=v;
         end;
      end;
   endsub;

   /*- Declare distortion functions ----------------------------------*/
   function VaR95(x)
            label="Value at risk at the 95% quantile";
      if x<.95 then do;
         return(0);
      end;
      else do;
         return(1);
      end;
   endsub;
   
   function CVaR95(x)
            label="Mean excess over the VaR at the 95% quantile";
      if x<.95 then do;
         return(0);
      end;
      else do;
         return((x-.95)/(1-.95));
      end;
   endsub;
run;

/*
options append=(cmplib=work.utils);   
*/

%let codeFilePath = c:/junk/dummy.txt;

%macro riskEval(inputDS, simCount, riskFunction, ouputDS);
   data _NULL_;
      rc = filename("tmpDummy","&codeFilePath");
   run;

   data _NULL_;
      set &inputDS. end = last;
      length str $128;
      file tmpDummy;

      str = '%Params_' || strip(SevModel) || '(' || strip(SevMedian) || ',' || strip(SevLikeMax) || ',' || strip(LikeMaxLevel) || ');';
      put str;

      put "do i = 1 to &simCount.;";
      put "   total_loss = 0;";

      str = '   simFreq = rand("POISSON", ' || FreqMean || ');';
      put str;

      put "   ";
      put "   do j = 1 to simFreq;";

      str =     '      %Sim' || strip(SevModel) || '( loss );';
      put str;

      put "      ";
      put "      total_loss = total_loss + loss;";
      put "   end;";
      put "   ;";

      str = '   index = ' || _N_ || ";";
      put str;
      
      str = '   losses[i,' || _N_ || '] = total_loss;";
      put str;

      put "   output;";
      put "end;";
      put "/*-------------------------------------*/;";

      if last then do;
         call symputx("modelCount", _N_ );
      end;
   run;

   data &ouputDS.;
      array losses(&simCount., &modelCount.) _temporary_;
      %include "&codeFilePath";
      /* Figure out the parameters                         */
      
      *&Model._Params(sevMedian, sevLikeMax, LikeMaxLevel);
      /* For lognormal, there is a condition that needs to */
      /* be met, otherwise, the input is meaningless.      */

      /* The condition is that:                            */
      /* Probit(MaxLevel)^2 > 2 x ln(SevMax/SevMean)       */
      /* In that case:                                     */ 
      /* sigma = Probit(MaxLevel) +                        */
      /*         + (-)Sqrt[                                */
      /*                   Probit(MaxLevel)^2              */
      /*                   - 2 ln(SevMax/SevMean)          */
      /*                     ]                             */
      /* mu = ln(SevMax) - sigma Probit(MaxLevel)          */
   run;

   proc sort data = &ouputDS.;
      by i index;
   run;

   proc transpose data = &ouputDS.(keep=loss i index) out=atr(drop=_NAME_) prefix=loss_;
      by i;
      id index;
   run;
%mend riskEval;

%macro Var_without(input, nreps);
      if last then do;      /* Calculate var without   */
         do col = 1 to &num_cells_P1;
            call orva_sort_mat( mat, col );
            risk  = 0;
            lastg = 0;

            do n    = 1 to &nreps;
               Fn   = n/&nreps;
               loss = mat[n, col];
               g = &distortion(Fn);
               g = &distortion(Fn, &distortionparam);

               

               Delta   = g-lastg;
               risk    = risk + loss * delta;
               lastg = g;
            end; /* do n = 1 to &nreps; */

         
      end; /* do col = 1 to &num_cells_P1; */    
   end;   /* if last ---------------*/
run;

%mend;

%VaR_Me_This_Mon_Homme(this_needs_to_be_VaRed, 10, .95, abcd);
