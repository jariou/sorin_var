option mprint;
data this_needs_to_be_VaRed;
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

%macro Var_without(input, nreps);
   data output; 
      set input end=last;
      array mat[&nreps, 1] _temporary_ ( &NUM_CELLS_NREPS * 0);

      length index 8;
      length key $8;
      length data $8;
      index = data;

      do col = 1 to &num_cells_P1;
         mat[_REP_, col] + loss;
      end;

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

               %if &useCapitalAtRisk = true %then %do;
                  g = g-Fn;
               %end;

               Delta   = g-lastg;
               risk    = risk + loss * delta;
               lastg = g;
            end; /* do n = 1 to &nreps; */

         /* Store the number away so we can compute the marginals*/
         risk_measures[col] = risk;
      end; /* do col = 1 to &num_cells_P1; */
      
      /* Risk measure from the total cell */
      total_risk = risk_measures[&num_cells_P1]; 

       /* compute marginal risks */
      total_marginal = 0;
      do col = 1 to &num_cells;
         risk_measures[col] = total_risk - risk_measures[col] ;
         total_marginal + risk_measures[col];
      end;

      /* Now compute porportions ---------*/
      do col = 1 to &num_cells;
         if total_marginal = 0 then do; 
            risk_measures[col] = 0;
         end;
         else do;
            risk_measures[col] = &totalproportion * risk_measures[col]/total_marginal;
         end;
      end;

      /* output proportions ---------------*/
      keep CCUID proportion;

      do col = 1 to &num_cells;
         CCUID = scan("&CELLS", col);
         proportion = risk_measures[col];
         output;
      end;
   end;   /* if last ---------------*/
run;

%mend;

%macro SimWBL(loss);
    &loss = rand("WEIBULL", Param_2, Param_1);
%mend SimWBL;

%macro SimLOGN(loss);
    &loss=exp(rand("NORMAL") * Param_2 + Param_1);
%mend SimLOGN;

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
   
   function _VaR_(x,param)
           label="Value at risk at the user given quantile";
   
      if param LE 0 or param GE 1 then return(.);
   
      if x < param then do;
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
   
   function CVaR(x,param)
            label="Mean excess over the VaR at the user provided quantile";
   
      if param LE 0 or param GE 1 then return(.);
   
      if x < param then do;
         return( 0);
      end;
      else do;
         return((x-param)/(1-param));
      end;
   endsub;
   

run;

/*
options append=(cmplib=work.utils);   
*/

%let codeFilePath = c:/junk/dummy.txt;

%macro VaR_Me_This_Mon_Homme(inputDS, simCount, varLevel, ouputDS);
   data _NULL_;
      rc = filename("tmpDummy","&codeFilePath");
   run;

   data _NULL_;
      set &inputDS.;
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
      /*
      put "   losses[i] = total_loss;";
      */
      put "   output;";
      put "end;";
      put "/*-------------------------------------*/;";
   run;

   data &ouputDS.;
      array losses(&simCount.) _temporary_;
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

   proc sort data=&ouputDS.;
      by i index;
   run;

   proc transpose data= &ouputDS.(keep=loss i index) out=atr(drop=_NAME_) prefix=loss_;
      by i;
      id index;
   run;
%mend VaR_Me_This_Mon_Homme;

%macro Params_LOGN(v1, v2, p);
    param_1 = log(&v1.);
    param_2 = log(&v2. / &v1.)/probit(&p.);
%mend Params_LOGN;

%macro Params_WBL(v1, v2, p);
    param_1 = log(log(.5) / log(1 - &p.)) / log(&v1. / &v2.);    
    param_2 = &v1. / (-log(.5))**(1/param_1);
%mend Params_WBL;

%VaR_Me_This_Mon_Homme(this_needs_to_be_VaRed, 10, .95, abcd);
