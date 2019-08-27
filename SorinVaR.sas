data this_needs_to_be_VaRed;
   Freq     = 2;
   SevMean  = 10000;
   SevMax   = 75000;
   MaxLevel = .9;
   output;

   Freq     = 1;
   SevMean  = 20000;
   SevMax   = 21000;
   MaxLevel = .95;
   output;

   Freq     = 5;
   SevMean  = 10000;
   SevMax   = 175000;
   MaxLevel = .90;
   output;
run;

%macro VaR_Me_This_Mon_Homme(inputDS, simCount, varLevel, ouputDS);
   data &ouputDS.;
      set inputDS.;

      /* Figure pout the parameters */
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
      
      do i = 1 to &simCount.;
         total_loss = 0;
         simFreq = rand("POISSON", Freq);
      
         do j = 1 to simFreq;
            total_loss = total_loss + RAND('WEIBULL', a, b);
         end;
      end;
   run;
%mend VaR_Me_This_Mon_Homme;


data _NULL_;
   do x= 1 to 30 do;
      y = probit(x/31);
      put y=;
   end;
run;
