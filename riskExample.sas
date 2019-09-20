
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

/* Instantiate risk functions               */
%riskFunctions;

/* Instantiate severity distribution models */
%riskSevDist;


/* Call the risk evaluation macro */
%riskEval(
          models, 
          10000, 
          VaR95, 
          riskResults
          );
