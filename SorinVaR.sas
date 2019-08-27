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

      for i = 1 to &simCount. do;
         simFreq = rand("POISSON", Freq);
         
         for j = 1 to simFreq do;
            rand(
         end;
      end;

%mend VaR_Me_This_Mon_Homme;
