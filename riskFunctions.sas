%macro riskFunctions(library, dataset);
   /* Some useful FCMP functions */
   proc fcmp outlib = &library..&dataset..risk;
      /*- Shell sort ------------------*/
      subroutine sort_matrix( mat[*,*], col);
         outargs mat;

         inc = 1;
         n   = dim1(mat);

         /*Determine the starting increment */
         do while ( inc <= n);
            inc = inc * 3;
            inc += 1;
         end;

         /* begin the main loop */
         do while (inc > 1);      
            inc = floor(inc / 3);

            do i = inc + 1 to n;
               v = mat[i, col];
               j = i;
               t =  mat[j - inc , col];

               do while (t > v);
                  mat[j, col] = mat[j - inc, col];
                  j = j - inc;

                  if j <= inc then goto next;

                  t =  mat[j - inc , col];
               end;

   next:       mat[j, col] = v;
            end;
         end;
      endsub;

      /*- Declare distortion functions ----------------------------------*/
      function VaR95(x)
               label = "Value at risk at the 95% quantile";
         if x<.95 then do;
            return(0);
         end;
         else do;
            return(1);
         end;
      endsub;
   
      function CVaR95(x)
               label = "Mean excess over the VaR at the 95% quantile";
         if x<.95 then do;
            return(0);
         end;
         else do;
            return((x - .95)/(1 - .95));
         end;
      endsub;
   run;
   
   options append =(cmplib = &library..&dataset.);   
%mend riskFunctions;
