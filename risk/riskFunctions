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
   
      function Dual_Block_Minima(x,param)
              label="Dual of mimimum among the given number of independent uniform draws";
   
         if (param GE 1) or (param LE 0) then return(.);
         
      	return(1-(1-x)**param);
      endsub;
   
      function Block_Maxima(x,param)
               label="Dual of maximum among the given number of independent uniform draws";
         if param LE 1 then return(.);
      	return( x**param);
      endsub;
   
      function Wang_Transform_Q(x,param)
         label="Wang Transform parameterized through the given quantile";
   
   	   /* Quantile should be positive to assure convex distortion */
      	if param LE 0 then return(.);
   
      	if x < 1.0 then do;
   	   	if x > 0 then do;
   		   	return(  probnorm(probit(x)-param));
      		end;
      		else do;
   	   		return(0);
   		   end;
      	end;
      	else do;
   	   	return(  1.0);
      	end;
      endsub;
   
      function Wang_Transform_P(x,param)
         label="Wang Transform, parameterized through probability";
   
   	   /* Test ensure the distortion is convex and well defined */
         if param LE 0.5 or param GE 1 then return(.);
   
      	if x < 1.0 then do;
   	   	if x > 0 then do;
   		   	return( probnorm(probit(x)-probit(param)));
      		end;
      		else do;
   	   		return(0);
   		   end;
      	end;
      	else do;
   	   	return( 1.0);
      	end;
      endsub;
   run;
   
   options append =(cmplib = &library..&dataset.);   
%mend riskFunctions;
