/*#############################################################################
  
  # BayClone2 :  A Bayesian Sequence and Copy Number Caller for Subclones Using 
  #              NGS Data
  # Copyright (C) <2014>  <Juhee Lee>
  # This file is part of BayClone2.

  #  BayClone2 is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.

  #  BayClone2 is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.

  #  You should have received a copy of the GNU General Public License
  #  along with BayClone2. If not, see <http://www.gnu.org/licenses/>.
  
#############################################################################*/

#include "stdio.h"
#include "R.h"
#include "Rmath.h"
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

//#include <stddef.h>
//#include <veclib/clapack.h>
//#include "CT_rho_fn.h"

#define MAX_LINE_LEN	50000
#define MAX_FLD_LEN		500
#define SAMPLE_OFFSET	9
#define FORMAT_INDX		8

//R CMD SHLIB intratumor.cpp
extern "C"
{

	int readNextLine (char *l, FILE *f)
	{
 		if (fgets (l, MAX_LINE_LEN, f) == NULL) 
		{
 			/* NULL.  EOF was read.  reset l, return */
 			l[0] = '\0';
 	   	 	return 1;
    	}
 		else
		{ 
			if (strlen (l) == MAX_LINE_LEN -1) 
			{
    			/* Line too long - do some malloc stuff.
    	 		 * For now, just terminate and return. */
     	 		error("Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
 	     		l[0] = '\0';
 	     		return 2;
    		}
			else    
			{
 				/* Proper line read. Return 0 */
 	    		return 0;
 	   		}
		}
	}

	void export_N_n(char** R_vcf_name,char** R_N_tot_file,char** R_n_alt_file)
	{
	
		//Rprintf("%s:%s:%s\n",*R_vcf_name,*R_N_tot_file,*R_n_alt_file);

		char vcf_name[100];
 		char N_tot_file[100];
 		char n_alt_file[100];

		FILE* fp;
		FILE* ft1,* ft2;
		char vcf_line[MAX_LINE_LEN], header_line[MAX_LINE_LEN];
		char* pch,*pch1,*pch2;
		int cnt1,s_cnt;
		int AD_fld;
		int GT_fld;
		char fmt_fld[MAX_FLD_LEN];
		int i,j;
		int e_flag;
		int N_tot,n_alt,N_minus_n;
		char**N_n_cnt_str;
		int *N_tot_samples,*n_alt_samples;


		strcpy(vcf_name,*R_vcf_name);
		strcpy(N_tot_file,*R_N_tot_file);
		strcpy(n_alt_file,*R_n_alt_file);

		fp = fopen(vcf_name,"r");		

		ft1=fopen(N_tot_file,"w");
 		ft2=fopen(n_alt_file,"w");


		if((NULL==fp)||(NULL==ft1)||(NULL==ft2))
		{
			error("file open error !\n");
			return;
		}
		else
		{
			while (readNextLine (vcf_line, fp) == 0)
			{
 				/* If it is a comment, keep track as header line */
 	     		if (vcf_line[0] == '#')
        			strcpy (header_line, vcf_line);
      			/* If not a comment, break out of this loop */
      			else
        			break;
    		}
			// we can now find number of samples
			cnt1=0;
			pch = strtok (header_line,"\t");
  			while (pch != NULL)
  			{
    			//Rprintf("%s\n",pch);
				pch = strtok (NULL, "\t");
				cnt1++;
  			}
			s_cnt = cnt1-SAMPLE_OFFSET;

			N_n_cnt_str=(char**)(calloc(s_cnt,sizeof(char*)));
			for(i=0;i<s_cnt;i++)
				N_n_cnt_str[i] = (char*)(calloc(1,MAX_FLD_LEN*sizeof(char)));

			N_tot_samples = (int*)(calloc(s_cnt,sizeof(int)));
			n_alt_samples = (int*)(calloc(s_cnt,sizeof(int)));

			cnt1=0;
			j=0;
			pch = strtok (vcf_line,"\t");
  			while (pch != NULL)
  			{
    			//Rprintf("%s\n",pch);
				if(cnt1==FORMAT_INDX)
				{
					strcpy(fmt_fld,pch);
				}
				else
				{
					if(cnt1>FORMAT_INDX)
					{
						strcpy(N_n_cnt_str[j],pch);
						j++;
					}
				}
				pch = strtok (NULL, "\t");
				cnt1++;
			}
			cnt1=0;
			pch = strtok (fmt_fld,":");
  			while (pch != NULL)
  			{
    			//Rprintf("%s\n",pch);
				if(strcmp("GT",pch)==0)
					GT_fld=cnt1;
				if(strcmp("AD",pch)==0)
					AD_fld=cnt1;

				pch = strtok (NULL, ":");
				cnt1++;
			}
			e_flag=0;
  			for(i=0;i<s_cnt;i++)
			{
				if(strcmp(N_n_cnt_str[i],"./.")==0 || strcmp(N_n_cnt_str[i],".|.")==0 )
				{
					e_flag=1;
					break;
				}
				else
				{
					j=0;
					pch1 = strtok(N_n_cnt_str[i],":");
 				    while (pch1 != NULL)
 				    {

 				    	if((j==GT_fld) && ((strcmp(pch1,"./.")==0) || (strcmp(pch1,".|.")==0)) )
 				    	{
 				    		e_flag = 1;
 				    		break;
 				    	}

						if(j==AD_fld)
						{
							sscanf(pch1,"%d,%d",&N_minus_n,&n_alt);
							N_tot_samples[i] = N_minus_n+n_alt;
							n_alt_samples[i] = n_alt;
							if(N_tot_samples[i] == 0)
							{
								e_flag = 1;
								break;
							}

						}
						pch1 = strtok (NULL, ":");
						j++;
					}
				}
				if(e_flag==1)
					break;
			} // end for

			if(e_flag==0)
			{
				for(i=0;i<s_cnt-1;i++)
				{
					fprintf(ft1,"%d\t",N_tot_samples[i]);
					fprintf(ft2,"%d\t",n_alt_samples[i]);
				}
				fprintf(ft1,"%d\n",N_tot_samples[i]);
				fprintf(ft2,"%d\n",n_alt_samples[i]);
			}
		
			while(readNextLine(vcf_line, fp) == 0)
			{
				cnt1=0;
				j=0;
				pch = strtok (vcf_line,"\t");
				while (pch != NULL)
				{
					//Rprintf("%s\n",pch);
					if(cnt1==FORMAT_INDX)
					{
						strcpy(fmt_fld,pch);
					}
					else
					{
						if(cnt1>FORMAT_INDX)
						{
							strcpy(N_n_cnt_str[j],pch);
							j++;
						}
					}
					pch = strtok (NULL, "\t");
						cnt1++;
				}
				cnt1=0;
				pch = strtok (fmt_fld,":");
				while (pch != NULL)
				{
					//Rprintf("%s\n",pch);
					if(strcmp("GT",pch)==0)
						GT_fld=cnt1;
					if(strcmp("AD",pch)==0)
						AD_fld=cnt1;

					pch = strtok (NULL, ":");
					cnt1++;
				}
				e_flag=0;
				for(i=0;i<s_cnt;i++)
				{
					if(strcmp(N_n_cnt_str[i],"./.")==0 || strcmp(N_n_cnt_str[i],".|.")==0 )
					{
						e_flag=1;
						break;
					}
					else
					{
						j=0;
						pch1 = strtok(N_n_cnt_str[i],":");
						while (pch1 != NULL)
						{
							if((j==GT_fld) && ((strcmp(pch1,"./.")==0) || (strcmp(pch1,".|.")==0)) )
							{
								e_flag = 1;
						 		break;
							}
							if(j==AD_fld)
							{
								sscanf(pch1,"%d,%d",&N_minus_n,&n_alt);
								N_tot_samples[i] = N_minus_n+n_alt;
								n_alt_samples[i] = n_alt;
								if(N_tot_samples[i] == 0)
								{
									e_flag=1;
									break;
								}
							}
							pch1 = strtok (NULL, ":");
							j++;
						}
						if(e_flag==1)
							break;
					}
				} // end for

				if(e_flag==0)
				{
					for(i=0;i<s_cnt-1;i++)
					{
						fprintf(ft1,"%d\t",N_tot_samples[i]);
						fprintf(ft2,"%d\t",n_alt_samples[i]);
					}
					fprintf(ft1,"%d\n",N_tot_samples[i]);
					fprintf(ft2,"%d\n",n_alt_samples[i]);
				}
			}
		
		}// else file open success

		for(i=0;i<s_cnt;i++)
			free(N_n_cnt_str[i]);
		free(N_n_cnt_str);
		free(N_tot_samples);
		free(n_alt_samples);


		fclose(fp);
		fclose(ft1);
		fclose(ft2);

		Rprintf("files are generated!\n");

	}

    //UPDATE W, M AND P WITH NEW THETA(W_STAR) FOR ONE TISSUE(t)
	void fn_theta_to_w_M_p_1(double *theta, double p0_z, double *L, double *Z, 
                             int t, int T, int S, int C, double *w, 
                             double *M, double *p)
	{
		int i_c, i_s;
		
		double p_tmp, M_tmp, theta_sum = 0.0;
		
		for (i_c=0; i_c < C; i_c++) {
			theta_sum = theta_sum + theta[i_c*T + t];
		}
		
		for (i_c=0; i_c < C; i_c++) {
			w[i_c] = theta[i_c*T + t]/theta_sum;
		}
		
		for (i_s=0; i_s < S; i_s++) {
			M_tmp = w[0]*L[i_s];

			for (i_c=1; i_c < C; i_c++) {
				M_tmp = M_tmp + w[i_c]*L[i_c*S + i_s];
			}
			M[i_s] = M_tmp;

			p_tmp = p0_z*w[0]*Z[i_s];
            
			for (i_c=1; i_c < C; i_c++) {
				p_tmp = p_tmp + w[i_c]*Z[i_c*S + i_s];
			}
			p[i_s] = p_tmp/M[i_s];
		}
		
		return;
	}

	    
    //CONVERT THETA TO W FOR ALL TISSUES
	void fn_theta_to_w(double *theta, int T, int C, double *w)
	{
		int i_c, i_t;
		double *theta_sum;
		theta_sum = new double[T];
		
        
		for (i_t=0; i_t < T; i_t++) {
			theta_sum[i_t]=0.0;
			for (i_c=0; i_c < C; i_c++) {
				theta_sum[i_t] = theta_sum[i_t] + theta[i_c*T + i_t];
			}
		}
		
		for (i_t=0; i_t < T; i_t++) {
			for (i_c=0; i_c < C; i_c++) {
				w[i_c*T + i_t] = theta[i_c*T + i_t]/theta_sum[i_t];
			}
		}
		
		delete[] theta_sum; theta_sum = NULL;
		return;
	}
    
    //UPDATE M AND P WITH NEW W AND PO_L FOR ALL TISSUES
	void fn_w_to_M_p(double p0_z, double *L, double *Z, int T, int S, int C, 
                     double *w, double *M, double *p)
	{
		//int S = *SS;
        //int C = *CC;
        //int T = *TT;
		
		int i_c, i_s, i_t;
		double M_tmp, p_tmp;
		
		for (i_t=0; i_t < T; i_t++) {
			for (i_s=0; i_s < S; i_s++) {
                //BACKGROUND
				M_tmp = w[i_t]*L[i_s];
                
                //FOR Z
				for (i_c=1; i_c < C; i_c++) {
					M_tmp = M_tmp + w[i_c*T + i_t]*L[i_c*S + i_s];
				}
				M[i_t*S + i_s] = M_tmp;

                //BACKGROUND
				p_tmp = p0_z*w[i_t]*Z[i_s];
                
                //FOR Z
				for (i_c=1; i_c < C; i_c++) {
					p_tmp = p_tmp + w[i_c*T + i_t]*Z[i_c*S + i_s];
				}
				p[i_t*S + i_s] = p_tmp/M[i_t*S + i_s];
			}
		}
		
		return;
	}
    

    //UPDATE P WITH NEW W AND P0_Z FOR ALL TISSUES
	void fn_w_to_p(double p0_z, double *Z, int T, int S, int C, double *w, 
                   double *M, double *p)
	{
		//int S = *SS;
        //int C = *CC;
        //int T = *TT;
		
		int i_c, i_s, i_t;
		double p_tmp;
		
		for (i_t=0; i_t < T; i_t++) {
			for (i_s=0; i_s < S; i_s++) {
                //BACKGROUND
				p_tmp = p0_z*w[i_t]*Z[i_s];
                
                //FOR Z
				for (i_c=1; i_c < C; i_c++) {
					p_tmp = p_tmp + w[i_c*T + i_t]*Z[i_c*S + i_s];
				}
				p[i_t*S + i_s] = p_tmp/M[i_t*S + i_s];
			}
		}
		
		return;
	}


	//UPDATE ONE THETA_TC AND SO RE-COMPUTE M_ST AND P_ST FOR ALL S
    void fn_update_one_theta(double *theta_cur, double *theta_pro, int T, int S, 
                             int C, int i_t, int i_c, double a, double phi_t, 
                             double p0_z, double *L, double *Z, double *w_cur, 
                             double *M_cur, double *p_cur, double *N, double *n, 
                             double *b)
    {
        double theta_tmp, alpha_pro, alpha_cur, u;
        
        //W, M AND P FOR ONE TISSUE SAMPLE
		double *w_pro; w_pro = new double[C];
		double *M_pro; M_pro = new double[S];
		double *p_pro; p_pro = new double[S];

        int ii_c, i_s;
        
        GetRNGstate();
        theta_tmp = exp(log(theta_cur[T*i_c + i_t]) + rnorm(0.0, 0.5));
        PutRNGstate();
        
        theta_pro[T*i_c + i_t] = theta_tmp;
        
        //printf("\n theta_tmp=%f", theta_tmp);
        //printf("\n theta_cur=%f, theta_pro=%f", theta_cur[T*i_c + i_t], 
        //       theta_pro[T*i_c + i_t]);
        
        //PRIOR (Gamma(a, 1)) + JACOBIAN
        alpha_cur = a*log(theta_cur[i_c*T + i_t]) - theta_cur[i_c*T + i_t]; 
        //PRIOR (Gamma(a, 1)) + JACOBIAN
        alpha_pro = a*log(theta_pro[i_c*T + i_t]) - theta_pro[i_c*T + i_t]; 
        
        //printf("\n alpha_cur=%f, alpha_pro=%f", alpha_cur, alpha_pro);
        //FOR PROPOSED THETA, COMPUTE W, M, P
        
        //THETA --> W, M AND P
        fn_theta_to_w_M_p_1(theta_pro, p0_z, L, Z, i_t, T, S, C, w_pro, 
                            M_pro, p_pro); 
        
        for (i_s=0; i_s < S; i_s++) {
            //LIKELIHOOD--POISSON
            alpha_cur = alpha_cur + N[i_t*S + i_s]*
                                    log(M_cur[i_t*S + i_s])-phi_t*
                                    M_cur[i_t*S + i_s]*b[i_t*S + i_s]/2.0;  
            //LIKELIHOOD--POISSON
            alpha_pro = alpha_pro + N[i_t*S + i_s]*log(M_pro[i_s]) - 
                                    phi_t*M_pro[i_s]*b[i_t*S + i_s]/2.0;  

            //LIKELIHOOD--BINOMIAL
            alpha_cur = alpha_cur + n[i_t*S + i_s]*log(p_cur[i_t*S + i_s]) + 
                                    (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                    log(1.0-p_cur[i_t*S + i_s]); 
            //LIKELIHOOD--BINOMIAL                        
            alpha_pro = alpha_pro + n[i_t*S + i_s]*log(p_pro[i_s]) + 
                                    (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                    log(1.0-p_pro[i_s]);  
        }
        
        //printf("\n alpha_cur=%f, alpha_pro=%f", alpha_cur, alpha_pro);
        
        //ACCEPT OR NOT
        GetRNGstate();
        u = runif(0.0, 1.0);
        PutRNGstate();
        
        //printf("\n log(u)=%f, (alpha_cur-alpha_pro)=%f", log(u), 
        //      alpha_pro-alpha_cur);
        
        if (log(u) < (alpha_pro-alpha_cur)) {
            theta_cur[T*i_c + i_t] = theta_tmp;
            
            for (ii_c = 0; ii_c < C; ii_c++) {
                w_cur[T*ii_c + i_t] = w_pro[ii_c];
            }
            
            for (i_s=0; i_s < S; i_s++) {
                M_cur[i_t*S + i_s] = M_pro[i_s];
                p_cur[i_t*S + i_s] = p_pro[i_s];
            }
        }else {
            theta_pro[T*i_c + i_t] = theta_cur[T*i_c + i_t];
        }
        
		delete[] w_pro; w_pro = NULL;
		delete[] M_pro; M_pro = NULL;
		delete[] p_pro; p_pro = NULL;
    }
    

    
  
    //W_STAR=THETA
	void fn_update_theta_p0(double d0, double d, double a0_z, double b0_z, 
                            double *Z, double *L, double *phi, int T, int S, 
                            int C, double *n, double *N, double *b, 
                            double *p0_z, double *theta_cur, double *w_cur, 
                            double *p_cur, double *M_cur)
	{
        double p0_z_cur, p0_z_pro;
        p0_z_cur = *p0_z;

        double *theta_pro; theta_pro = new double[C*T];
		double *p_pro; p_pro = new double[S*T];
        
        
		double alpha_cur, alpha_pro, u;
		
		int i_s, i_t, i_c, i;
		
		//COPY CURRENT VALUES TO THE PROPOSAL
		for (i = 0; i < (C*T); i++) {
			theta_pro[i] = theta_cur[i];
		}
		
        
        //UPDATE THETA (W_STAR) -- RECOMPUTE W, M AND P.
		for (i_t=0; i_t < T; i_t++) {
            i_c = 0;
            fn_update_one_theta(theta_cur, theta_pro, T, S, C, i_t, i_c, d0, 
                                phi[i_t], p0_z_cur, L, Z, w_cur, M_cur, 
                                p_cur, N, n, b);

			for (i_c=1; i_c < C; i_c++) {
                //printf("\n t=%d, c=%d", i_t, i_c);
                fn_update_one_theta(theta_cur, theta_pro, T, S, C, i_t, i_c, d, 
                                    phi[i_t], p0_z_cur, L, Z, w_cur, M_cur, 
                                    p_cur, N, n, b);
			}
		}
        
  		delete[] theta_pro; theta_pro = NULL;
        
        ///////////////////////////////////////////////////////////////////////
        //UPDATE P0_Z//////////////////////////////////////////////////////////
        GetRNGstate();
        p0_z_pro = exp(log(p0_z_cur) + rnorm(0.0, 0.2));
        PutRNGstate();
        
        //printf("\n q0_cur=%f, q0_pro=%f", p0_z_cur, p0_z_pro);

        //UPDATE P WITH P0 AND P0_CUR
        fn_w_to_p(p0_z_pro, Z, T, S, C, w_cur, M_cur, p_pro);

        //PRIOR (Beta(a0, b0)) + JACOBIAN
        alpha_cur = a0_z*log(p0_z_cur) + (b0_z-1.0)*log(1.0-p0_z_cur); 
        //PRIOR (Beta(a0, b0)) + JACOBIAN
        alpha_pro = a0_z*log(p0_z_pro) + (b0_z-1.0)*log(1.0-p0_z_pro); 
        
        //printf("\n alpha_cur=%f, alpha_pro=%f", alpha_cur, alpha_pro);
        
        for (i_s=0; i_s < S; i_s++) {
            for (i_t=0; i_t < T; i_t++) {
                //LIKELIHOOD --BINOMIAL
                alpha_cur = alpha_cur + n[i_t*S + i_s]*log(p_cur[i_t*S + i_s]) + 
                                        (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                        log(1.0-p_cur[i_t*S + i_s]);  
                //LIKELIHOOD --BINOMIAL                        
                alpha_pro = alpha_pro + n[i_t*S + i_s]*log(p_pro[i_t*S + i_s]) + 
                                        (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                        log(1.0-p_pro[i_t*S + i_s]);  
            }
        }
        
        //ACCEPT OR NOT
        GetRNGstate();
        u = runif(0.0, 1.0);
        PutRNGstate();
        
        //printf("\n log(u)=%f, (alpha_cur-alpha_pro)=%f", log(u), 
        //       alpha_pro-alpha_cur);
        
        if (log(u) < (alpha_pro-alpha_cur)) {
            *p0_z = p0_z_pro;
            for(i=0; i < S*T; i++)p_cur[i] = p_pro[i];
        }
        
		delete[] p_pro; p_pro = NULL;
	}
    
    //UPDATE ONE THETA_TC AND SO RE-COMPUTE M_ST AND P_ST FOR ALL S
    void fn_update_one_sample_theta(double *theta_cur, double *theta_pro, int T, 
                                    int S, int C, int i_t, double a, 
                                    double phi_t, double p0_z, double *L, 
                                    double *Z, double *w_cur, double *M_cur, 
                                    double *p_cur, double *N, double *n, 
                                    double *b)
    {
        double theta_tmp, alpha_pro, alpha_cur, u;
        
        //W, M AND P FOR ONE TISSUE SAMPLE
		double *w_pro; w_pro = new double[C];
		double *M_pro; M_pro = new double[S];
		double *p_pro; p_pro = new double[S];
        
        int i_c, i_s;
        
        alpha_cur = 0.0;
        alpha_pro = 0.0;
        
        for (i_c = 0; i_c < C; i_c++) {
            GetRNGstate();
            theta_tmp = exp(log(theta_cur[T*i_c + i_t]) + rnorm(0.0, 0.2));
            PutRNGstate();
            
            theta_pro[T*i_c + i_t] = theta_tmp;
            
            //printf("\n theta_tmp=%f", theta_tmp);
            //printf("\n theta_cur=%f, theta_pro=%f", theta_cur[T*i_c + i_t], 
            //       theta_pro[T*i_c + i_t]);
            
            //PRIOR (Gamma(a, 1)) + JACOBIAN
            alpha_cur = alpha_cur + a*log(theta_cur[i_c*T + i_t]) - 
                                    theta_cur[i_c*T + i_t]; 
            //PRIOR (Gamma(a, 1)) + JACOBIAN
            alpha_pro = alpha_pro + a*log(theta_pro[i_c*T + i_t]) - 
                                    theta_pro[i_c*T + i_t]; 
            
        }
        
        //printf("\n alpha_cur=%f, alpha_pro=%f", alpha_cur, alpha_pro);
        //FOR PROPOSED THETA, COMPUTE W, M, P
        fn_theta_to_w_M_p_1(theta_pro, p0_z, L, Z, i_t, T, S, C, w_pro, 
                            M_pro, p_pro); //THETA --> W, M AND P
        
        for (i_s=0; i_s < S; i_s++) {
            //LIKELIHOOD--POISSON
            alpha_cur = alpha_cur + N[i_t*S + i_s]*log(M_cur[i_t*S + i_s]) - 
                                    phi_t*M_cur[i_t*S + i_s]*
                                    b[i_t*S + i_s]/2.0;  
            //LIKELIHOOD--POISSON                        
            alpha_pro = alpha_pro + N[i_t*S + i_s]*log(M_pro[i_s]) - 
                                    phi_t*M_pro[i_s]*b[i_t*S + i_s]/2.0;  
            
            //LIKELIHOOD--BINOMIAL
            alpha_cur = alpha_cur + n[i_t*S + i_s]*log(p_cur[i_t*S + i_s]) + 
                                    (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                    log(1.0-p_cur[i_t*S + i_s]);
            //LIKELIHOOD--BINOMIAL
            alpha_pro = alpha_pro + n[i_t*S + i_s]*log(p_pro[i_s]) + 
                                    (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                    log(1.0-p_pro[i_s]);  
        }
        
        //printf("\n alpha_cur=%f, alpha_pro=%f", alpha_cur, alpha_pro);
        
        //ACCEPT OR NOT
        GetRNGstate();
        u = runif(0.0, 1.0);
        PutRNGstate();
        
        //printf("\n log(u)=%f, (alpha_cur-alpha_pro)=%f", log(u), 
        //       alpha_pro-alpha_cur);
        
        if (log(u) < (alpha_pro-alpha_cur)) {
            for (i_c = 0; i_c < C; i_c ++) {
                theta_cur[T*i_c + i_t] = theta_pro[T*i_c + i_t];
                w_cur[T*i_c + i_t] = w_pro[i_c];
            }

            for (i_s=0; i_s < S; i_s++) {
                M_cur[i_t*S + i_s] = M_pro[i_s];
                p_cur[i_t*S + i_s] = p_pro[i_s];
            }
            
            //printf("\n multi sample w accept");
        }else {
            for (i_c = 0; i_c < C; i_c ++) {
                theta_pro[T*i_c + i_t] = theta_cur[T*i_c + i_t];
            }
        }
        
		delete[] w_pro; w_pro = NULL;
		delete[] M_pro; M_pro = NULL;
		delete[] p_pro; p_pro = NULL;
    }
    
    //W_STAR=THETA
	void fn_update_multi_theta_p0(double d0, double d, double a0_z, double b0_z, 
                                  double *Z, double *L, double *phi,  int T, 
                                  int S, int C, double *n, double *N, double *b, 
                                  double *p0_z, double *theta_cur, 
                                  double *w_cur, double *p_cur, double *M_cur)
	{
        double p0_z_cur, p0_z_pro;
        p0_z_cur = *p0_z;
        
        double *theta_pro; theta_pro = new double[C*T];
		double *p_pro; p_pro = new double[S*T];
        
        
		double alpha_cur, alpha_pro, u;
		
		int i_s, i_t, i;
		
		//COPY CURRENT VALUES TO THE PROPOSAL
		for (i = 0; i < (C*T); i++) {
			theta_pro[i] = theta_cur[i];
		}
		
        
        //UPDATE THETA (W_STAR) -- RECOMPUTE W, M AND P.
		for (i_t=0; i_t < T; i_t++) {
            fn_update_one_sample_theta(theta_cur, theta_pro, T, S, C, i_t, d0, 
                                        phi[i_t], p0_z_cur, L, Z, w_cur, M_cur, 
                                        p_cur, N, n, b);
		}
        
  		delete[] theta_pro; theta_pro = NULL;
        
        //////////////////////////////////////////////////////////////////////
        //UPDATE P0_Z/////////////////////////////////////////////////////////
        GetRNGstate();
        p0_z_pro = exp(log(p0_z_cur) + rnorm(0.0, 0.1));
        PutRNGstate();
        
        //printf("\n q0_cur=%f, q0_pro=%f", p0_z_cur, p0_z_pro);
        
        //UPDATE P WITH P0 AND P0_CUR
        fn_w_to_p(p0_z_pro, Z, T, S, C, w_cur, M_cur, p_pro);
        
        //PRIOR (Beta(a0, b0)) + JACOBIAN
        alpha_cur = a0_z*log(p0_z_cur) + (b0_z-1.0)*log(1.0-p0_z_cur);
        //PRIOR (Beta(a0, b0)) + JACOBIAN
        alpha_pro = a0_z*log(p0_z_pro) + (b0_z-1.0)*log(1.0-p0_z_pro);
        
        //printf("\n alpha_cur=%f, alpha_pro=%f", alpha_cur, alpha_pro);
        
        for (i_s=0; i_s < S; i_s++) {
            for (i_t=0; i_t < T; i_t++) {
                //LIKELIHOOD --BINOMIAL
                alpha_cur = alpha_cur + n[i_t*S + i_s]*log(p_cur[i_t*S + i_s])+ 
                                        (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                        log(1.0-p_cur[i_t*S + i_s]);
                //LIKELIHOOD --BINOMIAL
                alpha_pro = alpha_pro + n[i_t*S + i_s]*log(p_pro[i_t*S + i_s])+ 
                                        (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                        log(1.0-p_pro[i_t*S + i_s]);  
            }
        }
        
        //ACCEPT OR NOT
        GetRNGstate();
        u = runif(0.0, 1.0);
        PutRNGstate();
        
        //printf("\n log(u)=%f, (alpha_cur-alpha_pro)=%f", log(u), 
        //       alpha_pro-alpha_cur);
        
        if (log(u) < (alpha_pro-alpha_cur)) {
            *p0_z = p0_z_pro;
            for(i=0; i < S*T; i++)p_cur[i] = p_pro[i];
        }
        
		delete[] p_pro; p_pro = NULL;
	}


    //UPDATE PI
    void fn_update_pi(double alpha, double beta, double *gam, int Q, int C, 
                      int S, double *A, double *pi_mat)
    {
        double *pi_star; pi_star = new double[Q];
        double pi_star_sum;
        int i_q;
        
        //CELL TYPE 0
        int i_c = 0;
        for (i_q=0; i_q <= Q; i_q++) {
            pi_mat[C*i_q + i_c] = 0.0;
        }
        
        pi_mat[C*2 + i_c] = 1.0;
        
        //NO CELL TYPE 0 FOR Q=2 (NO CNV)
        i_q = 2;
        for (i_c=1; i_c < C; i_c++) {
            GetRNGstate();
            pi_mat[C*i_q + i_c] = rbeta((beta + A[C*i_q + i_c]), 
                                        (alpha + (1.0*S - A[C*i_q + i_c])));
            PutRNGstate();
        }
        
        //NO CELL TYPE 0 FOR Q!=2 (CNV)
        for(i_c=1; i_c < C ; i_c++)
        {
            //PI_STAR
            pi_star_sum = 0.0;
            for (i_q=0; i_q < 2; i_q++) {
                GetRNGstate();
                pi_star[i_q] = rgamma(gam[i_q] + A[C*i_q + i_c], 1.0);
                PutRNGstate();
                pi_star_sum = pi_star_sum + pi_star[i_q];
                
                //printf("\n i_c=%d, i_q=%d, gam[%d]=%f, A[%d, %d]=%f, 
                //       pi_star=%f", i_c, i_q, i_q, gam[i_q], i_c, i_q, 
                //       A[C*i_q + i_c], pi_star[i_q]);
            }
            
            for (i_q=3; i_q < (Q+1); i_q++) {
                GetRNGstate();
                pi_star[i_q-1] = rgamma(gam[i_q-1] + A[C*i_q + i_c], 1.0);
                PutRNGstate();
                pi_star_sum = pi_star_sum + pi_star[i_q-1];
                //printf("\n i_c=%d, i_q=%d, gam[%d]=%f, A[%d, %d]=%f, 
                //       pi_star=%f", i_c, i_q, i_q, gam[i_q-1], i_c, i_q, 
                //       A[C*i_q + i_c], pi_star[i_q-1]);
            }

           // printf("\n pi_2=%f, pi_sum=%f", pi_mat[C*2 + i_c], pi_star_sum);
            //PI_ST
            i_q = 0;
            pi_mat[C*i_q + i_c] = pi_star[i_q]/pi_star_sum*
                                  (1.0-pi_mat[C*2 + i_c]);
            
            //printf("\n i_c=%d, i_q=%d, pi=%f", i_c, i_q, pi_mat[C*i_q + i_c]);

            i_q = 1;
            pi_mat[C*i_q + i_c] = pi_star[i_q]/pi_star_sum*
                                  (1.0-pi_mat[C*2 + i_c]);
            
            //printf("\n i_c=%d, i_q=%d, pi=%f", i_c, i_q, pi_mat[C*i_q + i_c]);
            for (i_q=3; i_q <= Q; i_q++) {
                pi_mat[C*i_q + i_c] = pi_star[i_q-1]/pi_star_sum*
                                      (1.0-pi_mat[C*2 + i_c]);
                //printf("\n i_c=%d, i_q=%d, pi=%f", i_c, i_q, 
                //       pi_mat[C*i_q + i_c]);
            }
        } //for(i.c in 1:CC)
        
        delete[] pi_star; pi_star = NULL;

    }

    
    
    //UPDATE PHI
    void fn_update_phi(int S, int T, double a, double b, double *M, 
                       double *N_sum, double *b_st, double *phi)
    {
        int i_t, i_s;
        double M_sum;
        

        for (i_t=0; i_t < T; i_t++) {
            M_sum = 0.0;
            
            for (i_s=0; i_s < S; i_s++) {
                M_sum = M_sum + M[i_t*S + i_s]*b_st[i_t*S + i_s];
            }

            GetRNGstate();
            //rgamma(a,b) for mean=ab
            phi[i_t] = rgamma(a+N_sum[i_t], 1.0/(b + M_sum/2.0));  
            PutRNGstate();
        }
    }
    
	double mymax(double a, double b)
	{
		return (a > b) ? a : b;
	}
	
	double mymin(double a, double b)
	{
		return (a < b) ? a : b;
	}

    
 
    
    //UPDATE L -- NEED TO RECOMPUTE A, M AND P, AS WELL
    void fn_update_multi_L(int C, int S, int T, int Q, double *n, double *N, 
                           double *b, double *ppi, double *phi, double *w, 
                           double *Z, double *A, double p_z0, double *M_cur, 
                           double *p_cur, double *L_cur)
    {
        double *w_Z; w_Z = new double[T*C];

        double *L_s_cur; L_s_cur = new double[C];
        double *L_s_pro; L_s_pro = new double[C];

        double *M_s_cur; M_s_cur = new double[T];
        double *M_s_pro; M_s_pro = new double[T];
        
        double *p_s_cur; p_s_cur = new double[T];
        double *p_s_pro; p_s_pro = new double[T];
       
        int i_t, i_c, i_s;
        int l_sc_cur, l_sc_pro;
        double Low, Upp;
        
        double M_tmp, p_tmp, u, prob_cur, prob_pro;
        
        i_c = 0;
        for (i_t = 0; i_t < T; i_t ++) {
            w_Z[T*i_c + i_t] = w[T*i_c + i_t]*(p_z0);
        }

        for (i_t = 0; i_t < T; i_t ++) {
            for (i_c = 1; i_c < C; i_c++) {
                w_Z[T*i_c + i_t] = w[T*i_c + i_t];
            }
        }

        for (i_s = 0; i_s < S; i_s++) {
            //COPY ROW S
            for (i_c=0; i_c < C; i_c++) {
                L_s_cur[i_c]=L_cur[i_c*S + i_s];
            }

            for (i_t=0; i_t < T; i_t++) {
                M_s_cur[i_t] = M_cur[i_t*S + i_s];
                p_s_cur[i_t] = p_cur[i_t*S + i_s];
            }
            
            
            //MAKE A PROPOSAL OF L_S AND COMPUTE M_ST AND P_ST ACCORDINGLY
            prob_cur = 0.0;
            prob_pro = 0.0;

            L_s_pro[0]=2.0; //CELL TYPE 0
            for (i_c=1; i_c < C; i_c++) {
                Low = mymax(L_s_cur[i_c]-1.0, Z[i_c*S + i_s]);
                Upp = mymin(L_s_cur[i_c]+1.0, Q);
                if (Low==Upp) {
                    L_s_pro[i_c] = Low;
                }else{
                    GetRNGstate();
                    u = runif(0.0, 1.0);
                    PutRNGstate();
                    
                    L_s_pro[i_c] = Low + floor(u*(Upp - Low + 1.0));
                }
                
                //PROB OF PROPOSAL Q(L_OLD -> L_NEW)
                prob_cur = prob_cur - log(Upp - Low + 1.0);  
            } //for (i_c=1; i_c < C; i_c++) {
            
            
            for (i_c=1; i_c < C; i_c++) {
                Low = mymax(L_s_pro[i_c]-1.0, Z[i_c*S + i_s]);
                Upp = mymin(L_s_pro[i_c]+1.0, Q);
                
                //PROB OF PROPOSAL Q(L_NEW -> L_OLD)
                prob_cur = prob_pro - log(Upp - Low + 1.0);   
            } //for (i_c=1; i_c < C; i_c++) {

            //COMPUTE M AND P WITH l_sc=i.q==NEW
            for (i_t=0; i_t < T; i_t++) {
                M_tmp = 0.0;
                for (i_c=0; i_c < C; i_c++) {
                    M_tmp = M_tmp + w[T*i_c + i_t]*L_s_pro[i_c];
                }
                M_s_pro[i_t] = M_tmp;
                
                p_tmp = 0.0;
                for (i_c=0; i_c < C; i_c++) {
                    p_tmp = p_tmp + w_Z[T*i_c + i_t]*Z[i_c*S + i_s];
                }
                p_s_pro[i_t] = p_tmp/M_tmp;
            }

            
            //COMPUTE PROBABILITY FOR CUR AND PRO
            for (i_c=0; i_c < C; i_c++) {
                l_sc_cur = L_s_cur[i_c];
                prob_cur = prob_cur + log(ppi[C*l_sc_cur + i_c]) - 
                                      log(l_sc_cur + 1.0);

                l_sc_pro = L_s_pro[i_c];
                prob_pro = prob_pro + log(ppi[C*l_sc_pro + i_c]) - 
                                      log(l_sc_pro + 1.0);
            }//for (i_c=0; i_c < C; i_c++) {
            
            for (i_t = 0; i_t < T; i_t ++) {
                //LIKELIHOOD--POISSON
                prob_cur = prob_cur + N[i_t*S + i_s]*log(M_s_cur[i_t]) - 
                                      phi[i_t]*M_s_cur[i_t]*b[i_t*S + i_s]/2.0; 
                //LIKELIHOOD--BINOMIAL                      
                prob_cur = prob_cur + n[i_t*S + i_s]*log(p_s_cur[i_t]) + 
                                      (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                      log(1.0-p_s_cur[i_t]); 

                //LIKELIHOOD--POISSON
                prob_pro = prob_pro + N[i_t*S + i_s]*log(M_s_pro[i_t]) - 
                                      phi[i_t]*M_s_pro[i_t]*b[i_t*S + i_s]/2.0;
                //LIKELIHOOD--BINOMIAL
                prob_pro = prob_pro + n[i_t*S + i_s]*log(p_s_pro[i_t]) + 
                                      (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                      log(1.0-p_s_pro[i_t]);
            }//for (i_t = 0; i_t < T; i_t ++) {
            
            //ACCEPT??
            GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
            
            if(log(u) < (prob_pro - prob_cur))
            {
                for (i_c=0; i_c < C; i_c++) {
                    l_sc_cur = L_s_cur[i_c];
                    A[C*l_sc_cur + i_c] = A[C*l_sc_cur + i_c] - 1.0;
                    
                    l_sc_pro = L_s_pro[i_c];
                    A[C*l_sc_pro + i_c] = A[C*l_sc_pro + i_c] + 1.0;

                    L_cur[i_c*S + i_s] = L_s_pro[i_c];
                }
                
                for (i_t=0; i_t < T; i_t++) {
                    M_cur[i_t*S + i_s] = M_s_pro[i_t];
                    p_cur[i_t*S + i_s] = p_s_pro[i_t];
                }
            }
        }//for (i_s = 0; i_s < S; i_s++) {
        
        
        delete[] w_Z; w_Z = NULL;
        delete[] L_s_cur; L_s_cur = NULL;
        delete[] L_s_pro; L_s_pro = NULL;
        delete[] M_s_cur; M_s_cur = NULL;
        delete[] M_s_pro; M_s_pro = NULL;
        delete[] p_s_cur; p_s_cur = NULL;
        delete[] p_s_pro; p_s_pro = NULL;

    }//void fn_update_multi_L
    


    
    //UPDATE Z -- NEED TO RECOMPUTE P
    void fn_update_multi_Z(int C, int S, int T, double *n, double *N, double *w, 
                           double *Z_cur, double p_z0, double *M, double *p_cur, 
                           double *L)
    {
        double *w_Z; w_Z = new double[T*C];
        
        double *Z_s_cur; Z_s_cur = new double[C];
        double *Z_s_pro; Z_s_pro = new double[C];
        
        double *p_s_cur; p_s_cur = new double[T];
        double *p_s_pro; p_s_pro = new double[T];
        
        
        int i_t, i_c, i_s;
        double Low, Upp;
        
        double p_tmp, u, prob_cur, prob_pro;
        
        i_c = 0;
        for (i_t = 0; i_t < T; i_t ++) {
            w_Z[T*i_c + i_t] = w[T*i_c + i_t]*(p_z0);
        }
        
        for (i_t = 0; i_t < T; i_t ++) {
            for (i_c = 1; i_c < C; i_c++) {
                w_Z[T*i_c + i_t] = w[T*i_c + i_t];
            }
        }
        
        for (i_s = 0; i_s < S; i_s++) {
            //COPY ROW S
            for (i_c=0; i_c < C; i_c++) {
                Z_s_cur[i_c]=Z_cur[i_c*S + i_s];
            }
            
            for (i_t=0; i_t < T; i_t++) {
                p_s_cur[i_t] = p_cur[i_t*S + i_s];
            }
            
            
            //MAKE A PROPOSAL OF L_S AND COMPUTE M_ST AND P_ST ACCORDINGLY
            prob_cur = 0.0;
            prob_pro = 0.0;

            Z_s_pro[0]=2.0; //CELL TYPE 0
            for (i_c=1; i_c < C; i_c++) {
                Low = mymax(Z_s_cur[i_c]-1.0, 0.0);
                Upp = mymin(L[i_c*S + i_s], Z_s_cur[i_c]+1.0);
                
                if (Low==Upp) {
                    Z_s_pro[i_c] = Low;
                }else{
                    GetRNGstate();
                    u = runif(0.0, 1.0);
                    PutRNGstate();
                    
                    Z_s_pro[i_c] = Low + floor(u*(Upp - Low + 1.0));
                }
                
                //PROPOSAL DIST Q(Z_OLD - > Z_NEW)
                prob_cur = prob_cur - log(Upp - Low + 1.0);  
            }//for (i_c=1; i_c < C; i_c++) {
            
            
            for (i_c=1; i_c < C; i_c++) {
                Low = mymax(Z_s_pro[i_c]-1.0, 0.0);
                Upp = mymin(L[i_c*S + i_s], Z_s_pro[i_c]+1.0);
                
                //PROPOSAL DIST Q(Z_NEW - > Z_OLD)
                prob_pro = prob_pro - log(Upp - Low + 1.0);  
            }//for (i_c=1; i_c < C; i_c++) {

            //COMPUTE M AND P WITH l_sc=i.q
            for (i_t=0; i_t < T; i_t++) {
                p_tmp = 0.0;
                for (i_c=0; i_c < C; i_c++) {
                    p_tmp = p_tmp + w_Z[T*i_c + i_t]*Z_s_pro[i_c];
                }
                p_s_pro[i_t] = p_tmp/M[i_t*S + i_s];
            }//for (i_t=0; i_t < T; i_t++) {
            
            
            //COMPUTE PROBABILITY FOR CUR AND PRO
            for (i_t = 0; i_t < T; i_t ++) {
                //LIKELIHOOD--BINOMIAL
                prob_cur = prob_cur + n[i_t*S + i_s]*log(p_s_cur[i_t]) + 
                                      (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                      log(1.0-p_s_cur[i_t]); 
                //LIKELIHOOD--BINOMIAL
                prob_pro = prob_pro + n[i_t*S + i_s]*log(p_s_pro[i_t]) + 
                                      (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                      log(1.0-p_s_pro[i_t]); 
            }//for (i_t = 0; i_t < T; i_t ++) {
            
            //ACCEPT??
            GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
            
            if(log(u) < (prob_pro - prob_cur))
            {
                for (i_c=0; i_c < C; i_c++) {
                    Z_cur[i_c*S + i_s] = Z_s_pro[i_c];
                }
                
                for (i_t=0; i_t < T; i_t++) {
                    p_cur[i_t*S + i_s] = p_s_pro[i_t];
                }
            }
        }//for (i_s = 0; i_s < S; i_s++) {
        
        delete[] w_Z; w_Z = NULL;
        delete[] Z_s_cur; Z_s_cur = NULL;
        delete[] Z_s_pro; Z_s_pro = NULL;
        delete[] p_s_cur; p_s_cur = NULL;
        delete[] p_s_pro; p_s_pro = NULL;

    }//void fn_update_multi_Z

    
    
    //UPDATE L-- NEED TO RECOMPUTE A, M AND P, AS WELL
    void fn_update_L(int C, int S, int T, int Q, double *n, double *N, 
                     double *b, double *ppi, double *phi, double *w, double *Z, 
                     double *A, double p_z0, double *M_cur, double *p_cur, 
                     double *L_cur)
    {
        double *w_Z; w_Z = new double[T*C];
        
        double *L_s; L_s = new double[C];
        
        double *M_Q; M_Q = new double[(Q+1)*T];  //Q+1 BY T MATRIX
        double *p_Q; p_Q = new double[(Q+1)*T];  //Q+1 BY T MATRIX

        
        double *prob_dist; prob_dist = new double[Q+1];
        int *l_set; l_set= new int[Q+1];
        
        int i_t, i_c, ii_c, i_s, i_q, i;
        int l_sc_cur, l_sc_new, z_sc;  //ALTHOUGH THESE ARE INTEGERS, FINE.
        int n_cnt;
        
        double M_tmp, p_tmp, u, prob_max, prob_sum, prob_tmp;
        
        i_c = 0;
        for (i_t = 0; i_t < T; i_t ++) {
            w_Z[T*i_c + i_t] = w[T*i_c + i_t]*(p_z0);
        }
        
        for (i_t = 0; i_t < T; i_t ++) {
            for (i_c = 1; i_c < C; i_c++) {
                w_Z[T*i_c + i_t] = w[T*i_c + i_t];
            }
        }
        
        
        for(i_s=0; i_s < S; i_s++)
        {
            //printf("\n\n i_s=%d", i_s);
            for(i_c=1; i_c < C; i_c++)
            {
                //printf("\n\n i_c=%d", i_c);
                //COPY ROW S
                for (ii_c=0; ii_c < C; ii_c++) {
                    L_s[ii_c]=L_cur[ii_c*S + i_s];
                }
                n_cnt = 0;

                l_sc_cur = L_s[i_c];
                z_sc = Z[i_c*S + i_s];
 
                //printf("\n  CURRENT l_sc=%d, z_sc=%d", l_sc_cur, z_sc);
                //l_sc == l_sc_cur
                i_q = l_sc_cur;
                l_set[n_cnt] = i_q;

                
                prob_tmp = log(ppi[C*i_q + i_c]) - log(i_q+1.0);
                //printf("\n prob=%f", prob_tmp);
                for (i_t=0; i_t < T; i_t++) {
                    //LIKELIHOOD--POISSON
                    prob_tmp = prob_tmp + N[i_t*S + i_s]*
                                          log(M_cur[i_t*S + i_s]) - phi[i_t]*
                                          M_cur[i_t*S + i_s]*b[i_t*S + i_s]/2.0;
                    //LIKELIHOOD--BINOMIAL
                    prob_tmp = prob_tmp + n[i_t*S + i_s]*
                                          log(p_cur[i_t*S + i_s]) + 
                                          (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                          log(1.0-p_cur[i_t*S + i_s]); 
                    
                    //printf("\n i_t=%d, M=%f, p=%f, prob=%f", i_t, 
                    //       M_cur[i_t*S + i_s], p_cur[i_t*S + i_s], prob_tmp);
                }
                //printf("\n i_q=%d, n_cnt=%d, l_set=%d, prob=%f", i_q, n_cnt, 
                //       l_set[n_cnt], prob_tmp);
                prob_dist[n_cnt] = prob_tmp;
                prob_max = prob_tmp;
                n_cnt = n_cnt + 1;
                
                for (i_q=z_sc; i_q <= Q; i_q++) {
                    if (!(i_q==l_sc_cur)) {
                        //printf("\n i_q=%d, n_cnt=%d", i_q, n_cnt);
                        l_set[n_cnt] = i_q;
                        L_s[i_c] = 1.0*i_q;
                        
                        //COMPUTE M AND P WITH l_sc=i.q
                        for(i_t=0; i_t < T; i_t++)
                        {
                            M_tmp = 0.0;
                            p_tmp = 0.0;
                            
                            //printf("\n i_t=%d", i_t);
                            for (ii_c=0; ii_c < C; ii_c++) {
                                M_tmp = M_tmp + w[T*ii_c + i_t]*L_s[ii_c];
                                p_tmp = p_tmp + w_Z[T*ii_c + i_t]*
                                                Z[ii_c*S + i_s];
                                //printf("\n ii_c=%d, w_L=%f, L=%f, w_Z=%f, 
                                //       Z=%f, M_tmp=%f, p_tmp=%f", ii_c, 
                                //       w_L[T*ii_c + i_t], L_s[ii_c], 
                                //       w_Z[T*ii_c + i_t], Z[ii_c*S + i_s], 
                                //       M_tmp, p_tmp);
                            }
                            M_Q[(Q+1)*i_t + n_cnt] = M_tmp;//M_ST WITH L_SC=I_Q
                            p_Q[(Q+1)*i_t + n_cnt] = p_tmp/M_tmp;
                        }
                        
                        prob_tmp = log(ppi[C*i_q + i_c]) - log(i_q+1.0);
                        //printf("\n prob=%f", prob_tmp);

                        for (i_t=0; i_t < T; i_t++) {
                            //LIKELIHOOD--POISSON
                            prob_tmp = prob_tmp + N[i_t*S + i_s]*
                                                  log(M_Q[(Q+1)*i_t + n_cnt]) - 
                                                  phi[i_t]*M_Q[(Q+1)*i_t + 
                                                  n_cnt]*b[i_t*S + i_s]/2.0;
                            //LIKELIHOOD--BINOMIAL
                            prob_tmp = prob_tmp+n[i_t*S + i_s]*
                                                log(p_Q[(Q+1)*i_t + n_cnt]) + 
                                                (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                                log(1.0-p_Q[(Q+1)*i_t + n_cnt]); 
                            
                            //printf("\n i_t=%d, M=%f, p=%f, prob=%f", i_t, 
                            //       M_Q[(Q+1)*i_t + n_cnt], 
                            //       p_Q[(Q+1)*i_t + n_cnt], prob_tmp);

                        }
                        prob_dist[n_cnt] = prob_tmp;
                        prob_max = mymax(prob_max,prob_tmp);
                        //printf("\n i_q=%d, n_cnt=%d, l_set=%d, prob=%f", 
                        //       i_q, n_cnt, l_set[n_cnt], prob_tmp);

                        n_cnt = n_cnt + 1;
                    } //if (!(i_q==l_sc_cur)) {
                } //for (i_q=z_sc; i_q <= Q; i_q++)
                
                
                prob_sum = 0.0;
                for (i=0; i < n_cnt; i++) {
                    prob_dist[i] = exp(prob_dist[i] - prob_max);
                    prob_sum = prob_sum + prob_dist[i];
                }
                
                prob_dist[0] = prob_dist[0]/prob_sum;
                //printf("\n prob[%d]=%f", 0, prob_dist[0]);
                for (i=1; i < n_cnt; i++) {
                    prob_dist[i] = prob_dist[i-1] + prob_dist[i]/prob_sum;
                    //printf("\n prob[%d]=%f", i, prob_dist[i]);
                }
                
                
                //SAMPLE L_SC 
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();

                //printf("\n u=%f", u);
                for (i=0; i < n_cnt; i++) {
                    if (u < prob_dist[i]) break;
                }
                
                l_sc_new = l_set[i];

                //printf("\n i_s=%d, i_c=%d, n_cnt=%d, i=%d, z_sc=%d, 
                //       l_sc_cur=%d, l_sc_new=%d", i_s, i_c, n_cnt, 
                //       i, z_sc, l_sc_cur, l_sc_new);
                //if (l_sc_new > Q) {
                //    printf("\n\n LLLLLLLLLLLLLLLLLLLLLOOOOOOOOOOOOOOOOOOOO");
                //}
               
                if(!(l_sc_cur == l_sc_new))
                {
                    //printf("\n A[%d]=%f, A[%d]=%f, ", C*l_sc_cur + i_c, 
                    //       A[C*l_sc_cur + i_c], C*l_sc_new + i_c, 
                    //       A[C*l_sc_new + i_c]);
                    
                    A[C*l_sc_cur + i_c] = A[C*l_sc_cur + i_c] - 1.0;
                    A[C*l_sc_new + i_c] = A[C*l_sc_new + i_c] + 1.0;
                    
                    //printf("\n A[%d]=%f, A[%d]=%f, ", C*l_sc_cur + i_c, 
                    //       A[C*l_sc_cur + i_c], C*l_sc_new + i_c, 
                    //       A[C*l_sc_new + i_c]);
                    L_cur[i_c*S + i_s] = 1.0*l_sc_new;
                    
                    for (i_t=0; i_t < T; i_t++) {
                        M_cur[i_t*S + i_s] = M_Q[(Q+1)*i_t + i];
                        p_cur[i_t*S + i_s] = p_Q[(Q+1)*i_t + i];
                    }
                } //if(!(l_sc_cur == l_sc_new))
                
                //printf("\n L[%d, %d]=%f", i_s, i_c, L_cur[i_c*S + i_s]);
              
            }//for(i_c=1; i_c < C; i_c++)
         
        } //for(i_s=0; i_s < S; i_s++)
        
         
        delete w_Z; w_Z = NULL;
        delete L_s; L_s = NULL;
        
        delete M_Q; M_Q = NULL;
        delete p_Q; p_Q = NULL;
        
        delete prob_dist; prob_dist = NULL;
        delete l_set; l_set = NULL;
        
    }

    
    //UPDATE L AND Z JOINTLY-- NEED TO RECOMPUTE A, M AND P, AS WELL
    void fn_update_L_Z(int C, int S, int T, int Q, double *n, double *N, 
                       double *b, double *ppi, double *phi, double *w, 
                       double *Z_cur, double *A, double p_z0, double *M_cur, 
                       double *p_cur, double *L_cur)
    {
        double *w_Z; w_Z = new double[T*C];
        
        double *L_s; L_s = new double[C];
        double *Z_s; Z_s = new double[C];
        
        double *M_s; M_s = new double[T];
        double *p_s; p_s = new double[T];
        
        
        int i_t, i_c, ii_c, i_s;
        //ALTHOUGH THESE ARE INTEGERS, FINE.
        int l_sc_cur, l_sc_new, z_sc_cur, z_sc_new;  
        
        double u, prob_cur, prob_pro, prob_max, prob_sum, M_tmp, p_tmp;
        
        i_c = 0;
        for (i_t = 0; i_t < T; i_t ++) {
            w_Z[T*i_c + i_t] = w[T*i_c + i_t]*(p_z0);
        }
        
        for (i_t = 0; i_t < T; i_t ++) {
            for (i_c = 1; i_c < C; i_c++) {
                w_Z[T*i_c + i_t] = w[T*i_c + i_t];
            }
        }
        
        
        for(i_s=0; i_s < S; i_s++)
        {
            //printf("\n\n i_s=%d", i_s);
            for(i_c=1; i_c < C; i_c++)
            {
                //printf("\n i_s=%d, i_c=%d", i_s, i_c);
                //COPY ROW S
                for (ii_c=0; ii_c < C; ii_c++) {
                    L_s[ii_c]=L_cur[ii_c*S + i_s];
                    Z_s[ii_c]=Z_cur[ii_c*S + i_s];
                }
                
                //CURRENT L_SC AND Z_SC
                l_sc_cur = L_s[i_c];
                z_sc_cur = Z_s[i_c];
                
                //printf("\n CURRENT l_sc=%d, z_sc=%d", l_sc_cur, z_sc_cur);

                
                //SAMPLE L_SC_PRO
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();

                l_sc_new = u*(Q+1);
                
                //printf("\n PROPOSED u=%f, l_sc=%d", u, l_sc_new);
                
                //SAMPLE Z_SC_PRO
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();
                
                z_sc_new = u*(l_sc_new+1);

                //printf("\n PROPOSED u=%f, z_sc=%d", u, z_sc_new);
                

                //printf("\n PROPOSED l_sc=%d, z_sc=%d", l_sc_new, z_sc_new);
                //p(l_sc | pi)p(z_sc|l_sc)/q(z_sc)--CURRENT
                prob_cur = log(ppi[C*l_sc_cur + i_c]); 
                //p(l_sc | pi)p(z_sc|l_sc)/q(z_sc)--PROPOSED
                prob_pro = log(ppi[C*l_sc_new + i_c]); 
                
                //printf("\n prob=%f", prob_tmp);
                for (i_t=0; i_t < T; i_t++) {
                    //LIKELIHOOD--POISSON
                    prob_cur = prob_cur + N[i_t*S + i_s]*
                                          log(M_cur[i_t*S + i_s]) - 
                                          phi[i_t]*M_cur[i_t*S + i_s]*
                                          b[i_t*S + i_s]/2.0;  
                    //LIKELIHOOD--BINOMIAL
                    prob_cur = prob_cur + n[i_t*S + i_s]*
                                          log(p_cur[i_t*S + i_s]) + 
                                          (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                          log(1.0-p_cur[i_t*S + i_s]); 
                    
                    //printf("\n i_t=%d, M=%f, p=%f, prob=%f", i_t, 
                    //       M_cur[i_t*S + i_s], p_cur[i_t*S + i_s], prob_tmp);
                }
                
                L_s[i_c] = 1.0*l_sc_new;
                Z_s[i_c] = 1.0*z_sc_new;
                
                //COMPUTE M AND P WITH l_sc=i.q
                for(i_t=0; i_t < T; i_t++)
                {
                    M_tmp = 0.0;
                    p_tmp = 0.0;
                    
                    //printf("\n i_t=%d", i_t);
                    for (ii_c=0; ii_c < C; ii_c++) {
                        M_tmp = M_tmp + w[T*ii_c + i_t]*L_s[ii_c];
                        p_tmp = p_tmp + w_Z[T*ii_c + i_t]*Z_s[ii_c];
                    }
                    M_s[i_t] = M_tmp;
                    p_s[i_t] = p_tmp/M_tmp;
                }
                
                for (i_t=0; i_t < T; i_t++) {
                    //LIKELIHOOD--POISSON
                    prob_pro = prob_pro + N[i_t*S + i_s]*log(M_s[i_t]) - 
                                          phi[i_t]*M_s[i_t]*b[i_t*S + i_s]/2.0;
                    //LIKELIHOOD--BINOMIAL
                    prob_pro = prob_pro + n[i_t*S + i_s]*log(p_s[i_t]) + 
                                          (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                          log(1.0-p_s[i_t]); 
                    
                    //printf("\n i_t=%d, M=%f, p=%f", i_t, M_s[i_t], p_s[i_t]);
                    
                }
                
                //printf("\n prob_cur=%f, prob_pro=%f", prob_cur, prob_pro);

                
                //SAMPLE!!
                prob_max = mymax(prob_cur,prob_pro);

                prob_cur = exp(prob_cur - prob_max);
                prob_pro = exp(prob_pro - prob_max);
                
                prob_sum = prob_cur + prob_pro;
                
                
                //SAMPLE L_SC
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();
                
                //printf("\n u=%f, prob_cur=%f, prob_pro=%f", u, prob_cur, 
                //       prob_pro);
                
                if( u < prob_pro/prob_sum)
                {
                    //printf("\n A[%d]=%f, A[%d]=%f, ", C*l_sc_cur + i_c, 
                    //       A[C*l_sc_cur + i_c], C*l_sc_new + i_c, 
                    //       A[C*l_sc_new + i_c]);
                    
                    A[C*l_sc_cur + i_c] = A[C*l_sc_cur + i_c] - 1.0;
                    A[C*l_sc_new + i_c] = A[C*l_sc_new + i_c] + 1.0;
                    
                    //printf("\n A[%d]=%f, A[%d]=%f, ", C*l_sc_cur + i_c, 
                    //       A[C*l_sc_cur + i_c], C*l_sc_new + i_c, 
                    //       A[C*l_sc_new + i_c]);
                    L_cur[i_c*S + i_s] = 1.0*l_sc_new;
                    Z_cur[i_c*S + i_s] = 1.0*z_sc_new;
                    
                    for (i_t=0; i_t < T; i_t++) {
                        M_cur[i_t*S + i_s] = M_s[i_t];
                        p_cur[i_t*S + i_s] = p_s[i_t];
                    }
                } //if(!(l_sc_cur == l_sc_new))
                
                //printf("\n L[%d, %d]=%f, Z[%d, %d]=%f", i_s, i_c, 
                //       L_cur[i_c*S + i_s], i_s, i_c, Z_cur[i_c*S + i_s]);
                //printf("\n Z[%d, %d]=%f", i_s, i_c, Z_cur[i_c*S + i_s]);
                
            }//for(i_c=1; i_c < C; i_c++)
            
        } //for(i_s=0; i_s < S; i_s++)
        
        
        delete w_Z; w_Z = NULL;
        delete L_s; L_s = NULL;
        delete Z_s; Z_s = NULL;
        
        delete M_s; M_s = NULL;
        delete p_s; p_s = NULL;
    }

    void pbeta_test(double *a, double *b, double *x, double*prob)
    {
        //pbeta(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
        *prob = pbeta(*x, *a, *b, 1, 1); 
        
        //printf("\n prob=%f", *prob);
        
    }
    
    //UPDATE Z-- NEED TO RECOMPUTE A, AND P, AS WELL
    void fn_update_Z(int C, int S, int T, int Q, double *n, double *N, 
                     double *w, double *Z_cur, double p_z0, double *M, 
                     double *p_cur, double *L)
    {
        double *w_Z; w_Z = new double[T*C];
        double *Z_s; Z_s = new double[C];
        double *p_Q; p_Q = new double[(Q+1)*T];
        double *prob_dist; prob_dist = new double[Q+1];
        int *z_set; z_set= new int[Q+1];
        
        int i_t, i_c, ii_c, i_s, i_q, i;
        double z_sc_cur, z_sc_new, l_sc;
        int n_cnt;
        
        double p_tmp, u, prob_max, prob_sum, prob_tmp;
        
        i_c = 0;
        for (i_t = 0; i_t < T; i_t ++) {
            w_Z[T*i_c + i_t] = w[T*i_c + i_t]*(p_z0);
        }
        
        for (i_t = 0; i_t < T; i_t ++) {
            for (i_c = 1; i_c < C; i_c++) {
                w_Z[T*i_c + i_t] = w[T*i_c + i_t];
            }
        }
        
        //printf("\n S=%d, C=%d, T=%d, Q=%d", S, C, T, Q);
        
        for(i_s=0; i_s < S; i_s++)
        {
            //printf("\n i_s=%d", i_s);
            for(i_c=1; i_c < C; i_c++)
            {
                //printf("\n i_c=%d", i_c);
                
                n_cnt = 0;
                
                //COPY ROW S
                for (ii_c=0; ii_c < C; ii_c++) {
                    Z_s[ii_c]=Z_cur[ii_c*S + i_s];
                }
                
                z_sc_cur = Z_s[i_c];
                l_sc = L[i_c*S + i_s];
                
                //printf("\n CURRENT z_sc=%f, l_sc=%f", z_sc_cur, l_sc);
                //l_sc == l_sc_cur
                i_q = z_sc_cur;
                z_set[n_cnt] = i_q;
                
                
                //printf("\n i_q=%d, n_cnt=%d", i_q, n_cnt);
                
                //COMPUTE M AND P WITH l_sc=i.q
                //for(i_t=0; i_t < T; i_t++)
                //{
                //    printf("\n i_t=%d, M=%f, p=%f", i_t, M[i_t*S + i_s], 
                //           p_cur[i_t*S + i_s]);
                //}
                
                
                prob_tmp = 0.0;
                for (i_t=0; i_t < T; i_t++) {
                    //LIKELIHOOD--BINOMIAL
                    prob_tmp = prob_tmp + n[i_t*S + i_s]*
                                        log(p_cur[i_t*S + i_s]) + 
                                        (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                        log(1.0-p_cur[i_t*S + i_s]); 
                }
                
                //printf("\n i_q=%d, n_cnt=%d, z_set=%d, prob=%f", i_q, n_cnt, 
                //       z_set[n_cnt], prob_tmp);
                prob_dist[n_cnt] = prob_tmp;
                prob_max = prob_tmp;
                n_cnt = n_cnt + 1;
                
                for (i_q=0; i_q <= l_sc; i_q++) {
                    if (!(i_q==z_sc_cur)) {
                        //printf("\n i_q=%d, n_cnt=%d", i_q, n_cnt);
                        z_set[n_cnt] = i_q;
                        Z_s[i_c] = 1.0*i_q;
                        
                        //COMPUTE M AND P WITH l_sc=i.q
                        for(i_t=0; i_t < T; i_t++)
                        {
                            p_tmp = 0.0;
                            for (ii_c=0; ii_c < C; ii_c++) {
                                p_tmp = p_tmp + w_Z[T*ii_c + i_t]*Z_s[ii_c];
                            }
                            p_Q[(Q+1)*i_t + n_cnt] = p_tmp/M[i_t*S + i_s];
                            
                            //printf("\n i_t=%d, p=%f", i_t, 
                            //       p_Q[(Q+1)*i_t + n_cnt]);
                        }
                        
                        prob_tmp = 0.0;
                        for (i_t=0; i_t < T; i_t++) {
                            //LIKELIHOOD--BINOMIAL
                            prob_tmp = prob_tmp + n[i_t*S + i_s]*
                                                log(p_Q[(Q+1)*i_t + n_cnt]) + 
                                                (N[i_t*S + i_s]-n[i_t*S + i_s])*
                                                log(1.0-p_Q[(Q+1)*i_t + n_cnt]); 
                        }
                        prob_dist[n_cnt] = prob_tmp;
                        prob_max = mymax(prob_max,prob_tmp);
                        //printf("\n i_q=%d, n_cnt=%d, z_set=%d, prob=%f", i_q, 
                        //       n_cnt, z_set[n_cnt], prob_tmp);
                        
                        n_cnt = n_cnt + 1;
                    } //if (!(i_q==l_sc_cur)) {
                } //for (i_q=z_sc; i_q <= Q; i_q++)
                
                prob_sum = 0.0;
                for (i=0; i < n_cnt; i++) {
                    prob_dist[i] = exp(prob_dist[i] - prob_max);
                    prob_sum = prob_sum + prob_dist[i];
                }
                
                prob_dist[0] = prob_dist[0]/prob_sum;
                //printf("\n prob[%d]=%f", 0, prob_dist[0]);
                for (i=1; i < n_cnt; i++) {
                    prob_dist[i] = prob_dist[i-1] + prob_dist[i]/prob_sum;
                    //printf("\n prob[%d]=%f", i, prob_dist[i]);
                }
                
                //SAMPLE L_SC
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();
                
                //printf("\n u=%f", u);
                for (i=0; i < n_cnt; i++) {
                    if (u < prob_dist[i]) break;
                }
                
                z_sc_new = z_set[i];
                
                //printf("\n i_s=%d, i_c=%d, n_cnt=%d, u=%f, i=%d, z_sc_cur=%f, 
                //       z_sc_new=%f", i_s, i_c, n_cnt, u, i, 
                //       z_sc_cur, z_sc_new);
                
                
                if(!(z_sc_cur == z_sc_new))
                {
                    Z_cur[i_c*S + i_s] = 1.0*z_sc_new;
                    
                    for (i_t=0; i_t < T; i_t++) {
                        p_cur[i_t*S + i_s] = p_Q[(Q+1)*i_t + i];
                    }
                } //if(!(l_sc_cur == l_sc_new))
                
                //printf("\n Z_cur[%d, %d]=%f", i_s, i_c, 
                //       Z_cur[i_c*S + i_s]);
                
            }//for(i_c=1; i_c < C; i_c++)
        } //for(i_s=0; i_s < S; i_s++)
        
        delete w_Z; w_Z = NULL;
        delete Z_s; Z_s = NULL;
        
        delete p_Q; p_Q = NULL;
        
        delete prob_dist; prob_dist = NULL;
        delete z_set; z_set = NULL;
        
    }

    
    
    void fn_eval_log_post(double *n, double *N, int *TT, int *SS, int *CC, 
                          int *QQ, double *pi_mat, double *A, double *gam, 
                          double *beta, double *alpha, double *L, double *a,
                          double *b, double *phi, double *d, double *w, 
                          double *p0_l, double *a0_l, double *b0_l, 
                          double *p0_z, double *a0_z, double *b0_z,
                          double *M, double *p, double *loglike_1)
	{
		int S = *SS;
		int T = *TT;
		int C = *CC;
        int Q = *QQ;
		int i_s, i_t, i_c, i_q;
        
		double pi_star, log_like = 0.0;
		
        
        //PI
        for(i_c=1; i_c < C ; i_c++)
        {
            i_q = 2;
            log_like = log_like + ((*beta) + A[C*i_q + i_c] - 1.0)*
                                log(pi_mat[C*i_q + i_c]) + 
                                ((*alpha) + (S - A[C*i_q + i_c]) - 1.0)*
                                log(1.0-pi_mat[C*i_q + i_c]);
            
            //PI_STAR
            for (i_q=0; i_q < 2; i_q++) {
                pi_star = pi_mat[C*i_q + i_c]/(1.0-pi_mat[C*2 + i_c]);
                log_like = log_like + (gam[i_q] + A[C*i_q + i_c]-1.0)*
                                        log(pi_star);
            }
            
            for (i_q=3; i_q < (Q+1); i_q++) {
                pi_star = pi_mat[C*i_q + i_c]/(1.0-pi_mat[C*2 + i_c]);
                log_like = log_like + (gam[i_q-1] + A[C*i_q + i_c]-1.0)*
                                        log(pi_star);
            }
        } //for(i.c in 1:CC)

        for (i_s=0; i_s < S; i_s++) {
            for (i_c=1; i_c < C; i_c++) {
                log_like = log_like - log(L[S*i_c + i_s]+1.0);
            }
        }

        for (i_t = 0; i_t < T; i_t ++) {
            log_like = log_like + (*a-1.0)*log(phi[i_t]) - (*b)*phi[i_t];
            for (i_c=1; i_c < C; i_c++) {
                log_like = log_like + (d[i_c]-1.0)*log(w[T*i_c + i_t]);
            }
        }
        
        log_like = log_like + (*a0_l-1.0)*log(*p0_l) + 
                              (*b0_l-1.0)*log(1.0-*p0_z) - 
                              pbeta(*p0_l, *a0_z, *b0_z, 1, 1) + 
                              (*a0_z-1.0)*log(*p0_z) + 
                              (*b0_z-1.0)*log(1.0-*p0_z);
        
 		for (i_s=0; i_s < S; i_s++) {
			for (i_t=0; i_t < T; i_t++) {
                //LIKELIHOOD--POISSON
                log_like = log_like + N[i_t*S + i_s]*log(M[i_t*S + i_s]) - 
                                      phi[i_t]*M[i_t*S + i_s]/2.0;
                //LIKELIHOOD--BINOMIAL
				log_like = log_like + n[i_t*S + i_s]*log(p[i_t*S + i_s]) + 
                                      (N[i_t*S + i_s] - n[i_t*S + i_s])*
                                      log(1.0 - p[i_t*S + i_s]);  
			}
		}
        
		*loglike_1 = log_like;
	}

    
    
    //UPDATE P FOR GIVEN THETAS AND EVALUATE THE LIKELIHOOD
    double fn_eval_like(double *n, double *N, double *b, double *p, double *phi, 
                        double *M, int T, int S)
    {
        double lambda, loglike = 0.0;
        int i_s, i_t;
        
        for (i_s=0; i_s < S; i_s++) {
            for (i_t=0; i_t < T; i_t++) {
                lambda = phi[i_t]*M[i_t*S + i_s]*b[i_t*S + i_s]/2.0;
                loglike = loglike + N[i_t*S + i_s]*log(lambda) - lambda;
                loglike = loglike + n[i_t*S + i_s]*log(p[i_t*S + i_s]) + 
                                    (N[i_t*S + i_s] - n[i_t*S + i_s])*
                                    log(1.0 - p[i_t*S + i_s]);
            }
        }
        if (isnan(loglike)) loglike= -INFINITY;
        
        return(loglike);
    }
    

 	void fn_CNV_MCMC_1(double *alpha, double *beta, double *gam, int *QQ, 
                       double *a, double *b, double *d0, double *d, 
                       double *aa0_z, double *bb0_z,int *CC, double *ppi, 
                       double *phi, double *L, double *Z, double *A, 
                       double *p0_z, double *th, double *w, double *M, 
                       double *p,int *SS, int *TT, double *n, double *N, 
                       double *N_sum, double *BB, int *NN_iter)
	{
        int C = *CC;
        int S = *SS;
        int T = *TT;
        int Q = *QQ;
        
		double a0_z = *aa0_z;
		double b0_z = *bb0_z;
        
        int i_iter;
        
        //BURN-IN
        for (i_iter=0; i_iter < *NN_iter; i_iter ++) {
            //printf("\n i.iter=%d", i_iter);
            //UPDATE Z-- NEED TO RECOMPUTE A, AND P, AS WELL
            fn_update_Z(C, S, T, Q, n, N, w, Z, *p0_z, M, p, L);
            
            //UPDATE L-- NEED TO RECOMPUTE A, M AND P, AS WELL
            fn_update_L(C, S, T, Q, n, N, BB, ppi, phi, w, Z, A, *p0_z, M, 
                        p, L);
            
            //UPDATE PHI
            fn_update_phi(S, T, *a, *b, M, N_sum, BB, phi);

            //W_STAR=THETA
            fn_update_theta_p0(*d0, *d, a0_z, b0_z, Z, L, phi, T, S, C, n, N, 
                               BB, p0_z, th, w, p, M);
            
            //UPDATE PI
            fn_update_pi(*alpha, *beta, gam, Q, C, S, A, ppi);
            
            //UPDATE Z -- NEED TO RECOMPUTE P
            fn_update_multi_Z(C, S, T, n, N, w, Z, *p0_z, M, p, L);
            
            //UPDATE L -- NEED TO RECOMPUTE A, M AND P, AS WELL
            fn_update_multi_L(C, S, T, Q, n, N, BB, ppi, phi, w, Z, A, *p0_z, 
                              M, p, L);

            //UPDATE PHI
            fn_update_phi(S, T, *a, *b, M, N_sum, BB, phi);
            
            //W_STAR=THETA
            fn_update_multi_theta_p0(*d0, *d, a0_z, b0_z, Z, L, phi, T, S, C, 
                                     n, N, BB, p0_z, th, w, p, M);
            
            //UPDATE PI
            fn_update_pi(*alpha, *beta, gam, Q, C, S, A, ppi);
            
            //UPDATE L AND Z JOINTLY
            fn_update_L_Z(C, S, T, Q, n, N, BB, ppi, phi, w, Z, A, *p0_z, M, 
                          p, L);
        }
  	}

    
    void fn_CNV_MCMC_2(double *alpha, double *beta, double *gam, int *QQ, 
                       double *a, double *b, double *d0, double *d, 
                       double *aa0_z, double *bb0_z, int *CC, double *ppi, 
                       double *phi, double *L, double *Z, double *A, 
                       double *p0_z, double *th, double *w, double *M, 
                       double *p, int *SS, int *TT, double *n_tr, double *N_tr, 
                       double *N_tr_sum, double *BB_tr, double *n_te, 
                       double *N_te, double *N_te_sum, double *BB_te, 
                       double *loglike, int *NN_iter)
	{
        int C = *CC;
        int S = *SS;
        int T = *TT;
        int Q = *QQ;
        
		double a0_z = *aa0_z;
		double b0_z = *bb0_z;
        
        int i_iter;
        
        
        for (i_iter=0; i_iter < (*NN_iter); i_iter++) {
            //TRAINING
            //UPDATE Z-- NEED TO RECOMPUTE A, AND P, AS WELL
            fn_update_Z(C, S, T, Q, n_tr, N_tr, w, Z, *p0_z, M, p, L);
            
            //UPDATE L-- NEED TO RECOMPUTE A, M AND P, AS WELL
            fn_update_L(C, S, T, Q, n_tr, N_tr, BB_tr, ppi, phi, w, Z, A, *p0_z,
                        M, p, L);
            
            //UPDATE PHI
            fn_update_phi(S, T, *a, *b, M, N_tr_sum, BB_tr, phi);
            
            //W_STAR=THETA
            fn_update_theta_p0(*d0, *d, a0_z, b0_z, Z, L, phi, T, S, C, n_tr, 
                               N_tr, BB_tr, p0_z, th, w, p, M);
            
            //UPDATE PI
            fn_update_pi(*alpha, *beta, gam, Q, C, S, A, ppi);
            
            //UPDATE Z -- NEED TO RECOMPUTE P
            fn_update_multi_Z(C, S, T, n_tr, N_tr, w, Z, *p0_z, M, p, L);
            
            //UPDATE L -- NEED TO RECOMPUTE A, M AND P, AS WELL
            fn_update_multi_L(C, S, T, Q, n_tr, N_tr, BB_tr, ppi, phi, w, Z, A, 
                              *p0_z, M, p, L);
            
            //UPDATE PHI
            fn_update_phi(S, T, *a, *b, M, N_tr_sum, BB_tr, phi);
            
            //W_STAR=THETA
            fn_update_multi_theta_p0(*d0, *d, a0_z, b0_z, Z, L, phi, T, S, C, 
                                     n_tr, N_tr, BB_tr, p0_z, th, w, p, M);
            
            //UPDATE PI
            fn_update_pi(*alpha, *beta, gam, Q, C, S, A, ppi);
            
            //UPDATE L AND Z JOINTLY
            fn_update_L_Z(C, S, T, Q, n_tr, N_tr, BB_tr, ppi, phi, w, Z, A, 
                          *p0_z, M, p, L);
        }
        
        *loglike = fn_eval_like(n_te, N_te, BB_te, p, phi, M, T, S);
    }

    // from BayClone2_post.cpp
    
    double myabs(double a)
    {
		return (a > 0.0) ? a : (-1.0*a);
	}

    void fn_sum_Z_1(int *SS, int *CC, int *ind_NN, int *ind_set, int *Z_1, 
                    int *Z_2, double *min_dist)
	{
		int S = *SS;
		int C = *CC;
        int ind_N = *ind_NN;
        
        int i1, i2, i_c, i_s, i_N;
        double diff_s, mmin_dist, tmp_min;
        
        double *D; D = new double[C*C];
        
        for (i1=1; i1 <= C; i1++) {
            for (i2=1; i2 <= C; i2++) {
                //DIFFERENCE IN Z_C
                diff_s = 0.0;
                for (i_s=0; i_s < S; i_s++) {
                    if (Z_1[i1*S + i_s] != Z_2[i2*S + i_s]) {
                        diff_s = diff_s + 1.0;
                    }
                }
                D[(i2-1)*C + (i1-1)] = diff_s;
                //printf("\n i1=%d, i2=%d, D=%f", i1, i2, D[(i2-1)*C + (i1-1)]);
            }
        }
        
        mmin_dist = INFINITY;
        
        for (i_N = 0; i_N < ind_N; i_N++) {
            tmp_min = 0.0;
            for (i_c=0; i_c < C; i_c++) {
                i1 = i_c;
                i2 = ind_set[ind_N*i_c + i_N];
                tmp_min = tmp_min + D[i1*C + i2];
                //printf("\n i_1=%d, i_2=%d, D=%f, tmp_min=%f", i1, i2, 
                        //D[i1*C + i2], tmp_min);
            }
            
            if (tmp_min < mmin_dist) {
                mmin_dist = tmp_min;
            }
            //printf("\n tmp_min=%f, min_dist=%f", tmp_min, mmin_dist);
        }
        
        *min_dist = mmin_dist;
        delete[] D; D = NULL;
    }

/*	void fn_sum_Z_1(int *SS, int *TT, int *CC, int *ind_NN, int *ind_set, 
                    int *Z_1, int *Z_2, double *w_1, double *w_2, 
                    double *min_dist)
	{
		int T = *TT;
		int S = *SS;
		int C = *CC;
        int ind_N = *ind_NN;
        
        int i1, i2, i_c, i_t, i_s, i_N;
        double ave_w, diff_s, mmin_dist, tmp_min;
        
        double *D; D = new double[C*C];
        
        for (i1=1; i1 <= C; i1++) {
            for (i2=1; i2 <= C; i2++) {
                
                //AVERAGE W
                ave_w = 0.0;
                for (i_t = 0; i_t < T; i_t ++) {
                    ave_w = ave_w + myabs(w_1[i1*T + i_t] - w_2[i2*T + i_t]);
                }
                
                //DIFFERENCE IN Z_C
                diff_s = 0.0;
                for (i_s=0; i_s < S; i_s++) {
                    if (Z_1[i1*S + i_s] != Z_2[i2*S + i_s]) {
                        diff_s = diff_s + 1.0;
                    }
                }
                D[(i2-1)*C + (i1-1)] = diff_s*ave_w/(1.0*T);
                //printf("\n i1=%d, i2=%d, D=%f", i1, i2, D[(i2-1)*C + (i1-1)]);
            }
        }
        
        mmin_dist = INFINITY;
 
        for (i_N = 0; i_N < ind_N; i_N++) {
            tmp_min = 0.0;
            for (i_c=0; i_c < C; i_c++) {
                i1 = i_c;
                i2 = ind_set[ind_N*i_c + i_N];
                tmp_min = tmp_min + D[i1*C + i2];
                //printf("\n i_1=%d, i_2=%d, D=%f, tmp_min=%f", i1, i2, 
                        //D[i1*C + i2], tmp_min);
            }
            
            if (tmp_min < mmin_dist) {
                mmin_dist = tmp_min;
            }
            //printf("\n tmp_min=%f, min_dist=%f", tmp_min, mmin_dist);
        }
        
        *min_dist = mmin_dist;
        delete[] D; D = NULL;
    }
 */
    


}




