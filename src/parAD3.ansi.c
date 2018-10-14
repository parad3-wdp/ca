#include <math.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>


#define no_argument 0
#define required_argument 1
#define optional_argument 2

#define NEARLY_EQ_TOL(a,b,tol) (((a)-(b))*((a)-(b))<=(tol))
#define NEARLY_BINARY(a,tol) (NEARLY_EQ_TOL((a),1.0,(tol)) || NEARLY_EQ_TOL((a),0.0,(tol)))
#define NEARLY_ZERO_TOL(a,tol) (((a)<=(tol)) && ((a)>=(-(tol))))

#define MAX_LINE_LENGTH 128*1024

#define DEBUG 0

#define true 1
#define false 0

enum OptimizationStatus {
    STATUS_OPTIMAL_INTEGER = 0,
    STATUS_OPTIMAL_FRACTIONAL,
    STATUS_INFEASIBLE,
    STATUS_UNSOLVED
};

typedef struct {
    double *first;
    int *second;

} intdouble;

typedef struct {
    int nvars;
    int nfactors;
    int totallinks;
    int maxfactordegree;
    double *logpotentials;
    int *variabledegree;
    int **variableslinks;
    int *factorsdegree;
    int *factorsnlinksaccum;
    int *link2var;
    int *link2factor;
    double *cached_q_alpha;
    double * cached_q_alpha_posterior;
    double *posteriors;
} factorgraphstruct;


int isFactorAtMostOne(char *str) {
    if (str==strstr(str,"ATMOSTONE")) return true;
    return false;
}

int getFactorDegree(char *str, int *pos) {
    char *cad=&str[10];
    *pos=strchr(cad,32)-cad+1+10;
    return atoi(cad);
}

int getFactorVar(char *str, int *pos) {
    char *cad=&str[*pos];
    *pos=strchr(cad,32)-cad+1+*pos;
    return atoi(cad);
}

int diff_ms(struct timeval t1, struct timeval t2) {
    return (((t1.tv_sec - t2.tv_sec) * 1000000) +
            (t1.tv_usec - t2.tv_usec))/1000;
}

void InsertionSort(double *a, int*b, int length) {
    int i, j;
    double tmpfirst;
    int tmpsecond;

    for (i = 1; i < length; i++) {
        j = i;
        while (j > 0 && a[j - 1] > a[j]) {
            tmpfirst = a[j];
            tmpsecond = b[j];
            a[j] = a[j - 1];
            b[j] = b[j - 1];
            a[j - 1] = tmpfirst;
            b[j - 1] = tmpsecond;
            j--;
        }
    }
}


int readFile(char *filein, factorgraphstruct* f) {
    int nvars;
    int nfactors;
    double *logpotentials;
    int *variabledegree;
    int *variabledegree_count;
    int **variableslinks;
    int *factorsdegree;
    int *factorsnlinksaccum;
    int **factorsvariables;
    int *link2var;
    int *link2factor;
    double *cached_q_alpha;
    double * cached_q_alpha_posterior;
    double * posteriors;
    char *tmp;


    // BEGIN READING INPUT FACTOR GRAPH FILE

    FILE* fin;
    fin=fopen(filein,"r");
    int currentline=0;
    if (!fin) {
		printf("File not found: %s\n", filein);
        return 0;
    }
    char readline[MAX_LINE_LENGTH];
    currentline++;
    fgets(readline,MAX_LINE_LENGTH,fin);
    nvars=atoi(readline);
    if (nvars<1) {
   		printf("Invaild number of variables: %s\n", readline);
        return 0;
    }
    currentline++;
    fgets(readline,MAX_LINE_LENGTH,fin);
    nfactors=atoi(readline);
    if (nfactors<1) {
       	printf("Invaild number of factors: %s\n", readline);
        return 0;
    }
	printf("nvars: %i nfactors: %i\n", nvars,nfactors);

    logpotentials = (double*) malloc(sizeof(double)*nvars);
    variabledegree = (int *) malloc(sizeof(int)*nvars);
    variabledegree_count = (int *) malloc(sizeof(int)*nvars);
    variableslinks= (int**) malloc(sizeof(int*)*nvars);

	printf("Reading variables\n");
    // Reading nvars
	int i,j;
    for(i=0; i<nvars; i++) {
        currentline++;
        fgets(readline,MAX_LINE_LENGTH,fin);
        logpotentials[i]=strtod(readline,&tmp);
        variabledegree[i]=0;
        variabledegree_count[i]=0;
    }
    // Reading factors
    factorsdegree = (int*) malloc(sizeof(int)*nfactors);
    factorsnlinksaccum = (int*) malloc(sizeof(int)*nfactors);
    factorsvariables = (int**) malloc(sizeof(int*)*nfactors);
    int maxfactordegree=0;

    int line_parsing_position;
    factorsnlinksaccum[0]=0;
    int totallinks=0;

	printf("Reading factors\n");
    for(i=0; i<nfactors; i++) {
        currentline++;
        fgets(readline,MAX_LINE_LENGTH,fin);
        if (!isFactorAtMostOne(readline)) {
			printf("Factor not supported, Line : %i\n ",currentline); 
        }
        factorsdegree[i]=getFactorDegree(readline,&line_parsing_position);
        if (factorsdegree[i]>maxfactordegree) maxfactordegree=factorsdegree[i];
        if (i) factorsnlinksaccum[i]=factorsnlinksaccum[i-1]+factorsdegree[i-1];
        factorsvariables[i]= (int *) malloc(sizeof(int)*factorsdegree[i]);
        totallinks+=factorsdegree[i];
        for (j=0; j<factorsdegree[i]; j++) {
            factorsvariables[i][j]=getFactorVar(readline,&line_parsing_position)-1;
        }
    }


    link2var= (int *) malloc(sizeof(int)*totallinks);
    link2factor=(int *) malloc(sizeof(int)*totallinks);
    cached_q_alpha= (double *) malloc(sizeof(double)*totallinks);
    cached_q_alpha_posterior = (double *) malloc(sizeof(double)*totallinks);


    for (i=0; i<nfactors; i++) {
        for (j=0; j<factorsdegree[i]; j++) {
            int link=factorsnlinksaccum[i]+j;
            int v=factorsvariables[i][j];
            link2var[link]=v;
            link2factor[link]=i;
            variabledegree[v]++;
            cached_q_alpha[link]=-1.0;
            cached_q_alpha_posterior[link]=-1.0;
        }
    }


    for (i=0; i<nvars; i++) {
        variableslinks[i]= (int *) malloc(sizeof(int)*variabledegree[i]);
    }

    for (i=0; i<nfactors; i++) {
        for (j=0; j<factorsdegree[i]; j++) {
            int v=factorsvariables[i][j];
            int link=factorsnlinksaccum[i]+j;
            variableslinks[v][variabledegree_count[v]]=link;
            variabledegree_count[v]++;
        }
    }


    free(variabledegree_count);
	free(factorsvariables);
    f->nvars=nvars;
    f->nfactors=nfactors;
    f->totallinks=totallinks;
    f->maxfactordegree=maxfactordegree;
    f->logpotentials=logpotentials;
    f->variabledegree=variabledegree;
    f->variableslinks=variableslinks;
    f->factorsdegree=factorsdegree;
    f->factorsnlinksaccum=factorsnlinksaccum;
    f->link2var=link2var;
    f->link2factor=link2factor;
    f->cached_q_alpha=cached_q_alpha;
    f->cached_q_alpha_posterior=cached_q_alpha_posterior;
    f->posteriors=(double*) malloc(sizeof(double)*nvars);
    return 1;
}




int RunAD3(double *value, double *upper_bound,factorgraphstruct* f,double ad3_eta_, double residual_threshold,int ad3_max_iterations_, int threads) {
    struct timeval start, end;
    gettimeofday(&start, NULL);
    double lower_bound =-1e100;
    int ad3_adapt_eta_=true;
    int verbosity_=2;
    double max_eta = 100.0;											 // Private for every thread - constant
    double min_eta = 1e-3;											 // Private for every thread - constant
    double gamma_primal = 100.0; // 10.0							 // Private for every thread - constant
    double gamma_dual = 10.0;										 // Private for every thread - constant
    double factor_step = 2.0;										 // Private for every thread - constant
    double tau = 1.0;												 // Private for every thread - constant
    int num_iterations_adapt_eta = 10; // 1							 // Private for every thread - constant
    int num_iterations_reset = 50;									 // Private for every thread - constant
    double cache_tolerance = 1e-12;								 	 // Private for every thread - constant
    int caching = true;  							    			 // Private for every thread - constant
    int skipfactorcheck=true;
    int nvars;														 // Private for every thread - constant
    int *variabledegree;											 // Private for every thread - constant
    int **variablefactors;										     // Not used, reserved
    int **variableslinks; 										     // Private for every thread - constant

    double *logpotentials;											 // Private for every thread - constant

    int nfactors;													 // Private for every thread - constant
    //int *factorstype;												 // Not used, reserved
    int *factorsdegree;												 // Private for every thread - constant
    int *factorsnlinksaccum;										 // Private for every thread - constant
    int maxfactordegree=0;											 // Private for every thread - constant
    double *logpotentials_fac;										 // Private for every thread - variable but only used for store intermediate data
    double *invdj;													 // Private for every thread - constant


    int totallinks;  												 // Private for every thread - constant
    intdouble *lastsort;	  									     // Private for every trhead if loop2 static


    // Privates


    int *link2var;											 // Private for every thread - constant
    int *link2factor;
    double dual_obj,dual_residual,primal_residual;                   // Shared
    double *vposteriors;											 // Shared
    double *cached_q_alpha;									 // Shared
    double *cached_q_alpha_posterior;							 // Shared
    double *lambdas,*maps_av,*maps_sum2,*diffmaps;			 // Shared
    int *variable_is_active;										 // Shared
    int *factor_is_active;
    int optimal = false;											 // Shared
    int reached_lower_bound = false;								 // Shared
    int eta_changed = true;										 // Shared
    double *maps_sum;												 // Shared

    double dual_obj_best = 1e100, primal_rel_obj_best = -1e100;		 // Shared
    double primal_obj_best = -1e100;								 // Shared
    int num_iterations_compute_dual = 50;					   		 // Shared
    double extra_score =0.0;										 // Shared
    double eta = ad3_eta_;											 // Shared
    int compute_dual = false;										 // Shared
    int compute_primal_rel = false;							 	 // Shared
    int  maxiterations= ad3_max_iterations_; 						 // Shared
    double *posteriors;
    int activefactors;
    int activevars;

	printf("starting AD3\n");
    nvars=f->nvars;
    nfactors=f->nfactors;
    totallinks=f->totallinks;
    maxfactordegree=f->maxfactordegree;
    logpotentials=f->logpotentials;
    variabledegree=f->variabledegree;
    variableslinks=f->variableslinks;
    factorsdegree=f->factorsdegree;
    factorsnlinksaccum=f->factorsnlinksaccum;
    link2var=f->link2var;
    link2factor=f->link2factor;
    cached_q_alpha=f->cached_q_alpha;
    cached_q_alpha_posterior=f->cached_q_alpha_posterior;
    posteriors=f->posteriors;


//    factor_is_active = new bool[nfactors];
	factor_is_active = (int *) malloc(sizeof(int)*nfactors);
    lastsort = (intdouble*) malloc (sizeof(intdouble)*nfactors);
    int i,j;
    for (i=0; i<nfactors; i++) {
        factor_is_active[i]=true;
        lastsort[i].first = (double*) malloc(sizeof(double)*factorsdegree[i]);
        lastsort[i].second = (int*) malloc(sizeof(int)*factorsdegree[i]);
        for (j=0; j<factorsdegree[i]; j++) {
            lastsort[i].second[j]=j;
        }
    }
    for (i=0; i < nvars; ++i) {
        if ((variabledegree[i] == 0) && (logpotentials[i] > 0))
            extra_score += logpotentials[i];
    }


 /*   lambdas = new double[totallinks];
    maps_sum2= new double[totallinks];
    maps_av = new double[nvars];
    diffmaps = new double[nvars];
    maps_sum = new double[nvars];
    variable_is_active = new bool[nvars];
    factor_is_active = new bool[nfactors];
    vposteriors=new double[nvars]; */
	
	lambdas = (double *) malloc (sizeof(double)*totallinks);
	maps_sum2 = (double *) malloc (sizeof(double)*totallinks);
	maps_av = (double *) malloc (sizeof(double)*nvars);
	diffmaps = (double *) malloc (sizeof(double)*nvars);
	maps_sum = (double *) malloc (sizeof(double)*nvars);
	
    variable_is_active =  (int *) malloc (sizeof(int)*nvars);
    factor_is_active =  (int *) malloc (sizeof(int)*nfactors);
	vposteriors=(double *) malloc (sizeof(double)*nvars);
	
    //double *savemaps = new double[totallinks];
	double *savemaps = (double *) malloc (sizeof(double)*totallinks);
	
    double *one = savemaps;
    double *two = cached_q_alpha_posterior;
    int originalmapping=true;

    for (i=0; i<totallinks; i++) {
        lambdas[i]=0.0;
        savemaps[i]=0.0;
    }
    for (i=0; i<nvars; i++) {
        maps_av[i]=0.5;
        maps_sum[i]=0;
    }
    for (i=0; i<nfactors; i++) {
        factor_is_active[i]=true;
    }

    int currentt;
    omp_set_num_threads(threads);

    #pragma omp parallel
    {


        int i,j,t;

        for (t = 0; t < maxiterations; ++t) {

            #pragma omp master
            {
                dual_residual = 0.0;
                primal_residual = 0.0;
            }

            #pragma omp barrier
            int recompute_everything=(((t % num_iterations_reset)==0) || eta_changed);

            if (eta_changed) {
                // LOOP 1
                // LOOP WAS VECTORIZED
				int link;
                #pragma omp  for		
                for (link = 0; link < totallinks; ++link) {
                    int k=link2var[link];
                    double val=logpotentials[k]/ (double) (variabledegree[k])+2.0*lambdas[link];
                    cached_q_alpha[link]=maps_av[k] + val / (2.0 * eta);

                }
            }
            double tau_projection=0.0;

            // LOOP 2
			int j;
            #pragma omp  for
            //  loop was not vectorized: not inner loop
            for (j = 0; j < nfactors; ++j) {
                //if (!factor_is_active[j]) continue;
                int skip=((!recompute_everything) && (!factor_is_active[j])) ;
                if (skip) continue;
                double s = 0.0;
                double s2= 0.0;
                tau_projection=0.0;
                int fd=factorsdegree[j];
                int linkbase=factorsnlinksaccum[j];
                int f,i;
                intdouble *ls=&lastsort[j];
				int link;
                // REMAINDER LOOP WAS VECTORIZED
				#pragma simd reduction (+:s)
				for (link=linkbase; link < linkbase+fd; ++link) {
                    cached_q_alpha_posterior[link] = (cached_q_alpha[link]<0.0)?0.0:cached_q_alpha[link];
                    s += cached_q_alpha_posterior[link];

                }
                if (s > 1.0) {
                    // LOOP WAS VECTORIZED
                    for (link=linkbase; link < linkbase+fd; ++link) {
                        s2 += cached_q_alpha[link];
                        ls->first[link-linkbase] = cached_q_alpha[linkbase+ls->second[link-linkbase]];
                    }
                    InsertionSort(ls->first,ls->second,fd);
                    for (i = 0; i < fd; i++) {
                        tau_projection = (s2 - 1.0)/(fd-i);
                        if (ls->first[i] > tau_projection) break;
                        s2 -= ls->first[i];
                    }
                    // LOOP WAS VECTORIZED
                    for (link=linkbase; link < linkbase+fd; link++) {
                        cached_q_alpha_posterior[link] = (cached_q_alpha[link]<tau_projection)?0.0:cached_q_alpha[link]-tau_projection;
                    }

                }
            }

            #pragma omp master
            {
                savemaps=(originalmapping)?two:one;
                cached_q_alpha_posterior=(originalmapping)?one:two;
                originalmapping=!originalmapping;
            }
            #pragma omp barrier

            // LOOP 3
            if (t == 0 || eta_changed) {
				int k;
                #pragma omp for
                //  LOOP WAS VECTORIZED
                for (k=0; k<nvars; k++) {
                    variable_is_active[k] = true;
                }
            } else {
					int k;
                    #pragma omp for
                    for (k=0; k<nvars; k++) {
                        variable_is_active[k] = false;
                        int l;
						for (l=0; l <variabledegree[k]; l++) {
                            int link=variableslinks[k][l];
                            double previous_cached_q_alpha=cached_q_alpha_posterior[link];
                            if (!(NEARLY_EQ_TOL(savemaps[link], maps_av[k], cache_tolerance))) {
                                variable_is_active[k] = true;
                                break;
                            }
                            if (!NEARLY_EQ_TOL(savemaps[link],previous_cached_q_alpha , cache_tolerance)) {
                                variable_is_active[k] = true;
                                break;
                            }
                            if (!NEARLY_BINARY(savemaps[link], cache_tolerance)) {
                                variable_is_active[k] = true;
                                break;
                            }
                        }

                    }
            }

            // LOOP 4			
		    int k;
            #pragma omp for
            for (k=0; k<nvars; k++) {
                double s=0.0;
                int l;
				for (l=0; l <variabledegree[k]; l++) {
                    int link=variableslinks[k][l];
                    maps_sum[k]+=savemaps[link]-cached_q_alpha_posterior[link];
                }
            }
            #pragma omp barrier

          
            // LOOP 5
			int i;
            #pragma omp for reduction(+:dual_residual)
            // loop was not vectorized: existence of vector dependence
            for (i = 0; i < nvars; ++i) {
                int variable_degree = variabledegree[i];
                if (!variable_is_active[i]) {
                    if (variable_degree == 0) {
                        maps_av[i] = (logpotentials[i] > 0)? 1.0 : 0.0;
                    }
                    continue;
                }
                double map_av_prev = maps_av[i];
                if (variable_degree == 0) {
                    maps_av[i] = (logpotentials[i] > 0)? 1.0 : 0.0;
                } else {
                    maps_av[i] = maps_sum[i] / (double) (variable_degree);
                }
                double diff = maps_av[i] - map_av_prev;
                diffmaps[i]=diff;
                // vector dependence: assumed ANTI dependence between dual_residual line *** and dual_residual line ***
                // vector dependence: assumed FLOW dependence between dual_residual line *** and dual_residual line ***
                dual_residual += variabledegree[i] * diff * diff;
            }

			// LOOP 6 : LOOP WAS VECTORIZED
            #pragma omp for
            for (k=0; k<nfactors; k++) factor_is_active[k] = false;
            // LOOP 7 : LOOP WAS VECTORIZED
			int link;
            #pragma omp for reduction(+:primal_residual)
             for (link=0; link<totallinks; link++) {

                int k=link2var[link];
                if (!variable_is_active[k]) continue;
                int f=link2factor[link];
                factor_is_active[f]=true;
                double diff_penalty = (savemaps[link] - maps_av[k]);
                cached_q_alpha[link]  += diffmaps[k] - tau * diff_penalty;

                lambdas[link] -= tau * eta * diff_penalty;
                primal_residual += diff_penalty * diff_penalty;
            }
            #pragma omp master
            {
                // END Iterate VARS
                primal_residual = sqrt(primal_residual / totallinks);
                dual_residual = sqrt(dual_residual / totallinks);

                // If primal residual is low enough or enough iterations
                // have passed, compute the dual.
                compute_dual = false;
                compute_primal_rel = false;
                // TODO: && dual_residual < residual_threshold?
                if (primal_residual < residual_threshold) {
                    compute_dual = true;
                    compute_primal_rel = true;
                } else if (t > 0 && 0 == (t % num_iterations_compute_dual)) {
                    compute_dual = true;
                }
                dual_obj = 1e100;
            }
            #pragma omp barrier


            if (compute_dual) {
                dual_obj = 0.0;
                // LOOP 8				
                #pragma omp  for reduction(+:dual_obj)

                for (j = 0; j < nfactors; ++j) {
                    double delta = 0.0;
                    double maxlp=-1e100;
                    // INNER LOOP 9: LOOP WAS VECTORIZED
                    int i;
					for (i = 0; i < factorsdegree[j]; ++i) {
                        int link = factorsnlinksaccum[j]+i;
                        int k=link2var[link];
                        int variable_degree = variabledegree[k];
                        double currentlp = logpotentials[k] / (double) (variabledegree[k])  + 2.0 * lambdas[link];
                        if (maxlp < currentlp) maxlp=currentlp;
                        delta -= lambdas[link];
                    }
                    double val = 0.0;
                    int all_zeros = true;
                    if (maxlp > 0.0) {
                        val += maxlp;
                        all_zeros = false;
                    }
                    dual_obj += val + delta;
                }
                #pragma omp master
                dual_obj += extra_score;
            }
            #pragma omp master
            {
                int dobreak =false;
                double primal_rel_obj = -1e100;
                if (compute_primal_rel) {
                    primal_rel_obj = 0.0;
                    // LOOP WAS VECTORIZED
                    int i;
					for (i = 0; i < nvars; ++i) {
                        primal_rel_obj += maps_av[i] * logpotentials[i];
                    }
                }
                // Compute primal objective.
                double primal_obj = -1e100;
                int compute_primal = false;
                if (compute_primal) {
                    //TODO
                    if (primal_obj > primal_obj_best) {
                        primal_obj_best = primal_obj;
                    }
                }
                if (dual_obj_best > dual_obj) {
                    dual_obj_best = dual_obj;
                    // LOOP WAS VECTORIZED
                    int i;
					for (i = 0; i < nvars; ++i) {
                        vposteriors[i] = maps_av[i];
                    }
                    if (dual_obj_best < lower_bound) {
                        reached_lower_bound = true;
                        dobreak=true;

                    }
                }
                if (!dobreak) {
                    if (primal_rel_obj_best < primal_rel_obj) {
                        primal_rel_obj_best = primal_rel_obj;
                    }
                    if (compute_dual) {
                        double mult=1.0;
                        gettimeofday(&end, NULL);
                        if (verbosity_ > 1) {
    							printf("Iteration = %i\tDual obj = %e\tPrimal rel obj = %e\tPrimal obj = %e\tDual residual = %e\tPrimal residual = %e\t",t,dual_obj,primal_rel_obj,primal_obj,dual_residual,primal_residual);
							printf("Best dual obj = %e\tBest primal rel obj = %e\ttBest primal obj = %e\tCached factors = 0\teta = %e\tChanged eta = ",dual_obj_best,primal_rel_obj_best,primal_obj_best,eta);
							if (eta_changed) printf("true"); else printf("false");
							printf("\tTime = %e sec.\n",((double) mult*diff_ms(end,start))/1000.0);
                        }
                    }
                    //double gap = dual_obj_best - primal_rel_obj_best;

                    // If both primal and dual residuals fall below a threshold,
                    // we are done. TODO: also use gap?
                    if (dual_residual < residual_threshold &&
                            primal_residual < residual_threshold) {
                        //LOOP WAS VECTORIZED
						int i;
                        for (i = 0; i < nvars; ++i) {
                            vposteriors[i] = maps_av[i];
                        }
                        optimal = true;
                        dobreak=true;
                    }
                }
                if (!dobreak) {

                    // Adjust the stepsize if residuals are very asymmetric.
                    eta_changed = false;
                    if (ad3_adapt_eta_ && 0 == (t % num_iterations_adapt_eta)) {
                        if (primal_residual > gamma_primal * dual_residual) {
                            if (eta < max_eta) {
                                eta *= factor_step;
                                eta_changed = true;
                            }
                        } else if (dual_residual > gamma_dual * primal_residual) {
                            if (eta > min_eta) {
                                eta /= factor_step;
                                eta_changed = true;
                            }
                        }

                    }
                } // From if !dobreak
                if (dobreak) maxiterations=0;
                else currentt=t+1;
            } // From pragma omp master
            #pragma omp barrier


        }// END Main loop
    } // From pragma omp parallel


    int fractional = false;
    *value = 0.0;
    for (i = 0; i < nvars; ++i) {
        if (!NEARLY_BINARY((vposteriors)[i], 1e-12)) fractional = true;
        *value += logpotentials[i] * (vposteriors)[i];
        posteriors[i]=vposteriors[i];
    }


    if (verbosity_ > 1) {
		printf("Solution value after %i iterations (AD3) = %e\n",currentt,*value);
    }
    *upper_bound = dual_obj_best;

    gettimeofday(&end, NULL);
    if (verbosity_ > 1) {
		printf("Took %e sec.\n",((double) diff_ms(end,start))/1000.0);
    }

    if (optimal) {
        if (!fractional) {
            if (verbosity_ > 1) {
				printf("Solution is integer\n");
            }
            return STATUS_OPTIMAL_INTEGER;
        } else {
            if (verbosity_ > 1) {
				printf("Solution is fractional\n");
            }
            return STATUS_OPTIMAL_FRACTIONAL;
        }
    } else {
        if (reached_lower_bound) {
            if (verbosity_ > 1) {
				printf("Reached lower bound: %e.\n",lower_bound);
            }
            return STATUS_INFEASIBLE;
        } else {
            if (verbosity_ > 1) {
				printf("Solution is only approximate.\n");
            }
            return STATUS_UNSOLVED;
        }
    }
}


int main(int argc, char * argv[])
{
    char message[1000];
	
	sprintf(message,"Usage: %s " \
                    "--file_graphs=[IN] --file_posteriors=[OUT] " \
                    "(--max_iterations=[NUM] --eta=[NUM]  " \
                    "--residual_threshold=[NUM] --threads=[NUM])", argv[0]);


    const struct option longopts[] =
    {
        {"file_graphs",     required_argument,        0, 'i'},
        {"file_posteriors", required_argument,        0, 'o'},
        {"max_iterations",  required_argument,  	  0, 'm'},
        {"eta",  			required_argument,  	  0, 'e'},
        {"residual_threshold",  required_argument,  	  0, 'r'},
        {"threads",  required_argument,  	  0, 't'},
        {0,0,0,0},
    };

    int index;
    int iarg=0;
    opterr=1;

    char filein[4351];
    char fileout[4351];
    char *tmp;
    int threads;



    int iset,oset,mset,eset,rset,tset;
    // Input and output filenames have no default values
    // Eta,max_iterations and residual threshold have.
    iset=oset=0;
    mset=eset=rset=tset=1;

    int max_iterations=1000;
    double eta=1000;
    double residual_threshold=1e-06;

    while(iarg != -1)
    {
        iarg = getopt_long(argc, argv, "iomert", longopts, &index);

        switch (iarg)
        {
        case 'i':
            memcpy(filein,optarg,strlen(optarg));
            iset=1;
            break;

        case 'o':
            oset=1;
            memcpy(fileout,optarg,strlen(optarg));
            break;

        case 'm':
            max_iterations=atoi(optarg);
            if (max_iterations<1) {
                mset=0;
				printf("Error determining max iterations %s\n", optarg);
            }
            break;

        case 'e':
            eta=strtod(optarg,&tmp);
            if (eta==0.0) {
                eset=0;
				printf("Invalid eta %s\n", optarg);
            }
            break;

        case 'r':
            residual_threshold=strtod(optarg,&tmp);
            if (residual_threshold==0.0) {
                rset=0;
				printf("Invalid residual_threshold %s\n", optarg);
            }
            break;
        case 't':
            threads=atoi(optarg);
            break;
        }
    }

    if (iset*oset*mset*eset*rset!=1) {
		printf("%s\n", message);
        return 0;
    }
    printf("Parameter summary: \n\n");
    printf("File graphs: %s\n",filein);
	printf("File posteriors: %s\n",fileout);
    printf("ETA:  %e\n",eta);
    printf("max iterations: %i\n",max_iterations);
    printf("residual threshold: %e\n",residual_threshold);
    printf("threads: %i\n\n",threads);

    // END PARSING PARAMETERS
    double value,upper_bound;
    factorgraphstruct fgdata;
    readFile(filein,&fgdata);
    RunAD3(&value,&upper_bound,&fgdata,eta,residual_threshold,max_iterations,threads);

	FILE * f;
	f=fopen(fileout,"w");
	if (f!=NULL) {
    	int i;
		for (i=0; i<fgdata.nvars; i++) {
			fprintf(f,"%f\n",fgdata.posteriors[i]);
    	}
		fclose(f);
	} else {
		printf("Error creating file %s\n",fileout);
	}
    return 0;
}
