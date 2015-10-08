#include<Rcpp.h>
#include <libgen.h>
#include <string.h>
#include "corpus.h"
#include "state.h"
#include "utils.h"
using namespace Rcpp;
using namespace std;


//[[Rcpp::export]]
Rcpp::List RunHLDA(List Documents, int max_iter = 500, int max_time = 3600){
    
  if(Documents.size()==0){
    Rcpp::stop("No Networks list was provided");
  }
  //if the paramters are not correct use the default ones
  if(max_iter <= 0 && max_time <=0){
    max_iter = 500;
    max_time = 3600; //One hour :)
  }
    
  // Using the HDP implementation by Blei's lab 
  // with adaptation to Rcpp and removal of the GSL dependency
  // The original source code can be found at 
  // https://github.com/Blei-Lab/hdp/tree/master/hdp-faster


  time_t t; time(&t);    
  // Data parameters.
  //char* train_data = NULL;
  
  // Model parameters.
  double eta = 0.01;
  double gamma = 1.0;
  double alpha = 1.0;
  double gamma_a = 1.0;
  double gamma_b = 1.0;
  double alpha_a = 1.0;
  double alpha_b = 1.0;
  int sample_hyper = 0;
  
  
    // Reading one of the training data.        
    Corpus* c_train = new Corpus();
    c_train->read_data(Documents);
    
    
    time_t start, current;
    int total_time = 0;
    int iter = 0;
    
    
    Rcout << "Initializing HDP" << std::endl;
        
    HDP* hdp = new HDP();
    hdp->init_hdp(eta, gamma, alpha, c_train->size_vocab_);
    
    Rcout << "Setting up HDP state" << std::endl;
    hdp->setup_doc_states(c_train->docs_);    
    
    hdp->iterate_gibbs_state(false, false);
      
    while ((max_iter == -1 || iter < max_iter) && (max_time == -1 || total_time < max_time)) {
      ++iter;
      time (&start);
      
      // Iterations.
      hdp->iterate_gibbs_state(true, true);
      // Scoring the documents.
      //double likelihood = hdp->log_likelihood(NULL);
      hdp->compact_hdp_state();

      if (sample_hyper) hdp->hyper_inference(gamma_a, gamma_b, alpha_a, alpha_b);
      
      // Record the time.
      time(&current);
      int elapse = (int) difftime(current, start);
      total_time += elapse;
            
    }
        
    Rcpp::List Results(hdp->save_state());
    
    return Results;
}


