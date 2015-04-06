#include "corpus.h"
#include <assert.h>
#include <stdio.h>

Corpus::Corpus() {
  num_docs_ = 0;
  size_vocab_ = 0;
  num_total_words_ = 0;
}

Corpus::~Corpus() {
  free_corpus();
}

void Corpus::free_corpus() {
  for (int i = 0; i < num_docs_; ++i) {
    Document* doc = docs_[i];
    delete doc;
  }
  docs_.clear();

  num_docs_ = 0;
  size_vocab_ = 0;
  num_total_words_ = 0;
}


void Corpus::read_data(Rcpp::List Documents){
  free_corpus();
  int length = 0, n = 0, nd = 0;  
  
  int nbDocs = Documents.size();
  
  Rcpp::NumericVector vocab;

  for(int i=0;i< nbDocs; i++){
  NumericMatrix d = Documents[i];  
	length = d.ncol();
  Document * doc = new Document(length);	
    for (n = 0; n < length; ++n) {
      //fscanf(fileptr, "%10d:%10d", &word, &count);
	        
      doc->words_[n] = d(0,n);
      doc->counts_[n] = d(1,n);
      doc->total_ += d(1,n);
	    vocab.push_back(d(0,n));
      /*if (d(0,n) >= nw)  // Here we suppose that the words are ordered by number, it consider them to be 
        nw = d(0,n) + 1;*/
    }
    num_total_words_ += doc->total_;
    doc->id_ = nd; 
    docs_.push_back(doc);
    nd++;
  }
 
  num_docs_ += nd;
  NumericVector terms= unique(vocab);
  size_vocab_ = terms.size();
  Rcout << "number of networks :" << nd << std::endl;
  Rcout << "number of terms :" << terms.size() << std::endl;
  Rcout << "number of total words :" << num_total_words_ << std::endl;  
}


int Corpus::max_corpus_length() const {
  int max_length = 0;

  for (int d = 0; d < num_docs_; d++) {
    if (docs_[d]->length_ > max_length)
        max_length = docs_[d]->length_;
  }
  return max_length;
}

