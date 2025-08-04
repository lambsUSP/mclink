// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
#include <vector>
#include <algorithm>


// Convert sparse column in a dense vector with: -1, 0 and 1
arma::Col<int> convert_sign_vector(const arma::sp_vec &v, arma::uword n_rows) {
  arma::Col<int> result(n_rows, arma::fill::zeros);
  
  for (arma::uword i = 0; i < v.n_nonzero; ++i) {
    arma::uword row = v.row_indices[i];
    double val = v.values[i];
    result(row) = (val < 0) ? -1 : (val > 0 ? 1 : 0);
  }
  
  return result;
}


// Compare two vectors and calculate the proportion of equalities ignoring -1
double hamming_proportion(const arma::Col<int> &a, const arma::Col<int> &b) {
  arma::uword n = a.n_elem;
  arma::uword matches = 0, valid = 0;
  
  for (arma::uword i = 0; i < n; ++i) {
    if (a[i] < 0 || b[i] < 0) continue;
    valid++;
    if (a[i] == b[i]) matches++;
  }
  
  return (valid > 0) ? static_cast<double>(matches) / valid : -1.0;
}


// [[Rcpp::export]]
arma::mat relative_diff(const arma::sp_mat &M) {
  arma::uword ncols = M.n_cols;
  arma::uword nrows = M.n_rows;
  arma::mat result(ncols, ncols);
  
  // Convert all sparse columns in dense columns 
  std::vector<arma::Col<int>> sign_cols(ncols);
  for (arma::uword i = 0; i < ncols; ++i) {
    sign_cols[i] = convert_sign_vector(M.col(i), nrows);
  }
  
  for (arma::uword i = 0; i < ncols; ++i) {
    for (arma::uword j = i + 1; j < ncols; ++j) {
      double dist = 1 - hamming_proportion(sign_cols[i], sign_cols[j]);
      result(i, j) = dist;
      result(j, i) = dist;
    }
  }
 
  return result;
}


std::vector<int> mc_sample_rows(const int nrows, const size_t sample_size = 1000, const int min_distance = 1000){
  int end = nrows;
  std::vector<int> rows(end);
  std::iota(rows.begin(), rows.end(), 0); // Fill with the elements {0,1,...,M.nrows}
  std::vector<int> sample;
  std::vector<bool> avaiable(end, true); 
  
  sample.reserve(sample_size); // Reserve memory to at least sample_size elements
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  while(sample.size() < sample_size && end > 0){
    std::uniform_int_distribution<int> unif(0, end - 1);
    
    int rpos = unif(gen);
    int row_sampled = rows[rpos];
    
    if(!avaiable[row_sampled]){
      std::swap(rows[rpos], rows[end-1]);
      --end;
      continue;
    }
    
    sample.push_back(row_sampled);
    
    for(int i = std::max(0, row_sampled - min_distance); i <= std::min(static_cast<int>(avaiable.size()), row_sampled + min_distance); ++i){
      avaiable[i] = false;
    }
    
    std::swap(rows[rpos], rows[end-1]);
    --end;
  }
  
  return sample;
}


// [[Rcpp::export]]
arma::sp_mat mc_sample_matrix(const arma::sp_mat &M, const size_t sample_size = 1000, const int min_distance = 1000) {
  std::vector<int> rows = mc_sample_rows(M.n_rows, sample_size, min_distance);
  arma::uword n_sample = rows.size();
  arma::uword n_cols = M.n_cols;
  
  std::vector<arma::uword> row_inds;
  std::vector<arma::uword> col_inds;
  std::vector<double> values;
  
  // Optimization: allocate a little more memory than necessary.
  row_inds.reserve(M.n_nonzero);  
  col_inds.reserve(M.n_nonzero);
  values.reserve(M.n_nonzero);
  
  for (arma::uword i = 0; i < n_sample; ++i) {
    arma::uword row = rows[i];
    for (arma::sp_mat::const_row_iterator it = M.begin_row(row); it != M.end_row(row); ++it) {
      row_inds.push_back(i);          
      col_inds.push_back(it.col());   
      values.push_back(*it);          
    }
  }
  
  arma::umat locations(2, values.size());
  
  std::copy(row_inds.begin(), row_inds.end(), locations.row(0).begin());
  std::copy(col_inds.begin(), col_inds.end(), locations.row(1).begin());
  
  // Construct the sparse matrix directly
  arma::sp_mat sample_matrix(
      locations,
      arma::vec(values),
      n_sample, n_cols,
      true, false
  );
  
  return sample_matrix;
}


// [[Rcpp::export]]
arma::sp_mat mc_shuffle_matrix(const arma::sp_mat &M) {
  const arma::uword n_rows = M.n_rows;
  const arma::uword n_cols = M.n_cols;
  
  std::vector<arma::uword> row_inds;
  std::vector<arma::uword> col_inds;
  std::vector<double> values;
  
  row_inds.reserve(M.n_nonzero);
  col_inds.reserve(M.n_nonzero);
  values.reserve(M.n_nonzero);
   
  std::random_device rd;
  std::mt19937 gen(rd());
   
  for (arma::uword i = 0; i < n_rows; ++i) {
    std::vector<arma::uword> original_cols;
    std::vector<double> original_vals;
     
    for (arma::sp_mat::const_row_iterator it = M.begin_row(i); it != M.end_row(i); ++it) {
      original_cols.push_back(it.col());
      original_vals.push_back(*it);
    }
     
    arma::uword nnz = original_cols.size();
    if (nnz == 0) continue;
     
    // Generate shuffled columns 
    std::vector<arma::uword> shuffled_cols(n_cols);
    std::iota(shuffled_cols.begin(), shuffled_cols.end(), 0);
    std::shuffle(shuffled_cols.begin(), shuffled_cols.end(), gen);
    shuffled_cols.resize(nnz);
     
    // Save new triple (line, shuffled column and original value)
    for (arma::uword j = 0; j < nnz; ++j) {
      row_inds.push_back(i);
      col_inds.push_back(shuffled_cols[j]);
      values.push_back(original_vals[j]);
    } 
  } 
  
  arma::umat locations(2, values.size());
  std::copy(row_inds.begin(), row_inds.end(), locations.row(0).begin());
  std::copy(col_inds.begin(), col_inds.end(), locations.row(1).begin());
   
  arma::vec val_vec(values);
   
  return arma::sp_mat(locations, val_vec, n_rows, n_cols);
} 


