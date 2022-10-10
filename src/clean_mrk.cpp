#include <Rcpp.h>

using namespace Rcpp;

//' Clean genetic markers
//' 
//' \code{clean_mrk} filters and/or recodes genetic markers such as single
//' nucleotide polymorphism markers.
//'
//' @param data Either a \code{matrix} or a \code{data.frame} containing marker
//'             genotype data. Markers must be arranged in rows and genotypes
//'             (i.e. lines, hybrids, clones or individuals) in columns. Marker 
//'             names can be provided as row names or in the first column (see 
//'             argument \code{row_names}).
//' @param na_freq Upper limit for the frequency of missing values.
//' @param rm_mono \code{TRUE} to remove all monomorphic markers.
//' @param maf Lower limit for the minor allele frequency.
//' @param recode \code{TRUE} to recode to numeric code.
//' @param rm_dupes \code{TRUE} to remove duplicate markers. Markers sharing
//'                 identical information with missing values for different 
//'                 genotypes are also considered duplicates.
//' @param keep Character vector containing markers to keep.
//' @param impute \code{TRUE} to impute data with mean (for numeric data) or
//'               mode (for categorical data).
//' @param homo Character vector defining homozygous marker genotypes. Default
//'             is c("A", "C", "G", "T") with "A" standing for homozygous for 
//'             nucleotide A.
//' @param hetero Character vector defining heterozygous marker genotypes.
//'               Default is ("R", "Y", "S", "W", "K", "M") with "R" standing 
//'               for a heterozygous genotype with nucleotides A and G according
//'               to IUPAC code.
//' @param homo_num_min Integer value used for recoding of genotypes homozygous
//'                     for minor allele.
//' @param homo_num_maj Integer value used for recoding of genotypes homozygous
//'                     for major allele.
//' @param hetero_num Integer value used for recoding of heterozygous genotypes.
//' @param row_names Must be \code{TRUE} if marker names are provided as row
//'                  names.
//' @param verbose \code{TRUE} to print results at each step.
//'  
//' @return \code{data.frame} with filtered marker genotypes.
//' 
//' @export
// [[Rcpp::export]]
DataFrame clean_mrk(RObject data,
                    Nullable<NumericVector> na_freq = R_NilValue,
                    bool rm_mono = false,
                    Nullable<NumericVector> maf = R_NilValue, 
                    bool recode = false,
                    bool rm_dupes = false,
                    Nullable<CharacterVector> keep = R_NilValue, 
                    bool impute = false,
                    CharacterVector homo = CharacterVector::create("A", "C", "G", "T"),
                    CharacterVector hetero = CharacterVector::create("R", "Y", "S", "W", "K", "M"),
                    int homo_num_min = -1,
                    int homo_num_maj = 1,
                    int hetero_num = 0,
                    bool row_names = false,
                    bool verbose = true) {
  
  // Prepare object templates
  CharacterMatrix mat;
  CharacterMatrix geno;
  NumericMatrix geno_num;
  CharacterVector mrk_names;
  CharacterVector column_names;
  List freq;
  DataFrame geno_out;
  
  // Prepare data matrix
  if (Rf_isMatrix(data)) {
    mat  = as<CharacterMatrix>(data);
  }
  if (is<DataFrame>(data)) {
    DataFrame df = as<DataFrame>(data);
    mat = CharacterMatrix(df.nrow(), df.size());
    for (int i = 0; i < df.size(); i++) {
      mat(_, i) = CharacterVector(df[i]);
    }
    CharacterVector row_names_df = df.attr("row.names");
    rownames(mat) = row_names_df;
    CharacterVector column_names_df = df.attr("names");
    colnames(mat) = column_names_df;
  }
  
  // Extract genotypes and names
  if (row_names) {
    mrk_names = rownames(mat);
    geno = mat;
  } else {
    mrk_names = mat(_, 0);
    geno = mat(_, Range(1, mat.ncol() - 1));
  }
  column_names = colnames(mat);
  
  // Remove markers with missing values
  if (na_freq.isNotNull()) {
    if (verbose) {
      Rcout << "  Removing markers with missing values: ";
    }
    NumericVector na_freq_nn(na_freq);
    int nrow_geno = geno.nrow();
    int ncol_geno = geno.ncol();
    LogicalVector na_index(nrow_geno);
    for (int i = 0; i < nrow_geno; i++) {
      na_index[i] = sum(is_na(geno(i, _))) > ncol_geno * na_freq_nn[0];
    }
    if (keep.isNotNull()) {
      CharacterVector keep_nn(keep);
      na_index = na_index & !in(mrk_names, keep_nn);
    }
    CharacterMatrix geno_clean(sum(!na_index), ncol_geno);
    for (int i = 0, j = 0; i < nrow_geno; i++) {
      if (!na_index[i]) {
        geno_clean(j, _) = geno(i, _);
        j = j + 1;
      }
    }
    geno = geno_clean;
    mrk_names = mrk_names[!na_index];
    if (verbose) {
      Rcout << sum(na_index) << "\n";
    }
  }
  
  // Calculate genotype frequencies if needed
  if (rm_mono | maf.isNotNull() | recode | impute) {
    int nrow_geno = geno.nrow();
    freq = List(nrow_geno);
    for (int i = 0; i < nrow_geno; i++) {
      CharacterVector geno_i = geno(i, _);
      geno_i = geno_i[!is_na(geno_i)];
      IntegerVector geno_i_tbl = table(geno_i);
      NumericVector freq_i = as<NumericVector>(geno_i_tbl) / geno_i.size();
      freq_i.names() = geno_i_tbl.names();
      freq[i] = freq_i;
    }
  }
  
  // Remove monomorphic markers
  if (rm_mono) {
    if (verbose) {
      Rcout << "  Removing monomorphic markers:         ";
    }
    int size_freq = freq.size();
    LogicalVector mono_index(size_freq);
    for (int i = 0; i < size_freq; i++) {
      NumericVector freq_i = freq[i];
      if(freq_i.size() == 1) {
        mono_index[i] = true;
      }
    }
    if (keep.isNotNull()) {
      CharacterVector keep_nn(keep);
      mono_index = mono_index & !in(mrk_names, keep_nn);
    }
    CharacterMatrix geno_clean(sum(!mono_index), geno.ncol());
    for (int i = 0, j = 0; i < geno.nrow(); i++) {
      if (!mono_index[i]) {
        geno_clean(j, _) = geno(i, _);
        j = j + 1;
      }
    }
    geno = geno_clean;
    mrk_names = mrk_names[!mono_index];
    freq = freq[!mono_index];
    if (verbose) {
      Rcout << sum(mono_index) << "\n";
    }
  }
  
  // Remove markers with low MAF
  if (maf.isNotNull()){
    if (verbose) {
      Rcout << "  Removing markers due to low MAF:      ";
    }
    NumericVector maf_nn(maf);
    int size_freq = freq.size();
    LogicalVector maf_index(size_freq);
    for (int i = 0; i < size_freq; i++) {
      NumericVector freq_i = freq[i];
      if (freq_i.size() > 1) {
        CharacterVector freq_i_names = freq_i.names();
        LogicalVector is_homo = in(freq_i_names, homo);
        NumericVector freq_i_homo = freq_i[is_homo];
        if (is_true(any(!is_homo))) {
          NumericVector freq_i_hetero = freq_i[!is_homo];
          freq_i_homo = freq_i_homo + 0.5 * as<double>(freq_i_hetero);
        }
        maf_index[i] = max(freq_i_homo) > 1 - maf_nn[0];
      }
    }
    if (keep.isNotNull()) {
      CharacterVector keep_nn(keep);
      maf_index = maf_index & !in(mrk_names, keep_nn);
    }
    CharacterMatrix geno_clean(sum(!maf_index), geno.ncol());
    for (int i = 0, j = 0; i < geno.nrow(); i++) {
      if (!maf_index[i]) {
        geno_clean(j, _) = geno(i, _);
        j = j + 1;
      }
    }
    geno = geno_clean;
    mrk_names = mrk_names[!maf_index];
    freq = freq[!maf_index];
    if (verbose) {
      Rcout << sum(maf_index) << "\n";  
    }
  }
  
  // Recode to numeric code
  if (recode) {
    if (verbose) {
      Rcout << "  Recoding markers to numeric code\n";
    }
    int nrow_geno = geno.nrow();
    int ncol_geno = geno.ncol();
    CharacterVector major(nrow_geno);
    for (int i = 0; i < nrow_geno; i++) {
      NumericVector freq_i = freq[i];
      if (freq_i.size() > 0) {
        CharacterVector freq_i_names = freq_i.names();
        if (is_true(any(in(freq_i_names, homo)))) {
          NumericVector freq_i_homo = freq_i[in(freq_i_names, homo)];
          CharacterVector freq_i_homo_names = freq_i_homo.names();
          int major_index  = which_max(freq_i_homo);
          major[i] = freq_i_homo_names[major_index];
        }
      } else {
        major[i] = NA_STRING;
      }
    }
    geno_num = NumericMatrix(nrow_geno, ncol_geno);
    for (int i = 0; i < nrow_geno; i++) {
      CharacterVector geno_i = geno(i, _);
      CharacterVector major_i = as<CharacterVector>(major[i]);
      LogicalVector minor_index = !in(geno_i, major_i);
      LogicalVector major_index = in(geno_i, major_i);
      LogicalVector hetero_index = in(geno_i, hetero);
      LogicalVector na_index = is_na(geno_i);
      NumericVector geno_num_i(ncol_geno);
      geno_num_i[minor_index] = homo_num_min;
      geno_num_i[major_index] = homo_num_maj;
      geno_num_i[hetero_index] = hetero_num;
      geno_num_i[na_index] = NA_REAL;
      geno_num(i, _) = geno_num_i;
    }
    
    // Remove duplicate markers
    if (rm_dupes) {
      if (verbose) {
        Rcout << "  Removing duplicate markers:           ";
      }
      int nrow_geno_num = geno_num.nrow();
      int ncol_geno_num = geno_num.ncol();
      LogicalMatrix match_mat(nrow_geno_num);
      for (int i = 0; i < nrow_geno_num; i++) {
        for (int j = i ; j < nrow_geno_num; j++) {
          match_mat(i, j) = match_mat(j, i) = is_true(
            all(
              (geno_num(i, _) == geno_num(j, _)) |
                is_na(geno_num(i, _)) |
                is_na(geno_num(j, _))
            )
          );
        }
      }
      CharacterVector dupes_keep(nrow_geno_num);
      for (int i = 0; i < nrow_geno_num; i++) {
        LogicalVector match_vec = match_mat(i, _);
        if (sum(match_vec) > 1) {
          CharacterVector group_names = mrk_names[match_vec];
          int size_group = sum(match_vec);
          CharacterMatrix geno_num_group(size_group, ncol_geno_num);
          for (int j = 0, k = 0; j < nrow_geno_num; j++) {
            if(match_vec[j]) {
              geno_num_group(k, _) = geno_num(j, _);
              k = k + 1;
            }
          }
          IntegerVector na_count(size_group);
          for (int j = 0; j < size_group; j++) {
            na_count[j] = sum(is_na(geno_num_group(j, _)));
          }
          dupes_keep[i] = group_names[which_min(na_count)];
        } else {
          dupes_keep[i] = mrk_names[i];
        }
      }
      LogicalVector dupe_index = !in(mrk_names, dupes_keep);
      if (keep.isNotNull()) {
        CharacterVector keep_nn(keep);
        dupe_index = dupe_index & !in(mrk_names, keep_nn);
      }
      NumericMatrix geno_num_clean(sum(!dupe_index), ncol_geno_num);
      for (int i = 0, j = 0; i < nrow_geno_num; i++) {
        if (!dupe_index[i]) {
          geno_num_clean(j, _) = geno_num(i, _);
          j = j + 1;
        }
      }
      geno_num = geno_num_clean;
      mrk_names = mrk_names[!dupe_index];
      freq = freq[!dupe_index];
      if (verbose) {
        Rcout << sum(dupe_index) << "\n";
      }
    }
  }
  
  // Impute missing values with mean or mode
  if (impute) {
    if (verbose) {
      Rcout << "  Imputing missing values\n";
    }
    if (recode) {
      for (int i = 0; i < geno_num.nrow(); i++) {
        NumericVector geno_num_i = geno_num(i, _);
        if (is_true(any(!is_na(geno_num_i)))) {
          geno_num_i[is_na(geno_num_i)] = mean(na_omit(geno_num_i));
          geno_num(i, _) = geno_num_i;
        }
      }
    } else {
      for (int i = 0; i < geno.nrow(); i++) {
        CharacterVector geno_i = geno(i, _);
        if (is_true(any(!is_na(geno_i)))) {
          NumericVector freq_i = freq[i];
          CharacterVector mode = freq_i.names();
          mode = mode[which_max(freq_i)];
          geno_i[is_na(geno_i)] = mode;
          geno(i, _) = geno_i;
        }
      }
    } 
  }
  
  // Create and return DataFrame
  if (row_names) {
    if (recode) {
      geno_out = DataFrame::create(Named("") = geno_num);
    } else {
      geno_out = DataFrame::create(Named("") = geno);
    }
    geno_out.attr("row.names") = mrk_names;
  } else {
    if (recode) {
      geno_out = DataFrame::create(Named("") = mrk_names,
                                   Named("") = geno_num);
    } else {
      geno_out = DataFrame::create(Named("") = mrk_names,
                                   Named("") = geno);
    }
  }
  geno_out.attr("names") = column_names;
  if (verbose) {
    Rcout << "  -------------------------------------\n" <<
      "  Remaining markers:                    " << geno_out.nrow();
  }
  return geno_out;
}