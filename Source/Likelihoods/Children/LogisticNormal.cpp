/*
 * LogisticNormal.cpp
 *
 *  Created on: Oct 26, 2016
 *      Author: Zaita
 */

#include <Likelihoods/Children/LogisticNormal.h>

#include "Utilities/doubleCompare.h"
#include "Utilities/Math.h"
#include "Utilities/RandomNumberGenerator.h"

namespace niwa {
namespace likelihoods {

using std::set;

namespace dc = niwa::utilities::doublecompare;
namespace math = niwa::utilities::math;

LogisticNormal::LogisticNormal(Model* model)  : Likelihood(model) {
  parameters_.Bind<string>(PARAM_LABEL, &label_, "Label for Logisitic Normal Likelihood", "");
  parameters_.Bind<string>(PARAM_TYPE, &type_, "Type of likelihood", "");
  covariance_table_ = new parameters::Table(PARAM_COVARIANCE_MATRIX);

  parameters_.Bind<float>(PARAM_RHO, &rho_, "The auto-correlation parameter $\rho$", "");
  parameters_.BindTable(PARAM_COVARIANCE_MATRIX, covariance_table_, "User defined Covariance matrix", "",false,true);
  parameters_.Bind<float>(PARAM_SIGMA, &sigma_, "Sigma parameter in the likelihood", "");
  parameters_.Bind<bool>(PARAM_ARMA, &arma_, "Defines if two rho parameters supplied then covar is assumed to have the correlation matrix of an ARMA(1,1) process", "");
  parameters_.Bind<unsigned>(PARAM_BIN_LABELS, &bins_, "If no covariance matrix parameter then list a vector of bin labels that the covariance matrix will be built for, can be ages or lengths.", "",false);
  parameters_.Bind<bool>(PARAM_SEXED, &sexed_, "Will the observation be split by sex?", "",false);
  parameters_.Bind<bool>(PARAM_ROBUST, &robust_, "Robustification term for zero observations", "",false);
  parameters_.Bind<bool>(PARAM_SEPERATE_BY_SEX, &sep_by_sex_, "If data is sexed, should the covariance matrix be seperated by sex?", "",false);
  parameters_.Bind<bool>(PARAM_SEX_LAG, &sex_lag_, "if T and data are sexed, then the AR(n) correlation structure ignores sex and sets lag = |i-j|+1, where i and j index the age or length classes in the data.  Ignored if data are not sexed.", "",false);

  RegisterAsAddressable(PARAM_SIGMA, &sigma_);
  RegisterAsAddressable(PARAM_RHO, &rho_);
}

LogisticNormal::~LogisticNormal() {
  delete covariance_table_;
}

void LogisticNormal::DoValidate() {
  LOG_TRACE();
  if(!arma_ & (rho_.size() > 1)) {
    if(rho_[1] <= -1 || rho_[1] >= (1- fabs(rho_[0])))
      LOG_ERROR_P(PARAM_RHO) << "Incorrect values for rho. Rho cannot be less than -1 or the first value cannot be greater than (1 - |first value|)";
  }

  if(rho_.size() > 2)
    LOG_ERROR_P(PARAM_RHO) << "Can only have one value or two";
  /*
   * Build our covariance matrix with user defined values
   */
  if (parameters_.GetTable(PARAM_COVARIANCE_MATRIX)->has_been_defined()) {
    LOG_FINEST() << "Converting user defined covariance matrix";
    vector<vector<string>>& covariance_values_data = covariance_table_->data();
    int rows = covariance_values_data.size();
    int cols = covariance_values_data[0].size();
    if(rows != cols) {
      LOG_ERROR_P(PARAM_COVARIANCE_MATRIX) << "Covariance matrix must be square (equal rows and columns). rows = " << rows << " cols = " << cols;
    }
    covariance_matrix_(rows,cols);
    // Covariance should be an age by age matrix or age x sex by age x sex matrix; if comp data given by sex.
    LOG_FINEST() << " has " << covariance_values_data.size() << " rows defined";
    unsigned col_iter = 0;
    for (vector<string>& covariance_values_data_line : covariance_values_data) {
      unsigned bin = 0;
      vector<float> temp;
      if (!utilities::To<unsigned>(covariance_values_data_line[0], bin)) {
        LOG_ERROR_P(PARAM_COVARIANCE_MATRIX) << " value " << covariance_values_data_line[0] << " could not be converted in to an unsigned integer. It should be the year for this line";
      } else {
        for (unsigned i = 1; i < covariance_values_data_line.size(); ++i) {
          float value = 0;
          if (!utilities::To<float>(covariance_values_data_line[i], value)) {
            LOG_ERROR_P(PARAM_SCANNED) << " value (" << covariance_values_data_line[i] << ") could not be converted to a float";
          }
          covariance_matrix_(i - 1, col_iter) = value;
        }
      }
      bins_.push_back(bin);
      ++col_iter;
    }
    n_bins_ = bins_.size();
    if (sexed_)
      unique_bins_ /= 2;
    else
      unique_bins_ = n_bins_;
  } else {
    LOG_FINEST() << "Calculating Covariance matrix with user specified parameters";
    unique_bins_ = bins_.size();
    if (sexed_)
      n_bins_ = unique_bins_ * 2;
    else
      n_bins_ = unique_bins_;

    //initialise covariance matrix
    LOG_FINEST() << "number of bins = " << n_bins_ << " number of unique bins = " << unique_bins_;
    // Create a covariance matrix with the user defined parameters
    calculate_covariance();
  }

  LOG_FINEST() << "Printing Covariance top left triangle matrix";
  for (unsigned k = 0; k < covariance_matrix_.size1(); ++k) {
    for (unsigned j = 0; j < covariance_matrix_.size2(); ++j) {
      LOG_FINEST() << "row = " << k << " col = " << j << " val = " << covariance_matrix_(j,k) << " ";
    }
  }
}


/**
 * Simulate observed values
 *
 * @param comparisons A collection of comparisons passed by the observation
 */
void LogisticNormal::SimulateObserved(map<unsigned, map<string, vector<observations::Comparison> > >& comparisons) {
  LOG_FINE();
  // Generate a multivariate variable  X
  vector<float> year_totals(comparisons.size(), 0.0);
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  unsigned year_storer = 0;

  if (!parameters_.GetTable(PARAM_COVARIANCE_MATRIX)->has_been_defined() && rho_[0] == 0.0) {
    LOG_FINEST() << "sigma_ = " << sigma_;
    // We just need to generate an independent multivariate normal distribtuion
    for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
      LOG_FINE() << "Simulating values for year: " << year_iterator->first;
      for (auto second_iter = year_iterator->second.begin(); second_iter != year_iterator->second.end(); ++second_iter) {
        LOG_FINE() << "Simulating values for cell: " << second_iter->first;
        for (observations::Comparison& comparison : second_iter->second) {
          comparison.simulated_ = exp(rng.normal(0.0,1.0) * sigma_ + log(comparison.expected_));
          LOG_FINEST() << "random deviate " << rng.normal(0.0,1.0) << " age = " << comparison.age_ << " simuiulated val = " << comparison.simulated_  << " expected = " << comparison.expected_ ;
          year_totals[year_storer] += comparison.simulated_;
        }
      }
      ++year_storer;
    }
  } else {
    if (!parameters_.GetTable(PARAM_COVARIANCE_MATRIX)->has_been_defined())
      calculate_covariance();

    LOG_FINEST() << "Printing Covariance top left triangle";
    for (unsigned k = 0; k < covariance_matrix_.size1(); ++k) {
      for (unsigned j = 0; j <= k; ++j) {
        LOG_FINEST() << "row = " << k << " col = " << j << " val = " << covariance_matrix_(j,k) << " ";
      }
    }
    if (!DoCholeskyDecmposition())
      LOG_FATAL()<< "Cholesky decomposition failed. Cannot continue Simulating from a logisitic-normal likelihood";
    // Calculate multivariate normal distribution
    year_storer = 0;
    vector<float> normals(covariance_matrix_.size1(), 0.0);
    float row_sum;
    for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
      LOG_FINE() << "Simulating values for year: " << year_iterator->first;
      for (auto second_iter = year_iterator->second.begin(); second_iter != year_iterator->second.end(); ++second_iter) {
        LOG_FINE() << "Simulating values for cell: " << second_iter->first;
        unsigned nbins = year_iterator->second.size();
        for (observations::Comparison& comparison : second_iter->second) {
          for (unsigned i = 0; i < nbins; ++i) {
            normals[i] = rng.normal();
          }
          row_sum = 0.0;
          for (unsigned j = 0; j < nbins; ++j) {
            row_sum += covariance_matrix_lt(j,0) * normals[j];
          }
          comparison.simulated_ = exp(row_sum + log(comparison.expected_));
          //LOG_FINEST() << " age = " << comparison.age_ << " simuiulated val = " << comparison.observed_  << " expected = " << comparison.expected_  << " multivariate offset = " << row_sum << " log expectations = " << log(comparison.expected_);
          year_totals[year_storer] += comparison.simulated_;
        }
      }
      ++year_storer;
    }
  }
  // Do the logistic transformation to get our desired values.
  year_storer = 0;
  for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
    LOG_FINEST() << "year = " << year_iterator->first;
    for (auto second_iter = year_iterator->second.begin(); second_iter != year_iterator->second.end(); ++second_iter) {
      LOG_FINE() << "Simulating values for cell: " << second_iter->first;
      for (observations::Comparison& comparison : second_iter->second) {
        comparison.simulated_ /= year_totals[year_storer];
        LOG_FINEST() << "Simulated val = " << comparison.simulated_ << " expected = " << comparison.expected_;
      }
    }
    ++year_storer;
  }
  LOG_FINEST() << "check out the totals";
  for(auto num :year_totals)
    LOG_FINEST() << num ;
}



void LogisticNormal::calculate_covariance() {
  LOG_TRACE();
  unsigned n_phi = rho_.size();
  vector<float> rhovec;
  vector<vector<float>> covar;
  covar.resize(n_bins_, vector<float>(n_bins_, 0.0));

  LOG_FINEST() << "covariance rows = " << covar.size() << " cols = " << covar[0].size();
  // initialise covar as the identity matrix
  for (unsigned diag = 0; diag < covar.size(); ++ diag)
    covar[diag][diag] = 1.0 * sigma_ * sigma_;

  if (rho_[0] == 0.0) {
    // Do nothing all zeros in the off diagonal
  } else if (sexed_ && sep_by_sex_) {

      LOG_FINEST() << "Calculating covar for sexed patrition with sepbysex = true, covar dims = " << covar.size()  << " " << covar[0].size();
      for (unsigned i = 0; i <= 1; ++i) {
        rhovec = GetRho(rho_,unique_bins_,arma_);
          for (unsigned diag = 0; diag < (unique_bins_ - 1); ++diag) {
            for (int row = 0; row < (int)unique_bins_; ++row) {
              for (int col = 0; col < (int)unique_bins_; ++col) {
                if (fabs(row - col) == diag + 1)
                  covar[row + (i * unique_bins_)][col + (i * unique_bins_)] = rhovec[diag] * sigma_ * sigma_;
              }
            }
          }
        }

    } else if (sexed_ && sex_lag_) {
      vector<int> binlab;
      int ages = 1;
      for (unsigned i = 1; i <= n_bins_; ++i, ++ages) {
        if (ages > (int)unique_bins_)
          ages = 1;
        binlab.push_back(ages);
        cout << ages << " ";
      }
      if (arma_ && n_phi == 2)
        rhovec = GetRho(rho_,unique_bins_,arma_);
      else if (!arma_ && n_phi == 2)
        rhovec = GetRho(rho_,unique_bins_ + 1,arma_);
        vector<int> col_vec, row_vec, offset_vec;

        for (unsigned row = 0; row < n_bins_; ++row) {
          for (unsigned col = 0; col < n_bins_; ++col) {
            col_vec.push_back(binlab[col]);
            row_vec.push_back(binlab[row]);
          }
        }
        // now calculate an offset vector;
        for (unsigned index = 0; index < col_vec.size(); ++index) {
          offset_vec.push_back(fabs(row_vec[index] - col_vec[index]) + 1.0);
          cout << offset_vec[index] << " ";
        }
        if (covar.size() * covar[0].size() != offset_vec.size())
          LOG_CODE_ERROR() << "covar.size() * covar[0].size() != offset_vec.size(). get_mat_index will fail as index vector needs to be same size as rows * col of the mat.";
        for (unsigned index = 1; index <  unique_bins_; ++index) {
          for (unsigned index2 = 0; index2 < offset_vec.size(); ++index2) {
            if ((int)index == offset_vec[index2]) {
              vector<int> indexes = math::get_mat_index((int)covar.size(),(int)covar[0].size(), index2);
              covar[indexes[0]][indexes[1]] = rhovec[index - 1];
              //cout << "row = " << indexes[0] << " col = " << indexes[1] << " index = " << index2 << endl;
            }
          }
        }


      // Add the identity mat
      for (unsigned row = 0; row < n_bins_; ++row)
        covar[row][row] = 1.0;
      // Store in the Covariance matrix
      for (unsigned row = 0; row < n_bins_; ++row) {
        for (unsigned col = 0; col < n_bins_; ++col) {
          covar[row][col] *= sigma_ * sigma_;
        }
      }
    } else {
      // Unisex or sexed but treated like a single covariance matrix
      rhovec = GetRho(rho_,unique_bins_,arma_);
      for (int diag = 0; diag < (int)unique_bins_; ++ diag) {
        for (int row = 0; row <  (int)n_bins_; ++ row) {
          for (int col = 0; col < (int)n_bins_; ++ col) {
            if (fabs(row - col) == diag + 1) {
              covar[row][col] = rhovec[diag] * sigma_ * sigma_;
            }
          }
        }
      }
    }
    // I am struggling to intialise this matrix with ints. so do the way I know will work
    ublas::matrix<float> temp_covar(covar.size(),covar.size(),0.0);
    for(unsigned i = 0; i < covar.size(); ++i){
      for(unsigned j = 0; j < covar[i].size(); ++j) {
        temp_covar(i,j) = covar[i][j] ;
      }
    }
    covariance_matrix_ = temp_covar;
    LOG_FINEST() << "Finished building covariance matrix";
}

vector<float> LogisticNormal::GetRho(vector<float>& Phi, unsigned nBin, bool ARMA) {
  LOG_TRACE();
  // declare all variables that will be used in this function
  vector<float> rhovec(nBin, 0.0);
  if (Phi.size() == 1) {
    LOG_FINEST() <<  "Found single Rho parameter = " <<Phi[0] << " number of bins = " << nBin;
    //calculation of AR(1) acf for  LN2
    for(unsigned i = 1; i <= nBin - 1; i++)
      rhovec[i - 1] = pow(Phi[0],i);
    LOG_FINEST() << "Finished building rho vec";
  } else {
    vector<float> ar, ma,final_acf,Cor;
    vector<vector<float> > A, ind;
    vector<float> psi, theta, Acf;
    // we are doing ARMAacf function
    unsigned p, q, r;
    if (ARMA) {
      q = 1;
      p = 1;
      ar.push_back(Phi[0]);
    } else {
      q = 0;
      p = 2;
      ar = Phi;
    }
    r = fmax(p, q + 1);
    if (p > 0) {
      if (r > 1) {
        if (r > p) {
          LOG_FINEST() << "calculating rho from an ARMA(1,1) process";
          for (unsigned i = 0; i < (r - p); ++i)
            ar.push_back(0.0);
          p = r;
        }
        LOG_FINEST() << "Structureing A";

        A.resize(p + 1, vector<float>(2 * p + 1, 0.0));
        ind.resize(2 * p + 1, vector<float>(p + 1, 0.0));
        for (int i = 0; i < (int)ind.size(); ++i) {
          for (int j = 0; j < (int)ind[i].size(); ++j) {
            ind[i][0] = i + 1;
            ind[i][1] = (i + 1) + (j + 1) - (i + 1);
          }
        }

        for (unsigned i = 0; i < A.size(); ++i) {
          A[i][i] = 1.0;
           A[i][i + 1] = -ar[0];
           A[i][i + 2] = -ar[1];
        }
        LOG_FINEST() << "Populate A. the second ar value" << ar[1];
        ublas::matrix<float> A_eig1(3,3,0.0);
        ublas::matrix<float> A_eig_inv(3,3,0.0);

        vector<float> rhs(3,0.0);
        // initialise rhs, which will be used to solve the following problem, that is Ax = b where b = rhs, so x = A^-1 b
        rhs[0] = 1.0;
        rhs[1] = 0.0;
        rhs[2] = 0.0;
        if (q > 0) {
          LOG_FINEST() << "Calculate ma coef";
          // deal with the ma coeffecient
          psi.push_back(1.0);
          psi.push_back(Phi[0] + Phi[1]);
          theta.push_back(1.0);
          theta.push_back(Phi[1]);
          for (unsigned i = 0; i <= q; ++i)
            theta.push_back(0.0);
          // Calculate rhs
          for (unsigned i = 0; i <= q; ++i) {
            float x1, x2;
            x1 = psi[0] * theta[i];
            x2 = psi[1] * theta[i + q];
            float val = 0.0;
            if (!utilities::To<float, float>(math::Sum({ x1, x2 }), val))
              LOG_CODE_ERROR() << " val " << math::Sum({ x1, x2 }) << " could not be converted in to a float";
            rhs[i] = val;
          }
          rhs[2] = 0.0;
        }
        LOG_FINEST() << "Calculate seq parameter";
        // Use the eigen library yo solve the inverse of for A with known vector B
        //vector<float> Ind;
        vector<unsigned> seq;
        for (unsigned i = 0; i <= p; ++i) {
          seq.push_back(p - i);
        }
        for (unsigned i = 0; i <= p; ++i) {
          for (unsigned j = 0; j <= p; ++j) {
            //LOG_FINEST() << ": i = " << i << " j = " << j << " i index = " << seq[i] << " j index = " << seq[j] << " mat value = " << A[seq[i]][seq[j]];
            float val = 0.0;
            if (j == 2) {
              if (!utilities::To<float, float>(A[i][j], val))
                LOG_CODE_ERROR() << "variable = " << val << " could not be converted in to a float";
              A_eig1(i,j) = val;
            } else {
              if (!utilities::To<float, float>(A[i][j] + A[i][2 * p  - j], val))
                LOG_CODE_ERROR() << "variable = " << val << " could not be converted in to a float";
              A_eig1(i,j) = val;
            }
          }
        }
        for (unsigned i = 0; i <= p; ++i) {
          for (unsigned j = 0; j <= p ; ++j) {
            A_eig1(i,j)=  A_eig1(seq[i],seq[j]);
          }
        }
        // the bodge
        A_eig1(1,2) = 0.0;


        LOG_FINEST() << "Check A mat that we are inverting\n" << A_eig1;

        // Now try mine
        bool inverted = InvertMatrix(A_eig1,A_eig_inv);
        if (!inverted) {
          LOG_FATAL() << "could not invert convariance matrix matrix, if it is a user supplied covariance, check the matrix is invertable, else it is a code error";
        }
        // Matrix multiplication to solve for x
        vector<float> result(A_eig1.size1(),0.0);

        for (unsigned i = 0; i < A_eig_inv.size1(); i++) {
         for (unsigned k = 0; k < A_eig_inv.size2(); k++) {
            result[i] += rhs[k] * A_eig_inv(i,k);
         }
        }
        LOG_FINEST() << "solution = ";
        for(auto num : result)
          LOG_FINEST() << num;

        for (unsigned i = 1; i <= 2; ++i) {
          final_acf.push_back(result[i] / result[0]);
        }
        LOG_FINEST() << "Final Acf";
        for (auto num : final_acf)
          LOG_FINEST() << num << " ";

        Cor = RecursiveFilter(ar, nBin, final_acf);

        // Add the initial coeffiecients to the beginning of the series.
        Cor.insert(Cor.begin(), final_acf[1]);
        Cor.insert(Cor.begin(), final_acf[0]);
        // Print results to screen
        vector<float>::const_iterator first = Cor.begin();
        vector<float>::const_iterator last = Cor.begin() + nBin;
        vector<float> corvec(first, last);
        rhovec = corvec;
        LOG_FINEST() << " Print rhovec";
        for (auto num : rhovec)
          LOG_FINEST()  << num << " ";

      }

    }
  }
return rhovec;
}

vector<float> LogisticNormal::RecursiveFilter(vector<float>& ar_coef, unsigned nBins, vector<float>& initial_vals) {
  LOG_TRACE();
  vector<float> store_vec(nBins + 1,0.0);

  if (ar_coef.size() > 2) {
    LOG_CODE_ERROR() <<  "RecursiveFilter(): has not been coded for more than 2 AR coeffiecients, ar_coef.size() > 2" << endl;
  }
  store_vec[0] = initial_vals[1];
  store_vec[1] = initial_vals[0];
  for (unsigned i = 1; i < nBins + 1; ++i) {
    if (i == 1) {
      store_vec[i] =   store_vec[i - 1] *ar_coef[0]  + store_vec[i] *  ar_coef[1];
    } else {
      store_vec[i] = store_vec[i - 1] *  ar_coef[0] + store_vec[i - 2] * ar_coef[1];
    }
    LOG_FINEST() << "value = " << store_vec[i];
  }
  // remove the first value
  store_vec.erase(store_vec.begin());
  return store_vec;
}

/**
 * Perform cholesky decomposition on our covariance
 * matrix before it's used in the simulations aspect
 * when simualting a multivariate normal distribution
 *
 * @return true on success, false on failure
 */

bool LogisticNormal::DoCholeskyDecmposition() {
  if (covariance_matrix_.size1() != covariance_matrix_.size2() )
      LOG_ERROR() << "Invalid covariance matrix (size1!=size2)";
    unsigned matrix_size1 = covariance_matrix_.size1();
    covariance_matrix_lt = covariance_matrix_;

    for (unsigned i = 0; i < matrix_size1; ++i) {
      for (unsigned j = 0; j < matrix_size1; ++j) {
        covariance_matrix_lt(i,j) = 0.0;
      }
    }

    for (unsigned i = 0; i < matrix_size1; ++i) {
      covariance_matrix_lt(i,i) = 1.0;
    }

    if (covariance_matrix_(0,0) < 0)
      return false;
    float sum = 0.0;

    covariance_matrix_lt(0,0) = sqrt(covariance_matrix_(0,0));

    for (unsigned i = 1; i < matrix_size1; ++i)
      covariance_matrix_lt(i,0) = covariance_matrix_(i,0)/covariance_matrix_lt(0,0);

    for (unsigned i = 1; i < matrix_size1; ++i) {
      sum = 0.0;
      for (unsigned j = 0; j < i; ++j)
        sum += covariance_matrix_lt(i,j) * covariance_matrix_lt(i,j);

      if (covariance_matrix_(i,i) <= sum)
        return false;
      covariance_matrix_lt(i,i) = sqrt(covariance_matrix_(i,i) - sum);
      for (unsigned j = i+1; j < matrix_size1; ++j) {
        sum = 0.0;
        for (unsigned k = 0; k < i; ++k)
          sum += covariance_matrix_lt(j,k) * covariance_matrix_lt(i,k);
        covariance_matrix_lt(j,i) = (covariance_matrix_(j,i) - sum) / covariance_matrix_lt(i,i);
      }
    }
    sum = 0.0;
    for (unsigned i = 0; i < (matrix_size1 - 1); ++i)
      sum += covariance_matrix_lt(matrix_size1 - 1,i) * covariance_matrix_lt(matrix_size1-1,i);
    if (covariance_matrix_(matrix_size1 - 1, matrix_size1 - 1) <= sum)
      return false;
    covariance_matrix_lt(matrix_size1 - 1, matrix_size1 - 1) = sqrt(covariance_matrix_(matrix_size1 - 1, matrix_size1 - 1) - sum);

   return true;
}

/*
 * Invert a square boost matrix
*/
bool LogisticNormal::InvertMatrix(const ublas::matrix<float>& input, ublas::matrix<float>& inverse) {
  typedef ublas::permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  ublas::matrix<float> A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = ublas::lu_factorize(A, pm);
  if (res != 0)
    return false;

  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<float> (A.size1()));

  // backsubstitute to get the inverse
  ublas::lu_substitute(A, pm, inverse);
  return true;
}

/*
 * Calcuate the determinant of a square boost matrix
*/
float LogisticNormal::det_fast(const ublas::matrix<float>& matrix) {
  // create a working copy of the input
 ublas::matrix<float> mLu(matrix);
 ublas::permutation_matrix<size_t> pivots(matrix.size1());

 auto isSingular = ublas::lu_factorize(mLu, pivots);

   if (isSingular)
       return static_cast<float>(0);

   float det = static_cast<float>(1);
   for (size_t i = 0; i < pivots.size(); ++i)
   {
       if (pivots(i) != i)
           det *= static_cast<float>(-1);

       det *= mLu(i, i);
   }
   return det;
}


} /* namespace likelihoods */
} /* namespace niwa */
