
// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  03/17/2006
#include "tier1NLP.h"
#include <coin/BonAmplInterface.hpp>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
#define KSCALE 1000.0

Tier1NLP::Tier1NLP(Index n, Index m, Index nnz_jac_g, Index nnz_h_lag,
    Map const * const ptrMap,
    ULCS1b const * const cSystem,  
    CORRE_MA_OPE const * const ptrGField,  
    double tier1TxTime
    ):
  m_numVariables(n+1),
  m_numConstraints(m),
  m_numNz_jac_g(nnz_jac_g),
  m_numNz_h_lag(nnz_h_lag),
  m_index_style(TNLP::FORTRAN_STYLE),
  m_tier1TxTime(tier1TxTime),
  m_cSystem(cSystem),
  m_ptrGField(ptrGField),
  m_ptrMap(ptrMap),
  printSol_(true)
{
  /* construct the std::vector of head index and variable index*/
  for (int i = 0; i < m_ptrMap->GetNumInitHeads(); ++i) {
    if (m_cSystem->vecHeadName.at(i) >= 0) {
      m_vecHeadTable.push_back(m_cSystem->vecHeadName.at(i)); 
    }
  }
  assert ( m_vecHeadTable.size() == m_numVariables);

  /* construct the std::vector of cluster entropy */
  std::list<std::list<int> >::const_iterator iterRow = m_cSystem->listCluMember->begin();
  std::vector<int> tmpIndicator(m_ptrMap->GetNumNodes(), 0);
  std::vector<double> tmpVariance(m_ptrMap->GetNumNodes(), 1.0);

  for (; iterRow != m_cSystem->listCluMember->end(); ++iterRow) {
    std::list<int>::const_iterator iterCol = iterRow->begin();
    if ( (*iterCol) >= 0) {
      for (; iterCol != iterRow->end(); ++iterCol) {
        tmpIndicator.at(*iterCol) = 1;
      }
    }
    m_vecClusterEntropy.push_back(m_ptrGField->GetJointEntropy(tmpIndicator, tmpVariance, 0, m_ptrMap->GetQBits()));
    std::fill(tmpIndicator.begin(), tmpIndicator.end(), 0);
  }
  for (int i = 0; i < m_vecClusterEntropy.size(); i++) {
    std::cout << i << ':' <<m_vecClusterEntropy.at(i) << ' ';
  }
  std::cout << endl;

  assert ( m_vecClusterEntropy.size() == m_numVariables);
}

double
Tier1NLP::GetMinimalPower()
{
  double totalPower = 0.0;
  for (int i = 0; i < m_vecHeadTable.size(); ++i) {
    totalPower += (pow(2.0, m_vecRate.at(i)/KSCALE/m_ptrMap->GetBandwidth()) - 1.0) * m_ptrMap->GetNoise() / m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i));
  }
  return totalPower;
}

bool 
Tier1NLP::get_variables_types(Index n, VariableType* var_types)
{
// cout << "=============================================== here get_variables_types ======================="<<endl;
  for (int i = 0; i < m_numVariables-1; ++i) {
    var_types[i] = CONTINUOUS;
  }
  var_types[m_numVariables-1] = BINARY;
  return true;
}


bool 
Tier1NLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
// cout << "=============================================== here get_variables_linearity ======================="<<endl;
  for (int i = 0; i < m_numVariables; ++i) {
    var_types[i] = Ipopt::TNLP::NON_LINEAR;
  }
  return true;
}


bool 
Tier1NLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
// cout << "=============================================== here get_constraints_linearity ======================="<<endl;
  assert (m==m_numConstraints);
  for (int i = 0; i < m_numConstraints; ++i) {
    const_types[i] = Ipopt::TNLP::NON_LINEAR;
  }
  return true;
}
bool 
Tier1NLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
                       Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
// cout << "=============================================== here get_nlp_info ======================="<<endl;
  n = m_numVariables;//number of variable
  m = m_numConstraints;//number of constraints
  nnz_jac_g = m_numNz_jac_g;//number of non zeroes in Jacobian
  nnz_h_lag = m_numNz_h_lag;//number of non zeroes in Hessian of Lagrangean
  index_style = TNLP::FORTRAN_STYLE;
  return true;
}

bool 
Tier1NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
// cout << "=============================================== here get_bounds_info ======================="<<endl;
  for (int i = 0; i < m_numVariables; ++i) {
    x_l[i] = 0.0;
    x_u[i] = DBL_MAX;
  }


  g_l[0] = 0.0;
  g_u[0] = m_tier1TxTime;
  

  return true;
}

bool 
Tier1NLP::get_starting_point(Index n, bool init_x, Number* x,
                             bool init_z, Number* z_L, Number* z_U,
                             Index m, bool init_lambda,
                             Number* lambda)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
// cout << "=============================================== here get_starting_point ======================="<<endl;

  
  assert(init_x);
  assert(!init_lambda);
  for (int i = 0; i < m_numVariables; ++i) {
    x[i] = 0.0;
  }
  return true;
}

bool 
Tier1NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n==m_numVariables);
  obj_value = 0.0;
  for (int i = 0; i < m_numVariables; ++i) {
    obj_value += (pow(2.0, x[i]/ KSCALE/m_ptrMap->GetBandwidth()) - 1.0 ) * m_ptrMap->GetNoise() / 
      m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i)); 
  }
  return true;
}

bool
Tier1NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n==m_numVariables);

// cout << "=============================================== here eval_grad_f ======================="<<endl;
  for (int i = 0; i < m_numVariables; ++i) {
    grad_f[i] = x[i] * pow(2.0, x[i]/KSCALE/m_ptrMap->GetBandwidth()) * m_ptrMap->GetNoise() /
     m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i)) / m_ptrMap->GetBandwidth()/KSCALE; 
    cout << m_ptrMap->GetNoise() / m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i)) << ' ' ;
  }
  cout << endl;
  return true;
}

bool
Tier1NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
// cout << "=============================================== here eval_g ======================="<<endl;
  g[0] = 0.0;
  for (int i = 0; i < m_numVariables; ++i) {
    g[0] += m_vecClusterEntropy.at(i)/x[i];
  }
  return true;
}

bool
Tier1NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                     Index m, Index nnz_jac, Index* iRow, Index *jCol,
                     Number* values)
{
  assert(n==m_numVariables);
  assert(nnz_jac == m_numNz_jac_g);
  assert(nnz_jac == m_numVariables);
// cout << "=============================================== here eval_jac_g ======================="<<endl;
  /* zero is included */
  if (values == NULL) {
    for (int i = 0; i < m_numVariables; ++i) {
      iRow[i] = 1;
      jCol[i] = i + 1;
    }
  }
  else {
    for (int i = 0; i < m_numVariables; ++i) {
      values[i] = -1.0 * m_vecClusterEntropy.at(i) * pow(x[i], -2.0); 
    }
  }
}

bool
Tier1NLP::eval_h(Index n, const Number* x, bool new_x,
                 Number obj_factor, Index m, const Number* lambda,
                 bool new_lambda, Index nele_hess, Index* iRow,
                 Index* jCol, Number* values)
{
  assert (n==m_numVariables);
  assert (m==m_numConstraints);
  assert(nele_hess==m_numNz_h_lag);
  assert(nele_hess== m_numVariables);
  if (values == NULL) {
    for (Index i = 0; i < m_numVariables; ++i) {
      iRow[i] = i+1;
      jCol[i] = i+1;
    }
  }
  else {
    for (Index i = 0; i < m_numVariables; ++i) {
      values[i] = (pow(2.0, x[i]/m_ptrMap->GetBandwidth())/KSCALE + pow(x[i], 2.0) * pow (2.0, x[i]/m_ptrMap->GetBandwidth()/KSCALE)/ 
          m_ptrMap->GetBandwidth()/KSCALE)  * m_ptrMap->GetNoise() /
        m_ptrMap->GetBandwidth()/KSCALE/m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i));
    }
    for (Index i = 0; i < m_numVariables; ++i) {
      values[i] += lambda[0] * 2 * m_vecClusterEntropy.at(i) * pow(x[i], -3.0); 
    }

  }
  return true;
}

void
Tier1NLP::finalize_solution(TMINLP::SolverReturn status,
                            Index n, const Number* x, Number obj_value)
{
#ifdef DEBUG
  std::cout<<"Problem status: "<<status<<std::endl;
  std::cout<<"Objective value: "<<obj_value<<std::endl;
  if(printSol_ && x != NULL){
    std::cout<<"Solution:"<<std::endl;
    for(int i = 0 ; i < n ; i++){
      std::cout<<"x["<<i<<"] = "<<x[i];
      if(i < n-1) std::cout<<", ";}
    std::cout<<std::endl;
  }
#endif

  /* Fill the solution here */
  if(x != NULL) {
    m_vecRate.resize(m_numVariables);
    for (int i = 0; i < m_numVariables; ++i) {
      m_vecRate.at(i) = x[i];
    }
    m_obj = obj_value;
  }
}
