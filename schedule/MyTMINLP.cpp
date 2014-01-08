// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  03/17/2006
#include "MyTMINLP.hpp"
#include <coin/BonAmplInterface.hpp>
#include <iostream>

using std::cout;
using std::endl;

MyTMINLP::MyTMINLP(Index n, Index m, Index nnz_jac_g, Index nnz_h_lag,
    const Eigen::MatrixXd& signma, const Eigen::MatrixXd constraints, const ClusterStructure* ptrCS):
  m_ptrCS(ptrCS),printSol_(true)
{
  m_numVariables = n;
  m_numConstraints = m;
  m_numNz_jac_g = nnz_jac_g;
  m_numNz_h_lag = nnz_h_lag;
  m_Constriants = constraints;
  m_Signma = signma;
  m_matI = Eigen::MatrixXd::Identity(n,n);

}

bool 
MyTMINLP::get_variables_types(Index n, VariableType* var_types)
{
//  cout << "=============================================== here get_variables_types ======================="<<endl;
  for (int i = 0; i < m_numVariables; ++i) {
    var_types[i] = BINARY;
  }
  return true;
}


bool 
MyTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
//  cout << "=============================================== here get_variables_linearity ======================="<<endl;
  for (int i = 0; i < m_numVariables; ++i) {
    var_types[i] = Ipopt::TNLP::LINEAR;
  }
  return true;
}


bool 
MyTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
  cout << "=============================================== here get_constraints_linearity ======================="<<endl;
  assert (m==m_numConstraints);
  for (int i = 0; i < m_numConstraints; ++i) {
    const_types[i] = Ipopt::TNLP::LINEAR;
  }
  return true;
}
bool 
MyTMINLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
                       Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
  cout << "=============================================== here get_nlp_info ======================="<<endl;
  n = m_numVariables;//number of variable
  m = m_numConstraints;//number of constraints
  nnz_jac_g = m_numNz_jac_g;//number of non zeroes in Jacobian
  nnz_h_lag = m_numNz_h_lag;//number of non zeroes in Hessian of Lagrangean
  index_style = TNLP::FORTRAN_STYLE;
  return true;
}

bool 
MyTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
  cout << "=============================================== here get_bounds_info ======================="<<endl;
  for (int i = 0; i < m_numVariables; ++i) {
    x_l[i] = 0.;
    x_u[i] = 1.;
  }

  for (int i = 0; i < m_ptrCS->GetNumHeads(); ++i) {
    g_l[i] = -DBL_MAX;
    g_u[i] = 0.;
  }
  for (int i = m_ptrCS->GetNumHeads(); i < m_numConstraints; ++i) {
    g_l[i] = -DBL_MAX;
    g_u[i] = 1.0;
  }
  return true;
}

bool 
MyTMINLP::get_starting_point(Index n, bool init_x, Number* x,
                             bool init_z, Number* z_L, Number* z_U,
                             Index m, bool init_lambda,
                             Number* lambda)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
  cout << "=============================================== here get_starting_point ======================="<<endl;

  
  assert(init_x);
  assert(!init_lambda);
  for (int i = 0; i < m_numVariables; ++i) {
    x[i] = 0.;
  }
  return true;
}

bool 
MyTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n==m_numVariables);
  Eigen::MatrixXd matX = Eigen::MatrixXd::Zero(m_ptrCS->GetNumNodes(),m_ptrCS->GetNumHeads());
  for (int i = 0; i < m_ptrCS->GetNumNodes(); ++i) {
    for (int j = 0; j < m_ptrCS->GetNumHeads(); ++j) {
      if (m_ptrCS->GetChIdxByName(i) == j) {
        matX(i,j) = x[i];
      }
    }
  }
  Eigen::MatrixXd TRM((matX.transpose() * m_Signma * matX).llt().matrixL());
  double sum = 0.0;
//  cout << "Det: " << (matX.transpose() * m_Signma * matX).determinant() << endl;
  Eigen::MatrixXd Diagonal(TRM.diagonal());
//  cout << Diagonal << endl;
  for (int i = 0; i < TRM.diagonal().rows(); ++i) {
    sum += 2.0*log2(Diagonal(i));
  }
  cout << sum << endl;
  obj_value = 1.*sum;

  return true;
}

bool
MyTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n==m_numVariables);

  Eigen::MatrixXd matX = Eigen::MatrixXd::Zero(m_ptrCS->GetNumNodes(),m_ptrCS->GetNumHeads());
  for (int i = 0; i < m_ptrCS->GetNumNodes(); ++i) {
    for (int j = 0; j < m_ptrCS->GetNumHeads(); ++j) {
      if (m_ptrCS->GetChIdxByName(i) == j) {
        matX(i,j) = x[i];
      }
    }
  }
  for (int i = 0; i < m_numVariables; ++i) {
//    grad_f[i] = 1.*m_Signma(i,i);
    grad_f[i] = x[i];
  }
  return true;
}

bool
MyTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
  Eigen::MatrixXd matX = Eigen::MatrixXd::Zero(m_ptrCS->GetNumNodes(),m_ptrCS->GetNumHeads());
  for (int i = 0; i < m_ptrCS->GetNumNodes(); ++i) {
    for (int j = 0; j < m_ptrCS->GetNumHeads(); ++j) {
      if (m_ptrCS->GetChIdxByName(i) == j) {
        matX(i,j) = x[i];
      }
    }
  }
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(m_ptrCS->GetNumHeads(), 1);
  Eigen::MatrixXd result(m_Constriants * matX * matOnes);
  //cout <<  m_Constriants * matX << endl << endl;
  for (int i = 0; i < m_ptrCS->GetNumHeads(); ++i) {
    g[i] = result(i);
//    cout << g[i] << ' ';
  }
//  cout << endl;
  for (int i = m_ptrCS->GetNumHeads(); i < m_numConstraints; ++i) {
    g[i] = matX.col(i-m_ptrCS->GetNumHeads()).sum();
    cout << matX.col(i-m_ptrCS->GetNumHeads()).sum() << endl; 
    cout << matX << endl;
  }
  return true;
}

bool
MyTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
                     Index m, Index nnz_jac, Index* iRow, Index *jCol,
                     Number* values)
{
  assert(n==m_numVariables);
  assert(nnz_jac == m_numNz_jac_g);
//  cout << "=============================================== here eval_jac_g ======================="<<endl;
  if(values == NULL) {
    for (int i = 0; i < m_numVariables; ++i) {
      for (int j = 0; j < m_numConstraints; ++j) {
        iRow[i*m_numConstraints+j] = j+1; 
        jCol[i*m_numConstraints+j] = i+1; 
//        cout << '(' << i*m_numConstraints+j <<"):" << iRow[i*m_numConstraints+j] << ',' <<jCol[i*m_numConstraints+j] << ' ';
      }
//      cout << endl;
    }
    return true;
  }
  else {
    for (int i = 0; i < m_numVariables; ++i) {
      for (int j = 0; j < m_ptrCS->GetNumHeads() ; ++j) {
        values[i*m_ptrCS->GetNumHeads()+j] = m_Constriants(j,i);
      }
    }
    for (int i = 0; i < m_numVariables; ++i) {
      for (int j = m_ptrCS->GetNumHeads(); j < m_numConstraints; ++j) {
        values[m_ptrCS->GetNumHeads()+i*m_ptrCS->GetNumHeads()+j] = 1.;
      }
    }
    return true;
  }
}

bool
MyTMINLP::eval_h(Index n, const Number* x, bool new_x,
                 Number obj_factor, Index m, const Number* lambda,
                 bool new_lambda, Index nele_hess, Index* iRow,
                 Index* jCol, Number* values)
{
  assert (n==m_numVariables);
  assert (m==m_numConstraints);
  assert(nele_hess==m_numNz_h_lag);
  return true;
}

void
MyTMINLP::finalize_solution(TMINLP::SolverReturn status,
                            Index n, const Number* x, Number obj_value)
{
  std::cout<<"Problem status: "<<status<<std::endl;
  std::cout<<"Objective value: "<<obj_value<<std::endl;
  if(printSol_ && x != NULL){
    std::cout<<"Solution:"<<std::endl;
    for(int i = 0 ; i < n ; i++){
      std::cout<<"x["<<i<<"] = "<<x[i];
      if(i < n-1) std::cout<<", ";}
    std::cout<<std::endl;
  }
}
