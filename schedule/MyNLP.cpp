// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyNLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#include "MyNLP.h"

#include <cassert>
#include <cfloat>
#include <cmath>

using namespace Ipopt;

/* Constructor. */

MyNLP::MyNLP(Index n, Index m, Index nnz_jac_g, Index nnz_h_lag, 
    Map const * const ptrMap,
    ULCS1b const * const cSystem,  
    CORRE_MA_OPE const * const ptrGField,  
    double tier1TxTime):
  m_numVariables(n),
  m_numConstraints(m),
  m_numNz_jac_g(nnz_jac_g),
  m_numNz_h_lag(nnz_h_lag),
  m_index_style(FORTRAN_STYLE),
  m_tier1TxTime(tier1TxTime),
  m_cSystem(cSystem),
  m_ptrGField(ptrGField),
  m_ptrMap(ptrMap)
{
  /* construct the vector of head index and variable index*/
  for (int i = 0; i < m_ptrMap->GetNumInitHeads(); ++i) {
    if (m_cSystem->vecHeadName.at(i) >= 0) {
      m_vecHeadTable.push_back(m_cSystem->vecHeadName.at(i)); 
    }
  }
  assert ( m_vecHeadTable.size() == m_numVariables);

  /* construct the vector of cluster entropy */
  list<list<int> >::const_iterator iterRow = m_cSystem->listCluMember->begin();
  vector<int> tmpIndicator(m_ptrMap->GetNumNodes(), 0);
  vector<double> tmpVariance(m_ptrMap->GetNumNodes(), 0.0);

  for (; iterRow != m_cSystem->listCluMember->end(); ++iterRow) {
    list<int>::const_iterator iterCol = iterRow->begin();
    if ( (*iterCol) >= 0) {
      for (; iterCol != iterRow->end(); ++iterCol) {
        tmpIndicator.at(*iterCol) = 1;
      }
    }
    m_vecClusterEntropy.push_back(m_ptrGField->GetJointEntropy(tmpIndicator, tmpVariance, 0, m_ptrMap->GetQBits()));
  }

  assert ( m_vecClusterEntropy.size() == m_numVariables);
}

MyNLP::~MyNLP()
{}

double
MyNLP::GetMinimalPower()
{
  double totalPower = 0.0;
  for (int i = 0; i < m_vecHeadTable.size(); ++i) {
    totalPower += (pow(2.0, m_vecRate.at(i)) - 1.0) * m_ptrMap->GetNoise() / m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i));
  }
  return totalPower;

}

bool MyNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{

 cout << "=============================================== here get_nlp_info ======================="<<endl;
  n = m_numVariables;

  m = m_numConstraints;

  nnz_jac_g = m_numNz_jac_g;

  nnz_h_lag = m_numNz_h_lag;

  index_style = m_index_style;

  return true;
}

bool MyNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
 cout << "=============================================== here get_bounds_info ======================="<<endl;
  assert(n == m_numVariables);
  assert(m == m_numConstraints);

  // x1 has a lower bound of -1 and an upper bound of 1
  for (int i = 0; i < m_numVariables; ++i) {
    x_l[i] = 0.0;
    x_u[i] = DBL_MAX;
  }


  // we have one equality constraint, so we set the bounds on this constraint
  // to be equal (and zero).
  g_l[0] = 0.0;
  g_u[0] = m_tier1TxTime;

  return true;
}

bool MyNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // we initialize x in bounds, in the upper right quadrant
  for (int i = 0; i < m_numVariables; ++i) {
    x[i] = 0.0;
  }

  return true;
}

bool MyNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  // return the value of the objective function
  assert(n == m_numVariables);
  obj_value = 0.0;
  for (int i = 0; i < m_numVariables; ++i) {
    obj_value += (pow(2.0, x[i]) - 1.0 ) * m_ptrMap->GetNoise() / 
      m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i)); 
  }
  return true;
}

bool MyNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)

  assert(n == m_numVariables);
  for (int i = 0; i < m_numVariables; ++i) {
    grad_f[i] = x[i] * pow(2.0, x[i]) * m_ptrMap->GetNoise() /
     m_ptrMap->GetGi0ByNode(m_vecHeadTable.at(i)); 
  }

  return true;
}

bool MyNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // return the value of the constraints: g(x)

  g[0] = 0.0;
  for (int i = 0; i < m_numVariables; ++i) {
    g[0] += m_vecClusterEntropy.at(i)/x[i];
  }

  return true;
}

bool MyNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  assert(n == m_numVariables);
  assert(m == m_numConstraints);
  if (values == NULL) {
    for (int i = 0; i < m_numVariables; ++i) {
      iRow[i] = 1;
      jCol[i] = i + 1;
    }
  }
  else {
    // return the values of the jacobian of the constraints
    Number x1 = x[0];

    // element at 1,1: grad_{x1} g_{1}(x)
    values[0] = -2.0 * x1;

    // element at 1,2: grad_{x1} g_{1}(x)
    values[1] = -1.0;
  }

  return true;
}

bool MyNLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    iRow[0] = 1;
    jCol[0] = 1;

    // element at 2,2: grad^2_{x2,x2} L(x,lambda)
    iRow[1] = 2;
    jCol[1] = 2;

    // Note: off-diagonal elements are zero for this problem
  }
  else {
    // return the values

    // element at 1,1: grad^2_{x1,x1} L(x,lambda)
    values[0] = -2.0 * lambda[0];

    // element at 2,2: grad^2_{x2,x2} L(x,lambda)
    values[1] = -2.0 * obj_factor;

    // Note: off-diagonal elements are zero for this problem
  }

  return true;
}

void MyNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
			      const IpoptData* ip_data,
			      IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution. Since the solution is displayed to the console,
  // we currently do nothing here.
  cout << "Status: " << status << endl;
  cout << "Obj_value: " << obj_value << endl;
}
