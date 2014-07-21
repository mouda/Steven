// (C) Copyright Carnegie Mellon University 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  03/17/2006
#include "minPowerImageMILP.hpp"
#include <coin/BonAmplInterface.hpp>
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

MinPowerImageMILP::MinPowerImageMILP(Index n, Index m, Index nnz_jac_g, Index nnz_h_lag,
    const double tau,
    const int tier2NumSlot,
    const double bandwidthKhz,    
    const ClusterStructure* ptrCS,
    const Map* ptrMap):
  m_txTimePerSlot(tau), 
  m_tier2NumSlot(tier2NumSlot),
  m_bandwidthKhz(bandwidthKhz),
  m_ptrCS(ptrCS), 
  m_ptrMap(ptrMap), 
  m_vecFi(ptrMap->GetNumNodes()),
  m_Gamma(1.0),
  m_Phi(10000.0),
  printSol_(true)
{
  m_numVariables = n;
  m_numConstraints = m;
  m_numNz_jac_g = nnz_jac_g;
  m_numNz_h_lag = nnz_h_lag;

  m_numNodes = m_ptrMap->GetNumNodes();
  m_numHeads =  m_ptrMap->GetNumInitHeads();
  m_maxPower = m_ptrMap->GetMaxPower();

//  fill(m_vecFi.begin(), m_vecFi.end(), pow(2.0, m_ptrMap->GetIdtEntropy()/m_txTimePerSlot/m_bandwidthKhz) - 1.0);
  for (int i = 0; i < m_ptrMap->GetNumNodes(); ++i) {
    m_vecFi.at(i) = pow(2.0, m_ptrMap->GetIdtEntropy(i)/m_txTimePerSlot/m_bandwidthKhz) - 1.0;
  }
  /* solution */
  
  // -------------------------------------------------------------------------- //
  // @Description: The constriant matrix 
  //
  //          |S| x N (Q_ij)     |S| x N (Y_ij)
  //         /              \   /              \
  //       ______________________________________
  //    / |                   |                  |
  // |S|  |                   |                  |
  //  x   |       Eq. 4.6     |                  |
  //  N   |                   |                  |
  //    \ |                   |                  |
  //      |___________________|__________________|
  //    / |                   |                  |
  // |S|  |       Eq. 4.9     |                  |
  //    \ |___________________|__________________|
  //      |                   |                  |
  //    / |       Eq. 4.10    |                  |
  // |H|  |                   |                  |
  //    \ |___________________|__________________|
  //      |                   |                  |
  //    / |                   |                  |
  // |H|  |       Eq. 4.13    |                  |
  //  x   |                   |                  |
  // |S|  |                   |                  |
  //    \ |___________________|__________________|
  //
  // @Provides: 
  // -------------------------------------------------------------------------- //
  
  /* Construct the constraints */
  m_matConstraints = Eigen::MatrixXd::Zero(m_numConstraints, m_numVariables);
  /* Interfereence constraints eq 4.6 */
  Eigen::IOFormat CleanFmt(2, 0, " ", "\n", "", "");
  for (int n_c = 0; n_c < m_tier2NumSlot; ++n_c) {
    for (int i_c = 0; i_c < m_numNodes; ++i_c) {

      /* row index of constraint matrix */
      int index_i = n_c * m_numNodes + i_c; 
      int headName = m_ptrCS->GetChNameByName(i_c);

      if ( headName < 0 || headName == i_c) continue; /* skip the un-supported node */

      for (int n_q = 0; n_q < m_tier2NumSlot; ++n_q) {
        for (int i_q = 0; i_q < m_numNodes; ++i_q) {

          /* column index of constraint matrix */
          int index_j = n_q * m_numNodes +  i_q;
          if ( n_c == n_q && i_c == i_q ) {
            m_matConstraints(index_i, index_j ) = m_ptrMap->GetGijByPair(i_q, headName);
          }
          else if (n_c == n_q && headName != m_ptrCS->GetChNameByName(i_q) ) {
            m_matConstraints(index_i, index_j) = -1 * m_Gamma * m_vecFi.at(i_q) * m_ptrMap->GetGijByPair(i_q, headName);  
          }
          else {
            m_matConstraints(index_i, index_j) = 0.0;
          }
        }
      }

      for (int n_y = 0; n_y < m_tier2NumSlot; ++n_y) {
          for (int i_y = 0; i_y < m_numNodes; ++i_y) {
          int index_j = m_tier2NumSlot * m_numNodes + n_y * m_numNodes +  i_y;
          if ( n_c == n_y && i_c == i_y ) {
            m_matConstraints(index_i, index_j) = -1.0 * m_Phi;
          }
          else {
            m_matConstraints(index_i, index_j) = 0.0;
          }
        }
      }
    }
  }


  /* eq 4.9 */
  for (int i_c = 0; i_c < m_numNodes; ++i_c) {
    /* row index of constraint matrix */
    int index_i = m_numNodes * m_tier2NumSlot + i_c;
    int headName = m_ptrCS->GetChNameByName(i_c);
    if ( headName < 0 || headName == i_c ) continue;/* skip the un-supported node */

    for (int n_y = 0; n_y < m_tier2NumSlot; ++n_y) {
      for (int i_y = 0; i_y < m_numNodes; ++i_y) {

        /* column index of constraint matrix */
        int index_j = m_tier2NumSlot * m_numNodes + n_y * m_numNodes +  i_y;
        if ( i_c == i_y ) {
          m_matConstraints(index_i, index_j) = 1.0;
        }
        else {
          m_matConstraints(index_i, index_j) = 0.0;
        }
      }
    }
  }

  /* eq 4.10 */
  for (int j_c = 0; j_c < m_numHeads; ++j_c) {
    /* row index of constraint matrix */
    int index_i = m_numNodes * m_tier2NumSlot + m_numNodes + j_c;

    if(m_ptrCS->GetVecHeadName().at(j_c) < 0 ) continue;
    for (int n_y = 0; n_y < m_tier2NumSlot; ++n_y) {
      for (int i_y = 0; i_y < m_numNodes; ++i_y) {

        /* column index of constraint matrix */
        int index_j = m_tier2NumSlot * m_numNodes + n_y * m_numNodes +  i_y;
        if (  j_c  ==  m_ptrCS->GetChIdxByName(i_y)) {
          m_matConstraints(index_i, index_j) = 1.0;
        }
        else {
          m_matConstraints(index_i, index_j) = 0.0;
        }
      }
    }
  }

  /* eq 4.13 */
  for (int n_c = 0; n_c < m_tier2NumSlot; ++n_c) {
    for (int j_c = 0; j_c < m_numHeads; ++j_c) {

      /* row index of constraint matrix */
      if (m_ptrCS->GetVecHeadName().at(j_c) < 0) continue;
      int index_i = m_numNodes * m_tier2NumSlot + m_numNodes + m_numHeads + n_c * m_numHeads + j_c;
      
      for (int n_y = 0; n_y < m_tier2NumSlot; ++n_y) {
        for (int i_y = 0; i_y < m_numNodes; ++i_y) {

          /* column index of constraint matrix */
          int index_j = m_tier2NumSlot * m_numNodes + n_y * m_numNodes +  i_y;
          if ( j_c   ==  m_ptrCS->GetChIdxByName(i_y) && n_c == n_y) {
            m_matConstraints(index_i, index_j) = 1.0;
          }
          else {
            m_matConstraints(index_i, index_j) = 0.0;
          }
        }
      }

    }
  }

//  fstream testFstream;
//  testFstream.open("TEST.out", ios::out);
//  testFstream << m_matConstraints.format(CleanFmt);
//  testFstream.close();
}

bool 
MinPowerImageMILP::get_variables_types(Index n, VariableType* var_types)
{
// cout << "=============================================== here get_variables_types ======================="<<endl;
  for (int i = 0; i < m_tier2NumSlot * m_numNodes; ++i) {
    var_types[i] = CONTINUOUS;
  }
  for (int i = m_tier2NumSlot * m_numNodes; i < m_numVariables; ++i) {
    var_types[i] = BINARY;
  }
  return true;
}


bool 
MinPowerImageMILP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
// cout << "=============================================== here get_variables_linearity ======================="<<endl;
  for (int i = 0; i < m_numVariables; ++i) {
    var_types[i] = Ipopt::TNLP::LINEAR;
  }
  return true;
}


bool 
MinPowerImageMILP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
// cout << "=============================================== here get_constraints_linearity ======================="<<endl;
  assert (m==m_numConstraints);
  for (int i = 0; i < m_numConstraints; ++i) {
    const_types[i] = Ipopt::TNLP::LINEAR;
  }
  return true;
}
bool 
MinPowerImageMILP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
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
MinPowerImageMILP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
// cout << "=============================================== here get_bounds_info ======================="<<endl;
  
  /* Bound of variables */
  for (int n_q = 0; n_q < m_tier2NumSlot; ++n_q) {
    for (int i_q = 0; i_q < m_numNodes; ++i_q) {
      int index_j = n_q * m_numNodes +  i_q; 
      if ( m_ptrCS->GetChIdxByName(i_q) >= 0 && m_ptrCS->GetChNameByName(i_q) != i_q ) {
        x_l[index_j] = 0.000001;
        x_u[index_j] = DBL_MAX;
      }
      else {
        x_l[index_j] = 0.0;
        x_u[index_j] = 0.0;
      }
    }
  }


  for (int n_y = 0; n_y < m_tier2NumSlot; ++n_y) {
    for (int i_y = 0; i_y < m_numNodes; ++i_y) {
      int index_j = m_tier2NumSlot * m_numNodes + n_y * m_numNodes +  i_y; 
      if ( m_ptrCS->GetChIdxByName(i_y) >= 0 && m_ptrCS->GetChNameByName(i_y) != i_y ) {
        x_l[index_j] = 0.0;
        x_u[index_j] = 1.0;
      }
      else {
        x_l[index_j] = 0.0;
        x_u[index_j] = 0.0;
      }
    }
  }

  /* Bound of constraints */
  for (int n_c = 0; n_c < m_tier2NumSlot; ++n_c) {
    for (int i_c = 0; i_c < m_numNodes; ++i_c) {

      /* row index of constraint matrix */
      int index_i = n_c * m_numNodes + i_c; 
      int headName = m_ptrCS->GetChNameByName(i_c);
      if ( headName < 0  || headName == i_c ) { /* skip the un-supported node */
        g_l[index_i] = -DBL_MAX;
        g_u[index_i] = DBL_MAX; 
      } 
      else {
        g_l[index_i] = 0.0;
        g_u[index_i] = DBL_MAX;
      }
    }
  }


  /* eq 4.9 */
  for (int i_c = 0; i_c < m_numNodes; ++i_c) {
    /* row index of constraint matrix */
    int index_i = m_numNodes * m_tier2NumSlot + i_c;
    int headName = m_ptrCS->GetChNameByName(i_c);
    if ( headName < 0 || headName == i_c) { /* skip the un-supported node */
      g_l[index_i] = -DBL_MAX;
      g_u[index_i] = DBL_MAX; 
    }
    else {
      g_l[index_i] = 0.0;
      g_u[index_i] = 1.0;
    }
  }

  /* eq 4.10 */
  list<list<int> >::const_iterator iterRow = m_ptrCS->GetListCluMemeber().begin();
  std::vector<int> vecCluSize(m_ptrMap->GetNumInitHeads(), -1);
  for (int j_c = 0 ; j_c < m_numHeads; ++j_c, ++iterRow) {
    /* row index of constraint matrix */
    int index_i = m_numNodes * m_tier2NumSlot + m_numNodes + j_c;
    if(m_ptrCS->GetVecHeadName().at(j_c) < 0 ) {
      g_l[index_i] = -DBL_MAX;
      g_u[index_i] = DBL_MAX;
    }
    else {
      //g_l[index_i] = static_cast<double>(iterRow->size() - 1); 
      g_l[index_i] = 0; 
      g_u[index_i] = static_cast<double>(iterRow->size() - 1);
      vecCluSize.at(j_c) = iterRow->size() - 1;
    }
  }

  /* eq 4.13 */
  for (int n_c = 0; n_c < m_tier2NumSlot; ++n_c) {
    for (int j_c = 0; j_c < m_numHeads; ++j_c) {

      /* row index of constraint matrix */
      cout << vecCluSize.at(j_c) << ' ';
      int index_i = m_numNodes * m_tier2NumSlot + m_numNodes + m_numHeads + n_c * m_numHeads + j_c;
      if (m_ptrCS->GetVecHeadName().at(j_c) < 0) {
        g_l[index_i] = -DBL_MAX;
        g_u[index_i] = DBL_MAX;
      }
      else {
        if (vecCluSize.at(j_c) > 0 ) {
          g_l[index_i] = 1.0;
          g_u[index_i] = 1.0;
          vecCluSize.at(j_c) = vecCluSize.at(j_c) -1;
        }
        else {
          g_l[index_i] = 0.0;
          g_u[index_i] = 0.0;
        }
      }

    }
    cout << endl;
  }

  return true;
}

bool 
MinPowerImageMILP::get_starting_point(Index n, bool init_x, Number* x,
                             bool init_z, Number* z_L, Number* z_U,
                             Index m, bool init_lambda,
                             Number* lambda)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
// cout << "=============================================== here get_starting_point ======================="<<endl;

  
  assert(init_x);
  assert(!init_lambda);
  for (int i = 0; i < m_tier2NumSlot * m_numNodes; ++i) {
    x[i] = 0.0;
  }
  for (int i = 0; i < m_numVariables; ++i) {
    x[i] = 0.0;
  }
  return true;
}

bool 
MinPowerImageMILP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n==m_numVariables);
  Number powerSum = 0.0;
  for (int i = 0; i < m_tier2NumSlot * m_numNodes; ++i) {
    powerSum += x[i];
  }
  obj_value = powerSum; 
  return true;
}

bool
MinPowerImageMILP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n==m_numVariables);

// cout << "=============================================== here eval_grad_f ======================="<<endl;
  for (int i = 0; i < m_tier2NumSlot * m_numNodes; ++i) {
    grad_f[i] = 1.0;
  }
  for (int i = m_tier2NumSlot * m_numNodes; i < 2 * m_tier2NumSlot * m_numNodes; ++i) {
    grad_f[i] = 0.0;
  }
  return true;
}

bool
MinPowerImageMILP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n==m_numVariables);
  assert(m==m_numConstraints);
// cout << "=============================================== here eval_g ======================="<<endl;
  Eigen::MatrixXd vecVariables = Eigen::MatrixXd::Zero(m_numVariables, 1);
  for (int i = 0; i < m_numVariables; ++i) {
    vecVariables(i,0) = x[i];
  }
  Eigen::MatrixXd vecResult = m_matConstraints * vecVariables;
  for (int n_c = 0; n_c < m_tier2NumSlot ; ++n_c) {
    for (int i_c = 0; i_c < m_numNodes; ++i_c) {
      int index_i = n_c * m_numNodes + i_c;
      g[index_i] = vecResult(index_i,0) + m_Phi - m_ptrMap->GetNoise() * m_vecFi.at(i_c) * m_Gamma;
    }
  }
  for (int index_i = m_tier2NumSlot * m_numNodes; index_i < m_numConstraints; ++index_i) {
    g[index_i] = vecResult(index_i,0);
  }
  return true;
}

bool
MinPowerImageMILP::eval_jac_g(Index n, const Number* x, bool new_x,
                     Index m, Index nnz_jac, Index* iRow, Index *jCol,
                     Number* values)
{
  assert(n==m_numVariables);
  assert(nnz_jac == m*n);
  assert(nnz_jac == m_numNz_jac_g);
// cout << "=============================================== here eval_jac_g ======================="<<endl;
  /* zero is included */
  if(values == NULL) {
    for (int i = 0; i < m_numConstraints; ++i) {
      for (int j = 0; j < m_numVariables; ++j) {
        iRow[i*m_numVariables+j] = i+1; 
        jCol[i*m_numVariables+j] = j+1; 
      }
    }
    return true;
  }
  else {
    for (int i = 0; i < m_numConstraints; ++i) {
      for (int j = 0; j < m_numVariables; ++j) {
        values[i*m_numVariables+j] = m_matConstraints(i,j); 
      }
    }
    return true;
  }
}

bool
MinPowerImageMILP::eval_h(Index n, const Number* x, bool new_x,
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
MinPowerImageMILP::finalize_solution(TMINLP::SolverReturn status,
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

  /* Fill the solution here */
  if(x != NULL) {
    m_vecSolution.resize(m_numVariables);
    for (int i = 0; i < m_numVariables; ++i) {
      m_vecSolution.at(i) = x[i];
    }
  }
}
