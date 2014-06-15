// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: MyTier1NLP.hpp 1861 2010-12-21 21:34:47Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#ifndef __MYNLP_HPP__
#define __MYNLP_HPP__

#include <coin/IpTNLP.hpp>
#include "clusterStructure.h"
#include "map.h"
#include "map.h"
#include "ULCS1b.h"
#include "CORRE_MA_OPE.h"

using namespace Ipopt;

/** C++ Example NLP for interfacing a problem with IPOPT.
 *  MyTier1NLP implements a C++ example showing how to interface with IPOPT
 *  through the TNLP interface. This example is designed to go along with
 *  the tutorial document (see Examples/CppTutorial/).
 *  This class implements the following NLP.
 *
 * min_x f(x) = -(x2-2)^2
 *  s.t.
 *       0 = x1^2 + x2 - 1
 *       -1 <= x1 <= 1
 *
 */
class MyTier1NLP : public TNLP
{
public:
  MyTier1NLP(Index n, Index m, Index nnz_jac_g, Index nnz_h_lag,
    Map const * const ptrMap,
    ULCS1b const * const cSystem,  
    CORRE_MA_OPE const * const ptrGField,  
    double tier1TxTime
    );

  /** default destructor */
  virtual ~MyTier1NLP();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  MyTier1NLP();
  MyTier1NLP(const MyTier1NLP&);
  MyTier1NLP& operator=(const MyTier1NLP&);
  //@}
  bool                        printSol_;
  int                         m_numNodes;
  int                         m_numHeads;
  double                      m_maxPower;

  Index                       m_numVariables;
  Index                       m_numConstraints;
  Index                       m_numNz_jac_g;
  Index                       m_numNz_h_lag;
  TNLP::IndexStyleEnum        m_index_style;
  Eigen::MatrixXd             m_Signma;
  Eigen::MatrixXd             m_Constriants;
  Map const * const           m_ptrMap;
  ULCS1b const * const        m_cSystem;   // system cluseter structure
  CORRE_MA_OPE const * const  m_ptrGField;

  double                      m_tier1TxTime;
  double                      m_obj;
  std::vector<double>         m_vecRate;
  std::vector<int>            m_vecHeadTable;
  std::vector<double>         m_vecClusterEntropy;
};


#endif
