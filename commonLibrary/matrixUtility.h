#include <vector>
#include <algorithm>//for sort
#include <climits>


using namespace std;

// -------------------------------------------------------------------------- //
// @Description: operation methods
// @Provides: 
// -------------------------------------------------------------------------- //

/* @brief    maxtrix products
 * @param    two std::vector matrix lhs, rhs 
 * @retval   std::vector matrix
 */

template <typename T> 
vector<vector<T> > multiply2D(  const vector<vector<T> > &lhs, const vector<vector<T> > &rhs) {
  int lhsRowSize = lhs.size();
  assert(lhsRowSize > 0);
  int lhsColSize = lhs[0].size(); 
  int rhsRowSize = rhs.size();
  assert(rhsRowSize > 0);
  T rhsColSize = rhs[0].size();
  vector<vector<T> > matReturn(lhsRowSize,vector<T>(rhsColSize));

  for (int i = 0; i < lhsRowSize; i++) {
    for (int j = 0; j < rhsColSize; j++) {
      int temp = 0;
      for (int k = 0; k < lhsColSize; k++) {
        temp += lhs[i][k]*rhs[k][j];
      }
      matReturn[i][j] = temp;
    }
  }

  return matReturn;
}

/* @brief    count specified element in the vector
 * @param    
 * @retval   
 */

template <typename T>
int countVectorElements( const vector<T>& rhs, const T& element ) {
  if (rhs.size() == 0) return -1;
  int count = 0;
  for (int i = 0; i < rhs.size(); ++i) {
    if (rhs[i] == element) ++count; 
  }
  return count; 
}

template <typename T>
int vecMinIdxNoConsiderZero( const vector<T> &rhs) {
  if (rhs.size() == 0) return -1;
  vector<T> tempVec(rhs);

  for (int i = 0; i < rhs.size(); i++) {
    if (tempVec[i] == 0) tempVec[i] = INT_MAX;  
  }
  return distance(tempVec.begin(),min_element(tempVec.begin(),tempVec.end()));
}





// -------------------------------------------------------------------------- //
// @Description: display methods
// @Provides: 
// -------------------------------------------------------------------------- //


template <typename T>
void mat2Ddisplay( const vector<vector<T> >&mat ){
  for (int i = 0; i < mat.size(); i++) {
    for (int j = 0; j < mat[0].size(); j++) {
      cout << mat[i][j] << ' ';
    }
    cout << endl;
  }
}

template <typename T>
void vec1DDisplay( const vector<T> &vec ) {
  for (int i = 0; i < vec.size(); i++) {
    cout << vec[i] << ' ';
  }
  cout << endl;
}

