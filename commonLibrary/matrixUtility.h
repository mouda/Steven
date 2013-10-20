#include <vector>
#include <algorithm>//for sort

using namespace std;

// -------------------------------------------------------------------------- //
// @Description: operation methods
// @Provides: 
// -------------------------------------------------------------------------- //

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

