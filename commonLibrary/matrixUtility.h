#include <vector>
#include <algorithm>//for sort

using namespace std;

vector<vector<int> > multiply2D(  const vector<vector<int> > &lhs, const vector<vector<int> > &rhs) {
  int lhsRowSize = lhs.size();
  assert(lhsRowSize > 0);
  int lhsColSize = lhs[0].size(); 
  int rhsRowSize = rhs.size();
  assert(rhsRowSize > 0);
  int rhsColSize = rhs[0].size();
  vector<vector<int> > matReturn(lhsRowSize,vector<int>(rhsColSize));

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

template <typename T>
void mat2Ddisplay( const vector<vector<T> >&mat ){
  for (int i = 0; i < mat.size(); i++) {
    for (int j = 0; j < mat[0].size(); j++) {
      cout << mat[i][j] << ' ';
    }
    cout << endl;
  }
}
