#include <iostream>
#include <vector>

using namespace std;

std::vector<std::vector<double>> inverse_matrix(std::vector<std::vector<double>> matrix, int dimension)
{
    std::vector<std::vector<double>> inverse_matrix(matrix.size(), std::vector<double> (matrix.size()));
    std::vector<std::vector<double>> lower(matrix.size(), std::vector<double> (matrix.size()));
    std::vector<std::vector<double>> upper(matrix.size(), std::vector<double> (matrix.size()));
    std::vector<std::vector<double>> Z(matrix.size(), std::vector<double> (matrix.size()));
    std::vector<std::vector<double>> I(matrix.size(), std::vector<double> (matrix.size()));

    for (int i=0; i<dimension; i++){
      for (int j=0; j<dimension; j++){
        if (i==j) I[i][i]=1;
        else I[i][j]=0;
        Z[i][j]=0;
      }
    }
  
    int i = 0, j = 0, k = 0;
    for (i = 0; i < dimension; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            if (j < i)
                lower[j][i] = 0;
            else
            {
                lower[j][i] = matrix[j][i];
                for (k = 0; k < i; k++)
                {
                    lower[j][i] = lower[j][i] - lower[j][k] * upper[k][i];
                }
            }
        }
        for (j = 0; j < dimension; j++)
        {
            if (j < i)
                upper[i][j] = 0;
            else if (j == i)
                upper[i][j] = 1;
            else
            {
                upper[i][j] = matrix[i][j] / lower[i][i];
                for (k = 0; k < i; k++)
                {
                    upper[i][j] = upper[i][j] - ((lower[i][k] * upper[k][j]) / lower[i][i]);
                }
            }
        }
    } 
  
    // compute z
    for(int col = 0; col < dimension; col++) {
        for(int row = 0; row < dimension; row++) {
            double sum = 0;
            for(int i = 0; i < dimension; i++) {
                if(i != row) {
                    sum += lower[row][i] * Z[i][col];
                }
            }
                Z[row][col] = (I[row][col] - sum)/lower[row][row];
                }
            }
;
    // compute inverse
    for(int col = 0; col < dimension; col++) {
        for(int row = dimension - 1; row >= 0; row--) {
          double sum = 0;
          for(int i = 0; i < dimension; i++) {
              if(i != row) {
                  sum += upper[row][i] * inverse_matrix[i][col];
              }
          }
            inverse_matrix[row][col] = (Z[row][col] - sum)/ upper[row][row];
        }
    }
  
    return inverse_matrix;
}

int main() {

  double dimension; 
  
    //cout<< "Input dimension of a matrix: ";
    cin>>dimension;
  
    //cout << "Input a matrix: " << endl;

    std::vector<std::vector<double>> mat{{1,2,3},{4,5,6}, {7,8,8}};
  std::vector<std::vector<double>> inv_mat = inverse_matrix(mat, dimension);
    
    cout << "Inverse matrix" << endl;
    for (int i=0; i<dimension; i++){
      for (int j=0; j<dimension; j++){
        std::cout<<inv_mat[i][j]<<" ";
      }
      std::cout<<"\n";
    }

    return 0;
}
