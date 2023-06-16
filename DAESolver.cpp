#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <stdexcept>
#include <functional>
#include <fstream>

/*
  DAE koji rjesavamo je
  y_1 + ni * t * y_2 = q(t)  algebarska
  dot(y_1) + ni * t * dot(y_2) + (1+ ni)y_2 = 0; diferencijalna
  dot(y_1) + ni * y_2 + ni * t * dot(y_2) = dot(q(t))
  */

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




void IspisiMatricu(const std::vector<std::vector<double>> mat) {
    for(std::size_t i = 0; i < mat.size(); i++) {
        for(std::size_t j = 0; j < mat.at(0).size(); j++) {
            std::cout << std::setw(8) << mat.at(i).at(j) << " ";
        }
        std::cout << "\n";
    }
}

std::vector<std::vector<double>>ProizvodMatrica(std::vector<std::vector<double>>a, std::vector<std::vector<double>>b){
  //pretpostavljam da matrice nisu grbave
  if(a.at(0).size() != b.size())throw std::domain_error("Nisu pogodne za mnozenje");

  //inicijalizacija
  std::vector<std::vector<double>>rjesenje(a.size(), std::vector<double>(b.at(0).size()));
  for(int i = 0; i < rjesenje.size(); i++){
    for(int j = 0; j < rjesenje.at(0).size(); j++)rjesenje.at(i).at(j) = 0;
  }

  //racunanje
  for(int i = 0; i < a.size(); i++){
    for(int j = 0; j < b.at(0).size(); j++){
      for(int k = 0; k < b.size(); k++)rjesenje.at(i).at(j)+=a.at(i).at(k) * b.at(k).at(j);
    }
  }
  return rjesenje;
}
std::vector<double>ProizvodMatriceIKolone(std::vector<std::vector<double>>a, std::vector<double>b){
  //pretpostavljam da matrice nisu grbave
  if(a.at(0).size() != b.size())throw std::domain_error("Nisu pogodne za mnozenje");

  //inicijalizacija
  std::vector<double>rjesenje(b.size());
  for(int i = 0; i < rjesenje.size(); i++){
    rjesenje.at(i) = 0;
  }

  //racunanje
  for(int i = 0; i < a.size(); i++){
      for(int k = 0; k < b.size(); k++)rjesenje.at(i)+=a.at(i).at(k) * b.at(k);
    }
  return rjesenje;
}
std::vector<std::vector<double>>MojDAESolverPrvogPrimjera(std::function<double(double)>q, double korak, double trenutak, double ni){
  std::vector<double> y1{};
  std::vector<double> y2{};
  std::vector<std::vector<double>> Y;
  y1.push_back(0);
  double delta_t = korak, t = 0;
  y2.push_back(-1);
  Y.push_back({0, -1});
  t = delta_t;
  for (int i = 0; i < trenutak / korak; i++) {
    /*
    y1 + ni*t*y2 = q(t)
    y1 - y1(t-1) /delta_t + ni *t * (y2 - y2(t-1))/delta_t  + y_2 (1+ni) = 0;    
      */
    
    std::vector<std::vector<double>>A{{1,ni * t},{1/delta_t, ni * t / delta_t + 1 + ni}};
    std::vector<double>b{q(t),Y.at(i).at(0) / delta_t + ni * t * Y.at(i).at(1) / delta_t};
    Y.push_back(ProizvodMatriceIKolone(inverse_matrix(A,A.size()), b));
    t += delta_t;
  }
  return Y;

  

  
}
int main() {
  double ni = 2.0;
  // koristimo q(t) kao sin(t)
  std::vector<double> y1{};
  std::vector<double> y2{};
  std::vector<std::vector<double>> Y;

  std::cout << "Unesite korak: ";
  double korak;
  std::cin >> korak;

  std::cout << "Do kojeg trenutka: ";
  double ts = 0;
  std::cin >> ts;
  Y = MojDAESolverPrvogPrimjera([](double x) -> double {
                       return sin(x);
                   }, korak, ts, ni);
  
  std::cout << "y1 je nakon " << ts << "s = " << Y.at(ts / korak - 1).at(0)
            << "dok je y2 nakon " << ts << "s = " << Y.at(ts / korak - 1).at(1);
  std::cout << std::endl;

  std::cout << "Tacni rezultati su: y1 = " << sin(ts) + ni * cos(ts) * ts;

  std::cout << "i y2 = ";
  std::cout << -cos(ts);
  std::cout << std::endl;
  std::cout << "Dakle napravljenje su greske redom "
            << Y.at(ts / korak - 1).at(0) - sin(ts) - ni * cos(ts) * ts << " i "
            << Y.at(ts / korak - 1).at(1) + cos(ts) << std::endl;
  std::cout << "Relativne greske u postotcima su: "
            << (Y.at(ts / korak - 1).at(0) - sin(ts) - ni * cos(ts) * ts) *
                   100 / (sin(ts) + ni * cos(ts) * ts)
            << "% i "
            << (Y.at(ts / korak - 1).at(1) + cos(ts)) * 100 / (-cos(ts)) << "%"
            << std::endl;
            std::ofstream out("y1.txt");
            std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
            std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(0)<<std::endl;

            std::ofstream out2("y2.txt");
            std::cout.rdbuf(out2.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(1)<<std::endl;



}
