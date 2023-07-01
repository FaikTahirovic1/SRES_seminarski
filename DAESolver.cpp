#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <iomanip>
#include <stdexcept>
#include <functional>
#include <fstream>

std::vector<double>Tdop{8.96, 6.0, 5.89};
std::vector<double>Tqop{0.31, 0.535, 0.6};
std::vector<double>Xd{0.146, 0.8958, 1.3125};
std::vector<double>Xdp{0.0608, 0.1198, 0.1813};
std::vector<double>Xq{0.0969, 0.8645, 1.2578};
std::vector<double>Xqp{0.0969, 0.1969, 0.2500};
std::vector<double>RS{0.0041,                 0.0026,                0.0035};
//std::vector<double>RS{0.1,                 0.1,                0.1};
double ws = 2*4*atan(1)*60;
std::vector<double>H{23.640,                 6.4000 ,               3.0100};
//std::vector<double>Dm{0.1*(2*23.64)/ws ,      0.2*(2*6.4)/ws    ,  0.3*(2*3.01)/ws};
std::vector<double>Dm{0.7 ,     0.7    , 0.7};
double TE = 0.3140;
double KE = 1.0;
double Ax = 0.0039;
double Bx = 1.5550;
double TF = 0.3500;
double KF = 0.0630;
double TA = 0.2000;
double KA = 20.000;
std::vector<double>Vg{1.04,1.025,1.025 }; 
std::vector<double>Vref{1.0934, 1.11812, 1.09769};

double TCH = 0.10;
double TSV = 0.05;
double RD = 0.05;

std::vector<double>TH{0.0,  0.161967,   0.0814153};


std::vector<double>PC{0.718633, 1.63659, 0.852446};


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
        }//std::cout<<"\n clan je "<<lower[i][i]<<"\n";
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
          }//if(fabs(upper[row][row])<0.1)std::cout<<"upper je "<<upper[row][row]<<"\n";
            inverse_matrix[row][col] = (Z[row][col] - sum)/ upper[row][row];
            
        }
    }
  
    return inverse_matrix;
}

bool DaLiJeKvadratna(std::vector<std::vector<double>>a){
  for (int i = 0; i < a.size(); i++)if(a.at(i).size() != a.size()) return false;
  return true;
}


void IspisiMatricu(const std::vector<std::vector<double>> mat) {
    for(std::size_t i = 0; i < mat.size(); i++) {
        for(std::size_t j = 0; j < mat.at(0).size(); j++) {
            std::cout << std::setw(8) << mat.at(i).at(j) << " ";
        }
        std::cout << "\nRed\n";
    }
}

void subMatrix(std::vector<std::vector<double>>mat, std::vector<std::vector<double>>&temp, int p, int q, int n) {
   int i = 0, j = 0;
   // filling the sub matrix
   for (int row = 0; row < n; row++) {
      for (int col = 0; col < n; col++) {
         // skipping if the current row or column is not equal to the current
         // element row and column
         std::cout<<"   hel   ";
         if (row != p && col != q) {
            temp[i][j++] = mat[row][col];
            if (j == n - 1) {
               j = 0;
               i++;
            }
         }std::cout<<"ovdje  "<<col;
      }std::cout<<"Da    \n";
   }
}

double determinantOfMatrix(std::vector<std::vector<double>>matrix, int n) {
   double determinant = 0;
   if (n == 1) {
      return matrix[0][0];
   }
   if (n == 2) {
      return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
   }
   std::vector<std::vector<double>> temp(3, std::vector<double>(3)); 
   int sign = 1;
   for (int i = 0; i < n; i++) {
      subMatrix(matrix, temp, 0, i, n);
      std::cout<<"ovdje";
      std::cout<<"derm zasad: "<<determinant<<"\n";
      determinant += sign * matrix[0][i] * determinantOfMatrix(temp, n - 1);
      sign = -sign;
   }
   return determinant;
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
std::vector<std::vector<double>>MojDAESolverPrvogPrimjera(double korak, double trenutak){
  std::vector<std::vector<double>> Y;
  double delta_t = korak, t = 0;
  std::vector<double>ypom;
  //init values gen1
  double Eqp10 = 1.05914;
  ypom.push_back(Eqp10);
  double Edp10 = 0;
  ypom.push_back(Edp10);
  double delta10 = 0.0614231;
  ypom.push_back(delta10);
  double w10 = 376.991;
  ypom.push_back(w10);
  double Efd10 = 1.08486;
  ypom.push_back(Efd10);
  double RF10 = 0.195275;
  ypom.push_back(RF10);
  double VR10 = 1.10772;
  ypom.push_back(VR10);
  double Id10 = 0.301852;
  ypom.push_back(Id10);
  double Iq10 = 0.671593;
  ypom.push_back(Iq10);
  double TM10 =  0.718633;
  ypom.push_back(TM10);
  double PSV10 = 0.718633;
  ypom.push_back(PSV10);


    //init values gen2
  double Eqp20 = 0.791927;
   ypom.push_back(Eqp20);
  double Edp20 = 0.623846;
  ypom.push_back(Edp20);
  double delta20 = 1.06445;
  ypom.push_back(delta20);
  double w20 = 376.991;
  ypom.push_back(w20);
  double Efd20 = 1.79169;
  ypom.push_back(Efd20);
  double RF20 = 0.322505;
  ypom.push_back(RF20);
  double VR20 = 1.90501;
  ypom.push_back(VR20);
  double Id20 = 1.28836;
  ypom.push_back(Id20);
  double Iq20 = 0.934461;
  ypom.push_back(Iq20);
  double TM20 =  1.63659;
  ypom.push_back(TM20);
  double PSV20 = 1.63659;
  ypom.push_back(PSV20);


      //init values gen3
  double Eqp30 = 0.770984;
   ypom.push_back(Eqp30);
  double Edp30 = 0.625046;
  ypom.push_back(Edp30);
  double delta30 = 0.943432;
  ypom.push_back(delta30);
  double w30 = 376.991;
  ypom.push_back(w30);
  double Efd30 = 1.40512;
  ypom.push_back(Efd30);
  double RF30 = 0.252921;
  ypom.push_back(RF30);
  double VR30 = 1.45383;
  ypom.push_back(VR30);
  double Id30 = 0.560582;
  ypom.push_back(Id30);
  double Iq30 = 0.620208;
  ypom.push_back(Iq30);
  double TM30 =  0.852446;
  ypom.push_back(TM30);
  double PSV30 = 0.852446;
  ypom.push_back(PSV30);
  Y.push_back(ypom);
  t = delta_t;
  for (int i = 0; i < trenutak / korak; i++) {


    if(i > 1/delta_t && i < 2/delta_t)PC.at(0) =0.718633 * 0.9;
    else PC.at(0) = 0.718633;



    std::vector<std::vector<double>>A{{Tdop.at(0)/delta_t + 1 , 0 , 0, 0, -1, 0, 0, Xd.at(0) - Xdp.at(0), 0, 0, 0,    0, 0, 0,0,0,0,0,0,0,0,0,       0, 0, 0,0,0,0,0,0,0,0,0},
    {0, Tqop.at(0)/delta_t + 1, 0, 0, 0, 0, 0, 0, Xq.at(0) - Xqp.at(0), 0 ,0,        0, 0, 0,0,0,0,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},
    {0,0,1/delta_t,-1,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,        0, 0, 0,0,0,0,0,0,0,0,0},
    {Y.at(i).at(8), Y.at(i).at(7),0, 2*H.at(0)/(ws * delta_t) + Dm.at(0) , 0,0,0, (Xqp.at(0) - Xdp.at(0)) * Y.at(i).at(8), 0,-1,0,    0,0,0,0,0,0,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},//prob
    {0,0,0,0, TE/delta_t + KE + Ax * exp(Bx * Y.at(i).at(4)), 0, -1, 0, 0, 0, 0,         0,0,0,0,0,0,0,0,0,0,0,        0,0,0,0,0,0,0,0,0,0,0,},//prob
    {0,0,0,0,-KF/TF,TF/delta_t + 1,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,KA * KF /TF,KA,TA/delta_t + 1, 0, 0, 0, 0,    0,0,0,0,0,0,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},
    
    

    {0, 1, 0,0,0,0,0,-RS.at(0), Xqp.at(0), 0, 0,     0,0,0,0,0,0,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},
    {1, 0, 0,0,0,0,0,Xdp.at(0), - RS.at(0), 0, 0,      0,0,0,0,0,0,0,0,0,0,0,        0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,1/(RD * ws),0,0,0,0,0,0,TSV/delta_t + 1,     0,0,0,0,0,0,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,TCH/delta_t + 1, -1,    0,0,0,0,0,0,0,0,0,0,0,       0,0,0,0,0,0,0,0,0,0,0},
    
    //gen2
    
    { 0, 0, 0,0,0,0,0,0,0,0,0,      Tdop.at(1)/delta_t + 1 , 0 , 0, 0, -1, 0, 0, Xd.at(1) - Xdp.at(1), 0, 0, 0,      0, 0, 0,0,0,0,0,0,0,0,0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, Tqop.at(1)/delta_t + 1, 0, 0, 0, 0, 0, 0, Xq.at(1) - Xqp.at(1), 0 ,0,         0, 0, 0,0,0,0,0,0,0,0,0},
    {0, 0, 0,0,0,0,0,0,0,0,0,     0,0,1/delta_t,-1,0,0,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,        Y.at(i).at(19), Y.at(i).at(18),0, 2*H.at(1)/(ws * delta_t) + Dm.at(1) , 0,0,0, (Xqp.at(1) - Xdp.at(1)) * Y.at(i).at(19), 0,-1,0,      0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,          0,0,0,0, TE/delta_t + KE + Ax * exp(Bx * Y.at(i).at(15)), 0, -1, 0, 0, 0, 0,       0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,           0,0,0,0,-KF/TF,TF/delta_t + 1,0,0,0,0,0,      0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,          0,0,0,0,KA * KF /TF,KA,TA/delta_t + 1, 0, 0, 0, 0,      0,0,0,0,0,0,0,0,0,0,0},
    
    
    {0,0,0,0,0,0,0,0,0,0,0,     0, 1, 0,0,0,0,0,-RS.at(1), Xqp.at(1), 0, 0,        0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,      1, 0, 0,0,0,0,0,Xdp.at(1), - RS.at(1), 0, 0,      0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,         0,0,0,1/(RD * ws),0,0,0,0,0,0,TSV/delta_t + 1,     0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,        0,0,0,0,0,0,0,0,0,TCH/delta_t + 1, -1,       0,0,0,0,0,0,0,0,0,0,0},
    
    
    //gen3
    
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     Tdop.at(2)/delta_t + 1 , 0 , 0, 0, -1, 0, 0, Xd.at(2) - Xdp.at(2), 0, 0, 0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     0, Tqop.at(2)/delta_t + 1, 0, 0, 0, 0, 0, 0, Xq.at(2) - Xqp.at(2), 0 ,0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     0,0,1/delta_t,-1,0,0,0,0,0,0,0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     Y.at(i).at(30), Y.at(i).at(29),0, 2*H.at(2)/(ws * delta_t) + Dm.at(2) , 0,0,0, (Xqp.at(2) - Xdp.at(2)) * Y.at(i).at(30), 0,-1,0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     0,0,0,0, TE/delta_t + KE + Ax * exp(Bx * Y.at(i).at(26)), 0, -1, 0, 0, 0, 0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,      0,0,0,0,-KF/TF,TF/delta_t + 1,0,0,0,0,0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     0,0,0,0,KA * KF /TF,KA,TA/delta_t + 1, 0, 0, 0, 0},
   
    
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     0, 1, 0,0,0,0,0,-RS.at(2), Xqp.at(2), 0, 0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     1, 0, 0,0,0,0,0,Xdp.at(2), - RS.at(2), 0, 0},
    {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     0,0,0,1/(RD * ws),0,0,0,0,0,0,TSV/delta_t + 1},
     {0, 0, 0,0,0,0,0,0,0,0,0,      0, 0, 0,0,0,0,0,0,0,0,0,     0,0,0,0,0,0,0,0,0,TCH/delta_t + 1, -1}}; 
    std::vector<double>b{Tdop.at(0) * Y.at(i).at(0) /delta_t,Tqop.at(0) * Y.at(i).at(1)/delta_t, Y.at(i).at(2)/delta_t - ws, 2 * H.at(0) /ws * Y.at(i).at(3) / delta_t + Dm.at(0) * ws, TE/delta_t * Y.at(i).at(4), TF * Y.at(i).at(5)/delta_t, TA * Y.at(i).at(6)/delta_t + KA *(Vref.at(0) - Vg.at(0)),    Vg.at(0) * sin(Y.at(i).at(2) - TH.at(0)),Vg.at(0) * cos(Y.at(i).at(2) - TH.at(0)),    TSV * Y.at(i).at(10)/delta_t + 1/RD  + PC.at(0),TCH * Y.at(i).at(9) / delta_t,
    Tdop.at(1) * Y.at(i).at(11) /delta_t,Tqop.at(1) * Y.at(i).at(12)/delta_t,Y.at(i).at(13)/delta_t - ws, 2 * H.at(1) /ws * Y.at(i).at(14) / delta_t + Dm.at(1) * ws, TE/delta_t * Y.at(i).at(15), TF * Y.at(i).at(16)/delta_t, TA * Y.at(i).at(17)/delta_t + KA *(Vref.at(1) - Vg.at(1)),     Vg.at(1) * sin(Y.at(i).at(13) - TH.at(1)),Vg.at(1) * cos(Y.at(i).at(13) - TH.at(1)),TSV * Y.at(i).at(21)/delta_t + 1/RD  + PC.at(1),TCH * Y.at(i).at(20) / delta_t, 
    Tdop.at(2) * Y.at(i).at(22) /delta_t,Tqop.at(2) * Y.at(i).at(23)/delta_t,Y.at(i).at(24)/delta_t - ws, 2 * H.at(2) /ws * Y.at(i).at(25) / delta_t + Dm.at(2) * ws, TE/delta_t * Y.at(i).at(26), TF * Y.at(i).at(27)/delta_t, TA * Y.at(i).at(28)/delta_t + KA *(Vref.at(2) - Vg.at(2)),    Vg.at(2) * sin(Y.at(i).at(24) - TH.at(2)),Vg.at(2) * cos(Y.at(i).at(24) - TH.at(2)),TSV * Y.at(i).at(32)/delta_t + 1/RD  + PC.at(2),TCH * Y.at(i).at(31) / delta_t};

    Y.push_back(ProizvodMatriceIKolone(inverse_matrix(A,A.size()), b));
    if(Y.at(i+1).at(6)>1.085)Y.at(i+1).at(6) = 1.085;
    if(Y.at(i+1).at(6)<1.073)Y.at(i+1).at(6) = 1.073;

    //std::cout<<"\n";

    //std::cout<<"\n nesto: "<<DaLiJeKvadratna(A);
    t += delta_t;
  }
  return Y;

  

  
}
int main() {

  std::vector<std::vector<double>> Y;

  std::cout << "Unesite korak: ";
  double korak;
  std::cin >> korak;

  std::cout << "Do kojeg trenutka: ";
  double ts = 0;
  std::cin >> ts;
  Y = MojDAESolverPrvogPrimjera(korak, ts);
  
  std::cout << "y1 je nakon " << ts << "s = " << Y.at(ts / korak - 1).at(0)
            << "dok je y2 nakon " << ts << "s = " << Y.at(ts / korak - 1).at(3);
  std::cout << std::endl;
  std::cout<<Dm.at(0);
  //for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(3)<<std::endl;

  
            std::ofstream out("Eqp1.txt");
            std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
            std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(0)<<std::endl;

            std::ofstream out2("w1.txt");
            std::cout.rdbuf(out2.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(3)<<std::endl;

            std::ofstream out3("w2.txt");
            std::cout.rdbuf(out3.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(14)<<std::endl;

            std::ofstream out4("w3.txt");
            std::cout.rdbuf(out4.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(25)<<std::endl;

            std::ofstream out5("Edp1.txt");
            std::cout.rdbuf(out5.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(1)<<std::endl;

            std::ofstream out6("delta1.txt");
            std::cout.rdbuf(out6.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(2)<<std::endl;


            std::ofstream out7("Efd1.txt");
            std::cout.rdbuf(out7.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(4)<<std::endl;

            std::ofstream out8("Vr1.txt");
            std::cout.rdbuf(out8.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(6)<<std::endl;

             std::ofstream out9("psv1.txt");
            std::cout.rdbuf(out9.rdbuf()); //redirect std::cout to out.txt!
            for(int i = 0; i < ts/korak;i++)std::cout<<Y.at(i).at(10)<<std::endl;





}
