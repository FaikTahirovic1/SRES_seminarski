#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <functional>
/*
  DAE koji rjesavamo je
  y_1 + ni * t * y_2 = q(t)  algebarska
  dot(y_1) + ni * t * dot(y_2) + (1+ ni)y_2 = 0; diferencijalna
  dot(y_1) + ni * y_2 + ni * t * dot(y_2) = dot(q(t))
  */

double dajDeterminantu(const std::vector<std::vector<double>> vect) {
    if(vect.size() != vect.at(0).size()) {
        throw std::runtime_error("Matrica nije kvadratna");
    } 
    int dimenzija = vect.size();

    if(dimenzija == 0) {
        return 1;
    }

    if(dimenzija == 1) {
        return vect.at(0).at(0);
    }

    if(dimenzija == 2) {
        return vect.at(0).at(0) * vect.at(1).at(1) - vect.at(0).at(1) * vect.at(1).at(0);
    }

    double rezultat = 0;
    int znak = 1;
    for(int i = 0; i < dimenzija; i++) {

        std::vector<std::vector<double>> subVect(dimenzija - 1, std::vector<double> (dimenzija - 1));
        for(int m = 1; m < dimenzija; m++) {
            int z = 0;
            for(int n = 0; n < dimenzija; n++) {
                if(n != i) {
                    subVect.at(m - 1).at(z) = vect.at(m).at(n);
                    z++;
                }
            }
        }

        rezultat = rezultat + znak * vect.at(0).at(i) * dajDeterminantu(subVect);
        znak = -znak;
    }

    return rezultat;
}

std::vector<std::vector<double>> TransponujMatricu(const std::vector<std::vector<double>> matrica1) {


    std::vector<std::vector<double>> rjesenje(matrica1.at(0).size(), std::vector<double> (matrica1.size()));

    //Filling rjesenje-matrica
    for(size_t i = 0; i < matrica1.size(); i++) {
        for(size_t j = 0; j < matrica1[0].size(); j++) {
            rjesenje.at(j).at(i) = matrica1.at(i).at(j);
        }
    }
    return rjesenje;
}

std::vector<std::vector<double>> DajKofaktorMatrice(const std::vector<std::vector<double>> vect) {
    if(vect.size() != vect.at(0).size()) {
        throw std::runtime_error("Matrica nije kvadratna");
    } 

    std::vector<std::vector<double>> rjesenje(vect.size(), std::vector<double> (vect.size()));
    std::vector<std::vector<double>> subVect(vect.size() - 1, std::vector<double> (vect.size() - 1));

    for(std::size_t i = 0; i < vect.size(); i++) {
        for(std::size_t j = 0; j < vect.at(0).size(); j++) {

            int p = 0;
            for(size_t x = 0; x < vect.size(); x++) {
                if(x == i) {
                    continue;
                }
                int q = 0;

                for(size_t y = 0; y < vect.size(); y++) {
                    if(y == j) {
                        continue;
                    }

                    subVect.at(p).at(q) = vect.at(x).at(y);
                    q++;
                }
                p++;
            }
            rjesenje.at(i).at(j) = pow(-1, i + j) * dajDeterminantu(subVect);
        }
    }
    return rjesenje;
}
std::vector<std::vector<double>> dajInverznu(const std::vector<std::vector<double>> vect) {
    if(dajDeterminantu(vect) == 0) {
        throw std::runtime_error("Determinant is 0");
    } 

    double d = 1.0/dajDeterminantu(vect);
    std::vector<std::vector<double>> rjesenje(vect.size(), std::vector<double> (vect.size()));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            rjesenje.at(i).at(j) = vect.at(i).at(j); 
        }
    }

    rjesenje = TransponujMatricu(DajKofaktorMatrice(rjesenje));

    for(size_t i = 0; i < vect.size(); i++) {
        for(size_t j = 0; j < vect.size(); j++) {
            rjesenje.at(i).at(j) *= d;
        }
    }

    return rjesenje;
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
  y2.push_back(1 / ni);
  Y.push_back({0, 1 / ni});
  t = delta_t;
  for (int i = 0; i < trenutak / korak; i++) {
    /*
    y1 + ni*t*y2 = q(t)
    y1 - y1(t-1) /delta_t + ni *t * (y2 - y2(t-1))/delta_t  + y_2 (1+ni) = 0;    
      */
    
    std::vector<std::vector<double>>A{{1,ni * t},{1/delta_t, ni * t / delta_t + 1 + ni}};
    std::vector<double>b{q(t),Y.at(i).at(0) / delta_t + ni * t * Y.at(i).at(1) / delta_t};
    Y.push_back(ProizvodMatriceIKolone(dajInverznu(A), b));
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
  std::cout << "Relativne greške u postotcima su: "
            << (Y.at(ts / korak - 1).at(0) - sin(ts) - ni * cos(ts) * ts) *
                   100 / (sin(ts) + ni * cos(ts) * ts)
            << "% i "
            << (Y.at(ts / korak - 1).at(1) + cos(ts)) * 100 / (-cos(ts)) << "%"
            << std::endl;

}
