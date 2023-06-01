#include <iostream>
#include <math.h>
#include <vector>
/*
  DAE koji rjesavamo je
  y_1 + ni * t * y_2 = q(t)  algebarska
  dot(y_1) + ni * t * dot(y_2) + (1+ ni)y_2 = 0; diferencijalna
  dot(y_1) + ni * y_2 + ni * t * dot(y_2) = dot(q(t))
  */

std::vector<double> RjesiSistem2x2(double a, double b, double c, double d,
                                   double e, double f) {
  /* ax + by = c
  dx +ey = f
  d*a*x  d*b*y = c * d
  a*d*x + a* e *y = a * f
  y(d * b - a * e) = c * d - a * f; bitno
  aex + bey = ce
  dbx + eby = fb
  x(ae - db) = ce - fb
*/

  return std::vector<double>{(c * e - f * b) / (a * e - d * b),
                             (c * d - a * f) / (d * b - a * e)};
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

  y1.push_back(0);
  double delta_t = korak, t = 0;
  y2.push_back(1 / ni);
  Y.push_back({0, 1 / ni});
  t = delta_t;

  std::cout << "Do kojeg trenutka: ";
  double ts = 0;
  std::cin >> ts;
  for (int i = 0; i < ts / korak; i++) {

    /*
    y1 + ni*t*y2 = q(t)
    y1 - y1(t-1) /delta_t + ni *t * (y2 - y2(t-1))/delta_t  + y_2 (1+ni) = 0;


      */
    Y.push_back(RjesiSistem2x2(
        1, ni * t, sin(t), 1 / delta_t, ni * t / delta_t + 1 + ni,
        Y.at(i).at(0) / delta_t + ni * t * Y.at(i).at(1) / delta_t));
    t += delta_t;
  }
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
