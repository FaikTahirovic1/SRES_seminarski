#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <complex>
#include <string>

using namespace std;
using namespace Eigen;
const double PI = std::atan(1.0) * 4.0;
const std::complex<double> i_u(0.0, 1.0);

std::complex<double> convertToComplex(const std::string& str) {
    std::size_t plusPos = str.find('+');
    std::size_t iPos = str.find('i');
    std::string realStr = str.substr(0, plusPos);
    std::string imagStr = str.substr(plusPos + 1, iPos - plusPos - 1);
    double real = std::stod(realStr);
    double imag = std::stod(imagStr);
    return std::complex<double>(real, imag);
}

void saveData(string fileName, MatrixXd  matrix)
{
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
	ofstream file(fileName);
	if (file.is_open())
	{
		file << matrix.format(CSVFormat);
		file.close();
	}
}

 MatrixXcd openDataComplex(string fileToOpen)
{
	vector<std::complex<double>> matrixEntries;
	ifstream matrixDataFile(fileToOpen);
	string matrixRowString;
	string matrixEntry;
	int matrixRowNumber = 0;
	while (getline(matrixDataFile, matrixRowString))
	{
		stringstream matrixRowStringStream(matrixRowString);

		while (getline(matrixRowStringStream, matrixEntry, ','))
		{
		 	matrixEntries.push_back(convertToComplex(matrixEntry));
		}
		matrixRowNumber++;
	}

	return Map<Matrix<std::complex<double>, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}

MatrixXd openData(string fileToOpen)
{
	vector<double> matrixEntries;
	ifstream matrixDataFile(fileToOpen);
	string matrixRowString;
	string matrixEntry;
	int matrixRowNumber = 0;


	while (getline(matrixDataFile, matrixRowString))
	{
		stringstream matrixRowStringStream(matrixRowString);

		while (getline(matrixRowStringStream, matrixEntry, ','))
		{
			matrixEntries.push_back(stod(matrixEntry));
		}
		matrixRowNumber++;
	}

	return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}

Eigen::MatrixXcd computeIphasor(MatrixXd PG, MatrixXd QG, MatrixXcd Vphasor, int m){
MatrixXcd Iphasor;
PG = PG(seq(0,m-1),0);
QG = QG(seq(0,m-1),0);
Iphasor = PG.transpose() + i_u*QG.transpose();
for(int i=0; i<m; i++){
    Vphasor(0,i) = 1./(Vphasor(0,i));
}
Iphasor = Iphasor.cwiseProduct(Vphasor);
Iphasor = Iphasor.conjugate();
return Iphasor;
}

Eigen::MatrixXd computeEm(Eigen::MatrixXcd E0, int m){
 Eigen::MatrixXd Em(1,m);
 for(int i =0; i<m; i++){
    Em(0,i) = std::abs(E0(0,i));
 }
return Em;
}

Eigen::MatrixXd computeD0(Eigen::MatrixXcd E0, int m){
 Eigen::MatrixXd D0(1,m);
 for(int i =0; i<m; i++){
    D0(0,i) = std::arg(E0(0,i));
 }
return D0;
}

Eigen::MatrixXd computeId0(Eigen::MatrixXcd Iphasor, Eigen::MatrixXd D0, int m){
 MatrixXcd mat(1,m);
 for(int i=0; i<m; i++){
    mat(0,i) = exp(-1.*i_u*( D0(0,i)-PI/2. ));
 }
 mat = mat.cwiseProduct(Iphasor);
 MatrixXd r(1,m);
 for(int i=0; i<m; i++){
    r(0,i) = std::real(mat(0,i));
 }
 return r;
}

Eigen::MatrixXd computeIq0(Eigen::MatrixXcd Iphasor, Eigen::MatrixXd D0, int m){
 MatrixXcd mat(1,m);
 for(int i=0; i<m; i++){
    mat(0,i) = exp(-1.*i_u*( D0(0,i)-PI/2. ));
 }
 mat = mat.cwiseProduct(Iphasor);
 MatrixXd r(1,m);
 for(int i=0; i<m; i++){
    r(0,i) = std::imag(mat(0,i));
 }
 return r;
}

int main(){
int m,n;
Eigen::MatrixXd H, Xd, Xdp, Xq, Xqp, Td0p, Tq0p,Rs, Dm,KA, TA, KE, TE, KF, TF, Ax, Bx,TCH, TSV, RD,MH;
m = 3; n = 9;
double ws = 2*PI*60;
Eigen::MatrixXd ws_vector = ws*Eigen::MatrixXd::Ones(1,m);
Eigen::MatrixXd w0 = ws_vector, Yabs, Yang;
double baseMVA = 100;
MatrixXcd Ybus;
MatrixXd bus, branch, gen, MD, ED, TD;
// ucitavanje podatak iz CSV fajlova u predvidjene matrice
// za ucitavanje CSV fajlova napisane su funkcije openDataComplex i openData
// u zavisnoti kojeg su tipa elementi matrica
Ybus = openDataComplex("Ybus_data.txt");
bus = openData("bus_data.txt");
branch = openData("branch_data.txt");
gen = openData("gen_data.txt");
MD = openData("Machine_data.txt");
ED = openData("Excitation_data.txt");
TD = openData("Turbine_data.txt");

Eigen::MatrixXd IC1, IC2,gen0, genP, IC3,genQ, IC4, IC5, IC6,IC(9,6),PL,QL,PG,QG,TH0,V0,VG0(1,m);
Eigen::MatrixXcd THG0(1,m);
IC1 = bus.col(7); IC2 = bus.col(8);
gen0 = Eigen::MatrixXd::Zero(n,1);
for(int i=0; i<m; i++){
    gen0(i,0) = gen(i,1);
}
genP = gen0; IC3 = genP;
IC3 = IC3/baseMVA;
gen0 = Eigen::MatrixXd::Zero(n,1);
for(int i=0; i<m; i++){
    gen0(i,0) = gen(i,2);
}
genQ = gen0; genQ = genQ + bus.col(5);
IC4 = genQ;
IC4 = IC4/baseMVA;
IC5 = bus.col(2); IC5 = IC5/baseMVA;
IC6 = bus.col(3); IC6 = IC6/baseMVA;
IC << IC1, IC2,IC3,IC4,IC5,IC6;
PL = IC.col(4); QL = IC.col(5);
PG = IC.col(2); QG = IC.col(3);
TH0 = IC.col(1)*PI/180.; TH0 = TH0.transpose();
V0 = IC.col(0); V0 = V0.transpose();
for(int i=0; i<m; i++){
    VG0(0,i) = V0(0,i);
    THG0(0,i) = TH0(0,i);
}

H = MD.row(0);
Xd = MD.row(1);
Xdp = MD.row(2);
Xq = MD.row(4);
Xqp = MD.row(5);
Td0p = MD.row(7);
Tq0p = MD.row(9);
Rs = MD.row(11);
Dm = MD.row(13);
KA = ED.row(0);
TA = ED.row(1);
KE = ED.row(2);
TE = ED.row(3);
KF = ED.row(4);
TF = ED.row(5);
Ax = ED.row(6);
Bx = ED.row(7);
TCH = TD.row(0);
TSV = TD.row(1);
RD = TD.row(2);
MH = (H*2)/ws;

// Racunanje pocetnih vrijednosti za direfencijalne i algebarske jednacine
Eigen::MatrixXcd Vphasor, Iphasor, E0;
Eigen::MatrixXd Em, D0, Id0, Iq0,Edp0, Eqp0, TM0, VR0, Efd0,RF0,Vref,PSV0,PC;
for(int i =0; i<m; i++) {
    THG0(0,i) = exp(i_u*THG0(0,i));
}
Vphasor = VG0.cwiseProduct(THG0);
Iphasor = computeIphasor(PG,QG,Vphasor,m);
E0 = Vphasor + (Rs + i_u*Xq).cwiseProduct(Iphasor);
Em = computeEm(E0,m);
D0 = computeD0(E0,m);
Id0 = computeId0(Iphasor, D0, m);
Iq0 = computeIq0(Iphasor, D0, m);
Edp0 = (Xq - Xqp).cwiseProduct(Iq0);
Eqp0 = Rs.cwiseProduct(Iq0) + Xdp.cwiseProduct(Id0)+V0(0,seq(0,m-1)).cwiseProduct( (MatrixXd)(((D0-TH0(0,seq(0,m-1))).array()).cos()) );
Efd0 = Eqp0 + (Xd - Xdp).cwiseProduct(Id0);
TM0 = Eqp0.cwiseProduct(Iq0) + Edp0.cwiseProduct(Id0) + ((Xqp-Xdp).cwiseProduct(Id0)).cwiseProduct(Iq0);
VR0 = (KE + Ax.cwiseProduct(MatrixXd(((Bx.cwiseProduct(Efd0)).array()).exp()))).cwiseProduct(Efd0);
RF0 = (MatrixXd(KF.array()/TF.array())).cwiseProduct(Efd0);
Vref = V0(0,seq(0,m-1)) + MatrixXd(VR0.array()/KA.array());
PSV0 = TM0;
PC = PSV0;

// Vektor pocetnih vrijednosti
Eigen::MatrixXd x0(1,11*m);
Eigen::MatrixXd a0(1,2*m+2*n);
x0(0,seq(0,m-1)) = Eqp0;
x0(0,seq(m,2*m-1)) = Edp0;
x0(0,seq(2*m,3*m-1)) = D0;
x0(0,seq(3*m,4*m-1)) = ws_vector;
x0(0,seq(4*m,5*m-1)) = Efd0;
x0(0,seq(5*m,6*m-1)) = RF0;
x0(0,seq(6*m,7*m-1)) = VR0;
x0(0,seq(7*m,8*m-1)) = TM0;
x0(0,seq(8*m,9*m-1)) = PSV0; // pocetne vrijednosti za varijable stanja

a0(0,seq(0,m-1)) = Id0;
a0(0,seq(m,2*m-1)) = Iq0;
a0(0,seq(2*m,2*m+n-1)) = V0;
a0(0,seq(2*m+n,2*m+2*n-1)) = TH0; // pocetne vrijednosti algebarskih jednacina

// racunanje matrice Yred
V0 = V0.transpose();
Eigen::MatrixXcd YL = (PL - i_u*QL);
YL = MatrixXcd( YL.array() / ((V0.array()).square())  );
Eigen::MatrixXcd YLdiag(n,n);
YLdiag.setZero();
for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        if(i==j) {
            YLdiag(i,j) = YL(i,0);
        }

    }
}

Eigen::MatrixXcd Y_Aug = Ybus + YLdiag;
Eigen::MatrixXcd Y11(m,m), Y12(m,n-m), Y21(n-m,m), Y22(n-m,n-m), Yred, YredInv;
Y11 = Y_Aug(seq(0,m-1), seq(0,m-1));
Y12 = Y_Aug(seq(0,m-1), seq(m,n-1));
Y21 = Y_Aug(seq(m,n-1), seq(0,m-1));
Y22 = Y_Aug(seq(m,n-1), seq(m,n-1));
Yred = Y11-Y12*Y22.inverse()*Y21;
YredInv = Yred.inverse();

Eigen::MatrixXd Eqp, Edp, Delta,w,Efd,RF,VR, TM, PSV, Id, Iq, V, TH;
Eqp = x0(0,seq(0,m-1));
Edp = x0(0,seq(m,2*m-1)) ;
Delta = x0(0,seq(2*m,3*m-1));
w = x0(0,seq(3*m,4*m-1));
Efd = x0(0,seq(4*m,5*m-1));
RF = x0(0,seq(5*m,6*m-1));
VR = x0(0,seq(6*m,7*m-1));
TM = x0(0,seq(7*m,8*m-1));
PSV = x0(0,seq(8*m,9*m-1));
Id = a0(0,seq(0,m-1));
Iq = a0(0,seq(m,2*m-1));
V = a0(0,seq(2*m,2*m+n-1));
TH = a0(0,seq(2*m+n,2*m+2*n-1));

std::cout << "Eqp: "<< Eqp <<std::endl;
std::cout << "Edp: "<< Edp <<std::endl;
std::cout << "Delta: "<< Delta <<std::endl;
std::cout << "w: "<< w <<std::endl;
std::cout << "Efd: "<< Efd <<std::endl;
std::cout << "RF: "<< RF <<std::endl;
std::cout << "VR: "<< VR <<std::endl;
std::cout << "TM: "<< TM <<std::endl;
std::cout << "PSV: "<< PSV <<std::endl;
std::cout << "Id: "<<Id <<std::endl;
std::cout << "Iq: "<< Iq <<std::endl;
std::cout << "V: "<< V <<std::endl;
std::cout << "TH: "<< TH <<std::endl;
std::cout << "PC: "<< PC <<std::endl;
std::cout << "Vref" << Vref << std::endl;
std::cout <<"PC: " << PC << std::endl;
std::cout <<"Yred: "<<Yred << std::endl;

return 0;
}

