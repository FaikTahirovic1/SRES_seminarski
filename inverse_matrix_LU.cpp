#include <iostream>

using namespace std;

int dimension;  //ostaje globalna varijabla, jer se nigdje ne prenosi u funkciji, moze se izmijeniti da se i prenosi kao parametar u funkciji logy

const int N = 100; //moze se povecati, na pocetku se alocira vektor ove duzine, a u mainu se trazi tacna dimenzija matrice

double *mat[N];
double *lower[N], *upper[N], *inverse[N], *Z[N], *I[N];;

void lu()
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < dimension; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            if (j < i)
                lower[j][i] = 0;
            else
            {
                lower[j][i] = mat[j][i];
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
                upper[i][j] = mat[i][j] / lower[i][i];
                for (k = 0; k < i; k++)
                {
                    upper[i][j] = upper[i][j] - ((lower[i][k] * upper[k][j]) / lower[i][i]);
                }
            }
        }
    }
}

void output(int *matrix[dimension])
{
    int i = 0, j = 0;
    for (i = 0; i < dimension; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            cout << matrix[i][j] << '\t';
        }
        cout << endl;
    }
}


void output(float *matrix[dimension])
{
    int i = 0, j = 0;
    for (i = 0; i < dimension; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            cout << matrix[i][j] << '\t';
        }
        cout << endl;
    }
}

void output(double *matrix[dimension])
{
    int i = 0, j = 0;
    for (i = 0; i < dimension; i++)
    {
        for (j = 0; j < dimension; j++)
        {
            cout << matrix[i][j] << '\t';
        }
        cout << endl;
    }
}

double compute_z(int col, int row) {
    double sum = 0;
    for(int i = 0; i < dimension; i++) {
        if(i != row) {
            sum += lower[row][i] * Z[i][col];
        }
    }

    double result = I[row][col] - sum;
    result = result / lower[row][row];

    return result;
}

double compute_inverse(int col, int row) {
    double sum = 0;
    for(int i = 0; i < dimension; i++) {
        if(i != row) {
            sum += upper[row][i] * inverse[i][col];
        }
    }

    double result = Z[row][col] - sum;
    result = result / upper[row][row];


    return result;
}


void inverse_matrix() {

    // compute z
    for(int col = 0; col < dimension; col++) {
        for(int row = 0; row < dimension; row++) {
            Z[row][col] = compute_z(col, row);
        }
    }
;
    // compute inverse
    for(int col = 0; col < dimension; col++) {
        for(int row = dimension - 1; row >= 0; row--) {
            inverse[row][col] = compute_inverse(col, row);
        }
    }
}

// Driver code
int main() {

    cout<< "Input dimension of a matrix: ";
    cin>>dimension;
  
    for (int i = 0; i < dimension; i++) {
        mat[i] = new double[dimension];
        lower[i] = new double[dimension];
        upper[i] = new double[dimension];
        inverse[i] = new double[dimension];
        Z[i] = new double[dimension];
        I[i] = new double[dimension];
    }
  
    cout << "Input a matrix: " << endl;

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            cin >> mat[i][j];
            Z[i][j] = 0;
            if(i == j)
                I[i][j] = 1;
        }
    }

    lu();
    inverse_matrix();
    
    cout << "Matrix" << endl;
    output(mat);
    cout << endl << "Lower " << endl;
    output(lower);
    cout << endl <<  "Upper " << endl;
    output(upper);
    cout << endl <<  "Z " << endl;
    output(Z);
    
    cout << endl << "inverse " << endl;
    output(inverse);

    return 0;
}
