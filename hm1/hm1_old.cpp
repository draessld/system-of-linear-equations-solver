#include <iostream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
using namespace std;

#define EPSILON 1e-5

int n;

//  power of p
double pow(double a, double p)
{
    double res = a;
    for (int i = 1; i < p; i++)
        res *= a;
    return a;
}

//  return lenght of vector
double vector_lenght(double *vector)
{
    double sum;
    for (int i = 0; i < n; i++)
        sum += pow(vector[i], 2);
    return pow(sum,0.5);
}

//  control print of matrix
void print_m(double **a)
{

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%lf ", a[i][j]);
        printf("\n");
    }
}

//  control print of vector
void print_v(double *b)
{
    for (int j = 0; j < n; j++)
        printf("%lf ", b[j]);
    printf("\n");
}



//  matice * sloupcovy vector
double *muptiply_m_v(double **a, double *b)
{
    double *res = new double[n];
    for (int i = 0; i < n; i++)
    {
        res[i] = 0;
        for (int j = 0; j < n; j++)
            res[i] += a[i][j] * b[j];
    }
    return res;
}

double *subtract_v_v(double *a, double *b)
{
    double *res = new double[n];
    for (int i = 0; i < n; i++)
        res[i] = a[i] - b[i];
    return res;
}

double *addition_v_v(double *a, double *b)
{
    double *res = new double[n];
    for (int i = 0; i < n; i++)
        res[i] = a[i] + b[i];
    return res;
}

//  sloupcovy vector * radkovy vector
double multiply_v_v(double *a, double *b)
{
    double sum;
    for (int i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

//  vector zvetseny o konstantu
double *multiply_v_c(double *a, double c)
{
    double *res = new double[n];
    for (int i = 0; i < n; i++)
        res[i] = a[i] * c;
    return res;
}

int main(int argc, char *argv[])
{
    // int opt;
    // while ((opt = getopt(argc, argv, ":ism")) != -1)
    // {
    //     switch (opt)
    //     {
    //     case 's':
    //         cout << "matrix will be saved as " << endl;

    //         break;
    //     case 'm':
    //         break;
    //     case 'f':
    //         cout << "input file is " << endl;
    //         break;
    //     default: /* '?' */
    //         fprintf(stderr, "Usage: %s [-f input_file] [-s matrix_format] [-m mathod]\n",
    //                 argv[0]);
    //         exit(EXIT_FAILURE);
    //     }
    // }

    /*  READ INPUTS */
    //  Matici a vektor prave strany nactete ze souboru.
    fstream file(argv[1]);

    //  read size of matrix n
    file >> n;

    //  initialize matrix
    double **A = new double *[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];

    //  read matrix
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            file >> A[i][j];

    //  initialize and read vector and start point
    double *b = new double[n];
    double *x = new double[n];
    for (int i = 0; i < n; i++)
    {
        file >> b[i];
        x[i] = 0;
    }

    double a;  //    step lenght

    double *gk = new double[n];
    double *g0 = new double[n];
    double e = 0;

    gk = subtract_v_v(b, muptiply_m_v(A, x));
    g0 = subtract_v_v(b, muptiply_m_v(A, x));

    do
    {
        a = multiply_v_v(gk, gk) / multiply_v_v(gk, muptiply_m_v(A, gk));
        x = addition_v_v(x, multiply_v_c(gk, a));
        gk = subtract_v_v(b, muptiply_m_v(A, x));
        print_v(x);
        e = (vector_lenght(gk) / vector_lenght(g0));
        cout << "e: " << e << endl;
        (e > EPSILON);
    } while (e > EPSILON);

    //  tidy up

    return 0;
}
