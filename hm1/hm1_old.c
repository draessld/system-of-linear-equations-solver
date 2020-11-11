#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#define EPSILON 1e-5

size_t n = 3;

//  return lenght of vector
double vector_lenght(double *vector)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += pow(vector[i], 2.);
    return sqrt(sum);
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

//  read matrix as 2d
int read_matrix(double **a, const char *filename)
{

    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            fscanf(pf, "%lf", a[i] + j);
    }

    fclose(pf);
    return 1;
}

//  read vector
int read_vector(double *a, const char *filename)
{
    FILE *pf;
    pf = fopen(filename, "r");
    if (pf == NULL)
        return 0;

    for (int i = 0; i < n; ++i)
        fscanf(pf, "%lf", &a[i]);

    fclose(pf);
    return 1;
}

//  matrix * column vector
void muptiply_m_v(double **a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
    {
        res[i] = 0;
        for (int j = 0; j < n; j++)
        {
            // printf("%lf + (%lf * %lf) = %lf + %lf  ",res[j],a[j][i], b[i],res[j], a[j][i] * b[i]);
            res[i] += (a[i][j] * b[j]);
            // printf(" = %lf; \n",res[j]);
        }
    }
}

//  substract two vectors
void subtract_v_v(double *a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] - b[i];
}

//  vector zvetseny o konstantu
void multiply_v_c(double *a, double c)
{
    for (int i = 0; i < n; i++)
        a[i] *= c;
}

//  create copy of vector
void copy_vector(double *s, double *r)
{
    for (int i = 0; i < n; i++)
        r[i] = s[i];
}

//  additive two vector
void addition_v_v(double *a, double *b, double *res)
{
    for (int i = 0; i < n; i++)
        res[i] = a[i] + b[i];
}

//  sloupcovy vector * radkovy vector
double multiply_v_v(double *a, double *b)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += a[i] * b[i];
    return sum;
}

//  gradient descent method
void gradient_descent_method(double **matrix, double *b, double *x)
{

    //  initialization vector x
    double pom[n];
    for (int i = 0; i < n; i++)
    {
        x[i] = 0;
        pom[i] = 0;
    }

    double *gk = (double *)malloc(n * sizeof(double));
    muptiply_m_v(matrix, x, pom);
    subtract_v_v(b, pom, pom);

    copy_vector(pom, gk);
    double lgo = vector_lenght(gk);
    double a = 0;   //  lenght of step
    double e = 0; //  difference

    //  lets do it! :))
    do
    {
        double y = multiply_v_v(gk, gk);
        // printf("y: %lf \n",y);
        muptiply_m_v(matrix, gk, pom);
        // print_v(pom);
        a = y / multiply_v_v(gk, pom);
        // printf("a: %lf \n", a);
        multiply_v_c(gk, a);
        // print_v(gk);
        addition_v_v(x, gk, x);
        muptiply_m_v(matrix, x, pom);
        subtract_v_v(b, pom, gk);
        // print_v(x);
        double lgk = vector_lenght(gk);
        if (lgo == 0)
            e = lgk;
        else
            e = lgk / lgo;
        // printf("e: %lf \n",e);
    } while (e > EPSILON);

    free(gk);
}

//  main
int main(int argc, char const *argv[])
{

    /*
    //  declare matrix
    double **matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        matrix[i] = (double *)malloc(n * sizeof(double));

    // declare vector
    double *vector;

    char c;
    int nl = 0; //  number of newlines

    //  read with unknown lenght
    while (fscanf(input, "%c", &c) != EOF)
    {
        printf("%c",c);
        if (c == '\n')
        {
            //  number of newline ++
            nl++;
            if (nl == 2)
            {
                //  read vector
                nl = -1;
                i = j = 0;
                //  set size of matrix
                n = sizeof(matrix[0]);
                vector = (double *)malloc(n * sizeof(double));
            }
            else
            {
                //  new line of matrix
                i++;
            }
        }
        else
        {
            //  check if realloc required
            if (j >= n)
            {
                n *= 2;
                matrix = (double **)realloc(matrix, n);
                for (int k = 0; k < n; k++)
                    matrix[k] = (double *)realloc(matrix, n);
            }
            if (nl == -1)
            {
                //  fill vector
                vector[i] = (double)c;
            }
            else
            {
                //  fill cell
                matrix[i][j] = (double)c;
                j++;
            }
        }
    }*/

    //  declare matrix
    double **matrix = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        matrix[i] = (double *)malloc(n * sizeof(double));

    // declare vector
    double *vector = (double *)malloc(n * sizeof(double));

    if (read_matrix(matrix, argv[1]) == 0)
    {
        printf("filename %s does not exist", argv[1]);
        return 1;
    }
    if (read_vector(vector, argv[2]) == 0)
    {
        printf("filename %s does not exist", argv[2]);
        return 1;
    }

    // print_m(matrix);
    // print_v(vector);

    double *output = (double *)malloc(n * sizeof(double));
    gradient_descent_method(matrix, vector, output);

    //  print result
    print_v(output);

    //  tidy up
    for (int i = 0; i < n; i++)
        free(matrix[i]);
    free(matrix);
    free(vector);
    free(output);

    return 0;
}
