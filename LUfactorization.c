#include<stdio.h>
#include<math.h>

#define MAX_SIZE 100 /* 矩阵最大维数 */
#define ZERO 0.000000001 /* 当一个正数小于ZERO就认为该数是0 */

typedef enum {false, true} bool;
bool Direct( int n, double a[][MAX_SIZE], double *b )
{
    double u[n][n], l[n][n];
    for(int i=0;i<n;i++)
        u[0][i] = a[0][i];
    for(int i=1;i<n;i++)
        l[i][0] = a[i][0]/u[0][0];
    for(int r=1;r<n;r++)
    {
        for(int i=r;i<n;i++)
        {
            double sum=0;
            for(int k=0;k<=r-1;k++)
                sum+=l[r][k]*u[k][i];
            u[r][i] = a[r][i] - sum;
        }
        for(int i=r+1;i<n&&r!=n-1;i++)
        {
            double sum=0;
            for(int k=0;k<=r-1;k++)
                sum+=l[i][k]*u[k][r];
            l[i][r]=(a[i][r] - sum)/u[r][r];
        }
    }
    for(int i=0;i<n;i++)
        if(u[i][i]==0)
            return false;
    double y[n],x[n];
    y[0] = b[0];
    for(int i=1;i<n;i++)
    {
        double sum=0;
        for(int k=0;k<=i-1;k++)
            sum += l[i][k]*y[k];
        y[i] = b[i] - sum;
    }
    x[n-1] = y[n-1]/u[n-1][n-1];
    for(int i=n-2;i>-1;i--)
    {
        double sum=0;
        for(int k=i+1;k<n;k++)
            sum += u[i][k]*x[k];
        x[i] = (y[i]-sum)/u[i][i];
    }
    for(int i=0;i<n;i++)
        b[i]=x[i];
    return true;
}
int main()
{
  int n, i, j;
  double a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE];
  while ( scanf("%d", &n) != EOF ) {
    for ( i=0; i<n; i++ ) {
      for ( j=0; j<n; j++ )
        scanf("%lf", &a[i][j]);
      scanf("%lf", &b[i]);
    }
    if ( Direct(n, a, b) ) {
      printf("Result of direct method:\n");
      for ( j=0; j<n; j++ )
        printf("%.8f\n", b[j]);
    }
    else
      printf("Doolittle factorization failed.\n");
    printf("\n");
  }
  return 0;
}
