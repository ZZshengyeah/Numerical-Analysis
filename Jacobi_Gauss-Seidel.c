#include<stdio.h>
#include<math.h>

#define MAX_SIZE 100 /* 矩阵最大维数 */
#define bound pow(2, 127) /* 判断叠代发散的边界值 */
#define ZERO 0.000000001 /* 当一个正数小于ZERO就认为该数是0 */

double norm(int n, double x[], double y[])
{
    double max=0;
    for(int i=0;i<n;i++)
    {
        double val=fabs(x[i]-y[i]);
        if(val>max)
            max=val;
    }
    return max;
}
void swap(int i,int j, int n, double a[][MAX_SIZE], double b[])
{
    double temp;
    for(int count=0;count<n;count++)
    {
        temp = a[i][count];
        a[i][count] = a[j][count];
        a[j][count] = temp;
        temp = b[i];
        b[i] = b[j];
        b[j] = temp;
    }
}
void change(int n, double a[][MAX_SIZE], double b[])
{
    for(int i=0;i<n;i++)
    {
        double max = a[i][i];
        for(int j=i+1;j<n;j++)
        {
            if(fabs(a[j][i])>fabs(max))
                swap(i,j,n,a,b);
        }
    }
}
int convergence(int n, double x[])
{
    for(int i=0;i<n;i++)
    {
        if(x[i]>bound || x[i]<-bound)
            return 0;
    }
    return 1;
}
int test(int n, double a[][MAX_SIZE])
{
    for(int i=0;i<n;i++)
    {
        int count=0;
        for(int j=0;j<n;j++)
        {
            if(a[j][i]==0)
                count++;
        }
        if(count == n)
            return 0;
    }
    return 1;
}
int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN )
{
    change(n,a,b);
    if(test(n, a)==0)
        return -1;
    double update_x[n];
    int k=1;
    while(1)
    {
        for(int i=0;i<n;i++)
        {
            double sum=0;
            for(int j=0;j<n;j++)
            {
                if(j==i)
                    continue;
                sum += a[i][j]*x[j];
            }
            update_x[i] = (b[i]-sum)/a[i][i];
        }
        if(norm(n,x,update_x)<=TOL)
        {
            for(int i=0;i<n;i++)
                x[i]=update_x[i];
            break;
        }
        if(convergence(n,update_x) == 0)
            return -2;
        for(int i=0;i<n;i++)
            x[i]=update_x[i];
        if(k>MAXN)
            return 0;
        k++;
    }
    return k;
}
int Gauss_Seidel ( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN )
{
    //在进行雅可比迭代时已经进行过选取主元
    if(test(n, a)==0)
        return -1;
    double update_x[n];
    int k=1;
    while(1)
    {
        for(int i=0;i<n;i++)
        {
            double sum1=0,sum2=0;
            for(int j=0;j<=i-1;j++)
                sum1 += a[i][j]*update_x[j];
            for(int j=i+1;j<n;j++)
                sum2 += a[i][j]*x[j];
            update_x[i] = (b[i]-sum1-sum2)/a[i][i];
        }
        if(norm(n,x,update_x)<=TOL)
        {
            for(int i=0;i<n;i++)
                x[i]=update_x[i];
            break;
        }
        if(convergence(n,update_x) == 0)
            return -2;
        for(int i=0;i<n;i++)
            x[i]=update_x[i];
        if(k>MAXN)
            return 0;
        k++;
    }
    return k;
}
int main()
{
  int n, MAXN, i, j, k;
  double a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE], x[MAX_SIZE];
  double TOL;

  while ( scanf("%d", &n) != EOF ) { /* 读取裁判测试用例 */
    for ( i=0; i<n; i++ ) {
      for ( j=0; j<n; j++ )
        scanf("%lf", &a[i][j]);
      scanf("%lf", &b[i]);
    }
    scanf("%lf %d", &TOL, &MAXN);

    /* 输出雅可比算法的结果 */
    printf("Result of Jacobi method:\n");
    for ( i=0; i<n; i++ )
      x[i] = 0.0;
    k = Jacobi( n, a, b, x, TOL, MAXN );
    switch ( k ) {
      case -2:
        printf("No convergence.\n");
        break;
      case -1:
        printf("Matrix has a zero column. No unique solution exists.\n");
        break;
      case 0:
        printf("Maximum number of iterations exceeded.\n");
        break;
      default:
        printf("no_iteration = %d\n", k);
        for ( j=0; j<n; j++ )
          printf("%.8f\n", x[j]);
        break;
    }

    /* 输出高斯-塞德尔算法的结果 */
    printf("Result of Gauss-Seidel method:\n");
    for ( i=0; i<n; i++ )
      x[i] = 0.0;
    k = Gauss_Seidel( n, a, b, x, TOL, MAXN );
    switch ( k ) {
      case -2:
        printf("No convergence.\n");
        break;
      case -1:
        printf("Matrix has a zero column. No unique solution exists.\n");
        break;
      case 0:
        printf("Maximum number of iterations exceeded.\n");
        break;
      default:
        printf("no_iteration = %d\n", k);
        for ( j=0; j<n; j++ )
          printf("%.8f\n", x[j]);
        break;
    }
    printf("\n");
  }

  return 0;
}

