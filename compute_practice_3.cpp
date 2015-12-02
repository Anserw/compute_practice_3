#include  <iostream> 
#include  <math.h> 
#include  <iomanip> 
 
using namespace std; 
const int MAT_SIZE = 4; 
const double PRICISION = 1E-12; 
const double SIGMA = 1E-7; 
const int MAX_ITERATION = 1000; 
 
//获得F的值的向量
void getFValue(double b[], const double* const offset, const double* const x)
{
    b[1] = -(0.5*cos(x[1])+x[2]+x[3]+x[4]-offset[1]); 
    b[2] = -(x[1]+0.5*sin(x[2])+x[3]+x[4]-offset[2]); 
    b[3] = -(0.5*x[1]+x[2]+cos(x[3])+x[4]-offset[3]); 
    b[4] = -(x[1]+0.5*x[2]+x[3]+sin(x[4])-offset[4]); 
}

//获得F的导数矩阵
void getFDerivative(double Fd[][MAT_SIZE+1], const double* const x)
{
    Fd[1][1] = -0.5*sin(x[1]); 
    Fd[1][2] = 1; 
    Fd[1][3] = 1; 
    Fd[1][4] = 1; 
 
    Fd[2][1] = 1; 
    Fd[2][2] = 0.5*cos(x[2]); 
    Fd[2][3] = 1; 
    Fd[2][4] = 1; 
 
    Fd[3][1] = 0.5; 
    Fd[3][2] = 1; 
    Fd[3][3] = -sin(x[3]); 
    Fd[3][4] = 1; 
 
    Fd[4][1] = 1; 
    Fd[4][2] = 0.5; 
    Fd[4][3] = 1; 
    Fd[4][4] = cos(x[4]); 
}
 
//解线性方程组
void solveEquations(double(*mat)[MAT_SIZE+1], double* b, double x[])
{
    for(int k = 1; k < MAT_SIZE; k++) {
        int i = k; 
        for(int j = k; j <= MAT_SIZE; j++) {
            if(mat[j][k] > mat[i][k])
                i = j; 
        }
        for(int j = k; j <= MAT_SIZE; j++) {
            double temp = mat[k][j]; 
            mat[k][j] = mat[i][j]; 
            mat[i][j] = temp; 
        }
        double temp = b[k]; 
        b[k] = b[i]; 
        b[i] = temp; 
        for(int i = k+1; i <= MAT_SIZE; i++) {
            double mik = mat[i][k]/mat[k][k]; 
            for(int j = k+1; j <= MAT_SIZE; j++)
                mat[i][j] = mat[i][j]-mik*mat[k][j]; 
            b[i] = b[i]-mik*b[k]; 
        }
    }
    x[MAT_SIZE] = b[MAT_SIZE]/mat[MAT_SIZE][MAT_SIZE]; 
    for(int k = MAT_SIZE-1; k >= 1; k--) {
        x[k] = 0; 
        for(int j = k+1; j <= MAT_SIZE; j++)
            x[k] -= mat[k][j]*x[j]; 
        x[k] += b[k]; 
        x[k] = x[k]/mat[k][k]; 
    }
}
 
//计算向量的范数, 此处采用无穷范数
double getNorm(const double* const vec)
{
    double temp = 0; 
    for(int i = 1; i <= MAT_SIZE; i++)
        if(temp < fabs(vec[i]))
            temp = fabs(vec[i]); 
    return temp; 
}
 
//Newton法解非线性方程组
void Newton(double x[], const double* const offset)
{
    char a; 
    double delta_x[MAT_SIZE+1]; 
    double tempx[MAT_SIZE+1]; 
    for(int i = 1; i <= MAT_SIZE; i++)
        x[i] = 1; 
    double FDerivative[MAT_SIZE+1][MAT_SIZE+1]; 
    double F[MAT_SIZE+1]; 
    for(int i = 1; i < MAX_ITERATION; i++) {
        getFDerivative(FDerivative, x); 
        getFValue(F, offset, x); 
        solveEquations(FDerivative, F, delta_x); 
        if(getNorm(delta_x)/getNorm(x) <= PRICISION)
            break; 
        else
            for(int i = 1; i <= MAT_SIZE; i++)
                x[i] += delta_x[i]; 
    }
}
 
//计算分片二次插值
void interpolate(const double t[][21], const double u[][21], double z[][21], const int xs, const int ys)
{
    double tt[6] = {0, 0.2, 0.4, 0.6, 0.8, 1}; 
    double uu[6] = {0, 0.4, 0.8, 1.2, 1.6, 2}; 
    double zz[6][6] = {{-0.5, -0.34, 0.14, 0.94, 2.06, 3.5}, 
                     {-0.42, -0.5, -0.26, 0.3, 1.18, 2.38}, 
                     {-0.18, -0.5, -0.5, -0.18, 0.46, 1.42}, 
                     {0.22, -0.34, -0.58, -0.5, -0.1, 0.62}, 
                     {0.78, -0.02, -0.5, -0.66, -0.5, -0.02}, 
                     {1.5, 0.46, -0.26, -0.66, -0.74, -0.5}}; 
    double h = 0.2, tao = 0.4; 
    int n = 5, m = 5; 
    for(int i = 0; i < xs; i++) {
        for(int j = 0; j < ys; j++) {
            int xp = 0, yp = 0; 
            if(t[i][j] <= tt[1]+h/2)
                xp = 1; 
            else if(t[i][j] > tt[n-1]-h/2)
                xp = n-1; 
            else {
                for(int q = 2; q <= n-2; q++)
                    if(t[i][j] > tt[q]-h/2 && t[i][j] <= tt[q]+h/2)
                        xp = q; 
            }
            if(u[i][j] <= uu[1]+tao/2)
                yp = 1; 
            else if(u[i][j] > uu[m-1]-tao/2)
                yp = n-1; 
            else {
                for(int q = 2; q <= m-2; q++)
                    if((u[i][j] > uu[q]-tao/2)&&(u[i][j] <= uu[q]+tao/2))
                        yp = q; 
            }
            z[i][j] = 0; 
            for(int k = xp-1; k <= xp+1; k++) {
                double lk = 1; 
                for(int ti = xp-1; ti <= xp+1; ti++) {
                    if(ti != k)
                        lk *= (t[i][j]-tt[ti])/(tt[k]-tt[ti]); 
                }
                for(int r = yp-1; r <= yp+1; r++) {
                    double lr = 1; 
                    for(int ti = yp-1; ti <= yp+1; ti++) {
                        if(ti != r)
                            lr *= (u[i][j]-uu[ti])/(uu[r]-uu[ti]); 
                    }
                    z[i][j] += lk*lr*zz[k][r]; 
                }
            }
        }
    } 
}
 
//计算两个矩阵的乘积
void multiplyMat(const double* const A, const double* const B, const int r, const int s, const int t, double* const result)
{
    for(int i = 0; i < r; i++) {
        for(int j = 0; j < t; j++) {
            result[i*t+j] = 0; 
            for(int k = 0; k < s; k++)
                result[i*t+j] += A[i*s+k]*B[k*t+j]; 
        }
    }
}
 
//计算矩阵的转置
void transposeMat(const double* const A, const int r, const int s, double* const result)
{
    for(int i = 0; i < r; i++)
        for(int j = 0; j < s; j++)
            result[j*r+i] = A[i*s+j]; 
}
 
//拷贝矩阵
void copyMat(const double* const A, const int r, const int s, const int width, double* const result)
{
    for(int i = 0; i < r; i++)
        for(int j = 0; j < s; j++)
            result[i*s+j] = A[i*width+j]; 
}
 
//计算矩阵的逆
void inverseMat(const double* const A, const int r, double* const result)
{
    double A_copy[r*r]; 
    copyMat(A, r, r, r, A_copy); 
    for(int i = 0; i < r; i++)
        for(int j = 0; j < r; j++) {
            if(i != j)
                result[i*r+j] = 0; 
            else
                result[i*r+j] = 1; 
        }
    for(int i = 0; i < r; i++) {
        int index = i; 
        for(int j = i+1; j < r; j++)
            if(fabs(A_copy[j*r+i]) > fabs(A_copy[index*r+i]))
                index = j; 
 
        if(index != i) {
            for(int j = 0; j < r; j++) {
                double temp = A_copy[index*r+j]; 
                A_copy[index*r+j] = A_copy[i*r+j]; 
                A_copy[i*r+j] = temp; 
                temp = result[index*r+j]; 
                result[index*r+j] = result[i*r+j]; 
                result[i*r+j] = temp; 
            }
        }
        double temp = A_copy[i*r+i]; 
        if(temp != 1) {
            for(int j = 0; j < r; j++) {
                A_copy[i*r+j] /= temp; 
                result[i*r+j] /= temp; 
            }
        }
        for(int j = 0; j < r; j++) {
            if(A_copy[j*r+i] != 0 && i != j) {
                temp = A_copy[j*r+i]; 
                for(int k = 0; k < r; k++) {
                    A_copy[j*r+k] -= temp*A_copy[i*r+k]; 
                    result[j*r+k] -= temp*result[i*r+k]; 
                }
            }
        } 
    }
}
 
//曲面拟合
double *fitSurface(double *z, int &kvalue)
{
    int xs = 11, ys = 21, num = 9; 
    double B[xs*num], G[ys*num], P[xs*ys]; 
    double B_temp[xs*ys], G_temp[xs*ys]; 
    double B_T[xs*ys], G_T[xs*ys]; 
    double BB[xs*ys], GG[xs*ys]; 
	double *C = new double[xs*ys];
    for(int i = 0; i < num; i++) {
        for(int j = 0; j < xs; j++)
            B[j*num+i] = pow(0.08*j, i); 
        for(int j = 0; j < ys; j++)
            G[j*num+i] = pow(0.5+0.05*j, i); 
    }
    double sigma = 0; 
	cout << endl << "选择过程的k和sigma值" << endl;
    for(int i = 0; i < num; i++) {
        sigma = 0; 
		//compute  (B_t * B)-1 * B_t * U * G * (G_t * G)-1
        copyMat(B, xs, i+1, num, B_temp); 
        transposeMat(B_temp, xs, i+1, B_T); 
        multiplyMat(B_T, B_temp, i+1, xs, i+1, BB); 
        inverseMat(BB, i+1, B_temp); 
        copyMat(G, ys, i+1, num, G_temp); 
        transposeMat(G_temp, ys, i+1, G_T); 
        multiplyMat(G_T, G_temp, i+1, ys, i+1, GG); 
        inverseMat(GG, i+1, G_T); 
        multiplyMat(B_temp, B_T, i+1, i+1, xs, BB); 
        multiplyMat(BB, z, i+1, xs, ys, GG); 
        multiplyMat(GG, G_temp, i+1, ys, i+1, BB); 
        multiplyMat(BB, G_T, i+1, i+1, i+1, C); 
        for(int j = 0; j < xs; j++) {
            for(int k = 0; k < ys; k++) {
                double temp = 0; 
                for(int p = 0; p < i+1; p++)
                    for(int q = 0; q < i+1; q++)
                        temp += C[p*(i+1)+q]*B[j*num+p]*G[k*num+q]; 
                P[j*ys+k] = temp; 
                sigma += (z[j*ys+k]-temp)*(z[j*ys+k]-temp); 
            }
        }
        cout << "k = " << i << setprecision(11) << setiosflags(ios::scientific|ios::uppercase) << " sigma = " << sigma << endl; 
        if(sigma < SIGMA) {
            kvalue = i; 
			cout << endl << "达到精度要求时的k和sigma值" << endl;
            cout << "k = " << i << setprecision(11) << setiosflags(ios::scientific|ios::uppercase) << " sigma = " << sigma << endl; 
			cout << endl << "p(x, y)中的系数Crs" << endl;
            for(int k = 0; k <= i; k++) {
                for(int j = 0; j <= i; j++) {
                    cout << "C[" << k << "][" << j << "] = " << C[k*(i+1)+j] << endl; 
                }
            }
            return C; 
        }
    }
    return NULL; 
}
   
int  main()
{
    double t[11][21], u[11][21], z[11][21], offset[MAT_SIZE+1], x[MAT_SIZE+1]; 
    int kvalue = 0; 
    for(int i = 0; i <= 10; i++) {
        for(int j = 0; j <= 20; j++) {
            offset[1] = 2.67+0.08*i; 
            offset[2] = 1.07+0.5+0.05*j; 
            offset[3] = 3.74+0.08*i; 
            offset[4] = 0.79+0.5+0.05*j; 
            Newton(x, offset); 
            t[i][j] = x[1]; 
            u[i][j] = x[2]; 
        }
    }
    interpolate(t, u, z, 11, 21); 
	cout << "数表: (xi, yi, f(xi, yi))" << endl;
    for(int i = 0; i <= 10; i++) {
        for(int j = 0; j <= 20; j++) {
            cout << resetiosflags(ios::scientific|ios::uppercase) << "x = " << 0.08*i << "   y = " << 0.5+0.05*j; 
            cout << setprecision(11) << setiosflags(ios::scientific|ios::uppercase) << "   f(x, y) = " << z[i][j] << endl; 
        }
    }
    double *C = fitSurface(z[0], kvalue); 
    for(int i = 1; i <= 8; i++) {
        for(int j = 1; j <= 5; j++) {
            offset[1] = 2.67+0.1*i; 
            offset[2] = 1.07+0.5+0.2*j; 
            offset[3] = 3.74+0.1*i; 
            offset[4] = 0.79+0.5+0.2*j; 
            Newton(x, offset); 
            t[i-1][j-1] = x[1]; 
            u[i-1][j-1] = x[2]; 
        }
    }
    interpolate(t, u, z, 8, 5); 
    double p[8][5]; 
    for(int i = 1; i <= 8; i++)
        for(int j = 1; j <= 5; j++) {
            double temp = 0; 
            for(int ii = 0; ii <= kvalue; ii++)
                for(int jj = 0; jj <= kvalue; jj++)
                    temp += C[ii*(kvalue+1)+jj]*pow(0.1*i, ii)*pow(0.5+0.2*j, jj); 
            p[i-1][j-1] = temp; 
        }
	cout << endl << "数表: (xi*, yi*, f(xi*, yi*), p(xi*, yi*))" << endl;
    for(int i = 0; i < 8; i++) {
        for(int j = 0; j < 5; j++) {
            cout << resetiosflags(ios::scientific|ios::uppercase); 
            cout << "x[" << i << "] = " << 0.1*i << " y[" << j << "] = " << 0.5+0.2*j << endl; 
            cout << setprecision(11) << setiosflags(ios::scientific|ios::uppercase); 
            cout << "p(x, y) = " << p[i][j] << "  f(x, y) = " << z[i][j] << endl; 
            cout << "delta = " << p[i][j]-z[i][j] << endl; 
        }
    }
    return 0; 
}
