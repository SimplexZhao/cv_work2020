//A class contains Mat Rotate, Multi & Inv.
#include<iomanip>
#include<iostream>
using namespace std;
class MatOperation
{
	public:    
	MatOperation(){};
    virtual ~MatOperation(){};
	
	// Generate R matrix, left=R*right
	// left:image right:world
	// flag==2, φ-ω-κ | flag==1, x-y-z| flag==0, z-y-x
	template<typename T>
	static void GetRotationMat(T* mat, T a1, T a2, T a3, int flag=1)
	{
		if(flag==0)
		{
			T r=a1;
			T p=a2;
			T y=a3;
			mat[0] = cos(r)*cos(y);
			mat[1] = -cos(r)*sin(y);
			mat[2] = -sin(r);
			mat[3] = cos(p)*sin(y)+sin(r)*sin(p)*cos(y);
			mat[4] = cos(p)*cos(y)-sin(r)*sin(p)*sin(y);
			mat[5] = cos(r)*sin(p);
			mat[6] = sin(r)*cos(y)*cos(p)-sin(p)*sin(y);
			mat[7] = -cos(p)*sin(r)*sin(y)-sin(p)*cos(y);
			mat[8] = cos(r)*cos(p);
		}
		if(flag==1)
		{
			T r=a1;
			T p=a2;
			T y=a3;
			mat[0] = cos(p)*cos(y);
			mat[1] = -cos(p)*sin(y);
			mat[2] = -sin(p);
			mat[3] = cos(r)*sin(y)+sin(p)*sin(r)*cos(y);
			mat[4] = cos(r)*cos(y)-sin(p)*sin(r)*sin(y);
			mat[5] = cos(p)*sin(r);
			mat[6] = sin(p)*cos(y)*cos(r)-sin(r)*sin(y);
			mat[7] = -cos(r)*sin(p)*sin(y)-sin(r)*cos(y);
			mat[8] = cos(p)*cos(r);
		}
		if(flag==2)
		{
			T phi=a1;  // a1 b1 c1
			T omega=a2;// a2 b2 c2
			T kappa=a3;// a3 b3 c3
			mat[0] = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa);
			mat[1] = cos(omega)*sin(kappa);
			mat[2] = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa);
			mat[3] = -cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa);
			mat[4] = cos(omega)*cos(kappa);
			mat[5] = -sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa);
			mat[6] = -sin(phi)*cos(omega);
			mat[7] = -sin(omega);
			mat[8] = cos(phi)*cos(omega);
		}
	}

	template<typename T>
	static T getunit(T* v, int len=3)
	{
		T m=0;
		for(int i=0;i<len;++i)m+=v[i]*v[i];
		m=sqrt(m);
		if(m<1e-8)return m;
		for(int i=0;i<len;++i)v[i]/=m;
		return m;
	}

	// 3D vector's cross-product
	template<typename T>
	static void Xprod_3d(T *v1, T *v2, T* res)
	{
		res[0]=v1[1]*v2[2]-v1[2]*v2[1];
		res[1]=-v1[0]*v2[2]+v1[2]*v2[0];
		res[2]=v1[0]*v2[1]-v1[1]*v2[0];
	}

	template<typename T>
	static T Dotprod_nd(T*v1, T*v2, int len)
	{
		if(!(len>0))return -0.2333;
		T result=0;
		for(int i=0;i<len;++i)
		{
			result+=v1[i]*v2[i];
		}
		return result;
	}

	// result=v1-v2
	template<typename T>
	static void Subtraction(T* result, T*v1, T*v2, int len)
	{
		if(!(len>0))return;
		for(int i=0;i<len;++i)
		{
			result[i]=v1[i]-v2[i];
		}
	}

	// Show a m*n matrix; for mx1 vectors n=1 by default
	template<typename T>
	static void ShowMat(T *p, int m, int n=1)
	{
		if(!m*n>0)return;
		int iter=0;
		cout<<endl<<"[";
		cout.flags(ios::fixed);
		cout.precision(6);
		for(int i=0;i<m-1;++i)
		{
			for(int j=0;j<n-1;++j)
			{
				cout<<p[iter]<<", ";
				++iter;
			}
			cout<<p[iter]<<";"<<endl;
			++iter;
		}
		for(int j=0;j<n-1;++j)
		{
			cout<<p[iter]<<", ";
			++iter;
		}
		cout<<p[iter]<<"]"<<endl;
	}

	template<typename T>
	static void Zero(T *p, int n)
	{
		if (n == 0)
			return;
		while (n--) *p++ = 0;
	}
	
	//src m*n, unable to rotate self.
	template<typename T>
	static bool Rotate(T *src, T *dst, int m, int n)
	{
		if (src == 0 || dst == 0)return 0;
		T *ps = src, *pd = dst;
		for (int i = 0; i < m; i++)
		{
			pd = dst + i;
			for (int j = 0; j < n; j++)
			{
				*pd = *ps;
				pd += m;
				++ps;
			}
		}
		pd = ps = 0; return 1;
	}

	///////////////////////////////////////////////
	//设A为m x l阶矩阵，B为l x n阶矩阵，C为m x n阶矩阵，计算 C=A x B的子程序为：
	template<typename T>
	static bool MatrixMulti(T * A, T * B, T * C, int M, int N, int L)
	{
		return 1;
	}

	//对称正定矩阵求逆 a为n*n阶对称正定矩阵，n为矩阵阶数
	template<typename T>
	static int MatrixInversion(T *a, int n)
	{

		return(2);
	}
};