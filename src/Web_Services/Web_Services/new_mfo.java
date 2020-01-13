package Web_Services;


import static Web_Services.Ga_Service_Selection.Cost_Params;
import static Web_Services.Ga_Service_Selection.Provider_Count;
import static Web_Services.Ga_Service_Selection.Provider_Max_Service;
import static Web_Services.Ga_Service_Selection.Providers_Resource;
import static Web_Services.Ga_Service_Selection.ServiceSet_Capability_Matrix;
import static Web_Services.Ga_Service_Selection.Resource_Cost;
import static Web_Services.Ga_Service_Selection.Weight_cost;
import static Web_Services.Ga_Service_Selection.Base_Option_Count;
import com.softtechdesign.ga.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import org.cloudbus.cloudsim.Log;


//abstract class f_xj {
//	abstract double func(double x[]);
//}
public class new_mfo {
	int n;
	int nd;
	int max_iter;
	double pa;
	double tol;
	//moth为飞蛾，flame为火焰
	double flame[][];
	double moth[][];
	double xmin[];
	double xmax[];
	//类f_xj生成一个实例ff，调用其中的func方法
	f_xj ff;

//iff为适应度值，由函数func生成，in为种群数量，ipa为，itol为，ixmin为下限数组，ixmax为上限数组
	// imax_iter为最大的迭代次数
	public new_mfo(f_xj iff, int in, double ipa, double itol, double ixmin[], double ixmax[],int imax_iter) {
	//n为矩阵的行数，nd为列数
		n = in;
		nd = ixmin.length;
		pa = ipa;
		tol = itol;
		max_iter = imax_iter;
		// flame为二维数组
		flame = new double[n][nd];
		moth = new double[n][nd];
		xmin = new double[nd];
		xmax = new double[nd];
		// ff为适应度值
		ff = iff;
//		double[] one = ones();
		for (int i = 0; i < nd; i++) {
			xmin[i] = ixmin[i];
			xmax[i] = ixmax[i];
		}
		System.out.println("hello mew_mfo");

	}

	// 数组a为随机数，长度与下限长度一致，为列数
	double[] randn() {
		double a[] = new double[nd];
		for (int i = 0; i < nd; i++) {
			a[i] = Math.random();
		}
		return a;
	}

	// 数组ones为全为1的数组，长度与下限长度一致，长度与列数一致
	double[] ones() {
		double v[] = new double[nd];
		for (int i = 0; i < nd; i++) {
			v[i] = 1.0;
		}
		return v;
	}

	// 种群初始化，（上限-下限）*随机数+下限，返回一个三维数组
	double[][][] initialize() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				moth[i][j] = xmin[j] + ((xmax[j] - xmin[j]) * Math.random());
			}
		}
		double[] fitness = new double[n];
		// 适应度值数组初始化为10的7次方
		for (int i = 0; i < n; i++) {
			fitness[i] = Math.pow(10, 7);
		}
		// dp为一个三维数组，第一层为初始化的种群，第二层为适应度值
		double dp[][][] = new double[2][n][nd];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				dp[0][i][j] = moth[i][j];
			}
		}
		for (int i = 0; i < n; i++) {
			dp[1][i][0] = fitness[i];
		}
		// 函数返回这个三维数组
		return dp;
	}



	// 函数getmival_index参数为一个数组，返回一个长度为2的数组
	//dep[0]存储最小值，dep[1]存储最小值的位置
	double[] getminval_index(double[] a) {
		double m = 0.0;
		double b[] = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			b[i] = a[i];
		}
		// 将minval复制为数组中的最小值
		double minval = a[0];
		for (int i = 0; i < a.length; i++) {
			if (a[i] < minval) {
				minval = a[i];
			}
		}
		// m为最小值的位置
		for (int i = 0; i < a.length; i++) {
			if (b[i] == minval) {
				m = i;
				break;
			}
		}
		;
		// dep为长度为2的数组，dep[0]存储最小值，dep[1]存储最小值的位置
		double[] dep = new double[2];
		dep[0] = minval;
		dep[1] = m;
		return dep;
	}

	// 返回为一个三维数组
	double[][][] get_best_flame(double[][] moth, double[][] flame, double[] fitness) {
		// 新种群的适应度值小于原来的适应度值，原来的适应度值更新为新的
		double fnew = 0.0;
		for (int j = 0; j < n; j++) {
			fnew = ff.func(flame[j]);
			if (fnew < fitness[j]) {
				fitness[j] = fnew;
				moth[j] = flame[j];
			}
		}
		// 得到fitness数组中的最小值和最小值的下标，fmin为最小值，index为下标
		double[] dep = getminval_index(fitness);
		double fmin = dep[0];
		int index = (int) dep[1];
		// best为一个长度为nd的数组，值为适应度值最低的那一行
		double[] best = new double[nd];
		for (int i = 0; i < nd; i++) {
			best[i] = moth[index][i];
		}
		// depp为三维数组，depp[0][0][0]为fitness数组中的最小值
		double depp[][][] = new double[4][n][nd];
		depp[0][0][0] = fmin;
		// depp第一层为fitness数组
		for (int i = 0; i < n; i++) {
			depp[1][i][0] = fitness[i];
		}
		// depp第二层为适应度值最好的一组，变成火焰
		for (int i = 0; i < nd; i++) {
			depp[2][0][i] = best[i];
		}
		// depp第三层为飞蛾
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				depp[3][i][j] = moth[i][j];
			}
		}
		return depp;
	}

	// Random permutation of matrix
	// 矩阵的随机置换，一行一行置换，顺序可以打得很乱
	double[][] randperm(double a[][]) {
		double y[][] = new double[a.length][a[0].length];
		double[] sw = new double[a[0].length];
		for (int i = 0; i < a.length; i++) {
			y[i] = a[i];
		}
		int m = a.length;
		int rp = 0;
		for (int i = 0; i < m - 1; i++) {
			rp = (int) ((m - i) * Math.random() + i);
			sw = y[i];
			y[i] = y[rp];
			y[rp] = sw;
		}
		return y;
	}

	// 得到飞蛾的位置更新，返回一个一维数组，并进行位置的修正
	double[][] get_moths(double[][] moth, double flame[],int iter,int max_iter) {
		double[] s = new double[nd];	
		//b为常数1
		double b=1.0;
		//t为[-1，1]之间的随机数
		double t = Math.random()*2-1;
		
		// 循环所有种群
		for (int j = 0; j < n; j++) {
			s = moth[j];
			for (int i = 0; i < nd; i++)
			{
				//加入权重w
				double w = Math.sin((Math.PI*iter)/(2*max_iter)+Math.PI)+1;
				//螺旋线函数更新公式
				s[i]=Math.abs((s[i] - flame[i]))*Math.exp(b*t)*Math.cos(2*Math.PI*t)+w*flame[i];
			}
	//修正位置，将位置固定在上限和下限中
			moth[j] = simplebounds(s);
		}
		return moth;
	}

	// 返回s,防止越界
	double[] simplebounds(double s[]) {
		for (int i = 0; i < nd; i++) {
			if (s[i] < xmin[i]) {
				s[i] = xmin[i];
			}
			if (s[i] > xmax[i]) {
				s[i] = xmax[i];
			}
		}
		return s;
	}

	// 飞蛾随机迁移，得到一个新的矩阵newmoth
	double[][] emptynests(double[][] moth) {
		double[][] K = new double[n][nd];
		double[][] randmat = new double[n][nd];
		double[][] moth1 = new double[n][nd];
		double[][] moth2 = new double[n][nd];
		double[][] moth3 = new double[n][nd];
		double[][] moth4 = new double[n][nd];
		double[][] stepsize = new double[n][nd];
		double[][] newmoth = new double[n][nd];
		// 二维数组K要么取0，要么取1
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				randmat[i][j] = Math.random();
				if (randmat[i][j] < pa) {
					K[i][j] = 0.0;
				}
				if (randmat[i][j] > pa) {
					K[i][j] = 1.0;
				}
			}
		}
		// moth1和moth2均为参数的值
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				moth1[i][j] = moth[i][j];
				moth2[i][j] = moth[i][j];
			}
		}
		// moth1和moth2随机置换之后，moth1-moth2
		moth3 = Matrix.substract(randperm(moth1), randperm(moth2));
		// 二维数组stepsize跟moth3有关
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				stepsize[i][j] = moth3[i][j] * Math.random();
			}
		}
		// 矩阵stepsize和矩阵K点积得到moth4
		moth4 = Matrix.dotproduct(stepsize, K);
		// 矩阵newmoth为原矩阵moth和moth4相加
		newmoth = Matrix.add(flame, moth4);
		return newmoth;
	}

	void solution() {
		double moth[][]=new double[n][nd];
		double flame[]=new double[n];
		double fitness[] = new double[n];
		
//		double bestnest[] = new double[nd];
		double best[] = new double[nd];
		double newmoth[][] = new double[n][nd];
		// 数组d1为初始化之后的矩阵
		double d1[][][] = initialize();
		// moth为d1的第0层
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				moth[i][j] = d1[0][i][j];
			}
		}
		// fitness为d1的第一层
		for (int i = 0; i < n; i++) {
			fitness[i] = d1[1][i][0];
		}
		// d2为一个三维矩阵
		double d2[][][] = get_best_flame(moth, moth, fitness);
		// fmin为fitness中的最小值
		double fmin = d2[0][0][0];
		// 赋值fitness，bestnest，nestt
		for (int i = 0; i < n; i++) {
			fitness[i] = d2[1][i][0];
		}
		for (int i = 0; i < nd; i++) {
			flame[i] = d2[2][0][i];
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
			moth[i][j] = d2[3][i][j];
			}
		}
		// iter为迭代步数
		int iter = 0;
		double fnew = 0.0;
		// 得到一个新的矩阵moth矩阵
		moth = emptynests(moth);
		// System.out.println(Matrix.toString(nestt));
		// 当最小的适应度值fmin大于tol，tol为参数
		//并且当前的迭代次数要小于最大的迭代次数
		while (fmin > tol&&iter<max_iter) { 
			//飞蛾进行螺旋线飞行
			newmoth = get_moths(moth, flame,iter,max_iter);
			// d3为三维数组，并且fnew为最小的适应度值
			double d3[][][] = get_best_flame(moth,newmoth, fitness);
			fnew = d3[0][0][0];
			//
			for (int i = 0; i < n; i++) {
				fitness[i] = d3[1][i][0];
			}
			//
			for (int i = 0; i < nd; i++) {
				best[i] = d3[2][0][i];
			}
			//
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < nd; j++) {
					moth[i][j] = d3[3][i][j];
				}
			}
			// 随机迁移，得到一个新的矩阵moth
			moth = emptynests(moth);
			// System.out.println(Matrix.toString(nestt));
			double d4[][][] = get_best_flame(moth, newmoth, fitness);
			// 更新之后的最小适应度值为fnew
			fnew = d4[0][0][0];
			for (int i = 0; i < n; i++) {
				fitness[i] = d4[1][i][0];
			}
			for (int i = 0; i < nd; i++) {
				best[i] = d4[2][0][i];
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < nd; j++) {
					moth[i][j] = d4[3][i][j];
				}
			}
			iter += 1;
			// 选出最小的适应度值，最好的种群为flame
			if (fnew < fmin) {
				fmin = fnew;
				flame = best;
			}
		}
		for (int i = 0; i < flame.length; i++) {
			System.out.println("x[" + i + "] = " + flame[i]);
		}
	}

}
