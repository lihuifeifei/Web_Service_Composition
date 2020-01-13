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
public class mfo {
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
	public mfo(f_xj iff, int in, double itol, double ixmin[], double ixmax[],int imax_iter) {
	//n为矩阵的行数，nd为列数
		n = in;
		nd = ixmin.length;
//		pa = ipa;
		tol = itol;
		max_iter = imax_iter;
		// flame为二维数组
		flame = new double[n][nd];
		moth = new double[n][nd];
		xmin = new double[nd];
		xmax = new double[nd];
		// ff为适应度值
		ff = iff;
		double[] one = ones();
		for (int i = 0; i < nd; i++) {
			xmin[i] = ixmin[i];
			xmax[i] = ixmax[i];
		}
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
	double[][][] get_best_flame(double[][] flame, double[][] newflame, double[] fitness) {
		// 新种群的适应度值小于原来的适应度值，原来的适应度值更新为新的
		double fnew = 0.0;
		for (int j = 0; j < n; j++) {
			fnew = ff.func(newflame[j]);
			if (fnew < fitness[j]) {
				fitness[j] = fnew;
				flame[j] = newflame[j];
			}
		}
		// 得到fitness数组中的最小值和最小值的下标，fmin为最小值，index为下标
		double[] dep = getminval_index(fitness);
		double fmin = dep[0];
		int index = (int) dep[1];
		// best为一个长度为列数的数组，值为适应度值最低的那一行
		double[] best = new double[nd];
		for (int i = 0; i < nd; i++) {
			best[i] = flame[index][i];
		}
		// depp为三维数组，depp[0][0][0]为fitness数组中的最小值
		double depp[][][] = new double[4][n][nd];
		depp[0][0][0] = fmin;
		// depp第一层为fitness数组
		for (int i = 0; i < n; i++) {
			depp[1][i][0] = fitness[i];
		}
		// depp第二层为best数组
		for (int i = 0; i < nd; i++) {
			depp[2][0][i] = best[i];
		}
		// depp第三层为flame二维数组
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				depp[3][i][j] = flame[i][j];
			}
		}
		return depp;
	}


	// 得到最好的飞蛾位置，返回一个一维数组
	double[][] get_moths(double[][] flame, double best[]) {
		double[] s = new double[nd];	
		//b为常数1
		double b=1.0;
		//t为[-1，1]之间的随机数
		double t = Math.random()*2-1;
		
		// 循环所有种群
		for (int j = 0; j < n; j++) {
			s = flame[j];
			for (int i = 0; i < nd; i++)
			{
				//螺旋线函数更新公式
				s[i]=Math.abs((s[i] - best[i]))*Math.exp(b*t)*Math.cos(2*Math.PI*t)+best[i];
			}
			flame[j] = s;
		}
		return flame;
	}



	void solution() {
		double nestt[][] = new double[n][nd];
		double fitness[] = new double[n];
		double bestnest[] = new double[nd];
		double best[] = new double[nd];
		double newnest[][] = new double[n][nd];
		// 数组d1为初始化之后的矩阵
		double d1[][][] = initialize();
		// nestt为d1的第0层
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				nestt[i][j] = d1[0][i][j];
			}
		}
		// fitness为d1的第一层
		for (int i = 0; i < n; i++) {
			fitness[i] = d1[1][i][0];
		}
		// d2为一个三维矩阵
		double d2[][][] = get_best_flame(nestt, nestt, fitness);
		// fmin为fitness中的最小值
		double fmin = d2[0][0][0];
		// 赋值fitness，bestnest，nestt
		for (int i = 0; i < n; i++) {
			fitness[i] = d2[1][i][0];
		}
		for (int i = 0; i < nd; i++) {
			bestnest[i] = d2[2][0][i];
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				nestt[i][j] = d2[3][i][j];
			}
		}
		// iter为迭代步数
		int iter = 0;
		double fnew = 0.0;
		// System.out.println(Matrix.toString(nestt));
		// 当最小的适应度值fmin大于tol，tol为参数
		//并且当前的迭代次数要小于最大的迭代次数
		while (fmin > tol&&iter<max_iter) { // newnest为一维数组
			newnest = get_moths(flame, bestnest);
			// d3为三维数组，并且fnew为最小的适应度值
			double d3[][][] = get_best_flame(nestt, newnest, fitness);
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
					nestt[i][j] = d3[3][i][j];
				}
			}
			// n为参数，
			iter += 1;
			// System.out.println(Matrix.toString(nestt));
			double d4[][][] = get_best_flame(nestt, newnest, fitness);
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
					nestt[i][j] = d4[3][i][j];
				}
			}
			iter += 1;
			// 选出最小的适应度值，最好的种群为bestnest
			if (fnew < fmin) {
				fmin = fnew;
				bestnest = best;
			}
		}
		System.out.println("最小适应度值为：Optimum value of a function = " + fmin);
		for (int i = 0; i < bestnest.length; i++) {
			System.out.println("x[" + i + "] = " + bestnest[i]);
		}
	}

}
