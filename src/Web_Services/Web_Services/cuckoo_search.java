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

abstract class f_xj {
	abstract double func(double x[]);
}

public class cuckoo_search {
	int n;
	int nd;
	double pa;
	double tol;
	double nest[][];
	double xmin[];
	double xmax[];
	f_xj ff;

//iff为适应度值，由函数func生成，in为种群数量，ipa为，itol为，ixmin为下限数组，ixmax为上限数组
	public cuckoo_search(f_xj iff, int in, double ipa, double itol, double ixmin[], double ixmax[]) {
		n = in;
		nd = ixmin.length;
		pa = ipa;
		tol = itol;
		// nest为二维数组
		nest = new double[n][nd];
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

	// 数组a为随机数，长度与下限长度一致
	double[] randn() {
		double a[] = new double[nd];
		for (int i = 0; i < nd; i++) {
			a[i] = Math.random();
		}
		return a;
	}

	// 数组ones为全为1的数组，长度与下限长度一致
	double[] ones() {
		double v[] = new double[nd];
		for (int i = 0; i < nd; i++) {
			v[i] = 1.0;
		}
		return v;
	}

	// 种群初始化，（上限-下限）*随机数+下限
	double[][][] initialize() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				nest[i][j] = xmin[j] + ((xmax[j] - xmin[j]) * Math.random());
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
				dp[0][i][j] = nest[i][j];
			}
		}
		for (int i = 0; i < n; i++) {
			dp[1][i][0] = fitness[i];
		}
		// 函数返回这个三维数组
		return dp;
	}

	// logGamma返回一个数
	double logGamma(double x) {
		double tmp = (x - 0.5) * Math.log(x + 4.5) - (x + 4.5);
		double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1) + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
				+ 0.00120858003 / (x + 4) - 0.00000536382 / (x + 5);
		return tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
	}

	// gamma返回一个e的次方的数
	double gamma(double x) {
		return Math.exp(logGamma(x));
	}

	// 函数getmival_index参数为一个数组，返回一个长度为2的数组
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
	double[][][] get_best_nest(double[][] nest, double[][] newnest, double[] fitness) {
		// 新种群的适应度值小于原来的适应度值，原来的适应度值更新为新的
		double fnew = 0.0;
		for (int j = 0; j < n; j++) {
			fnew = ff.func(newnest[j]);
			if (fnew < fitness[j]) {
				fitness[j] = fnew;
				nest[j] = newnest[j];
			}
		}
		// 得到fitness数组中的最小值和最小值的下标，fmin为最小值，index为下标
		double[] dep = getminval_index(fitness);
		double fmin = dep[0];
		int index = (int) dep[1];
		// best为一个长度为nd的数组，值为适应度值最低的那一行
		double[] best = new double[nd];
		for (int i = 0; i < nd; i++) {
			best[i] = nest[index][i];
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
		// depp第三层为nest二维数组
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				depp[3][i][j] = nest[i][j];
			}
		}
		return depp;
	}

	// Random permutation of matrix
	// 矩阵的随机置换，一行一行置换
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

	// 返回一个一位数组
	double[][] get_cuckoos(double[][] nest, double best[]) {
		double[] s = new double[nd];
		double[] u = new double[nd];
		double[] v = new double[nd];
		double[] step = new double[nd];
		double[] stepsize = new double[nd];
		// 此处得到一个sigma，为一个指数的数
		double beta = 1.5;
		double sigma = Math.pow((gamma(1 + beta) * Math.sin(3.1415 * beta / 2.0)
				/ (gamma((1.0 + beta) / 2.0) * beta * Math.pow(2.0, (beta - 1.0) / 2.0))), (1.0 / beta));
		double[] sk = new double[nd];
		// rndm为一个随机数
		Random rndm = new Random();
		// 循环所有种群
		for (int j = 0; j < n; j++) {
			s = nest[j];
			for (int i = 0; i < nd; i++)
			// u[i]和v[i]都是符合高斯分布的数组
			{
				u[i] = (rndm.nextGaussian()) * sigma;
				v[i] = (rndm.nextGaussian());
				step[i] = u[i] / Math.pow((Math.abs(v[i])), (1 / beta));
				stepsize[i] = 0.01 * step[i] * (s[i] - best[i]);
			}
			// 得到的s为一位数组
			for (int k = 0; k < nd; k++) {
				sk[k] = stepsize[k] * rndm.nextGaussian();
				s[k] = sk[k] + s[k];
			}
			nest[j] = simplebounds(s);
		}
		return nest;
	}

	// 返回s要么是数组的上限要么是下限
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

	// 得到一个新的矩阵newnest
	double[][] emptynests(double[][] nest) {
		double[][] K = new double[n][nd];
		double[][] randmat = new double[n][nd];
		double[] perm1 = new double[n];
		double[] perm2 = new double[n];
		double[][] nest1 = new double[n][nd];
		double[][] nest2 = new double[n][nd];
		double[][] nest3 = new double[n][nd];
		double[][] nest4 = new double[n][nd];
		double[][] stepsize = new double[n][nd];
		double[][] randmat2 = new double[n][nd];
		double[][] newnest = new double[n][nd];
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
		// nest1和nest2为参数的nest
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				nest1[i][j] = nest[i][j];
				nest2[i][j] = nest[i][j];
			}
		}
		// 稀疏矩阵运算器，得到nest3为二维矩阵
		nest3 = Matrix.substract(randperm(nest1), randperm(nest2));
		// 二维数组stepsize跟nest3有关
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nd; j++) {
				stepsize[i][j] = nest3[i][j] * Math.random();
			}
		}
		// 矩阵stepsize和矩阵K点积得到nest4
		nest4 = Matrix.dotproduct(stepsize, K);
		// newnest跟矩阵nest和nest4有关
		newnest = Matrix.add(nest, nest4);
		return newnest;
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
		double d2[][][] = get_best_nest(nestt, nestt, fitness);
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
		// nestt得到一个新的矩阵nestt
		nestt = emptynests(nestt);
		// System.out.println(Matrix.toString(nestt));
		// 当最小的适应度值fmin大于tol，tol为参数
		while (fmin > tol) { // newnest为一维数组
			newnest = get_cuckoos(nest, bestnest);
			// d3为三维数组，并且fnew为最小的适应度值
			double d3[][][] = get_best_nest(nestt, newnest, fitness);
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
			iter += n;
			// 得到一个新的矩阵newnest
			newnest = emptynests(nestt);
			// System.out.println(Matrix.toString(nestt));
			double d4[][][] = get_best_nest(nestt, newnest, fitness);
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
			iter += n;
			// 选出最小的适应度值，最好的种群为bestnest
			if (fnew < fmin) {
				fmin = fnew;
				bestnest = best;
			}
		}
		for (int i = 0; i < bestnest.length; i++) {
			System.out.println("x[" + i + "] = " + bestnest[i]);
		}
	}
}