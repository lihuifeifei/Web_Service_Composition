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

class funfit extends f_xj {
	public double[] ServiceSet_Cost_Temp;
	public static double[] a;
	public String[] cuckoo;
	double[] Best;
	int co = 0, co1 = 0;
	public int z;
	long best_cuckoo = 3000;
	long current_cuckoo = 18500;

//初始化序列
	void setInitialSequence() {
		cuckoo = new String[3000];
		a = new double[3000];
		z = 0;
		ServiceSet_Cost_Temp = new double[20000];
		Weight_cost = 10;
		for (int i = 0; i < Provider_Count; i++) {
			for (int re = 0; re < Providers_Resource[i]; re++) {
				for (int k = 0; k < 4; k++) {
					ServiceSet_Cost_Temp[z] += ServiceSet_Capability_Matrix[i][re][k];

				}
				z++;
			}
		}
		Best = new double[3000];
		Base_Option_Count = 7;
	}

	static int y = 0;
	// rn为随机数
	static Random rn = new Random();

	double func(double x[]) {
		setInitialSequence();
		for (int i = 0; i < 50; i++)
			a[i] = rn.nextDouble() * (best_cuckoo) + current_cuckoo;
		co1 = 0;
		String s = "";
		int counter = 0;
		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i) == '1') {
				counter++;
			}
		}
		String hok = "";
		double cost = 0;
		for (int i = 0; i < s.length(); i++) {
			int o = rn.nextInt(z);
			cost += ServiceSet_Cost_Temp[o];
			hok += o;
		}
		cuckoo[co1] = hok;
		Best[co1] = Math.abs(cost);
		co1++;
		co++;
		return Math.abs(cost);
	}

	public static double[] removeDuplicates(double[] arr) {
		Set<Double> alreadyPresent = new HashSet<>();
		double[] whitelist = new double[arr.length];
		int i = 0;
		for (double element : arr) {
			if (alreadyPresent.add(element)) {
				whitelist[i++] = element;
			}
		}
		y = i;
		return Arrays.copyOf(whitelist, i);
	}

	public double[] get_all_fit() {
		double[] temp = new double[50];
		temp = a;
		// 冒泡排序，temp最后为一个从大到小的数组
		for (int i = 0; i < 50 - 1; i++)
			for (int j = i; j < 50; j++) {
				if (temp[j] > temp[i]) {
					double k = temp[i];
					temp[i] = temp[j];
					temp[j] = k;
				}
			}
		try {
			File file = new File("fitness.txt");
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			for (int i = 0; i < 50; i++) {
//                    		System.out.println("适应度值为:"+String.format("%.4f",Math.abs(temp[i])));		

				bw.write(String.format("%.4f", Math.abs(temp[i])) + "\r\n");
			}
			bw.close();
		} catch (IOException e) {
		}
		return temp;
	}

}
