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

public class Genetic extends GAString {
	public double[] ServiceSet_Cost_Temp;
	public static double[] a;
	public String[] Chrom;
	double[] Best;
	int co = 0, co1 = 0;
	public int z;

	void setInitialSequence() {
		Chrom = new String[3000];
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

	public Genetic(int num_chr, int pop_num, double cross, int ge, double mut, int prob) throws GAException {
		super(num_chr * 2, // chromosome dim (number of genes)
				pop_num, // population of chromosomes
				cross, // crossover probability
				6, // random selection chance % (regardless of fitness)
				ge, // stop after this many generations
				0, // num prelim runs (to build good breeding stock for final--full run)
				10, // max prelim generations
				mut, // chromosome mutation prob.
				prob, "01", Crossover.ctRoulette, // crossover type
				true); // compute statistics?
		setInitialSequence();
	}

	static Random rn = new Random();

	protected double getFitness(int iChromIndex) {
		if (co == 0)
			setInitialSequence();
		co1 = 0;
		String s = this.getChromosome(iChromIndex).getGenesAsStr();
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
		Chrom[co1] = hok;
		Best[co1] = Math.abs(cost);
		a[co] = Math.abs(cost);
		co1++;
		co++;
		return Math.abs(cost);
	}

	String Best_Chrom() {
		for (int i = 0; i < co1; i++)
			for (int j = i; j < co1; j++) {
				if (Best[j] > Best[i]) {
					String e = Chrom[i];
					double k = Best[i];
					Chrom[i] = Chrom[j];
					Best[i] = Best[j];
					Chrom[j] = e;
					Best[j] = k;
				}
			}
		return Chrom[0];
	}

	static int y = 0;

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
		removeDuplicates(a);
		double[] temp = new double[y];
		temp = removeDuplicates(a);
		for (int i = 0; i < y - 1; i++)
			for (int j = i; j < y; j++) {
				if (temp[j] > temp[i]) {
					double k = temp[i];
					temp[i] = temp[j];
					temp[j] = k;
				}
			}
		if (y < this.populationDim) {
			for (int i = y; i < populationDim; i++)
				temp[i] = temp[y - 1];
		}
		try {
			File file = new File("fitness.txt");
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			for (int i = 0; i < this.populationDim; i++) {
				bw.write(String.format("%.4f", Math.abs(temp[i])) + "\r\n");
			}
			bw.close();
		} catch (IOException e) {
		}
		return temp;
	}

}
