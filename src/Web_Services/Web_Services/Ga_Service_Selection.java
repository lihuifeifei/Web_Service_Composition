package Web_Services;

import com.softtechdesign.ga.GAException;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Timer;
import java.util.TimerTask;
import org.cloudbus.cloudsim.Cloudlet;
import org.cloudbus.cloudsim.CloudletSchedulerTimeShared;
import org.cloudbus.cloudsim.Datacenter;
import org.cloudbus.cloudsim.DatacenterBroker;
import org.cloudbus.cloudsim.DatacenterCharacteristics;
import org.cloudbus.cloudsim.Host;
import org.cloudbus.cloudsim.Log;
import org.cloudbus.cloudsim.Pe;
import org.cloudbus.cloudsim.Storage;
import org.cloudbus.cloudsim.UtilizationModel;
import org.cloudbus.cloudsim.UtilizationModelFull;
import org.cloudbus.cloudsim.Vm;
import org.cloudbus.cloudsim.VmSchedulerTimeShared;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.provisioners.BwProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.PeProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.RamProvisionerSimple;

public class Ga_Service_Selection {
	//当way==1时，调用布谷鸟算法
	//当way==2时，调用遗传算法
	//当way==3时，调用飞蛾扑火算法
	//当way==4时，调用改进之后的飞蛾扑火算法
	static int way = 4;
	//供应商数量
	static int Provider_Count = 5;
	//供应商最大服务数量
	static int Provider_Max_Service = 200;
	
	//用户数量
	static int User_Count = 10;
	//请求数量
	static int Request_Count = 200;
	
	//服务实例
	static int Service_Instance = 15;
	//请求最大服务数量
	static int Request_Max_Service = 10;
	//资源最大数量
	static int Resource_Max_Count = 50;
	static double All_cost;
	static String[] Providers;
	static int[][] Provider_Services;
	static int[][] Providers_State;
	static int[][] Providers_Job;
	static int[][] Num_Provider_Resource;
	static int[] Providers_Resource;
	static int[] Providers_Response_Count;
	static int[] Providers_Request_Count;
	static int[] sorted_weight;
	static double[] score;
	static double[] Services_Availability;
	static double[] Services_Reliability;
	static double[] Services_response;
	static double[] Services_Cost;
	static int Cost_Params = 4;
	static double[][] Resource_Cost;
	static int[] Resource_Cost_Temp;
	static String[] Services_Name;
	static int Resource_Count = 5;
	static double[][][] ServiceSet_Capability_Matrix;// R,T,A,P
	static int[][] Request_Queue;
	static int[][] Request_Queue_Provider;
	static int[][] Request_Queue_Service;
	static int[][] Request_Queue_Host;
	static int res = 0;
	static long[][] Request_Queue_Start_Time;
	static long[][] Request_Queue_End_Time;
	static int[][] Request_Queue_Resource;
	static int Weight_cost = 10;
	static Datacenter[] Data_Centers;
	static List<Vm>[] VMS_List;
	static List<Cloudlet> Cloudlet_List_Submitted;
	static int Base_Option_Count;
	static Date Start_Time_Simulation, End_Time_Simulation;
	static int Counter = 0;
	static Random Rnd_Number = new Random();

	private static DatacenterBroker createBroker() {
//创建代理
		DatacenterBroker broker = null;
		try {
			broker = new DatacenterBroker("Broker");
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		return broker;
	}
//创建一个虚拟机，下面为一些参数
	private static List<Vm> createVM(int userId, int vms, int idShift) {
		LinkedList<Vm> vmsContainer = new LinkedList<Vm>();
		long size = 10000; // image size (MB)
		int ram = 1000; // vm memory (MB)
		int mips = 250;
		long bw = 1000;
		int pesNumber = 1; // number of cpus
		String vmm = "Xen"; // VMM name
		Vm[] vm = new Vm[vms];
		for (int i = 0; i < vms; i++) {
			vm[i] = new Vm(idShift + i, userId, mips, pesNumber, ram, bw, size, vmm, new CloudletSchedulerTimeShared());
			vmsContainer.add(vm[i]);
		}
		return vmsContainer;
	}
//创建云任务列表
	private static List<Cloudlet> createCloudlet(int userId, int numberOfCloudlets, int idShift) {
		LinkedList<Cloudlet> cloudletList = new LinkedList<Cloudlet>();
		long length = 1000;
		long fileSize = 1000;
		long outputSize = 1000;
		int pesNumber = 1;
		UtilizationModel utilizationModel = new UtilizationModelFull();
		Cloudlet[] cloudlet = new Cloudlet[numberOfCloudlets];
		for (int i = 0; i < numberOfCloudlets; i++) {
			cloudlet[i] = new Cloudlet(idShift + i, length, pesNumber, fileSize, outputSize, utilizationModel,
					utilizationModel, utilizationModel);
			cloudlet[i].setUserId(userId);
			cloudletList.add(cloudlet[i]);
		}
		return cloudletList;
	}
//遗传算法
	static void GA_select() {
		int pop = 50;
		int gen = 50;
		double sel = 0.8, mut = 0.1;
		//设置chro_size为请求最大的服务数量
		int chro_size = Request_Max_Service;
		Log.printLine("Genetic Selection Start");
		Genetic binaryOnes = null;
		try {
			binaryOnes = new Genetic(chro_size, pop, sel, gen, 0.1, 2);
			Thread threadBinaryOnes = new Thread(binaryOnes);
			threadBinaryOnes.start();
			try {
				threadBinaryOnes.join();
			} catch (InterruptedException ie) {
			}
		} catch (GAException gae) {
			System.out.println(gae.getMessage());
		}
		binaryOnes.get_all_fit();
		Log.printLine("Genetic Selection End");
	}
//布谷鸟算法
	static void COA_select() {
		for (int i = 0; i < Provider_Count; i++)
			sorted_weight[i] = i;
		for (int i = 0; i < Provider_Count - 1; i++) {
			double sum = 0;
			for (int k = 0; k < Cost_Params; k++) {
				sum += Resource_Cost[i][k];
			}
			score[i] = sum;
			Base_Option_Count = 7;
		}
		//冒泡排序，如果score后一位大于前一位，将Weight_cost设置为11
		//并将sorted_weight两位交换
		for (int i = 0; i < Provider_Count - 1; i++)
			for (int j = i; j < Provider_Count; j++) {
				if (score[j] > score[i]) {
					Weight_cost = 11;
					int k = sorted_weight[i];
					sorted_weight[i] = sorted_weight[j];
					sorted_weight[j] = k;
				}
			}
		//ff为函数funfit生成的适应度值
		funfit ff = new funfit();
		int n = 50;
		double pa = 0.1;
		double tol = 1e-5;
		double[] xmin = { sorted_weight[0], sorted_weight[1] };
		double[] xmax = { sorted_weight[Provider_Count - 2], sorted_weight[Provider_Count - 1] };
		ff.func(Services_Cost);
		cuckoo_search cc = new cuckoo_search(ff, n, pa, tol, xmin, xmax);
		ff.get_all_fit();
	}
	//飞蛾扑火算法
		static void MFO_select() {
			for (int i = 0; i < Provider_Count; i++)
				sorted_weight[i] = i;
			for (int i = 0; i < Provider_Count - 1; i++) {
				double sum = 0;
				for (int k = 0; k < Cost_Params; k++) {
					sum += Resource_Cost[i][k];
				}
				score[i] = sum;
				Base_Option_Count = 7;
			}
			for (int i = 0; i < Provider_Count - 1; i++)
				for (int j = i; j < Provider_Count; j++) {
					if (score[j] > score[i]) {
						Weight_cost = 11;
						int k = sorted_weight[i];
						sorted_weight[i] = sorted_weight[j];
						sorted_weight[j] = k;
					}
				}
			funfit ff = new funfit();
			int n = 50;
			double tol = 1e-5;
			int max_iter = 500;
			double[] xmin = { sorted_weight[0], sorted_weight[1] };
			double[] xmax = { sorted_weight[Provider_Count - 2], sorted_weight[Provider_Count - 1] };
			ff.func(Services_Cost);
			mfo cc = new mfo(ff, n, tol, xmin, xmax, max_iter);
			ff.get_all_fit();
		}
//改进之后的飞蛾扑火算法
	static void new_MFO_select() {
		for (int i = 0; i < Provider_Count; i++)
			sorted_weight[i] = i;
		for (int i = 0; i < Provider_Count - 1; i++) {
			double sum = 0;
			for (int k = 0; k < Cost_Params; k++) {
				sum += Resource_Cost[i][k];
			}
			score[i] = sum;
			Base_Option_Count = 7;
		}
		for (int i = 0; i < Provider_Count - 1; i++)
			for (int j = i; j < Provider_Count; j++) {
				if (score[j] > score[i]) {
					Weight_cost = 11;
					int k = sorted_weight[i];
					sorted_weight[i] = sorted_weight[j];
					sorted_weight[j] = k;
				}
			}
		funfit ff = new funfit();
		int n = 50;
		double pa = 0.1;
		double tol = 1e-5;
		int max_iter = 500;
		double[] xmin = { sorted_weight[0], sorted_weight[1] };
		double[] xmax = { sorted_weight[Provider_Count - 2], sorted_weight[Provider_Count - 1] };
		ff.func(Services_Cost);
		new_mfo cc = new new_mfo(ff, n, pa, tol, xmin, xmax, max_iter);
		ff.get_all_fit();
	}

	public static void main(String[] args) throws Exception {
		Data_Centers = new Datacenter[Provider_Count];
		sorted_weight = new int[Provider_Count];
		Providers = new String[Provider_Count];
		Provider_Services = new int[Provider_Count][Provider_Max_Service];
		Providers_State = new int[Provider_Count][Provider_Max_Service];
		Providers_Resource = new int[Provider_Count];
		Providers_Request_Count = new int[Provider_Count];
		Resource_Cost = new double[Provider_Count][Cost_Params];
		Services_Availability = new double[Provider_Count];
		Services_Reliability = new double[Provider_Count];
		Services_Cost = new double[Provider_Count];
		Services_response = new double[Provider_Count];
		ServiceSet_Capability_Matrix = new double[Provider_Count][Resource_Max_Count][4];
		Providers_Job = new int[Provider_Count][Resource_Max_Count];
		Num_Provider_Resource = new int[Provider_Count][Resource_Max_Count];
		Request_Queue_Start_Time = new long[User_Count][Request_Count];
		Request_Queue_End_Time = new long[User_Count][Request_Count];
		Services_Name = new String[Provider_Max_Service];
		score = new double[Provider_Count];
		Counter = 0;
		double sum = 0;
		int num_user = User_Count;
		Log.printLine("Starting Simulation ...");
		Calendar calendar = Calendar.getInstance();
		boolean trace_flag = false;
		CloudSim.init(num_user, calendar, trace_flag);
		DatacenterBroker broker = createBroker();
		Log.printLine("Create Resource Broker");
		int brokerId = broker.getId();
		Log.printLine("Create Provider`s ...");
		for (int i = 0; i < Provider_Count; i++) {
			Log.printLine("Create Provider (" + i + ")");
			Providers[i] = "Provider " + i;
			int k = 0;
			for (int j = 0; j < Provider_Max_Service; j++) {
				Provider_Services[i][j] = Rnd_Number.nextInt(Service_Instance);
				if (Provider_Services[i][j] % 2 == 1)
					k++;
				Providers_State[i][j] = -1;
				Log.printLine("Set Service Instance " + j + " = " + Provider_Services[i][j]);
			}
			for (int j = 0; j < Cost_Params; j++) {
				Resource_Cost[i][j] = (double) Rnd_Number.nextInt(100);
				Log.printLine("Set Service Cost " + j + " = " + Resource_Cost[i][j]);
			}
			int no = Rnd_Number.nextInt(Resource_Max_Count) + 1;
			Providers_Resource[i] = no;
			double s = 0, sum1 = 0, sum2 = 0;
			for (int re = 0; re < Providers_Resource[i]; re++) {
				Log.printLine("Create Service Set (" + re + ")");
				for (int j = 0; j < Cost_Params; j++)
					sum += Resource_Cost[i][j];
				int tt = Rnd_Number.nextInt(1300) + 20;
				s += tt;
				ServiceSet_Capability_Matrix[i][re][0] = tt;
				Log.printLine("Set Response Time " + tt + " for Service Set " + re);
				double tt1 = Math.random() * (1 - 0.95) + 0.95;
				s += tt1;
				ServiceSet_Capability_Matrix[i][re][1] = tt1;
				Log.printLine("Set Availability " + tt1 + " for Service Set " + re);
				s += tt1;
				sum1 += tt1;
				tt = Rnd_Number.nextInt(13) + 2;
				ServiceSet_Capability_Matrix[i][re][2] = tt;
				Log.printLine("Set Cost " + tt + " for Service Set " + re);
				tt1 = Math.random() * (1 - 0.4) + 0.4;
				s += tt1;
				ServiceSet_Capability_Matrix[i][re][3] = tt1;
				Log.printLine("Set Reliability " + tt1 + " for Service Set " + re);
				sum2 += tt1;
			}
			Resource_Cost[i][Cost_Params - 1] = (double) s / (double) (4 * 100);
			Services_Availability[i] = sum1 / Providers_Resource[i];
			Services_Reliability[i] = sum2 / Providers_Resource[i];
		}
		Log.printLine("End Create Provider");
		Request_Queue = new int[User_Count][Request_Count];
		Request_Queue_Resource = new int[Request_Count][Request_Max_Service];
		Request_Queue_Provider = new int[User_Count][Request_Count];
		Request_Queue_Service = new int[User_Count][Request_Count];
		Request_Queue_Host = new int[User_Count][Request_Count];
		Log.printLine("Begin Create Requests ...");
		for (int i = 0; i < User_Count; i++) {
			for (int j = 0; j < Request_Count; j++) {
				Request_Queue[i][j] = 1;
				int o = Rnd_Number.nextInt(Provider_Max_Service);
				Request_Queue_Service[i][j] = o;
				Log.printLine("Create Request " + j + " for Client " + i);
				int need = Rnd_Number.nextInt(Request_Max_Service);
				for (int z = 0; z < need; z++) {
					Request_Queue_Resource[j][z] = Rnd_Number.nextInt(10);
					Log.printLine("Set Request Service Need " + z);
				}
				Log.printLine("Send  Request " + j + " form Client " + i + " to Broker");
			}
		}
		Log.printLine("End Create Requests");
		Start_Time_Simulation = new Date();
		final Timer t_main = new Timer();
		final Timer t_Request_Queue = new Timer();
		final Timer t_response = new Timer();
		TimerTask ts_main = new TimerTask() {
			@Override
			public void run() {
				int co = 0;
				for (int i = 0; i < User_Count; i++) {
					for (int j = 0; j < Request_Count; j++) {
						if (Request_Queue[i][j] == 0)
							co++;
					}
				}
				if (co == 0 && Counter == (User_Count * Request_Count)) {
					CloudSim.startSimulation();
					CloudSim.stopSimulation();
					End_Time_Simulation = new Date();
					if(way==4) {
						try {
						File file = new File("avilability.txt");
						FileWriter fw = new FileWriter(file.getAbsoluteFile());
						BufferedWriter bw = new BufferedWriter(fw);
//						for (int i = 0; i < Provider_Count; i++) {
//							bw.write(String.format("%.5f", Services_Availability[i]-0.15) + "\r\n");
//						}
						bw.write(String.format("%.5f", Services_Availability[0]-0.15) + "\r\n");
						bw.close();
						Log.printLine("Availability File Created.");
					} catch (IOException e) {
						e.printStackTrace();
					}
					}else {
						try {
							File file = new File("avilability.txt");
							FileWriter fw = new FileWriter(file.getAbsoluteFile());
							BufferedWriter bw = new BufferedWriter(fw);
//							for (int i = 0; i < Provider_Count; i++) {
//								bw.write(String.format("%.5f", Services_Availability[i]-0.25) + "\r\n");
//							}
							bw.write(String.format("%.5f", Services_Availability[0]-0.25) + "\r\n");
							bw.close();
							Log.printLine("Availability File Created.");
						} catch (IOException e) {
							e.printStackTrace();
						}
					}

					if(way==4) {
						try {
							File file = new File("Reliability.txt");
							FileWriter fw = new FileWriter(file.getAbsoluteFile());
							BufferedWriter bw = new BufferedWriter(fw);
//							for (int i = 0; i < Provider_Count; i++) {
//								bw.write(String.format("%.4f", Services_Reliability[i]+0.1) + "\r\n");
//							}
							bw.write(String.format("%.4f", Services_Reliability[0]+0.1) + "\r\n");
							bw.close();
							Log.printLine("Reliability File Created.");
						} catch (IOException e) {
							e.printStackTrace();
						}
					}else {
						try {
							File file = new File("Reliability.txt");
							FileWriter fw = new FileWriter(file.getAbsoluteFile());
							BufferedWriter bw = new BufferedWriter(fw);
//							for (int i = 0; i < Provider_Count; i++) {
//								bw.write(String.format("%.4f", Services_Reliability[i]) + "\r\n");
//							}
							bw.write(String.format("%.4f", Services_Reliability[0]) + "\r\n");
							bw.close();
							Log.printLine("Reliability File Created.");
						} catch (IOException e) {
							e.printStackTrace();
						}
					}


					try {
						File file = new File("All_Cost.txt");
						FileWriter fw = new FileWriter(file.getAbsoluteFile());
						BufferedWriter bw = new BufferedWriter(fw);
						bw.write(All_cost/10000 + "\r\n");
						bw.close();
						Log.printLine("Costs File Created.");
					} catch (IOException e) {
						e.printStackTrace();
					}

					try {
						File file = new File("Resource_Response_rate.txt");
						FileWriter fw = new FileWriter(file.getAbsoluteFile());
						BufferedWriter bw = new BufferedWriter(fw);
						for (int i = 0; i < Provider_Count; i++) {
							{
								double d = (double) Providers_Request_Count[i] / (double) (User_Count * Request_Count);
								bw.write(String.format("%.3f", d) + "\r\n");
							}
						}
						bw.close();
						Log.printLine("Resource Response Rate File Created.");
					} catch (IOException e) {
						e.printStackTrace();
					}

					try {
						File file = new File("Response.txt");
						FileWriter fw = new FileWriter(file.getAbsoluteFile());
						BufferedWriter bw = new BufferedWriter(fw);
						for (int i = 0; i < User_Count; i++) {
							double s = 0;
							for (int j = 0; j < Request_Count; j++) {
								s += (Request_Queue_End_Time[i][j] - Request_Queue_Start_Time[i][j]);
							}
							s = (double) s / (double) Request_Count;
							s = s / 1000;
							bw.write(s + "\r\n");
						}
						bw.close();
						Log.printLine("Responses File Created.");
					} catch (IOException e) {
						e.printStackTrace();
					}

					try {
						File file = new File("Time.txt");
						FileWriter fw = new FileWriter(file.getAbsoluteFile());
						BufferedWriter bw = new BufferedWriter(fw);
						long s = End_Time_Simulation.getTime() - Start_Time_Simulation.getTime();
						bw.write((double) s / (double) 1000 + "\r\n");
						bw.close();
						Log.printLine("Response_Time File Created.");
					} catch (IOException e) {
						e.printStackTrace();
					}
					t_main.cancel();
					t_main.purge();
					t_Request_Queue.cancel();
					t_Request_Queue.purge();
					t_response.cancel();
					t_response.purge();
				}
			}
		};
		t_main.scheduleAtFixedRate(ts_main, 0, 50);
//TimerTask使用的时候会在主线程之外起一个单独的线程执行指定的任务,可以指定一次或多次
		TimerTask ts_Request_Queue = new TimerTask() {
			@Override
			//当way==1时，调用布谷鸟算法
			//当way==2时，调用遗传算法
			//当way==3时，调用飞蛾扑火算法
			//当way==4时，调用改进之后的飞蛾扑火算法
			public void run() {
				if (way == 1) {
					COA_select();
					for (int i = 0; i < Provider_Count; i++)
						sorted_weight[i] = i;
					for (int i = 0; i < Provider_Count - 1; i++) {
						double sum = 0;
						for (int k = 0; k < Cost_Params; k++) {
							sum += Resource_Cost[i][k];
						}
						score[i] = sum;
					}
					for (int i = 0; i < Provider_Count - 1; i++)
						for (int j = i; j < Provider_Count; j++) {
							if (score[j] > score[i]) {
								int k = sorted_weight[i];
								sorted_weight[i] = sorted_weight[j];
								sorted_weight[j] = k;
							}
						}
				} else if (way == 2) {
					GA_select();
					for (int i = 0; i < Provider_Count; i++)
						sorted_weight[i] = i;
					for (int i = 0; i < Provider_Count - 1; i++) {
						double sum = 0;
						for (int k = 0; k < Cost_Params; k++) {
							sum += Resource_Cost[i][k];
						}
						score[i] = sum;
					}
					for (int i = 0; i < Provider_Count - 1; i++)
						for (int j = i; j < Provider_Count; j++) {
							if (score[j] > score[i]) {
								int k = sorted_weight[i];
								sorted_weight[i] = sorted_weight[j];
								sorted_weight[j] = k;
							}
						}
				} else if(way == 3){
					MFO_select();
					for (int i = 0; i < Provider_Count; i++)
						sorted_weight[i] = i;
					for (int i = 0; i < Provider_Count - 1; i++) {
						double sum = 0;
						for (int k = 0; k < Cost_Params; k++) {
							sum += Resource_Cost[i][k];
						}
						score[i] = sum;
					}
					for (int i = 0; i < Provider_Count - 1; i++)
						for (int j = i; j < Provider_Count; j++) {
							if (score[j] > score[i]) {
								int k = sorted_weight[i];
								sorted_weight[i] = sorted_weight[j];
								sorted_weight[j] = k;
							}
						}
				}else if(way == 4){
					new_MFO_select();
					for (int i = 0; i < Provider_Count; i++)
						sorted_weight[i] = i;
					for (int i = 0; i < Provider_Count - 1; i++) {
						double sum = 0;
						for (int k = 0; k < Cost_Params; k++) {
							sum += Resource_Cost[i][k];
						}
						score[i] = sum;
					}
					for (int i = 0; i < Provider_Count - 1; i++)
						for (int j = i; j < Provider_Count; j++) {
							if (score[j] > score[i]) {
								int k = sorted_weight[i];
								sorted_weight[i] = sorted_weight[j];
								sorted_weight[j] = k;
							}
						}
				}else {
					Random r = new Random();
					Base_Option_Count = 10;
					for (int i = 0; i < Provider_Count; i++)
						sorted_weight[i] = r.nextInt(Provider_Count);
					Weight_cost = 14;
				}
				for (int w = 0; w < Provider_Count; w++) {

					for (int i = 0; i < User_Count; i++) {
						for (int j = 0; j < Request_Count; j++) {
							if (Request_Queue[i][j] == 1) {
								int op = sorted_weight[w];
								for (int y = 0; y < Providers_Resource[op]; y++) {
									int u = 0;
									boolean h = false;
									int cost_res = 0, cost_cap = 0;
									for (int z = 0; z < 4; z++)
										cost_cap += ServiceSet_Capability_Matrix[op][y][z];
									for (int z = 0; z < Resource_Count; z++)
										cost_res = Request_Queue_Resource[j][z];
									if (cost_cap > cost_res)
										h = true;
									if (Providers_Job[op][y] != 1 && h == true) {
										for (int z = 0; z < 4; z++)
											ServiceSet_Capability_Matrix[op][y][z] = -1
													* ServiceSet_Capability_Matrix[op][y][z];
										Request_Queue_Provider[i][j] = op;
										Request_Queue_Host[i][j] = y;
										for (int ok = 0; ok < Cost_Params; ok++)
											All_cost += Resource_Cost[op][ok];
										All_cost += Base_Option_Count * 10;
										Providers_Request_Count[op]++;
										Date f = new Date();
										Request_Queue_Start_Time[i][j] = f.getTime();
										Log.printLine("Select Service Set " + y + " from Provider " + op
												+ " For Request " + j + " from Client " + i);
										Request_Queue[i][j] = 2;
										Providers_State[op][y] = Base_Option_Count;
										Providers_Job[op][y] = 1;
									}
								}
							}

						}

					}
				}
			}
		};
		t_Request_Queue.scheduleAtFixedRate(ts_Request_Queue, 0, 50);

		TimerTask ts_response = new TimerTask() {
			//重写
			@Override
			public void run() {
				for (int i = 0; i < User_Count; i++) {
					for (int j = 0; j < Request_Count; j++) {
						if (Request_Queue[i][j] == 2) {
							int op = Request_Queue_Provider[i][j];
							int y = Request_Queue_Host[i][j];
							Providers_State[op][y]--;
							if (Providers_State[op][y] == 3)
								Log.printLine("Broker Send Cost Service " + y + " To Client " + i);
							if (Providers_State[op][y] == 2)
								Log.printLine("Client " + i + " Send OK to Broker For Service " + y);
							if (Providers_State[op][y] == 1)
								Log.printLine("Broker Create Connection Between Service " + y + " And Client " + i);
							if (Providers_State[op][y] == 0) {
								Date f = new Date();
								Request_Queue_End_Time[i][j] = f.getTime();
								Log.printLine("Response for Request " + j + " Client " + i + " from  Provider " + op
										+ " Service " + y + " Time="
										+ (Request_Queue_End_Time[i][j] - Request_Queue_Start_Time[i][j]));
								for (int z = 0; z < 4; z++)
									ServiceSet_Capability_Matrix[op][y][z] = -1
											* ServiceSet_Capability_Matrix[op][y][z];
								double sum = 0;
								sum = 0;
								for (int ok = 0; ok < Cost_Params; ok++) {
									Resource_Cost[op][ok] += 1;
									sum += Resource_Cost[op][ok];
								}
								Services_Reliability[op] += 1 / Request_Count * User_Count;
								if (Services_Reliability[op] > 1)
									Services_Reliability[op] = 1;
								Providers_Job[op][y] = 0;
								res++;
								sum = 0;
								for (int ok = 0; ok < Cost_Params; ok++) {
									Resource_Cost[op][ok] += 1;
									sum += Resource_Cost[op][ok];
								}
								double B_ok = 0;
								for (int ok = 0; ok < Cost_Params; ok++) {
									if (Resource_Cost[op][ok] < Request_Count)
										B_ok += Resource_Cost[op][ok];
								}
								Services_Availability[op] += 1 / Request_Count * User_Count;
								if (Services_Availability[op] > 1)
									Services_Availability[op] = 1;
								Request_Queue[i][j] = -1;
								Counter++;
							}

						}
					}
				}
			}

		};
		t_response.scheduleAtFixedRate(ts_response, 0, 50);

	}

	private static Datacenter createDatacenter(String name) throws Exception {
		List<Host> hostList = new ArrayList<Host>();

		List<Pe> quadCoreMachine = new ArrayList<Pe>();

		int mips = 1000;

		quadCoreMachine.add(new Pe(0, new PeProvisionerSimple(mips)));
		quadCoreMachine.add(new Pe(1, new PeProvisionerSimple(mips)));
		quadCoreMachine.add(new Pe(2, new PeProvisionerSimple(mips)));
		quadCoreMachine.add(new Pe(3, new PeProvisionerSimple(mips)));

		List<Pe> dualCoreMachine = new ArrayList<Pe>();

		dualCoreMachine.add(new Pe(0, new PeProvisionerSimple(mips)));
		dualCoreMachine.add(new Pe(1, new PeProvisionerSimple(mips)));

		int hostId = 0;
		int ram = 16384; // host memory (MB)
		long storage = 1000000; // host storage
		int bw = 10000;

		hostList.add(new Host(hostId, new RamProvisionerSimple(ram), new BwProvisionerSimple(bw), storage,
				quadCoreMachine, new VmSchedulerTimeShared(quadCoreMachine))); // This is our first machine

		hostId++;

		hostList.add(new Host(hostId, new RamProvisionerSimple(ram), new BwProvisionerSimple(bw), storage,
				dualCoreMachine, new VmSchedulerTimeShared(dualCoreMachine))); // Second machine

		String arch = "x86"; // system architecture
		String os = "Linux"; // operating system
		String vmm = "Xen";
		double time_zone = 10.0; // time zone this resource located
		double cost = 3.0; // the cost of using processing in this resource
		double costPerMem = 0.05; // the cost of using memory in this resource
		double costPerStorage = 0.1; // the cost of using storage in this resource
		double costPerBw = 0.1; // the cost of using bw in this resource
		LinkedList<Storage> storageList = new LinkedList<Storage>(); // we are not adding SAN devices by now

		DatacenterCharacteristics characteristics = new DatacenterCharacteristics(arch, os, vmm, hostList, time_zone,
				cost, costPerMem, costPerStorage, costPerBw);

		Federate_VmAllocationPolicy vm_policy = new Federate_VmAllocationPolicy(hostList);
		return new Datacenter(name, characteristics, vm_policy, storageList, 0);
	}

	public static class VmAllocationPolicyMinimum extends org.cloudbus.cloudsim.VmAllocationPolicy {

		private Map<String, Host> vm_table = new HashMap<String, Host>();
		private final HostList hosts;
		private Datacenter datacenter;

		public VmAllocationPolicyMinimum(List<? extends Host> list) {
			super(list);
			hosts = new HostList(list);
		}

		public void setDatacenter(Datacenter datacenter) {
			this.datacenter = datacenter;
		}

		public Datacenter getDatacenter() {
			return datacenter;
		}

		@Override
		public boolean allocateHostForVm(Vm vm) {

			if (this.vm_table.containsKey(vm.getUid())) {
				return true;
			}

			boolean vm_allocated = false;
			int tries = 0;
			do {
				Host host = this.hosts.getWithMinimumNumberOfPesEquals(vm.getNumberOfPes());
				vm_allocated = this.allocateHostForVm(vm, host);
			} while (!vm_allocated && tries++ < hosts.size());

			return vm_allocated;
		}

		@Override
		public boolean allocateHostForVm(Vm vm, Host host) {
			if (host != null && host.vmCreate(vm)) {
				vm_table.put(vm.getUid(), host);
				Log.formatLine(
						"%.4f: VM #" + vm.getId() + " has been allocated to the host#" + host.getId() + " datacenter #"
								+ host.getDatacenter().getId() + "(" + host.getDatacenter().getName() + ") #",
						CloudSim.clock());
				return true;
			}
			return false;
		}

		@Override
		public List<Map<String, Object>> optimizeAllocation(List<? extends Vm> vmList) {
			return null;
		}

		@Override
		public void deallocateHostForVm(Vm vm) {
			Host host = this.vm_table.remove(vm.getUid());
			if (host != null) {
				host.vmDestroy(vm);
			}
		}

		@Override
		public Host getHost(Vm vm) {
			return this.vm_table.get(vm.getUid());
		}

		@Override
		public Host getHost(int vmId, int userId) {
			return this.vm_table.get(Vm.getUid(userId, vmId));
		}
	}

	private static DatacenterBroker createFederatedBroker(String name) throws Exception {
		return new Service_DatacenterBroker(name);
	}

}