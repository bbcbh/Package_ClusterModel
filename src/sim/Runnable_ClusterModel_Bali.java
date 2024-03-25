package sim;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.regex.Pattern;

import relationship.ContactMap;

public class Runnable_ClusterModel_Bali extends Runnable_ClusterModel_MultiTransmission {
	
	// 0: TP, 1: NG, 2: CT, 3: HIV
	private static final int num_inf = 4; 
	// 0: ANY, 1: URETHAL, 2: RECTAL
	private static final int num_site = 3;
	 //0 = ANY, 1 = ANAL
	private static final int num_act = 2;
	
	
	public static final Pattern PROP_TYPE_PATTERN = Pattern
			.compile("Bali_Model");
	
	// Note: Excl time column
	private static final int[] COL_SEL_INF_GENDER = new int[] { 2, 6 };
	private static final int[] COL_SEL_INF_SITE = new int[] { 0, 5, 6, 7 };
	private static final int[] COL_SEL_INF_GENDER_SITE = new int[] { 8, 25, 26, 27 };
	private static final int[] COL_SEL_INF_GENDER_SITE_AT = new int[] { 64, 65, 192, 194, 196, 198, 200, 202, 204,
			206 };
	

	public Runnable_ClusterModel_Bali(long cMap_seed, long sim_seed, int[] pop_composition, ContactMap base_cMap,
			int numTimeStepsPerSnap, int numSnap) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap, num_inf, num_site, num_act);		
	}
	
	@SuppressWarnings("unchecked")
	protected void postSimulation() {
		String key, fileName;
		HashMap<Integer, int[]> countMap;
		String filePrefix = this.getRunnableId() == null ? "" : this.getRunnableId();

		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE) != 0) {
			key = String.format(SIM_OUTPUT_KEY_CUMUL_TREATMENT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE);

			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);

		}
		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE) != 0) {
			key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);

			key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE_SITE,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);

			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_SITE,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d_Site_%d", new int[] { NUM_INF, NUM_GENDER, NUM_SITE },
					COL_SEL_INF_GENDER_SITE);

		}

		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE) != 0) {

			key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_GENDER_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(
					filePrefix + "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);

			key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_SITE_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(
					filePrefix + "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Site_%d", new int[] { NUM_INF, NUM_SITE }, COL_SEL_INF_SITE);

			key = String.format(SIM_OUTPUT_KEY_INFECTED_AT_GENDER_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(
					filePrefix + "Infected_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE, cMAP_SEED,
					sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d_Infected_SiteInc_%d",
					new int[] { NUM_INF, NUM_GENDER, 1 << (NUM_SITE + 1) }, COL_SEL_INF_GENDER_SITE_AT);

			key = String.format(SIM_OUTPUT_KEY_INFECTED_SITE_STAGE_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			HashMap<Integer, HashMap<String, Integer>> infected_site_stage_count = (HashMap<Integer, HashMap<String, Integer>>) sim_output
					.get(key);

			fileName = String.format(
					filePrefix + "Infected_All_Stages_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE,
					cMAP_SEED, sIM_SEED);
			try {
				PrintWriter pWri = new PrintWriter(new FileWriter(new java.io.File(baseDir, fileName)));
				Integer[] timeArr = infected_site_stage_count.keySet().toArray(new Integer[0]);
				Arrays.sort(timeArr);
				ArrayList<String[]> col_names = new ArrayList<>();

				for (Integer time : timeArr) {
					HashMap<String, Integer> ent = infected_site_stage_count.get(time);
					for (String col_name : ent.keySet()) {
						String[] cs = col_name.split(",");
						int pt = Collections.binarySearch(col_names, cs, new Comparator<String[]>() {
							@Override
							public int compare(String[] s1, String[] s2) {
								int res = 0;
								for (int i = 0; i < s1.length; i++) {
									res = Integer.compare(Integer.parseInt(s1[i]), Integer.parseInt(s2[i]));
									if (res != 0) {
										return res;
									}
								}
								return res;
							}
						});
						if (pt < 0) {
							col_names.add(~pt, cs);
						}
					}

				}

				// Header
				pWri.print("Time");
				for (String[] col_name_split : col_names) {
					pWri.print(',');
					pWri.printf("Inf_%s_Site_%s_Stage_%s", col_name_split[0], col_name_split[1], col_name_split[2]);
				}
				pWri.println();

				for (Integer time : timeArr) {
					HashMap<String, Integer> ent = infected_site_stage_count.get(time);
					pWri.print(time);
					for (int i = 0; i < col_names.size(); i++) {
						pWri.print(',');
						String ckey = String.format("%s,%s,%s", col_names.get(i)[0], col_names.get(i)[1],
								col_names.get(i)[2]);
						Integer val = ent.get(ckey);
						if (val == null) {
							val = 0;
						}
						pWri.print(val);
					}
					pWri.println();
				}
				pWri.close();

			} catch (IOException e) {
				e.printStackTrace(System.err);
			}

		}

		if (print_progress != null && runnableId != null) {
			try {
				print_progress.printf("Post simulation file generation for Thread <%s> completed. Timestamp = %tc.\n",
						runnableId, System.currentTimeMillis());
			} catch (Exception ex) {
				System.err.printf("Post simulation file generation for Thread <%s> completed.\n", runnableId);
			}
		}
	}

}
