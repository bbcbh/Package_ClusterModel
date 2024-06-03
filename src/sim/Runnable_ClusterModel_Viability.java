package sim;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_Viability extends Runnable_ClusterModel_MultiTransmission {

	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("ClusterModel_Viability");

	private static final int num_inf = 3; // TP, NG and CT
	private static final int num_site = 4;
	private static final int num_act = 5;

	// FILENAME_PREVALENCE_PERSON, FILENAME_CUMUL_INCIDENCE_PERSON
	private static final int[] COL_SEL_INF_GENDER = new int[] { 2, 6, 10 };
	// FILENAME_CUMUL_INCIDENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE = new int[] { 8, 25, 26, 27, 41, 42, 43 };
	// "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_SITE = new int[] { 0, 5, 6, 7, 9, 10, 11 };
	// "Infected_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE_AT = null; // new int[] { 64, 65, 192, 194, 196, 198, 200, 202,
																	// 204 160, 162, 164, 166, 224, 225 };

	protected static final int STAGE_ID_JUST_TREATED = -2;
	protected static final int STAGE_ID_NON_VIABLE = 3;

	protected float[][] prob_non_viabile_from_treatment; // new float[num_inf][num_site];
	protected float[][] dur_adj_non_viable_from_treatment; // new float[num_inf][num_site];

	protected RandomGenerator rng_viability;
	protected HashMap<String, Integer> pre_treatment_duration;
	protected int[] cumul_treatment_non_viable = new int[NUM_GENDER * NUM_INF];

	private static final String SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE = "SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE";

	public static final String PROP_PROB_NON_VIABLE_TREATMENT = "PROP_PROB_NON_VIABLE_TREATMENT";
	public static final String PROP_DUR_ADJ_NON_VIABLE_TREATMENT = "PROP_DUR_ADJ_NON_VIABLE_TREATMENT";

	public Runnable_ClusterModel_Viability(long cMap_seed, long sim_seed, ContactMap base_cMap, Properties prop) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);
		rng_viability = new MersenneTwisterRandomGenerator(sim_seed);

		String defaultStr = Arrays.deepToString(new float[num_inf][num_site]);
		prob_non_viabile_from_treatment = (float[][]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PROB_NON_VIABLE_TREATMENT, defaultStr), float[][].class);
		dur_adj_non_viable_from_treatment = (float[][]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_DUR_ADJ_NON_VIABLE_TREATMENT, defaultStr), float[][].class);

		for (int i = 0; i < num_inf; i++) {
			for (int s = 0; s < num_site; s++) {
				prob_non_viabile_from_treatment[i][s] = 0;
				dur_adj_non_viable_from_treatment[i][s] = 1;
			}
		}

		pre_treatment_duration = new HashMap<>();
		Arrays.fill(cumul_treatment_non_viable, 0);
	}

	@Override
	public ArrayList<Integer> loadOptParamter(String[] parameter_settings, double[] point, int[][] seedInfectNum,
			boolean display_only) {

		ArrayList<String> common_parameter_name = new ArrayList<>();
		ArrayList<Double> common_parameter_value = new ArrayList<>();

		for (int i = 0; i < parameter_settings.length; i++) {
			if (parameter_settings[i].startsWith(PROP_PROB_NON_VIABLE_TREATMENT)
					|| parameter_settings[i].startsWith(PROP_DUR_ADJ_NON_VIABLE_TREATMENT)) {

				String header = parameter_settings[i].startsWith(PROP_PROB_NON_VIABLE_TREATMENT)
						? PROP_PROB_NON_VIABLE_TREATMENT
						: PROP_DUR_ADJ_NON_VIABLE_TREATMENT;

				Matcher m = Pattern.compile(header + "_(\\d+)_(\\d+)").matcher(parameter_settings[i]);

				if (m.matches()) {
					int inf = Integer.parseInt(m.group(1));
					int site = Integer.parseInt(m.group(2));

					if (header.equals(PROP_PROB_NON_VIABLE_TREATMENT)) {
						prob_non_viabile_from_treatment[inf][site] = (float) point[i];
					} else {
						dur_adj_non_viable_from_treatment[inf][site] = (float) point[i];
					}

				}

			} else {
				common_parameter_name.add(parameter_settings[i]);
				common_parameter_value.add(point[i]);
			}

		}

		Double[] common_parameter_val_obj = common_parameter_value.toArray(new Double[common_parameter_value.size()]);

		double[] common_parameter_val = new double[common_parameter_value.size()];
		for (int i = 0; i < common_parameter_val.length; i++) {
			common_parameter_val[i] = common_parameter_val_obj[i].doubleValue();
		}

		return super.loadOptParamter(parameter_settings, point, seedInfectNum, display_only);
	}

	@Override
	protected void applyTreatment(int currentTime, int infId, int pid, int[][] inf_stage) {
		int[][] infection_switch = map_infection_stage_switch.get(pid);

		boolean hasNonViableOnly = true;
		for (int s = 0; s < num_site; s++) {
			if (inf_stage[infId][s] != AbstractIndividualInterface.INFECT_S
					&& prob_non_viabile_from_treatment[infId][s] > 0) {
				String key = String.format("%d,%d,%d", pid, infId, s);
				pre_treatment_duration.put(key, infection_switch[infId][s]);
			}
			hasNonViableOnly &= inf_stage[infId][s] == STAGE_ID_NON_VIABLE;
		}

		if (hasNonViableOnly) {
			int gI = getGenderType(pid);
			cumul_treatment_non_viable[gI * NUM_GENDER + infId]++;
		}

		super.applyTreatment(currentTime, infId, pid, inf_stage);
	}

	@Override
	protected int[] handleNoNextStage(Integer pid, int infection_id, int site_id, int current_infection_stage,
			int current_time) {

		int[] outcome = null; // {next_stage, duration}

		if (current_infection_stage == STAGE_ID_JUST_TREATED
				&& prob_non_viabile_from_treatment[infection_id][site_id] > 0) {
			if (rng_viability.nextFloat() < prob_non_viabile_from_treatment[infection_id][site_id]) {
				int next_stage = STAGE_ID_NON_VIABLE;
				String key = String.format("%d,%d,%d", pid, infection_id, site_id);
				Integer current_infection_state_switch_at = pre_treatment_duration.remove(key);
				if (current_infection_state_switch_at != null) {
					int non_viable_until = current_time + Math.round((current_infection_state_switch_at - current_time)
							* dur_adj_non_viable_from_treatment[infection_id][site_id]);
					outcome = new int[] { next_stage, non_viable_until };
				}
			}
		}

		if (outcome == null) {
			return super.handleNoNextStage(pid, infection_id, site_id, current_infection_stage, current_time);
		} else {
			return outcome;
		}
	}

	@Override
	protected void postTimeStep(int currentTime) {
		super.postTimeStep(currentTime);
		if (currentTime % nUM_TIME_STEPS_PER_SNAP == 0) {
			@SuppressWarnings("unchecked")
			HashMap<Integer, int[]> countMap = (HashMap<Integer, int[]>) sim_output
					.get(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE);
			if (countMap == null) {
				countMap = new HashMap<>();
				sim_output.put(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE, countMap);
			}
			countMap.put(currentTime, Arrays.copyOf(cumul_treatment_non_viable, cumul_treatment_non_viable.length));
		}
	}

	@SuppressWarnings("unchecked")
	protected void postSimulation() {
		String key, fileName;
		HashMap<Integer, int[]> countMap;
		String filePrefix = this.getRunnableId() == null ? "" : this.getRunnableId();

		// Non-viable treatment
		countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE);
		if (countMap != null) {
			fileName = String.format(filePrefix + "Treatment_Non_Viable_%d_%d.csv", cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);
		}

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
