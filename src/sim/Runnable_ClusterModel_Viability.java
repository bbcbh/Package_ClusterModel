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

	protected static final int[][] STAGE_ID_NON_VIABLE = new int[][] { null, new int[] { 3, 3, 3, 3 },
			new int[] { 3, 3, 3, 3 } };

	protected float[][] prob_non_viabile_from_treatment; // new float[num_inf][num_site];
	protected float[][] dur_adj_non_viable_from_treatment; // new float[num_inf][num_site];

	protected float[][][] prob_non_viabile_from_transmission; // new float[num_inf][num_site_from][num_site_to];
	protected float[][][] dur_adj_non_viable_from_transmission; // new float[num_inf][num_site]{parameter};

	protected RandomGenerator rng_viability;
	protected int[] cumul_treatment_non_viable = new int[NUM_GENDER * NUM_INF];
	protected int[] cumul_treatment_non_viable_site = new int[NUM_GENDER * NUM_INF * NUM_SITE];
	protected int[] cumul_treatment_non_infected_site = new int[NUM_GENDER * NUM_INF * NUM_SITE];

	private final String[] optParameter_optionsHeader;
	private final Object[] optParameter_optionsTar;

	private static final String SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE = "SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE";
	private static final String SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE_SITE = "SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE_SITE";
	private static final String SIM_OUTPUT_KEY_TREATMENT_NON_INFECTED_SITE = "SIM_OUTPUT_KEY_TREATMENT_NON_INFECTED_SITE";

	public static final String PROP_PROB_NON_VIABLE_TREATMENT = "PROP_PROB_NON_VIABLE_TREATMENT";
	public static final String PROP_DUR_ADJ_NON_VIABLE_TREATMENT = "PROP_DUR_ADJ_NON_VIABLE_TREATMENT";

	public static final String PROP_PROB_NON_VIABLE_TRANSMISSION = "PROP_PROB_NON_VIABLE_TRANSMISSION";
	public static final String PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION = "PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION";

	public Runnable_ClusterModel_Viability(long cMap_seed, long sim_seed, ContactMap base_cMap, Properties prop) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);
		rng_viability = new MersenneTwisterRandomGenerator(sim_seed);

		String defaultStr = Arrays.deepToString(new float[num_inf][num_site]);

		prob_non_viabile_from_treatment = (float[][]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PROB_NON_VIABLE_TREATMENT, defaultStr), float[][].class);
		dur_adj_non_viable_from_treatment = (float[][]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_DUR_ADJ_NON_VIABLE_TREATMENT, defaultStr), float[][].class);

		defaultStr = Arrays.deepToString(new float[num_inf][num_site][num_site]);
		prob_non_viabile_from_transmission = (float[][][]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PROB_NON_VIABLE_TRANSMISSION, defaultStr), float[][][].class);
		defaultStr = Arrays.deepToString(new float[num_inf][num_site][2]);
		dur_adj_non_viable_from_transmission = (float[][][]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION, defaultStr), float[][][].class);

		for (int i = 0; i < num_inf; i++) {
			for (int s = 0; s < num_site; s++) {
				if (prop.getProperty(PROP_PROB_NON_VIABLE_TREATMENT) == null) {
					prob_non_viabile_from_treatment[i][s] = 0;
				}
				if (prop.getProperty(PROP_DUR_ADJ_NON_VIABLE_TREATMENT) == null) {
					dur_adj_non_viable_from_treatment[i][s] = 1;
				}
				if (prop.getProperty(PROP_PROB_NON_VIABLE_TRANSMISSION) == null) {
					Arrays.fill(prob_non_viabile_from_transmission[i][s], 0);
				}
				if (prop.getProperty(PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION) == null) {
					Arrays.fill(dur_adj_non_viable_from_transmission[i][s], 0);
				}
			}
		}

		Arrays.fill(cumul_treatment_non_viable, 0);
		Arrays.fill(cumul_treatment_non_viable_site, 0);

		optParameter_optionsHeader = new String[] { PROP_PROB_NON_VIABLE_TREATMENT, PROP_DUR_ADJ_NON_VIABLE_TREATMENT,
				PROP_PROB_NON_VIABLE_TRANSMISSION, PROP_DUR_ADJ_NON_VIABLE_TRANSMISSION };

		optParameter_optionsTar = new Object[] { prob_non_viabile_from_treatment, dur_adj_non_viable_from_treatment,
				prob_non_viabile_from_transmission, dur_adj_non_viable_from_transmission };
	}

	@Override
	public ArrayList<Integer> loadOptParameter(String[] parameter_settings, double[] point, int[][] seedInfectNum,
			boolean display_only) {

		ArrayList<String> common_parameter_name = new ArrayList<>();
		ArrayList<Double> common_parameter_value = new ArrayList<>();

		for (int i = 0; i < parameter_settings.length; i++) {
			String header = null;
			Object target = null;
			for (int j = 0; j < optParameter_optionsHeader.length; j++) {
				if (parameter_settings[i].startsWith(optParameter_optionsHeader[j])) {
					header = optParameter_optionsHeader[j];
					target = optParameter_optionsTar[j];
				}
			}
			if (header != null && target != null) {
				Matcher m;
				if (target instanceof float[][]) {
					m = Pattern.compile(header + "_(\\d+)_(\\d+)").matcher(parameter_settings[i]);
					if (m.matches()) {
						int inf = Integer.parseInt(m.group(1));
						int site = Integer.parseInt(m.group(2));
						((float[][]) target)[inf][site] = (float) point[i];
					}
				} else if (target instanceof float[][][]) {
					m = Pattern.compile(header + "_(\\d+)_(\\d+)_(\\d+)").matcher(parameter_settings[i]);
					if (m.matches()) {
						((float[][][]) target)[Integer.parseInt(m.group(1))][Integer.parseInt(m.group(2))][Integer
								.parseInt(m.group(3))] = (float) point[i];
					}

				} else {
					System.err.printf("Warning: Opt parameter for %s type mismatch. Value ignored.\n", header);
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

		return super.loadOptParameter(common_parameter_name.toArray(new String[common_parameter_name.size()]),
				common_parameter_val, seedInfectNum, display_only);
	}

	@Override
	protected void simulate_transmission_failed_act(int currentTime, int inf_id, Integer pid_inf_src, int pid_inf_tar,
			int src_site, int tar_site) {
		// Non-viable from transmission
		super.simulate_transmission_failed_act(currentTime, inf_id, pid_inf_src, pid_inf_tar, src_site, tar_site);
		if (prob_non_viabile_from_transmission[inf_id][src_site][tar_site] > 0) {
			if (rng_viability.nextFloat() < prob_non_viabile_from_transmission[inf_id][src_site][tar_site]) {
				int non_via_duration_range = Math.round(dur_adj_non_viable_from_transmission[inf_id][tar_site][1]
						- dur_adj_non_viable_from_transmission[inf_id][tar_site][0]);
				int non_viable_until = currentTime
						+ Math.round(dur_adj_non_viable_from_transmission[inf_id][tar_site][0])
						+ (non_via_duration_range > 0 ? rng_viability.nextInt(non_via_duration_range) : 0);

				if (non_viable_until > currentTime) {
					int[][] inf_stage = map_currrent_infection_stage.get(pid_inf_tar);
					if (inf_stage == null) {
						inf_stage = new int[NUM_INF][NUM_SITE];
						for (int[] stage_by_infection : inf_stage) {
							Arrays.fill(stage_by_infection, AbstractIndividualInterface.INFECT_S);
						}
						map_currrent_infection_stage.put(pid_inf_tar, inf_stage);
					}
					int[][] infection_state_switch = map_infection_stage_switch.get(pid_inf_tar);
					if (infection_state_switch == null) {
						infection_state_switch = new int[NUM_INF][NUM_SITE];
						map_infection_stage_switch.put(pid_inf_tar, infection_state_switch);
					}

					inf_stage[inf_id][tar_site] = STAGE_ID_NON_VIABLE[inf_id][tar_site];
					updateInfectionStage(pid_inf_tar, inf_id, tar_site, STAGE_ID_NON_VIABLE[inf_id][tar_site],
							currentTime, inf_stage, infection_state_switch, non_viable_until - currentTime);
				}

			}

		}
	}

	@Override
	protected boolean isValidInfectionTargetSite(int inf_id, int tar_site, int[][] tar_infection_stages) {
		return super.isValidInfectionTargetSite(inf_id, tar_site, tar_infection_stages)
				|| ((STAGE_ID_NON_VIABLE[inf_id] != null)
						&& (tar_infection_stages[inf_id][tar_site] == STAGE_ID_NON_VIABLE[inf_id][tar_site]));
	}

	@Override
	protected void applyTreatment(int currentTime, int infId, int pid, int[][] inf_stage) {
		int[][] infection_switch = map_infection_stage_switch.get(pid);
		boolean hasInfectiousPreTreatment = false;
		int[] becomeNonViableUtil = new int[num_site];
		for (int siteId = 0; siteId < num_site; siteId++) {
			boolean is_infectious_stage = (inf_stage[infId][siteId] >= 0)
					&& (lookupTable_infection_infectious_stages[infId][siteId] & (1 << inf_stage[infId][siteId])) != 0;
			hasInfectiousPreTreatment |= is_infectious_stage;
			if (is_infectious_stage && prob_non_viabile_from_treatment[infId][siteId] > 0) {
				// Check if possible to turn viable to non-viable through treatment
				int preTreatmentDur = infection_switch[infId][siteId];
				if (rng_viability.nextFloat() < prob_non_viabile_from_treatment[infId][siteId]) {
					int non_viable_until = currentTime + Math
							.round((preTreatmentDur - currentTime) * dur_adj_non_viable_from_treatment[infId][siteId]);
					becomeNonViableUtil[siteId] = non_viable_until;
				}
			}
			if (!is_infectious_stage) {
				if (inf_stage[infId][siteId] >= 0) {
					cumul_treatment_non_viable_site[infId * NUM_GENDER * NUM_SITE + getGenderType(pid) * NUM_SITE
							+ siteId]++;
				} else {
					cumul_treatment_non_infected_site[infId * NUM_GENDER * NUM_SITE + getGenderType(pid) * NUM_SITE
							+ siteId]++;
				}
			}
		}

		if (!hasInfectiousPreTreatment) {
			cumul_treatment_non_viable[infId * NUM_GENDER + getGenderType(pid)]++;
		}else {
			super.applyTreatment(currentTime, infId, pid, inf_stage);
		}

		// Set non-viable infection post treatment
		for (int siteId = 0; siteId < num_site; siteId++) {
			if (becomeNonViableUtil[siteId] > 0) {
				inf_stage[infId][siteId] = STAGE_ID_NON_VIABLE[infId][siteId];
				updateInfectStageChangeSchedule(pid, infId, siteId, becomeNonViableUtil[siteId], currentTime + 1);
				updateInfectionStage(pid, infId, siteId, STAGE_ID_NON_VIABLE[infId][siteId], currentTime, inf_stage,
						infection_switch, becomeNonViableUtil[siteId] - currentTime);
			}
		}

	}

	@SuppressWarnings("unchecked")
	@Override
	protected void postTimeStep(int currentTime) {
		super.postTimeStep(currentTime);
		if (currentTime % nUM_TIME_STEPS_PER_SNAP == 0) {
			HashMap<Integer, int[]> countMap;
			countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE);
			if (countMap == null) {
				countMap = new HashMap<>();
				sim_output.put(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE, countMap);
			}
			countMap.put(currentTime, Arrays.copyOf(cumul_treatment_non_viable, cumul_treatment_non_viable.length));
			countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE_SITE);
			if (countMap == null) {
				countMap = new HashMap<>();
				sim_output.put(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE_SITE, countMap);
			}
			countMap.put(currentTime,
					Arrays.copyOf(cumul_treatment_non_viable_site, cumul_treatment_non_viable_site.length));

			countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_TREATMENT_NON_INFECTED_SITE);
			if (countMap == null) {
				countMap = new HashMap<>();
				sim_output.put(SIM_OUTPUT_KEY_TREATMENT_NON_INFECTED_SITE, countMap);
			}
			countMap.put(currentTime,
					Arrays.copyOf(cumul_treatment_non_infected_site, cumul_treatment_non_infected_site.length));

		}
	}

	@Override
	@SuppressWarnings("unchecked")
	protected void postSimulation() {
		super.postSimulation();
		
		String key, fileName;
		HashMap<Integer, int[]> countMap;
		String filePrefix = this.getRunnableId() == null ? "" : this.getRunnableId();

		// Non-viable treatment
		countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE);
		if (countMap != null) {
			fileName = String.format(filePrefix + "Treatment_Non_Viable_Person_%d_%d.csv", cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);
		}

		countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_TREATMENT_NON_VIABLE_SITE);
		if (countMap != null) {
			fileName = String.format(filePrefix + "Treatment_Non_Viable_Site_%d_%d.csv", cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d_Site_%d", new int[] { NUM_INF, NUM_GENDER, NUM_SITE },
					COL_SEL_INF_GENDER_SITE);
		}

		countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_TREATMENT_NON_INFECTED_SITE);
		if (countMap != null) {
			fileName = String.format(filePrefix + "Treatment_Non_Infected_Site_%d_%d.csv", cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d_Site_%d", new int[] { NUM_INF, NUM_GENDER, NUM_SITE },
					COL_SEL_INF_GENDER_SITE);
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
					for (String[] col_name : col_names) {
						pWri.print(',');
						String ckey = String.format("%s,%s,%s", col_name[0], col_name[1], col_name[2]);
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
