package sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.regex.Pattern;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_Bali extends Runnable_ClusterModel_MultiTransmission {

	// 0: TP, 1: NG, 2: CT, 3: HIV
	private static final int num_inf = 4;
	// 0: ANY, 1: URETHAL, 2: RECTAL
	private static final int num_site = 3;
	// 0 = ANY, 1 = ANAL
	private static final int num_act = 2;

	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("Bali_Model");

	// Note: Excl time column
	private static final int[] COL_SEL_INF_GENDER = new int[] { 2, 6, 10, 14 };
	private static final int[] COL_SEL_INF_SITE = new int[] { 0, 4, 5, 7, 8, 9 };
	private static final int[] COL_SEL_INF_GENDER_SITE = new int[] { 6, 19, 20, 31, 32, 42 };
	private static final int[] COL_SEL_INF_GENDER_SITE_AT = new int[] { 32, 33, 96, 98, 100, 102, 160, 162, 164, 166,
			224, 225 };

	final RandomGenerator rng_test;

	// Parameter - to be set
	public static final String DX_SWITCH_FILENAME = "dx_Bali.prop";
	protected final int TEST_RATE_ADJ = -1;

	protected int dxSwitchAt = -1;
	protected int dxSwitchGrp = 1 << 6;
	protected float dxSwitchProb = 1f;
	protected int[] dxSwitchTestTimeRange = new int[] { 360, 180 };
	protected float switchPrEPUptake = 0.0f;
	protected int switchPrEPAt = -1;

	protected final int INF_ID_HIV = 3;
	protected final int STAGE_PREP = 4;

	public Runnable_ClusterModel_Bali(long cMap_seed, long sim_seed, int[] pop_composition, ContactMap base_cMap,
			int numTimeStepsPerSnap, int numSnap) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap, num_inf, num_site,
				num_act);
		rng_test = new MersenneTwisterRandomGenerator(getSim_seed());
	}

	@Override
	public void initialse() {
		super.initialse();
		File dxSwitchFile = new File(baseDir, DX_SWITCH_FILENAME);
		if (dxSwitchFile.exists()) {
			try {
				Properties dxProp = new Properties();
				FileInputStream f = new FileInputStream(dxSwitchFile);
				dxProp.loadFromXML(f);
				f.close();
				dxSwitchAt = Integer.parseInt(dxProp.getProperty("DX_SWITCH_TIME", Integer.toString(dxSwitchAt)));
				dxSwitchGrp = Integer.parseInt(dxProp.getProperty("DX_SWITCH_GRP", Integer.toString(dxSwitchGrp)));
				dxSwitchProb = Float.parseFloat(dxProp.getProperty("DX_SWITCH_PROB", Float.toString(dxSwitchProb)));
				dxSwitchTestTimeRange = (int[]) PropValUtils.propStrToObject(
						dxProp.getProperty("DX_SWITCH_RANGE", Arrays.toString(dxSwitchTestTimeRange)), int[].class);
				switchPrEPAt = Integer.parseInt(dxProp.getProperty("PREP_SWITCH_AT", Integer.toString(switchPrEPAt)));
				switchPrEPUptake = Float
						.parseFloat(dxProp.getProperty("PREP_UPTAKE", Float.toString(switchPrEPUptake)));

			} catch (Exception e) {
				e.printStackTrace(System.err);
			}
		}
	}

	@Override
	protected void postTimeStep(int currentTime) {
		super.postTimeStep(currentTime);
		if (currentTime == dxSwitchAt) {
			for (Integer personId : test_rate_index_map.keySet()) {
				HashMap<Integer, Integer> test_rate_index = test_rate_index_map.get(personId);
				if ((test_rate_index.get(0).intValue() << dxSwitchGrp != 0) && rng_test.nextFloat() < dxSwitchProb) {
					test_rate_index.put(TEST_RATE_ADJ, 0);
					test_rate_index.remove(0); // Replace 1st schedule
					scheduleNextTest(personId, currentTime);
				}
			}
		}
	}

	@Override
	protected void testPerson(int currentTime, int pid, int infIncl, int siteIncl, int[][] cumul_treatment_by_person) {
		int[][] inf_stage = map_currrent_infection_stage.get(pid);
		// Negative HIV DX
		if (switchPrEPAt > 0 && currentTime >= switchPrEPAt && (infIncl & 1 << INF_ID_HIV) != 0
				&& (inf_stage == null || inf_stage[INF_ID_HIV][0] == AbstractIndividualInterface.INFECT_S)) {
			if (rng_test.nextFloat() < switchPrEPUptake) {
				int[][] current_stage_arr = map_currrent_infection_stage.get(pid);
				if (current_stage_arr == null) {
					current_stage_arr = new int[NUM_INF][NUM_SITE];
					for (int[] stage_by_infection : current_stage_arr) {
						Arrays.fill(stage_by_infection, AbstractIndividualInterface.INFECT_S);
					}
					map_currrent_infection_stage.put(pid, current_stage_arr);
				}
				int[][] infection_state_switch = map_infection_stage_switch.get(pid);
				if (infection_state_switch == null) {
					infection_state_switch = new int[NUM_INF][NUM_SITE];
					map_infection_stage_switch.put(pid, infection_state_switch);
				}

				current_stage_arr[INF_ID_HIV][0] = STAGE_PREP;
				int PrEP_switch_time = (int) Math
						.round(currentTime + dist_stage_period[INF_ID_HIV][0][STAGE_PREP].sample());

				updateInfectStageChangeSchedule(pid, INF_ID_HIV, 0, PrEP_switch_time,
						infection_state_switch[INF_ID_HIV][0]);
			}
		} else {
			super.testPerson(currentTime, pid, infIncl, siteIncl, cumul_treatment_by_person);
		}

	}

	@Override
	public void scheduleNextTest(Integer personId, int lastTestTime) {
		HashMap<Integer, Integer> test_rate_index = test_rate_index_map.get(personId);
		// Overwrite existing test rate
		if (test_rate_index != null && test_rate_index.containsKey(TEST_RATE_ADJ)) {
			int test_pt = test_rate_index.get(TEST_RATE_ADJ);

			int nextTestAfter = (int) (dxSwitchTestTimeRange[test_pt + 1]
					+ RNG.nextInt((int) dxSwitchTestTimeRange[test_pt] - (int) dxSwitchTestTimeRange[test_pt + 1]));

			int nextTestDate = lastTestTime + nextTestAfter;

			ArrayList<int[]> day_sch = schedule_testing.get(nextTestDate);

			if (day_sch == null) {
				day_sch = new ArrayList<>();
				schedule_testing.put(nextTestDate, day_sch);
			}

			int[] test_pair = new int[] { personId, 7, 7 };
			int pt_t = Collections.binarySearch(day_sch, test_pair, new Comparator<int[]>() {
				@Override
				public int compare(int[] o1, int[] o2) {
					int res = 0;
					int pt = 0;
					while (res == 0 && pt < o1.length) {
						res = Integer.compare(o1[pt], o2[pt]);
						pt++;
					}
					return res;
				}
			});

			if (pt_t < 0) {
				day_sch.add(~pt_t, test_pair);
			} else {
				int[] org_pair = day_sch.get(pt_t);
				org_pair[1] |= 7;
			}

		}
		super.scheduleNextTest(personId, lastTestTime);

	}

	@SuppressWarnings("unchecked")
	@Override
	protected void postSimulation() {
		super.postSimulation();
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
