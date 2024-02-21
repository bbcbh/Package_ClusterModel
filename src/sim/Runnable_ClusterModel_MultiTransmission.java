package sim;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.RealDistribution;

import person.AbstractIndividualInterface;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_MultiTransmission extends Abstract_Runnable_ClusterModel_Transmission {

	public Object[] runnable_fields = {
			// 0: RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ
			// double{{ACT_TYPE, GENDER_INCLUDE_INDEX_FROM, GENDER_INCLUDE_INDEX_TO,
			// SITE_P1, SITE_P2, SITE_1_RISK_GRP_INCLUDE_INDEX, ACT_PER_DAY},...}
			// or
			// double{{ACT_TYPE, GENDER_INCLUDE_INDEX_FROM, GENDER_INCLUDE_INDEX_TO,
			// SITE_P1, SITE_P2, SITE_1_RISK_GRP_INCLUDE_INDEX, ACT_PER_DAY,
			// CONDOM_EFFICACY, USAGE_REG, USAGE_CASUAL},...}
			new double[][] {},
			// 1: RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE
			// double{{INF_ID, STATE_INCLUDE_INDEX, SITE_FROM, SITE_TO,
			// TRANSMISSION_PARAM_0,
			// TRANSMISSION_PARAM_1...},... }
			new double[][] {},
			// 2: RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD
			// double{{INF_ID, SITE, STAGE_INCLUDE_INDEX}, ...}
			new double[][] {},
			// 3: RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD
			// double{{INF_ID, SITE, STAGE_ID_FROM, STAGE_ID_TO, SWTICH_STAGE_PROB,
			// DURATION_STAGE_ID_TO_PARAM_0 ...},...}
			new double[][] {},
			// 4: RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			// double{{INF_ID, SITE, STAGE_ID, PROB_PARAM},...}
			new double[][] {},
			// 5: RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			// double{{GENDER_INCLUDE_INDEX, RISK_GRP_DEF_ID, NUM_CASUAL_PARTNER_LOWER,
			// NUM_CASUAL_PARTNER_UPPER}, ... }
			new double[][] {},
			// 6: RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES
			// double{{GENDER_INCLUDE_INDEX, INF_INCLUDE_INDEX, SITE_INCLUDE_INDEX,
			// RISK_GRP_INCLUDE_INDEX, TESTING_RATE_PARAM...}, ...}
			new double[][] {},
			// 7: RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM;
			// double{{GENDER_INCLUDE_INDEX, SITE, DAILY_SOUGHT_TEST_RATE}, ...}
			new double[][] {},
			// 8 RUNNABLE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES
			// double{{INF_INDEX, SITE_INDEX, TEST_ACCURACY,
			// TARGET_STAGE_INC, TREATMENT_SUC_STAGE}, ...}
			new double[][] {}, };

	// Fixed input value
	protected final int NUM_INF;
	protected final int NUM_SITE;
	protected final int NUM_ACT;
	protected final int NUM_GENDER;

	// For RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ
	public static final int FIELD_ACT_FREQ_ACT_TYPE = 0;
	public static final int FIELD_ACT_FREQ_GENDER_INCLUDE_INDEX_FROM = FIELD_ACT_FREQ_ACT_TYPE + 1;
	public static final int FIELD_ACT_FREQ_GENDER_INCLUDE_INDEX_TO = FIELD_ACT_FREQ_GENDER_INCLUDE_INDEX_FROM + 1;
	public static final int FIELD_ACT_FREQ_SITE_P1 = FIELD_ACT_FREQ_GENDER_INCLUDE_INDEX_TO + 1;
	public static final int FIELD_ACT_FREQ_SITE_P2 = FIELD_ACT_FREQ_SITE_P1 + 1;
	public static final int FIELD_ACT_FREQ_P1_RISK_GRP_INCLUDE_INDEX = FIELD_ACT_FREQ_SITE_P2 + 1;
	public static final int FIELD_ACT_FREQ_ACT_PER_DAY = FIELD_ACT_FREQ_P1_RISK_GRP_INCLUDE_INDEX + 1;
	public static final int FIELD_ACT_FREQ_CONDOM_EFFICACY = FIELD_ACT_FREQ_ACT_PER_DAY + 1;
	public static final int FIELD_ACT_FREQ_USAGE_REG = FIELD_ACT_FREQ_CONDOM_EFFICACY + 1;
	public static final int FIELD_ACT_FREQ_USAGE_CASUAL = FIELD_ACT_FREQ_USAGE_REG + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE
	public static final int FIELD_TRANSMISSION_RATE_INF_ID = 0;
	public static final int FIELD_TRANSMISSION_RATE_STATE_INCLUDE_INDEX = FIELD_TRANSMISSION_RATE_INF_ID + 1;
	public static final int FIELD_TRANSMISSION_RATE_SITE_FROM = FIELD_TRANSMISSION_RATE_STATE_INCLUDE_INDEX + 1;
	public static final int FIELD_TRANSMISSION_RATE_SITE_TO = FIELD_TRANSMISSION_RATE_SITE_FROM + 1;
	public static final int FIELD_TRANSMISSION_RATE_TRANS_PARAM_START = FIELD_TRANSMISSION_RATE_SITE_TO + 1;
	// RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD
	public static final int FIELD_INFECTIOUS_PERIOD_INF_ID = 0;
	public static final int FIELD_INFECTIOUS_PERIOD_SITE = FIELD_INFECTIOUS_PERIOD_INF_ID + 1;
	public static final int FIELD_INFECTIOUS_PERIOD_STAGE_INCLUDE_INDEX = FIELD_INFECTIOUS_PERIOD_SITE + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD
	public static final int FIELD_STAGE_PERIOD_INF_ID = 0;
	public static final int FIELD_STAGE_PERIOD_SITE = FIELD_STAGE_PERIOD_INF_ID + 1;
	public static final int FIELD_STAGE_PERIOD_STATE_FROM = FIELD_STAGE_PERIOD_SITE + 1;
	public static final int FIELD_STAGE_PERIOD_STATE_TO = FIELD_STAGE_PERIOD_STATE_FROM + 1;
	public static final int FIELD_STAGE_PERIOD_SWTICH_STAGE_PROB = FIELD_STAGE_PERIOD_STATE_TO + 1;
	public static final int FIELD_STAGE_PERIOD_DURATION_PARAM_START = FIELD_STAGE_PERIOD_SWTICH_STAGE_PROB + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
	public static final int FIELD_SYM_RATE_INF_ID = 0;
	public static final int FIELD_SYM_RATE_SITE = FIELD_SYM_RATE_INF_ID + 1;
	public static final int FIELD_SYM_RATE_STATE_ID = FIELD_SYM_RATE_SITE + 1;
	public static final int FIELD_SYM_RATE_PROB = FIELD_SYM_RATE_STATE_ID + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS
	public static final int FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_GENDER_INCLUDE_INDEX = 0;
	public static final int FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_RISK_GRP_DEF_ID = FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_GENDER_INCLUDE_INDEX
			+ 1;
	public static final int FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_NUM_CASUAL_PARTNER_LOWER = FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_RISK_GRP_DEF_ID
			+ 1;
	public static final int FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_NUM_CASUAL_PARTNER_UPPER = FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_NUM_CASUAL_PARTNER_LOWER
			+ 1;
	// For RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES
	public static final int FIELD_TESTING_RATE_BY_RISK_CATEGORIES_GENDER_INCLUDE_INDEX = 0;
	public static final int FIELD_TESTING_RATE_BY_RISK_CATEGORIES_INF_INCLUDE_INDEX = FIELD_TESTING_RATE_BY_RISK_CATEGORIES_GENDER_INCLUDE_INDEX
			+ 1;
	public static final int FIELD_TESTING_RATE_BY_RISK_CATEGORIES_SITE_INCLUDE_INDEX = FIELD_TESTING_RATE_BY_RISK_CATEGORIES_INF_INCLUDE_INDEX
			+ 1;
	public static final int FIELD_TESTING_RATE_BY_RISK_CATEGORIES_RISK_GRP_INCLUDE_INDEX = FIELD_TESTING_RATE_BY_RISK_CATEGORIES_SITE_INCLUDE_INDEX
			+ 1;
	public static final int FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START = FIELD_TESTING_RATE_BY_RISK_CATEGORIES_RISK_GRP_INCLUDE_INDEX
			+ 1;
	// For RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
	public static final int FIELD_SOUGHT_TEST_PERIOD_BY_SYM_INF_ID = 0;
	public static final int FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SITE = FIELD_SOUGHT_TEST_PERIOD_BY_SYM_INF_ID + 1;
	public static final int FIELD_SOUGHT_TEST_PERIOD_BY_SYM_DAILY_SOUGHT_TEST_RATE = FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SITE
			+ 1;
	// For RUNNALBE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES
	public static final int FIELD_DX_TEST_PROPERTIES_INF_ID = 0;
	public static final int FIELD_DX_TEST_PROPERTIES_SITE = FIELD_DX_TEST_PROPERTIES_INF_ID + 1;
	public static final int FIELD_DX_TEST_PROPERTIES_ACCURACY = FIELD_DX_TEST_PROPERTIES_SITE + 1;
	public static final int FIELD_DX_TEST_PROPERTIES_STAGE_INC = FIELD_DX_TEST_PROPERTIES_ACCURACY + 1;
	public static final int FIELD_DX_TEST_PROPERTIES_TREATMENT_SUC_STAGE = FIELD_DX_TEST_PROPERTIES_STAGE_INC + 1;

	// Transient field
	protected transient HashMap<Integer, double[][][][]> map_trans_prob; // Key=PID,V=double[INF_ID][SITE_FROM][SITE_TO][STAGE_ID]
	protected transient HashMap<Integer, int[][]> map_currrent_infection_stage; // Key=PID,V=int[INF_ID][SITE]{infection_stage}
	protected transient HashMap<Integer, int[][]> map_infection_stage_switch; // Key=PID,V=int[INF_ID][SITE]{switch_time_at_site_0};
	protected transient HashMap<String, ArrayList<Integer>> map_currently_infectious; // Key=Inf_ID,SiteID,V=ArrayList
																						// of
																						// pid with infectious

	// For schedule, key are day and entries are list of person id, ordered
	protected transient HashMap<Integer, ArrayList<int[]>> schedule_testing; // V=int[]{PID,INF_INCLUDE_INDEX,SITE_INCLUDE_INDEX}
	protected transient HashMap<Integer, ArrayList<ArrayList<ArrayList<Integer>>>> schedule_stage_change;
	// V=ArrayList<Inf><SITE>{PIds of those need to change}

	// Helper objects set by fields

	protected transient HashMap<String, double[]> lookupTable_infection_stage_path; // Key="INF_ID,SITE_ID,STAGE_ID",V={Cumul_Prob_0,...
																					// State_ID_0...}
	protected transient int[][] lookupTable_infection_infectious_stages; // int[INF_ID][SITE]=STAGE_INCLUDE_ID
	protected transient HashMap<String, double[]> lookupTable_testing_properties; // Key="INF_ID,SITE_ID,
																					// V=TEST_ACCURACY_INFO

	protected transient double[][][][] table_act_frequency; // double[ACT_ID][G_TO][G_FROM]=fieldEntry

	protected transient RealDistribution[][][] dist_sym_rate; // RealDistribution[INF_ID][SITE][STAGE_ID]
	protected transient RealDistribution[][][][] dist_tranmissionMatrix; // RealDistribution[INF_ID][SITE_FROM][SITE_TO][STAGE_ID]
	protected transient RealDistribution[][][] dist_stage_period; // ReadDistribution[INF_ID][SITE][STAGE_ID]

	protected static final int STAGE_ID_JUST_INFECTED = Integer.MIN_VALUE;
	
	protected static final String SIM_OUTPUT_KEY_INFECTIOUS_GENDER_COUNT = "Output_%d";
	protected static final String SIM_OUTPUT_KEY_INFECTIOUS_SITE_COUNT = "Output_%d_Infectious_Site";
	protected static final String SIM_OUTPUT_KEY_INFECTED_SITE_STAGE_COUNT =  "Output_%d_Infected_Site_Stage";
	protected static final String SIM_OUTPUT_KEY_CUMUL_TREATMENT = "Output_%d";
	protected static final String SIM_OUTPUT_KEY_CUMUL_INCIDENCE = "Output_%d";
	protected static final String SIM_OUTPUT_KEY_CUMUL_INCIDENCE_SITE = "Output_%d_S";
	
	
	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("MultiTransmission_(\\d+)_(\\d+)_(\\d+)");
	
	

	public Runnable_ClusterModel_MultiTransmission(long cMap_seed, long sim_seed, int[] pop_composition,
			ContactMap base_cMap, int numTimeStepsPerSnap, int numSnap, int num_inf, int num_site, int num_act) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap);
		NUM_INF = num_inf;
		NUM_SITE = num_site;
		NUM_ACT = num_act;
		NUM_GENDER = pop_composition.length;

		// Initiate transient field, lookup table etc
		map_trans_prob = new HashMap<>();
		map_currrent_infection_stage = new HashMap<>();
		lookupTable_infection_stage_path = new HashMap<>();
		lookupTable_infection_infectious_stages = new int[NUM_INF][NUM_SITE];
		lookupTable_testing_properties = new HashMap<>();
		map_infection_stage_switch = new HashMap<>();
		map_currently_infectious = new HashMap<>();

		schedule_testing = new HashMap<>();
		schedule_stage_change = new HashMap<>();

		table_act_frequency = new double[NUM_ACT][NUM_GENDER][NUM_GENDER][];

		dist_tranmissionMatrix = new RealDistribution[NUM_INF][NUM_SITE][NUM_SITE][];

		dist_stage_period = new RealDistribution[NUM_INF][NUM_SITE][];
		dist_sym_rate = new RealDistribution[NUM_INF][NUM_SITE][];
	}

	public void refreshField(int fieldId, boolean clearAll) {
		double[][] field = null;

		try {
			field = (double[][]) getRunnable_fields()[fieldId];
		} catch (ClassCastException ex) {
			field = null;
		}

		// Reset runnable field

		switch (fieldId) {
		case RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ:
			if (clearAll) {
				for (double[][][] actType : table_act_frequency) {
					for (double[][] actType_gIf : actType) {
						Arrays.fill(actType_gIf, null);
					}
				}
			}
			for (double[] ent : field) {
				int act_type = (int) ent[FIELD_ACT_FREQ_ACT_TYPE];
				int gI_from = (int) ent[FIELD_ACT_FREQ_GENDER_INCLUDE_INDEX_FROM];
				int gI_to = (int) ent[FIELD_ACT_FREQ_GENDER_INCLUDE_INDEX_TO];
				for (int gf = 0; gf < NUM_GENDER; gf++) {
					if ((gI_from & 1 << gf) != 0) {
						for (int gt = 0; gt < NUM_GENDER; gt++) {
							if ((gI_to & 1 << gt) != 0) {
								table_act_frequency[act_type][gf][gt] = ent;
							}
						}
					}
				}
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE:
			if (clearAll) {
				for (RealDistribution[][][] inf : dist_tranmissionMatrix) {
					for (RealDistribution[][] inf_sf : inf) {
						Arrays.fill(inf_sf, null);
					}
				}

			}
			for (double[] ent : field) {
				int inf_id = (int) ent[FIELD_TRANSMISSION_RATE_INF_ID];
				int state_include = (int) ent[FIELD_TRANSMISSION_RATE_STATE_INCLUDE_INDEX];
				int site_from = (int) ent[FIELD_TRANSMISSION_RATE_SITE_FROM];
				int site_to = (int) ent[FIELD_TRANSMISSION_RATE_SITE_TO];
				double[] param = Arrays.copyOfRange(ent, FIELD_TRANSMISSION_RATE_TRANS_PARAM_START, ent.length);

				int comparator = state_include;
				int statePt = 0;

				while (comparator != 0) {
					if ((comparator & 1) != 0) {
						if (dist_tranmissionMatrix[inf_id][site_from][site_to] == null) {
							dist_tranmissionMatrix[inf_id][site_from][site_to] = new RealDistribution[statePt + 1];
						} else if (!(statePt < dist_tranmissionMatrix[inf_id][site_from][site_to].length)) {
							dist_tranmissionMatrix[inf_id][site_from][site_to] = Arrays
									.copyOf(dist_tranmissionMatrix[inf_id][site_from][site_to], statePt + 1);
						}
						if (param.length == 1) {
							dist_tranmissionMatrix[inf_id][site_from][site_to][statePt] = generateNonDistribution(
									param);
						} else {
							// Beta distribution
							dist_tranmissionMatrix[inf_id][site_from][site_to][statePt] = generateBetaDistribution(
									param);
						}
					}
					comparator = comparator >> 1;
					statePt++;
				}

			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD:
			if (clearAll) {
				for (int[] inf : lookupTable_infection_infectious_stages) {
					Arrays.fill(inf, 0);
				}
			}
			for (double[] ent : field) {
				int inf_id = (int) ent[FIELD_INFECTIOUS_PERIOD_INF_ID];
				int site = (int) ent[FIELD_INFECTIOUS_PERIOD_SITE];
				lookupTable_infection_infectious_stages[inf_id][site] = (int) ent[FIELD_INFECTIOUS_PERIOD_STAGE_INCLUDE_INDEX];
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD:
			if (clearAll) {
				for (RealDistribution[][] inf : dist_stage_period) {
					Arrays.fill(inf, null);

				}
				lookupTable_infection_stage_path.clear();
			}

			HashMap<String, HashMap<Integer, Double>> map_infection_stage_switch; // Key="INF_ID,SITE_ID,STAGE_ID",V=Map(Stage->Prob)
			HashMap<Integer, Double> map_target_stage;
			map_infection_stage_switch = new HashMap<>();
			for (double[] ent : field) {
				int inf_id = (int) ent[FIELD_STAGE_PERIOD_INF_ID];
				int site = (int) ent[FIELD_STAGE_PERIOD_SITE];
				int stage_id_from = (int) ent[FIELD_STAGE_PERIOD_STATE_FROM];
				int stage_id_to = (int) ent[FIELD_STAGE_PERIOD_STATE_TO];
				double switch_prob = ent[FIELD_STAGE_PERIOD_SWTICH_STAGE_PROB];
				double[] param = Arrays.copyOfRange(ent, FIELD_STAGE_PERIOD_DURATION_PARAM_START, ent.length);
				String key = String.format("%d,%d,%d", inf_id, site, stage_id_from);
				map_target_stage = map_infection_stage_switch.get(key);
				if (map_target_stage == null) {
					map_target_stage = new HashMap<>();
					map_infection_stage_switch.put(key, map_target_stage);
				}
				map_target_stage.put(stage_id_to, switch_prob);
				if (dist_stage_period[inf_id][site] == null) {
					dist_stage_period[inf_id][site] = new RealDistribution[stage_id_to + 1];
				} else if (!(stage_id_to < dist_stage_period[inf_id][site].length)) {
					dist_stage_period[inf_id][site] = Arrays.copyOf(dist_stage_period[inf_id][site], stage_id_to + 1);
				}
				if (param.length == 1) {
					dist_stage_period[inf_id][site][stage_id_to] = generateNonDistribution(param);
				} else {
					if (param[1] < 0) {
						// Uniform distribution
						double[] param_u = Arrays.copyOf(param, param.length);
						param_u[1] = Math.abs(param_u[1]);
						dist_stage_period[inf_id][site][stage_id_to] = generateUniformDistribution(param_u);
					} else {
						// Gamma Distribution
						dist_stage_period[inf_id][site][stage_id_to] = generateGammaDistribution(param);
					}
				}
			}
			for (String key : map_infection_stage_switch.keySet()) {
				map_target_stage = map_infection_stage_switch.get(key);
				double[] cumul_prob = new double[map_target_stage.size() * 2];
				Integer[] target_stages = map_target_stage.keySet().toArray(new Integer[map_target_stage.size()]);
				Arrays.sort(target_stages);
				double prob_sum = 0;
				for (int i = 0; i < target_stages.length; i++) {
					prob_sum += map_target_stage.get(target_stages[i]).doubleValue();
				}
				double cumul_p = 0;
				for (int i = 0; i < target_stages.length; i++) {
					cumul_prob[i] = cumul_p + map_target_stage.get(target_stages[i]).doubleValue() / prob_sum;
					cumul_p += cumul_prob[i];
					cumul_prob[i + target_stages.length] = target_stages[i].intValue();
				}
				lookupTable_infection_stage_path.put(key, cumul_prob);
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_SYM_RATE:
			if (clearAll) {
				for (RealDistribution[][] inf : dist_sym_rate) {
					Arrays.fill(inf, null);
				}
			}
			for (double[] ent : field) {
				int inf_Id = (int) ent[FIELD_SYM_RATE_INF_ID];
				int site = (int) ent[FIELD_SYM_RATE_SITE];
				int state = (int) ent[FIELD_SYM_RATE_STATE_ID];
				double[] param = Arrays.copyOfRange(ent, FIELD_SYM_RATE_PROB, ent.length);

				if (dist_sym_rate[inf_Id][site] == null) {
					dist_sym_rate[inf_Id][site] = new RealDistribution[state + 1];
				} else if (dist_sym_rate[inf_Id][site].length < state) {
					dist_sym_rate[inf_Id][site] = Arrays.copyOf(dist_sym_rate[inf_Id][site], state + 1);
				}

				if (param.length == 1) {
					dist_sym_rate[inf_Id][site][state] = generateNonDistribution(param);
				} else {
					dist_sym_rate[inf_Id][site][state] = generateBetaDistribution(param);
				}
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES:
			if (clearAll) {
				lookupTable_testing_properties.clear();
			}
			for (double[] ent : field) {
				int inf_Id = (int) ent[FIELD_DX_TEST_PROPERTIES_INF_ID];
				int site = (int) ent[FIELD_DX_TEST_PROPERTIES_SITE];
				String key = String.format("%d,%d", inf_Id, site);
				lookupTable_testing_properties.put(key, ent);
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS:
		case RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES:
		case RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM:
			// Do nothing
			break;
		default:
			System.err.printf("Warning: refreshField(fieldId=%d) not defined.\n", fieldId);
			if (field != null) {
				System.err.println("Existing field:");
				for (double[] ent : field) {
					System.err.println(Arrays.toString(ent));
				}
			}

		}

	}

	@Override
	public Object[] getRunnable_fields() {
		return runnable_fields;
	}

	@Override
	public void initialse() {
		map_trans_prob.clear();
		map_currrent_infection_stage.clear();
		map_infection_stage_switch.clear();

		map_currently_infectious.clear();

		schedule_testing.clear();
		schedule_stage_change.clear();

		for (int r = 0; r < runnable_fields.length; r++) {
			refreshField(r, true);
		}

		sim_output = new HashMap<>();
	}

	@Override
	public void run() {

		// Check if there is a switch in prop
		Integer[] switchTime = propSwitch_map.keySet().toArray(new Integer[propSwitch_map.size()]);
		Arrays.sort(switchTime);
		int switchTimeIndex = 0;

		// Current contact map
		ContactMap cMap = new ContactMap();
		Integer[][] edges_array = getEdgesArrayFromBaseConctactMap();
		int edges_array_pt = 0;
		HashMap<Integer, ArrayList<Integer[]>> removeEdges = new HashMap<>();
		edges_array_pt = initaliseCMap(cMap, edges_array, edges_array_pt, firstSeedTime, removeEdges);

		// Pre allocate risk categories (mainly form MSM)
		setPreAllocatedRiskFromFile();

		int startTime = firstSeedTime;
		// Schedule testing
		for (Integer personId : bASE_CONTACT_MAP.vertexSet()) {
			scheduleNextTest(personId, startTime);
		} // End of schedule test

		int snap_index = 0;
		boolean hasInfectious = hasInfectious();

		ArrayList<int[]> infected_today = new ArrayList<>();

		int[][][] cumul_incidence_by_site = new int[NUM_INF][NUM_GENDER][NUM_SITE];
		int[][] cumul_incidence_by_person = new int[NUM_INF][NUM_GENDER];
		int[][] cumul_treatment_by_person = new int[NUM_INF][NUM_GENDER];

		HashMap<String, int[]> acted_today = new HashMap<>();
		// K="PID,PID", V=int[act] =
		// if SITE_P1 != SITE_P2
		// = 0 if no acts
		// = 1 if p1s1 -> p2s2, e.g. p1 insertive, or if the site is the same
		// = 2 if p2s1 -> p1s2 e.g. p2 insertive

		for (int currentTime = startTime; currentTime < startTime + nUM_TIME_STEPS_PER_SNAP * nUM_SNAP
				&& hasInfectious; currentTime++) {

			infected_today.clear();

			if (switchTimeIndex < switchTime.length && switchTime[switchTimeIndex] == currentTime) {
				HashMap<Integer, String> switch_ent = propSwitch_map.get(currentTime);
				for (Integer switch_index : switch_ent.keySet()) {
					String str_obj = switch_ent.get(switch_index);
					int fieldId = switch_index - RUNNABLE_OFFSET;
					getRunnable_fields()[fieldId] = PropValUtils.propStrToObject(str_obj,
							getRunnable_fields()[fieldId].getClass());
					refreshField(fieldId, false);
				}
				switchTimeIndex++;
			}

			edges_array_pt = updateCMap(cMap, currentTime, edges_array, edges_array_pt, removeEdges);

			// Check for stage change
			ArrayList<ArrayList<ArrayList<Integer>>> state_change_today = schedule_stage_change.remove(currentTime);
			if (state_change_today != null) {
				for (int inf_id = 0; inf_id < NUM_INF; inf_id++) {
					if (state_change_today.get(inf_id) != null && state_change_today.get(inf_id).size() > 0) {
						for (int site_id = 0; site_id < NUM_SITE; site_id++) {
							if (state_change_today.get(inf_id).get(site_id) != null
									&& state_change_today.get(inf_id).get(site_id).size() > 0) {
								ArrayList<Integer> state_change_pids = state_change_today.get(inf_id).get(site_id);

								// Update infection state
								for (Integer pid : state_change_pids) {
									int[][] current_stage_arr = map_currrent_infection_stage.get(pid);
									int[][] infection_state_switch = map_infection_stage_switch.get(pid);
									updateInfectionStage(pid, inf_id, site_id, current_stage_arr[inf_id][site_id],
											currentTime, current_stage_arr, infection_state_switch, 0);
								}
							}
						}
					}
				}
			}

			// Transmission
			String[] currently_infectious_keys = map_currently_infectious.keySet()
					.toArray(new String[map_currently_infectious.size()]);
			Arrays.sort(currently_infectious_keys);

			acted_today.clear();

			for (String key : currently_infectious_keys) {
				ArrayList<Integer> currenty_infectious_ent = map_currently_infectious.get(key);

				String[] key_s = key.split(",");
				int inf_id = Integer.parseInt(key_s[0]);
				int infected_site_id = Integer.parseInt(key_s[1]);

				for (Integer pid_inf_src : currenty_infectious_ent) {

					if (cMap.containsVertex(pid_inf_src)) {
						int g_s = getGenderType(pid_inf_src);
						Integer[][] edges = cMap.edgesOf(pid_inf_src).toArray(new Integer[0][]);

						for (Integer[] e : edges) {
							int pid_inf_tar = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1].equals(pid_inf_src)
									? e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]
									: e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];

							int g_t = getGenderType(pid_inf_tar);

							int[] partners = new int[] { pid_inf_src, pid_inf_tar };
							Arrays.sort(partners);

							int inf_src_pt = partners[0] == pid_inf_src ? 0 : 1;

							int[] hasActed = acted_today.get(Arrays.toString(partners));

							if (hasActed == null) {
								hasActed = new int[NUM_ACT];
								int anyAct = 0;
								for (int a = 0; a < NUM_ACT; a++) {
									double[] fieldEntry = table_act_frequency[a][g_s][g_t];
									if (fieldEntry[FIELD_ACT_FREQ_ACT_PER_DAY] >= 0) {
										if (RNG.nextDouble() < fieldEntry[FIELD_ACT_FREQ_ACT_PER_DAY]) {
											if (fieldEntry[FIELD_ACT_FREQ_SITE_P1] == fieldEntry[FIELD_ACT_FREQ_SITE_P2]) {
												if (!((FIELD_ACT_FREQ_CONDOM_EFFICACY < fieldEntry.length) && (RNG
														.nextDouble() < fieldEntry[FIELD_ACT_FREQ_CONDOM_EFFICACY]
																* (e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] > 1
																		? fieldEntry[FIELD_ACT_FREQ_USAGE_REG]
																		: fieldEntry[FIELD_ACT_FREQ_USAGE_CASUAL])))) {
													hasActed[a] = 1;
												}
											} else {
												// Check for receptiveness
												int p_1_risk_grp_inc = (int) fieldEntry[FIELD_ACT_FREQ_P1_RISK_GRP_INCLUDE_INDEX];

												// Valid risk group of the insertive partner
												int valid_p1 = 0;
												for (int p = 0; p < partners.length; p++) {
													int risk_p = risk_cat_map.get(partners[p]);
													if ((p_1_risk_grp_inc & 1 << risk_p) != 0) {
														valid_p1 |= 1 << p;
													}
												}
												if (!((FIELD_ACT_FREQ_CONDOM_EFFICACY < fieldEntry.length) && (RNG
														.nextDouble() < fieldEntry[FIELD_ACT_FREQ_CONDOM_EFFICACY]
																* (e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] > 1
																		? fieldEntry[FIELD_ACT_FREQ_USAGE_REG]
																		: fieldEntry[FIELD_ACT_FREQ_USAGE_CASUAL])))) {
													switch (valid_p1) {
													case 1:
														hasActed[a] = 1;
														break;
													case 2:
														hasActed[a] = 2;
														break;
													default:
														hasActed[a] = RNG.nextInt(2) + 1;
													}
												}
											}
										}
									}
									anyAct += hasActed[a];
								}

								for (int a = 0; a < NUM_ACT; a++) {
									double[] fieldEntry = table_act_frequency[a][g_s][g_t];
									if (fieldEntry[FIELD_ACT_FREQ_ACT_PER_DAY] < 0) {
										hasActed[a] = (anyAct > 0) ? 1 : 0;
									}
								}
								acted_today.put(Arrays.toString(partners), hasActed);
							}

							for (int a = 0; a < NUM_ACT; a++) {
								double[] fieldEntry = table_act_frequency[a][g_s][g_t];
								for (int s = 0; s < 2; s++) { // if s==0 p_1=site_1, else p_1=site_2
									if ((hasActed[a] & 1 << s) != 0) {
										int p1_site, p2_site;
										p1_site = (int) (s == 0 ? fieldEntry[FIELD_ACT_FREQ_SITE_P1]
												: fieldEntry[FIELD_ACT_FREQ_SITE_P2]);
										p2_site = (int) (s == 0 ? fieldEntry[FIELD_ACT_FREQ_SITE_P2]
												: fieldEntry[FIELD_ACT_FREQ_SITE_P1]);
										int src_site = inf_src_pt == 0 ? p1_site : p2_site;
										if (src_site == infected_site_id) {
											int tar_site = inf_src_pt == 0 ? p2_site : p1_site;
											int[][] tar_infection_stages = map_currrent_infection_stage
													.get(pid_inf_tar);

											if (tar_infection_stages == null
													|| tar_infection_stages[inf_id][tar_site] == AbstractIndividualInterface.INFECT_S) {

												int[][] src_infection_stages = map_currrent_infection_stage
														.get(pid_inf_src);
												int src_stage = src_infection_stages[inf_id][src_site];
												double[][][][] trans_prob = map_trans_prob.get(pid_inf_src);

												if (trans_prob[inf_id][src_site][tar_site][src_stage] > 0) {
													if (RNG.nextDouble() < trans_prob[inf_id][src_site][tar_site][src_stage]) {
														// Transmission successful
														cumul_incidence_by_site[inf_id][g_t][tar_site]++;

														boolean newInfection = true;
														if (tar_infection_stages == null) {
															tar_infection_stages = new int[NUM_INF][NUM_SITE];
															for (int[] stage_by_infection : tar_infection_stages) {
																Arrays.fill(stage_by_infection,
																		AbstractIndividualInterface.INFECT_S);
															}
															map_currrent_infection_stage.put(pid_inf_tar,
																	tar_infection_stages);
														} else {
															for (int siteId = 0; siteId < NUM_SITE
																	&& newInfection; siteId++) {
																newInfection &= tar_infection_stages[inf_id][siteId] == AbstractIndividualInterface.INFECT_S;
															}
														}
														if (newInfection) {
															cumul_incidence_by_person[inf_id][g_t]++;
														}

														tar_infection_stages[inf_id][tar_site] = STAGE_ID_JUST_INFECTED;
														int[][] tar_infection_state_switch = map_infection_stage_switch
																.get(pid_inf_tar);
														if (tar_infection_state_switch == null) {
															tar_infection_state_switch = new int[NUM_INF][NUM_SITE];
															map_infection_stage_switch.put(pid_inf_tar,
																	tar_infection_state_switch);
														}

														updateInfectStageChangeSchedule(pid_inf_tar, inf_id, tar_site,
																currentTime + 1,
																tar_infection_state_switch[inf_id][tar_site]);

													} // End of successful transmission
												} // End of possible transmission check
											} // End of checking target site as susceptible
										} // End of infectious src site
									} // End of specifying site for p1 and p2
								} // End of determining site of p1 and p2 through act
							} // End of checking all acts
						} // End of checking all edge of infectious
					} // End of checking infectious in cMap
				} // End of checking one infectious
			} // End loop for all infectious infection-site combinations

			// Testing
			ArrayList<int[]> testToday = schedule_testing.remove(currentTime);
			if (testToday != null) {
				for (int[] testing_stat : testToday) {
					// Testing stat: int[]{PID,INF_INCLUDE_INDEX,SITE_INCLUDE_INDEX}
					int pid = testing_stat[0];
					int infIncl = testing_stat[1];
					int siteIncl = testing_stat[2];
					for (int infId = 0; infId < NUM_INF; infId++) {
						if ((infIncl & 1 << infId) != 0) {

							boolean applyTreatment = false;
							double[] test_properies;
							int tested_stage_inc;
							int[][] inf_stage = null;

							for (int siteId = 0; siteId < NUM_SITE && !applyTreatment; siteId++) {
								if ((siteIncl & 1 << siteId) != 0) {
									// Test for the site
									test_properies = lookupTable_testing_properties
											.get(String.format("%d,%d", infId, siteId));
									if (test_properies != null) {
										inf_stage = map_currrent_infection_stage.get(pid);
										if (inf_stage != null && inf_stage[infId][siteId] >= 0) {
											tested_stage_inc = (int) test_properies[FIELD_DX_TEST_PROPERTIES_STAGE_INC];
											if ((tested_stage_inc & 1 << inf_stage[infId][siteId]) != 0) {
												double testSensitivity = test_properies[FIELD_DX_TEST_PROPERTIES_ACCURACY];
												if (testSensitivity > 0) {
													applyTreatment |= RNG.nextDouble() < testSensitivity;
												}
											}
										}
									}
								} // End of testing pid for infection and site
							}
							if (applyTreatment) {
								cumul_treatment_by_person[infId][getGenderType(pid)]++;
								int[][] infection_switch = map_infection_stage_switch.get(pid);
								for (int siteId = 0; siteId < NUM_SITE; siteId++) {
									test_properies = lookupTable_testing_properties
											.get(String.format("%d,%d", infId, siteId));
									if (test_properies != null) {
										tested_stage_inc = (int) test_properies[FIELD_DX_TEST_PROPERTIES_STAGE_INC];
										if ((tested_stage_inc & 1 << inf_stage[infId][siteId]) != 0) {
											int nextStage = (int) test_properies[FIELD_DX_TEST_PROPERTIES_TREATMENT_SUC_STAGE];
											inf_stage[infId][siteId] = nextStage;
											updateInfectStageChangeSchedule(pid, infId, siteId, currentTime + 1,
													infection_switch[infId][siteId]);
										}
									}
								}
							} // End of apply for treatment
						} // End of testing for infection
					} // Check if infection is tested
				} // End of testing for single individual
			} // End of all scheduled testing for today

			// Storing of outputs
			if (snap_index == 0) {
				String key;

				if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE) != 0) {
					key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_GENDER_COUNT,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> current_infected_gender_count = (HashMap<Integer, int[]>) sim_output
							.get(key);
					if (current_infected_gender_count == null) {
						current_infected_gender_count = new HashMap<>();
						sim_output.put(key, current_infected_gender_count);
					}

					key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_SITE_COUNT,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> current_infected_site_count = (HashMap<Integer, int[]>) sim_output.get(key);
					if (current_infected_site_count == null) {
						current_infected_site_count = new HashMap<>();
						sim_output.put(key, current_infected_site_count);
					}

					int[] num_inf_count = new int[NUM_INF * NUM_SITE];
					int[] num_inf_by_gender = new int[NUM_INF * NUM_GENDER];

					current_infected_gender_count.put(currentTime, num_inf_by_gender);
					current_infected_site_count.put(currentTime, num_inf_count);

					int pt = 0;
					for (int i = 0; i < NUM_INF; i++) {
						ArrayList<Integer> alreadyCounted = new ArrayList<>();
						for (int s = 0; s < NUM_SITE; s++) {
							String infKey = String.format("%d,%d", i, s);
							ArrayList<Integer> arr = map_currently_infectious.get(infKey);
							if (arr != null) {
								num_inf_count[pt] = arr.size();
								for (Integer pid : arr) {
									int pt_inf = Collections.binarySearch(alreadyCounted, pid);
									if (pt_inf < 0) {
										num_inf_by_gender[i* NUM_GENDER + getGenderType(pid)]++;
										alreadyCounted.add(~pt_inf, pid);
									}

								}

							}
							pt++;
						}
					}

					key = String.format(SIM_OUTPUT_KEY_INFECTED_SITE_STAGE_COUNT,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);

					@SuppressWarnings("unchecked")
					HashMap<Integer, HashMap<String, Integer>> current_infected_site_stage_count = (HashMap<Integer, HashMap<String, Integer>>) sim_output
							.get(key);

					if (current_infected_site_stage_count == null) {
						current_infected_site_stage_count = new HashMap<>();
						sim_output.put(key, current_infected_site_stage_count);
					}

					HashMap<String, Integer> infect_site_stage_count_current = new HashMap<>();
					current_infected_site_stage_count.put(currentTime, infect_site_stage_count_current);

					for (Integer pid : map_currrent_infection_stage.keySet()) {
						int[][] inf_stage = map_currrent_infection_stage.get(pid);
						for (int i = 0; i < NUM_INF; i++) {
							for (int s = 0; s < NUM_SITE; s++) {
								String count_key = String.format("%d,%d,%d", i, s, inf_stage[i][s]);
								Integer ent = infect_site_stage_count_current.get(count_key);
								if (ent == null) {
									ent = 0;
								}
								infect_site_stage_count_current.put(count_key, ent + 1);
							}
						}
					}
				}

				if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE) != 0) {
					// int[][] cumul_treatment_by_person = new int[NUM_INF][NUM_GENDER];
					key = String.format(SIM_OUTPUT_KEY_CUMUL_TREATMENT,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE);
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> cumulative_treatment = (HashMap<Integer, int[]>) sim_output.get(key);
					if (cumulative_treatment == null) {
						cumulative_treatment = new HashMap<>();
						sim_output.put(key, cumulative_treatment);
					}
					int[] ent = new int[NUM_INF * NUM_GENDER];
					cumulative_treatment.put(currentTime, ent);

					int pt = 0;
					for (int i = 0; i < NUM_INF; i++) {
						for (int g = 0; g < NUM_GENDER; g++) {
							ent[pt] = cumul_treatment_by_person[i][g];
							pt++;
						}
					}

				}

				if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE) != 0) {
					key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);

					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> cumulative_incidence_by_person = (HashMap<Integer, int[]>) sim_output
							.get(key);
					if (cumulative_incidence_by_person == null) {
						cumulative_incidence_by_person = new HashMap<>();
						sim_output.put(key, cumulative_incidence_by_person);
					}
					int[] ent = new int[NUM_INF * NUM_GENDER];
					cumulative_incidence_by_person.put(currentTime, ent);

					int pt = 0;
					for (int i = 0; i < NUM_INF; i++) {
						for (int g = 0; g < NUM_GENDER; g++) {
							ent[pt] = cumul_incidence_by_person[i][g];
							pt++;
						}
					}
					key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE_SITE,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);

					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> cumulative_incidence_by_site = (HashMap<Integer, int[]>) sim_output
							.get(key);
					if (cumulative_incidence_by_site == null) {
						cumulative_incidence_by_site = new HashMap<>();
						sim_output.put(key, cumulative_incidence_by_site);
					}
					ent = new int[NUM_INF * NUM_GENDER * NUM_SITE];
					cumulative_incidence_by_site.put(currentTime, ent);
					pt = 0;
					for (int i = 0; i < NUM_INF; i++) {
						for (int g = 0; g < NUM_GENDER; g++) {
							for (int s = 0; s < NUM_SITE; s++) {
								ent[pt] = cumul_incidence_by_site[i][g][s];
								pt++;
							}
						}
					}
				}
				if (print_progress != null && runnableId != null) {
					try {
						print_progress.printf("Thread <%s>: t = %d . Timestamp = %tc.\n", runnableId, currentTime,
								System.currentTimeMillis());
					} catch (Exception ex) {
						System.err.printf("Thread <%s>: t = %d .\n", runnableId, currentTime);
					}
				}
			}

			snap_index = (snap_index + 1) % nUM_TIME_STEPS_PER_SNAP;
			hasInfectious = hasInfectious();

		} // End of time step

		if (runnableId != null) {
			System.out.printf("Thread <%s> completed.\n", runnableId);
		}

		// Post simulation
		postSimulation();

	}

	@SuppressWarnings("unchecked")
	protected void postSimulation() {
		String key, fileName;
		HashMap<Integer, int[]> countMap;
		String filePrefix = this.getRunnableId() == null ? "" : this.getRunnableId();
		
		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE) != 0) {
			key = String.format(SIM_OUTPUT_KEY_CUMUL_TREATMENT, Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE);

			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER });

		}
		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE) != 0) {
			key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE, Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER });

			key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE_SITE, Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);

			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_SITE,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d_Site_%d", new int[] { NUM_INF, NUM_GENDER, NUM_SITE });

		}

		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE) != 0) {
			
			key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_GENDER_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER });			
			
			
			key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_SITE_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Site_%d", new int[] { NUM_INF, NUM_SITE });

			// TODO: To check timing 
			
			/*
			key = String.format(SIM_OUTPUT_KEY_INFECTED_SITE_STAGE_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			
			HashMap<Integer, HashMap<String, Integer>> infect_site_stage_count = (HashMap<Integer, HashMap<String, Integer>>) sim_output
					.get(key);

			Comparator<String> stateKeyComp = new Comparator<String>() {
				@Override
				public int compare(String o1, String o2) {
					String[] s1 = o1.split(",");
					String[] s2 = o2.split(",");
					int res = 0;
					int pt = 0;
					while (res == 0 && pt < Math.min(s1.length, s2.length)) {
						res = Integer.compare(Integer.parseInt(s1[pt]), Integer.parseInt(s2[pt]));
					}
					return res;
				}
			};

			ArrayList<String> headerKey = new ArrayList<>();
			Integer[] timeArr = infect_site_stage_count.keySet().toArray(new Integer[0]);
			Arrays.sort(timeArr);

			// Check max state
			for (Integer time : timeArr) {
				HashMap<String, Integer> infect_site_stage_count_current = infect_site_stage_count.get(time);
				String[] state_key_arr = infect_site_stage_count_current.keySet().toArray(new String[0]);
				for (String state_key : state_key_arr) {
					int pt = Collections.binarySearch(headerKey, state_key, stateKeyComp);
					if (pt < 0) {
						headerKey.add(~pt, state_key);
					}
				}
			}
			// Print entry
			PrintWriter pWri;
			StringWriter sWri = null;
			try {
				FileOutputStream fOut = new FileOutputStream(new File(baseDir,
						String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE,
								cMAP_SEED, sIM_SEED)),true);
				pWri = new PrintWriter(fOut);
			} catch (FileNotFoundException e) {
				e.printStackTrace(System.err);
				sWri = new StringWriter();
				pWri = new PrintWriter(sWri);
			}

			pWri.println("Infection-Site-Stage-Count");
			pWri.print("Time");
			for (String header : headerKey) {
				pWri.print(',');
				pWri.print(header);
			}
			pWri.println();

			for (Integer time : timeArr) {
				pWri.print(time);
				HashMap<String, Integer> infect_site_stage_count_current = infect_site_stage_count.get(time);
				for (String header : headerKey) {
					Integer count = infect_site_stage_count_current.get(header);
					if (count == null) {
						count = 0;
					}
					pWri.print(',');
					pWri.print(count);
				}
				pWri.println();
			}

			pWri.close();

			if (sWri != null) {
				System.out.println(sWri.toString());
			}
			
			*/

		}

		
		
		if (print_progress != null && runnableId != null) {
			try {
				print_progress.printf("Post simulation file generation for Thread <%s> completed. Timestamp = %tc.\n", runnableId,
						System.currentTimeMillis());
			} catch (Exception ex) {
				System.err.printf("Post simulation file generation for Thread <%s> completed.\n", runnableId);
			}
		}
	}

	public void printCountMap(HashMap<Integer, int[]> countMap, String fileName, String headerFormat, int[] dimension) {

		PrintWriter pWri;
		StringWriter s = null;
		try {
			pWri = new PrintWriter(new File(baseDir, fileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace(System.err);
			s = new StringWriter();
			pWri = new PrintWriter(s);
		}

		// Header
		pWri.print("Time");
		recursiveHeaderGeneration(dimension, 0, headerFormat, pWri);
		pWri.println();

		// Entry
		Integer[] timePt = countMap.keySet().toArray(new Integer[countMap.size()]);
		Arrays.sort(timePt);
		for (Integer time : timePt) {
			int[] ent = countMap.get(time);
			pWri.print(time);
			for (int i = 0; i < ent.length; i++) {
				pWri.print(',');
				pWri.print(ent[i]);
			}
			pWri.println();
		}

		pWri.close();

		if (s != null) {
			System.out.println(s.toString());
		}

	}

	private void recursiveHeaderGeneration(int[] dimension, int dPt, String seedStr, PrintWriter pWri) {
		for (int i = 0; i < dimension[dPt]; i++) {
			if (dPt == dimension.length - 1) {
				pWri.print(',');
				pWri.print(seedStr.replaceFirst("%d", Integer.toString(i)));
			} else {
				recursiveHeaderGeneration(dimension, dPt + 1, seedStr.replaceFirst("%d", Integer.toString(i)), pWri);
			}
		}
	}

	public void scheduleNextTest(Integer personId, int lastTestTime) {
		int genderType = getGenderType(personId);
		int riskCat;
		if (risk_cat_map.containsKey(personId)) {
			riskCat = risk_cat_map.get(personId);
		} else {
			riskCat = -1;
			// Based on number of casual partners
			Set<Integer[]> edges = bASE_CONTACT_MAP.edgesOf(personId);
			int numCasual = 0;
			int firstPartnerTime = Integer.MAX_VALUE;
			int lastPartnerTime = 0;
			for (Integer[] e : edges) {
				if (e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] >= firstSeedTime
						&& e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] <= 1) {
					firstPartnerTime = Math.min(firstPartnerTime,
							e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
					lastPartnerTime = Math.max(lastPartnerTime,
							e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
					numCasual++;
				}
			}
			double numCasualPerYear = (((double) AbstractIndividualInterface.ONE_YEAR_INT) * numCasual)
					/ (lastPartnerTime - firstPartnerTime);

			double[][] riskCatDefs = (double[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS];

			for (int i = 0; i < riskCatDefs.length && riskCat == -1; i++) {
				int gIncl = (int) riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_GENDER_INCLUDE_INDEX];
				if ((gIncl & 1 << genderType) != 0) {
					if (riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_NUM_CASUAL_PARTNER_LOWER] < numCasualPerYear
							&& numCasualPerYear < riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_NUM_CASUAL_PARTNER_UPPER]) {
						riskCat = (int) riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_RISK_GRP_DEF_ID];
					}
				}

			}
			risk_cat_map.put(personId, riskCat);
		}

		if (riskCat != -1) {
			double[][] testRateDefs = (double[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES];
			for (int i = 0; i < testRateDefs.length; i++) {
				double[] testRateDef = testRateDefs[i];
				int gIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_GENDER_INCLUDE_INDEX];
				int iIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_INF_INCLUDE_INDEX];
				int rIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_RISK_GRP_INCLUDE_INDEX];
				if (((gIncl & (1 << genderType)) & (rIncl & (1 << riskCat))) != 0) {
					int sIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_SITE_INCLUDE_INDEX];
					// FORMAT: {Cumul_Prob_0, Cummul_Prob_1, .... test_gap_time_0, test_gap_time_1,
					// test_gap_time_2..}
					int cumul_end_index = FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START
							+ ((testRateDef.length - FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START) - 1)
									/ 2;
					int pt = Arrays.binarySearch(testRateDef,
							FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START, cumul_end_index,
							RNG.nextDouble());
					if (pt < 0) {
						pt = ~pt;
					}
					int nextTestAfter = (int) (testRateDef[cumul_end_index
							- FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START + pt + 1]
							+ RNG.nextInt((int) testRateDef[cumul_end_index
									- FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START + pt]
									- (int) testRateDef[cumul_end_index
											- FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START + pt + 1]));

					int nextTestDate = lastTestTime + nextTestAfter;

					ArrayList<int[]> day_sch = schedule_testing.get(nextTestDate);

					if (day_sch == null) {
						day_sch = new ArrayList<>();
						schedule_testing.put(nextTestDate, day_sch);
					}

					int[] test_pair = new int[] { personId, iIncl, sIncl };

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
						org_pair[1] |= iIncl;
					}
				}
			}
		}
	}

	protected boolean hasInfectious() {
		boolean hasInfected = false;
		for (String k : map_currently_infectious.keySet()) {
			if (map_currently_infectious.get(k).size() != 0) {
				return true;
			}
		}

		return hasInfected;
	}

	@Override
	public void allocateSeedInfection(int[][] num_infected, int time) {
		// num_infected = int[infection_id]{GENDER_INC_INDEX_0, SITE_INDEX_0,
		// Number_INF_0,...}

		int lastPid = cUMULATIVE_POP_COMPOSITION[cUMULATIVE_POP_COMPOSITION.length - 1];
		firstSeedTime = Math.min(firstSeedTime, time);

		int inf_id = 0;
		ArrayList<Integer> candidate = new ArrayList<>();
		for (int[] inf_setting : num_infected) {
			int pt = 0;
			while (pt < inf_setting.length) {
				int genderInclude = inf_setting[pt];
				pt++;
				int site_index = inf_setting[pt];
				pt++;
				int num_inf = inf_setting[pt];
				pt++;

				int validStage = lookupTable_infection_infectious_stages[inf_id][site_index];

				if (validStage == 0) {
					System.err.printf("Warning: no infectious stage find for Inf_#%d at site #%d.\n", inf_id,
							site_index);
				}

				int stage_pt = 0;
				ArrayList<Integer> validInfectiousStages = new ArrayList<>();
				while (validStage != 0) {
					if ((validStage & 1) > 0) {
						validInfectiousStages.add(stage_pt);
					}
					validStage = validStage >> 1;
					stage_pt++;
				}

				if (validInfectiousStages.size() > 0) {
					// List out all valid candidate
					candidate.clear();
					for (int pid = 1; pid <= lastPid; pid++) {
						if ((genderInclude & 1 << getGenderType(pid)) > 0) {
							candidate.add(pid);
						}
					}

					int counter = 0;
					for (Integer pid : candidate) {
						if (RNG.nextInt(candidate.size() - counter) < num_inf) {
							int[] cumulProb = new int[validInfectiousStages.size()];
							for (int i = 0; i < cumulProb.length; i++) {
								if (i > 0) {
									cumulProb[i] = cumulProb[i - 1];
								}
								cumulProb[i] += (int) Math.ceil(
										dist_stage_period[inf_id][site_index][validInfectiousStages.get(i)].sample());
							}

							int stageDur = RNG.nextInt(cumulProb[cumulProb.length - 1]);
							int stage_sel = Arrays.binarySearch(cumulProb, stageDur);
							if (stage_sel < 0) {
								stage_sel = ~stage_sel;
							}
							validStage = validInfectiousStages.get(stage_sel);
							if (stage_sel > 0) {
								stageDur -= cumulProb[stage_sel - 1];
							}

							addInfectious(pid, inf_id, site_index, validStage, time, stageDur);
							num_inf--;
						}
						counter++;
					}
				}
			}
			inf_id++;
		}
	}

	@Override
	public int addInfectious(Integer infectedPId, int infectionId, int site_id, int stage_id, int infectious_time,
			int state_duration_adj) {

		double[][][][] trans_prob = map_trans_prob.get(infectedPId);

		if (trans_prob == null) {
			trans_prob = new double[NUM_INF][NUM_SITE][NUM_SITE][];
			for (int inf = 0; inf < trans_prob.length; inf++) {
				for (int s_f = 0; s_f < trans_prob[inf].length; s_f++) {
					for (int s_t = 0; s_t < trans_prob[inf][s_f].length; s_t++) {
						if (dist_tranmissionMatrix[inf][s_f][s_t] != null) {
							trans_prob[inf][s_f][s_t] = new double[dist_tranmissionMatrix[inf][s_f][s_t].length];
							for (int stateId = 0; stateId < trans_prob[inf][s_f][s_t].length; stateId++) {
								if (dist_tranmissionMatrix[inf][s_f][s_t][stateId] != null) {
									trans_prob[inf][s_f][s_t][stateId] = dist_tranmissionMatrix[inf][s_f][s_t][stateId]
											.sample();
								}
							}
						}
					}
				}

			}
			map_trans_prob.put(infectedPId, trans_prob);
		}

		int[][] current_stage_arr = map_currrent_infection_stage.get(infectedPId);
		if (current_stage_arr == null) {
			current_stage_arr = new int[NUM_INF][NUM_SITE];
			for (int[] stage_by_infection : current_stage_arr) {
				Arrays.fill(stage_by_infection, AbstractIndividualInterface.INFECT_S);
			}
			map_currrent_infection_stage.put(infectedPId, current_stage_arr);
		}
		int[][] infection_state_switch = map_infection_stage_switch.get(infectedPId);
		if (infection_state_switch == null) {
			infection_state_switch = new int[NUM_INF][NUM_SITE];
			map_infection_stage_switch.put(infectedPId, infection_state_switch);
		}

		int pre_ent = infection_state_switch[infectionId][site_id];

		ArrayList<ArrayList<ArrayList<Integer>>> schedule_inf;
		ArrayList<ArrayList<Integer>> schedule_inf_site;
		ArrayList<Integer> schedule_inf_site_arr;
		int pt;

		// Remove previous entry of the site if needed.
		if (pre_ent > 0) {
			schedule_inf = schedule_stage_change.get(pre_ent);
			if (schedule_inf != null) {
				schedule_inf_site = schedule_inf.get(infectionId);
				if (schedule_inf_site != null) {
					schedule_inf_site_arr = schedule_inf_site.get(site_id);
					if (schedule_inf_site_arr != null) {
						pt = Collections.binarySearch(schedule_inf_site_arr, infectedPId);
						if (pt >= 0) {
							schedule_inf_site_arr.remove(pt);
						}
					}
				}
			}
		}

		// Set infection state

		int current_state = stage_id;
		int current_infect_switch_time = updateInfectionStage(infectedPId, infectionId, site_id, current_state,
				infectious_time, current_stage_arr, infection_state_switch, state_duration_adj);

		return current_infect_switch_time;
	}

	// If state_duration_adj > 0, simple duration update, otherwise switch to next
	// stage and resample
	private int updateInfectionStage(Integer pid, int infection_id, int site_id, int current_state, int current_time,
			int[][] current_stage_arr, int[][] infection_state_switch, int state_duration_adj) {
		int pt;
		double state_duration;
		int infect_switch_time = current_time;

		if (current_state == STAGE_ID_JUST_INFECTED) {
			// Always to first stage
			current_state = 0;
			infect_switch_time = addInfectious(pid, infection_id, site_id, current_state, current_time,
					(int) dist_stage_period[infection_id][site_id][current_state].sample());
		} else {
			if (state_duration_adj <= 0) {
				while (infect_switch_time == current_time) {
					String key = String.format("%d,%d,%d", infection_id, site_id, current_state);
					double[] nextProb = lookupTable_infection_stage_path.get(key);
					if (nextProb == null) {
						// No next stage - return to susceptible
						current_state = AbstractIndividualInterface.INFECT_S;
						infect_switch_time = -1;
						break;
					} else {
						int state_pt = Arrays.binarySearch(nextProb, 0, nextProb.length / 2, RNG.nextDouble());
						if (state_pt <= 0) {
							state_pt = ~state_pt;
						}
						current_state = (int) nextProb[nextProb.length / 2 + state_pt];
						state_duration = dist_stage_period[infection_id][site_id][current_state].sample();
						infect_switch_time = (int) Math.round(current_time + state_duration);
					}
				}
			} else {
				infect_switch_time = current_time + state_duration_adj;
			}
		}

		// Update state_switch map
		current_stage_arr[infection_id][site_id] = current_state;
		infection_state_switch[infection_id][site_id] = infect_switch_time;

		if (infect_switch_time > current_time) {
			// Update schedule
			updateInfectStageChangeSchedule(pid, infection_id, site_id, infect_switch_time, -1);
		}

		// Update infectious status
		String key = String.format("%d,%d", infection_id, site_id);
		ArrayList<Integer> currenty_infectious_ent = map_currently_infectious.get(key);

		if (((1 << current_state) & lookupTable_infection_infectious_stages[infection_id][site_id]) != 0) {
			// Infectious stage
			if (currenty_infectious_ent == null) {
				currenty_infectious_ent = new ArrayList<>();
				map_currently_infectious.put(key, currenty_infectious_ent);
			}
			pt = Collections.binarySearch(currenty_infectious_ent, pid);
			if (pt < 0) {
				currenty_infectious_ent.add(~pt, pid);
			}
		} else {
			// Non infectious stage
			if (currenty_infectious_ent != null) {
				pt = Collections.binarySearch(currenty_infectious_ent, pid);
				if (pt >= 0) {
					currenty_infectious_ent.remove(pt);
				}
			}
		}
		return infect_switch_time;
	}

	private void updateInfectStageChangeSchedule(Integer pid, int infection_id, int site_id, int infect_switch_time,
			int org_infection_switch_time) {
		ArrayList<ArrayList<ArrayList<Integer>>> schedule_inf;
		ArrayList<ArrayList<Integer>> schedule_inf_site;
		ArrayList<Integer> schedule_inf_site_arr;
		int pt;

		if (org_infection_switch_time > 0) {
			schedule_inf = schedule_stage_change.get(org_infection_switch_time);
			if (schedule_inf != null) { // Remove original schedule
				schedule_inf_site_arr = schedule_inf.get(infection_id).get(site_id);
				pt = Collections.binarySearch(schedule_inf_site_arr, pid);
				if (pt > 0) {
					schedule_inf_site_arr.remove(pt);
				}
			}
		}
		schedule_inf = schedule_stage_change.get(infect_switch_time);
		if (schedule_inf == null) {
			schedule_inf = new ArrayList<>(NUM_INF);
			for (int i = 0; i < NUM_INF; i++) {
				schedule_inf_site = new ArrayList<>(NUM_SITE);
				for (int s = 0; s < NUM_SITE; s++) {
					schedule_inf_site.add(new ArrayList<>());
				}
				schedule_inf.add(schedule_inf_site);
			}
			schedule_stage_change.put(infect_switch_time, schedule_inf);
		}

		schedule_inf_site_arr = schedule_inf.get(infection_id).get(site_id);
		pt = Collections.binarySearch(schedule_inf_site_arr, pid);
		if (pt < 0) {
			schedule_inf_site_arr.add(~pt, pid);
		}
	}

}
