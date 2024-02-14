package sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.RealDistribution;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_MultiTransmission extends Abstract_Runnable_ClusterModel_Transmission {

	public Object[] runnable_fields = {
			// 0: RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ
			// double{{ACT_TYPE, GENDER_INCLUDE_INDEX_FROM, GENDER_INCLUDE_INDEX_TO,
			// ACT_PER_DAY},...}
			// or
			// double{{ACT_TYPE, GENDER_INCLUDE_INDEX_FROM, GENDER_INCLUDE_INDEX_TO,
			// ACT_PER_DAY, CONDOM_EFFICACY, USAGE_REG, USAGE_CASUAL},...}
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
			// 8 RUNNABLE_FIELD_TRANSMISSION_DX_TEST_ACCURACY
			// double{{GENDER_INCLUDE_INDEX, INF_INCLUDE_INDEX, SITE, TEST_ACCURACY}, ...}
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
	public static final int FIELD_ACT_FREQ_ACT_PER_DAY = FIELD_ACT_FREQ_GENDER_INCLUDE_INDEX_TO + 1;
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
	// For RUNNALBE_FIELD_TRANSMISSION_DX_TEST_ACCURACY
	public static final int FIELD_DX_TEST_ACCURACY_INF_ID = 0;
	public static final int FIELD_DX_TEST_ACCURACY_SITE = FIELD_DX_TEST_ACCURACY_INF_ID + 1;
	public static final int FIELD_TEST_ACCURACY = FIELD_DX_TEST_ACCURACY_SITE + 1;

	// Transient field
	protected transient HashMap<Integer, double[][][][]> map_trans_prob; // Key=PID,V=double[INF_ID][SITE_FROM][SITE_TO]
	protected transient HashMap<Integer, int[][]> map_currrent_infection_stage; // Key=PID,V=int[INF_ID][SITE]{infection_stage}
	protected transient HashMap<Integer, int[][]> map_infection_stage_switch; // Key=PID,V=int[inf_id]{switch_time_at_site_0};
	protected transient HashMap<String, ArrayList<Integer>> map_currently_infectious; // Key=Inf_ID_Site,V=ArrayList of
																						// pid with infectious

	// For schedule, key are day and entries are list of person id, ordered
	protected transient HashMap<Integer, ArrayList<int[]>> schedule_testing; // V=PID,INF_INCLUDE_INDEX
	protected transient HashMap<Integer, ArrayList<ArrayList<ArrayList<Integer>>>> schedule_stage_change;
	// V=ArrayList<Inf><SITE>{PIds of those need to change}

	// Helper objects set by fields

	protected transient HashMap<String, double[]> lookupTable_infection_stage_path; // Key="INF_ID,SITE_ID,STAGE_ID",V={Cumul_Prob_0,...
																					// State_ID_0...}
	protected transient int[][] lookupTable_infection_infectious_stages; // int[INF_ID][SITE]=STAGE_INCLUDE_ID

	protected transient double[][][][] table_act_frequency; // double[ACT_ID][G_TO][G_FROM]=fieldEntry

	protected transient RealDistribution[][][] dist_sym_rate; // RealDistribution[INF_ID][SITE][STAGE_ID]
	protected transient RealDistribution[][][][] dist_tranmissionMatrix; // RealDistribution[INF_ID][SITE_FROM][SITE_TO][STAGE_ID]
	protected transient RealDistribution[][][] dist_stage_period; // ReadDistribution[INF_ID][SITE][STAGE_ID]

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
		case RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS:
		case RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES:
		case RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM:
		case RUNNABLE_FIELD_TRANSMISSION_DX_TEST_ACCURACY:
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
						// FORMAT: {Cumul_Prob_0, Cummul_Prob_1, .... test_gap_time_0, test_gap_time_1,
						// test_gap_time_2..}
						int cumul_end_index = FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START
								+ ((testRateDef.length - FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START)
										- 1) / 2;
						int pt = Arrays.binarySearch(testRateDef, RNG.nextInt(),
								FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START, cumul_end_index);
						if (pt < 0) {
							pt = ~pt;
						}
						int nextTestAfter = (int) (testRateDef[cumul_end_index
								- FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START + pt + 1]
								+ RNG.nextInt((int) testRateDef[cumul_end_index
										- FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START + pt]
										- (int) testRateDef[cumul_end_index
												- FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START + pt
												+ 1]));

						int nextTestDate = startTime + nextTestAfter;

						ArrayList<int[]> day_sch = schedule_testing.get(nextTestDate);

						if (day_sch == null) {
							day_sch = new ArrayList<>();
							schedule_testing.put(nextTestDate, day_sch);
						}

						int[] test_pair = new int[] { personId, iIncl };

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
		} // End of schedule test

		int snap_index = 0;
		boolean hasInfectious = hasInfectious();

		ArrayList<int[]> infected_today = new ArrayList<>();

		int[][][] cumul_incidence_by_site = new int[NUM_INF][NUM_GENDER][NUM_SITE];
		int[][] cumul_incidence_by_person = new int[NUM_INF][NUM_GENDER];
		int[][] cumul_treatment_by_person = new int[NUM_INF][NUM_GENDER];
		int[][] cumul_positive_dx_test_by_person = new int[NUM_INF][NUM_GENDER];
		int[][] cumul_positive_dx_sought_by_person = new int[NUM_INF][NUM_GENDER];
		HashMap<String, boolean[]> acted_today = new HashMap<>(); // K="PID,PID", V=boolean[] = has_acted

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

				// Transmission
				String[] currently_infectious_keys = map_currently_infectious.keySet()
						.toArray(new String[map_currently_infectious.size()]);
				Arrays.sort(currently_infectious_keys);

				for (String key : currently_infectious_keys) {
					ArrayList<Integer> currenty_infectious_ent = map_currently_infectious.get(key);
					ArrayList<Integer> infectious_next_step = new ArrayList<>();

					String[] key_s = key.split(",");
					int inf_id = Integer.parseInt(key_s[0]);
					int site_id = Integer.parseInt(key_s[1]);

					acted_today.clear();

					for (Integer pid_inf_src : currenty_infectious_ent) {

						if (cMap.containsVertex(pid_inf_src)) {
							int g_s = getGenderType(pid_inf_src);
							Integer[][] edges = cMap.edgesOf(pid_inf_src).toArray(new Integer[0][]);

							for (Integer[] e : edges) {
								int pid_inf_tar = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]
										.equals(pid_inf_src) ? e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]
												: e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];

								int g_t = getGenderType(pid_inf_tar);

								int[] partners = new int[] { pid_inf_src, pid_inf_tar };
								Arrays.sort(partners);

								boolean[] hasActed = acted_today.get(Arrays.toString(partners));

								if (hasActed == null) {
									hasActed = new boolean[NUM_ACT];
									for (int a = 0; a < NUM_ACT; a++) {
										double[] fieldEntry = table_act_frequency[a][g_s][g_t];
										hasActed[a] = RNG.nextDouble() < fieldEntry[0];
										if (hasActed[a] && fieldEntry.length > 1) {
											// Partnership specific condom usage											
											hasActed[a] = RNG.nextDouble() >= (
													e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] > 1 ? 
														fieldEntry[FIELD_ACT_FREQ_USAGE_REG] : fieldEntry[FIELD_ACT_FREQ_USAGE_CASUAL]);
										}
									}
									acted_today.put(Arrays.toString(partners), hasActed);
								}

								// TODO: Determine if acts has occurred.
								
								
								
								

							}

						}

					}

					for (Integer pid_inf_next : infectious_next_step) {
						int pt = Collections.binarySearch(currenty_infectious_ent, pid_inf_next);
						if (pt < 0) {
							currenty_infectious_ent.add(~pt, pid_inf_next);
						}
					}
				}

			}

		} // End of time step

		// TODO: To be implement/Check

		System.out.println("Debug: Run completed.");

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
							validStage = validInfectiousStages.get(RNG.nextInt(validInfectiousStages.size()));
							addInfectious(pid, inf_id, site_index, validStage, time, 2);
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

	private int updateInfectionStage(Integer pid, int infection_id, int site_id, int current_state, int current_time,
			int[][] current_stage_arr, int[][] infection_state_switch, int state_duration_adj) {
		ArrayList<ArrayList<ArrayList<Integer>>> schedule_inf;
		ArrayList<ArrayList<Integer>> schedule_inf_site;
		ArrayList<Integer> schedule_inf_site_arr;
		int pt;
		double state_duration;
		int infect_switch_time = current_time;
		while (infect_switch_time == current_time) {
			String key = String.format("%d,%d,%d", infection_id, site_id, current_state);
			double[] nextProb = lookupTable_infection_stage_path.get(key);
			if (nextProb == null) {
				break;
			} else {
				int state_pt = Arrays.binarySearch(nextProb, 0, nextProb.length / 2, RNG.nextDouble());
				if (state_pt <= 0) {
					state_pt = ~state_pt;
				}
				current_state = (int) nextProb[nextProb.length / 2 + state_pt];
				state_duration = dist_stage_period[infection_id][site_id][current_state].sample();
				if (state_duration_adj != 0) {
					// Initial offset
					state_duration /= state_duration_adj;
				}

				infect_switch_time = (int) Math.round(current_time + state_duration);
			}
		}

		// Update state_switch map
		current_stage_arr[infection_id][site_id] = current_state;
		infection_state_switch[infection_id][site_id] = infect_switch_time;

		// Update schedule
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

}
