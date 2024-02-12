package sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.RealDistribution;

import relationship.ContactMap;

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
			// double{{INF_ID, SITE, INFECTIOUS_PERIOD_PARAM_0, ...}, ...}
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
	public static final int FIELD_INFECTIOUS_PERIOD_PARAM_START = FIELD_INFECTIOUS_PERIOD_SITE;
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
	protected transient HashMap<Integer, int[][]> map_currrent_infectState; // Key=PID,V=int[INF_ID][SITE]{infection_state}
	protected transient HashMap<Integer, int[][]> map_infection_state_switch; // Key=PID,V=int[inf_id]{switch_time_at_site_0};
	protected transient HashMap<String, ArrayList<Integer>> map_currently_infectious; // Key=Inf_ID_Site,V=ArrayList of
																						// pid with infectious

	// For schedule, key are day and entries are list of person id, ordered
	protected transient HashMap<Integer, ArrayList<Integer>> schedule_testing;
	protected transient HashMap<Integer, ArrayList<ArrayList<ArrayList<Integer>>>> schedule_state_change;
	// V=ArrayList<Inf><SITE>{PIds of those need to change}

	// Helper objects set by fields

	protected transient HashMap<String, double[]> lookupTable_infection_stage_path; // Key="INF_ID,SITE_ID,STAGE_ID",
	// V={Cumul_Prob_0,... State_ID_0...}

	protected transient double[][][][] table_act_frequency; // double[ACT_ID][G_TO][G_FROM]{fieldEntry}

	protected transient RealDistribution[][][] dist_sym_rate; // RealDistribution[INF_ID][SITE][STAGE_ID]
	protected transient RealDistribution[][][][] dist_tranmissionMatrix; // RealDistribution[INF_ID][SITE_FROM][SITE_TO][STAGE_ID]
	protected transient RealDistribution[][] dist_infectious_period; // RealDistribution[INF_ID][SITE]
	protected transient RealDistribution[][][] dist_stages_period; // ReadDistribution[INF_ID][SITE][STAGE_ID]

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
		map_currrent_infectState = new HashMap<>();
		lookupTable_infection_stage_path = new HashMap<>();
		map_infection_state_switch = new HashMap<>();
		map_currently_infectious = new HashMap<>();

		schedule_testing = new HashMap<>();
		schedule_state_change = new HashMap<>();

		table_act_frequency = new double[NUM_ACT][NUM_GENDER][NUM_GENDER][];

		dist_tranmissionMatrix = new RealDistribution[NUM_INF][NUM_SITE][NUM_SITE][];
		dist_infectious_period = new RealDistribution[NUM_INF][NUM_SITE];
		dist_stages_period = new RealDistribution[NUM_INF][NUM_SITE][];
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
								if (fieldId == RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ) {
									table_act_frequency[act_type][gf][gt] = ent;
								}
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

				while (comparator > 0) {
					if (comparator % 1 > 0) {
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
					comparator = comparator / 2;
					statePt++;
				}

			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD:
			if (clearAll) {
				for (RealDistribution[] inf : dist_infectious_period) {
					Arrays.fill(inf, null);
				}
			}
			for (double[] ent : field) {
				int inf_id = (int) ent[FIELD_INFECTIOUS_PERIOD_INF_ID];
				int site = (int) ent[FIELD_INFECTIOUS_PERIOD_SITE];
				double[] param = Arrays.copyOfRange(ent, FIELD_INFECTIOUS_PERIOD_PARAM_START, ent.length);
				if (param.length > 0) { // If param.length == 0, duration of infectious is determined elsewhere.
					if (param.length == 1) {
						dist_infectious_period[inf_id][site] = generateNonDistribution(param);
					} else {
						if (param[1] < 0) {
							// Uniform distribution
							double[] param_u = Arrays.copyOf(param, param.length);
							param_u[1] = Math.abs(param_u[1]);
							dist_infectious_period[inf_id][site] = generateUniformDistribution(param_u);
						} else {
							dist_infectious_period[inf_id][site] = generateGammaDistribution(param);
						}
					}
				}
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD:
			if (clearAll) {
				for (RealDistribution[][] inf : dist_stages_period) {
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
				if (dist_stages_period[inf_id][site] == null) {
					dist_stages_period[inf_id][site] = new RealDistribution[stage_id_to + 1];
				} else if (!(stage_id_to < dist_stages_period[inf_id][site].length)) {
					dist_stages_period[inf_id][site] = Arrays.copyOf(dist_stages_period[inf_id][site], stage_id_to + 1);
				}
				if (param.length == 1) {
					dist_stages_period[inf_id][site][stage_id_to] = generateNonDistribution(param);
				} else {
					if (param[1] < 0) {
						// Uniform distribution
						double[] param_u = Arrays.copyOf(param, param.length);
						param_u[1] = Math.abs(param_u[1]);
						dist_stages_period[inf_id][site][stage_id_to] = generateUniformDistribution(param_u);
					} else {
						// Gamma Distribution
						dist_stages_period[inf_id][site][stage_id_to] = generateGammaDistribution(param);
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
		map_currrent_infectState.clear();
		map_infection_state_switch.clear();

		map_currently_infectious.clear();

		schedule_testing.clear();
		schedule_state_change.clear();

		for (int r = 0; r < runnable_fields.length; r++) {
			refreshField(r, true);
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub

	}

	@Override
	public void allocateSeedInfection(int[][] num_infected, int time) {
		// num_infected = int[infection_id]{GENDER_INC_INDEX_0, SITE_INDEX_0,
		// Number_INF_0,...}
		int lastPid = cUMULATIVE_POP_COMPOSITION[cUMULATIVE_POP_COMPOSITION.length - 1];

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
						addInfectious(pid, inf_id, site_index, 0, time, -1);
						num_inf--;
					}
					counter++;
				}
			}
			inf_id++;
		}
	}

	@Override
	public int addInfectious(Integer infectedId, int infId, int site, int state, int infectious_time, int recoveredAt) {

		double[][][][] trans_prob = map_trans_prob.get(infectedId);

		if (trans_prob == null) {
			trans_prob = new double[NUM_INF][NUM_SITE][NUM_SITE][];
			for (int inf = 0; inf < trans_prob.length; inf++) {
				for (int s_f = 0; s_f < trans_prob[inf].length; s_f++) {
					for (int s_t = 0; s_t < trans_prob[inf][s_f].length; s_t++) {
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
			map_trans_prob.put(infectedId, trans_prob);
		}
		
		int[][] current_state_arr = map_currrent_infectState.get(infectedId);
		if(current_state_arr == null) {
			current_state_arr = new int[NUM_INF][NUM_SITE];
			map_currrent_infectState.put(infectedId, current_state_arr);
		}
		int[][] infection_state_switch = map_infection_state_switch.get(infectedId);
		if (infection_state_switch == null) {
			infection_state_switch = new int[NUM_INF][NUM_SITE];
			map_infection_state_switch.put(infectedId, infection_state_switch);
		}

		int pre_ent = infection_state_switch[infId][site];

		ArrayList<ArrayList<ArrayList<Integer>>> schedule_inf;
		ArrayList<ArrayList<Integer>> schedule_inf_site;
		ArrayList<Integer> schedule_inf_site_arr;

		// Remove previous entry of the site if needed.
		if (pre_ent > 0) {
			schedule_inf = schedule_state_change.get(pre_ent);
			if (schedule_inf != null) {
				schedule_inf_site = schedule_inf.get(infId);
				if (schedule_inf_site != null) {
					schedule_inf_site_arr = schedule_inf_site.get(site);
					if (schedule_inf_site_arr != null) {
						int pt = Collections.binarySearch(schedule_inf_site_arr, infectedId);
						if (pt >= 0) {
							schedule_inf_site_arr.remove(pt);
						}
					}
				}
			}
		}

		// Set infection state

		int current_state = state;
		int current_infect_switch_time = (int) Math.round(infectious_time + dist_stages_period[infId][site][current_state].sample());

		while (current_infect_switch_time == infectious_time) {
			String key = String.format("%d,%d,%d", infId, site, current_state);
			double[] nextProb = lookupTable_infection_stage_path.get(key);
			if (nextProb == null) {
				break;
			} else {
				int state_pt = Arrays.binarySearch(nextProb, 0, nextProb.length/2, RNG.nextDouble());
				if(state_pt <= 0) {
					state_pt = ~state_pt;
				}
				current_state = (int) nextProb[nextProb.length/2 + state_pt];				
				current_infect_switch_time = (int) Math.round(infectious_time + dist_infectious_period[infId][current_state].sample());
			}
		}
		
		
		// Update state_switch map
		current_state_arr[infId][site] = current_state;		
		infection_state_switch[infId][site] = current_infect_switch_time;

		// Update schedule
		schedule_inf = schedule_state_change.get(current_infect_switch_time);
		if (schedule_inf == null) {
			schedule_inf = new ArrayList<>(NUM_INF);
			for (int i = 0; i < NUM_INF; i++) {
				schedule_inf_site = new ArrayList<>(NUM_SITE);
				for (int s = 0; s < NUM_SITE; s++) {
					schedule_inf.add(new ArrayList<>());
				}
				schedule_inf.add(schedule_inf_site);
			}
			schedule_state_change.put(current_infect_switch_time, schedule_inf);
		}
		schedule_inf_site_arr = schedule_inf.get(infId).get(site);
		int pt = Collections.binarySearch(schedule_inf_site_arr, infectedId);
		if (pt < 0) {
			schedule_inf_site_arr.add(pt, infectedId);
		}

		// TODO: Update list_currently_infectious if need
		
		
		

		return current_infect_switch_time;
	}

}