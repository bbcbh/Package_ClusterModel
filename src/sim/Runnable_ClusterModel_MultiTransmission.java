package sim;

import java.util.ArrayList;
import java.util.Arrays;
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
			// double{{INF_ID, SITE_FROM, SITE_TO, TRANSMISSION_PARAM_0,
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
	public static final int FIELD_TRANSMISSION_RATE_SITE_FROM = FIELD_TRANSMISSION_RATE_INF_ID + 1;
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
	protected transient HashMap<Integer, double[][][]> map_trans_prob; // Key=PID, V=double[INF_ID][SITE_FROM][SITE_TO]
	protected transient HashMap<Integer, int[][]> map_currrent_infectState; // Key=PID,
																			// V=int[INF_ID][SITE]{infection_state}
	protected transient ArrayList<Integer>[][] list_currently_infectious; // ArrayList<Integer>[INF_ID][SITE]

	// For schedule, key are day and entries are list of person id, ordered
	protected transient HashMap<Integer, ArrayList<Integer>> schedule_testing;

	// Helper objects set by fields

	protected transient HashMap<String, double[]> map_infection_stage_path; // Key="INF_ID,SITE_ID,STAGE_ID",
																			// V={Cumul_Prob_0,... State_ID_0...}

	protected transient double[][][][] table_act_frequency; // double[ACT_ID][G_TO][G_FROM]{fieldEntry}

	protected transient RealDistribution[][][] dist_sym_rate; // RealDistribution[INF_ID][SITE][STAGE_ID]
	protected transient RealDistribution[][][] dist_tranmissionMatrix; // RealDistribution[INF_ID][SITE_FROM][SITE_TO]
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
		schedule_testing = new HashMap<>();

		table_act_frequency = new double[NUM_ACT][NUM_GENDER][NUM_GENDER][];
		map_infection_stage_path = new HashMap<>();

		dist_tranmissionMatrix = new RealDistribution[NUM_INF][NUM_SITE][NUM_SITE];
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
						for (double[] actType_gIf_gIt : actType_gIf) {
							Arrays.fill(actType_gIf_gIt, 0);
						}
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
				for (RealDistribution[][] inf : dist_tranmissionMatrix) {
					for (RealDistribution[] inf_sf : inf) {
						Arrays.fill(inf_sf, null);
					}
				}

			}
			for (double[] ent : field) {
				int inf_id = (int) ent[FIELD_TRANSMISSION_RATE_INF_ID];
				int site_from = (int) ent[FIELD_TRANSMISSION_RATE_SITE_FROM];
				int site_to = (int) ent[FIELD_TRANSMISSION_RATE_SITE_TO];
				double[] param = Arrays.copyOfRange(ent, FIELD_TRANSMISSION_RATE_TRANS_PARAM_START, ent.length);
				if (param.length == 1) {
					dist_tranmissionMatrix[inf_id][site_from][site_to] = generateNonDistribution(param);
				} else {
					// Beta distribution
					dist_tranmissionMatrix[inf_id][site_from][site_to] = generateBetaDistribution(param);
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
				if (param.length == 1) {
					dist_infectious_period[inf_id][site] = generateNonDistribution(param);
				} else {
					// Gamma distribution
					dist_infectious_period[inf_id][site] = generateGammaDistribution(param);
				}
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD:
			if (clearAll) {
				for (RealDistribution[][] inf : dist_stages_period) {
					for (RealDistribution[] inf_site : inf) {
						Arrays.fill(inf_site, null);
					}
				}
				map_infection_stage_path.clear();
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
				} else if (dist_stages_period[inf_id][site].length < stage_id_to) {
					dist_stages_period[inf_id][site] = Arrays.copyOf(dist_stages_period[inf_id][site], stage_id_to + 1);
				}
				if (param.length == 1) {
					dist_stages_period[inf_id][site][stage_id_to] = generateNonDistribution(param);
				} else {
					// Gamma distribution
					dist_stages_period[inf_id][site][stage_id_to] = generateGammaDistribution(param);
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
				map_infection_stage_path.put(key, cumul_prob);
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_SYM_RATE:
			if(clearAll) {
				for(RealDistribution[][] inf : dist_sym_rate) {
					for (RealDistribution[] inf_site : inf) {
						Arrays.fill(inf_site, null);
					}
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
		schedule_testing.clear();

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
		// num_infected = int[infection_id]{GENDER_INC_INDEX_0, SITE_INDEX_0, Number_INF_0,...}		
		
		
		// TODO Auto-generated method stub
		
	}

	@Override
	public int addInfectious(Integer infectedId, int infId, int site, int infectious_time, int recoveredAt) {
		
		
		// TODO Auto-generated method stub
		return 0;
	}

	
}
