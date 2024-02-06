package sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import relationship.ContactMap;

public class Runnable_ClusterModel_MultiTransmission extends Abstract_Runnable_ClusterModel_Transmission {

	public Object[] runnable_fields = {
			// 0: RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ
			// double{{ACT_TYPE, GENDER_INCLUDE_INDEX_FROM, GENDER_INCLUDE_INDEX_TO},...}
			// or
			// double{{ACT_TYPE, GENDER_INCLUDE_INDEX_FROM, GENDER_INCLUDE_INDEX_TO,
			// CONDOM_EFFICACY, USAGE_REG, USAGE_CASUAL},...}
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
			// DURATION_PARAM_0 ...},...}
			new double[][] {},
			// 4: RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			// float{{INF_ID, SITE, STAGE_ID, RATE},...}
			new double[][] {},
			// 5: RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			// double{{GENDER_INCLUDE_INDEX, RISK_GRP_DEF_ID, NUM_CASUAL_PARTNER_LOWER,
			// NUM_CASUAL_PARTNER_UPPER}, ... }
			new double[][] {},
			// 6: RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES
			// double{{GENDER_INCLUDE_INDEX, INF_INCLUDE_INDEX, SITE_INCLUDE_INDEX,
			// RISK_GRP_TEST_ID, TESTING_RATE_PARAM...}, ...}
			new double[][] {},
			// 7: RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM;
			// double{{GENDER_INCLUDE_INDEX, SITE, DAILY_SOUGHT_TEST_RATE}, ...}
			new double[][] {},
			// 8 RUNNALBE_FIELD_TRANSMISSION_DX_TEST_ACCURACY
			// double{{GENDER_INCLUDE_INDEX, INF_INCLUDE_INDEX, SITE, TEST_ACCURACY}, ...}
			new double[][] {}, };

	protected transient HashMap<Integer, double[][][]> map_trans_prob; // Key = PID, V =
																		// double[INF_ID][SITE_FROM][SITE_TO]
	protected transient HashMap<Integer, int[][]> map_currrent_infectState; // Key = PID, V =
																			// int[INF_ID][SITE]{infection_state}

	protected transient ArrayList<Integer>[][] arr_currently_infectious; // ArrayList<Integer>[INF_ID][SITE]

	// Fixed input value
	protected final int NUM_INF;
	protected final int NUM_ACT;
	protected final int NUM_GENDER;

	// For RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ
	public static final int FIELD_INF_ID = 0;
	public static final int FIELD_SITE = FIELD_INF_ID + 1;
	public static final int FIELD_GENDER_INCLUDE_INDEX = 0;
	public static final int FIELD_INF_INCLUDE_INDEX = FIELD_GENDER_INCLUDE_INDEX + 1;
	public static final int FIELD_SITE_INCLUDE_INDEX = FIELD_INF_INCLUDE_INDEX + 1;

	// For RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ only
	public static final int FIELD_ACT_TYPE = 0;
	public static final int FIELD_GENDER_INCLUDE_INDEX_FROM = FIELD_ACT_TYPE + 1;
	public static final int FIELD_GENDER_INCLUDE_INDEX_TO = FIELD_GENDER_INCLUDE_INDEX_FROM + 1;
	public static final int FIELD_CONDOM_EFFICACY = FIELD_GENDER_INCLUDE_INDEX_TO + 1;
	public static final int FIELD_USAGE_REG = FIELD_CONDOM_EFFICACY + 1;
	public static final int FIELD_USAGE_CASUAL = FIELD_USAGE_REG + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD only
	// For RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD only
	public static final int FIELD_STATE_FROM = FIELD_SITE + 1;
	public static final int FIELD_STATE_TO = FIELD_STATE_FROM + 1;
	public static final int FIELD_SWTICH_STAGE_PROB = FIELD_STATE_TO + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_SYM_RATE only
	public static final int FIELD_STATE_ID = FIELD_SITE + 1;
	public static final int FIELD_SYM_RATE = FIELD_STATE_ID + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS only
	public static final int FIELD_RISK_GRP_DEF_ID = FIELD_GENDER_INCLUDE_INDEX + 1;
	public static final int FIELD_NUM_CASUAL_PARTNER_LOWER = FIELD_RISK_GRP_DEF_ID + 1;
	public static final int FIELD_NUM_CASUAL_PARTNER_UPPER = FIELD_NUM_CASUAL_PARTNER_LOWER + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES only
	public static final int FIELD_RISK_GRP_TEST_ID = FIELD_SITE_INCLUDE_INDEX + 1;
	// For RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM only
	public static final int FIELD_DAILY_SOUGHT_TEST_RATE = FIELD_SITE + 1;
	// For RUNNALBE_FIELD_TRANSMISSION_DX_TEST_ACCURACY only
	public static final int FIELD_TEST_ACCURACY = FIELD_SITE + 1;

	public Runnable_ClusterModel_MultiTransmission(long cMap_seed, long sim_seed, int[] pop_composition,
			ContactMap base_cMap, int numTimeStepsPerSnap, int numSnap, int num_inf, int num_act) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap);
		NUM_INF = num_inf;
		NUM_ACT = num_act;
		NUM_GENDER = pop_composition.length;

		// Initiate transient field, lookup table etc
		map_trans_prob = new HashMap<>();
		map_currrent_infectState = new HashMap<>();

	}

	public void initialse() {
		for (int r = 0; r < runnable_fields.length; r++) {
			refreshField(r);
		}
	}

	public void refreshField(int fieldId) {
		double[][] field = null;

		try {
			field = (double[][]) getRunnable_fields()[fieldId];
		} catch (ClassCastException ex) {
			field = null;
		}

		// Reset runnable field
		// TODO: to be implemented
		switch (fieldId) {
		case RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ:

			break;
		default:
			System.err.printf("Warning: refreshField(fieldId=%d) not defined.\n", fieldId);
			if (field != null) {
				System.err.println("Existing field:");
				for(double[] ent : field) {
					System.err.println(Arrays.toString(ent));
				}
			}

		}

	}

	@Override
	public void run() {
		// TODO Auto-generated method stub

	}

	@Override
	public Object[] getRunnable_fields() {
		return runnable_fields;
	}

}
