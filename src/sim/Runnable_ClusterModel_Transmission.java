package sim;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import population.person.Person_Bridging_Pop;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_Transmission extends Abstract_Runnable_ClusterModel {

	public static final int ACT_INDEX_GENITAL = 0;
	public static final int ACT_INDEX_ANAL = ACT_INDEX_GENITAL + 1;
	public static final int ACT_INDEX_FELLATIO = ACT_INDEX_ANAL + 1;
	public static final int ACT_INDEX_RIMMING = ACT_INDEX_FELLATIO + 1;
	public static final int ACT_INDEX_KISSING = ACT_INDEX_RIMMING + 1;

	public static final int SITE_VAGINA = 0;
	public static final int SITE_PENIS = SITE_VAGINA + 1;
	public static final int SITE_RECTUM = SITE_PENIS + 1;
	public static final int SITE_OROPHARYNX = SITE_RECTUM + 1;
	public static final int LENGTH_SITE = SITE_OROPHARYNX + 1;

	protected static final int TEST_OUTCOME_TREATMENT_APPLIED = 0;
	protected static final int TEST_OUTCOME_TREATMENT_APPLIED_ON_TRUE_INFECTION = TEST_OUTCOME_TREATMENT_APPLIED + 1;
	protected static final int TEST_OUTCOME_POSITIVE_TEST = TEST_OUTCOME_TREATMENT_APPLIED_ON_TRUE_INFECTION + 1;
	protected static final int TEST_OUTCOME_POSITIVE_TEST_ON_TRUE_INFECTION = TEST_OUTCOME_POSITIVE_TEST + 1;

	protected int simSetting = 1;

	// Transmission
	private double[] DEFAULT_TRANS_V2P = new double[] { 0.4, 0.10 };
	private double[] DEFAULT_TRANS_P2V = new double[] { 0.2, 0.05 };
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private double[] DEFAULT_TRANS_P2R = new double[] { 0.63, 0 };
	private double[] DEFAULT_TRANS_R2P = new double[] { 0.02, 0 };
	private double[] DEFAULT_TRANS_P2O = new double[] { 0.44, 0 };
	private double[] DEFAULT_TRANS_O2P = new double[] { 0.01, 0 };

	// Duration
	private double[] DEFAULT_INFECTIOUS_PERIOD_VAGINA = new double[] { 15 * 7, 5 * 7 };
	private double[] DEFAULT_INFECTIOUS_PERIOD_PENIS = new double[] { 15 * 7, 5 * 7 };
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private double[] DEFAULT_INFECTIOUS_PERIOD_RECTUM = new double[] { 307.2, 6.08 };
	private double[] DEFAULT_INFECTIOUS_PERIOD_OROPHARYNX = new double[] { 80.4, 5.67 };

	// Incubation
	private double[] DEFAULT_INCUBATION_RANGE = new double[] { 3, 6 }; // 3 - 5 days

	// Sought test period
	// From Qibin's paper of Gamma(3, 0.86)
	private double[] DEFAULT_SYM_TEST = new double[] { 3 * 0.86, Math.sqrt(3 * 0.86 * 0.86) };

	// Act frequency
	// ASHR2: Those who were in heterosexual relationships had had
	// sex on average 1.84 times a week in the past 4 weeks,
	private float DEFAULT_ACT_GENITAL_FREQ = 1.84f / 7;
	// ASHR2 : only 1% of men and 0.4% of women had anal intercourse
	// the last time they had sex with a partner of the other sex
	private float DEFAULT_ACT_ANAL_FREQ_HETRO = 0.1f / 100;
	// ASHR2: oral sex was reported in only approximately one in four encounters
	private float DEFAULT_ACT_FELLATIO_FREQ_HETRO = 0.25f;
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private float DEFAULT_ACT_ANAL_FREQ_MSM = ((1.6f + 2.4f)) / 2 / 7;
	private float DEFAULT_ACT_FELLATIO_FREQ_MSM = ((1.6f + 2.4f) / 2) / 7;
	private float DEFAULT_ACT_RIMMING_FREQ_MSM = ((1.2f + 1.8f) / 2) / 7;
	private float DEFAULT_ACT_KISSING_FREQ_MSM = ((2.4f + 3.6f) / 2) / 7;

	private float[] DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM = new float[] { 20 };
	private float[] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_LOW_RISK = new float[] { 0.375f, 1, 720, 360, 90 };
	private float[] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_HIGH_RISK = new float[] { 0.05f, 0.49f, 0.71f, 0.99f, 720,
			360, 180, 120, 90 };
	private float[][] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM = new float[][] {
			DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_LOW_RISK, DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_HIGH_RISK };

	private float[] DEFAULT_SYM_RATE_BY_SITE_MALE = new float[] { Float.NaN, 1.0f, 0.12f, 0 };
	private float[] DEFAULT_SYM_RATE_BY_SITE_FEMALE = new float[] { 0.4f, Float.NaN, 0.12f, 0 };

	// DX and viability test
	// float[dx_sensitivity, dx_specificity][gender][site]
	private float[][][] DEFAULT_TEST_ACCURACY = new float[][][] {
			// Dx sensitivity - i.e. +ive identified correctly
			new float[][] { new float[] { 1f, 1f, 1f, 1f }, new float[] { 1f, 1f, 1f, 1f },
					new float[] { 1f, 1f, 1f, 1f }, new float[] { 1f, 1f, 1f, 1f } },
			// Dx specificity - i.e. -ive identified correctly
			new float[][] { new float[] { 1f, 1f, 1f, 1f }, new float[] { 1f, 1f, 1f, 1f },
					new float[] { 1f, 1f, 1f, 1f }, new float[] { 1f, 1f, 1f, 1f } },

	};

	public final static int DX_SENSITIVITY_INDEX = 0;
	public final static int DX_SPECIFICITY_INDEX = DX_SENSITIVITY_INDEX + 1;
	public final static int DX_TREATMENT_RATE_INDEX = DX_SPECIFICITY_INDEX + 1;

	// float[gender][site][setting]
	private float[][][] DEFAULT_NON_VIABLE_INFECTION_SETTING = new float[Population_Bridging.LENGTH_GENDER][LENGTH_SITE][LENGTH_NON_VIABILITY_SETTING];

	public static final int NON_VIABILITY_CONTACT_INDUCED_PROB = 0;
	public static final int NON_VIABILITY_CONTACT_INDUCED_DURATION_MEAN = NON_VIABILITY_CONTACT_INDUCED_PROB + 1;
	public static final int NON_VIABILITY_CONTACT_INDUCED_DURATION_SD = NON_VIABILITY_CONTACT_INDUCED_DURATION_MEAN + 1;
	public static final int NON_VIABILITY_TREATMENT_INDUCED_PROB = NON_VIABILITY_CONTACT_INDUCED_DURATION_SD + 1;
	public static final int NON_VIABILITY_TREATMENT_INDUCED_INFECTION_REDUCTION_ADJ = NON_VIABILITY_TREATMENT_INDUCED_PROB
			+ 1;
	public static final int NON_VIABILITY_RECOVERY_INDUDCED_PROB = NON_VIABILITY_TREATMENT_INDUCED_INFECTION_REDUCTION_ADJ
			+ 1;
	public static final int NON_VIABILITY_RECOVERY_INDUDCED_MEAN = NON_VIABILITY_RECOVERY_INDUDCED_PROB + 1;
	public static final int NON_VIABILITY_RECOVERY_INDUDCED_SD = NON_VIABILITY_RECOVERY_INDUDCED_MEAN + 1;
	public static final int LENGTH_NON_VIABILITY_SETTING = NON_VIABILITY_RECOVERY_INDUDCED_SD + 1;

	// double[gender][mean, sd]
	private double[][] DEFAULT_ANTIBIOTIC_DURATION = new double[][] { new double[] { 7, 7 }, new double[] { 7, 7 },
			new double[] { 7, 7 }, new double[] { 7, 7 }, };

	// float[gender][site][effect]
	private float[][][] DEFAULT_VACCINE_PROPERTIES = new float[Population_Bridging.LENGTH_GENDER][LENGTH_SITE][LENGTH_VACCINE_PROPERTIES];

	public static final int VACCINE_PROPERTIES_MEAN_VAC_DURATION = 0;
	public static final int VACCINE_PROPERTIES_INF_DURATION_ADJUST = VACCINE_PROPERTIES_MEAN_VAC_DURATION + 1;
	public static final int VACCINE_PROPERTIES_SUSCEPTIBILITY_ADJUST = VACCINE_PROPERTIES_INF_DURATION_ADJUST + 1;
	public static final int VACCINE_PROPERTIES_TRANSMISSION_ADJ = VACCINE_PROPERTIES_SUSCEPTIBILITY_ADJUST + 1;
	public static final int VACCINE_PROPERTIES_SYMPTOM_ADJ = VACCINE_PROPERTIES_TRANSMISSION_ADJ + 1;
	public static final int LENGTH_VACCINE_PROPERTIES = VACCINE_PROPERTIES_SYMPTOM_ADJ + 1;

	private float[][] DEFAULT_VACCINATION_SETTING = new float[Population_Bridging.LENGTH_GENDER][LENGTH_VACCINATION_SETTING];

	public static final int VACCINATION_SETTING_BOOSTER_PERIOD = 0;
	public static final int VACCINATION_SETTING_BOOSTER_LIMIT = VACCINATION_SETTING_BOOSTER_PERIOD + 1;
	public static final int VACCINATION_SETTING_RATE_PER_TEST = VACCINATION_SETTING_BOOSTER_LIMIT + 1;
	public static final int LENGTH_VACCINATION_SETTING = VACCINATION_SETTING_RATE_PER_TEST + 1;

	public static final int RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ = 0;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE = RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ + 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD = RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_INCUBATION_PERIOD = RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SYM_RATE = RUNNABLE_FIELD_TRANSMISSION_INCUBATION_PERIOD + 1;

	public static final int RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS = RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES = RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM = RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_DX_TEST_ACCURACY = RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_NON_VIABLE_INFECTION_SETTING = RUNNABLE_FIELD_TRANSMISSION_DX_TEST_ACCURACY
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_ANTIBIOTIC_DURATION = RUNNABLE_FIELD_TRANSMISSION_NON_VIABLE_INFECTION_SETTING
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_VACCINE_PROPERTIES = RUNNABLE_FIELD_TRANSMISSION_ANTIBIOTIC_DURATION
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_VACCINE_SETTING = RUNNABLE_FIELD_TRANSMISSION_VACCINE_PROPERTIES
			+ 1;
	public static final int LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD = RUNNABLE_FIELD_TRANSMISSION_VACCINE_SETTING + 1;

	public Object[] runnable_fields = {
			// RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ
			// double[ACT_TYPE][GENDER_FROM][GENDER_TO]
			new float[][][] {
					// ACT_INDEX_GENITAL
					new float[][] { new float[] { 0, DEFAULT_ACT_GENITAL_FREQ, 0, DEFAULT_ACT_GENITAL_FREQ },
							new float[] { DEFAULT_ACT_GENITAL_FREQ, 0, 0, 0 }, new float[] { 0, 0, 0, 0 },
							new float[] { DEFAULT_ACT_GENITAL_FREQ, 0, 0, 0 }, },
					// ACT_INDEX_ANAL
					new float[][] { new float[] { 0, DEFAULT_ACT_ANAL_FREQ_HETRO, 0, DEFAULT_ACT_ANAL_FREQ_HETRO },
							new float[] { DEFAULT_ACT_ANAL_FREQ_HETRO, DEFAULT_ACT_ANAL_FREQ_HETRO, 0, 0 },
							new float[] { 0, 0, DEFAULT_ACT_ANAL_FREQ_MSM, DEFAULT_ACT_ANAL_FREQ_MSM },
							new float[] { DEFAULT_ACT_ANAL_FREQ_HETRO, 0, DEFAULT_ACT_ANAL_FREQ_MSM,
									DEFAULT_ACT_ANAL_FREQ_MSM }, },
					// ACT_INDEX_FELLATIO
					new float[][] {
							new float[] { 0, DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, DEFAULT_ACT_FELLATIO_FREQ_HETRO },
							new float[] { DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, 0, 0 },
							new float[] { 0, 0, DEFAULT_ACT_FELLATIO_FREQ_MSM, DEFAULT_ACT_FELLATIO_FREQ_MSM },
							new float[] { DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, DEFAULT_ACT_FELLATIO_FREQ_MSM,
									DEFAULT_ACT_FELLATIO_FREQ_MSM }, },
					// ACT_INDEX_RIMMING (OROPHARYNX <-> RECTUM)
					new float[][] { new float[] { 0, 0, 0, 0 }, new float[] { 0, 0, 0, 0 },
							new float[] { 0, 0, DEFAULT_ACT_RIMMING_FREQ_MSM, DEFAULT_ACT_RIMMING_FREQ_MSM },
							new float[] { 0, 0, DEFAULT_ACT_RIMMING_FREQ_MSM, DEFAULT_ACT_RIMMING_FREQ_MSM }, },

					// ACT_INDEX_KISSING (OROPHARYNX <-> OROPHARYNX)
					new float[][] { new float[] { 0, 0, 0, 0 }, new float[] { 0, 0, 0, 0 },
							new float[] { 0, 0, DEFAULT_ACT_KISSING_FREQ_MSM, DEFAULT_ACT_KISSING_FREQ_MSM },
							new float[] { 0, 0, DEFAULT_ACT_KISSING_FREQ_MSM, DEFAULT_ACT_KISSING_FREQ_MSM }, },

			},
			// RUNNABLE_FIELD_TRANSMISSION_MAP_TRANSMISSION_RATE
			// double[SITE_FROM][SITE_TO][] =
			new double[][][] {
					// FROM SITE_VAGINA
					new double[][] { null, DEFAULT_TRANS_V2P, null, null },
					// FROM SITE_PENIS
					new double[][] { DEFAULT_TRANS_P2V, null, DEFAULT_TRANS_P2R, DEFAULT_TRANS_P2O },
					// FROM SITE_RECTUM
					new double[][] { null, DEFAULT_TRANS_R2P, null, null },
					// FROM SITE_OROPHARYNX
					new double[][] { null, DEFAULT_TRANS_O2P, null, null }, },

			// RUNNABLE_FIELD_INFECTIOUS_PERIOD
			// double[SITE]
			new double[][] { DEFAULT_INFECTIOUS_PERIOD_VAGINA, DEFAULT_INFECTIOUS_PERIOD_PENIS,
					DEFAULT_INFECTIOUS_PERIOD_RECTUM, DEFAULT_INFECTIOUS_PERIOD_OROPHARYNX, },
			// RUNNABLE_FIELD_INCUBATION_PERIOD
			// double[SITE]
			new double[][] { DEFAULT_INCUBATION_RANGE, DEFAULT_INCUBATION_RANGE, DEFAULT_INCUBATION_RANGE,
					DEFAULT_INCUBATION_RANGE },
			// RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			new float[][] { DEFAULT_SYM_RATE_BY_SITE_FEMALE, DEFAULT_SYM_RATE_BY_SITE_MALE,
					DEFAULT_SYM_RATE_BY_SITE_MALE, DEFAULT_SYM_RATE_BY_SITE_MALE },
			// RUNNABLE_FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			new float[][] { null, null, DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM,
					DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM },
			// RUNNABLE_FIELD_TESTING_RATE_BY_RISK_CATEGORIES
			new float[][][] { null, null, DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM,
					DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM },
			// RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
			// Format: double[] {mean, sd}
			// Alt Format: double[] { mean_gender_0, sd_gender_0, ....}
			DEFAULT_SYM_TEST,
			// RUNNABLE_FIELD_TRANSMISSION_TEST_ACCURACY
			// Format: float[dx_sensitivity, dx_specificity][gender][site]
			// or float[dx_sensitivity, dx_specificity, treatment_rate][gender][site]

			DEFAULT_TEST_ACCURACY,
			// RUNNABLE_FIELD_TREATMENT_INDUCED_NON_VIABILITY
			DEFAULT_NON_VIABLE_INFECTION_SETTING,
			// RUNNABLE_FIELD_TRANSMISSION_ANTIBIOTIC_DURATION
			DEFAULT_ANTIBIOTIC_DURATION,
			// RUNNABLE_FIELD_TRANSMISSION_VACCINE_EFFECT
			DEFAULT_VACCINE_PROPERTIES,
			// RUNNABLE_FIELD_TRANSMISSION_VACCINE_SETTING
			DEFAULT_VACCINATION_SETTING, };

	// FIELD_POP_COMPOSITION
	// int[] {NUM_FEMALE, NUM_MALE, NUM_MSMO, NUM_MSMW}
	final int[] cumulative_pop_composition;
	final ContactMap BASE_CONTACT_MAP;
	final int NUM_TIME_STEPS_PER_SNAP;
	final int NUM_SNAP;
	protected final long sim_seed;
	protected final long cMap_seed;

	protected RandomGenerator RNG;

	protected transient RealDistribution[][] tranmissionMatrix = new RealDistribution[LENGTH_SITE][LENGTH_SITE];
	protected transient RealDistribution[] infectious_period = new RealDistribution[LENGTH_SITE];
	protected transient RealDistribution[] incubation_period = new RealDistribution[LENGTH_SITE];
	protected transient RealDistribution[] sym_test_period_by_gender = new RealDistribution[Population_Bridging.LENGTH_GENDER];

	protected transient ArrayList<Integer>[] currently_infectious;

	// For schedule, key are day and entries are list of person id, ordered
	protected transient HashMap<Integer, ArrayList<Integer>>[] schedule_incubation;
	protected transient HashMap<Integer, ArrayList<Integer>>[] schedule_recovery;
	protected transient HashMap<Integer, ArrayList<Integer>> schedule_testing;
	protected transient HashMap<Integer, Integer[]> mapping_infection_schedule;

	protected transient HashMap<Integer, double[][]> trans_prob;

	protected transient HashMap<Integer, Integer> risk_cat_map;
	protected transient int firstSeedTime = Integer.MAX_VALUE;
	protected transient HashMap<String, Object> sim_output = null;

	protected transient HashMap<Integer, Integer> has_non_viable_bacteria_until;
	// Key = time, V = Map(K = run_parameter_id, V = PropVal_Str
	protected HashMap<Integer, HashMap<Integer, String>> propSwitch_map;

	// For infection tracking
	// Key = pid, V = [site][infection_start_time_1,
	// infection_end_time_1...];
	protected transient HashMap<Integer, ArrayList<Integer>[]> infection_history;

	// For antibiotic tracking
	protected transient RealDistribution[][] non_viable_inf_by_contact_duration = new RealDistribution[Population_Bridging.LENGTH_GENDER][LENGTH_SITE];
	protected transient RealDistribution[][] non_viable_inf_by_recovery_duration = new RealDistribution[Population_Bridging.LENGTH_GENDER][LENGTH_SITE];
	protected transient RealDistribution[] antibotic_duration = new RealDistribution[Population_Bridging.LENGTH_GENDER];
	// if pid < 0, it is over-treatment
	protected transient ArrayList<Integer> currently_has_antibiotic;
	// Key = date, V = pid (or -pid if overtreatment)
	protected transient HashMap<Integer, ArrayList<Integer>> schedule_antibiotic_clearance;

	// For vaccine tracking
	protected transient ArrayList<Integer> currently_vaccinated;
	// Key = pid, V = int{expiry_at_site_0, ... }>
	protected transient HashMap<Integer, int[]> vaccine_expiry_by_indivdual;
	// Key = day, V = pid of those who receive booster shot
	protected transient HashMap<Integer, ArrayList<Integer>> schedule_vaccination;

	protected transient RealDistribution[][] vacc_dur = new RealDistribution[Population_Bridging.LENGTH_GENDER][LENGTH_SITE];

	// HashMap<Integer, int[][]> with K = time, V= int[gender][site]
	public static final String SIM_OUTPUT_INFECTIOUS_COUNT = "SIM_OUTPUT_INFECTIOUS_COUNT";
	// HashMap<Integer, int[][]> with K = time, V= int[gender][site]
	public static final String SIM_OUTPUT_CUMUL_INCIDENCE = "SIM_OUTPUT_CUMUL_INCIDENCE";

	// HashMap<Integer, int[]> with K = time, V= int[gender]
	public static final String SIM_OUTPUT_INFECTIOUS_COUNT_BY_PERSON = "SIM_OUTPUT_INFECTIOUS_COUNT_BY_PERSON";
	// HashMap<Integer, int[]> with K = time, V= int[gender]
	public static final String SIM_OUTPUT_CUMUL_INCIDENCE_BY_PERSON = "SIM_OUTPUT_CUMUL_INCIDENCE_BY_PERSON";
	// HashMap<Integer, int[]> with K = time, V= int[gender_apply_treatment,
	// gender_true_infection]
	public static final String SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON = "SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON";
	// HashMap<Integer, int[]> with K = time, V= int[gender_positive_dx,
	// gender_true_infection]
	public static final String SIM_OUTPUT_CUMUL_POS_DX_BY_PERSON = "SIM_OUTPUT_CUMUL_POS_DX_BY_PERSON";
	// HashMap<Integer, int[]> with K = time, V= int[gender_positive_dx,
	// gender_true_infection]
	public static final String SIM_OUTPUT_CUMUL_POS_DX_SOUGHT_BY_PERSON = "SIM_OUTPUT_CUMUL_POS_DX_SOUGHT_BY_PERSON";
	

	// HashMap<Integer, int[][]>
	// with K = time, V= int[gender]{proper,over treatment} measured in person-day
	public static final String SIM_OUTPUT_CUMUL_ANTIBOTIC_USAGE = "SIM_OUTPUT_CUMUL_ANTIBOTIC_USAGE";
	// Key = pid, V = [site][infection_start_time_1,
	// infection_end_time_1...];
	public static final String SIM_OUTPUT_INFECTION_HISTORY = "SIM_OUTPUT_INFECTION_HISTORY";

	// HashMap<Integer, int[][][]>
	// with K = time, V= int[gender][site]{valid, expired}
	public static final String SIM_OUTPUT_VACCINE_COVERAGE = "SIM_OUTPUT_VACCINE_COVERAGE";

	private ArrayList<Integer[]> edges_list;
	private static final int RUNNABLE_OFFSET = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
			+ +Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
			+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD;

	public Runnable_ClusterModel_Transmission(long cMap_seed, long sim_seed, int[] POP_COMPOSITION,
			ContactMap BASE_CONTACT_MAP, int NUM_TIME_STEPS_PER_SNAP, int NUM_SNAP) {
		super();

		this.cumulative_pop_composition = new int[POP_COMPOSITION.length];
		int offset = 0;

		for (int g = 0; g < this.cumulative_pop_composition.length; g++) {
			this.cumulative_pop_composition[g] = offset + POP_COMPOSITION[g];
			offset += POP_COMPOSITION[g];
		}

		this.BASE_CONTACT_MAP = BASE_CONTACT_MAP;
		this.NUM_TIME_STEPS_PER_SNAP = NUM_TIME_STEPS_PER_SNAP;
		this.NUM_SNAP = NUM_SNAP;
		this.cMap_seed = cMap_seed;
		this.sim_seed = sim_seed;

	}

	public void setPropSwitch_map(HashMap<Integer, HashMap<Integer, String>> propSwitch_map) {
		this.propSwitch_map = propSwitch_map;
	}

	public int getSimSetting() {
		return simSetting;
	}

	public void setSimSetting(int simSetting) {
		this.simSetting = simSetting;
	}

	public void setEdges_list(ArrayList<Integer[]> edges_list) {
		this.edges_list = edges_list;
	}

	@Override
	public Object[] getRunnable_fields() {
		return runnable_fields;
	}

	@SuppressWarnings("unchecked")
	public void initialse() {
		RNG = new MersenneTwisterRandomGenerator(sim_seed);

		// Transmission
		double[][][] tranParm = (double[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE];
		for (int sf = 0; sf < tranmissionMatrix.length; sf++) {
			for (int st = 0; st < tranmissionMatrix[sf].length; st++) {
				double[] param = tranParm[sf][st];
				if (param != null) {
					tranmissionMatrix[sf][st] = generateBetaDistribution(param);
				}
			}
		}

		// Duration & incubation
		double[][] durParam = (double[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD];
		double[][] incParam = (double[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_INCUBATION_PERIOD];
		for (int s = 0; s < infectious_period.length; s++) {
			infectious_period[s] = generateGammaDistribution(durParam[s]);
			incubation_period[s] = new UniformRealDistribution(RNG, incParam[s][0], incParam[s][1]);
		}

		sym_test_period_by_gender = new RealDistribution[Population_Bridging.LENGTH_GENDER];
		double[] sought_test_param = (double[]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM];

		if (sought_test_param.length == 2) {
			RealDistribution dist = generateGammaDistribution(sought_test_param);
			for (int g = 0; g < sym_test_period_by_gender.length; g++) {
				sym_test_period_by_gender[g] = dist;
			}
		} else {
			int offset = 0;
			for (int g = 0; g < sym_test_period_by_gender.length; g++) {
				sym_test_period_by_gender[g] = generateGammaDistribution(
						Arrays.copyOfRange(sought_test_param, offset, offset + 2));
				offset += 2;
			}
		}

		// Lists
		currently_infectious = new ArrayList[LENGTH_SITE];
		schedule_incubation = new HashMap[LENGTH_SITE];
		schedule_recovery = new HashMap[LENGTH_SITE];

		for (int s = 0; s < LENGTH_SITE; s++) {
			currently_infectious[s] = new ArrayList<>();
			schedule_incubation[s] = new HashMap<>();
			schedule_recovery[s] = new HashMap<>();

		}

		trans_prob = new HashMap<>();
		sim_output = new HashMap<>();
		schedule_testing = new HashMap<>();
		mapping_infection_schedule = new HashMap<>();
		risk_cat_map = new HashMap<>();

		// For infection tracking
		// Key = pid, V = [infection_start_time_1, infection_end_time_1...];
		infection_history = new HashMap<>();

		// Antibiotic tracking
		currently_has_antibiotic = new ArrayList<>();
		schedule_antibiotic_clearance = new HashMap<>();
		has_non_viable_bacteria_until = new HashMap<>();

		float[][][] non_viable_infection_setting = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_NON_VIABLE_INFECTION_SETTING]);
		double[][] antibotic_dur = (double[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_ANTIBIOTIC_DURATION];

		for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
			for (int s = 0; s < LENGTH_SITE; s++) {
				if (non_viable_infection_setting[g][s][NON_VIABILITY_CONTACT_INDUCED_PROB] > 0) {
					non_viable_inf_by_contact_duration[g][s] = generateGammaDistribution(new double[] {
							non_viable_infection_setting[g][s][NON_VIABILITY_CONTACT_INDUCED_DURATION_MEAN],
							non_viable_infection_setting[g][s][NON_VIABILITY_CONTACT_INDUCED_DURATION_SD] });
				}
				if (non_viable_infection_setting[g][s][NON_VIABILITY_RECOVERY_INDUDCED_PROB] > 0) {
					non_viable_inf_by_recovery_duration[g][s] = generateGammaDistribution(
							new double[] { non_viable_infection_setting[g][s][NON_VIABILITY_RECOVERY_INDUDCED_MEAN],
									non_viable_infection_setting[g][s][NON_VIABILITY_RECOVERY_INDUDCED_SD] });

				}
			}

			antibotic_duration[g] = generateGammaDistribution(antibotic_dur[g]);
		}

		// Runnable properties switch
		propSwitch_map = new HashMap<>();

		// Vaccines
		currently_vaccinated = new ArrayList<>();
		vaccine_expiry_by_indivdual = new HashMap<>();

		for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
			for (int s = 0; s < LENGTH_SITE; s++) {
				float[] vacc_setting = ((float[][][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_VACCINE_PROPERTIES])[g][s];
				if (vacc_setting[VACCINE_PROPERTIES_MEAN_VAC_DURATION] > 0) {
					vacc_dur[g][s] = new ExponentialDistribution(vacc_setting[VACCINE_PROPERTIES_MEAN_VAC_DURATION]);
				} else {
					vacc_dur[g][s] = null;
				}
			}
		}
		schedule_vaccination = new HashMap<>();

	}

	public HashMap<String, Object> getSim_output() {
		return sim_output;
	}

	public void scheduleNextTest(Integer personId, int lastTestTime) {
		int genderType = getGenderType(personId);

		float[][] testRate = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES])[genderType];
		if (testRate != null) {
			int riskCat = Math.max(0, getRiskCategories(personId, genderType));
			float[] testRateByCat = testRate[riskCat];
			int divder = (testRateByCat.length - 1) / 2;

			int pI;

			float p = RNG.nextFloat();
			pI = Arrays.binarySearch(testRateByCat, 0, divder, p);
			if (pI < 0) {
				pI = ~pI;
			}

			if (pI < divder) {

				double testGapTime = (testRateByCat[divder + pI] + testRateByCat[divder + pI + 1]) / 2;
				testGapTime *= 1 + RNG.nextGaussian() / 10;

				Integer testTime = (int) Math.round(lastTestTime + (testGapTime));

				ArrayList<Integer> testEnt = schedule_testing.get(testTime);
				if (testEnt == null) {
					testEnt = new ArrayList<>();
					schedule_testing.put(testTime, testEnt);
				}
				int key = Collections.binarySearch(testEnt, personId);
				if (key < 0) {
					testEnt.add(~key, personId);
				}
			}

		}

	}

	private int getRiskCategories(Integer personId, int genderType) {

		if (risk_cat_map.containsKey(personId)) {
			return risk_cat_map.get(personId);
		} else {
			int riskCat = -1;
			float[] riskCatList = ((float[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS])[genderType];
			if (riskCatList != null) {
				int numCasual = 0;
				Set<Integer[]> edges = BASE_CONTACT_MAP.edgesOf(personId);
				for (Integer[] e : edges) {
					if (e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] >= firstSeedTime
							&& e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] < firstSeedTime
									+ NUM_TIME_STEPS_PER_SNAP
							&& e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] <= 1) {
						numCasual++;
					}
				}

				float numCasual1Year = ((float) AbstractIndividualInterface.ONE_YEAR_INT) * numCasual
						/ NUM_TIME_STEPS_PER_SNAP;
				riskCat = Arrays.binarySearch(riskCatList, numCasual1Year);
				if (riskCat < 0) {
					riskCat = ~riskCat;
				}

			}

			risk_cat_map.put(personId, riskCat);

			return riskCat;
		}

	}

	public void allocateSeedInfection(int[][] num_infectioned_by_gender_site, int time) {

		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] pid_collection = new ArrayList[Population_Bridging.LENGTH_GENDER];
		for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
			pid_collection[g] = new ArrayList<>();
		}

		for (Integer v : BASE_CONTACT_MAP.vertexSet()) {
			pid_collection[getGenderType(v)].add(v);
		}

		for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
			Integer[] candidate = pid_collection[g].toArray(new Integer[pid_collection[g].size()]);
			for (int s = 0; s < LENGTH_SITE; s++) {
				int numInfect = num_infectioned_by_gender_site[g][s];
				for (int p = 0; p < candidate.length && numInfect > 0; p++) {
					if (RNG.nextInt(candidate.length - p) < numInfect) {
						addInfectious(candidate[p], s, time, time + (int) Math.round(infectious_period[s].sample()));
						numInfect--;
					}
				}
			}
		}
	}

	@SuppressWarnings("unchecked")
	public int addInfectious(Integer infectedId, int site, int infectious_time, int recoveredAt) {
		int key = Collections.binarySearch(currently_infectious[site], infectedId);

		if (key < 0) {
			currently_infectious[site].add(~key, infectedId);
			firstSeedTime = Math.min(firstSeedTime, infectious_time);

			// Recovery
			ArrayList<Integer> sch = schedule_recovery[site].get(recoveredAt);
			if (sch == null) {
				sch = new ArrayList<>();
				schedule_recovery[site].put(recoveredAt, sch);
			}
			int keyR = Collections.binarySearch(sch, infectedId);
			if (keyR < 0) {
				sch.add(~keyR, infectedId);
				updateScheduleMap(infectedId, site, null);
				updateScheduleMap(infectedId, LENGTH_SITE + site, recoveredAt);
			}

			int gender = getGenderType(infectedId);

			// Transmission probability
			if (!trans_prob.containsKey(infectedId)) {
				double[][] trans = new double[LENGTH_SITE][LENGTH_SITE];
				for (int sf = 0; sf < LENGTH_SITE; sf++) {
					boolean sample = (gender == Person_Bridging_Pop.GENDER_TYPE_FEMALE) ? sf != SITE_PENIS
							: sf != SITE_VAGINA;
					if (sample) {
						for (int st = 0; st < LENGTH_SITE; st++) {
							if (tranmissionMatrix[sf][st] != null) {
								trans[sf][st] = tranmissionMatrix[sf][st].sample();
							}
						}
					}
				}
				trans_prob.put(infectedId, trans);
			}

			// Schedule sym test if needed.
			float sym_rate = ((float[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_SYM_RATE])[gender][site];

			// Effect of vaccine reducing symptoms
			int[] vacc_expiry = vaccine_expiry_by_indivdual.get(infectedId);
			if (vacc_expiry != null) {
				float[][][] vaccine_effect = (float[][][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_VACCINE_PROPERTIES];
				if (infectious_time < vacc_expiry[site]) {
					sym_rate *= vaccine_effect[gender][site][VACCINE_PROPERTIES_SYMPTOM_ADJ];
				}
			}

			if (sym_rate > 0) {
				if (RNG.nextFloat() < sym_rate) {
					int sym_test_day = (int) Math.round(infectious_time + sym_test_period_by_gender[gender].sample());
					sch = schedule_testing.get(sym_test_day);

					if (sch == null) {
						sch = new ArrayList<>();
						schedule_testing.put(sym_test_day, sch);
					}

					int keyS = Collections.binarySearch(sch, infectedId);
					if (keyS < 0) {
						keyS = Collections.binarySearch(sch, -infectedId);
						if (keyS < 0) {
							sch.add(~keyS, -infectedId);
						}
					}

				}

			}

			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_INFECTION_HISTORY) > 0) {
				ArrayList<Integer>[] hist = infection_history.get(infectedId);
				if (hist == null) {
					hist = new ArrayList[LENGTH_SITE];
					for (int s = 0; s < LENGTH_SITE; s++) {
						hist[s] = new ArrayList<>();
					}
					infection_history.put(infectedId, hist);
				}
				hist[site].add(infectious_time);
			}

		}
		return key;
	}

	public int removeInfected(Integer infectedId, int site, int recoverTime) {
		int key = Collections.binarySearch(currently_infectious[site], infectedId);
		if (key >= 0) {
			currently_infectious[site].remove(key);
			updateScheduleMap(infectedId, LENGTH_SITE + site, null);
			if (key >= 0) {
				if ((simSetting
						& 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_INFECTION_HISTORY) > 0) {
					ArrayList<Integer>[] hist = infection_history.get(infectedId);
					hist[site].add(recoverTime);
				}
			}

			int gender = getGenderType(infectedId);

			float[] non_viable_infection_setting_by_gender_site = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_NON_VIABLE_INFECTION_SETTING])[gender][site];

			if (non_viable_infection_setting_by_gender_site[NON_VIABILITY_RECOVERY_INDUDCED_PROB] > 0) {
				if (has_non_viable_bacteria_until.getOrDefault(infectedId, recoverTime) <= recoverTime) {
					if (RNG.nextFloat() < non_viable_infection_setting_by_gender_site[NON_VIABILITY_RECOVERY_INDUDCED_PROB]) {
						double hasNABacteria = non_viable_inf_by_recovery_duration[gender][site].sample();
						has_non_viable_bacteria_until.put(infectedId, (int) Math.round(recoverTime + hasNABacteria));
					}
				}
			}

		}
		return key;
	}

	public int getGenderType(Integer personId) {
		return getGenderType(personId, cumulative_pop_composition);
	}

	public static int getGenderType(Integer personId, int[] cumul_pop_comp) {
		int index = Arrays.binarySearch(cumul_pop_comp, personId);

		if (index < 0) {
			return ~index;
		} else {
			return index;
		}
	}

	@Override
	public void run() {
		int startTime = firstSeedTime;
		float[][][] vaccine_effect = (float[][][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_VACCINE_PROPERTIES];

		if (startTime < Integer.MAX_VALUE) {

			Object[] simulation_store = preSimulation();

			Integer[] switchTime = propSwitch_map.keySet().toArray(new Integer[propSwitch_map.size()]);

			Arrays.sort(switchTime);
			int switchTimeIndex = 0;

			// Current contact map
			ContactMap cMap = new ContactMap();

			if (edges_list == null) {

				try {
					edges_list = generateMapEdgeArray(BASE_CONTACT_MAP).call();
				} catch (Exception e) {
					e.printStackTrace(System.err);
					System.err.println("Error in generating edge list from BASE_CONTACT_MAP. Exiting...");
					edges_list = new ArrayList<>();
					System.exit(-1);

				}
			}
			Integer[][] edges_array = edges_list.toArray(new Integer[edges_list.size()][]);
			int edges_array_pt = 0;

			HashMap<Integer, ArrayList<Integer[]>> removeEdges = new HashMap<>();
			ArrayList<Integer[]> toRemove;

			// Skip invalid edges
			while (edges_array_pt < edges_array.length
					&& edges_array[edges_array_pt][Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] < startTime) {
				Integer[] edge = edges_array[edges_array_pt];
				int edge_start_time = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
				int expireAt = edge_start_time + edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];

				if (expireAt > startTime) {
					toRemove = removeEdges.get(expireAt);
					if (toRemove == null) {
						toRemove = new ArrayList<>();
						removeEdges.put(expireAt, toRemove);
					}
					toRemove.add(edge);

					for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
							Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
						if (!cMap.containsVertex(edge[index])) {
							cMap.addVertex(edge[index]);
						}
					}
					cMap.addEdge(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
							edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], edge);
				}

				edges_array_pt++;
			}

			// Schedule testing
			for (Integer personId : BASE_CONTACT_MAP.vertexSet()) {
				scheduleNextTest(personId, startTime);
			}
			int snap_index = 0;

			boolean hasInfected = hasInfectedInPop();

			ArrayList<Integer> infected_today = new ArrayList<>();

			int[][] cumul_incidence = new int[Population_Bridging.LENGTH_GENDER][LENGTH_SITE];
			int[] cumul_incidence_by_person = new int[Population_Bridging.LENGTH_GENDER];
			int[][] cumul_antibiotic_use = new int[Population_Bridging.LENGTH_GENDER][2];
			int[] cumul_treatment_by_person = new int[Population_Bridging.LENGTH_GENDER * 2];
			int[] cumul_positive_dx_by_person = new int[Population_Bridging.LENGTH_GENDER * 2];
			int[] cumul_positive_dx_sought_by_person = new int[Population_Bridging.LENGTH_GENDER * 2];

			for (int currentTime = startTime; currentTime < startTime + NUM_TIME_STEPS_PER_SNAP * NUM_SNAP
					&& hasInfected; currentTime++) {

				infected_today.clear();

				if (switchTimeIndex < switchTime.length && switchTime[switchTimeIndex] == currentTime) {
					HashMap<Integer, String> switch_ent = propSwitch_map.get(currentTime);
					for (Integer switch_index : switch_ent.keySet()) {
						String str_obj = switch_ent.get(switch_index);
						getRunnable_fields()[switch_index - RUNNABLE_OFFSET] = PropValUtils.propStrToObject(str_obj,
								getRunnable_fields()[switch_index - RUNNABLE_OFFSET].getClass());
					}
					switchTimeIndex++;
				}

				// Remove expired edges
				toRemove = removeEdges.get(currentTime);
				if (toRemove != null) {
					for (Integer[] edge : toRemove) {
						cMap.removeEdge(edge);
					}
				}

				// Add new edges and update removal schedule
				while (edges_array_pt < edges_array.length
						&& edges_array[edges_array_pt][Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] <= currentTime) {

					Integer[] edge = edges_array[edges_array_pt];
					Integer expireAt = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
							+ Math.max(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION], 1);

					toRemove = removeEdges.get(expireAt);
					if (toRemove == null) {
						toRemove = new ArrayList<>();
						removeEdges.put(expireAt, toRemove);
					}
					toRemove.add(edge);

					for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
							Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
						if (!cMap.containsVertex(edge[index])) {
							cMap.addVertex(edge[index]);
						}
					}

					cMap.addEdge(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
							edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], edge);

					edges_array_pt++;
				}

				for (int site_src = 0; site_src < LENGTH_SITE; site_src++) {
					// Update infectious
					ArrayList<Integer> becomeInfectiousToday = schedule_incubation[site_src].remove(currentTime);

					if (becomeInfectiousToday != null) {
						for (Integer toInfectiousId : becomeInfectiousToday) {
							int recoveredAt = (int) Math.round(infectious_period[site_src].sample()) + currentTime;

							// Reduction of infection duration due to vaccine
							int[] vac_expiry = vaccine_expiry_by_indivdual.get(toInfectiousId);
							if (vac_expiry != null) {
								if (currentTime < vac_expiry[site_src]) {
									recoveredAt = currentTime
											+ Math.round((recoveredAt - currentTime) * vaccine_effect[getGenderType(
													toInfectiousId)][site_src][VACCINE_PROPERTIES_INF_DURATION_ADJUST]);
								}
							}

							addInfectious(toInfectiousId, site_src, currentTime, recoveredAt);
						}
					}

					// Update recovery
					ArrayList<Integer> recoveredToday = schedule_recovery[site_src].remove(currentTime);
					if (recoveredToday != null) {
						for (Integer toRecoveredId : recoveredToday) {
							removeInfected(toRecoveredId, site_src, currentTime);

						}
					}

					// Transmission
					for (Integer infectious : currently_infectious[site_src]) {
						if (cMap.containsVertex(infectious)) {
							Set<Integer[]> edges = cMap.edgesOf(infectious);
							double[][] trans = trans_prob.get(infectious);

							int[] vaccine_expiry_src = vaccine_expiry_by_indivdual.get(infectious);

							for (Integer[] e : edges) {

								int partner = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1].equals(infectious)
										? e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]
										: e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];

								int g_s = getGenderType(infectious);
								int g_t = getGenderType(partner);

								int[] valid_target = g_t == Population_Bridging.GENDER_FEMALE
										? new int[] { SITE_VAGINA, SITE_RECTUM, SITE_OROPHARYNX }
										: new int[] { SITE_PENIS, SITE_RECTUM, SITE_OROPHARYNX };

								for (int site_target : valid_target) {
									if (trans[site_src][site_target] != 0) {

										// Transmission is possible
										if (Collections.binarySearch(currently_infectious[site_target], partner) < 0) {

											int actType;

											// Determine act type
											switch (site_src) {
											case SITE_VAGINA:
												actType = ACT_INDEX_GENITAL;
												break;
											case SITE_OROPHARYNX:
												switch (site_target) {
												case SITE_RECTUM:
													actType = ACT_INDEX_RIMMING;
													break;
												case SITE_OROPHARYNX:
													actType = ACT_INDEX_KISSING;
													break;
												default:
													actType = ACT_INDEX_FELLATIO;
												}
												break;
											case SITE_RECTUM:
												switch (site_target) {
												case SITE_OROPHARYNX:
													actType = ACT_INDEX_RIMMING;
													break;
												default:
													actType = ACT_INDEX_ANAL;
												}
												break;
											default:
												switch (site_target) {
												case SITE_VAGINA:
													actType = ACT_INDEX_GENITAL;
													break;
												case SITE_OROPHARYNX:
													actType = ACT_INDEX_FELLATIO;
													break;
												case SITE_RECTUM:
													actType = ACT_INDEX_ANAL;
													break;
												default:
													actType = -1;
												}
											}
											boolean transmission_possible = actType != -1;

											if (transmission_possible) {

												int[] vaccine_expiry_target = vaccine_expiry_by_indivdual.get(partner);

												float actProb = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ])[actType][g_s][g_t];
												float transProb = (float) trans[site_src][site_target];
												transmission_possible &= actProb > 0;

												if (transmission_possible) {
													transmission_possible &= RNG.nextFloat() < actProb;

													if (transmission_possible) {

														// Reduction of transProb due to vaccine
														if (vaccine_expiry_src != null) {
															if (currentTime < vaccine_expiry_src[site_src]) {
																transProb *= vaccine_effect[g_s][site_src][VACCINE_PROPERTIES_TRANSMISSION_ADJ];
															}
														}
														if (vaccine_expiry_target != null) {
															if (currentTime < vaccine_expiry_target[site_target]) {
																transProb *= vaccine_effect[g_t][site_target][VACCINE_PROPERTIES_SUSCEPTIBILITY_ADJUST];
															}
														}

														transmission_possible &= RNG.nextFloat() < transProb;

														if (!transmission_possible) {
															// Passing on non-viable bacteria

															if (has_non_viable_bacteria_until.getOrDefault(partner,
																	currentTime) <= currentTime) {
																// Acted but no transmission
																float[] non_viable_infection_setting_by_gender_site = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_NON_VIABLE_INFECTION_SETTING])[g_t][site_target];

																if (non_viable_infection_setting_by_gender_site[NON_VIABILITY_CONTACT_INDUCED_PROB] > 0) {
																	if (RNG.nextFloat() < non_viable_infection_setting_by_gender_site[NON_VIABILITY_CONTACT_INDUCED_PROB]) {
																		double hasNABacteria = non_viable_inf_by_contact_duration[g_t][site_target]
																				.sample();
																		has_non_viable_bacteria_until.put(partner,
																				(int) Math.round(
																						currentTime + hasNABacteria));

																	}
																}
															}

														} else {
															cumul_incidence[g_t][site_target]++;
															int k = Collections.binarySearch(infected_today, partner);
															if (k < 0) {
																cumul_incidence_by_person[g_t]++;
																infected_today.add(~k, partner);
															}
															transmission_success(currentTime, infectious, partner,
																	site_target, actType, simulation_store);
														}

													}
												}
											}

										}
									}
								}
							}
						}
					}
				}

				// Testing
				ArrayList<Integer> testToday = schedule_testing.remove(currentTime);
				if (testToday != null) {
					for (Integer tId : testToday) {
						int testOutcome = testPerson(currentTime, tId);

						int gI = getGenderType(Math.abs(tId));
						if ((testOutcome & 1 << TEST_OUTCOME_TREATMENT_APPLIED) != 0) {
							cumul_treatment_by_person[gI]++;
						}
						if ((testOutcome & 1 << TEST_OUTCOME_TREATMENT_APPLIED_ON_TRUE_INFECTION) != 0) {
							cumul_treatment_by_person[Population_Bridging.LENGTH_GENDER + gI]++;
						}
						if ((testOutcome & 1 << TEST_OUTCOME_POSITIVE_TEST) != 0) {
							if (tId < 0) {
								cumul_positive_dx_sought_by_person[gI]++;
							} else {
								cumul_positive_dx_by_person[gI]++;
							}
						}

						if ((testOutcome & 1 << TEST_OUTCOME_POSITIVE_TEST_ON_TRUE_INFECTION) != 0) {
							if (tId < 0) {
								cumul_positive_dx_sought_by_person[Population_Bridging.LENGTH_GENDER + gI]++;
							} else {
								cumul_positive_dx_by_person[Population_Bridging.LENGTH_GENDER + gI]++;
							}
						}

					}
				}

				// Antibiotic clearance
				if ((simSetting
						& 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE) != 0) {
					ArrayList<Integer> antibiotic_clearance_ent = schedule_antibiotic_clearance.remove(currentTime);
					if (antibiotic_clearance_ent != null) {
						for (Integer antibiotic_clearance : antibiotic_clearance_ent) {
							int key = Collections.binarySearch(currently_has_antibiotic, antibiotic_clearance);
							if (key >= 0) {
								currently_has_antibiotic.remove(key);
							}
						}
					}

					for (Integer antibiotic_user : currently_has_antibiotic) {
						int usage_key = antibiotic_user > 0 ? 0 : 1;
						cumul_antibiotic_use[getGenderType(Math.abs(antibiotic_user))][usage_key]++;
					}

				}

				// Vaccine or vaccine booster
				ArrayList<Integer> vaccineToday = schedule_vaccination.remove(currentTime);
				if (vaccineToday != null) {
					for (Integer vId : vaccineToday) {
						vaccine_person(vId, currentTime);
					}
				}

				// Storing snapshot infected in sim_output

				if (snap_index == 0) {
					@SuppressWarnings("unchecked")
					// K = time, V= int [gender][site]
					HashMap<Integer, int[][]> infectious_count_map = (HashMap<Integer, int[][]>) sim_output
							.get(SIM_OUTPUT_INFECTIOUS_COUNT);

					if (infectious_count_map == null) {
						infectious_count_map = new HashMap<>();
						sim_output.put(SIM_OUTPUT_INFECTIOUS_COUNT, infectious_count_map);
					}

					@SuppressWarnings("unchecked")
					// K = time, V= int [gender]
					HashMap<Integer, int[]> infectious_count_map_by_person = (HashMap<Integer, int[]>) sim_output
							.get(SIM_OUTPUT_INFECTIOUS_COUNT_BY_PERSON);

					if (infectious_count_map_by_person == null) {
						infectious_count_map_by_person = new HashMap<>();
						sim_output.put(SIM_OUTPUT_INFECTIOUS_COUNT_BY_PERSON, infectious_count_map_by_person);
					}

					int[][] infectious_count = new int[cumulative_pop_composition.length][LENGTH_SITE];
					int[] infectious_count_by_person = new int[cumulative_pop_composition.length];

					ArrayList<Integer> infected_person = new ArrayList<>();
					for (int site = 0; site < LENGTH_SITE; site++) {
						for (Integer infected_id : currently_infectious[site]) {
							int gender_type = getGenderType(infected_id);
							infectious_count[gender_type][site]++;
							int k = Collections.binarySearch(infected_person, infected_id);
							if (k < 0) {
								infectious_count_by_person[gender_type]++;
								infected_person.add(~k, infected_id);
							}
						}
					}
					infectious_count_map.put(currentTime, infectious_count);
					infectious_count_map_by_person.put(currentTime, infectious_count_by_person);

					@SuppressWarnings("unchecked")
					HashMap<Integer, int[][]> cumul_incidence_map = (HashMap<Integer, int[][]>) sim_output
							.get(SIM_OUTPUT_CUMUL_INCIDENCE);
					if (cumul_incidence_map == null) {
						cumul_incidence_map = new HashMap<>();
						sim_output.put(SIM_OUTPUT_CUMUL_INCIDENCE, cumul_incidence_map);
					}

					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> cumul_incidence_map_by_person = (HashMap<Integer, int[]>) sim_output
							.get(SIM_OUTPUT_CUMUL_INCIDENCE_BY_PERSON);
					if (cumul_incidence_map_by_person == null) {
						cumul_incidence_map_by_person = new HashMap<>();
						sim_output.put(SIM_OUTPUT_CUMUL_INCIDENCE_BY_PERSON, cumul_incidence_map_by_person);
					}

					int[][] incidence_snap = new int[cumul_incidence.length][];

					for (int g = 0; g < cumul_incidence.length; g++) {
						incidence_snap[g] = Arrays.copyOf(cumul_incidence[g], cumul_incidence[g].length);
					}

					cumul_incidence_map.put(currentTime, incidence_snap);
					cumul_incidence_map_by_person.put(currentTime,
							Arrays.copyOf(cumul_incidence_by_person, cumul_incidence_by_person.length));

					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> cumul_treatment_map_by_person = (HashMap<Integer, int[]>) sim_output
							.get(SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON);
					if (cumul_treatment_map_by_person == null) {
						cumul_treatment_map_by_person = new HashMap<>();
						sim_output.put(SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON, cumul_treatment_map_by_person);
					}
					cumul_treatment_map_by_person.put(currentTime,
							Arrays.copyOf(cumul_treatment_by_person, cumul_treatment_by_person.length));										
					
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> cumul_positive_dx_map_by_person = (HashMap<Integer, int[]>) sim_output
							.get(SIM_OUTPUT_CUMUL_POS_DX_BY_PERSON);
					if (cumul_positive_dx_map_by_person == null) {
						cumul_positive_dx_map_by_person = new HashMap<>();
						sim_output.put(SIM_OUTPUT_CUMUL_POS_DX_BY_PERSON, cumul_positive_dx_map_by_person);
					}
					cumul_positive_dx_map_by_person.put(currentTime,
							Arrays.copyOf(cumul_positive_dx_by_person, cumul_positive_dx_by_person.length));
					
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> cumul_positive_dx_sought_map_by_person = (HashMap<Integer, int[]>) sim_output
							.get(SIM_OUTPUT_CUMUL_POS_DX_SOUGHT_BY_PERSON);
					if (cumul_positive_dx_sought_map_by_person == null) {
						cumul_positive_dx_sought_map_by_person = new HashMap<>();
						sim_output.put(SIM_OUTPUT_CUMUL_POS_DX_SOUGHT_BY_PERSON, cumul_positive_dx_sought_map_by_person);
					}
					cumul_positive_dx_map_by_person.put(currentTime,
							Arrays.copyOf(cumul_positive_dx_sought_by_person, cumul_positive_dx_sought_by_person.length));
					
					

					if ((simSetting
							& 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE) != 0) {

						@SuppressWarnings("unchecked")
						HashMap<Integer, int[][]> cumul_antibiotic_map = (HashMap<Integer, int[][]>) sim_output
								.get(SIM_OUTPUT_CUMUL_ANTIBOTIC_USAGE);

						if (cumul_antibiotic_map == null) {
							cumul_antibiotic_map = new HashMap<>();
							sim_output.put(SIM_OUTPUT_CUMUL_ANTIBOTIC_USAGE, cumul_antibiotic_map);
						}

						int[][] antibiotic_usage_snap = new int[cumul_antibiotic_use.length][];
						for (int g = 0; g < cumul_incidence.length; g++) {
							antibiotic_usage_snap[g] = Arrays.copyOf(cumul_antibiotic_use[g],
									cumul_antibiotic_use[g].length);
						}
						cumul_antibiotic_map.put(currentTime, antibiotic_usage_snap);

					}

					if ((simSetting
							& 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_INFECTION_HISTORY) != 0) {
						sim_output.put(SIM_OUTPUT_INFECTION_HISTORY, infection_history);

					}

					if ((simSetting
							& 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_VACCINE_COVERAGE) != 0) {
						@SuppressWarnings("unchecked")
						HashMap<Integer, int[][][]> vaccine_coverage_map = (HashMap<Integer, int[][][]>) sim_output
								.get(SIM_OUTPUT_VACCINE_COVERAGE);
						if (vaccine_coverage_map == null) {
							vaccine_coverage_map = new HashMap<>();
							sim_output.put(SIM_OUTPUT_VACCINE_COVERAGE, vaccine_coverage_map);
						}

						int[][][] vaccine_coverage_ent = new int[Population_Bridging.LENGTH_GENDER][LENGTH_SITE][2];

						for (Integer pid : currently_vaccinated) {
							int g = getGenderType(pid);
							int[] vaccine_expiry = vaccine_expiry_by_indivdual.get(pid);
							for (int s = 0; s < LENGTH_SITE; s++) {
								if (currentTime < vaccine_expiry[s]) {
									vaccine_coverage_ent[g][s][0]++;
								} else {
									vaccine_coverage_ent[g][s][1]++;
								}

							}
						}
						vaccine_coverage_map.put(currentTime, vaccine_coverage_ent);

					}

				}

				snap_index = (snap_index + 1) % NUM_TIME_STEPS_PER_SNAP;

				hasInfected = hasInfectedInPop();

			}
			// End of simulations
			if (runnableId != null) {
				System.out.println(String.format("Thread <%s> completed.", runnableId));
			}
			postSimulation(simulation_store);

		}

	}

	protected int testPerson(int currentTime, Integer testing_pid_signed) {
		Integer test_pid = Math.abs(testing_pid_signed);

		int gender = getGenderType(test_pid);

		boolean positve_dx = false;
		boolean positiv_dx_on_infectious = false;
		boolean applyTreatment = false;
		boolean applyTreatment_on_infectious = false;

		// float[dx_sensitivity, dx_specificity][gender][site]
		// or float[dx_sensitivity, dx_specificity, treatment_rate][gender][site]
		float[][][] test_accuracy = (float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_DX_TEST_ACCURACY];

		// Treatment due to DX
		int[] infectious_key_index = new int[LENGTH_SITE];
		Arrays.fill(infectious_key_index, -1);

		for (int site = 0; site < LENGTH_SITE; site++) {
			infectious_key_index[site] = Collections.binarySearch(currently_infectious[site], test_pid);
			if (infectious_key_index[site] >= 0) {
				// True infectious

				// Positive test
				boolean posDx = test_accuracy[DX_SENSITIVITY_INDEX][gender][site] >= 1;

				if (!posDx) {
					posDx = RNG.nextFloat() < test_accuracy[DX_SENSITIVITY_INDEX][gender][site];
				}

				positve_dx |= posDx;
				positiv_dx_on_infectious |= posDx;

				// Apply treatment
				if (posDx) {
					boolean applyTr = !(DX_TREATMENT_RATE_INDEX < test_accuracy.length);
					if (!applyTr) {
						applyTr = test_accuracy[DX_TREATMENT_RATE_INDEX][gender][site] >= 1;
						if (!applyTr) {
							applyTr = RNG.nextFloat() < test_accuracy[DX_TREATMENT_RATE_INDEX][gender][site];
						}
					}
					if (applyTr) {
						applyTreatment |= applyTr;
						applyTreatment_on_infectious |= applyTr;
					}
				}

			} else if (!applyTreatment
					&& has_non_viable_bacteria_until.getOrDefault(test_pid, currentTime) > currentTime) {
				// Overtreatment due to non-viable
				boolean posDx = test_accuracy[DX_SPECIFICITY_INDEX][gender][site] < 1;
				if (posDx) {
					posDx = RNG.nextFloat() >= test_accuracy[DX_SPECIFICITY_INDEX][gender][site];
					positve_dx |= posDx;
					if (posDx) {
						boolean applyTr = !(DX_TREATMENT_RATE_INDEX < test_accuracy.length);
						if (!applyTr) {
							applyTr = test_accuracy[DX_TREATMENT_RATE_INDEX][gender][site] >= 1;
							if (!applyTr) {
								applyTr = RNG.nextFloat() < test_accuracy[DX_TREATMENT_RATE_INDEX][gender][site];
							}
						}
						if (applyTr) {
							applyTreatment |= applyTr;
						}
					}
				}

			}
		}

		if (applyTreatment) {

			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE) != 0) {
				int key;
				key = Collections.binarySearch(currently_has_antibiotic,
						applyTreatment_on_infectious ? test_pid : -test_pid);
				if (key < 0) {
					currently_has_antibiotic.add(~key, applyTreatment_on_infectious ? test_pid : -test_pid);
					int antibiotic_flush_at = (int) Math.round(currentTime + antibotic_duration[gender].sample());
					ArrayList<Integer> sch_antibiotic_flush = schedule_antibiotic_clearance.get(antibiotic_flush_at);
					if (sch_antibiotic_flush == null) {
						sch_antibiotic_flush = new ArrayList<>();
						schedule_antibiotic_clearance.put(antibiotic_flush_at, sch_antibiotic_flush);
					}
					sch_antibiotic_flush.add(applyTreatment_on_infectious ? test_pid : -test_pid);

				}

			}

			// infection_schMap = {incubation by site_0, ... recovery_by_site_0 ...}
			Integer[] infection_schMap = mapping_infection_schedule.get(test_pid);

			// float[gender][site][setting]
			float[][] non_viable_infection_setting_by_gender = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_NON_VIABLE_INFECTION_SETTING])[gender];

			// Assemble all non-viable bacteria were clear at treatment first
			if (has_non_viable_bacteria_until.containsKey(test_pid)) {
				has_non_viable_bacteria_until.put(test_pid, currentTime);
			}

			for (int site = 0; site < LENGTH_SITE; site++) {
				if (infectious_key_index[site] >= 0) {
					// Remove from infectious
					currently_infectious[site].remove(infectious_key_index[site]);
					// Treatment induced non-viability
					if (non_viable_infection_setting_by_gender[site][NON_VIABILITY_TREATMENT_INDUCED_PROB] > 0) {
						if (RNG.nextFloat() < non_viable_infection_setting_by_gender[site][NON_VIABILITY_TREATMENT_INDUCED_PROB]) {
							Integer recoverDate_org = infection_schMap[LENGTH_SITE + site];
							Integer has_non_viable = has_non_viable_bacteria_until.getOrDefault(test_pid, currentTime);

							int non_viable_dur = Math.max(has_non_viable, currentTime + Math.round((recoverDate_org
									- currentTime)
									* non_viable_infection_setting_by_gender[site][NON_VIABILITY_TREATMENT_INDUCED_INFECTION_REDUCTION_ADJ]));

							has_non_viable_bacteria_until.put(test_pid, non_viable_dur);

						}
					}
				}
			}

			// Remove from incubation and recovery

			if (infection_schMap != null) {
				for (int i = 0; i < infection_schMap.length; i++) {
					if (infection_schMap[i] != null) {
						Integer dateEnt = infection_schMap[i];
						int site = i % LENGTH_SITE;
						ArrayList<Integer> schArr;
						if (i < LENGTH_SITE) {
							schArr = schedule_incubation[site].get(dateEnt);
						} else {
							schArr = schedule_recovery[site].get(dateEnt);
						}
						int key = Collections.binarySearch(schArr, test_pid);
						if (key >= 0) {
							schArr.remove(key);
						}

					}
				}
				Arrays.fill(infection_schMap, null);
			}

		}

		if (testing_pid_signed < 0) {
			// Schedule next test
			scheduleNextTest(test_pid, currentTime);
		}

		int res = 0;

		if (positve_dx) {
			res |= 1 << TEST_OUTCOME_POSITIVE_TEST;
		}
		if (positiv_dx_on_infectious) {
			res |= 1 << TEST_OUTCOME_POSITIVE_TEST_ON_TRUE_INFECTION;
		}
		if (applyTreatment) {
			res |= 1 << TEST_OUTCOME_TREATMENT_APPLIED;
		}
		if (applyTreatment_on_infectious) {
			res |= 1 << TEST_OUTCOME_TREATMENT_APPLIED_ON_TRUE_INFECTION;
		}

		// Vaccination by test
		float[] vacc_setting = ((float[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_VACCINE_SETTING])[gender];
		float vacc_coverage_by_test = vacc_setting[VACCINATION_SETTING_RATE_PER_TEST];

		int key = Collections.binarySearch(currently_vaccinated, test_pid);

		// Only vaccinate possible if they have not done so before
		if (key < 0) {
			if (vacc_coverage_by_test > 0 && RNG.nextFloat() < vacc_coverage_by_test) {
				vaccine_person(test_pid, currentTime);
			}
		}

		return res;
	}

	protected void vaccine_person(Integer vaccinate_pid, int currentTime) {
		int key = Collections.binarySearch(currently_vaccinated, vaccinate_pid);
		int gender = getGenderType(vaccinate_pid);
		float[] vacc_setting = ((float[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_VACCINE_SETTING])[gender];

		if (key < 0) {
			// Schedule booster on first run
			currently_vaccinated.add(~key, vaccinate_pid);
			if (vacc_setting[VACCINATION_SETTING_BOOSTER_PERIOD] > 0) {
				if (vacc_setting[VACCINATION_SETTING_BOOSTER_LIMIT] > 0) {
					for (int b = 0; b < vacc_setting[VACCINATION_SETTING_BOOSTER_LIMIT]; b++) {
						int booster_date = Math
								.round(currentTime + (b + 1) * vacc_setting[VACCINATION_SETTING_BOOSTER_PERIOD]);
						ArrayList<Integer> booster_candidate = schedule_vaccination.get(booster_date);
						if (booster_candidate == null) {
							booster_candidate = new ArrayList<>();
							schedule_vaccination.put(booster_date, booster_candidate);
						}
						int b_key = Collections.binarySearch(booster_candidate, vaccinate_pid);
						if (b_key < 0) {
							booster_candidate.add(~b_key, vaccinate_pid);
						}
					}
				}
			}
		}

		if (vacc_setting[VACCINATION_SETTING_BOOSTER_PERIOD] > 0
				&& vacc_setting[VACCINATION_SETTING_BOOSTER_LIMIT] < 0) { // Inf. booster
			int booster_date = Math.round(currentTime + vacc_setting[VACCINATION_SETTING_BOOSTER_PERIOD]);
			ArrayList<Integer> booster_candidate = schedule_vaccination.get(booster_date);
			if (booster_candidate == null) {
				booster_candidate = new ArrayList<>();
				schedule_vaccination.put(booster_date, booster_candidate);
			}
			int b_key = Collections.binarySearch(booster_candidate, vaccinate_pid);
			if (b_key < 0) {
				booster_candidate.add(~b_key, vaccinate_pid);
			}
		}

		int[] vacc_expiry = vaccine_expiry_by_indivdual.get(vaccinate_pid);
		if (vacc_expiry == null) {
			vacc_expiry = new int[LENGTH_SITE];
			vaccine_expiry_by_indivdual.put(vaccinate_pid, vacc_expiry);
		}
		for (int site = 0; site < LENGTH_SITE; site++) {
			if (vacc_dur[gender][site] != null) {
				vacc_expiry[site] = (int) Math.round(currentTime + vacc_dur[gender][site].sample());
			} else {
				vacc_expiry[site] = 0;
			}

		}
	}

	protected boolean hasInfectedInPop() {
		boolean hasInfected = false;
		for (int site_src = 0; site_src < LENGTH_SITE && !hasInfected; site_src++) {
			hasInfected |= !(schedule_incubation[site_src].isEmpty() && currently_infectious[site_src].isEmpty());
		}
		return hasInfected;
	}

	protected Object[] preSimulation() {
		// Do nothing by default
		return new Object[0];

	}

	/**
	 * Procedure that are called after simulations. *
	 * 
	 */
	@SuppressWarnings("unchecked")
	protected void postSimulation(Object[] simulation_store) {
		PrintWriter pWri;
		HashMap<Integer, int[][]> count_map;
		HashMap<Integer, int[]> count_map_by_person;
		StringBuilder str = null;

		try {
			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE) != 0) {
				count_map = (HashMap<Integer, int[][]>) sim_output.get(SIM_OUTPUT_INFECTIOUS_COUNT);
				str = printCountMap(count_map, "Gender_%d_Site_%d");
				pWri = new PrintWriter(new File(baseDir, String
						.format(Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();

				count_map_by_person = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_INFECTIOUS_COUNT_BY_PERSON);
				str = printCountMap(count_map_by_person, Population_Bridging.LENGTH_GENDER, "Gender_%d");
				pWri = new PrintWriter(new File(baseDir, String
						.format(Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();

			}
			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE) != 0) {
				count_map = (HashMap<Integer, int[][]>) sim_output.get(SIM_OUTPUT_CUMUL_INCIDENCE);
				str = printCountMap(count_map, "Gender_%d_Site_%d");
				pWri = new PrintWriter(new File(baseDir, String.format(
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_SITE, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();

				count_map_by_person = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_CUMUL_INCIDENCE_BY_PERSON);
				str = printCountMap(count_map_by_person, Population_Bridging.LENGTH_GENDER, "Gender_%d");
				pWri = new PrintWriter(new File(baseDir, String.format(
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();
			}
			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE) != 0) {
				
				count_map_by_person = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_CUMUL_POS_DX_BY_PERSON);				
				str = printCountMap(count_map_by_person, Population_Bridging.LENGTH_GENDER * 2,
						new String[] { "Total_Positive_DX_Gender_%d", "True_Positive_DX_Gender_%d" });
				pWri = new PrintWriter(new File(baseDir, String.format(
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_PERSON, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();
				
				count_map_by_person = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_CUMUL_POS_DX_SOUGHT_BY_PERSON);				
				str = printCountMap(count_map_by_person, Population_Bridging.LENGTH_GENDER * 2,
						new String[] { "Total_Positive_DX_SOUGHT_Gender_%d", "True_Positive_DX_SOUGHT_Gender_%d" });
				pWri = new PrintWriter(new File(baseDir, String.format(
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();
				

				count_map_by_person = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON);
				str = printCountMap(count_map_by_person, Population_Bridging.LENGTH_GENDER * 2,
						new String[] { "Total_Treatment_Gender_%d", "True_Treatment_Gender_%d" });
				pWri = new PrintWriter(new File(baseDir, String.format(
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();

			}

			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_INFECTION_HISTORY) != 0) {
				Integer[] pids = infection_history.keySet().toArray(new Integer[infection_history.size()]);
				Arrays.sort(pids);
				pWri = new PrintWriter(new File(baseDir, String
						.format(Simulation_ClusterModelTransmission.FILENAME_INFECTION_HISTORY, cMap_seed, sim_seed)));
				for (Integer pid : pids) {
					ArrayList<Integer>[] hist = infection_history.get(pid);
					for (int site = 0; site < hist.length; site++) {
						pWri.print(pid.toString());
						pWri.print(',');
						pWri.print(site);
						for (Integer timeEnt : hist[site]) {
							pWri.print(',');
							pWri.print(timeEnt);
						}
						pWri.println();
					}
				}

				pWri.close();

			}

			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE) != 0) {
				count_map = (HashMap<Integer, int[][]>) sim_output.get(SIM_OUTPUT_CUMUL_ANTIBOTIC_USAGE);
				str = printCountMap(count_map, new int[] { Population_Bridging.LENGTH_GENDER, 2 },
						"Gender_%d_Usage_%d"); // Proper, Over treatment
				pWri = new PrintWriter(new File(baseDir, String.format(
						Simulation_ClusterModelTransmission.FILENAME_CUMUL_ANTIBIOTIC_USAGE, cMap_seed, sim_seed)));
				pWri.println(str.toString());
				pWri.close();
			}

			if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TRACK_VACCINE_COVERAGE) != 0) {
				HashMap<Integer, int[][][]> count_map_vacc = (HashMap<Integer, int[][][]>) sim_output
						.get(SIM_OUTPUT_VACCINE_COVERAGE);
				pWri = new PrintWriter(new File(baseDir, String
						.format(Simulation_ClusterModelTransmission.FILENAME_VACCINE_COVERAGE, cMap_seed, sim_seed)));
				Integer[] time_array = count_map_vacc.keySet().toArray(new Integer[count_map_vacc.size()]);
				Arrays.sort(time_array);
				str = new StringBuilder();

				str.append("Time");
				for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
					for (int s = 0; s < LENGTH_SITE; s++) {
						for (String v : new String[] { "valid", "expired" }) {
							str.append(',');
							str.append(g);
							str.append('_');
							str.append(s);
							str.append('_');
							str.append(v);
						}
					}
				}
				str.append('\n');
				for (Integer t : time_array) {
					int[][][] ent = count_map_vacc.get(t);
					str.append(t);
					for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
						for (int s = 0; s < LENGTH_SITE; s++) {
							for (int v = 0; v < ent[g][s].length; v++) {
								str.append(',');
								str.append(ent[g][s][v]);
							}
						}
					}
					str.append('\n');
				}
				pWri.println(str.toString());
				pWri.close();
			}

		} catch (Exception e) {
			e.printStackTrace(System.err);
			if (str != null) {
				System.err.println("Outcome so far");
				System.err.println(str.toString());
			}

		}

	}

	protected StringBuilder printCountMap(HashMap<Integer, int[]> count_map_by_person, int dimension,
			String headingFormat) {
		return printCountMap(count_map_by_person, dimension, new String[] { headingFormat });
	}

	protected StringBuilder printCountMap(HashMap<Integer, int[]> count_map_by_person, int dimension,
			String[] headingFormat) {
		StringBuilder str;
		Integer[] time_array;
		str = new StringBuilder();

		time_array = count_map_by_person.keySet().toArray(new Integer[count_map_by_person.size()]);
		Arrays.sort(time_array);
		str.append("Time");

		int numPerSet = dimension / headingFormat.length;

		for (int g = 0; g < dimension; g++) {
			str.append(',');
			str.append(String.format(headingFormat[g / numPerSet], g % numPerSet));

		}
		str.append('\n');
		for (Integer time : time_array) {
			str.append(time);
			int[] ent = count_map_by_person.get(time);
			for (int g = 0; g < dimension; g++) {
				str.append(',');
				str.append(ent[g]);

			}
			str.append('\n');
		}
		return str;
	}

	private StringBuilder printCountMap(HashMap<Integer, int[][]> count_map, String headingFormat) {
		return printCountMap(count_map, new int[] { Population_Bridging.LENGTH_GENDER, LENGTH_SITE }, headingFormat);
	}

	private StringBuilder printCountMap(HashMap<Integer, int[][]> count_map, int[] dimension, String headingFormat) {
		// K = time, V= int [gender][site] (for example)
		Integer[] time_array;
		StringBuilder str;
		time_array = count_map.keySet().toArray(new Integer[count_map.size()]);
		Arrays.sort(time_array);

		str = new StringBuilder();
		str.append("Time");

		for (int g = 0; g < dimension[0]; g++) {
			for (int s = 0; s < dimension[1]; s++) {
				str.append(',');
				str.append(String.format(headingFormat, g, s));
			}
		}
		str.append('\n');
		for (Integer time : time_array) {
			str.append(time);
			int[][] ent = count_map.get(time);
			for (int g = 0; g < dimension[0]; g++) {
				for (int s = 0; s < dimension[1]; s++) {
					str.append(',');
					str.append(ent[g][s]);
				}
			}
			str.append('\n');
		}
		return str;
	}

	protected void transmission_success(int currentTime, Integer infectious, int partner, int site_target, int actType,
			Object[] simulation_store) {
		Integer incubation_end_at = currentTime + (int) incubation_period[site_target].sample();
		ArrayList<Integer> ent = schedule_incubation[site_target].get(incubation_end_at);
		if (ent == null) {
			ent = new ArrayList<>();
			schedule_incubation[site_target].put(incubation_end_at, ent);
		}
		int key = Collections.binarySearch(ent, partner);
		if (key < 0) {
			ent.add(~key, partner);
			updateScheduleMap(partner, site_target, incubation_end_at);
		}

	}

	protected void updateScheduleMap(int personId, int schMap_index, Integer schMap_ent) {

		Integer[] schMap = mapping_infection_schedule.get(personId);

		// key = person_id, ent = {incubation by site_0, ... recovery_by_site_0 ...}

		if (schMap == null) {
			schMap = new Integer[2 * LENGTH_SITE];
			Arrays.fill(schMap, null);
			mapping_infection_schedule.put(personId, schMap);
		}
		schMap[schMap_index] = schMap_ent;
	}

	protected final AbstractRealDistribution generateNonDistribution(double[] input) {
		return new AbstractRealDistribution(RNG) {
			private static final long serialVersionUID = -4946118496555960005L;

			@Override
			public boolean isSupportUpperBoundInclusive() {
				return true;
			}

			@Override
			public boolean isSupportLowerBoundInclusive() {
				return true;
			}

			@Override
			public boolean isSupportConnected() {
				return true;
			}

			@Override
			public double getSupportUpperBound() {
				return input[0];
			}

			@Override
			public double getSupportLowerBound() {
				return input[0];
			}

			@Override
			public double getNumericalVariance() {
				return 0;
			}

			@Override
			public double getNumericalMean() {
				return input[0];
			}

			@Override
			public double density(double x) {
				return x == input[0] ? 1 : 0;
			}

			@Override
			public double cumulativeProbability(double x) {
				return x < input[0] ? 0 : 1;
			}

			@Override
			public double sample() {
				return input[0];
			}

		};

	}

	protected AbstractRealDistribution generateGammaDistribution(double[] input) {

		if (input[1] != 0) {
			// For Gamma distribution
			// GammaDistribution(RandomGenerator rng, double shape, double scale)
			// shape = mean / scale i.e. mean / (var / mean)
			// scale = var / mean
			double[] res = new double[2];
			double var = input[1] * input[1];
			// scale
			res[1] = var / input[0];
			// shape
			res[0] = input[0] / res[1];
			return new GammaDistribution(RNG, res[0], res[1]);
		} else {
			return generateNonDistribution(input);

		}
	}

	protected AbstractRealDistribution generateBetaDistribution(double[] input) {
		if (input[1] != 0) {

			// For Beta distribution,
			// alpha = mean*(mean*(1-mean)/variance - 1)
			// beta = (1-mean)*(mean*(1-mean)/variance - 1)
			double[] res = new double[2];
			double var = input[1] * input[1];
			double rP = input[0] * (1 - input[0]) / var - 1;
			// alpha
			res[0] = rP * input[0];
			// beta
			res[1] = rP * (1 - input[0]);
			return new BetaDistribution(RNG, res[0], res[1]);
		} else {
			return generateNonDistribution(input);
		}
	}

}
