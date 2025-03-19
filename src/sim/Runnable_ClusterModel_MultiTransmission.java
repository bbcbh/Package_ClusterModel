package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.RealDistribution;

import optimisation.Optimisation_Factory;
import person.AbstractIndividualInterface;
import population.Population_Bridging;
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
			// double{
			// if STAGE_ID_TO >= 0
			// {INF_ID, SITE, STAGE_ID_FROM, STAGE_ID_TO, SWTICH_STAGE_PROB,
			// DURATION_STAGE_ID_TO_PARAM_0 ...},
			// else (instantaneous stage switch - e.g. infected in last turn)
			// {INF_ID, SITE, RISK_INC_INDEX, STAGE_ID_TO,
			// SWITCH_STAGE_1, SWITCH_STAGE_2, ... SWTICH_STAGE_PROB_1, SWTICH_STAGE_PROB_2,
			// ... ...}
			/// ...}
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
			// double{{GENDER_INCLUDE_INDEX, SEEK_TEST_MEAN, SEEK_TEST_SD}, ...}
			// or double{{GENDER_INCLUDE_INDEX, SEEK_TEST_MEAN, SEEK_TEST_SD,
			// WILL_SEEK_TREATMENT}, ...}
			new double[][] {},
			// 8 RUNNABLE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES
			// double{{INF_INDEX, SITE_INDEX, TEST_ACCURACY,
			// TARGET_STAGE_INC_1, TREATMENT_SUC_STAGE_1,TARGET_STAGE_INC_2,
			// TREATMENT_SUC_STAGE_2..}, ...}
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

	public static final int FIELD_STAGE_PERIOD_INSTANT_RISKGRP_INC = FIELD_STAGE_PERIOD_SITE + 1;
	public static final int FIELD_STAGE_PERIOD_INSTANT_STATE_TO = FIELD_STAGE_PERIOD_STATE_FROM + 1;
	public static final int FIELD_STAGE_PERIOD_INSTANT_STAGE_OPTIONS_START = FIELD_STAGE_PERIOD_INSTANT_STATE_TO + 1;

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
	public static final int FIELD_SOUGHT_TEST_PERIOD_BY_SYM_GENDER_INC = 0;
	public static final int FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SEEK_TEST_MEAN = FIELD_SOUGHT_TEST_PERIOD_BY_SYM_GENDER_INC
			+ 1;
	public static final int FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SEEK_TEST_SD = FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SEEK_TEST_MEAN
			+ 1;
	public static final int FIELD_SOUGHT_TEST_RATE = FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SEEK_TEST_SD + 1;

	// For RUNNALBE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES
	public static final int FIELD_DX_TEST_PROPERTIES_INF_ID = 0;
	public static final int FIELD_DX_TEST_PROPERTIES_SITE = FIELD_DX_TEST_PROPERTIES_INF_ID + 1;
	public static final int FIELD_DX_TEST_PROPERTIES_ACCURACY_START = FIELD_DX_TEST_PROPERTIES_SITE + 1;

	protected HashMap<Integer, double[][][][]> map_trans_prob; // Key=PID,V=double[INF_ID][SITE_FROM][SITE_TO][STAGE_ID]
	protected HashMap<Integer, int[][]> map_currrent_infection_stage; // Key=PID,V=int[INF_ID][SITE]{infection_stage}
	protected HashMap<Integer, int[][]> map_infection_stage_switch; // Key=PID,V=int[INF_ID][SITE]{switch_time_at_site_0};
	// Key=Inf_ID,SiteID,V=ArrayList of pid with infectious
	protected transient HashMap<String, ArrayList<Integer>> map_currently_infectious;

	// For schedule, key are day and entries are list of person id, ordered
	protected HashMap<Integer, ArrayList<int[]>> schedule_testing; // V=int[]{PID,INF_INCLUDE_INDEX,SITE_INCLUDE_INDEX}
	protected HashMap<Integer, ArrayList<ArrayList<ArrayList<Integer>>>> schedule_stage_change;
	// V=ArrayList<Inf><SITE>{PIds of those need to change}

	protected int[][][] cumul_incidence_by_site;
	protected int[][] cumul_incidence_by_person;
	protected int[][] cumul_treatment_by_person;
	protected int lastStateSwitch = -1;

	protected HashMap<Integer, HashMap<Integer, Integer>> test_rate_index_map; // KEY = PID, V=Map<TestRateId, Pt_ID>

	// Helper objects set by fields - Not exported
	protected transient HashMap<String, double[]> lookupTable_infection_stage_path; // Key="INF_ID,SITE_ID,STAGE_ID",V={Cumul_Prob_0,...State_ID_0...}
	protected transient int[][] lookupTable_infection_infectious_stages; // int[INF_ID][SITE]=STAGE_INCLUDE_ID
	protected transient HashMap<String, double[]> lookupTable_test_treatment_properties; // Key="INF_ID,SITE_ID,V=TEST_ACCURACY_INFO
	protected transient HashMap<String, ArrayList<double[]>> lookupTable_instant_stage_switch; // KEY="Inf_ID,SITE_ID,STAGE_ID",V
																								// = fieldEntry

	protected transient double[][][][] table_act_frequency; // double[ACT_ID][G_TO][G_FROM]=fieldEntry

	protected transient RealDistribution[][][] dist_sym_rate; // RealDistribution[INF_ID][SITE][STAGE_ID]
	protected transient RealDistribution[][][][] dist_tranmissionMatrix; // RealDistribution[INF_ID][SITE_FROM][SITE_TO][STAGE_ID]
	protected transient RealDistribution[][][] dist_stage_period; // RealDistribution[INF_ID][SITE][STAGE_ID]
	protected transient RealDistribution[] dist_seek_test_period; // RealDistribution[GENDER_ID]
	protected transient double[] seek_test_rate; // double[GENDER_ID]

	// Output setting
	protected static final int STAGE_ID_JUST_INFECTED = Integer.MIN_VALUE;

	public static final String SIM_OUTPUT_KEY_INFECTIOUS_GENDER_COUNT = "Output_%d_Gender_Count";
	public static final String SIM_OUTPUT_KEY_INFECTED_AT_GENDER_COUNT = "Output_%d_Infected_At_By_Gender";
	public static final String SIM_OUTPUT_KEY_INFECTIOUS_SITE_COUNT = "Output_%d_Infectious_Site";
	public static final String SIM_OUTPUT_KEY_INFECTED_SITE_STAGE_COUNT = "Output_%d_Infected_Site_Stage";

	protected static final String SIM_OUTPUT_KEY_CUMUL_TREATMENT = "Output_%d";
	protected static final String SIM_OUTPUT_KEY_CUMUL_INCIDENCE = "Output_%d";
	protected static final String SIM_OUTPUT_KEY_CUMUL_INCIDENCE_SITE = "Output_%d_S";

	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("MultiTransmission_(\\d+)_(\\d+)_(\\d+)");

	public Runnable_ClusterModel_MultiTransmission(long cMap_seed, long sim_seed, ContactMap base_cMap, Properties prop,
			int num_inf, int num_site, int num_act) {
		this(cMap_seed, sim_seed,
				(int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey), int[].class), base_cMap,
				Integer.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])),
				Integer.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])),
				num_inf, num_site, num_act);
	}

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
		map_infection_stage_switch = new HashMap<>();
		map_currently_infectious = new HashMap<>();

		schedule_testing = new HashMap<>();
		schedule_stage_change = new HashMap<>();

		lookupTable_infection_stage_path = new HashMap<>();
		lookupTable_infection_infectious_stages = new int[NUM_INF][NUM_SITE];
		lookupTable_test_treatment_properties = new HashMap<>();
		lookupTable_instant_stage_switch = new HashMap<>();

		table_act_frequency = new double[NUM_ACT][NUM_GENDER][NUM_GENDER][];

		dist_tranmissionMatrix = new RealDistribution[NUM_INF][NUM_SITE][NUM_SITE][];

		dist_stage_period = new RealDistribution[NUM_INF][NUM_SITE][];
		dist_sym_rate = new RealDistribution[NUM_INF][NUM_SITE][];
		dist_seek_test_period = new RealDistribution[NUM_GENDER];
		seek_test_rate = new double[NUM_GENDER];
		Arrays.fill(seek_test_rate, 1);
	}

	protected static String[] exportFileFormat = new String[] {
			Abstract_Runnable_ClusterModel_Transmission.exportFileFormat[Abstract_Runnable_ClusterModel_Transmission.EXPORT_RNG],
			Abstract_Runnable_ClusterModel_Transmission.exportFileFormat[Abstract_Runnable_ClusterModel_Transmission.EXPORT_SIMOUTPUT],
			"Export_map_trans_prob_%d_%d_%d.map", "Export_map_currrent_infection_stage_%d_%d_%d.map",
			"Export_map_infection_stage_switch_%d_%d_%d.map", "Export_schedule_testing_%d_%d_%d.map",
			"Export_schedule_stage_change_%d_%d_%d.map", "Export_test_rate_index_map_%d_%d_%d.map",
			"Export_cumul_incidence_by_site_%d_%d_%d.arr", "Export_cumul_incidence_by_person_%d_%d_%d.arr",
			"Export_cumul_treatment_by_person_%d_%d_%d.arr"

	};

	protected final int EXPORT_MAP_TRANS_PROB = EXPORT_SIMOUTPUT + 1;
	protected final int EXPORT_MAP_INF_STAGE = EXPORT_MAP_TRANS_PROB + 1;
	protected final int EXPORT_MAP_INF_SWICH = EXPORT_MAP_INF_STAGE + 1;
	protected final int EXPORT_SCHEDULE_TESTING = EXPORT_MAP_INF_SWICH + 1;
	protected final int EXPORT_SCHEDULE_STAGE_CHANGE = EXPORT_SCHEDULE_TESTING + 1;
	protected final int EXPORT_TEST_RATE_INDEX = EXPORT_SCHEDULE_STAGE_CHANGE + 1;
	protected final int EXPORT_CUMUL_INC_SITE = EXPORT_TEST_RATE_INDEX + 1;
	protected final int EXPORT_CUMUL_INC_PERSON = EXPORT_CUMUL_INC_SITE + 1;
	protected final int EXPORT_CUMUL_TREATMENT_PERSON = EXPORT_CUMUL_INC_PERSON + 1;

	protected final int[] EXPORT_INDEX_ARRAY = new int[] { EXPORT_MAP_TRANS_PROB, EXPORT_MAP_INF_STAGE,
			EXPORT_MAP_INF_SWICH, EXPORT_SCHEDULE_TESTING, EXPORT_SCHEDULE_STAGE_CHANGE, EXPORT_TEST_RATE_INDEX,
			EXPORT_CUMUL_INC_SITE, EXPORT_CUMUL_INC_PERSON, EXPORT_CUMUL_TREATMENT_PERSON };

	protected int[] exportTime = new int[0];

	@Override
	public void exportRunnableTransmission(int time_pt) throws IOException {
		super.exportRunnableTransmission(time_pt);
		File csvFile;
		PrintWriter pWri;

		// Zip file
		try {
			long tic = System.currentTimeMillis();

			for (int eI : EXPORT_INDEX_ARRAY) {
				csvFile = new File(baseDir, String.format(exportFileFormat[eI], cMAP_SEED, sIM_SEED, time_pt));
				pWri = new PrintWriter(csvFile);
				switch (eI) {
				case EXPORT_MAP_TRANS_PROB:
					for (Integer pid : map_trans_prob.keySet()) {
						double[][][][] ent = map_trans_prob.get(pid);
						pWri.print(pid);
						pWri.print(':');
						pWri.print(util.PropValUtils.objectToPropStr(ent, ent.getClass()));
						pWri.println();
					}
					break;
				case EXPORT_MAP_INF_STAGE:
				case EXPORT_MAP_INF_SWICH:
					HashMap<Integer, int[][]> map;
					switch (eI) {
					case EXPORT_MAP_INF_STAGE:
						map = map_currrent_infection_stage;
						break;
					case EXPORT_MAP_INF_SWICH:
						map = map_infection_stage_switch;
						break;
					default:
						System.err.printf("exportRunnableTransmission: export file index of %d not defined.\n", eI);
						map = null;
					}
					if (map != null) {
						for (Integer pid : map.keySet()) {
							int[][] ent = map.get(pid);
							pWri.print(pid);
							pWri.print(':');
							pWri.print(util.PropValUtils.objectToPropStr(ent, ent.getClass()));
							pWri.println();
						}
					}
					break;
				case EXPORT_SCHEDULE_TESTING:
					for (Integer time : schedule_testing.keySet()) {
						ArrayList<int[]> entArr = schedule_testing.get(time);
						pWri.print(time);
						for (int[] ent : entArr) {
							pWri.print(':');
							pWri.print(util.PropValUtils.objectToPropStr(ent, int[].class));
						}
						pWri.println();
					}
					break;
				case EXPORT_SCHEDULE_STAGE_CHANGE:
					for (Integer time : schedule_stage_change.keySet()) {
						ArrayList<ArrayList<ArrayList<Integer>>> ent = schedule_stage_change.get(time);
						pWri.print(time);
						for (int i = 0; i < NUM_INF; i++) {
							for (int s = 0; s < NUM_SITE; s++) {
								ArrayList<Integer> pids = ent.get(i).get(s);
								if (pids.size() > 0) {
									pWri.printf(":%d,%d", i, s);
									for (Integer pid : pids) {
										pWri.print(',');
										pWri.print(pid);
									}
								}
							}
						}
						pWri.println();
					}
					break;
				case EXPORT_TEST_RATE_INDEX:
					for (Integer pid : test_rate_index_map.keySet()) {
						HashMap<Integer, Integer> ent = test_rate_index_map.get(pid);
						pWri.print(pid);
						for (Entry<Integer, Integer> entry : ent.entrySet()) {
							pWri.printf(":%d,%d", entry.getKey(), entry.getValue());
						}
						pWri.println();
					}
					break;
				case EXPORT_CUMUL_INC_SITE:
					pWri.println(util.PropValUtils.objectToPropStr(cumul_incidence_by_site, int[][][].class));
					break;
				case EXPORT_CUMUL_INC_PERSON:
					pWri.println(util.PropValUtils.objectToPropStr(cumul_incidence_by_person, int[][].class));
					break;
				case EXPORT_CUMUL_TREATMENT_PERSON:
					pWri.println(util.PropValUtils.objectToPropStr(cumul_treatment_by_person, int[][].class));
					break;
				default:
					System.err.printf("exportRunnableTransmission: export file index of %d not defined.\n", eI);
				}
				pWri.close();

				File archiveFile = new File(baseDir, csvFile.getName() + ".7z");
				util.Util_7Z_CSV_Entry_Extract_Callable.zipFile(new File[] { csvFile }, archiveFile, false);
				Files.delete(csvFile.toPath());
			}

			// Delete old zip only when all zip okay
			for (String srcFormat : exportFileFormat) {
				File archiveFile = new File(baseDir, String.format(srcFormat, cMAP_SEED, sIM_SEED, time_pt) + ".7z");
				String pat_format = srcFormat.replaceFirst("%d", Long.toString(cMAP_SEED));
				pat_format = pat_format.replaceFirst("%d", Long.toString(sIM_SEED));
				pat_format = pat_format.replaceFirst("%d", "(\\\\d+)");
				Pattern pattern_previous_zips = Pattern.compile(pat_format + ".7z");
				File[] prev_zips = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pattern_previous_zips.matcher(pathname.getName()).matches()
								&& !pathname.equals(archiveFile);
					}
				});
				for (File f : prev_zips) {
					Files.delete(f.toPath());
				}
			}

			if (print_progress != null && runnableId != null) {
				try {
					print_progress.printf("Thread <%s>: Population stat exported at T = %d . Time req. = %.3fs\n",
							runnableId, time_pt, (System.currentTimeMillis() - tic) / 1000.0);
				} catch (Exception e) {
					e.printStackTrace(System.err);
				}
			}

		} catch (IOException e) {
			e.printStackTrace(System.err);
		}

	}

	@Override
	public void importRunnableTransmission(int time_pt) throws IOException, ClassNotFoundException {
		super.importRunnableTransmission(time_pt);
		File csvFile;
		BufferedReader reader;
		String line;

		for (int eI : EXPORT_INDEX_ARRAY) {
			long tic = System.currentTimeMillis();
			csvFile = new File(baseDir, String.format(exportFileFormat[eI], cMAP_SEED, sIM_SEED, time_pt));
			reader = new BufferedReader(new FileReader(csvFile));
			while ((line = reader.readLine()) != null) {
				String[] lineSp = line.split(":");
				int key = -1;
				if (lineSp.length > 1) {
					key = Integer.parseInt(lineSp[0]);
				}
				switch (eI) {
				case EXPORT_MAP_TRANS_PROB:
					map_trans_prob.put(key,
							(double[][][][]) util.PropValUtils.propStrToObject(lineSp[1], double[][][][].class));
					break;
				case EXPORT_MAP_INF_STAGE:
					map_currrent_infection_stage.put(key,
							(int[][]) util.PropValUtils.propStrToObject(lineSp[1], int[][].class));
					break;
				case EXPORT_MAP_INF_SWICH:
					map_infection_stage_switch.put(key,
							(int[][]) util.PropValUtils.propStrToObject(lineSp[1], int[][].class));
					break;
				case EXPORT_SCHEDULE_TESTING:
					ArrayList<int[]> ent = new ArrayList<>();
					for (int i = 1; i < lineSp.length; i++) {
						ent.add((int[]) util.PropValUtils.propStrToObject(lineSp[i], int[].class));
					}
					schedule_testing.put(key, ent);
					break;
				case EXPORT_SCHEDULE_STAGE_CHANGE:
					ArrayList<ArrayList<ArrayList<Integer>>> schedule_inf = new ArrayList<>(NUM_INF);
					for (int i = 0; i < NUM_INF; i++) {
						ArrayList<ArrayList<Integer>> schedule_inf_site = new ArrayList<>(NUM_SITE);
						for (int s = 0; s < NUM_SITE; s++) {
							schedule_inf_site.add(new ArrayList<>());
						}
						schedule_inf.add(schedule_inf_site);
					}
					for (int i = 1; i < lineSp.length; i++) {
						String[] intEnt = lineSp[i].split(",");
						if (intEnt.length > 2) {
							int inf_id = Integer.parseInt(intEnt[0]);
							int site_id = Integer.parseInt(intEnt[1]);
							for (int pI = 2; pI < intEnt.length; pI++) {
								schedule_inf.get(inf_id).get(site_id).add(Integer.parseInt(intEnt[pI]));
							}
						}
					}
					schedule_stage_change.put(key, schedule_inf);
					lastStateSwitch = Math.max(lastStateSwitch, key);
					break;
				case EXPORT_TEST_RATE_INDEX:
					HashMap<Integer, Integer> entry = new HashMap<>();
					for (int i = 1; i < lineSp.length; i++) {
						String[] intEnt = lineSp[i].split(",");
						if (intEnt.length > 2) {
							entry.put(Integer.parseInt(intEnt[0]), Integer.parseInt(intEnt[0]));
						}
					}
					test_rate_index_map.put(key, entry);
					break;
				case EXPORT_CUMUL_INC_SITE:
					cumul_incidence_by_site = (int[][][]) util.PropValUtils.propStrToObject(line, int[][][].class);
					break;
				case EXPORT_CUMUL_INC_PERSON:
					cumul_incidence_by_person = (int[][]) util.PropValUtils.propStrToObject(line, int[][].class);
					break;
				case EXPORT_CUMUL_TREATMENT_PERSON:
					cumul_treatment_by_person = (int[][]) util.PropValUtils.propStrToObject(line, int[][].class);
					break;
				default:
					System.err.printf("importRunnableTransmission: file index of %d not defined.\n", eI);
				}
			}
			reader.close();
			
			if (print_progress != null && runnableId != null) {
				print_progress.printf("Thread <%s>: Import of %s completed. Time = %.3fs\n", runnableId, 
						csvFile.getName(), (System.currentTimeMillis() - tic)/1000.0) ;
			}
		}

		// Fill others

		// map_currently_infectious
		if (map_currently_infectious == null) {
			map_currently_infectious = new HashMap<>();
		}
		for (int pid : map_currrent_infection_stage.keySet()) {
			int[][] ent = map_currrent_infection_stage.get(pid);
			for (int infection_id = 0; infection_id < NUM_INF; infection_id++) {
				for (int site_id = 0; site_id < NUM_SITE; site_id++) {
					int current_infection_stage = ent[infection_id][site_id];
					if (((1 << current_infection_stage)
							& lookupTable_infection_infectious_stages[infection_id][site_id]) != 0) {
						String key = String.format("%d,%d", infection_id, site_id);
						ArrayList<Integer> currenty_infectious_ent = map_currently_infectious.get(key);
						if (currenty_infectious_ent == null) {
							currenty_infectious_ent = new ArrayList<>();
							map_currently_infectious.put(key, currenty_infectious_ent);
						}
						int pt = Collections.binarySearch(currenty_infectious_ent, pid);
						if (pt < 0) {
							currenty_infectious_ent.add(~pt, pid);
						}
					}

				}
			}
		}
		firstSeedTime = time_pt+1; // Next time step
	}

	@Override
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

			random.RandomGenerator trans_RNG = RNG;

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
									trans_RNG, param);
						} else {
							// Beta distribution
							dist_tranmissionMatrix[inf_id][site_from][site_to][statePt] = generateBetaDistribution(
									trans_RNG, param);
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
				lookupTable_instant_stage_switch.clear();
			}

			HashMap<String, HashMap<Integer, Double>> map_infection_stage_switch; // Key="INF_ID,SITE_ID,STAGE_ID",V=Map(Stage->Prob)
			HashMap<Integer, Double> map_target_stage;
			map_infection_stage_switch = new HashMap<>();

			random.RandomGenerator stage_dist = RNG;

			for (double[] ent : field) {
				int inf_id = (int) ent[FIELD_STAGE_PERIOD_INF_ID];
				int site = (int) ent[FIELD_STAGE_PERIOD_SITE];
				int stage_id_to = (int) ent[FIELD_STAGE_PERIOD_STATE_TO];

				if (stage_id_to < 0) {
					// Instantaneous stage switch
					String key = String.format("%d,%d,%d", inf_id, site, stage_id_to);
					ArrayList<double[]> arrEnt = lookupTable_instant_stage_switch.get(key);
					if (arrEnt == null) {
						arrEnt = new ArrayList<>();
						lookupTable_instant_stage_switch.put(key, arrEnt);
					}
					arrEnt.add(ent);
				} else {
					int stage_id_from = (int) ent[FIELD_STAGE_PERIOD_STATE_FROM];
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
						dist_stage_period[inf_id][site] = Arrays.copyOf(dist_stage_period[inf_id][site],
								stage_id_to + 1);
					}
					if (param.length == 1 || Double.isNaN(param[1])) {
						dist_stage_period[inf_id][site][stage_id_to] = generateNonDistribution(stage_dist, param);
					} else {
						if (param[1] < 0) {
							// Uniform distribution
							double[] param_u = Arrays.copyOf(param, param.length);
							param_u[1] = Math.abs(param_u[1]);
							dist_stage_period[inf_id][site][stage_id_to] = generateUniformDistribution(stage_dist,
									param_u);
						} else {
							// Gamma Distribution
							dist_stage_period[inf_id][site][stage_id_to] = generateGammaDistribution(stage_dist, param);
						}
					}
				}
			}
			for (String key : map_infection_stage_switch.keySet()) {
				map_target_stage = map_infection_stage_switch.get(key);
				double[] cumul_prob = new double[map_target_stage.size() * 2];
				Integer[] target_stages = map_target_stage.keySet().toArray(new Integer[map_target_stage.size()]);
				Arrays.sort(target_stages);
				double prob_sum = 0;
				for (Integer target_stage : target_stages) {
					prob_sum += map_target_stage.get(target_stage).doubleValue();
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
			random.RandomGenerator sym_rate_rng = RNG;
			for (double[] ent : field) {
				int inf_Id = (int) ent[FIELD_SYM_RATE_INF_ID];
				int site = (int) ent[FIELD_SYM_RATE_SITE];
				int state = (int) ent[FIELD_SYM_RATE_STATE_ID];
				double[] param = Arrays.copyOfRange(ent, FIELD_SYM_RATE_PROB, ent.length);

				if (dist_sym_rate[inf_Id][site] == null) {
					dist_sym_rate[inf_Id][site] = new RealDistribution[state + 1];
				} else if (!(state < dist_sym_rate[inf_Id][site].length)) {
					dist_sym_rate[inf_Id][site] = Arrays.copyOf(dist_sym_rate[inf_Id][site], state + 1);
				}
				if (param.length == 1) {
					dist_sym_rate[inf_Id][site][state] = generateNonDistribution(sym_rate_rng, param);
				} else {
					dist_sym_rate[inf_Id][site][state] = generateBetaDistribution(sym_rate_rng, param);
				}
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES:
			if (clearAll) {
				lookupTable_test_treatment_properties.clear();
			}
			for (double[] ent : field) {
				int inf_Id = (int) ent[FIELD_DX_TEST_PROPERTIES_INF_ID];
				int site = (int) ent[FIELD_DX_TEST_PROPERTIES_SITE];
				String key = String.format("%d,%d", inf_Id, site);
				lookupTable_test_treatment_properties.put(key, ent);
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM:

			random.RandomGenerator test_period_rng = RNG;
			for (double[] ent : field) {
				int gender_inc = (int) ent[FIELD_SOUGHT_TEST_PERIOD_BY_SYM_GENDER_INC];
				for (int g = 0; g < dist_seek_test_period.length; g++) {
					if ((gender_inc & 1 << g) != 0) {
						dist_seek_test_period[g] = generateGammaDistribution(test_period_rng,
								new double[] { ent[FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SEEK_TEST_MEAN],
										ent[FIELD_SOUGHT_TEST_PERIOD_BY_SYM_SEEK_TEST_SD] });
						if (ent.length > FIELD_SOUGHT_TEST_RATE) {
							seek_test_rate[g] = ent[FIELD_SOUGHT_TEST_RATE];
						}
					}
				}
			}
			break;
		case RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS:
		case RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES:

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

		sim_output = new HashMap<>();

		if (test_rate_index_map == null) {
			test_rate_index_map = new HashMap<>();
		}

		// Check export time
		if (baseProp.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_POP_EXPORT_AT])) {
			exportTime = (int[]) util.PropValUtils.propStrToObject(
					baseProp.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_POP_EXPORT_AT]),
					int[].class);
		}
		// Import previous exported entries
		importExportedTransmissionStates(exportFileFormat);

		for (int r = 0; r < runnable_fields.length; r++) {
			refreshField(r, true);
		}

		// Non exported
		if (risk_cat_map == null) {
			risk_cat_map = new HashMap<>();
		}
		if (propSwitch_map == null) {
			propSwitch_map = new HashMap<>();
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

		if (importedAtTime < 0) {
			// Schedule testing if not already set
			if (bASE_CONTACT_MAP != null) {
				for (Integer personId : bASE_CONTACT_MAP.vertexSet()) {
					scheduleNextTest(personId, startTime);
				}
			}
			if (pop_stat != null) {
				// Schedule test for pop
				for (Integer pid : pop_stat.keySet()) {
					String[] popEnt = pop_stat.get(pid);
					scheduleNextTest(pid, Integer
							.parseInt(popEnt[Abstract_Runnable_ClusterModel_Transmission.POP_INDEX_ENTER_POP_AT]));
				}
			}
			// End of schedule test

			cumul_incidence_by_site = new int[NUM_INF][NUM_GENDER][NUM_SITE];
			cumul_incidence_by_person = new int[NUM_INF][NUM_GENDER];
			cumul_treatment_by_person = new int[NUM_INF][NUM_GENDER];
		}

		int snap_index = 0;
		boolean hasInfectious = hasInfectious();

		HashMap<String, int[]> acted_today = new HashMap<>();
		// K="PID,PID", V=int[act] =
		// if SITE_P1 != SITE_P2
		// = 0 if no acts
		// = 1 if p1s1 -> p2s2, e.g. p1 insertive, or if the site is the same
		// = 2 if p2s1 -> p1s2 e.g. p2 insertive

		// For timing stat
		long tic = System.currentTimeMillis();
		long debug_time_test = 0, debug_time_inf = 0, debug_time_cMap = 0;
		int debug_num_test = 0, debug_num_edge_add = 0, debug_num_new_inf = 0;

		for (int currentTime = startTime; currentTime < startTime + nUM_TIME_STEPS_PER_SNAP * nUM_SNAP
				&& (currentTime < lastStateSwitch || hasInfectious); currentTime++) {

			if (switchTimeIndex < switchTime.length && switchTime[switchTimeIndex] == currentTime) {
				HashMap<Integer, String> switch_ent = propSwitch_map.get(currentTime);
				for (Integer switch_index : switch_ent.keySet()) {
					String str_obj = switch_ent.get(switch_index);

					int fieldId = switch_index - RUNNABLE_OFFSET;
					if (fieldId < 0) {
						loadNonRunnableFieldSetting(switch_index, str_obj, currentTime);
					} else {
						getRunnable_fields()[fieldId] = PropValUtils.propStrToObject(str_obj,
								getRunnable_fields()[fieldId].getClass());
						refreshField(fieldId, false);
					}
				}
				switchTimeIndex++;
			}

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

			ArrayList<Integer> currenty_infectious_ent = new ArrayList<Integer>();
			for (String key : currently_infectious_keys) {
				currenty_infectious_ent.addAll(map_currently_infectious.get(key));
			}
			Collections.sort(currenty_infectious_ent);

			// Update cMap
			tic = System.currentTimeMillis();
			int[] updatecMap_Stat = updateCMap(cMap, currentTime, edges_array, edges_array_pt, removeEdges,
					currenty_infectious_ent);
			edges_array_pt = updatecMap_Stat[0];
			debug_num_edge_add += updatecMap_Stat[1];
			debug_time_cMap += System.currentTimeMillis() - tic;

			acted_today.clear();

			tic = System.currentTimeMillis();

			for (String key : currently_infectious_keys) {
				currenty_infectious_ent = map_currently_infectious.get(key);

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
												hasActed[a] = 1;
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

											if (isValidInfectionTargetSite(inf_id, tar_site, tar_infection_stages)) {

												double transProb = getTransmissionProb(currentTime, inf_id, pid_inf_src,
														pid_inf_tar,
														e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION], a,
														src_site, tar_site);

												if (transProb > 0) {
													if (RNG.nextDouble() < transProb) {
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
																newInfection &= isValidInfectionTargetSite(inf_id,
																		siteId, tar_infection_stages);
																// tar_infection_stages[inf_id][siteId] ==
																// AbstractIndividualInterface.INFECT_S;
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
														simulate_transmission_success_act(currentTime, inf_id,
																pid_inf_src, pid_inf_tar, src_site, tar_site);
														debug_num_new_inf++;
													} else {
														simulate_transmission_failed_act(currentTime, inf_id,
																pid_inf_src, pid_inf_tar, src_site, tar_site);

													} // End of determining transmission
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

			debug_time_inf += System.currentTimeMillis() - tic;

			// Simulate acts among non-infectious
			simulate_non_infectious_act(currentTime, cMap, acted_today);

			// Testing
			tic = System.currentTimeMillis();
			ArrayList<int[]> testToday = schedule_testing.remove(currentTime);
			if (testToday != null) {
				for (int[] testing_stat : testToday) {
					// Testing stat: int[]{PID,INF_INCLUDE_INDEX,SITE_INCLUDE_INDEX}
					// symptomatic or one off test if PID < 0
					int pid = Math.abs(testing_stat[0]);
					int infIncl = testing_stat[1];
					int siteIncl = testing_stat[2];
					testPerson(currentTime, pid, infIncl, siteIncl, cumul_treatment_by_person);
					// Schedule next test
					if (testing_stat[0] > 0) {
						scheduleNextTest(pid, currentTime, infIncl, siteIncl);
					}

				} // End of testing for single individual
			} // End of all scheduled testing for today

			debug_time_test += System.currentTimeMillis() - tic;
			if (testToday != null) {
				debug_num_test += testToday.size();
			}

			// Storing of outputs
			if (snap_index == 0) {

				String key;

				if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE) != 0) {
					key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_GENDER_COUNT,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> current_infectious_gender_count = (HashMap<Integer, int[]>) sim_output
							.get(key);
					if (current_infectious_gender_count == null) {
						current_infectious_gender_count = new HashMap<>();
						sim_output.put(key, current_infectious_gender_count);
					}

					key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_SITE_COUNT,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> current_infectious_site_count = (HashMap<Integer, int[]>) sim_output
							.get(key);
					if (current_infectious_site_count == null) {
						current_infectious_site_count = new HashMap<>();
						sim_output.put(key, current_infectious_site_count);
					}

					key = String.format(SIM_OUTPUT_KEY_INFECTED_AT_GENDER_COUNT,
							Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
					@SuppressWarnings("unchecked")
					HashMap<Integer, int[]> current_infected_count = (HashMap<Integer, int[]>) sim_output.get(key);
					if (current_infected_count == null) {
						current_infected_count = new HashMap<>();
						sim_output.put(key, current_infected_count);
					}

					int[] num_infectious_site_count = new int[NUM_INF * NUM_SITE];
					int[] num_infectious_person_count = new int[NUM_INF * NUM_GENDER];
					int[] num_infected_by_gender_site_at = new int[NUM_INF * NUM_GENDER * (1 << (NUM_SITE + 1))];

					current_infectious_gender_count.put(currentTime, num_infectious_person_count);
					current_infectious_site_count.put(currentTime, num_infectious_site_count);
					current_infected_count.put(currentTime, num_infected_by_gender_site_at);

					for (Integer pid : getCurrentPopulationPId(currentTime)) {
						int[][] current_stage_arr = map_currrent_infection_stage.get(pid);
						int gI = getGenderType(pid);
						for (int i = 0; i < NUM_INF; i++) {
							int infected_site = 0;
							int infectious_site = 0;

							for (int s = 0; s < NUM_SITE; s++) {
								int infection_stage = lookupTable_infection_infectious_stages[i][s];
								if (current_stage_arr != null) {
									if (current_stage_arr[i][s] != AbstractIndividualInterface.INFECT_S) {
										if (current_stage_arr[i][s] >= 0
												&& ((1 << current_stage_arr[i][s]) & infection_stage) != 0) {
											// Infectious
											infectious_site |= 1 << s;
											num_infectious_site_count[i * NUM_SITE + s]++;
										}

										infected_site |= 1 << s;
									}
								}
							}

							if (infectious_site != 0) {
								num_infectious_person_count[i * NUM_GENDER + getGenderType(pid)]++;
							}

							int pt_infected_site_at = i * NUM_GENDER * (1 << (NUM_SITE + 1))
									+ gI * (1 << (NUM_SITE + 1)) + infected_site;
							num_infected_by_gender_site_at[pt_infected_site_at]++;
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
								if (inf_stage[i][s] != AbstractIndividualInterface.INFECT_S) {
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

						print_progress.printf(
								"T=%d, CMAP={%d,%d}" + ", # edge_added = %d, time_edge_added = %.3f"
										+ ", # inf added = %d,  time inf = %.3f" + ", # test = %d, time_test = %.3fs\n",
								currentTime, cMap.vertexSet().size(), cMap.edgeSet().size(), debug_num_edge_add,
								debug_time_cMap / 1000f, debug_num_new_inf, debug_time_inf / 1000f, debug_num_test,
								debug_time_test / 1000f);

					} catch (Exception ex) {
						System.err.printf("Thread <%s>: t = %d .\n", runnableId, currentTime);
					}

					debug_time_test = 0;
					debug_time_inf = 0;
					debug_time_cMap = 0;
					debug_num_test = 0;
					debug_num_edge_add = 0;
					debug_num_new_inf = 0;
				}
			}

			snap_index = (snap_index + 1) % nUM_TIME_STEPS_PER_SNAP;
			hasInfectious = hasInfectious();

			postTimeStep(currentTime);

		} // End of time step

		if (runnableId != null) {
			System.out.printf("Thread <%s> completed.\n", runnableId);
		}

		// Post simulation
		postSimulation();

	}

	protected void simulate_transmission_failed_act(int currentTime, int inf_id, Integer pid_inf_src, int pid_inf_tar,
			int src_site, int tar_site) {
		// Do nothing by default

	}

	protected void simulate_transmission_success_act(int currentTime, int inf_id, Integer pid_inf_src, int pid_inf_tar,
			int src_site, int tar_site) {
		// Do nothing by default

	}

	protected boolean isValidInfectionTargetSite(int inf_id, int tar_site, int[][] tar_infection_stages) {
		return tar_infection_stages == null
				|| tar_infection_stages[inf_id][tar_site] == AbstractIndividualInterface.INFECT_S;
	}

	protected void testPerson(int currentTime, int pid, int infIncl, int siteIncl, int[][] cumul_treatment_by_person) {
		for (int infId = 0; infId < NUM_INF; infId++) {
			if ((infIncl & 1 << infId) != 0) {
				boolean applyTreatment = false;
				double[] test_properies;
				int tested_stage_inc;
				int[][] inf_stage = null;

				for (int siteId = 0; siteId < NUM_SITE && !applyTreatment; siteId++) {
					if ((siteIncl & 1 << siteId) != 0) {
						// Test for the site
						test_properies = lookupTable_test_treatment_properties
								.get(String.format("%d,%d", infId, siteId));
						if (test_properies != null) {
							inf_stage = map_currrent_infection_stage.get(pid);
							if (inf_stage != null && inf_stage[infId][siteId] >= 0) {
								double testSensitivity = 0;
								int stage_pt = FIELD_DX_TEST_PROPERTIES_ACCURACY_START;

								while (testSensitivity == 0 && stage_pt < test_properies.length) {
									// TEST_ACCURACY_1, TARGET_STAGE_INC_1, TREATMENT_SUC_STAGE_1 ..
									tested_stage_inc = (int) test_properies[stage_pt + 1];
									if ((tested_stage_inc & 1 << inf_stage[infId][siteId]) != 0) {
										testSensitivity = test_properies[stage_pt];
									}
									stage_pt += 3;
								}
								if (testSensitivity > 0) {
									applyTreatment |= RNG.nextDouble() < testSensitivity;
								}
							}
						}
					} // End of testing pid for infection and site
				}
				if (applyTreatment) {
					cumul_treatment_by_person[infId][getGenderType(pid)]++;
					applyTreatment(currentTime, infId, pid, inf_stage);
				}
			} // End of testing for infection
		} // Check if infection is tested
	}

	protected void postTimeStep(int currentTime) {
		if (exportTime.length > 0) {
			if (currentTime > 0 && currentTime % exportTime[0] == 0) {
				try {
					exportRunnableTransmission(currentTime);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	protected double getTransmissionProb(int currentTime, int inf_id, int pid_inf_src, int pid_inf_tar,
			int partnershiptDur, int actType, int src_site, int tar_site) {

		double[] actFieldEntry = table_act_frequency[actType][getGenderType(pid_inf_src)][getGenderType(pid_inf_tar)];

		int[][] src_infection_stages = map_currrent_infection_stage.get(pid_inf_src);
		int src_stage = src_infection_stages[inf_id][src_site];
		double[][][][] trans_prob = map_trans_prob.get(pid_inf_src);

		double prob = trans_prob[inf_id][src_site][tar_site][src_stage];

		if (FIELD_ACT_FREQ_CONDOM_EFFICACY > actFieldEntry.length) {
			double condomUsage = (partnershiptDur > 1 ? actFieldEntry[FIELD_ACT_FREQ_USAGE_REG]
					: actFieldEntry[FIELD_ACT_FREQ_USAGE_CASUAL]);
			if (RNG.nextDouble() < condomUsage) {
				prob *= (1 - actFieldEntry[FIELD_ACT_FREQ_CONDOM_EFFICACY]);
			}
		}

		return prob;
	}

	protected void simulate_non_infectious_act(int currentTime, ContactMap cMaop, HashMap<String, int[]> acted_today) {
		// Do noting on its own, but can be overwritten by subclasses
	}

	protected void applyTreatment(int currentTime, int infId, int pid, int[][] inf_stage) {
		double[] test_treatment_properies;
		int tested_stage_inc;
		int[][] infection_switch = map_infection_stage_switch.get(pid);

		for (int siteId = 0; siteId < NUM_SITE; siteId++) {
			test_treatment_properies = lookupTable_test_treatment_properties.get(String.format("%d,%d", infId, siteId));

			if (test_treatment_properies != null) {

				// TEST_ACCURACY_1, TARGET_STAGE_INC_1, TREATMENT_SUC_STAGE_1 ...
				int stage_inc_pt = FIELD_DX_TEST_PROPERTIES_ACCURACY_START + 1;

				while (stage_inc_pt < test_treatment_properies.length) {
					tested_stage_inc = (int) test_treatment_properies[stage_inc_pt];
					if ((tested_stage_inc & 1 << inf_stage[infId][siteId]) != 0) {
						int nextStage = (int) test_treatment_properies[stage_inc_pt + 1];
						inf_stage[infId][siteId] = nextStage;
						updateInfectStageChangeSchedule(pid, infId, siteId, currentTime + 1,
								infection_switch[infId][siteId]);
					}
					stage_inc_pt += 3;
				}
			}
		}
	}

	@Override
	protected void postSimulation() {
		super.postSimulation();
		// System.out.printf("Warning. %s.postSimulation() not set or overwritten.\n",
		// this.getClass().getName());
	}

	public void printCountMap(HashMap<Integer, int[]> countMap, String fileName, String headerFormat, int[] dimension) {
		printCountMap(countMap, fileName, headerFormat, dimension, null);
	}

	public void printCountMap(HashMap<Integer, int[]> countMap, String fileName, String headerFormat, int[] dimension,
			int[] col_to_print) {

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
		StringBuilder header = new StringBuilder();
		recursiveHeaderGeneration(dimension, 0, headerFormat, header);
		if (col_to_print != null) {
			StringBuilder sel_header = new StringBuilder();
			String[] ent = header.toString().split(",");
			for (int col_sel : col_to_print) { // Offset of 1
				if (col_sel + 1 < ent.length) {
					sel_header.append(',');
					sel_header.append(ent[col_sel + 1]);
				}
			}

			header = sel_header;
		}

		pWri.println(header.toString());

		// Entry
		Integer[] timePt = countMap.keySet().toArray(new Integer[countMap.size()]);
		Arrays.sort(timePt);
		for (Integer time : timePt) {
			int[] ent = countMap.get(time);
			pWri.print(time);
			for (int i = 0; i < ent.length; i++) {
				if (col_to_print == null || Arrays.binarySearch(col_to_print, i) >= 0) {
					pWri.print(',');
					pWri.print(ent[i]);
				}
			}
			pWri.println();
		}

		pWri.close();

		if (s != null) {
			System.out.println(s.toString());
		}

	}

	private void recursiveHeaderGeneration(int[] dimension, int dPt, String seedStr, StringBuilder strBuilder) {
		for (int i = 0; i < dimension[dPt]; i++) {
			if (dPt == dimension.length - 1) {
				strBuilder.append(',');
				strBuilder.append(seedStr.replaceFirst("%d", Integer.toString(i)));
			} else {
				recursiveHeaderGeneration(dimension, dPt + 1, seedStr.replaceFirst("%d", Integer.toString(i)),
						strBuilder);
			}
		}

	}

	public void test_rate_index_map_reset() {
		test_rate_index_map.clear();
	}

	public void scheduleNextTest(Integer personId, int lastTestTime) {
		scheduleNextTest(personId, lastTestTime, -1, -1);
	}

	public void scheduleNextTest(Integer personId, int lastTestTime, int last_test_infIncl, int last_test_siteIncl) {
		int riskCat;
		if (risk_cat_map.containsKey(personId)) {
			riskCat = risk_cat_map.get(personId);
		} else {
			if (bASE_CONTACT_MAP != null) {
				riskCat = determineRiskGrpByCasualPartners(personId);
			} else {
				riskCat = 0;
			}
			risk_cat_map.put(personId, riskCat);
		}

		if (riskCat != -1) {
			double[][] testRateDefs = (double[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES];
			int genderType = getGenderType(personId);

			HashMap<Integer, Integer> past_test_pt = test_rate_index_map.get(personId);

			for (int i = 0; i < testRateDefs.length; i++) {
				double[] testRateDef = testRateDefs[i];
				int gIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_GENDER_INCLUDE_INDEX];
				int iIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_INF_INCLUDE_INDEX];
				int rIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_RISK_GRP_INCLUDE_INDEX];
				int sIncl = (int) testRateDef[FIELD_TESTING_RATE_BY_RISK_CATEGORIES_SITE_INCLUDE_INDEX];

				if (((gIncl & (1 << genderType)) != 0) && ((rIncl & (1 << riskCat)) != 0)
						&& (last_test_infIncl < 0 || last_test_infIncl == iIncl)
						&& (last_test_siteIncl < 0 || last_test_siteIncl == sIncl)) {

					int test_pt;

					if (past_test_pt == null || past_test_pt.get(i) == null) {
						int pt;

						// FORMAT: {Cumul_Prob_0, Cummul_Prob_1, .... test_gap_time_0, test_gap_time_1,
						// test_gap_time_2..}
						int cumul_end_index = FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START
								+ ((testRateDef.length - FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START)
										- 1) / 2;

						if (past_test_pt == null) {
							past_test_pt = new HashMap<>();
							test_rate_index_map.put(personId, past_test_pt);
						}
						pt = Arrays.binarySearch(testRateDef,
								FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START, cumul_end_index,
								RNG.nextDouble());
						if (pt < 0) {
							pt = ~pt;
						}

						test_pt = cumul_end_index - FIELD_TESTING_RATE_BY_RISK_CATEGORIES_TEST_RATE_PARAM_START + pt;
						past_test_pt.put(i, test_pt);
					} else {
						test_pt = past_test_pt.get(i);
					}

					if (!Double.isInfinite(testRateDef[test_pt])) {
						int nextTestAfter = (int) (testRateDef[test_pt + 1]
								+ RNG.nextInt((int) testRateDef[test_pt] - (int) testRateDef[test_pt + 1]));

						int nextTestDate = lastTestTime + nextTestAfter;

						if (nextTestDate < exitPopAt(personId)) {
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
		}

	}

	protected int determineRiskGrpByCasualPartners(Integer personId) {
		int riskCat;
		riskCat = -1;
		int genderType = getGenderType(personId);
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
		double numCasualPerYear;

		if (lastPartnerTime > firstPartnerTime) {
			numCasualPerYear = (((double) AbstractIndividualInterface.ONE_YEAR_INT) * numCasual)
					/ (lastPartnerTime - firstPartnerTime);
		} else {
			numCasualPerYear = (((double) AbstractIndividualInterface.ONE_YEAR_INT) * numCasual)
					/ (nUM_TIME_STEPS_PER_SNAP * nUM_SNAP);
		}

		double[][] riskCatDefs = (double[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS];

		for (int i = 0; i < riskCatDefs.length && riskCat == -1; i++) {
			int gIncl = (int) riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_GENDER_INCLUDE_INDEX];
			if ((gIncl & 1 << genderType) != 0) {
				if (riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_NUM_CASUAL_PARTNER_LOWER] <= numCasualPerYear
						&& numCasualPerYear < riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_NUM_CASUAL_PARTNER_UPPER]) {
					riskCat = (int) riskCatDefs[i][FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS_RISK_GRP_DEF_ID];
				}
			}

		}
		return riskCat;
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
		// or
		// num_infected = int[infection_id]{-MIN_NUM_PARTNER_IN_LAST_12_MONTHS,
		// SITE_INDEX_0,
		// Number_INF_0,...}
		firstSeedTime = Math.min(firstSeedTime, time);

		int inf_id = 0;
		ArrayList<Integer> candidate = new ArrayList<>();

		ContactMap cMap = null;

		cMapLoop: for (int[] inf_setting : num_infected) {
			int pt = 0;
			while (pt < inf_setting.length) {
				if (inf_setting[pt] < 0) {
					if (bASE_CONTACT_MAP == null) {
						System.err.println("allocateSeedInfection: MIN_NUM_PARTNER_IN_LAST_12_MONTHS setting "
								+ "not yet support for null base contact map in this version");

					} else {
						cMap = Abstract_Runnable_ClusterModel.generateContactMapAcrossTimeRange(bASE_CONTACT_MAP,
								new int[] { time - AbstractIndividualInterface.ONE_YEAR_INT, time });
					}

					break cMapLoop;
				}
				pt += 3;
			}
		}

		Integer[] all_vertex = getCurrentPopulationPId(time); // bASE_CONTACT_MAP.vertexSet().toArray(new
																// Integer[bASE_CONTACT_MAP.vertexSet().size()]);
		Arrays.sort(all_vertex);

		for (int[] inf_setting : num_infected) {
			int pt = 0;
			while (pt < inf_setting.length) {
				int includeIndex = inf_setting[pt];
				pt++;
				int site_index = inf_setting[pt];
				pt++;
				int num_inf = inf_setting[pt];
				pt++;

				int validStage = 1; // lookupTable_infection_infectious_stages[inf_id][site_index];
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
					for (Integer pid : all_vertex) {
						int[][] current_stage_arr = map_currrent_infection_stage.get(pid);

						if (isValidInfectionTargetSite(inf_id, site_index, current_stage_arr)) {
							if (includeIndex > 0) { // By gender group
								if ((includeIndex & 1 << getGenderType(pid)) > 0) {
									candidate.add(pid);
								}
							} else { // By Risk Group
								int numPartLimit = -includeIndex;
								if (cMap != null && cMap.containsVertex(pid)) {
									// pWri.printf("%d,%d\n", pid, cMap.degreeOf(pid));
									if (cMap.degreeOf(pid) >= numPartLimit) {
										candidate.add(pid);
									}
								}
							}
						}
					}

					int counter = 0;
					int num_inf_added = 0;

					for (Integer pid : candidate) {
						if (RNG.nextInt(candidate.size() - counter) < (num_inf - num_inf_added)) {
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
							num_inf_added++;

						}
						counter++;
					}

					if (print_progress != null && num_inf_added > 0) {
						print_progress.printf("Time = %d. Inf #%d candidate = %d. # added = %d.\n", time, inf_id,
								candidate.size(), num_inf_added);
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

		int current_infect_switch_time = updateInfectionStage(infectedPId, infectionId, site_id, stage_id,
				infectious_time, current_stage_arr, infection_state_switch, state_duration_adj);

		return current_infect_switch_time;
	}

	// If state_duration_adj > 0, simple duration update, otherwise switch to next
	// stage and resample
	protected int updateInfectionStage(Integer pid, int infection_id, int site_id, int current_infection_stage,
			int current_time, int[][] current_stage_arr, int[][] infection_state_switch, int state_duration_preset) {
		int pt;
		double state_duration;
		int infect_switch_time = current_time;

		if (current_infection_stage == STAGE_ID_JUST_INFECTED) {
			String key = String.format("%d,%d,%d", infection_id, site_id, STAGE_ID_JUST_INFECTED);
			ArrayList<double[]> optionsArr = lookupTable_instant_stage_switch.get(key);
			if (optionsArr != null) {
				int riskGrp = risk_cat_map.get(pid);
				for (double[] options : optionsArr) {
					if (current_infection_stage != STAGE_ID_JUST_INFECTED) {
						break;
					}
					int riskGrpInc = (int) options[FIELD_STAGE_PERIOD_INSTANT_RISKGRP_INC];
					if ((riskGrpInc & 1 << riskGrp) != 0) {
						current_infection_stage = (int) options[FIELD_STAGE_PERIOD_INSTANT_STAGE_OPTIONS_START];
						if (FIELD_STAGE_PERIOD_INSTANT_STAGE_OPTIONS_START + 1 < options.length) {
							// More than one option
							int numPossibleStage = (options.length - FIELD_STAGE_PERIOD_INSTANT_STAGE_OPTIONS_START)
									/ 2;
							int opPt = Arrays.binarySearch(options,
									FIELD_STAGE_PERIOD_INSTANT_STAGE_OPTIONS_START + numPossibleStage, options.length,
									RNG.nextDouble());
							if (opPt < 0) {
								opPt = ~opPt;
							}
							current_infection_stage = (int) options[opPt - numPossibleStage];
						}

					}
				}
			}
			if (current_infection_stage == STAGE_ID_JUST_INFECTED) {
				// Default to first stage
				current_infection_stage = 0;
			}

			infect_switch_time = addInfectious(pid, infection_id, site_id, current_infection_stage, current_time,
					(int) dist_stage_period[infection_id][site_id][current_infection_stage].sample());

		} else {
			if (state_duration_preset <= 0) {
				while (infect_switch_time == current_time) {
					String key = String.format("%d,%d,%d", infection_id, site_id, current_infection_stage);
					double[] nextProb = lookupTable_infection_stage_path.get(key);
					if (nextProb == null) {
						// No next stage - return to susceptible
						int[] infect_switch_outcome = handleNoNextStage(pid, infection_id, site_id,
								current_infection_stage, current_time);

						current_infection_stage = infect_switch_outcome[0];
						infect_switch_time = infect_switch_outcome[1];

						break;
					} else {
						int state_pt = Arrays.binarySearch(nextProb, 0, nextProb.length / 2, RNG.nextDouble());
						if (state_pt <= 0) {
							state_pt = ~state_pt;
						}
						current_infection_stage = (int) nextProb[nextProb.length / 2 + state_pt];
						state_duration = dist_stage_period[infection_id][site_id][current_infection_stage].sample();
						infect_switch_time = (int) Math.round(current_time + state_duration);
					}
				}
			} else {
				infect_switch_time = current_time + state_duration_preset;
			}
		}

		// Update state_switch map
		current_stage_arr[infection_id][site_id] = current_infection_stage;

		if (infection_state_switch[infection_id][site_id] != infect_switch_time) {
			int org_switch = infection_state_switch[infection_id][site_id];
			infection_state_switch[infection_id][site_id] = infect_switch_time;
			// Update schedule
			updateInfectStageChangeSchedule(pid, infection_id, site_id, infect_switch_time, org_switch);
		}

		// Update infectious status
		String key = String.format("%d,%d", infection_id, site_id);
		ArrayList<Integer> currenty_infectious_ent = map_currently_infectious.get(key);

		if (((1 << current_infection_stage) & lookupTable_infection_infectious_stages[infection_id][site_id]) != 0) {
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

		// Check if the state has symptoms
		if (dist_sym_rate[infection_id][site_id] != null) {
			boolean hasSym = false;
			if (current_infection_stage >= 0 && current_infection_stage < dist_sym_rate[infection_id][site_id].length) {
				if (dist_sym_rate[infection_id][site_id][current_infection_stage] != null) {
					hasSym = RNG.nextDouble() < dist_sym_rate[infection_id][site_id][current_infection_stage].sample();
				}
			}
			// Check if seek treatment
			boolean seekTest = seek_test_rate[getGenderType(pid)] >= 1;
			if (seek_test_rate[getGenderType(pid)] < 1) {
				seekTest = RNG.nextDouble() < seek_test_rate[getGenderType(pid)];
			}
			if (hasSym & seekTest) {
				int seekTestAfter = (int) Math.max(1, Math.round(dist_seek_test_period[getGenderType(pid)].sample()));
				int nextTestDate = current_time + seekTestAfter;
				ArrayList<int[]> day_sch = schedule_testing.get(nextTestDate);
				if (day_sch == null) {
					day_sch = new ArrayList<>();
					schedule_testing.put(nextTestDate, day_sch);
				}
				int[] test_pair = new int[] { -pid, 1 << infection_id, 1 << site_id }; // -ive pid for symptomatic test

				int pt_t = Collections.binarySearch(day_sch, test_pair, new Comparator<int[]>() {
					@Override
					public int compare(int[] o1, int[] o2) {
						int res = 0;
						int pt = 0;
						while (res == 0 && pt < o1.length) {
							res = Integer.compare(Math.abs(o1[pt]), Math.abs(o2[pt]));
							pt++;
						}
						return res;
					}
				});
				if (pt_t < 0) {
					day_sch.add(~pt_t, test_pair);
				}
			}
		}

		return infect_switch_time;
	}

	protected int[] handleNoNextStage(Integer pid, int infection_id, int site_id, int current_infection_stage,
			int current_time) {
		int[] infect_switch_outcome; // {next_stage, duration}
		infect_switch_outcome = new int[] { AbstractIndividualInterface.INFECT_S, -1 };
		return infect_switch_outcome;
	}

	protected void updateInfectStageChangeSchedule(Integer pid, int infection_id, int site_id, int infect_switch_time,
			int org_infection_switch_time) {
		ArrayList<ArrayList<ArrayList<Integer>>> schedule_inf;
		ArrayList<ArrayList<Integer>> schedule_inf_site;
		ArrayList<Integer> schedule_inf_site_arr;
		int pt;

		if (org_infection_switch_time != infect_switch_time) {

			if (org_infection_switch_time > 0) {
				schedule_inf = schedule_stage_change.get(org_infection_switch_time);
				if (schedule_inf != null) { // Remove original schedule
					schedule_inf_site_arr = schedule_inf.get(infection_id).get(site_id);
					pt = Collections.binarySearch(schedule_inf_site_arr, pid);
					if (pt >= 0) {
						schedule_inf_site_arr.remove(pt);
					}
				}
			}

			if (infect_switch_time > 0) { // Never switch if < 0
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
					lastStateSwitch = Math.max(lastStateSwitch, infect_switch_time);

				}

				schedule_inf_site_arr = schedule_inf.get(infection_id).get(site_id);
				pt = Collections.binarySearch(schedule_inf_site_arr, pid);
				if (pt < 0) {
					schedule_inf_site_arr.add(~pt, pid);
				}
			}
		}
	}

	@Override
	public ArrayList<Integer> loadOptParameter(String[] parameter_settings, double[] point, int[][] seedInfectNum,
			boolean display_only) {

		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
				+ Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD;

		ArrayList<Integer> field_to_update = new ArrayList<>();
		for (int param_arr_index = 0; param_arr_index < parameter_settings.length; param_arr_index++) {
			String param_setting = parameter_settings[param_arr_index];
			param_setting = param_setting.replaceAll("\\s", "");
			String[] param_setting_arr = param_setting.split("_");
			int param_name_index = Integer.parseInt(param_setting_arr[0]);
			int field_id = param_name_index - Optimisation_Factory.RUNNABLE_OFFSET;

			if (field_id >= 0) { // Runnable field
				Object val = getRunnable_fields()[field_id];
				if (val != null) {
					int setting_level = 1;
					Optimisation_Factory.recursiveRunnableFieldReplace(val, param_arr_index, point, param_setting_arr,
							setting_level);
				}

				int pt = Collections.binarySearch(field_to_update, field_id);
				if (pt < 0) {
					field_to_update.add(~pt, field_id);
				}
			} else if (param_name_index == sim_offset + Simulation_ClusterModelTransmission.SIM_FIELD_SEED_INFECTION) {
				// Seed infection replacement
				int i_lvl = 1;
				int j_lvl = 2;
				for (int i = 0; i < seedInfectNum.length; i++) {
					int inc_index = Integer.parseInt(param_setting_arr[i_lvl]);
					if ((inc_index & 1 << i) != 0) {
						for (int j = 0; j < seedInfectNum[i].length; j++) {
							inc_index = Integer.parseInt(param_setting_arr[j_lvl]);
							if ((inc_index & 1 << j) != 0) {
								seedInfectNum[i][j] = (int) Math.round(point[param_arr_index]);
							}
						}
					}
				}
				System.out.printf("Seed infection replacement to %s\n", Arrays.deepToString(seedInfectNum));
			} else {
				String prop_str = (String) getSim_prop()
						.get(String.format("POP_PROP_INIT_PREFIX_%d", param_name_index));
				String prop_class = (String) getSim_prop()
						.get(String.format("POP_PROP_INIT_PREFIX_CLASS_%d", param_name_index));
				try {
					Object val = util.PropValUtils.propStrToObject(prop_str, Class.forName(prop_class));
					int setting_level = 1;
					Optimisation_Factory.recursiveRunnableFieldReplace(val, param_arr_index, point, param_setting_arr,
							setting_level);

					getSim_prop().put(String.format("POP_PROP_INIT_PREFIX_%d", param_name_index),
							util.PropValUtils.objectToPropStr(val, Class.forName(prop_class)));

				} catch (ClassNotFoundException ex) {
					ex.printStackTrace(System.err);

				}
			}
		}

		return field_to_update;
	}

}
