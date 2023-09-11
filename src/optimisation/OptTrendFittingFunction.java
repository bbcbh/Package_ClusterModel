package optimisation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import population.Population_Bridging;
import random.RandomGenerator;
import relationship.ContactMap;
import sim.Runnable_ClusterModel_ContactMap_Generation;
import sim.Runnable_ClusterModel_Transmission;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelGeneration;
import util.PropValUtils;

public class OptTrendFittingFunction extends OptFittingFunction {

	private static Pattern OPT_TREND_TYPE_FORMAT_BY_SITE = Pattern.compile("(.*)_(\\d+)");

	private final ContactMap[] baseCMaps;
	private final long[] baseCMapSeeds;
	private long[] sim_seeds;
	private final Properties prop;
	private final int num_threads;
	private final File baseDir;
	private final HashMap<String, double[][]> target_trend_collection;
	private final double[][] target_trend_time_range;
	private double[] bestResidue_by_runnable;

	// Arguments for OptTrendFittingFunction.calculate_residue_opt_trend
	public static final String ARGS_PROGRESS_DISP = "ARGS_PROGRESS_DISP";
	public static final String ARGS_OPT_METHOD = "ARGS_OPT_METHOD";
	public static final String ARGS_TAR_TRENDS_TIMERANGE = "ARGS_TAR_TRENDS_TIMERANGE";
	public static final String ARGS_TAR_TRENDS_COLLECTIONS = "ARGS_TAR_TRENDS_COLLECTIONS";
	public static final String ARGS_PROP = "ARGS_PROP";
	public static final String ARGS_BOUNDARIES = "ARGS_BOUNDARIES";
	public static final String ARGS_INIT_PARAM = "ARGS_INIT_PARAM";
	public static final String ARGS_NUM_EVAL = "ARGS_NUM_EVAL";
	public static final String ARGS_BASEDIR = "ARGS_BASEDIR";
	public static final String ARGS_SIM_SEED = "ARGS_SIM_SEED";
	public static final String ARGS_CMAP_SEED = "ARGS_CMAP_SEED";
	public static final String ARGS_CMAP = "ARGS_CMAP";
	public static final String ARGS_PREV_RESULTS = "ARGS_PREV_RESULTS";

	public static final String POP_PROP_OPT_PARAM_FIT_SETTING = "POP_PROP_OPT_PARAM_FIT_SETTING";
	// POP_PROP_OPT_PARAM_FIT_SETTING
	// Format: String[] { popPropInitPrefix_IncIndices_... , ...}
	// Examples:
	// 16_2_1_1 = TRANS_P2V
	// 16_1_2_1 = TRANS_V2P
	// 16_2_4_1 = TRANS_P2R
	// 16_4_2_1 = TRANS_R2P
	// 16_2_8_1 = TRANS_P2O
	// 16_8_2_1 = TRANS_O2P
	// 16_4_8_1 = TRANS_R2O
	// 16_8_4_1 = TRANS_O2R
	// 16_8_8_1 = TRANS_O2O
	// 17_1_1 = MEAN_DUR_V
	// 17_2_1 = MEAN_DUR_P
	// 17_4_1 = MEAN_DUR_R
	// 17_8_1 = MEAN_DUR_O
	// 19_14_2 = % Sym. for all male at P
	// 22_1 = Mean period sym hetro female of seeking treatment
	// 22_4 = Mean period sym hetro male of seeking treatment

	public static final Pattern POP_PROP_OPT_PARAM_FIT_SETTING_DIFF_FORMAT = Pattern.compile("Diff(\\d+)");
	// Eg. 18_Diff1 means Param[18] = opt. parameter value + value of opt parameter
	// with index 1

	// Opt trend output keys
	public static final String OPT_TREND_OUTPUT_RESULT_DISP = "OPT_TREND_OUTPUT_RESULT_DISP";
	public static final String OPT_TREND_OUTPUT_COUNT_BY_PERSON = "OPT_TREND_OUTPUT_COUNT_BY_PERSON";
	public static final String OPT_TREND_OUTPUT_COUNT_BY_SITE = "OPT_TREND_OUTPUT_COUNT_BY_SITE";
	public static final String OPT_TREND_OUTPUT_RUNNABLE = "OPT_TREND_OUTPUT_RUNNABLE";
	public static final String OPT_TREND_OUTPUT_BEST_RESIDUE = "OPT_TREND_OUTPUT_BEST_RESIDUE";
	public static final String OPT_TREND_CALLABLE_OUTPUT_BEST_SO_FAR = "OPT_TREND_CALLABLE_OUTPUT_BEST_SO_FAR";
	public static final String OPT_TREND_CALLABLE_OUTPUT_RESULT_KEY = "OPT_TREND_CALLABLE_OUTPUT_RESULT_KEY";

	// Opt trend count map index
	private static final int OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_PERSON = 0;
	private static final int OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_PERSON = OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_PERSON
			+ 1;
	private static final int OPT_TREND_COUNT_MAP_CUMUL_POS_DX_BY_PERSON = OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_PERSON
			+ 1;
	private static final int OPT_TREND_COUNT_MAP_CUMUL_POS_DX_SOUGHT_BY_PERSON = OPT_TREND_COUNT_MAP_CUMUL_POS_DX_BY_PERSON
			+ 1;
	private static final int LENGTH_OPT_TREND_COUNT_MAP_BY_PERSON = OPT_TREND_COUNT_MAP_CUMUL_POS_DX_SOUGHT_BY_PERSON
			+ 1;

	private static final int OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_SITE = 0;
	private static final int OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_SITE = OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_SITE
			+ 1;
	private static final int LENGTH_OPT_TREND_COUNT_MAP_BY_SITE = OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_SITE + 1;

	// Fields for input trend file
	public static final String OPT_TREND_MAP_KEY_FORMAT = "%s,%s,%d,%s,%d";
	public static final String OPT_TREND_INPUT_KEY_TYPE = "Type";
	public static final String OPT_TREND_INPUT_KEY_WEIGHT = "Weight";
	public static final String OPT_TREND_INPUT_KEY_FITFROM = "FitFrom";
	public static final String OPT_TREND_INPUT_TARGET_GRP = "Target_Grp";

	public static final String OPT_TREND_CSV_RANGE = "CSV_RANGE";
	private static final int OPT_TREND_MAP_KEY_PATH = 0;
	private static final int OPT_TREND_MAP_KEY_TYPE = OPT_TREND_MAP_KEY_PATH + 1;
	private static final int OPT_TREND_MAP_KEY_TARGET_GRP = OPT_TREND_MAP_KEY_TYPE + 1;
	private static final int OPT_TREND_MAP_KEY_WEIGHT = OPT_TREND_MAP_KEY_TARGET_GRP + 1;
	private static final int OPT_TREND_MAP_KEY_FITFROM = OPT_TREND_MAP_KEY_WEIGHT + 1;

	// Fields for output file
	public static final String OPT_TREND_FILE_NAME_TREND_OUTPUT = "Opt_trend_%d_%d_%d.txt";
	public static final String OPT_TREND_FILE_NAME_BEST_SO_FAR = "BestFit_%d_%d.txt";

	public static final String OPT_TREND_INPUT_TYPE_NUMINF = "NumInf";
	public static final String OPT_TREND_INPUT_TYPE_INCID = "CumulIncid";
	public static final String OPT_TREND_INPUT_TYPE_DX = "CumulDX";

	public static final String OPT_TREND_OUTPUT_PREFIX_CMAP = "CMAP    = ";
	public static final String OPT_TREND_OUTPUT_PREFIX_SIMSEED = "SimSeed = ";
	public static final String OPT_TREND_OUTPUT_PREFIX_PARAM = "Param   = ";
	public static final String OPT_TREND_OUTPUT_PREFIX_RESIDUE = "Residue = ";
	public static final String OPT_TREND_OUTPUT_PREFIX_OFFSET = "Offset  = ";

	public static final String OPT_SUMMARY_FILE = "Opt_Summary.csv";
	public static final String OPT_SUMMARY_TREND_FILE = "Opt_Summary_Trend_%d.csv";

	public OptTrendFittingFunction(File baseDir, Properties prop, ContactMap[] baseCMaps, long[] baseCMapSeeds,
			long[] sim_seeds, HashMap<String, double[][]> target_trend_collection, double[][] trend_time_range,
			int num_threads) {
		this.num_threads = num_threads;
		this.baseDir = baseDir;
		this.prop = prop;

		this.baseCMaps = baseCMaps;
		this.baseCMapSeeds = baseCMapSeeds;
		this.target_trend_collection = target_trend_collection;
		this.target_trend_time_range = trend_time_range;
		this.sim_seeds = sim_seeds;

		this.bestResidue_by_runnable = null;
	}

	public OptTrendFittingFunction(File baseDir, Properties prop, ContactMap[] baseCMaps, long[] baseCMapSeeds,
			int num_sim_per_map, RandomGenerator sim_seed_rng, HashMap<String, double[][]> target_trend_collection,
			double[][] trend_time_range, int num_threads) {

		this(baseDir, prop, baseCMaps, baseCMapSeeds, null, target_trend_collection, trend_time_range, num_threads);
		this.sim_seeds = new long[num_sim_per_map];
		if (sim_seed_rng != null) {
			for (int i = 0; i < sim_seeds.length; i++) {
				sim_seeds[i] = sim_seed_rng.nextLong();
			}
		}
	}

	@Override
	public long[] getSim_seeds() {
		return sim_seeds;
	}

	@Override
	public long[] getCMap_seeds() {
		return baseCMapSeeds;
	}

	@Override
	public Properties getProperties() {
		return prop;
	}

	@Override
	public double[] getBestResidue_by_runnable() {
		return bestResidue_by_runnable;
	}

	@Override
	public double value(double[] point) {
		long tic = System.currentTimeMillis();
		double best_fitting_sq_sum;

		ContactMap[] cMaps = baseCMaps;
		long[] cMap_seeds = baseCMapSeeds;

		if (baseCMaps.length == 0) {
			cMaps = new ContactMap[] { null }; // Special case for opt parameter value disp.
			cMap_seeds = new long[] { 0 };
		}

		HashMap<String, Object> cal_resiude_arg = new HashMap<>();
		HashMap<String, Object> cal_residue_output = new HashMap<>();

		cal_resiude_arg.put(OptTrendFittingFunction.ARGS_CMAP, cMaps);
		cal_resiude_arg.put(OptTrendFittingFunction.ARGS_CMAP_SEED, cMap_seeds);
		cal_resiude_arg.put(OptTrendFittingFunction.ARGS_BASEDIR, baseDir);
		cal_resiude_arg.put(OptTrendFittingFunction.ARGS_SIM_SEED, sim_seeds);
		cal_resiude_arg.put(OptTrendFittingFunction.ARGS_PROP, prop);
		cal_resiude_arg.put(OptTrendFittingFunction.ARGS_TAR_TRENDS_COLLECTIONS, target_trend_collection);
		cal_resiude_arg.put(OptTrendFittingFunction.ARGS_TAR_TRENDS_TIMERANGE, target_trend_time_range);

		if (prop != null && prop.containsKey(ARGS_PREV_RESULTS)) {
			cal_resiude_arg.put(ARGS_PREV_RESULTS, prop.get(ARGS_PREV_RESULTS));
		}

		bestResidue_by_runnable = OptTrendFittingFunction.calculate_residue_opt_trend(point, cal_resiude_arg,
				cal_residue_output, num_threads);

		best_fitting_sq_sum = 0;

		for (double r : bestResidue_by_runnable) {
			best_fitting_sq_sum += r;
		}

		String outMsg;

		StringBuilder param_disp = new StringBuilder();
		for (double pt : point) {
			if (param_disp.length() != 0) {
				param_disp.append(',');
			}
			param_disp.append(String.format("%.5f", pt));
		}

		if (bestResidue_by_runnable.length == 1) {
			outMsg = String.format("P = [%s], V = %.2e, map_seed= %d, sim_seed = %d, Time req = %.3fs\n",
					param_disp.toString(), best_fitting_sq_sum, cMap_seeds[0], sim_seeds[0],
					(System.currentTimeMillis() - tic) / 1000f);

		} else {
			outMsg = String.format("P = [%s], V = %f, Time req = %.3fs\n", param_disp.toString(), best_fitting_sq_sum,
					(System.currentTimeMillis() - tic) / 1000f);
		}

		try {
			File opt_output_file = new File(baseDir, Optimisation_Factory.FILENAME_OPT_RESULT);
			FileWriter fWri = new FileWriter(opt_output_file, true);
			PrintWriter pWri = new PrintWriter(fWri);
			pWri.print(outMsg);
			pWri.close();
			fWri.close();

		} catch (IOException ex) {
			ex.printStackTrace(System.err);
		}

		System.out.println(outMsg);

		return best_fitting_sq_sum;
	}

	public static void extractBestOptrendResults(File basedir, Pattern subDirPattern) throws IOException {

		final Pattern bestFileFile_pattern = Pattern
				.compile(OPT_TREND_FILE_NAME_BEST_SO_FAR.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));
		File[] subDirs = basedir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && subDirPattern.matcher(pathname.getName()).matches();
			}
		});

		// double residue, long cMapSeed, long simSeed, long offset, String[] param_str,
		// String[] trends
		ArrayList<Object[]> resultsCollections = new ArrayList<>();
		final Comparator<Object[]> resultCollectionsComp = new Comparator<Object[]>() {
			@Override
			public int compare(Object[] o1, Object[] o2) {
				int r = 0;
				for (int p = 0; (r == 0) && (p < o1.length); p++) {
					if (o1[p] instanceof Double) {
						r = Double.compare((Double) o1[p], (Double) o2[p]);
					} else if (o1[p] instanceof Long) {
						r = Long.compare((Long) o1[p], (Long) o2[p]);
					} else {
						String[] s1 = (String[]) o1[p];
						String[] s2 = (String[]) o2[p];
						for (int s = 0; ((r == 0) && s < s1.length); s++) {
							r = s1[s].compareTo(s2[s]);
						}
					}
				}
				return r;
			}
		};

		for (File subDir : subDirs) {
			File[] bestFitFiles = subDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return bestFileFile_pattern.matcher(pathname.getName()).matches();
				}
			});
			for (File bestFitFile : bestFitFiles) {

				BufferedReader reader = new BufferedReader(new FileReader(bestFitFile));
				String line;

				Matcher m = bestFileFile_pattern.matcher(bestFitFile.getName());
				m.matches();

				long cMapSeed = Long.parseLong(m.group(1));
				long simSeed = Long.parseLong(m.group(2));

				while ((line = reader.readLine()) != null) {
					if (line.startsWith(OPT_TREND_OUTPUT_PREFIX_PARAM)) {
						String[] param_str = line
								.substring(OPT_TREND_OUTPUT_PREFIX_PARAM.length() + 1, line.length() - 1).split(",");

						line = reader.readLine();
						double residue = Double.parseDouble(line.substring(OPT_TREND_OUTPUT_PREFIX_RESIDUE.length()));

						line = reader.readLine();
						long offset = Integer.parseInt(line.substring(OPT_TREND_OUTPUT_PREFIX_OFFSET.length()));

						ArrayList<String> trend_arr = new ArrayList<>();

						// TrendEntry
						while ((line = reader.readLine()) != null) {
							if(!line.isBlank()) {
								trend_arr.add(line);
							}
						}
						String[] trends = trend_arr.toArray(new String[trend_arr.size()]);
						Object[] val = new Object[] { residue, cMapSeed, simSeed, offset, param_str, trends };

						int key = Collections.binarySearch(resultsCollections, val, resultCollectionsComp);
						resultsCollections.add(~key, val);
					}
				}
				reader.close();
			}
		}

		// Print of results
		PrintWriter pWri_summary = new PrintWriter(new File(basedir, OPT_SUMMARY_FILE));
		pWri_summary.println("Residue,CMapSeed,SimSeed,Param");

		PrintWriter[] pWri_trend = null;

		for (Object[] val : resultsCollections) {
			pWri_summary.printf("%f,%d,%d", (Double) val[0], (Long) val[1], (Long) val[2]);
			String[] param = (String[]) val[4];
			for (int s = 0; s < param.length; s++) {
				pWri_summary.print(',');
				pWri_summary.print(param[s]);
			}
			pWri_summary.println();

			String[] trends = (String[]) val[5];
			int offset = ((Long) val[3]).intValue();

			if (pWri_trend == null) {
				String[] firstLine = trends[0].split(",");
				pWri_trend = new PrintWriter[firstLine.length];
				for (int w = 1; w < pWri_trend.length; w++) {
					pWri_trend[w] = new PrintWriter(new File(basedir, String.format(OPT_SUMMARY_TREND_FILE, w)));
				}
			}

			StringBuilder timeline = new StringBuilder();

			for (int t = 0; t < trends.length; t++) {
				String[] ent = trends[t].split(",");
				int time_adj = Integer.parseInt(ent[0]) - offset;
				if (t != 0) {
					timeline.append(',');
				}
				timeline.append(time_adj);
			}

			for (int v = 1; v < pWri_trend.length; v++) {
				pWri_trend[v].println(timeline.toString());
			}

			for (int t = 0; t < trends.length; t++) {
				String[] ent = trends[t].split(",");
				for (int v = 1; v < ent.length; v++) {
					if (t != 0) {
						pWri_trend[v].print(',');
					}
					pWri_trend[v].print(ent[v]);
				}
			}

			for (int v = 1; v < pWri_trend.length; v++) {
				pWri_trend[v].println();
			}

		}

		pWri_summary.close();
		for (PrintWriter pWri : pWri_trend) {
			if (pWri != null) {
				pWri.close();
			}
		}
		
		System.out.printf("Export opt trend result completed at %s./n", basedir.getAbsolutePath());

	}

	public static HashMap<String, double[][]> loadTrendCSV(Properties prop) throws FileNotFoundException, IOException {
		HashMap<String, double[][]> target_trend_collection = new HashMap<>();
		if (prop.getProperty(Optimisation_Factory.POP_PROP_OPT_TARGET) != null) {
			double[] time_range = new double[] { Double.NaN, Double.NaN };

			String[] target_trent_csv = prop.getProperty(Optimisation_Factory.POP_PROP_OPT_TARGET).replaceAll("\n", "")
					.split(",");
			for (String csv_path : target_trent_csv) {
				File csv_file = new File(csv_path);
				ArrayList<Double> fit_t = new ArrayList<>();
				ArrayList<Double> fit_y = new ArrayList<>();
				String line, type = null;
				int tar_grp = -1;
				int fitFrom = -1;
				double weight = -1;
				BufferedReader csv_reader = new BufferedReader(new FileReader(csv_file));
				while ((line = csv_reader.readLine()) != null) {
					String[] ent = line.split(",");
					if (OPT_TREND_INPUT_KEY_TYPE.equals(ent[0])) {
						type = ent[1];
					} else if (OPT_TREND_INPUT_KEY_FITFROM.equals(ent[0])) {
						fitFrom = Integer.parseInt(ent[1]);
					} else if (OPT_TREND_INPUT_TARGET_GRP.equals(ent[0])) {
						tar_grp = Integer.parseInt(ent[1]);
					} else if (OPT_TREND_INPUT_KEY_WEIGHT.equals(ent[0])) {
						weight = Double.parseDouble(ent[1]);
					} else {
						double t = Double.parseDouble(ent[0]);
						double y = Double.parseDouble(ent[1]);
						int key = Collections.binarySearch(fit_t, t);
						if (key < 0) {
							fit_t.add(~key, t);
							fit_y.add(~key, y);
						}
						if (Double.isNaN(time_range[0])) {
							time_range[0] = t;
							time_range[1] = t;

						} else {
							time_range[0] = Math.min(time_range[0], t);
							time_range[1] = Math.max(time_range[1], t);
						}
					}
				}
				csv_reader.close();

				if (fit_t.size() > 0 && fit_y.size() == fit_t.size()) {
					String mapKey = String.format(OPT_TREND_MAP_KEY_FORMAT, csv_path, type, tar_grp,
							Double.toString(weight), fitFrom);
					double[][] mapEntry = new double[2][fit_t.size()];
					for (int i = 0; i < fit_t.size(); i++) {
						mapEntry[0][i] = fit_t.get(i);
						mapEntry[1][i] = fit_y.get(i);
					}
					target_trend_collection.put(mapKey, mapEntry);
				}

			}
			if (target_trend_collection.size() > 0) {
				target_trend_collection.put(OPT_TREND_CSV_RANGE, new double[][] { time_range });
			}

		}
		return target_trend_collection;
	}

	@SuppressWarnings("unchecked")
	public static double[] calculate_residue_opt_trend(double[] point, HashMap<String, Object> args,
			HashMap<String, Object> outputMap, final int NUM_THREADS) {

		ContactMap[] cMap = (ContactMap[]) args.get(OptTrendFittingFunction.ARGS_CMAP);
		long[] cMap_seed = (long[]) args.get(OptTrendFittingFunction.ARGS_CMAP_SEED);
		long[] sim_seeds = (long[]) args.get(OptTrendFittingFunction.ARGS_SIM_SEED);
		int NUM_SIM_PER_MAP = sim_seeds.length;

		double[] bestResidue_by_runnable;
		bestResidue_by_runnable = new double[cMap.length * NUM_SIM_PER_MAP];
		Arrays.fill(bestResidue_by_runnable, Double.NaN);

		// From args
		File baseDir = (File) args.get(OptTrendFittingFunction.ARGS_BASEDIR);
		Properties prop = (Properties) args.get(OptTrendFittingFunction.ARGS_PROP);
		HashMap<String, double[][]> target_trend_collection = (HashMap<String, double[][]>) args
				.get(OptTrendFittingFunction.ARGS_TAR_TRENDS_COLLECTIONS);
		double[][] time_range = (double[][]) args.get(OptTrendFittingFunction.ARGS_TAR_TRENDS_TIMERANGE);
		Number[][] result_lookup = (Number[][]) args.get(ARGS_PREV_RESULTS);

		// From properties
		int[] POP_COMPOSITION = new int[] { 500000, 500000, 20000, 20000 };
		String popCompositionKey = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX
				+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
		if (prop.containsKey(popCompositionKey)) {
			POP_COMPOSITION = (int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey), int[].class);
		}

		int NUM_SNAP = 1;
		if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])) {
			NUM_SNAP = Integer
					.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
		}
		int NUM_TIME_STEPS_PER_SNAP = 1;
		if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
			NUM_TIME_STEPS_PER_SNAP = Integer
					.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ]));
		}
		int[][] SEED_INFECTION = null;
		if (prop.containsKey(Optimisation_Factory.PROP_SEED_INFECTION)) {
			SEED_INFECTION = (int[][]) PropValUtils
					.propStrToObject(prop.getProperty(Optimisation_Factory.PROP_SEED_INFECTION), int[][].class);
		}
		int START_TIME = 365;
		String contactMapRangeKey = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX
				+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
						+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
						+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);
		if (prop.containsKey(contactMapRangeKey)) {
			START_TIME = ((int[]) PropValUtils.propStrToObject(prop.getProperty(contactMapRangeKey), int[].class))[0];
		}

		Runnable_ClusterModel_Transmission[] runnable = new Runnable_ClusterModel_Transmission[bestResidue_by_runnable.length];
		String[] res_disp_all = new String[runnable.length];

		int[] bestMatchStart_by_runnable = new int[runnable.length];
		Arrays.fill(bestMatchStart_by_runnable, -1);

		StringBuilder param_str = new StringBuilder();
		for (double pt : point) {
			if (param_str.length() != 0) {
				param_str.append(',');
			}
			param_str.append(String.format("%f", pt));
		}

		int rId = 0;
		int cMap_id = 0;
		ExecutorService exec = null;

		for (ContactMap c : cMap) {
			for (int r = 0; r < NUM_SIM_PER_MAP; r++) {
				long sim_seed = sim_seeds[r];
				if (result_lookup != null) {
					Number[] val = new Number[point.length + 3];
					val[0] = cMap_seed[cMap_id];
					val[1] = sim_seed;
					for (int p = 2; p < val.length - 1; p++) {
						val[p] = point[p - 2];
					}
					int k = Arrays.binarySearch(result_lookup, val,
							Optimisation_Factory.COMPARATOR_RESULT_LOOKUP_PARAM);
					if (k >= 0) {
						bestResidue_by_runnable[rId] = (double) val[val.length - 1];
						runnable[rId] = null;
					}
				}
				if (Double.isNaN(bestResidue_by_runnable[rId])) {
					runnable[rId] = new Runnable_ClusterModel_Transmission(cMap_seed[cMap_id], sim_seed,
							POP_COMPOSITION, c, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP);
					runnable[rId].setBaseDir(baseDir);

					for (int i = Optimisation_Factory.RUNNABLE_OFFSET; i < Optimisation_Factory.RUNNABLE_OFFSET
							+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; i++) {

						String key = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX + Integer.toString(i);
						if (prop.containsKey(key)) {
							runnable[rId].getRunnable_fields()[i - Optimisation_Factory.RUNNABLE_OFFSET] = PropValUtils
									.propStrToObject(prop.getProperty(key),
											runnable[rId].getRunnable_fields()[i - Optimisation_Factory.RUNNABLE_OFFSET]
													.getClass());
						}
					}
					runnable[rId].setSimSetting(1); // No output
					Optimisation_Factory.setOptParamInRunnable(runnable[rId], prop, point, c == null);
					runnable[rId].initialse();
					runnable[rId].allocateSeedInfection(SEED_INFECTION, START_TIME);
				}
				rId++;
			}
			cMap_id++;
		}

		if (rId == 1 || NUM_THREADS <= 1) {
			for (int r = 0; r < rId; r++) {
				if (runnable[r] != null) {
					runnable[r].run();
				}
			}
		} else {
			exec = Executors.newFixedThreadPool(NUM_THREADS);
			for (int r = 0; r < rId; r++) {
				if (runnable[r] != null) {
					exec.submit(runnable[r]);
				}
			}
			exec.shutdown();
			try {
				if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
					throw new InterruptedException("Time out");
				}
			} catch (InterruptedException e) {
				e.printStackTrace(System.err);
				for (int r = 0; r < rId; r++) {
					if (runnable[r] != null) {
						runnable[r].run();
					}
				}
			}
		}
		HashMap<Integer, int[][]>[] countMapBySite = new HashMap[OptTrendFittingFunction.LENGTH_OPT_TREND_COUNT_MAP_BY_SITE];
		HashMap<Integer, int[]>[] countMapByPerson = new HashMap[OptTrendFittingFunction.LENGTH_OPT_TREND_COUNT_MAP_BY_PERSON];

		for (int r = 0; r < rId; r++) {

			if (runnable[r] != null) {

				countMapBySite[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_SITE] = (HashMap<Integer, int[][]>) runnable[r]
						.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT);
				countMapBySite[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_SITE] = (HashMap<Integer, int[][]>) runnable[r]
						.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_INCIDENCE);

				countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
						.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT_BY_PERSON);
				countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
						.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_INCIDENCE_BY_PERSON);
				countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_POS_DX_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
						.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_POS_DX_BY_PERSON);
				countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_POS_DX_SOUGHT_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
						.getSim_output()
						.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_POS_DX_SOUGHT_BY_PERSON);

				Integer[] simTime = null;

				for (HashMap<Integer, int[][]> map : countMapBySite) {
					if (map != null) {
						if (simTime == null || map.keySet().size() > simTime.length) {
							simTime = map.keySet().toArray(new Integer[map.size()]);
						}
					}
				}
				for (HashMap<Integer, int[]> map : countMapByPerson) {
					if (map != null) {
						if (simTime == null || map.keySet().size() > simTime.length) {
							simTime = map.keySet().toArray(new Integer[map.size()]);
						}
					}
				}

				if (simTime == null) {
					System.err.printf(
							"Simulation results not found for simulation with CMAP_SEED = %d and SIM_SEED = %d.\n",
							runnable[r].getcMap_seed(), runnable[r].getSim_seed());
				} else {
					Arrays.sort(simTime);
					int num_target_trend = target_trend_collection.size();
					String[] trend_target_key = target_trend_collection.keySet().toArray(new String[num_target_trend]);

					Arrays.sort(trend_target_key);
					double[] t_values = new double[simTime.length];
					for (int i = 0; i < simTime.length; i++) {
						t_values[i] = simTime[i];
					}

					// Trim simTime to valid starting point

					int minFitFrom = simTime[0];
					String[][] trend_target_key_split = new String[trend_target_key.length][];
					for (int trend_target_pt = 0; trend_target_pt < num_target_trend; trend_target_pt++) {
						trend_target_key_split[trend_target_pt] = trend_target_key[trend_target_pt].split(",");
						minFitFrom = Math.max(minFitFrom, Integer.parseInt(
								trend_target_key_split[trend_target_pt][OptTrendFittingFunction.OPT_TREND_MAP_KEY_FITFROM]));
					}

					if (minFitFrom > simTime[0]) {
						int trimFrom = Arrays.binarySearch(simTime, minFitFrom);
						if (trimFrom < 0) {
							trimFrom = ~trimFrom - 1;
						}
						simTime = Arrays.copyOfRange(simTime, Math.max(0, trimFrom), simTime.length);
					}

					double[][][] tar_values = new double[num_target_trend][][];
					double[] weight = new double[num_target_trend];
					UnivariateFunction[] interpolation = new PolynomialSplineFunction[num_target_trend];
					UnivariateInterpolator interpolator = new LinearInterpolator();

					StringBuilder[] str_disp = new StringBuilder[simTime.length];

					// Fill interpolation, tar_values and weight for trend_target
					for (int trend_target_pt = 0; trend_target_pt < num_target_trend; trend_target_pt++) {
						String[] trend_keys = trend_target_key_split[trend_target_pt];
						String type_key = trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE];
						int site = -1;
						int incl_grp = Integer
								.parseInt(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_TARGET_GRP]);

						double[] y_values = new double[t_values.length];

						tar_values[trend_target_pt] = target_trend_collection.get(trend_target_key[trend_target_pt]);
						weight[trend_target_pt] = Double
								.parseDouble(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_WEIGHT]);

						Matcher m = OPT_TREND_TYPE_FORMAT_BY_SITE
								.matcher(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE]);
						if (m.find()) {
							site = Integer.parseInt(m.group(2));
							type_key = m.group(1);
						}

						for (int t_pt = 0; t_pt < simTime.length; t_pt++) {
							Integer time = simTime[t_pt];

							if (str_disp[t_pt] == null) {
								str_disp[t_pt] = new StringBuilder();
								str_disp[t_pt].append(time);
							}

							double newVal = 0;

							if (OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_NUMINF.equals(type_key)
									|| OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_INCID.equals(type_key)) {
								if (site < 0) {
									// V= int[gender]
									int[] ent = OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_NUMINF.equals(type_key)
											? countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_PERSON]
													.get(time)
											: countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_PERSON]
													.get(time);
									if (ent != null) {
										for (int g = 0; g < ent.length; g++) {
											if ((incl_grp & 1 << g) != 0) {
												newVal += ent[g];
											}
										}
									}
								} else {
									// V= int[gender][site]
									int[][] ent = OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_NUMINF.equals(type_key)
											? countMapBySite[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_SITE]
													.get(time)
											: countMapBySite[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_SITE]
													.get(time);

									if (ent != null) {
										for (int g = 0; g < ent.length; g++) {
											if ((incl_grp & 1 << g) != 0) {
												newVal += ent[g][site];
											}
										}
									}
								}
							} else if (OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_DX.equals(type_key)) {
								// V= int[gender_positive_dx, gender_true_infection]
								int[][] ent_collection = new int[][] {
										countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_POS_DX_BY_PERSON]
												.get(time),
										countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_POS_DX_SOUGHT_BY_PERSON]
												.get(time),

								};
								for (int[] ent : ent_collection) {
									if (ent != null) {
										for (int g = 0; g < ent.length / 2; g++) {
											if ((incl_grp & 1 << g) != 0) {
												newVal += ent[g];
											}
										}
									}
								}
							}
							y_values[t_pt] = newVal;
							str_disp[t_pt].append(',');
							str_disp[t_pt].append(newVal);

						}

						interpolation[trend_target_pt] = interpolator.interpolate(t_values, y_values);

					}
					// Calculate best fit for all target
					boolean hasCompleteRun = simTime[0] < simTime[simTime.length - 1]
							- (time_range[0][1] - time_range[0][0]);

					if (hasCompleteRun) {
						for (int match_start_time = simTime[0]; match_start_time < simTime[simTime.length - 1]
								- (time_range[0][1] - time_range[0][0]); match_start_time++) {
							double residue = 0;

							for (int trend_target_pt = 0; trend_target_pt < num_target_trend; trend_target_pt++) {
								double offset = 0;
								if (trend_target_key_split[trend_target_pt][OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE]
										.startsWith("Cumul")) {
									offset = interpolation[trend_target_pt].value(match_start_time);
								}
								for (int i = 0; i < tar_values[trend_target_pt][0].length; i++) {
									double model_adj_tar_t = match_start_time + tar_values[trend_target_pt][0][i];
									double target_y = tar_values[trend_target_pt][1][i];
									double model_y = interpolation[trend_target_pt].value(model_adj_tar_t);

									residue += weight[trend_target_pt] * Math.pow((model_y - offset) - target_y, 2);
								}

							}
							// System.out.printf("Start_time = %d, R = %f\n", match_start_time,
							// residue);
							if (bestResidue_by_runnable[r] > residue) {
								bestResidue_by_runnable[r] = residue;
								bestMatchStart_by_runnable[r] = match_start_time;
							}
						}
					} else {
						// Has no complete run - e.g. due to extinction
						double residue = 0;
						int match_start_time = simTime[0];
						for (int trend_target_pt = 0; trend_target_pt < num_target_trend; trend_target_pt++) {
							double offset = 0;
							double no_match_val = 0;
							if (trend_target_key_split[trend_target_pt][OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE]
									.startsWith("Cumul")) {
								offset = interpolation[trend_target_pt].value(match_start_time);
								no_match_val = interpolation[trend_target_pt].value(simTime[simTime.length - 1]);

							}
							for (int i = 0; i < tar_values[trend_target_pt][0].length; i++) {
								double model_adj_tar_t = match_start_time + tar_values[trend_target_pt][0][i];
								double target_y = tar_values[trend_target_pt][1][i];
								double model_y;

								if (model_adj_tar_t <= simTime[simTime.length - 1]) {
									model_y = interpolation[trend_target_pt].value(model_adj_tar_t);
								} else {
									model_y = no_match_val;
								}

								residue += weight[trend_target_pt] * Math.pow((model_y - offset) - target_y, 2);
							}

						}
						// System.out.printf("Start_time = %d, R = %f\n", match_start_time,
						// residue);
						if (residue < bestResidue_by_runnable[r]) {
							bestResidue_by_runnable[r] = residue;
							bestMatchStart_by_runnable[r] = match_start_time;
						}

					}

					// Display trends
					try {
						File opt_output_file;

						opt_output_file = new File(baseDir, String.format(OPT_TREND_FILE_NAME_TREND_OUTPUT,
								runnable[r].getcMap_seed(), runnable[r].getSim_seed(), r));

						boolean newFile = !opt_output_file.exists();

						FileWriter fWri = new FileWriter(opt_output_file, true);
						PrintWriter pWri = new PrintWriter(fWri);

						if (newFile) {
							pWri.printf("Fitting %s target trend(s):\n", num_target_trend);
							for (String sKey : trend_target_key) {
								pWri.println(sKey);
							}
							pWri.println();
						}

						pWri.printf("%s%d\n", OPT_TREND_OUTPUT_PREFIX_CMAP, runnable[r].getcMap_seed());
						pWri.printf("%s%d\n", OPT_TREND_OUTPUT_PREFIX_SIMSEED, runnable[r].getSim_seed());
						pWri.printf("%s[%s]\n", OPT_TREND_OUTPUT_PREFIX_PARAM, param_str.toString());
						pWri.printf("%s%f\n", OPT_TREND_OUTPUT_PREFIX_RESIDUE, bestResidue_by_runnable[r]);
						pWri.printf("%s%d\n", OPT_TREND_OUTPUT_PREFIX_OFFSET, bestMatchStart_by_runnable[r]);

						for (StringBuilder s : str_disp) {
							pWri.println(s.toString());
						}

						pWri.println();

						pWri.close();
						fWri.close();

					} catch (IOException e) {
						e.printStackTrace(System.err);
					}

					StringBuilder res_disp = new StringBuilder();
					res_disp.append(String.format("Fitting %s target trend(s):\n", num_target_trend));
					for (String sKey : trend_target_key) {
						res_disp.append(sKey);
						res_disp.append('\n');
					}
					res_disp.append('\n');
					res_disp.append(String.format("%s[%s]\n", OPT_TREND_OUTPUT_PREFIX_PARAM, param_str.toString()));
					res_disp.append(
							String.format("%s%f\n", OPT_TREND_OUTPUT_PREFIX_RESIDUE, bestResidue_by_runnable[r]));
					res_disp.append(
							String.format("%s%d\n", OPT_TREND_OUTPUT_PREFIX_OFFSET, bestMatchStart_by_runnable[r]));
					for (StringBuilder s : str_disp) {
						res_disp.append(s.toString());
						res_disp.append('\n');
					}
					res_disp_all[r] = res_disp.toString();

				} // if (simTime != null) {

			} // if(runnable[r] != null) {

		}
		outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_RUNNABLE, runnable);
		outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_RESULT_DISP, res_disp_all);
		outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_BEST_RESIDUE, bestResidue_by_runnable);
		outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_COUNT_BY_SITE, countMapBySite);
		outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_COUNT_BY_PERSON, countMapByPerson);

		return bestResidue_by_runnable;
	}

}