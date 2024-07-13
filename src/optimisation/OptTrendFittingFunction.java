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
import org.apache.commons.math3.exception.OutOfRangeException;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import random.RandomGenerator;
import relationship.ContactMap;
import sim.Abstract_Runnable_ClusterModel_Transmission;
import sim.Runnable_ClusterModel_Bali;
import sim.Runnable_ClusterModel_ContactMap_Generation;
import sim.Runnable_ClusterModel_MultiTransmission;
import sim.Runnable_ClusterModel_Transmission;
import sim.Runnable_ClusterModel_Viability;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelGeneration;
import sim.Simulation_ClusterModelTransmission;
import util.PropValUtils;

public class OptTrendFittingFunction extends OptFittingFunction {

	public static Pattern OPT_TREND_TYPE_FORMAT_BY_SITE = Pattern.compile("(.*)_(\\d+)");

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
	public static final String ARGS_VERBOSE = "ARGS_VERBOSE";
	public static final String ARGS_SEEDLIST = "ARGS_SEEDLIST";

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
	public static final String POP_PROP_OPT_PARAM_FIT_SETTING = "POP_PROP_OPT_PARAM_FIT_SETTING";

	// POP_PROP_OPT_PARAM_TRANSFORM
	// Format: String[][] {popPropInitPrefix_IncIndices,
	// popPropInitPrefix_IncIndices_src_1, ratio_1 ...}
	// e.g. [A, *B, 0.2, C, 0.1, Const ] => A = A * 0.2 B + 0.1 C + Const
	public static final String POP_PROP_OPT_PARAM_TRANSFORM = "POP_PROP_OPT_PARAM_TRANSFORM";

	public static final Pattern POP_PROP_OPT_PARAM_FIT_SETTING_DIFF_FORMAT = Pattern.compile("Diff(\\d+)");
	// Eg. 18_Diff1 means Param[18] = opt. parameter value + value of opt parameter
	// with index 1

	// Opt trend output keys
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
	public static final int OPT_TREND_MAP_KEY_PATH = 0;
	public static final int OPT_TREND_MAP_KEY_TYPE = OPT_TREND_MAP_KEY_PATH + 1;
	public static final int OPT_TREND_MAP_KEY_TARGET_GRP = OPT_TREND_MAP_KEY_TYPE + 1;
	public static final int OPT_TREND_MAP_KEY_WEIGHT = OPT_TREND_MAP_KEY_TARGET_GRP + 1;
	public static final int OPT_TREND_MAP_KEY_FITFROM = OPT_TREND_MAP_KEY_WEIGHT + 1;

	// Fields for output file
	public static final String OPT_TREND_FILE_NAME_TREND_OUTPUT = "Opt_trend_%d_%d_%d.txt";

	public static final String OPT_TREND_INPUT_TYPE_NUMINF = "NumInf";
	public static final String OPT_TREND_INPUT_TYPE_INCID = "CumulIncid";
	public static final String OPT_TREND_INPUT_TYPE_DX = "CumulDX";
	public static final String OPT_TREND_INPUT_TYPR_DX_TESTING = "CumulDX_Test";

	public static final String OPT_TREND_OUTPUT_PREFIX_OFFSET = "Offset  = ";

	public static final String OPT_SUMMARY_FILE = "Opt_Summary.csv";
	public static final String OPT_SUMMARY_UNIQUE_FILE = "Opt_Summary_Unique.csv";
	public static final String OPT_SUMMARY_TREND_FILE = "Opt_Summary_Trend_%d.csv";
	public static final String OPT_SUMMARY_TREND_UNIQUE_FILE = "Opt_Summary_Trend_%d_Unique.csv";

	private static final Pattern[] pattern_multi_trans = new Pattern[] {
			Runnable_ClusterModel_MultiTransmission.PROP_TYPE_PATTERN,
			Runnable_ClusterModel_Viability.PROP_TYPE_PATTERN, Runnable_ClusterModel_Bali.PROP_TYPE_PATTERN };

	private static final int pattern_multi_trans_index_base = 0;
	private static final int pattern_multi_trans_index_viability = pattern_multi_trans_index_base + 1;
	private static final int pattern_multi_trans_index_Bali = pattern_multi_trans_index_viability + 1;

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

		if (prop != null) {
			if (prop.containsKey(ARGS_PREV_RESULTS)) {
				cal_resiude_arg.put(ARGS_PREV_RESULTS, prop.get(ARGS_PREV_RESULTS));
			}
			if (prop.containsKey(ARGS_VERBOSE)) {
				cal_resiude_arg.put(ARGS_VERBOSE, prop.get(ARGS_VERBOSE));
			}
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

	public static void extractUniqueTrendResults(File version_baseDir, File[] trend_file)
			throws FileNotFoundException, IOException {
		UnivariateInterpolator interpolator = new LinearInterpolator();

		String[] param_summary = util.Util_7Z_CSV_Entry_Extract_Callable
				.extracted_lines_from_text(new File(version_baseDir, OPT_SUMMARY_FILE));

		String[] param_summary_unique = util.Util_7Z_CSV_Entry_Extract_Callable
				.extracted_lines_from_text(new File(version_baseDir, OPT_SUMMARY_UNIQUE_FILE));

		for (int trendPt = 0; trendPt < trend_file.length; trendPt++) {
			// Key = Parameter, Value = Sq Sum
			HashMap<String, ArrayList<Double>> trendValueCollection = new HashMap<>();
			String[] trend_values = util.Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(
					new File(version_baseDir, String.format(OPT_SUMMARY_TREND_FILE, trendPt + 1)));

			String[] trend_setting = util.Util_7Z_CSV_Entry_Extract_Callable
					.extracted_lines_from_text(trend_file[trendPt]);

			int fitFrom = 0;
			float weight = 0;
			ArrayList<Double> t_val = new ArrayList<>();
			ArrayList<Double> y_val = new ArrayList<>();
			for (String line : trend_setting) {
				String[] line_sp = line.split(",");
				if (line_sp[0].equals(OPT_TREND_INPUT_KEY_FITFROM)) {
					fitFrom = Integer.parseInt(line_sp[1]);
				} else if (line_sp[0].equals(OPT_TREND_INPUT_KEY_WEIGHT)) {
					weight = Float.parseFloat(line_sp[1]);
				} else {
					if (Pattern.matches("^\\d+$", line_sp[0])) {
						t_val.add(Double.parseDouble(line_sp[0]));
						y_val.add(Double.parseDouble(line_sp[1]));
					}
				}
			}
			double[][] trend_tar_val = new double[2][t_val.size()];
			for (int i = 0; i < t_val.size(); i++) {
				trend_tar_val[0][i] = t_val.get(i) + fitFrom;
				trend_tar_val[1][i] = y_val.get(i);
			}

			UnivariateFunction interpolation_target = interpolator.interpolate(trend_tar_val[0], trend_tar_val[1]);
			String[] row, param_str;
			double[][] trend_sim_val;
			String param_key;
			int num_param = 0;

			for (int p = 1; p < param_summary.length; p++) {
				row = param_summary[p].split(",");
				param_str = Arrays.copyOfRange(row, 4, row.length);
				num_param = param_str.length;
				param_key = Arrays.toString(param_str);
				ArrayList<Double> sq_sum = trendValueCollection.get(param_key);

				if (sq_sum == null) {
					sq_sum = new ArrayList<>();
					trendValueCollection.put(param_key, sq_sum);
				}
				row = trend_values[p].split(",");

				trend_sim_val = new double[2][row.length / 2];
				for (int c = 0; c < row.length / 2; c++) {
					trend_sim_val[0][c] = Double.parseDouble(row[c]);
					trend_sim_val[1][c] = Double.parseDouble(row[row.length / 2 + c]);
				}

				UnivariateFunction interpolation_sim = interpolator.interpolate(trend_sim_val[0], trend_sim_val[1]);

				double sq_sum_val = 0;

				for (int t = (int) trend_tar_val[0][0]; t <= trend_tar_val[0][trend_tar_val[0].length
						- 1]; t += AbstractIndividualInterface.ONE_YEAR_INT) {
					double sim_val;
					try {
						sim_val = interpolation_sim.value(t);
						if (weight < 0) {
							sim_val -= interpolation_sim.value(t - AbstractIndividualInterface.ONE_YEAR_INT);
						}
					} catch (OutOfRangeException ex) {

						sim_val = 0;
					}

					sq_sum_val += Math.pow(sim_val - interpolation_target.value(t), 2);
				}
				sq_sum.add(sq_sum_val);
			}

			// Printing of result
			PrintWriter trend_unique_summary = new PrintWriter(
					new File(version_baseDir, String.format(OPT_SUMMARY_TREND_UNIQUE_FILE, trendPt + 1)));

			trend_unique_summary.print("Param");
			for (int i = 0; i < num_param; i++) {
				trend_unique_summary.print(",");
			}
			trend_unique_summary.println("Sq_Sum");

			for (int p = 1; p < param_summary_unique.length; p++) {
				row = param_summary_unique[p].split(",");
				param_str = Arrays.copyOfRange(row, 4, row.length);
				param_key = Arrays.toString(param_str);

				trend_unique_summary.print(param_key.substring(1, param_key.length() - 1));
				ArrayList<Double> sq_sum = trendValueCollection.get(param_key);
				for (double sq_sum_val : sq_sum) {
					trend_unique_summary.print(",");
					trend_unique_summary.print(sq_sum_val);
				}
				trend_unique_summary.println();
			}
			trend_unique_summary.close();

		}
	}

	public static void combineOptTrendResults(ArrayList<File> baseDirsArr, File combinedBaseDir)
			throws FileNotFoundException, IOException {

		File combinedSummary = new File(combinedBaseDir, OPT_SUMMARY_FILE);
		File[] combinedTrend = new File[0];

		ArrayList<Number[]> summary_store = new ArrayList<>();
		HashMap<String, String> trend_mapping = new HashMap<>();
		String trend_mapping_key_pattern = "%d_%d"; // Id, trend number

		int mapping_key = 0;
		Pattern pattern_trend_file = Pattern.compile(OPT_SUMMARY_TREND_FILE.replaceAll("%d", "(\\\\d+)"));

		for (File baseDir : baseDirsArr) {
			File summaryFile = new File(baseDir, OPT_SUMMARY_FILE);
			File[] trendFiles = baseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pattern_trend_file.matcher(pathname.getName()).matches();
				}
			});
			Arrays.sort(trendFiles, new Comparator<File>() {
				@Override
				public int compare(File o1, File o2) {
					Matcher m1 = pattern_trend_file.matcher(o1.getName());
					Matcher m2 = pattern_trend_file.matcher(o2.getName());
					m1.matches();
					m2.matches();
					return Integer.compare(Integer.parseInt(m1.group(1)), Integer.parseInt(m2.group(1)));
				}

			});

			if (combinedTrend.length == 0) {
				combinedTrend = new File[trendFiles.length];
				for (int i = 0; i < combinedTrend.length; i++) {
					combinedTrend[i] = new File(combinedBaseDir, trendFiles[i].getName());
				}
			}

			String line;
			BufferedReader reader_summary = new BufferedReader(new FileReader(summaryFile));
			BufferedReader[] reader_trends = new BufferedReader[trendFiles.length];

			// Skip first line
			for (int t = 0; t < reader_trends.length; t++) {
				reader_trends[t] = new BufferedReader(new FileReader(trendFiles[t]));
				reader_trends[t].readLine();
			}
			reader_summary.readLine();

			while ((line = reader_summary.readLine()) != null) {
				String[] ent = line.split(",");
				Number[] summary_store_ent = new Number[ent.length + 1];
				for (int i = 0; i < ent.length; i++) {
					switch (i) {
					case 0:
						summary_store_ent[i] = Double.parseDouble(ent[i]);
						break;
					case 1:
					case 2:
						summary_store_ent[i] = Long.parseLong(ent[i]);
						break;
					default:
						try {
							summary_store_ent[i] = Float.parseFloat(ent[i]);
						} catch (NumberFormatException ex) {
							summary_store_ent[i] = Float.NaN;
						}
					}
				}
				summary_store_ent[ent.length] = mapping_key;

				for (int t = 0; t < reader_trends.length; t++) {
					String trend_mapping_key = String.format(trend_mapping_key_pattern, mapping_key, t);
					trend_mapping.put(trend_mapping_key, reader_trends[t].readLine());
				}
				summary_store.add(summary_store_ent);
				mapping_key++;

			}

			reader_summary.close();
			for (BufferedReader reader_trend : reader_trends) {
				reader_trend.close();
			}
		}

		// Sort results based on residue

		Collections.sort(summary_store, new Comparator<Number[]>() {
			@Override
			public int compare(Number[] o1, Number[] o2) {
				int res = 0;
				int pt = 0;
				while (res == 0 && pt < o1.length) {
					if (o1[pt] instanceof Long) {
						res = Long.compare(o1[pt].longValue(), o2[pt].longValue());
					} else {
						res = Double.compare(o1[pt].doubleValue(), o2[pt].doubleValue());
					}
					pt++;
				}
				return res;
			}
		});

		PrintWriter pri_summary = new PrintWriter(combinedSummary);
		PrintWriter[] pri_trends = new PrintWriter[combinedTrend.length];

		// Print heading
		pri_summary.println("Residue,CMapSeed,SimSeed,Param,...Key_id");
		for (int p = 0; p < pri_trends.length; p++) {
			pri_trends[p] = new PrintWriter(combinedTrend[p]);
			pri_trends[p].println("Opt_trend time-value pairing");
		}

		for (Number[] ent : summary_store) {
			for (int i = 0; i < ent.length; i++) {
				if (i > 0) {
					pri_summary.print(',');
				}
				pri_summary.print(ent[i].toString());
			}
			pri_summary.println();
			for (int p = 0; p < pri_trends.length; p++) {
				String trend_mapping_key = String.format(trend_mapping_key_pattern, ent[ent.length - 1], p);
				pri_trends[p].println(trend_mapping.get(trend_mapping_key));
			}
		}

		pri_summary.close();
		for (PrintWriter p : pri_trends) {
			p.close();
		}

	}

	public static void extractBestOptrendResults(File basedir, Pattern subDirPattern) throws IOException {

		final Pattern bestFileFile_pattern = Pattern
				.compile(OPT_TREND_FILE_NAME_TREND_OUTPUT.replaceAll("%d", "(-{0,1}\\\\d+)"));

		File[] subDirs;
		if (subDirPattern != null) {
			subDirs = basedir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isDirectory() && subDirPattern.matcher(pathname.getName()).matches();
				}
			});
		} else {
			subDirs = new File[] { basedir };
		}

		// double residue, long cMapSeed, long simSeed, long offset, String[] param_str,
		// String[] trends, File trendFile
		ArrayList<Object[]> resultsTrendCollections = new ArrayList<>();
		final Comparator<Object[]> resultTrendCollectionsComp = new Comparator<>() {
			@Override
			public int compare(Object[] o1, Object[] o2) {
				int r = 0;
				for (int p = 0; (r == 0) && (p < o1.length); p++) {
					if (o1[p] instanceof Double) {
						r = Double.compare((Double) o1[p], (Double) o2[p]);
					} else if (o1[p] instanceof Long) {
						r = Long.compare((Long) o1[p], (Long) o2[p]);

					} else if (o1[p] instanceof File) {
						r = ((File) o1[p]).compareTo((File) o2[p]);
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

			File[] trendFiles = subDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return bestFileFile_pattern.matcher(pathname.getName()).matches();
				}
			});
			for (File trendFile : trendFiles) {

				BufferedReader reader = new BufferedReader(new FileReader(trendFile));
				String line;

				Matcher m = bestFileFile_pattern.matcher(trendFile.getName());
				m.matches();

				long cMapSeed = Long.parseLong(m.group(1));
				long simSeed = Long.parseLong(m.group(2));

				while ((line = reader.readLine()) != null) {
					if (line.startsWith(Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM)) {
						String[] param_str = line
								.substring(Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM.length() + 1, line.length() - 1)
								.split(",");

						line = reader.readLine();
						double residue = Double
								.parseDouble(line.substring(Optimisation_Factory.OPT_OUTPUT_PREFIX_RESIDUE.length()));

						line = reader.readLine();
						long offset = 0;
						if (line.startsWith(OPT_TREND_OUTPUT_PREFIX_OFFSET)) {
							offset = Integer.parseInt(line.substring(OPT_TREND_OUTPUT_PREFIX_OFFSET.length()));
							line = reader.readLine();
						}

						ArrayList<String> trend_arr = new ArrayList<>();

						// TrendEntry
						while (line != null) {
							if (line.trim().length() > 0) {
								trend_arr.add(line);
							} else {
								break;
							}
							line = reader.readLine();
						}
						String[] trends = trend_arr.toArray(new String[trend_arr.size()]);
						Object[] val = new Object[] { residue, cMapSeed, simSeed, offset, param_str, trends,
								trendFile };

						int key = Collections.binarySearch(resultsTrendCollections, val, resultTrendCollectionsComp);
						if (key < 0) {
							resultsTrendCollections.add(~key, val);
						}
					}
				}
				reader.close();
			}
		}

		// Print of results
		PrintWriter pWri_summary = new PrintWriter(new File(basedir, OPT_SUMMARY_FILE));
		pWri_summary.println("Residue,CMapSeed,SimSeed,Dir,Param");

		PrintWriter pWri_summary_unique = new PrintWriter(new File(basedir, OPT_SUMMARY_UNIQUE_FILE));
		pWri_summary_unique.println("Residue,CMapSeed,SimSeed,Dir,Param");
		ArrayList<String> printed_entries = new ArrayList<>();

		PrintWriter[] pWri_trend = null;

		for (Object[] val : resultsTrendCollections) {
			pWri_summary.printf("%f,%d,%d,%s", val[0], val[1], val[2],
					((File) val[val.length - 1]).getParentFile().getName());
			String[] param = (String[]) val[4];
			for (String element : param) {
				pWri_summary.print(',');
				pWri_summary.print(element);
			}
			pWri_summary.println();

			String paramStr = Arrays.toString(param);
			int r = Collections.binarySearch(printed_entries, paramStr);
			if (r < 0) {
				pWri_summary_unique.printf("%f,%d,%d,%s", val[0], val[1], val[2],
						((File) val[val.length - 1]).getParentFile().getName());
				for (String element : param) {
					pWri_summary_unique.print(',');
					pWri_summary_unique.print(element);
				}
				pWri_summary_unique.println();
				printed_entries.add(~r, paramStr);
			}

			String[] trends = (String[]) val[5];
			int offset = ((Long) val[3]).intValue();

			if (pWri_trend == null) {
				String[] firstLine = trends[0].split(",");
				pWri_trend = new PrintWriter[firstLine.length];
				for (int w = 1; w < pWri_trend.length; w++) {
					pWri_trend[w] = new PrintWriter(new File(basedir, String.format(OPT_SUMMARY_TREND_FILE, w)));
					pWri_trend[w].println("Opt_trend time-value pairing");
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
				pWri_trend[v].print(timeline.toString());
			}

			for (String trend : trends) {
				String[] ent = trend.split(",");
				for (int v = 1; v < ent.length; v++) {
					pWri_trend[v].print(',');
					pWri_trend[v].print(ent[v]);
				}
			}

			for (int v = 1; v < pWri_trend.length; v++) {
				pWri_trend[v].println();
			}

		}

		pWri_summary.close();
		pWri_summary_unique.close();
		if (pWri_trend != null) {

			for (PrintWriter pWri : pWri_trend) {
				if (pWri != null) {
					pWri.close();
				}
			}
		}

		System.out.printf("Export opt trend result completed at %s.\n", basedir.getAbsolutePath());

	}

	public static HashMap<String, double[][]> loadTrendCSV(Properties prop) throws FileNotFoundException, IOException {

		if (prop.getProperty(Optimisation_Factory.POP_PROP_OPT_TARGET) != null) {
			String[] target_trend_csv_path = prop.getProperty(Optimisation_Factory.POP_PROP_OPT_TARGET)
					.replaceAll("\n", "").split(",");
			int pt = 0;
			File[] target_trend_csv = new File[target_trend_csv_path.length];
			for (String csv_path : target_trend_csv_path) {
				target_trend_csv[pt] = new File(csv_path);
				pt++;
			}
			return loadTrendCSV(target_trend_csv);
		} else {
			System.err.printf("%s not found in Properties.\n", Optimisation_Factory.POP_PROP_OPT_TARGET);
			return null;
		}

	}

	public static HashMap<String, double[][]> loadTrendCSV(File[] target_trent_csv_path)
			throws FileNotFoundException, IOException {
		HashMap<String, double[][]> target_trend_collection = new HashMap<>();

		double[] time_range = new double[] { Double.NaN, Double.NaN };
		for (File csv_file : target_trent_csv_path) {
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
				String mapKey = String.format(OPT_TREND_MAP_KEY_FORMAT, csv_file.getAbsolutePath(), type, tar_grp,
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
		return target_trend_collection;
	}

	@SuppressWarnings("unchecked")
	public static double[] calculate_residue_opt_trend(double[] point, HashMap<String, Object> args,
			HashMap<String, Object> outputMap, final int NUM_THREADS) {

		ContactMap[] cMapArr = (ContactMap[]) args.get(OptTrendFittingFunction.ARGS_CMAP);
		long[] cMap_seed = (long[]) args.get(OptTrendFittingFunction.ARGS_CMAP_SEED);
		long[] sim_seeds = (long[]) args.get(OptTrendFittingFunction.ARGS_SIM_SEED);
		int NUM_SIM_PER_MAP = sim_seeds.length;

		boolean verbose = args.containsKey(OptTrendFittingFunction.ARGS_VERBOSE) ? true : false;

		double[] bestResidue_by_runnable;
		bestResidue_by_runnable = new double[cMapArr.length * NUM_SIM_PER_MAP];
		Arrays.fill(bestResidue_by_runnable, Double.NaN);

		// From args
		File baseDir = (File) args.get(OptTrendFittingFunction.ARGS_BASEDIR);
		Properties prop = (Properties) args.get(OptTrendFittingFunction.ARGS_PROP);
		HashMap<String, double[][]> target_trend_collection = (HashMap<String, double[][]>) args
				.get(OptTrendFittingFunction.ARGS_TAR_TRENDS_COLLECTIONS);
		Number[][] result_lookup = (Number[][]) args.get(ARGS_PREV_RESULTS);

		// Setting initial value from properties
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

		int simSetting = 1;

		if (prop.containsKey(Simulation_ClusterModelTransmission.PROP_SIM_SETTING)) {
			simSetting = Integer.parseInt(prop.getProperty(Simulation_ClusterModelTransmission.PROP_SIM_SETTING));
		}

		Abstract_Runnable_ClusterModel_Transmission[] runnable = new Abstract_Runnable_ClusterModel_Transmission[bestResidue_by_runnable.length];

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

		String popType = (String) prop.get(SimulationInterface.PROP_NAME[SimulationInterface.PROP_POP_TYPE]);

		int multipTrans_index = -1;
		for (int i = 0; i < pattern_multi_trans.length || multipTrans_index < 0; i++) {
			if (pattern_multi_trans[i].matcher(popType).matches()) {
				multipTrans_index = i;
			}
		}

		HashMap<String, ArrayList<String[]>> seedListParameter = new HashMap<>();
		File seedFile = (File) args.get(OptTrendFittingFunction.ARGS_SEEDLIST);
		if (seedFile != null && seedFile.exists()) {
			try {
				BufferedReader reader = new BufferedReader(new FileReader(seedFile));
				String[] headerline = reader.readLine().split(","); // Header
				String line;
				while ((line = reader.readLine()) != null) {
					String[] line_sp = line.split(",");
					String key = String.format("%d,%d", Long.parseLong(line_sp[0]), Long.parseLong(line_sp[1]));
					ArrayList<String[]> ent = seedListParameter.get(key);
					if (ent == null) {
						ent = new ArrayList<>();
						ent.add(Arrays.copyOfRange(headerline, 2, headerline.length));
						seedListParameter.put(key, ent);
					}
					ent.add(Arrays.copyOfRange(line_sp, 2, headerline.length));
				}
				reader.close();
			} catch (IOException e) {
				System.err.printf("Warning: Error in reading seed file %s.", seedFile.getAbsolutePath());
			}

		}

		for (ContactMap cMap : cMapArr) {
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
						Number[] found_res = result_lookup[k];
						bestResidue_by_runnable[rId] = (double) found_res[found_res.length - 1];
						runnable[rId] = null;
					}
				}
				if (Double.isNaN(bestResidue_by_runnable[rId])) {

					if (verbose) {
						System.out.printf("[%d,%d]->%s: starting simulations.\n", cMap_seed[cMap_id], sim_seed,
								Arrays.toString(point));
					}

					if (multipTrans_index >= 0) {
						switch (multipTrans_index) {
						case pattern_multi_trans_index_viability:
							runnable[rId] = new Runnable_ClusterModel_Viability(cMap_seed[cMap_id], sim_seed, cMap, prop) {
								@Override
								protected void postSimulation() {
									// Do nothing
								}
							};
							break;
						case pattern_multi_trans_index_Bali:
							runnable[rId] = new Runnable_ClusterModel_Bali(cMap_seed[cMap_id], sim_seed,
									POP_COMPOSITION, cMap, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP) {
								@Override
								protected void postSimulation() {
									// Do nothing
								}						};
							break;
						case pattern_multi_trans_index_base:
						default:
							Matcher m = Runnable_ClusterModel_MultiTransmission.PROP_TYPE_PATTERN.matcher(popType);
							runnable[rId] = new Runnable_ClusterModel_MultiTransmission(cMap_seed[cMap_id], sim_seed,
									POP_COMPOSITION, cMap, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP, Integer.parseInt(m.group(1)),
									Integer.parseInt(m.group(2)), Integer.parseInt(m.group(3))) {
								@Override
								protected void postSimulation() {
									// Do nothing
								}
							};

						}

						if (runnable[rId] == null) {
							System.err.printf("Error: Null runnable for Optimisation of PROP_POP_TYPE of %s\n.",
									popType);
							System.exit(-1);
						}

					} else {
						runnable[rId] = new Runnable_ClusterModel_Transmission(cMap_seed[cMap_id], sim_seed,
								POP_COMPOSITION, cMap, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP);
					}

					runnable[rId].setBaseDir(baseDir);

					for (int i = Optimisation_Factory.RUNNABLE_OFFSET; i < Optimisation_Factory.RUNNABLE_OFFSET
							+ runnable[rId].getRunnable_fields().length; i++) {
						String key = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX + Integer.toString(i);
						if (prop.containsKey(key)) {
							runnable[rId].getRunnable_fields()[i - Optimisation_Factory.RUNNABLE_OFFSET] = PropValUtils
									.propStrToObject(prop.getProperty(key),
											runnable[rId].getRunnable_fields()[i - Optimisation_Factory.RUNNABLE_OFFSET]
													.getClass());
						}
					}

					runnable[rId].setSimSetting(simSetting);
					Optimisation_Factory.setOptParamInRunnable(runnable[rId], prop, point, SEED_INFECTION, cMap == null);

					ArrayList<String[]> seedListParam = seedListParameter
							.get(String.format("%d,%d", cMap_seed[cMap_id], sim_seed));

					if (seedListParam != null && seedListParam.size() > 1) {
						String[] param_def = seedListParam.get(0);
						String[] param_val_str = seedListParam.remove(1);
						double[] param_val = new double[param_val_str.length];
						for (int i = 0; i < param_val.length; i++) {
							param_val[i] = Double.parseDouble(param_val_str[i]);
						}
						runnable[rId].loadOptParameter(param_def, param_val, SEED_INFECTION, cMap == null);
					}

					runnable[rId].initialse();

					File riskGrpDir = baseDir;
					File pre_allocate_risk_file = new File(baseDir, String.format(
							Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, cMap_seed[cMap_id]));

					if (!pre_allocate_risk_file.exists()) {
						// Try loading one in cMap folder
						if (prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC) != null) {
							riskGrpDir = new File(
									prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC));
							if (!riskGrpDir.exists() || !riskGrpDir.isDirectory()) {
								riskGrpDir = baseDir;
							}

							pre_allocate_risk_file = new File(riskGrpDir,
									String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP,
											cMap_seed[cMap_id]));
						}
					}

					if (pre_allocate_risk_file.exists()) {
						ArrayList<Number[]> riskGrpArr = new ArrayList<>();
						try {
							Simulation_ClusterModelTransmission.loadPreallocateRiskGrp(riskGrpArr, riskGrpDir,
									cMap_seed[cMap_id]);
						} catch (NumberFormatException | IOException e) {
							e.printStackTrace(System.err);
						}
						runnable[rId].fillRiskCatMap(riskGrpArr);
					}

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

		if (multipTrans_index >= 0) {
			// Set up fitting target
			int num_target_trend = target_trend_collection.size();
			String[] trend_target_key = target_trend_collection.keySet().toArray(new String[num_target_trend]);
			Arrays.sort(trend_target_key);

			UnivariateInterpolator interpolator = new LinearInterpolator();
			int minFitFrom = 0;
			double[] weight = new double[num_target_trend];
			int[] first_target_time = new int[num_target_trend];
			int[] last_target_time = new int[num_target_trend];
			UnivariateFunction[] interpolation = new PolynomialSplineFunction[num_target_trend];
			String[][] trend_target_key_split = new String[num_target_trend][];
			Pattern pattern_trend_type_composite = Pattern.compile("\\[(.*) (.*) (\\d+)\\]");
			Pattern pattern_trend_type = Pattern.compile("(.*)_C(-?\\d+)");

			for (int trend_target_pt = 0; trend_target_pt < num_target_trend; trend_target_pt++) {
				trend_target_key_split[trend_target_pt] = trend_target_key[trend_target_pt].split(",");
				String[] trend_keys = trend_target_key_split[trend_target_pt];
				minFitFrom = Math.max(minFitFrom,
						Integer.parseInt(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_FITFROM]));
				double[][] tar_values = target_trend_collection.get(trend_target_key[trend_target_pt]);
				interpolation[trend_target_pt] = interpolator.interpolate(tar_values[0], tar_values[1]);
				first_target_time[trend_target_pt] = (int) tar_values[0][0];
				last_target_time[trend_target_pt] = (int) tar_values[0][tar_values[0].length - 1];
				weight[trend_target_pt] = Double
						.parseDouble(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_WEIGHT]);

			}

			for (int r = 0; r < rId; r++) {
				if (runnable[r] != null) {
					if (verbose) {
						System.out.printf("[%d,%d]->%s: completed.\n", runnable[r].getcMap_seed(),
								runnable[r].getSim_seed(), Arrays.toString(point));
					}
					bestResidue_by_runnable[r] = 0;
					StringBuilder[] str_disp = null;

					for (int trend_target_pt = 0; trend_target_pt < num_target_trend; trend_target_pt++) {
						String trend_type = trend_target_key_split[trend_target_pt][OptTrendFittingFunction.OPT_TREND_MAP_KEY_TYPE];

						Matcher trend_matcher_composite = pattern_trend_type_composite.matcher(trend_type);
						String[] trend_type_arr;
						int trend_operation = -1;

						if (trend_matcher_composite.matches()) {
							trend_type_arr = new String[] { trend_matcher_composite.group(1),
									trend_matcher_composite.group(2) };
							trend_operation = Integer.parseInt(trend_matcher_composite.group(3));
						} else {
							trend_type_arr = new String[] { trend_type };
						}

						Integer[] sim_time = null;
						double[][] model_val = new double[2][];

						for (String trend_key : trend_type_arr) {
							Matcher trend_matcher = pattern_trend_type.matcher(trend_key);
							if (trend_matcher.matches()) {
								HashMap<Integer, int[]> countMap = (HashMap<Integer, int[]>) runnable[r].getSim_output()
										.get(trend_matcher.group(1));
								int col_number = Integer.parseInt(trend_matcher.group(2));

								if (sim_time == null) {
									sim_time = countMap.keySet().toArray(new Integer[countMap.size()]);
									Arrays.sort(sim_time);
									if (str_disp == null) {
										str_disp = new StringBuilder[sim_time.length];
									}
									model_val = new double[2][sim_time.length];
									for (int mP = 0; mP < model_val.length; mP++) {
										model_val[mP] = new double[sim_time.length];
										Arrays.fill(model_val[mP], Double.NaN);
									}
								}
								for (int i = 0; i < str_disp.length; i++) {
									model_val[0][i] = sim_time[i];
									double single_model_val = countMap.get(sim_time[i])[Math.abs(col_number)];
									double offset = 0;
									if (col_number < 0) {
										int pt = Arrays.binarySearch(sim_time, first_target_time[trend_target_pt]);
										if (pt < 0) {
											System.err.printf(
													"Warning: offset time %d not found in countMap of time range = %s. Average of neighbouring value used.\n",
													first_target_time[trend_target_pt], Arrays.deepToString(sim_time));
											offset = ((countMap.get(sim_time[~pt - 1]))[Math.abs(col_number)]
													+ (countMap.get(sim_time[~pt]))[Math.abs(col_number)]) / 2;
										} else {
											offset = (countMap.get(sim_time[pt]))[Math.abs(col_number)];
										}
									}

									if (weight[trend_target_pt] < 0) {
										if (i > 0 && countMap.containsKey((sim_time[i - 1]))) {
											offset += countMap.get(sim_time[i - 1])[Math.abs(col_number)];
										} else {
											offset = Double.NaN;
										}
									}

									single_model_val -= offset;

									if (Double.isNaN(model_val[1][i])) {
										model_val[1][i] = single_model_val;
									} else {
										switch (trend_operation) {
										case 0:
											model_val[1][i] += single_model_val;
											break;
										case 1:
											model_val[1][i] -= single_model_val;
											break;
										case 2:
											model_val[1][i] *= single_model_val;
											break;
										case 3:
											model_val[1][i] /= single_model_val;
											break;
										default:
											model_val[1][i] = single_model_val;

										}
									}
								}

							} else {
								System.err.printf("Warning: Ill formed trend type %s, fitting ignored.\n", trend_type);
							}
						}

						try {
							// Remove N/A model pair
							int first_valid_row = 0;

							// Print output
							for (int i = 0; i < str_disp.length; i++) {
								if (str_disp[i] == null) {
									str_disp[i] = new StringBuilder();
									str_disp[i].append(model_val[0][i]);
								}
								str_disp[i].append(',');
								str_disp[i].append(model_val[1][i]);
								if (Double.isNaN(model_val[1][i])) {
									first_valid_row++;
								}
							}
							if (first_valid_row != 0) {
								model_val[0] = Arrays.copyOfRange(model_val[0], first_valid_row, model_val[0].length);
								model_val[1] = Arrays.copyOfRange(model_val[1], first_valid_row, model_val[1].length);
							}
							UnivariateFunction model_interpolation = interpolator.interpolate(model_val[0],
									model_val[1]);
							int sample_time = Math.max(minFitFrom, first_target_time[trend_target_pt]);
							while (sample_time <= last_target_time[trend_target_pt]) {
								double sim_y;
								if (sample_time >= model_val[0][0]
										&& sample_time <= model_val[0][model_val[0].length - 1]) {
									sim_y = model_interpolation.value(sample_time);
								} else {
									// Outside model range - assume to be zero
									sim_y = 0;
								}
								bestResidue_by_runnable[r] += Math.abs(weight[trend_target_pt])
										* Math.pow((sim_y) - interpolation[trend_target_pt].value(sample_time), 2);
								sample_time += AbstractIndividualInterface.ONE_YEAR_INT;
							}

						} catch (NullPointerException | ArrayIndexOutOfBoundsException ex) {
							System.err.printf("Warning: Exception encountered for trend type %s, fitting ignored.\n",
									trend_type);
							ex.printStackTrace(System.err);
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

						pWri.printf("%s%d\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_CMAP, runnable[r].getcMap_seed());
						pWri.printf("%s%d\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_SIMSEED,
								runnable[r].getSim_seed());
						pWri.printf("%s[%s]\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM, param_str.toString());
						pWri.printf("%s%f\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_RESIDUE,
								bestResidue_by_runnable[r]);

						if (str_disp != null) {
							for (StringBuilder s : str_disp) {
								if (s != null) {
									pWri.println(s.toString());
								}
							}
						}
						pWri.println();

						pWri.close();
						fWri.close();

						if (verbose) {
							System.out.printf("[%d,%d]->%s: result exported to %s.\n", runnable[r].getcMap_seed(),
									runnable[r].getSim_seed(), Arrays.toString(point), opt_output_file.getName());
						}
					} catch (IOException e) {
						e.printStackTrace(System.err);
					}
				}
			}
		} else {
			HashMap<Integer, int[][]>[] countMapBySite = new HashMap[OptTrendFittingFunction.LENGTH_OPT_TREND_COUNT_MAP_BY_SITE];
			HashMap<Integer, int[]>[] countMapByPerson = new HashMap[OptTrendFittingFunction.LENGTH_OPT_TREND_COUNT_MAP_BY_PERSON];

			for (int r = 0; r < rId; r++) {

				if (runnable[r] != null) {

					if (verbose) {
						System.out.printf("[%d,%d]->%s: completed.\n", runnable[r].getcMap_seed(),
								runnable[r].getSim_seed(), Arrays.toString(point));
					}

					countMapBySite[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_SITE] = (HashMap<Integer, int[][]>) runnable[r]
							.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT);
					countMapBySite[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_SITE] = (HashMap<Integer, int[][]>) runnable[r]
							.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_INCIDENCE);

					countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_INFECTION_COUNT_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
							.getSim_output()
							.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT_BY_PERSON);
					countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_INCIDENCE_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
							.getSim_output()
							.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_INCIDENCE_BY_PERSON);
					countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_POS_DX_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
							.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_POS_DX_BY_PERSON);
					countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_POS_DX_SOUGHT_BY_PERSON] = (HashMap<Integer, int[]>) runnable[r]
							.getSim_output()
							.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_POS_DX_SOUGHT_BY_PERSON);

					Integer[] simTime = null;

					bestResidue_by_runnable[r] = 0;

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

						bestResidue_by_runnable[r] = Double.NaN;
					} else {
						Arrays.sort(simTime);
						int num_target_trend = target_trend_collection.size();
						String[] trend_target_key = target_trend_collection.keySet()
								.toArray(new String[num_target_trend]);

						Arrays.sort(trend_target_key);
						double[] model_t = new double[simTime.length];
						for (int i = 0; i < simTime.length; i++) {
							model_t[i] = simTime[i];
						}

						// Trim simTime to valid starting point

						String[][] trend_target_key_split = new String[trend_target_key.length][];
						for (int trend_target_pt = 0; trend_target_pt < num_target_trend; trend_target_pt++) {
							trend_target_key_split[trend_target_pt] = trend_target_key[trend_target_pt].split(",");
						}

						double[][][] tar_values = new double[num_target_trend][][];
						double[] weight = new double[num_target_trend];
						int[] fitFrom = new int[num_target_trend];
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

							double[] model_y = new double[model_t.length];

							tar_values[trend_target_pt] = target_trend_collection
									.get(trend_target_key[trend_target_pt]);
							weight[trend_target_pt] = Double
									.parseDouble(trend_keys[OptTrendFittingFunction.OPT_TREND_MAP_KEY_WEIGHT]);

							fitFrom[trend_target_pt] = Integer.parseInt(
									trend_target_key_split[trend_target_pt][OptTrendFittingFunction.OPT_TREND_MAP_KEY_FITFROM]);

							// Trend value with fitFrom adjustment
							double[] target_t = Arrays.copyOf(tar_values[trend_target_pt][0],
									tar_values[trend_target_pt][0].length);
							double[] target_y = Arrays.copyOf(tar_values[trend_target_pt][1],
									tar_values[trend_target_pt][1].length);

							for (int tt = 0; tt < target_t.length; tt++) {
								target_t[tt] += fitFrom[trend_target_pt];
							}
							interpolation[trend_target_pt] = interpolator.interpolate(target_t, target_y);

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
										int[][] ent = OptTrendFittingFunction.OPT_TREND_INPUT_TYPE_NUMINF
												.equals(type_key)
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
								} else if (OptTrendFittingFunction.OPT_TREND_INPUT_TYPR_DX_TESTING.equals(type_key)) {
									// V= int[gender_positive_dx, gender_true_infection]
									int[][] ent_collection = new int[][] {
											countMapByPerson[OptTrendFittingFunction.OPT_TREND_COUNT_MAP_CUMUL_POS_DX_BY_PERSON]
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
								model_y[t_pt] = newVal;
								str_disp[t_pt].append(',');
								str_disp[t_pt].append(newVal);

							}

							double sampleTime = target_t[0];

							UnivariateFunction model_value_interpol = interpolator.interpolate(model_t, model_y);

							while (sampleTime <= target_t[target_t.length - 1]) {
								double trend_y_val = interpolation[trend_target_pt].value(sampleTime);
								double model_y_val = Double.POSITIVE_INFINITY;

								if (model_t[0] <= sampleTime && sampleTime <= model_t[model_t.length - 1]) {
									model_y_val = model_value_interpol.value(sampleTime);
									if (weight[trend_target_pt] < 0) {
										// Offset from previous time step
										if (sampleTime - AbstractIndividualInterface.ONE_YEAR_INT >= model_t[0]) {
											model_y_val = model_y_val - model_value_interpol
													.value(sampleTime - AbstractIndividualInterface.ONE_YEAR_INT);
										} else {
											System.err.printf(
													"Warning: Negative weight (%f) for prior to time %d not available.\n",
													weight[trend_target_pt], model_t[0]);
											model_y_val = Double.POSITIVE_INFINITY;
										}
									}
								} else {
									// Out of range
									model_y_val = 0;
								}

								bestResidue_by_runnable[r] += Math.abs(weight[trend_target_pt])
										* Math.pow(model_y_val - trend_y_val, 2);

								sampleTime += AbstractIndividualInterface.ONE_YEAR_INT;
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

							pWri.printf("%s%d\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_CMAP,
									runnable[r].getcMap_seed());
							pWri.printf("%s%d\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_SIMSEED,
									runnable[r].getSim_seed());
							pWri.printf("%s[%s]\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM, param_str.toString());
							pWri.printf("%s%f\n", Optimisation_Factory.OPT_OUTPUT_PREFIX_RESIDUE,
									bestResidue_by_runnable[r]);

							for (StringBuilder s : str_disp) {
								pWri.println(s.toString());
							}

							pWri.println();

							pWri.close();
							fWri.close();

							if (verbose) {
								System.out.printf("[%d,%d]->%s: result exported to %s.\n", runnable[r].getcMap_seed(),
										runnable[r].getSim_seed(), Arrays.toString(point), opt_output_file.getName());
							}
						} catch (IOException e) {
							e.printStackTrace(System.err);
						}
					} // if (simTime != null) {
				} // if(runnable[r] != null) {
			} // for (int r = 0; r < rId; r++) {
			outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_COUNT_BY_SITE, countMapBySite);
			outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_COUNT_BY_PERSON, countMapByPerson);
		} // if(isMultiTrans){...} else{ }

		outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_RUNNABLE, runnable);
		outputMap.put(OptTrendFittingFunction.OPT_TREND_OUTPUT_BEST_RESIDUE, bestResidue_by_runnable);

		return bestResidue_by_runnable;
	}

}