package optimisation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.InvalidPropertiesFormatException;
import java.util.Properties;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import sim.Abstract_Runnable_ClusterModel;
import sim.Abstract_Runnable_ClusterModel_Transmission;
import sim.Runnable_ClusterModel_ContactMap_Generation;
import sim.Runnable_ClusterModel_Transmission;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelGeneration;
import sim.Simulation_ClusterModelTransmission;
import smile.math.kernel.GaussianKernel;
import smile.regression.GaussianProcessRegression;
import smile.stat.distribution.GaussianDistribution;
import util.ArrayUtilsRandomGenerator;
import util.PropValUtils;

public class Optimisation_Factory {

	public static final String OPT_OUTPUT_PREFIX_CMAP = "CMAP    = ";
	public static final String OPT_OUTPUT_PREFIX_SIMSEED = "SimSeed = ";
	public static final String OPT_OUTPUT_PREFIX_PARAM = "Param   = ";
	public static final String OPT_OUTPUT_PREFIX_RESIDUE = "Residue = ";

	static class OptTrendFittingCallable_FS implements Callable<HashMap<String, Object>> {

		HashMap<String, Object> opt_outputs;
		HashMap<String, Object> args;

		public static final String OPT_TREND_CALLABLE_OUTPUT_RESULT_KEY_FORMAT = "%d_%d";

		File baseDir;
		int num_eval = 100;
		double[] init_param;
		double[][] boundaries;
		int optMethod;
		int contact_map_start_time = 365;

		double best_so_far = Double.POSITIVE_INFINITY;

		public OptTrendFittingCallable_FS(HashMap<String, Object> args) {
			this.args = args;
			opt_outputs = new HashMap<>();
			try {
				this.baseDir = (File) args.get(OptTrendFittingFunction.ARGS_BASEDIR);
				this.num_eval = (int) args.get(OptTrendFittingFunction.ARGS_NUM_EVAL);
				this.init_param = (double[]) args.get(OptTrendFittingFunction.ARGS_INIT_PARAM);
				this.boundaries = (double[][]) args.get(OptTrendFittingFunction.ARGS_BOUNDARIES);
				this.optMethod = (int) args.get(OptTrendFittingFunction.ARGS_OPT_METHOD);

			} catch (ClassCastException e) {
				e.printStackTrace(System.err);
			}
		}

		@Override
		public HashMap<String, Object> call() throws Exception {
			long[] cMap_seed = (long[]) args.get(OptTrendFittingFunction.ARGS_CMAP_SEED);
			long[] sim_seed = (long[]) args.get(OptTrendFittingFunction.ARGS_SIM_SEED);

			opt_outputs.put(OptTrendFittingFunction.OPT_TREND_CALLABLE_OUTPUT_RESULT_KEY,
					String.format(OPT_TREND_CALLABLE_OUTPUT_RESULT_KEY_FORMAT, cMap_seed[0], sim_seed[0]));

			OptFittingFunction func = new OptFittingFunction() {

				@Override
				public long[] getCMap_seeds() {
					return cMap_seed;
				}

				@Override
				public long[] getSim_seeds() {
					return sim_seed;
				}

				@Override
				public double[] getBestResidue_by_runnable() {
					return (double[]) opt_outputs.get(OptTrendFittingFunction.OPT_TREND_OUTPUT_BEST_RESIDUE);
				}

				@Override
				public Properties getProperties() {
					return (Properties) args.get(OptTrendFittingFunction.ARGS_PROP);
				}

				@Override
				public double value(double[] point) {
					double[] bestResidue = OptTrendFittingFunction.calculate_residue_opt_trend(point, args, opt_outputs,
							1);
					if (bestResidue[0] < best_so_far) {
						best_so_far = bestResidue[0];
						opt_outputs.put(OptTrendFittingFunction.OPT_TREND_CALLABLE_OUTPUT_BEST_SO_FAR, best_so_far);
					}
					return bestResidue[0];
				}

			};

			switch (optMethod) {
			case OPT_METHOD_BAYESIAN_FS:
				runBayesianOpt(func, init_param, boundaries, num_eval, baseDir, sim_seed[0],
						(boolean) args.get(OptTrendFittingFunction.ARGS_PROGRESS_DISP),
						String.format("%d_%d_", cMap_seed[0], sim_seed[0]));
				break;
			default:
				runSimplex(func, init_param, boundaries, contact_map_start_time);
			}

			return opt_outputs;
		}

	}

	public static final class Comparator_ResultLookup implements Comparator<Number[]> {
		int[] range;
		final double PARAM_TOL = 0.0000005;

		public Comparator_ResultLookup(int[] range) {
			this.range = range;
		}

		@Override
		public int compare(Number[] o1, Number[] o2) {
			// CMap_Seed, Sim_Seed, Parameters ... , Residue
			int r = 0;
			int pt = range[0];

			int maxLength = Math.min(o1.length, o2.length);
			if (range[1] < 0) {
				maxLength = maxLength + range[1];
			}
			while (r == 0 && pt < maxLength) {
				if (o1[pt] instanceof Long) {
					r = Long.compare((Long) o1[pt], (Long) o2[pt]);
				} else {
					if (Math.abs((Double) o1[pt] - (Double) o2[pt]) < PARAM_TOL) {
						r = 0;
					} else {
						r = Double.compare((Double) o1[pt], (Double) o2[pt]);
					}
				}
				pt++;
			}
			return r;
		}

	}

	public static final Comparator_ResultLookup COMPARATOR_RESULT_LOOKUP_ALL = new Comparator_ResultLookup(
			new int[] { 0, Integer.MAX_VALUE });
	public static final Comparator_ResultLookup COMPARATOR_RESULT_LOOKUP_PARAM = new Comparator_ResultLookup(
			new int[] { 0, -1 });

	private static final String POP_PROP_INIT_PREFIX = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX;
	private static final String FILENAME_FORMAT_ALL_CMAP = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP;
	private static final String CMAP_REGEX_STR = FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");

	public static final String POP_PROP_OPT_TARGET = "POP_PROP_OPT_TARGET";

	// Format (for stable fit):
	// float[][] opt_target = new float[NUM_DATA_TO_FIT]
	// { DATA_FITTING_TYPE, GROUPS_TO_INCLUDE, OPT_WEIGHTING, OPT_TARGET_VALUES_0,
	// ...}

	private static final int OPT_TARGET_STABLE_FITTING_TYPE_NUM_INFECTED_BY_SITE = 0;
	private static final int OPT_TARGET_STABLE_FITTING_TYPE_NOTIFICATIONS_BY_PERSON = 1;

	private static final int OPT_TARGET_STABLE_DATA_FITTING_TYPE = 0;
	private static final int OPT_TARGET_STABLE_GROUPS_TO_INCLUDE = OPT_TARGET_STABLE_DATA_FITTING_TYPE + 1;
	private static final int OPT_TARGET_STABLE_OPT_WEIGHTING = OPT_TARGET_STABLE_GROUPS_TO_INCLUDE + 1;
	private static final int OPT_TARGET_STABLE_TARGET_VALUES = OPT_TARGET_STABLE_OPT_WEIGHTING + 1;

	public static final String PROP_SEED_INFECTION = POP_PROP_INIT_PREFIX + Integer.toString(
			Population_Bridging.LENGTH_FIELDS_BRIDGING_POP + Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
					+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
					+ Simulation_ClusterModelTransmission.SIM_FIELD_SEED_INFECTION); // "POP_PROP_INIT_PREFIX_14";

	public static final int RUNNABLE_OFFSET = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
			+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
			+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD;

	final static String FILENAME_OPT_RESULT = "Opt_res.txt";

	public static final int OPT_METHOD_SIMPLEX = 0;
	public static final int OPT_METHOD_BAYESIAN = OPT_METHOD_SIMPLEX + 1;
	public static final int OPT_METHOD_SIMPLEX_FS = OPT_METHOD_BAYESIAN + 1;
	public static final int OPT_METHOD_BAYESIAN_FS = OPT_METHOD_SIMPLEX_FS + 1;

	public final static String FILENAME_OPT_BAYESIAN_OBS = "Opt_Bayesian_Obs.csv";
	public final static String FILENAME_OPT_BAYESIAN_OBS_BOUNDED = "Opt_Bayesian_Obs_Bounded.csv";

	private static final Pattern PATTERN_OPT_TREND_TXT = Pattern
			.compile(OptTrendFittingFunction.OPT_TREND_FILE_NAME_TREND_OUTPUT.replaceAll("%d", "(-{0,1}\\\\d+)"));

	public static void trend_fit_Simplex(String[] args) throws FileNotFoundException, IOException {
		trend_fit_general(args, OPT_METHOD_SIMPLEX);
	}

	public static void trend_fit_Bayesian(String[] args) throws FileNotFoundException, IOException {
		trend_fit_general(args, OPT_METHOD_BAYESIAN);
	}

	public static void trend_fit_Simplex_fs(String[] args) throws FileNotFoundException, IOException {
		trend_fit_general(args, OPT_METHOD_SIMPLEX_FS);
	}

	public static void trend_fit_Bayesian_fs(String[] args) throws FileNotFoundException, IOException {
		trend_fit_general(args, OPT_METHOD_BAYESIAN_FS);
	}

	private static void trend_fit_general(String[] args, int optMethod) throws FileNotFoundException, IOException {

		int minArgLength = 4;
		final String FLAG_nEval = "-nEval=";
		final String FLAG_verbose = "-verbose";
		final String FLAG_forceInit = "-forceInit";

		final String USAGE_INFO = String.format(
				"Usage: PROP_FILE_DIRECTORY INIT_PARAM_VALUE (double[]) BOUNDARIES (double[][]) "
						+ "RESULT_LIST_FILENAME <optional: %sNUM_EVAL (int)>  <optional: %s)> <optional %s>",
				FLAG_nEval, FLAG_verbose, FLAG_forceInit);

		if (args.length < minArgLength) {
			System.out.println(USAGE_INFO);
			System.exit(0);
		}
		int numEval = 100;
		boolean verbose = false;
		boolean forceInit = false;

		// Read input argument
		File baseDir = new File(args[0]);
		double[] init_param_default = (double[]) PropValUtils.propStrToObject(args[1], double[].class);
		double[][] boundaries = (double[][]) PropValUtils.propStrToObject(args[2], double[][].class);

		// Usage: PROP_FILE_DIRECTORY INIT_PARAM_VALUE (double[]) BOUNDARIES
		// (double[][]) RESULT_LIST_FILENAME
		// <optional: -nEval=NUM_EVAL (int)> <optional: -verbose)>";

		for (int a = minArgLength; a < args.length; a++) {
			if (args[a].startsWith(FLAG_nEval)) {
				numEval = (int) Integer.parseInt(args[a].substring(FLAG_nEval.length()));
			}
			verbose |= FLAG_verbose.equals(args[a]);
			forceInit |= FLAG_forceInit.equals(args[a]);
		}

		File resList = new File(baseDir, args[3]);

		File propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);

		if (propFile.exists()) {

			FileInputStream fIS = new FileInputStream(propFile);
			Properties prop = new Properties();
			prop.loadFromXML(fIS);
			fIS.close();

			System.out.println(String.format("Properties file < %s > loaded.", propFile.getAbsolutePath()));

			long seed = System.currentTimeMillis();
			int numThreads = Runtime.getRuntime().availableProcessors();
			int numSimPerMap = 1;

			// Load trend CSV
			// Key: Path,type,tar_grp,weight
			HashMap<String, double[][]> target_trend_collection = OptTrendFittingFunction.loadTrendCSV(prop);

			if (target_trend_collection.size() > 0) {
				// ContactMap Dir
				File contactMapDir = baseDir;
				if (prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC) != null) {
					contactMapDir = new File(
							prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC));
					if (!contactMapDir.exists() || !contactMapDir.isDirectory()) {
						contactMapDir = baseDir;
					}
				}

				if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
					seed = Long.parseLong(
							prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
				}

				if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL])) {
					numThreads = Integer.parseInt(
							prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL]));
				}

				if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET])) {
					numSimPerMap = Integer.parseInt(
							prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET]));
				}

				File[] preGenClusterFiles;
				ArrayList<Long> sim_seeds = new ArrayList<>();
				ArrayList<Long> cMap_seeds = new ArrayList<>();

				if (resList.exists()) {
					BufferedReader reader = new BufferedReader(new FileReader(resList));
					reader.readLine(); // Header line

					String line;
					while ((line = reader.readLine()) != null) {
						String[] ent = line.split(",");
						if (ent.length > 2) {
							if (ent[0].length() > 0) {
								cMap_seeds.add(Long.parseLong(ent[0]));
							}
							if (ent[1].length() > 0) {
								sim_seeds.add(Long.parseLong(ent[1]));
							}
						}
					}
					reader.close();

					Long[] cMap_seeds_arr = cMap_seeds.toArray(new Long[cMap_seeds.size()]);
					Arrays.sort(cMap_seeds_arr);

					preGenClusterFiles = contactMapDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							Matcher m = Pattern.compile(CMAP_REGEX_STR).matcher(pathname.getName());
							int inList = -1;
							if (m.matches()) {
								Long ent = Long.parseLong(m.group(1));
								inList = Arrays.binarySearch(cMap_seeds_arr, ent);
							}

							return pathname.isFile() && inList >= 0;

						}
					});

				} else {
					preGenClusterFiles = contactMapDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pathname.isFile() && Pattern.matches(CMAP_REGEX_STR, pathname.getName());

						}
					});
				}

				long tic = System.currentTimeMillis();

				ContactMap[] baseCMaps = new ContactMap[preGenClusterFiles.length];
				long[] baseCMapSeeds = new long[baseCMaps.length];
				RandomGenerator rng = new MersenneTwisterRandomGenerator(seed);
				double[][] target_trend_time_range = target_trend_collection
						.remove(OptTrendFittingFunction.OPT_TREND_CSV_RANGE);

				final String popCompositionKey = POP_PROP_INIT_PREFIX
						+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
				int[] pop_composition = (int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey),
						int[].class);
				int[] cumulative_pop_composition = new int[pop_composition.length];
				int pop_offset = 0;
				for (int g = 0; g < cumulative_pop_composition.length; g++) {
					cumulative_pop_composition[g] = pop_offset + pop_composition[g];
					pop_offset += pop_composition[g];
				}

				final String riskCatListAllKey = POP_PROP_INIT_PREFIX
						+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
								+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
								+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
								+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD
								+ Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS);
				float[][] riskCatListAll = (float[][]) PropValUtils.propStrToObject(prop.getProperty(riskCatListAllKey),
						float[][].class);

				final String time_rangeKey = POP_PROP_INIT_PREFIX
						+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
								+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
								+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);

				int[] map_time_range = (int[]) PropValUtils.propStrToObject(prop.getProperty(time_rangeKey),
						int[].class);

				int cMap_count = extractContactMap(baseCMaps, baseCMapSeeds, preGenClusterFiles,
						Math.min(numThreads, preGenClusterFiles.length));

				if (cMap_count > 0) {
					System.out.printf("%d ContactMap(s) from %s loaded. Time req. = %.3fs\n", cMap_count,
							contactMapDir.getAbsolutePath(), (System.currentTimeMillis() - tic) / 1000f);
				} else {
					System.out.println("No contact map defined. Printing out runnable fields only;");
					Abstract_Runnable_ClusterModel_Transmission runnable = new Runnable_ClusterModel_Transmission(0, 0,
							pop_composition, null, 0, 0);

					for (int i = Optimisation_Factory.RUNNABLE_OFFSET; i < Optimisation_Factory.RUNNABLE_OFFSET
							+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; i++) {

						String key = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX + Integer.toString(i);
						if (prop.containsKey(key)) {
							runnable.getRunnable_fields()[i - Optimisation_Factory.RUNNABLE_OFFSET] = PropValUtils
									.propStrToObject(prop.getProperty(key),
											runnable.getRunnable_fields()[i - Optimisation_Factory.RUNNABLE_OFFSET]
													.getClass());
						}
					}

					setOptParamInRunnable(runnable, prop, init_param_default, true);
					System.exit(1);

				}

				// Generate pre_allocated prealloactedRiskGrpArr

				for (int c = 0; c < baseCMaps.length; c++) {
					ArrayList<Number[]> riskGrpArr = new ArrayList<>();
					File pre_allocate_risk_file = new File(baseDir, String.format(
							Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, baseCMapSeeds[c]));
					if (!pre_allocate_risk_file.exists()) {
						Simulation_ClusterModelTransmission.fillRiskGrpArrByCasualPartnership(riskGrpArr, baseCMaps[c],
								cumulative_pop_composition, riskCatListAll, map_time_range);

						Simulation_ClusterModelTransmission.reallocateRiskGrp(riskGrpArr, baseCMapSeeds[c],
								cumulative_pop_composition, riskCatListAll, baseDir, seed);

					}
				}

				// Loading of previous optTrend text output files
				int num_param = init_param_default.length;
				Number[][] results_lookup = extractPreviousOptTrend(baseDir, num_param);

				prop.put(OptTrendFittingFunction.ARGS_PREV_RESULTS, results_lookup);
				if (verbose) {
					prop.put(OptTrendFittingFunction.ARGS_VERBOSE, true);
				}
				double[] init_param;

				switch (optMethod) {
				case OPT_METHOD_SIMPLEX:
				case OPT_METHOD_BAYESIAN:
					init_param = Arrays.copyOf(init_param_default, init_param_default.length);
					if (!forceInit && (results_lookup != null && results_lookup.length > 0)) {
						Arrays.sort(results_lookup, new Comparator<Number[]>() {
							@Override
							public int compare(Number[] o1, Number[] o2) {
								return Double.compare((Double) o1[o1.length], (Double) o2[o2.length]);
							}
						});
						for (int i = 0; i < init_param.length; i++) {
							init_param[i] = (double) results_lookup[0][i + 2];
						}
					}

					OptFittingFunction opt_trend_obj_func;

					if (resList.exists()) {
						HashMap<Long, ContactMap> cMap_lookup = new HashMap<>();
						for (int i = 0; i < baseCMapSeeds.length; i++) {
							cMap_lookup.put(baseCMapSeeds[i], baseCMaps[i]);
						}

						int cMapPt = 0;
						long[] cMapSeeds_long = new long[cMap_seeds.size()];
						ContactMap[] cMaps = new ContactMap[cMap_seeds.size()];

						for (Long cMapSeed : cMap_seeds) {
							if (cMap_lookup.containsKey(cMapSeed)) {
								cMapSeeds_long[cMapPt] = cMapSeed;
								cMaps[cMapPt] = cMap_lookup.get(cMapSeed);
								cMapPt++;
							}
						}

						long[] sim_seeds_long = new long[sim_seeds.size()];
						int sSeedPt = 0;
						for (Long ss : sim_seeds) {
							sim_seeds_long[sSeedPt] = ss;
							sSeedPt++;
						}
						opt_trend_obj_func = new OptTrendFittingFunction(baseDir, prop, Arrays.copyOf(cMaps, cMapPt),
								Arrays.copyOf(cMapSeeds_long, cMapPt), sim_seeds_long, target_trend_collection,
								target_trend_time_range, numThreads);

					} else {
						opt_trend_obj_func = new OptTrendFittingFunction(baseDir, prop, baseCMaps, baseCMapSeeds,
								numSimPerMap, rng, target_trend_collection, target_trend_time_range, numThreads);
						long[] gen_sim_seed = opt_trend_obj_func.getSim_seeds();

						PrintWriter pWri = new PrintWriter(resList);
						pWri.println("CMAP_SEED,SIM_SEED");

						for (int i = 0; i < Integer.max(baseCMapSeeds.length, gen_sim_seed.length); i++) {
							if (i < baseCMapSeeds.length) {
								pWri.print(baseCMapSeeds[i]);
							}
							pWri.print(',');
							if (i < gen_sim_seed.length) {
								pWri.print(gen_sim_seed[i]);
							}
						}
						pWri.close();
					}

					switch (optMethod) {
					case OPT_METHOD_SIMPLEX:
						runSimplex(opt_trend_obj_func, init_param, boundaries, numEval);
						break;
					case OPT_METHOD_BAYESIAN:
						runBayesianOpt(opt_trend_obj_func, init_param, boundaries, numEval, baseDir, seed,
								numThreads == 1);
						break;
					default:
						System.err.printf("OptMethod #%d not implemented.", optMethod);
					}
					break;

				case OPT_METHOD_BAYESIAN_FS:
				case OPT_METHOD_SIMPLEX_FS:
					OptTrendFittingCallable_FS[] opt_callable = new OptTrendFittingCallable_FS[baseCMaps.length
							* numSimPerMap];

					Number[][] best_result_collection = new Number[opt_callable.length][];

					int cId = 0;

					// Import result if exist
					if (resList.exists()) {
						BufferedReader reader = new BufferedReader(new FileReader(resList));
						reader.readLine(); // Header line
						String line;
						while ((line = reader.readLine()) != null) {
							String[] ent = line.split(",");
							if (ent.length > 3) {
								best_result_collection[cId] = new Number[ent.length];
								best_result_collection[cId][0] = Long.parseLong(ent[0]); // CMap_Seed
								best_result_collection[cId][1] = Long.parseLong(ent[1]); // Sim_Seed
								for (int c = 2; c < ent.length; c++) {
									best_result_collection[cId][c] = Double.parseDouble(ent[c]); // Best fit parameter
								}
								best_result_collection[cId][ent.length - 1] = Double.parseDouble(ent[ent.length - 1]); // Residue;
								cId++;
							}
						}
						reader.close();
					} else {
						// New result list
						while (cId < opt_callable.length) {
							for (int mapId = 0; mapId < baseCMaps.length; mapId++) {
								long cMap_seed = baseCMapSeeds[mapId];
								for (int s = 0; s < numSimPerMap; s++) {
									best_result_collection[cId] = new Number[3 + init_param_default.length];
									best_result_collection[cId][0] = cMap_seed;
									best_result_collection[cId][1] = rng.nextLong();
									best_result_collection[cId][best_result_collection[cId].length - 1] = Double.NaN;
								}
								cId++;
							}
						}
					}

					// Fill results with result lookup (if needed)
					if (results_lookup != null) {
						for (Number[] best_res : best_result_collection) {
							Number[] searchkey_min = new Number[best_res.length];
							Number[] searchkey_max = new Number[best_res.length];
							Arrays.fill(searchkey_min, Double.NEGATIVE_INFINITY);
							Arrays.fill(searchkey_max, Double.POSITIVE_INFINITY);
							searchkey_min[0] = best_res[0];
							searchkey_min[1] = best_res[1];
							searchkey_max[0] = best_res[0];
							searchkey_max[1] = best_res[1];

							int k_min = Arrays.binarySearch(results_lookup, searchkey_min,
									COMPARATOR_RESULT_LOOKUP_ALL);
							int k_max = Arrays.binarySearch(results_lookup, searchkey_max,
									COMPARATOR_RESULT_LOOKUP_ALL);

							if (k_min < 0) {
								k_min = ~k_min;
							}
							if (k_max < 0) {
								k_max = ~k_max;
							}

							for (int check_k = k_min; check_k < k_max; check_k++) {
								double best_res_residue_so_far = best_res[best_res.length - 1].doubleValue();
								if (results_lookup[check_k][results_lookup[check_k].length - 1]
										.doubleValue() < best_res_residue_so_far) {
									for (int i = 2; i < best_res.length; i++) {
										best_res[i] = results_lookup[check_k][i];
									}
								}
							}
						}
					}

					// Export current result collection
					exportResultCollection(best_result_collection, resList);
					cId = 0;

					Arrays.sort(best_result_collection, new Comparator<Number[]>() {
						@Override
						public int compare(Number[] o1, Number[] o2) {
							int r = -Double.compare((Double) o1[o1.length - 1], (Double) o2[o2.length - 1]);
							if (r == 0) {
								return COMPARATOR_RESULT_LOOKUP_ALL.compare(o1, o2);
							}
							return r;
						}
					});

					for (Number[] row : best_result_collection) {
						long cMap_seed = row[0].longValue();
						int cMap_index = Arrays.binarySearch(baseCMapSeeds, cMap_seed);
						if (cMap_index >= 0) {
							ContactMap cMap = baseCMaps[cMap_index];
							long sim_seed = row[1].longValue();
							// Start optimisation from previous min value
							init_param = Arrays.copyOf(init_param_default, init_param_default.length);

							if (row.length == init_param.length) {
								for (int i = 0; i < init_param.length; i++) {
									if (row[i + 2] != null && !((Double) row[i + 2]).isNaN()) {
										init_param[i] = ((Double) row[i + 2]).doubleValue();
									}
								}
							} else {
								System.out.printf(
										"Previous parameter set of length %d != %d required. Previous parameter set not loaded.\n",
										row.length, init_param.length);
							}
							HashMap<String, Object> arg = new HashMap<>();
							arg.put(OptTrendFittingFunction.ARGS_CMAP, new ContactMap[] { cMap });
							arg.put(OptTrendFittingFunction.ARGS_CMAP_SEED, new long[] { cMap_seed });
							arg.put(OptTrendFittingFunction.ARGS_SIM_SEED, new long[] { sim_seed });
							arg.put(OptTrendFittingFunction.ARGS_BASEDIR, baseDir);
							arg.put(OptTrendFittingFunction.ARGS_NUM_EVAL, numEval);

							arg.put(OptTrendFittingFunction.ARGS_INIT_PARAM, init_param);
							arg.put(OptTrendFittingFunction.ARGS_BOUNDARIES, boundaries);
							arg.put(OptTrendFittingFunction.ARGS_TAR_TRENDS_COLLECTIONS, target_trend_collection);
							arg.put(OptTrendFittingFunction.ARGS_TAR_TRENDS_TIMERANGE, target_trend_time_range);
							arg.put(OptTrendFittingFunction.ARGS_PROP, prop);
							arg.put(OptTrendFittingFunction.ARGS_OPT_METHOD, optMethod);
							arg.put(OptTrendFittingFunction.ARGS_PROGRESS_DISP, numThreads <= 1);
							arg.put(OptTrendFittingFunction.ARGS_PREV_RESULTS, results_lookup);
							if (verbose) {
								arg.put(OptTrendFittingFunction.ARGS_VERBOSE, true);
							}

							opt_callable[cId] = new OptTrendFittingCallable_FS(arg);
							cId++;
						}

					}

					cId = 0;
					@SuppressWarnings("unchecked")
					HashMap<String, Object>[] opt_outputs = new HashMap[opt_callable.length];

					if (numThreads <= 1 || opt_callable.length == 1) {
						while (cId < opt_callable.length) {
							if (opt_callable[cId] != null) {
								try {
									opt_outputs[cId] = opt_callable[cId].call();
									best_result_collection[cId][0] = (Double) opt_outputs[cId]
											.get(OptTrendFittingFunction.OPT_TREND_CALLABLE_OUTPUT_BEST_SO_FAR);
									exportResultCollection(best_result_collection, resList);
								} catch (Exception e) {
									e.printStackTrace(System.err);
								}
							}
							cId++;
						}
					} else {
						ExecutorService exec = null;
						@SuppressWarnings("unchecked")
						Future<HashMap<String, Object>>[] output_future = new Future[opt_callable.length];
						int inPoolFrom = 0;

						while (cId < opt_callable.length) {
							if (exec == null) {
								exec = Executors.newFixedThreadPool(numThreads);
							}
							if (opt_callable[cId] != null) {
								output_future[cId] = exec.submit(opt_callable[cId]);
								inPoolFrom = cId;
							}
							cId++;

							if ((cId - inPoolFrom) == numThreads) {
								executeOptTrendFittingCallable(exec, opt_callable, output_future,
										best_result_collection, resList, inPoolFrom, cId);
								exec = null;
							}
						}
						if (exec != null) {
							executeOptTrendFittingCallable(exec, opt_callable, output_future, best_result_collection,
									resList, inPoolFrom, cId);
							exec = null;

						}
					}
					// Handle opt_outputs

					break;
				default:
					System.err.printf("OptMethod #%d not implemented.", optMethod);
				}
			} else { // End of if (target_trend_collection.size() > 0)
				System.out.println("Target trend CSV not load (e.g. wrong path?).");
			}
		} else { // End of propFile.exists()
			System.out.printf("Properties file < %s > NOT found.\n", propFile.getAbsolutePath());
		}

	}

	private static Number[][] extractPreviousOptTrend(File baseDir, int num_param_input)
			throws FileNotFoundException, IOException {
		// Add opt_trend from previous collection
		File[] optTrend_txt_files = baseDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return PATTERN_OPT_TREND_TXT.matcher(pathname.getName()).matches();
			}
		});

		ArrayList<Number[]> results_lookup_arr = new ArrayList<>();

		int num_param = num_param_input;

		for (File optTrend_txt_file : optTrend_txt_files) {
			BufferedReader reader = new BufferedReader(new FileReader(optTrend_txt_file));
			String line;

			if (num_param < 0) {
				while ((line = reader.readLine()) != null && num_param > 0) {
					if (line.startsWith(Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM)) {
						num_param = (line.substring(Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM.length() + 1,
								line.length() - 1)).split(",").length;
						break;
					}
				}
				reader.close();
				reader = new BufferedReader(new FileReader(optTrend_txt_file));
			}

			while ((line = reader.readLine()) != null && num_param > 0) {
				if (line.startsWith(Optimisation_Factory.OPT_OUTPUT_PREFIX_CMAP)) {
					// cMap_Seed, sim_seed, unbounded_parameters, residue
					Number[] val = new Number[num_param + 3];
					int pt = 0;
					val[pt] = Long.parseLong(line.substring(Optimisation_Factory.OPT_OUTPUT_PREFIX_CMAP.length()));
					pt++;
					line = reader.readLine();
					val[pt] = Long.parseLong(line.substring(Optimisation_Factory.OPT_OUTPUT_PREFIX_SIMSEED.length()));
					pt++;
					line = reader.readLine();
					String[] param_ent = (line.substring(Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM.length() + 1,
							line.length() - 1)).split(",");
					for (int p = 0; p < param_ent.length; p++) {
						val[pt] = Double.parseDouble(param_ent[p]);
						pt++;
					}

					line = reader.readLine();
					val[pt] = Double
							.parseDouble(line.substring(Optimisation_Factory.OPT_OUTPUT_PREFIX_RESIDUE.length()));
					pt++;

					int key = Collections.binarySearch(results_lookup_arr, val, COMPARATOR_RESULT_LOOKUP_ALL);
					if (key < 0) {
						results_lookup_arr.add(~key, val);
					}
				}
			}
			reader.close();
		}

		Number[][] results_lookup = results_lookup_arr.toArray(new Number[results_lookup_arr.size()][]);
		return results_lookup;
	}

	private static void executeOptTrendFittingCallable(ExecutorService exec, OptTrendFittingCallable_FS[] opt_callable,
			Future<HashMap<String, Object>>[] output_future, Number[][] best_result_collection, File resList,
			int inPoolFrom, int inPoolTo) {
		exec.shutdown();
		try {
			if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
				throw new InterruptedException("Time out");
			}
		} catch (InterruptedException e) {
			e.printStackTrace(System.err);
		}
		for (int c = inPoolFrom; c < Math.min(inPoolTo, opt_callable.length); c++) {
			if (opt_callable[c] != null) {
				try {
					HashMap<String, Object> opt_output = output_future[c].get();
					best_result_collection[c][0] = (Double) opt_output
							.get(OptTrendFittingFunction.OPT_TREND_CALLABLE_OUTPUT_BEST_SO_FAR);
				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace(System.err);
				}
			}
		}

		try {
			exportResultCollection(best_result_collection, resList);
		} catch (FileNotFoundException e) {
			e.printStackTrace(System.err);
		}

	}

	private static void exportResultCollection(Number[][] result_collection, File resList)
			throws FileNotFoundException {

		Arrays.sort(result_collection, new Comparator<Number[]>() {
			@Override
			public int compare(Number[] o1, Number[] o2) {
				int res = -Double.compare((Double) o1[o1.length - 1], (Double) o2[o2.length - 1]);
				if (res == 0) {
					res = Long.compare((Long) o1[1], (Long) o2[1]); // cMap_seed
				}
				return 0;
			}
		});

		PrintWriter pWri = new PrintWriter(resList);
		pWri.println("CMAP_SEED,SIM_SEED,PARAM");

		for (Number[] row : result_collection) {
			StringBuilder line = null;
			for (Number n : row) {
				if (line == null) {
					line = new StringBuilder();

				} else {
					line.append(',');
				}
				if (n != null) {
					line.append(n.toString());
				} else {
					line.append(Double.NaN);
				}
			}
			pWri.println(line.toString());
		}
		pWri.close();
	}

	private static void runBayesianOpt(OptFittingFunction opt_trend_obj_func, double[] init_param,
			double[][] boundaries, int numEval, File baseDir, long seed, boolean disp)
			throws FileNotFoundException, IOException {
		runBayesianOpt(opt_trend_obj_func, init_param, boundaries, numEval, baseDir, seed, disp, "");
	}

	private static void runBayesianOpt(OptFittingFunction opt_trend_obj_func, double[] init_param,
			double[][] boundaries, int numEval, File baseDir, long seed, boolean disp, String prefix)
			throws FileNotFoundException, IOException {
		double sigma = 0.06;
		double xi = 0.01;
		double noise = 0.01;
		int maxSize = 10000;

		ArrayList<double[]> observations = new ArrayList<>(); // cMap_Seed, sim_seed, unbounded_parameters, residue
		int nextRowToPrint = 0;
		int numEval_current = 0;
		double[] bestRow = null;

		// Observation CSV
		File obsCSV = new File(baseDir, prefix + FILENAME_OPT_BAYESIAN_OBS);
		if (obsCSV.exists()) {
			BufferedReader reader = new BufferedReader(new FileReader(obsCSV));
			String line;
			while ((line = reader.readLine()) != null) {
				String[] val_str = line.split(",");
				double[] val = new double[val_str.length];
				for (int i = 0; i < val.length; i++) {
					val[i] = Double.parseDouble(val_str[i]);
				}
				observations.add(val);
				nextRowToPrint++;
			}
			reader.close();
		}

		OptFittingFunctionWrapper opt_trend_obj_func_wrapper = new OptFittingFunctionWrapper(opt_trend_obj_func,
				boundaries[0], boundaries[1]);

		if (observations.size() < 2 * init_param.length) {
			addBayesianObservationCollection(observations, opt_trend_obj_func_wrapper, init_param, disp);
			numEval_current++;
			// Add initial observations
			for (int i = 0; i < init_param.length; i++) {
				double[] adj_param = Arrays.copyOf(init_param, init_param.length);
				adj_param[i] = (init_param[i] + boundaries[0][i]) / 2;
				addBayesianObservationCollection(observations, opt_trend_obj_func_wrapper, adj_param, disp);
				numEval_current++;
				adj_param[i] = (init_param[i] + boundaries[1][i]) / 2;
				addBayesianObservationCollection(observations, opt_trend_obj_func_wrapper, adj_param, disp);
				numEval_current++;
			}
			nextRowToPrint = exportBayesianObservationsToFile(observations, obsCSV, opt_trend_obj_func_wrapper, 0,
					prefix);
		}

		int numEval_base = numEval_current;

		while (numEval_current < numEval) {
			// Surrogate function

			double[][] x = new double[observations.size()][]; // unbounded
			double[] y = new double[observations.size()];
			GaussianProcessRegression<double[]> surrogate;

			double best_y = Double.NEGATIVE_INFINITY; // Max_y
			double[] best_x = null; // unbounded
			int pt = 0;

			for (double[] arr : observations) {
				x[pt] = Arrays.copyOfRange(arr, 2, arr.length - 1);
				y[pt] = arr[arr.length - 1];

				if (best_y < y[pt]) {
					best_y = y[pt];
					best_x = x[pt];
					bestRow = arr;
				}

				pt++;
			}

			if (observations.size() < maxSize) {
				surrogate = GaussianProcessRegression.fit(x, y, new GaussianKernel(sigma), noise);
			} else {
				RandomGenerator rngPick = new MersenneTwisterRandomGenerator(seed);
				double[][] t = new double[maxSize][];
				int t_pt = 0;
				for (int i = 0; i < x.length && t_pt < t.length; i++) {
					if (rngPick.nextInt(x.length - i) < maxSize - t_pt) {
						t[t_pt] = Arrays.copyOf(x[i], x[i].length);
						t_pt++;
					}
				}
				surrogate = GaussianProcessRegression.fit(x, y, t, new GaussianKernel(sigma), noise);
			}

			// Expected improvement
			// See Brochu et. al 2010, A Tutorial on Bayesian Optimization of
			// Expensive Cost Functions, with Application to
			// Active User Modeling and Hierarchical
			// Reinforcement Learning (Eq4)
			final double BEST_Y = best_y;
			final GaussianDistribution norm = new GaussianDistribution(0, 1);
			final double ADJ_XI = Math.pow(xi, numEval_current - numEval_base + 1);
			MultivariateFunction expected_improvement = new MultivariateFunction() {
				@Override
				public double value(double[] point) {
					double[] est = new double[2];
					surrogate.predict(point, est);

					double eI = 0;

					if (est[1] != 0) {
						double z = (est[0] - BEST_Y - ADJ_XI) / est[1];
						eI = (est[0] - BEST_Y - ADJ_XI) * norm.cdf(z) + est[1] * norm.p(z);
						// System.out.printf("CDF(z)=%f, PDF(z)=%f",norm.cdf(z), norm.p(z));
						// System.out.printf("Unbounbed param = %s, EI = %f\n", Arrays.toString(point),
						// eI);
					}

					return eI;
				}
			};

			// Use simplex to find the maximum
			double[] next_x; // unbounded
			final double RELATIVE_TOLERANCE = 1e-5;
			final double ABSOLUTE_TOLERANCE = 1e-10;
			final int SIMPLEX_MAX_EVAL = 100;

			ObjectiveFunction objFunc = new ObjectiveFunction(expected_improvement);

			InitialGuess initial_guess;
			final NelderMeadSimplex simplex;

			initial_guess = new InitialGuess(best_x);

			simplex = new NelderMeadSimplex(best_x.length);

			SimplexOptimizer optimizer = new SimplexOptimizer(RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);

			try {
				PointValuePair pV;
				pV = optimizer.optimize(objFunc, simplex, GoalType.MAXIMIZE, initial_guess,
						new MaxEval(SIMPLEX_MAX_EVAL));
				next_x = pV.getPoint();

			} catch (org.apache.commons.math3.exception.TooManyEvaluationsException ex) {
				PointValuePair[] res = simplex.getPoints();
				Arrays.sort(res, new Comparator<PointValuePair>() {
					@Override
					public int compare(PointValuePair o1, PointValuePair o2) {
						return Double.compare(o1.getValue(), o2.getValue());
					}

				});
				next_x = res[0].getPoint();
			}

			addBayesianObservationCollection(observations, opt_trend_obj_func_wrapper,
					opt_trend_obj_func_wrapper.unboundedToBounded(next_x), disp);
			nextRowToPrint = exportBayesianObservationsToFile(observations, obsCSV, opt_trend_obj_func_wrapper,
					nextRowToPrint, prefix);
			numEval_current++;
		}

		if (numEval_current >= numEval) {
			System.out.printf(prefix + " Max number of evaluations (%d) reached.\n", numEval);
		}

		if (bestRow != null) {
			System.out.printf(prefix + " Best obs. = %s\n", Arrays.toString(bestRow));
			double[] best_param_unbound = Arrays.copyOfRange(bestRow, 2, bestRow.length - 1);
			System.out.printf(prefix + " Best param (bounded) = %s\n",
					Arrays.toString(opt_trend_obj_func_wrapper.unboundedToBounded(best_param_unbound)));

		}
	}

	private static int exportBayesianObservationsToFile(ArrayList<double[]> observations, File obsCSV,
			MultivariateFunctionMappingAdapter wrapper, int fromRow, String prefix) throws IOException {

		PrintWriter pWri = new PrintWriter(new FileWriter(obsCSV, true));

		File dispFile = new File(obsCSV.getParent(), prefix + FILENAME_OPT_BAYESIAN_OBS_BOUNDED);

		PrintWriter pWri_Disp = new PrintWriter(new FileWriter(dispFile, true));

		if (!dispFile.exists()) {
			pWri_Disp.println("CMAP_Seed, Sim_Seed,Bounded_Obs");
		}

		int nextRow = fromRow;
		for (int r = fromRow; r < observations.size(); r++) {
			double[] arr = observations.get(r);
			double[] x = Arrays.copyOfRange(arr, 2, arr.length - 1);
			double[] bounded = wrapper.unboundedToBounded(x);

			pWri.print(arr[0]); // CMAP_Seed;
			pWri.print(',');
			pWri.print(arr[1]); // Sim_Seed;

			pWri_Disp.print(arr[0]); // CMAP_Seed;
			pWri_Disp.print(',');
			pWri_Disp.print(arr[1]); // Sim_Seed;

			for (double a : x) {
				pWri.print(',');
				pWri.print(a);
			}
			for (double a : bounded) {
				pWri_Disp.print(',');
				pWri_Disp.print(a);
			}
			pWri.print(',');
			pWri.print(arr[arr.length - 1]);

			pWri_Disp.print(',');
			pWri_Disp.print(arr[arr.length - 1]);

			pWri.println();
			pWri_Disp.println();
			nextRow++;
		}
		pWri.close();
		pWri_Disp.close();

		return nextRow;
	}

	private static void addBayesianObservationCollection(ArrayList<double[]> observationsCollection,
			OptFittingFunctionWrapper func, double[] param, boolean disp_val) {

		double[] unbound_param = func.boundedToUnbounded(param);

		double[] val = new double[param.length + 3];

		double best_r = func.value(unbound_param);
		double[] r_by_runnable = func.getBoundedFunc().getBestResidue_by_runnable();
		long[] cMap_seeds = func.getBoundedFunc().getCMap_seeds();
		long[] sim_seeds = func.getBoundedFunc().getSim_seeds();

		for (int r = 0; r < r_by_runnable.length; r++) {
			Arrays.fill(val, Double.NaN);
			val[0] = cMap_seeds[r];
			val[1] = sim_seeds[r];
			int pt = 2;
			for (double p : unbound_param) {
				val[pt] = p;
				pt++;
			}
			if (Double.isNaN(-r_by_runnable[r])) {
				val[pt] = -best_r;
			} else {
				val[pt] = -r_by_runnable[r]; // Maximum
			}
			observationsCollection.add(val);
			if (disp_val) {
				System.out.printf("Added Obs : %s -> %f\n", Arrays.toString(param), val[pt]);
			}
		}

	}

	public static void stable_prevalence_by_tranmission_fit_GA(String[] args)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException, InterruptedException {
		final String USAGE_INFO = "Usage: PROP_FILE_DIRECTORY INIT_PARAM_VALUE (double[]) BOUNDARIES (double[][]) GA_POP_SIZE (int) "
				+ "<-ta TOURNAMENT_ARITY (int)> <-mr MUTATION_RATE (float)> <-mg MAX_GENERATION (int)> "
				+ "<-rss RNG_SEED_SKIP (int)> " + "<-nib NUM_IN_BATCH (int)> " + "<-exportAll true|false> "
				+ "<-useInitalValues false|true> " + "<-useCMapMapping true|false> " + "<-useCMapEdgeList true|false>"
				+ "<-singleExecutor true|false> " + "<" + Simulation_ClusterModelTransmission.LAUNCH_ARGS_PRINT_PROGRESS
				+ " false|true>";

		if (args.length < 3) {
			System.out.println(USAGE_INFO);
			System.exit(0);
		} else {
			File baseDir = new File(args[0]);

			final double[] init_transmissionProb = (double[]) PropValUtils.propStrToObject(args[1], double[].class);
			final double[][] boundaries = (double[][]) PropValUtils.propStrToObject(args[2], double[][].class);
			final int GA_MAX_POP_SIZE = (int) Integer.parseInt(args[3]);
			final SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

			float max_num_generation = Float.POSITIVE_INFINITY;
			int tournament_arity = GA_MAX_POP_SIZE / 2;
			float mutation_rate = 0.1f;
			boolean exportAll = true;
			boolean useInitValue = false;
			boolean useCMapMapping = true;
			boolean useCMapEdgeList = true;
			boolean useSingleExecutor = true;
			boolean printProgress = false;
			int rngSeedSkip = 0;
			int numInBatch = 1;

			for (int a = 4; a < args.length; a += 2) {
				if ("-ta".equals(args[a])) {
					tournament_arity = Integer.parseInt(args[a + 1]);
				}
				if ("-mr".equals(args[a])) {
					mutation_rate = Float.parseFloat(args[a + 1]);
				}
				if ("-mg".equals(args[a])) {
					max_num_generation = Float.parseFloat(args[a + 1]);
				}
				if ("-rss".equals(args[a])) {
					rngSeedSkip = Integer.parseInt(args[a + 1]);
				}
				if ("-nib".equals(args[a])) {
					numInBatch = Integer.parseInt(args[a + 1]);
				}
				if ("-exportAll".equals(args[a])) {
					exportAll = Boolean.parseBoolean(args[a + 1]);
				}
				if ("-useInitalValues".equals(args[a])) {
					useInitValue = Boolean.parseBoolean(args[a + 1]);
				}
				if ("-useCMapMapping".equals(args[a])) {
					useCMapMapping = Boolean.parseBoolean(args[a + 1]);
				}
				if ("-useCMapEdgeList".equals(args[a])) {
					useCMapEdgeList = Boolean.parseBoolean(args[a + 1]);
				}
				if ("-singleExecutor".equals(args[a])) {
					useSingleExecutor = Boolean.parseBoolean(args[a + 1]);
				}
				if (Simulation_ClusterModelTransmission.LAUNCH_ARGS_PRINT_PROGRESS.equals(args[a])) {
					printProgress = Boolean.parseBoolean(args[a + 1]);
				}

			}

			final File propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);

			final String PROP_SEED_INFECTION = POP_PROP_INIT_PREFIX
					+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
							+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
							+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
							+ Simulation_ClusterModelTransmission.SIM_FIELD_SEED_INFECTION); // "POP_PROP_INIT_PREFIX_14";

			final int GA_ENT_FITNESS = 0;
			final int GA_ENT_CMAP_SEED = GA_ENT_FITNESS + 1;
			final int GA_ENT_SIM_SEED = GA_ENT_CMAP_SEED + 1;
			final int GA_ENT_PARM_START = GA_ENT_SIM_SEED + 1;

			final double ABSOLUTE_TOLERANCE = 1e-10;
			final File GA_ALL_FILE = new File(baseDir, "OptRes_GA_All.csv");
			final File GA_NEXT_POP_FILE = new File(baseDir, "OptRes_Next_GA_POP.csv");

			ArrayList<Number[]> ga_population = new ArrayList<>(GA_MAX_POP_SIZE + 1);

			Comparator<Number[]> ga_population_cmp = new Comparator<Number[]>() {
				@Override
				public int compare(Number[] o1, Number[] o2) {
					int res = Double.compare((Double) o1[GA_ENT_FITNESS], (Double) o2[GA_ENT_FITNESS]);
					for (int i = 0; i < o1.length && res == 0; i++) {
						if (o1[i] instanceof Long) {
							res = Long.compare((Long) o1[i], (Long) o2[i]);
						} else {
							res = Double.compare((Double) o1[i], (Double) o2[i]);
						}
					}
					return res;
				}

			};

			if (propFile.exists()) {
				FileInputStream fIS = new FileInputStream(propFile);
				Properties prop = new Properties();
				prop.loadFromXML(fIS);
				fIS.close();

				System.out.println(String.format("Properties file < %s > loaded.", propFile.getAbsolutePath()));

				long seed = System.currentTimeMillis();
				int numSnap = 1;
				int num_time_steps_per_snap = 1;
				int[] pop_composition = new int[] { 500000, 500000, 20000, 20000 };
				int numThreads = Runtime.getRuntime().availableProcessors();
				int[][] seed_infection = null;
				int contact_map_start_time = 365;
				float[][] opt_target = new float[0][];

				File contactMapDir = baseDir;

				if (prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC) != null) {
					contactMapDir = new File(
							prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC));
					if (!contactMapDir.exists() || !contactMapDir.isDirectory()) {
						contactMapDir = baseDir;
					}
				}

				if (prop.containsKey(POP_PROP_OPT_TARGET)) {
					opt_target = (float[][]) PropValUtils.propStrToObject(prop.getProperty(POP_PROP_OPT_TARGET),
							float[][].class);
				}

				if (prop.containsKey(PROP_SEED_INFECTION)) {
					seed_infection = (int[][]) PropValUtils.propStrToObject(prop.getProperty(PROP_SEED_INFECTION),
							int[][].class);
				}
				if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
					seed = Long.parseLong(
							prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
				}
				if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])) {
					numSnap = Integer.parseInt(
							prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
				}
				if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
					num_time_steps_per_snap = Integer.parseInt(
							prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ]));
				}
				if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL])) {
					numThreads = Integer.parseInt(
							prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL]));
				}
				String popCompositionKey = POP_PROP_INIT_PREFIX
						+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
				if (prop.containsKey(popCompositionKey)) {
					pop_composition = (int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey),
							int[].class);
				}
				String contactMapRangeKey = POP_PROP_INIT_PREFIX
						+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
								+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
								+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);
				if (prop.containsKey(contactMapRangeKey)) {
					contact_map_start_time = ((int[]) PropValUtils.propStrToObject(prop.getProperty(contactMapRangeKey),
							int[].class))[0];
				}

				RandomGenerator RNG = new MersenneTwisterRandomGenerator(seed);

				for (int s = 0; s < rngSeedSkip; s++) {
					RNG = new MersenneTwisterRandomGenerator(RNG.nextLong());
				}

				System.out.printf("# Sim Seed = %d. # processors available = %d. # threads used = %d.\n", seed,
						Runtime.getRuntime().availableProcessors(), numThreads);

				// Check for contact cluster generated
				final Pattern pattern_baseCMap_filename = Pattern.compile(CMAP_REGEX_STR);

				File[] preGenClusterFiles = contactMapDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isFile() && Pattern.matches(CMAP_REGEX_STR, pathname.getName());

					}
				});

				ConcurrentHashMap<Long, ContactMap> cMap_mapping = null;
				ConcurrentHashMap<Long, ArrayList<Integer[]>> cMap_edgelist_mapping = null;
				if (useCMapMapping) {
					cMap_mapping = new ConcurrentHashMap<>(preGenClusterFiles.length);
				}
				if (useCMapEdgeList) {
					cMap_edgelist_mapping = new ConcurrentHashMap<>(preGenClusterFiles.length);
				}

				long[] BASE_CONTACT_MAP_SEED = new long[preGenClusterFiles.length];

				for (int mapPt = 0; mapPt < preGenClusterFiles.length; mapPt++) {
					File cMap_file = preGenClusterFiles[mapPt];
					Matcher m = pattern_baseCMap_filename.matcher(cMap_file.getName());
					m.matches();
					BASE_CONTACT_MAP_SEED[mapPt] = Long.parseLong(m.group(1));

				}

				if (GA_NEXT_POP_FILE.exists()) {
					// Fill previous pop file
					BufferedReader reader = new BufferedReader(new FileReader(GA_NEXT_POP_FILE));
					String line;
					while ((line = reader.readLine()) != null) {
						String[] ent = line.split(",");
						Number[] entArr = new Number[ent.length];
						for (int i = 0; i < entArr.length; i++) {
							if (i == GA_ENT_CMAP_SEED || i == GA_ENT_SIM_SEED) {
								entArr[i] = Long.parseLong(ent[i]);
							} else {
								entArr[i] = Double.parseDouble(ent[i]);
							}
						}
						ga_population.add(entArr);
					}
					reader.close();

					// Check if there is replacement entries from GA_ALL_FILE
					if (GA_ALL_FILE.exists()) {
						Comparator<Number[]> ga_population_equ_cmp = new Comparator<Number[]>() {
							@Override
							public int compare(Number[] o1, Number[] o2) {
								int res = 0;
								for (int i = GA_ENT_CMAP_SEED; i < o1.length && res == 0; i++) {
									if (o1[i] instanceof Long) {
										res = Long.compare((Long) o1[i], (Long) o2[i]);
									} else {
										res = Double.compare((Double) o1[i], (Double) o2[i]);
									}
								}
								return res;
							}

						};

						ArrayList<Number[]> ga_all_array = new ArrayList<>();
						reader = new BufferedReader(new FileReader(GA_ALL_FILE));
						while ((line = reader.readLine()) != null) {
							String[] ent = line.split(",");
							Number[] entArr = new Number[ent.length];
							for (int i = 0; i < entArr.length; i++) {
								if (i == GA_ENT_CMAP_SEED || i == GA_ENT_SIM_SEED) {
									entArr[i] = Long.parseLong(ent[i]);
									// Dummy parameters to skip RNG value
									if (i == GA_ENT_CMAP_SEED) {
										RNG.nextInt(BASE_CONTACT_MAP_SEED.length);
									} else {
										RNG.nextLong();
									}
								} else {
									entArr[i] = Double.parseDouble(ent[i]);
									if (i != GA_ENT_FITNESS) {
										// Dummy parameters to skip RNG value
										RNG.nextDouble();
									}
								}
							}
							int key = Collections.binarySearch(ga_all_array, entArr, ga_population_equ_cmp);
							if (key < 0) {
								ga_all_array.add(~key, entArr);
							}

						}
						reader.close();

						// Check for updated GA population, and generate new one if needed.
						boolean updateGA = false;
						for (Number[] ent : ga_population) {
							if (Double.isNaN(ent[GA_ENT_FITNESS].doubleValue())) {
								int key = Collections.binarySearch(ga_all_array, ent, ga_population_equ_cmp);
								if (key >= 0) {
									updateGA = true;
									ent[GA_ENT_FITNESS] = ga_all_array.get(key)[GA_ENT_FITNESS];
								}
							}
						}
						if (updateGA && (rngSeedSkip == 0 || numInBatch == 1)) {
							PrintWriter pWri = new PrintWriter(GA_NEXT_POP_FILE);
							for (Number[] ga_ent : ga_population) {
								StringBuilder linebuilder = new StringBuilder();
								for (Number val : ga_ent) {
									if (linebuilder.length() != 0) {
										linebuilder.append(',');
									}
									linebuilder.append(val.toString());
								}
								pWri.println(linebuilder.toString());
							}
							pWri.close();
						}

					}

				} else if (GA_ALL_FILE.exists()) {
					BufferedReader reader = new BufferedReader(new FileReader(GA_ALL_FILE));
					String line;
					while ((line = reader.readLine()) != null) {
						String[] ent = line.split(",");
						Number[] entArr = new Number[ent.length];
						for (int i = 0; i < entArr.length; i++) {
							if (i == GA_ENT_CMAP_SEED || i == GA_ENT_SIM_SEED) {
								entArr[i] = Long.parseLong(ent[i]);
								// Dummy parameters to skip RNG value
								if (i == GA_ENT_CMAP_SEED) {
									RNG.nextInt(BASE_CONTACT_MAP_SEED.length);
								} else {
									RNG.nextLong();
								}
							} else {
								entArr[i] = Double.parseDouble(ent[i]);
								if (i != GA_ENT_FITNESS) {
									// Dummy parameters to skip RNG value
									RNG.nextDouble();
								}
							}
						}
						int key = Collections.binarySearch(ga_population, entArr, ga_population_cmp);
						if (key < 0) {
							ga_population.add(~key, entArr);
						}
						if (ga_population.size() > GA_MAX_POP_SIZE) {
							ga_population.remove(GA_MAX_POP_SIZE);
						}
					}
					reader.close();
				}

				// Populate GA_Population
				int prefill_size = ga_population.size();

				// Pre-fill with initial values

				if (useInitValue || prefill_size == 0) {
					Number[] init_ent = new Number[init_transmissionProb.length + 3];
					init_ent[GA_ENT_FITNESS] = Double.NaN;
					init_ent[GA_ENT_CMAP_SEED] = BASE_CONTACT_MAP_SEED[RNG.nextInt(BASE_CONTACT_MAP_SEED.length)];
					init_ent[GA_ENT_SIM_SEED] = RNG.nextLong();
					for (int p = 0; p < init_transmissionProb.length; p++) {
						init_ent[GA_ENT_PARM_START + p] = init_transmissionProb[p];
					}
					if (prefill_size < GA_MAX_POP_SIZE) {
						ga_population.add(prefill_size, init_ent);
					} else {
						ga_population.set(GA_MAX_POP_SIZE - 1, init_ent);
					}
				}

				// Populate the rest of GA_Population if needed
				prefill_size = ga_population.size();
				for (int g = prefill_size; g < GA_MAX_POP_SIZE; g++) {
					Number[] added_ent = new Number[init_transmissionProb.length + 3];
					added_ent[GA_ENT_FITNESS] = Double.NaN;
					added_ent[GA_ENT_CMAP_SEED] = BASE_CONTACT_MAP_SEED[RNG.nextInt(BASE_CONTACT_MAP_SEED.length)];
					added_ent[GA_ENT_SIM_SEED] = RNG.nextLong();
					for (int p = 0; p < init_transmissionProb.length; p++) {
						added_ent[GA_ENT_PARM_START + p] = boundaries[0][p]
								+ RNG.nextFloat() * (boundaries[1][p] - boundaries[0][p]);
					}
					ga_population.add(g, added_ent);
				}

				boolean exiting = false;
				int num_gen = 0;

				while (!exiting) {
					// GA evolve
					if (num_gen > max_num_generation) {
						System.out.printf("Maximum number of generations (%d) reached.\n", (int) max_num_generation);
						exiting = true;
					} else {
						// Populate GA_POPULATION fitness

						ExecutorService exec = Executors.newFixedThreadPool(numThreads);
						int toBeGeneratedCount = 0;
						int threadCount = 0;

						for (Number[] ga_ent : ga_population) {
							if (Double.isNaN(ga_ent[GA_ENT_FITNESS].doubleValue())) {
								int offset = toBeGeneratedCount - rngSeedSkip;
								if (offset >= 0 && offset % numInBatch == 0) {

									final int[] POP_COMPOSITION = pop_composition;
									final int NUM_TIME_STEPS_PER_SNAP = num_time_steps_per_snap;
									final int NUM_SNAP = numSnap;
									final int runnnable_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
											+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
											+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
											+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD;
									final Properties PROP = prop;
									final int[][] SEED_INFECTION = seed_infection;
									final int START_TIME = contact_map_start_time;
									final float[][] OPT_TARGET = opt_target;
									final boolean EXPORT_ALL = exportAll;
									final File CMAP_DIR = contactMapDir;
									final ConcurrentHashMap<Long, ContactMap> C_MAP_MAPPING = cMap_mapping;
									final ConcurrentHashMap<Long, ArrayList<Integer[]>> C_MAP_EDGES_LIST_MAPPING = cMap_edgelist_mapping;
									final boolean PRINT_PROGRESS = printProgress;

									Runnable fitness_thread = new Runnable() {
										@Override
										public void run() {
											long cMap_seed = (long) ga_ent[GA_ENT_CMAP_SEED];
											ContactMap cmap = null;

											if (C_MAP_MAPPING != null) {
												cmap = C_MAP_MAPPING.get(cMap_seed);
											}

											if (cmap == null) {
												String cmap_match_str = FILENAME_FORMAT_ALL_CMAP
														.replaceFirst("%d", Long.toString(cMap_seed))
														.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");
												File[] cmap_match_file = CMAP_DIR.listFiles(new FileFilter() {
													@Override
													public boolean accept(File pathname) {
														return pathname.isFile()
																&& Pattern.matches(cmap_match_str, pathname.getName());
													}
												});

												if (cmap_match_file.length > 0) {
													Callable<ContactMap> cm_read_callable = Abstract_Runnable_ClusterModel
															.generateContactMapCallable(cmap_match_file[0]);
													try {
														cmap = cm_read_callable.call();
													} catch (Exception e) {
														e.printStackTrace(System.err);
													}
												}

												if (C_MAP_MAPPING != null) {
													C_MAP_MAPPING.put(cMap_seed, cmap);
												}
											}

											if (cmap != null) {

												// Setting up runnable
												final String popCompositionKey = POP_PROP_INIT_PREFIX
														+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
												int[] pop_composition = (int[]) PropValUtils.propStrToObject(
														prop.getProperty(popCompositionKey), int[].class);
												int[] cumulative_pop_composition = new int[pop_composition.length];
												int pop_offset = 0;
												for (int g = 0; g < cumulative_pop_composition.length; g++) {
													cumulative_pop_composition[g] = pop_offset + pop_composition[g];
													pop_offset += pop_composition[g];
												}

												final String riskCatListAllKey = POP_PROP_INIT_PREFIX + Integer
														.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
																+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
																+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
																+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD
																+ Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS);
												float[][] riskCatListAll = (float[][]) PropValUtils.propStrToObject(
														prop.getProperty(riskCatListAllKey), float[][].class);

												File pre_allocate_risk_file = new File(baseDir, String.format(
														Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP,
														(long) ga_ent[GA_ENT_CMAP_SEED]));
												ArrayList<Number[]> riskGrpArr = new ArrayList<>();

												long seed = Long.parseLong(prop.getProperty(
														SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));

												boolean reallocate = true;

												if (!pre_allocate_risk_file.exists()) {

													final String time_rangeKey = POP_PROP_INIT_PREFIX + Integer
															.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
																	+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
																	+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);

													int[] map_time_range = (int[]) PropValUtils.propStrToObject(
															prop.getProperty(time_rangeKey), int[].class);

													Simulation_ClusterModelTransmission
															.fillRiskGrpArrByCasualPartnership(riskGrpArr, cmap,
																	cumulative_pop_composition, riskCatListAll,
																	map_time_range);

												} else {
													try {
														reallocate = Simulation_ClusterModelTransmission
																.loadPreallocateRiskGrp(riskGrpArr, baseDir,
																		(long) ga_ent[GA_ENT_CMAP_SEED]);
													} catch (Exception e) {
														e.printStackTrace(System.err);
													}
												}

												if (reallocate) {
													Simulation_ClusterModelTransmission.reallocateRiskGrp(riskGrpArr,
															(long) ga_ent[GA_ENT_CMAP_SEED], cumulative_pop_composition,
															riskCatListAll, baseDir, seed);
												}

												Runnable_ClusterModel_Transmission runnable = new Runnable_ClusterModel_Transmission(
														(long) ga_ent[GA_ENT_CMAP_SEED], (long) ga_ent[GA_ENT_SIM_SEED],
														POP_COMPOSITION, cmap, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP);
												runnable.setBaseDir(baseDir);

												if (PRINT_PROGRESS) {
													runnable.setPrint_progress(System.out);
													runnable.setRunnableId(
															String.format("%d,%d", (long) ga_ent[GA_ENT_CMAP_SEED],
																	(long) ga_ent[GA_ENT_SIM_SEED]));
												}

												if (C_MAP_EDGES_LIST_MAPPING != null) {
													ArrayList<Integer[]> cMap_edges = C_MAP_EDGES_LIST_MAPPING
															.get(cMap_seed);
													if (cMap_edges == null) {
														try {
															Callable<ArrayList<Integer[]>> callable = Runnable_ClusterModel_Transmission
																	.generateMapEdgeArray(cmap);
															cMap_edges = callable.call();
														} catch (Exception e) {
															e.printStackTrace(System.err);
														}
														C_MAP_EDGES_LIST_MAPPING.put(cMap_seed, cMap_edges);
													}
													runnable.setEdges_list(cMap_edges);
												}

												for (int i = runnnable_offset; i < runnnable_offset
														+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; i++) {

													String key = POP_PROP_INIT_PREFIX + Integer.toString(i);
													if (PROP.containsKey(key)) {
														runnable.getRunnable_fields()[i
																- runnnable_offset] = PropValUtils.propStrToObject(
																		PROP.getProperty(key),
																		runnable.getRunnable_fields()[i
																				- runnnable_offset].getClass());
													}
												}

												runnable.setSimSetting(1); // No output
												Number[] param_number = Arrays.copyOfRange(ga_ent, GA_ENT_PARM_START,
														ga_ent.length);
												double[] param_double = new double[param_number.length];
												for (int i = 0; i < param_number.length; i++) {
													param_double[i] = param_number[i].doubleValue();
												}

												setOptParamInRunnable(runnable, PROP, param_double, false);
												runnable.initialse();
												runnable.fillRiskCatMap(riskGrpArr);
												runnable.allocateSeedInfection(SEED_INFECTION, START_TIME);

												// Run simulation
												runnable.run();

												// Extract results
												int start_k = 2;

												Integer[] keys = new Integer[NUM_SNAP];
												keys[0] = NUM_TIME_STEPS_PER_SNAP;
												for (int k = 1; k < keys.length; k++) {
													keys[k] = keys[k - 1] + NUM_TIME_STEPS_PER_SNAP;
												}

												@SuppressWarnings("unchecked")
												HashMap<Integer, int[][]> infectious_count_map = (HashMap<Integer, int[][]>) runnable
														.getSim_output()
														.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT);

												@SuppressWarnings("unchecked")
												HashMap<Integer, int[]> cumul_treatment_map = (HashMap<Integer, int[]>) runnable
														.getSim_output()
														.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON);

												double sqSum = calculateOptFitness(param_double, OPT_TARGET,
														POP_COMPOSITION, infectious_count_map, cumul_treatment_map,
														keys, start_k,
														String.format("CM_Seed = %d, sim_seed = %d",
																(long) ga_ent[GA_ENT_CMAP_SEED],
																(long) ga_ent[GA_ENT_SIM_SEED]),
														null);

												ga_ent[GA_ENT_FITNESS] = sqSum;

												if (EXPORT_ALL) {
													StringBuilder ga_ent_disp = new StringBuilder();
													for (Number val : ga_ent) {
														if (ga_ent_disp.length() != 0) {
															ga_ent_disp.append(',');
														}
														ga_ent_disp.append(val);
													}
													synchronized (GA_ALL_FILE) {
														try {
															final int max_retryAttempt = 5;
															int attempt = 0;
															while (!GA_ALL_FILE.canWrite()
																	&& attempt < max_retryAttempt) {
																try {
																	Thread.sleep(new Random().nextInt(1000));
																} catch (InterruptedException e) {
																	e.printStackTrace(System.err);
																}
																attempt++;
															}
															FileWriter export_all_fWri = new FileWriter(GA_ALL_FILE,
																	true);
															PrintWriter export_all_pWri = new PrintWriter(
																	export_all_fWri);
															export_all_pWri.println(ga_ent_disp.toString());
															export_all_pWri.close();
															export_all_fWri.close();
														} catch (IOException e) {
															e.printStackTrace(System.err);
														}

													}
												}
											}

										}
									};

									exec.submit(fitness_thread);
									threadCount++;

									if (!useSingleExecutor) {
										if (threadCount == numThreads) {
											exec.shutdown();
											if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
												System.err.println("Thread time-out!");
											}
											exec = Executors.newFixedThreadPool(numThreads);
											threadCount = 0;
										}
									}

								}
								toBeGeneratedCount++;
							}
						}

						if (threadCount != 0) {
							exec.shutdown();
							if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
								System.err.println("Thread time-out!");
							}
						}

						// Sort GA_POPULATION
						Collections.sort(ga_population, ga_population_cmp);

						if (ga_population.get(GA_MAX_POP_SIZE - 1)[GA_ENT_FITNESS].doubleValue()
								- ga_population.get(0)[GA_ENT_FITNESS].doubleValue() < ABSOLUTE_TOLERANCE) {
							System.out.printf(
									"The difference between best and worst fit in GA_POPULATION is less than ABSOLUTE_TOLERANCE of %.f.\n",
									ABSOLUTE_TOLERANCE);
							exiting = true;

						} else {

							// Fill the rest with next generation
							int next_gen_start = tournament_arity;

							for (int g = next_gen_start; g < GA_MAX_POP_SIZE; g++) {
								Number[] ga_ent = ga_population.get(g);
								ga_ent[GA_ENT_FITNESS] = Double.NaN;
								ga_ent[GA_ENT_CMAP_SEED] = BASE_CONTACT_MAP_SEED[RNG
										.nextInt(BASE_CONTACT_MAP_SEED.length)];
								ga_ent[GA_ENT_SIM_SEED] = RNG.nextLong();

								int[] sel_indices = ArrayUtilsRandomGenerator.randomSelectIndex(2, 0, tournament_arity,
										RNG);

								for (int p = GA_ENT_PARM_START; p < ga_ent.length; p++) {
									if (RNG.nextFloat() < mutation_rate) {
										ga_ent[p] = boundaries[0][p - GA_ENT_PARM_START]
												+ RNG.nextFloat() * (boundaries[1][p - GA_ENT_PARM_START]
														- boundaries[0][p - GA_ENT_PARM_START]);
									} else {
										// Best 2 as crossover candidate
										ga_ent[p] = ga_population.get(sel_indices[RNG.nextInt(2)])[p];
									}
								}
							}
						}

						// Export GA_POPULATION
						if (rngSeedSkip == 0 || numInBatch == 1) {
							PrintWriter pWri = new PrintWriter(GA_NEXT_POP_FILE);
							for (Number[] ga_ent : ga_population) {
								StringBuilder line = new StringBuilder();
								for (Number val : ga_ent) {
									if (line.length() != 0) {
										line.append(',');
									}
									line.append(val.toString());
								}
								pWri.println(line.toString());
							}
							pWri.close();

							System.out.printf("Generation #%d formed at %s and exported to %s.\n", num_gen,
									dateFormat.format(new Date(System.currentTimeMillis())),
									GA_NEXT_POP_FILE.getAbsolutePath());
						}

						num_gen++;
					}
				}

			}

		}

	}

	public static void stable_prevalence_by_tranmission_fit_Simplex(String[] args)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException {

		final String FLAG_nEval = "-nEval=";
		final String FLAG_resList = "-resList=";
		final String USAGE_INFO = String
				.format("Usage: PROP_FILE_DIRECTORY INIT_PARAM_VALUE (double[]) BOUNDARIES (double[][])  "
						+ "<optional: %sRESULT_LIST_FILENAME> <optional: %sNUM_EVAL (int)>", FLAG_resList, FLAG_nEval);

		int numEval = 100;
		int minArgLength = 3;

		if (args.length < minArgLength) {
			System.out.println(USAGE_INFO);
			System.exit(0);
		} else {
			File baseDir = new File(args[0]);
			double[] init_transmissionProb = (double[]) PropValUtils.propStrToObject(args[1], double[].class);
			double[][] boundaries = (double[][]) PropValUtils.propStrToObject(args[2], double[][].class);
			File resultList = null;

			for (int a = minArgLength; a < args.length; a++) {
				if (args[a].startsWith(FLAG_nEval)) {
					numEval = (int) Integer.parseInt(args[a].substring(FLAG_nEval.length()));
				}
				if (args[a].startsWith(FLAG_resList)) {
					resultList = new File(baseDir, (args[a].substring(FLAG_resList.length())));
				}
			}
			stable_prevalence_by_tranmission_fit_Simplex(baseDir, init_transmissionProb, boundaries, numEval,
					resultList);
		}

	}

	public static void stable_prevalence_by_tranmission_fit_Simplex(File baseDir, final double[] init_param,
			final double[][] boundaries, int numEval, File resultList)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException {

		final String oPT_RESULT_FILENAME_FORMAT = "Opt_result_%d_%d.txt";

		final File propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);

		if (propFile.exists()) {
			FileInputStream fIS = new FileInputStream(propFile);
			Properties prop = new Properties();
			prop.loadFromXML(fIS);
			fIS.close();

			System.out.println(String.format("Properties file < %s > loaded.", propFile.getAbsolutePath()));

			long seed = System.currentTimeMillis();
			int numSnap = 1;
			int num_time_steps_per_snap = 1;
			int[] pop_composition = new int[] { 500000, 500000, 20000, 20000 };
			int numThreads = Runtime.getRuntime().availableProcessors();
			int[][] seed_infection = null;
			int contact_map_start_time = 365;
			int numSimPerMap = 1;

			float[][] opt_target = new float[0][];

			File contactMapDir = baseDir;

			if (prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC) != null) {
				contactMapDir = new File(prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC));
				if (!contactMapDir.exists() || !contactMapDir.isDirectory()) {
					contactMapDir = baseDir;
				}
			}

			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET])) {
				numSimPerMap = Integer.parseInt(
						prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET]));
			}

			if (prop.containsKey(POP_PROP_OPT_TARGET)) {
				opt_target = (float[][]) PropValUtils.propStrToObject(prop.getProperty(POP_PROP_OPT_TARGET),
						float[][].class);
			}

			if (prop.containsKey(PROP_SEED_INFECTION)) {
				seed_infection = (int[][]) PropValUtils.propStrToObject(prop.getProperty(PROP_SEED_INFECTION),
						int[][].class);
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
				seed = Long
						.parseLong(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])) {
				numSnap = Integer
						.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
				num_time_steps_per_snap = Integer
						.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ]));
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL])) {
				numThreads = Integer.parseInt(
						prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL]));
			}
			String popCompositionKey = POP_PROP_INIT_PREFIX
					+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
			if (prop.containsKey(popCompositionKey)) {
				pop_composition = (int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey),
						int[].class);
			}
			String contactMapRangeKey = POP_PROP_INIT_PREFIX
					+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
							+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
							+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);
			if (prop.containsKey(contactMapRangeKey)) {
				contact_map_start_time = ((int[]) PropValUtils.propStrToObject(prop.getProperty(contactMapRangeKey),
						int[].class))[0];
			}

			if (opt_target == null || opt_target.length == 0) {
				System.out.println("OPT_TARGET missing. Printing out runnable fields instead.");

				Abstract_Runnable_ClusterModel_Transmission runnable = new Runnable_ClusterModel_Transmission(-1, -1,
						pop_composition, null, num_time_steps_per_snap, numSnap);
				runnable.setBaseDir(baseDir);

				// Set fields based on prop file
				for (int i = RUNNABLE_OFFSET; i < RUNNABLE_OFFSET
						+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; i++) {
					String key = POP_PROP_INIT_PREFIX + Integer.toString(i);
					if (prop.containsKey(key)) {
						runnable.getRunnable_fields()[i - RUNNABLE_OFFSET] = PropValUtils.propStrToObject(
								prop.getProperty(key), runnable.getRunnable_fields()[i - RUNNABLE_OFFSET].getClass());
					}
				}

				runnable.setSimSetting(1); // No output
				setOptParamInRunnable(runnable, prop, init_param, true);
			}

			// Check for contact cluster generated

			File[] pregen_cMap_files = contactMapDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isFile() && Pattern.matches(CMAP_REGEX_STR, pathname.getName());

				}
			});

			// Parameters to be pass into anonymous functions
			final ContactMap[] bASE_CONTACT_MAP = new ContactMap[pregen_cMap_files.length];
			final long[] bASE_CONTACT_MAP_SEED = new long[bASE_CONTACT_MAP.length];
			final int nUM_THREADS = numThreads;
			final int[] pOP_COMPOSITION = pop_composition;
			final int nUM_TIME_STEPS_PER_SNAP = num_time_steps_per_snap;
			final int nUM_SNAP = numSnap;
			final int sTART_TIME = contact_map_start_time;
			final int[][] sEED_INFECTION = seed_infection;
			final float[][] oPT_TARGET = opt_target;

			// Loading of cMap;
			long tic = System.currentTimeMillis();
			int cMap_count = extractContactMap(bASE_CONTACT_MAP, bASE_CONTACT_MAP_SEED, pregen_cMap_files, nUM_THREADS);
			System.out.printf("%d ContactMap(s) from %s loaded. Time req. = %.3fs\n", cMap_count,
					contactMapDir.getAbsolutePath(), (System.currentTimeMillis() - tic) / 1000f);

			if (resultList != null) {
				ArrayList<Number[]> resultArr = new ArrayList<>();
				final String rESULT_KEY_FORMAT = "%d_%d"; // Key: "[CMAP_SEED]_[SIM_SEED]"
				HashMap<String, Number[]> bestResultMap = new HashMap<>();

				// Generate or reading of result list
				if (!resultList.exists()) {
					RandomGenerator rng = new MersenneTwisterRandomGenerator(seed);
					long[] sim_seeds = new long[numSimPerMap];
					for (int i = 0; i < sim_seeds.length; i++) {
						sim_seeds[i] = rng.nextLong();
					}
					for (long cMap_seed : bASE_CONTACT_MAP_SEED) {
						for (int i = 0; i < sim_seeds.length; i++) {
							Number[] resultArr_ent = new Number[2 + init_param.length + 1]; // cMap_seed, sim_seed,
																							// param, residue
							Arrays.fill(resultArr_ent, Double.NaN);
							resultArr_ent[0] = cMap_seed;
							resultArr_ent[1] = sim_seeds[i];
							resultArr_ent[resultArr_ent.length - 1] = Double.POSITIVE_INFINITY;
							resultArr.add(resultArr_ent);
							bestResultMap.put(String.format(rESULT_KEY_FORMAT, cMap_seed, sim_seeds[i]), resultArr_ent);
						}
					}

				} else {
					BufferedReader reader = new BufferedReader(new FileReader(resultList));
					reader.readLine(); // Header line

					String line;
					while ((line = reader.readLine()) != null) {
						String[] ent = line.split(",");
						long cMap_seed = Long.parseLong(ent[0]);
						Number[] resultArr_ent = new Number[ent.length];
						Arrays.fill(resultArr_ent, Double.NaN);
						resultArr_ent[0] = cMap_seed;
						resultArr_ent[1] = Long.parseLong(ent[1]);
						for (int i = 2; i < resultArr_ent.length; i++) {
							resultArr_ent[i] = Double.parseDouble(ent[i]);
						}

						resultArr.add(resultArr_ent);
						bestResultMap.put(String.format(rESULT_KEY_FORMAT, cMap_seed, resultArr_ent[1].longValue()),
								resultArr_ent);
					}
					reader.close();
				}

				// Reading past result and extract the best one
				Pattern pattern_previous_opt_out_file = Pattern
						.compile(oPT_RESULT_FILENAME_FORMAT.replaceAll("%d", "(-{0,1}\\\\d+)"));
				File[] prev_output_files = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						Matcher m = pattern_previous_opt_out_file.matcher(pathname.getName());
						return m.matches();
					}
				});

				HashMap<String, ArrayList<Number[]>> all_result_map = new HashMap<>();

				for (File prev_output_file : prev_output_files) {
					BufferedReader reader = new BufferedReader(new FileReader(prev_output_file));
					String line;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(OPT_OUTPUT_PREFIX_CMAP)) {
							try {
								long cMap_seed = Long.parseLong(line.substring(OPT_OUTPUT_PREFIX_CMAP.length()));
								line = reader.readLine();
								long sim_seed = Long.parseLong(line.substring(OPT_OUTPUT_PREFIX_SIMSEED.length()));
								line = reader.readLine();
								String[] paramString = (line.substring(
										Optimisation_Factory.OPT_OUTPUT_PREFIX_PARAM.length() + 1, line.length() - 1))
										.split(",");
								line = reader.readLine();
								double residue = Double.parseDouble(line.substring(OPT_OUTPUT_PREFIX_RESIDUE.length()));

								String key = String.format(rESULT_KEY_FORMAT, cMap_seed, sim_seed);

								ArrayList<Number[]> res_store = all_result_map.get(key);
								if (res_store == null) {
									res_store = new ArrayList<>();
									all_result_map.put(key, res_store);
								}
								Number[] readEnt = new Number[2 + paramString.length + 1];
								readEnt[0] = cMap_seed;
								readEnt[1] = sim_seed;
								readEnt[readEnt.length - 1] = residue;
								for (int i = 0; i < paramString.length; i++) {
									readEnt[i + 2] = Double.parseDouble(paramString[i]);
								}
								res_store.add(readEnt);

								Number[] ent = bestResultMap.get(key);
								if (ent != null) {
									double res_current = ent[ent.length - 1].doubleValue();
									if (Double.isNaN(res_current) || residue < res_current) {
										for (int i = 2; i < ent.length; i++) {
											ent[i] = readEnt[i];
										}
									}
								}
							} catch (IOException | NullPointerException ex) {
								ex.printStackTrace(System.err);
							}

						}
					}
					reader.close();
				}

				Number[][] result_collection_best_so_far = resultArr.toArray(new Number[resultArr.size()][]);
				exportResultCollection(result_collection_best_so_far, resultList);

				HashMap<Long, ContactMap> cMap_mapping = new HashMap<>();
				for (int c = 0; c < bASE_CONTACT_MAP_SEED.length; c++) {
					cMap_mapping.put(bASE_CONTACT_MAP_SEED[c], bASE_CONTACT_MAP[c]);
				}

				ExecutorService exec = null;
				int inExec = 0;

				for (Number[] result_single_param : result_collection_best_so_far) {
					Long cMap_seed = (Long) result_single_param[0];
					Long sim_seed = (Long) result_single_param[1];
					ArrayList<Number[]> all_result = all_result_map
							.get(String.format(rESULT_KEY_FORMAT, cMap_seed, sim_seed));
					Number[][] all_result_arr;

					if (all_result != null) {
						all_result_arr = all_result.toArray(new Number[all_result.size()][]);
						Arrays.sort(all_result_arr, COMPARATOR_RESULT_LOOKUP_ALL);
					} else {
						all_result_arr = new Number[0][];
					}

					MultivariateFunction func = new MultivariateFunction() {
						@Override
						public double value(double[] point) {
							// Look up previous results
							Number[] val = new Number[point.length + 3];
							val[0] = cMap_seed;
							val[1] = sim_seed;
							for (int p = 2; p < val.length - 1; p++) {
								val[p] = point[p - 2];
							}
							int all_res_k = Arrays.binarySearch(all_result_arr, val,
									Optimisation_Factory.COMPARATOR_RESULT_LOOKUP_PARAM);

							if (all_res_k >= 0) {
								Number[] found_res = all_result_arr[all_res_k];
								return found_res[found_res.length - 1].doubleValue();
							} else {

								Runnable_ClusterModel_Transmission runnable = new Runnable_ClusterModel_Transmission(
										cMap_seed, sim_seed, pOP_COMPOSITION, cMap_mapping.get(cMap_seed),
										nUM_TIME_STEPS_PER_SNAP, nUM_SNAP);
								runnable.setBaseDir(baseDir);

								// Set fields based on prop file
								for (int i = RUNNABLE_OFFSET; i < RUNNABLE_OFFSET
										+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; i++) {
									String key = POP_PROP_INIT_PREFIX + Integer.toString(i);
									if (prop.containsKey(key)) {
										runnable.getRunnable_fields()[i - RUNNABLE_OFFSET] = PropValUtils
												.propStrToObject(prop.getProperty(key),
														runnable.getRunnable_fields()[i - RUNNABLE_OFFSET].getClass());
									}
								}

								runnable.setSimSetting(1); // No output
								setOptParamInRunnable(runnable, prop, point, false);
								runnable.initialse();

								// Set pre-allocated risk group
								File pre_allocate_risk_file = new File(baseDir, String.format(
										Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, cMap_seed));
								if (pre_allocate_risk_file.exists()) {
									ArrayList<Number[]> riskGrpArr = new ArrayList<>();
									try {
										Simulation_ClusterModelTransmission.loadPreallocateRiskGrp(riskGrpArr, baseDir,
												cMap_seed);
									} catch (Exception e) {
										e.printStackTrace(System.err);
									}
									runnable.fillRiskCatMap(riskGrpArr);
								}
								runnable.allocateSeedInfection(sEED_INFECTION, sTART_TIME);

								runnable.run();
								int start_k = 2;

								Integer[] keys = new Integer[nUM_SNAP];
								keys[0] = nUM_TIME_STEPS_PER_SNAP;
								for (int k = 1; k < keys.length; k++) {
									keys[k] = keys[k - 1] + nUM_TIME_STEPS_PER_SNAP;
								}

								@SuppressWarnings("unchecked")
								HashMap<Integer, int[][]> infectious_count_map = (HashMap<Integer, int[][]>) runnable
										.getSim_output()
										.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT);

								@SuppressWarnings("unchecked")
								HashMap<Integer, int[]> cumul_treatment_map = (HashMap<Integer, int[]>) runnable
										.getSim_output()
										.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON);

								StringBuilder str_disp = new StringBuilder();

								double sqSum = calculateOptFitness(point, oPT_TARGET, pOP_COMPOSITION,
										infectious_count_map, cumul_treatment_map, keys, start_k,
										String.format("CM_Seed = %d, sim_seed = %d", cMap_seed, sim_seed), str_disp);

								// Print simulation result
								try {
									File opt_output_file = new File(baseDir,
											String.format(oPT_RESULT_FILENAME_FORMAT, cMap_seed, sim_seed));
									boolean newFile = !opt_output_file.exists();
									PrintWriter pWri = new PrintWriter(new FileWriter(opt_output_file, true));

									if (newFile) {
										pWri.printf("Opt target (%s in total):\n", oPT_TARGET.length);
										for (float[] opt_target_ent : oPT_TARGET) {
											pWri.println(Arrays.toString(opt_target_ent));
										}
										pWri.println();
									}

									// Param value
									StringBuilder param_str = new StringBuilder();
									for (double pt : point) {
										if (param_str.length() != 0) {
											param_str.append(',');
										}
										param_str.append(String.format("%f", pt));
									}
									pWri.printf("%s%d\n", OPT_OUTPUT_PREFIX_CMAP, cMap_seed);
									pWri.printf("%s%d\n", OPT_OUTPUT_PREFIX_SIMSEED, sim_seed);
									pWri.printf("%s[%s]\n", OPT_OUTPUT_PREFIX_PARAM, param_str);
									pWri.printf("%s%f\n", OPT_OUTPUT_PREFIX_RESIDUE, sqSum);
									pWri.println(str_disp.toString());
									pWri.close();

								} catch (IOException e) {
									e.printStackTrace(System.err);
									System.out.println(str_disp.toString());
								}
								return sqSum;

							}
						}
					};

					double[] start_param = Arrays.copyOf(init_param, init_param.length);
					for (int i = 2; i < result_single_param.length - 1; i++) {
						if (!((Double) result_single_param[i]).isNaN()) {
							start_param[i - 2] = result_single_param[i].doubleValue();
						}
					}

					// Run optimisation thread
					Runnable opt_runnable = new Runnable() {
						@Override
						public void run() {
							runSimplex(func, start_param, boundaries, numEval);
						}
					};

					if (exec == null) {
						exec = Executors.newFixedThreadPool(nUM_THREADS);
					}
					exec.submit(opt_runnable);
					inExec++;

					if (inExec == nUM_THREADS) {
						exec.shutdown();
						try {
							if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
								throw new InterruptedException("Time out");
							}
						} catch (InterruptedException e) {
							e.printStackTrace(System.err);
							opt_runnable.run();
						}
						exec = null;
						inExec = 0;
					}
				}

				// Run eRemainder to exec
				if (exec != null) {
					exec.shutdown();
					try {
						if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
							throw new InterruptedException("Time out");
						}
					} catch (InterruptedException e) {
						e.printStackTrace(System.err);
					}
					exec = null;
					inExec = 0;
				}

			} else {
				// Backward compatibility
				System.out.println("Optimistion using multiple seeds under same cMap - might remove support in future");

				final RandomGenerator RNG;

				final Properties PROP;

				RNG = new MersenneTwisterRandomGenerator(seed);

				PROP = prop;

				MultivariateFunction func = new MultivariateFunction() {
					@Override
					public double value(double[] point) {
						long sim_seed = RNG.nextLong();

						ContactMap[] cMap = bASE_CONTACT_MAP;
						long[] cMap_seed = bASE_CONTACT_MAP_SEED;

						if (bASE_CONTACT_MAP.length == 0) {
							cMap = new ContactMap[] { null }; // Special case for opt parameter value disp.
							cMap_seed = new long[] { 0 };
						}

						Runnable_ClusterModel_Transmission[] runnable = new Runnable_ClusterModel_Transmission[cMap.length];
						int rId = 0;

						ExecutorService exec = null;

						long tic = System.currentTimeMillis();

						for (ContactMap c : cMap) {

							final String popCompositionKey = POP_PROP_INIT_PREFIX
									+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
							int[] pop_composition = (int[]) PropValUtils
									.propStrToObject(prop.getProperty(popCompositionKey), int[].class);
							int[] cumulative_pop_composition = new int[pop_composition.length];
							int pop_offset = 0;
							for (int g = 0; g < cumulative_pop_composition.length; g++) {
								cumulative_pop_composition[g] = pop_offset + pop_composition[g];
								pop_offset += pop_composition[g];
							}

							final String riskCatListAllKey = POP_PROP_INIT_PREFIX
									+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
											+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
											+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
											+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD
											+ Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS);
							float[][] riskCatListAll = (float[][]) PropValUtils
									.propStrToObject(prop.getProperty(riskCatListAllKey), float[][].class);

							File pre_allocate_risk_file = new File(baseDir,
									String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP,
											cMap_seed[rId]));
							ArrayList<Number[]> riskGrpArr = new ArrayList<>();

							long seed = Long.parseLong(
									prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));

							boolean reallocate = true;

							if (!pre_allocate_risk_file.exists()) {

								final String time_rangeKey = POP_PROP_INIT_PREFIX
										+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
												+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
												+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);

								int[] map_time_range = (int[]) PropValUtils
										.propStrToObject(prop.getProperty(time_rangeKey), int[].class);

								Simulation_ClusterModelTransmission.fillRiskGrpArrByCasualPartnership(riskGrpArr, c,
										cumulative_pop_composition, riskCatListAll, map_time_range);

							} else {
								try {
									reallocate = Simulation_ClusterModelTransmission.loadPreallocateRiskGrp(riskGrpArr,
											baseDir, cMap_seed[rId]);
								} catch (Exception e) {
									e.printStackTrace(System.err);
								}
							}

							if (reallocate) {
								Simulation_ClusterModelTransmission.reallocateRiskGrp(riskGrpArr, cMap_seed[rId],
										cumulative_pop_composition, riskCatListAll, baseDir, seed);
							}

							runnable[rId] = new Runnable_ClusterModel_Transmission(cMap_seed[rId], sim_seed,
									pOP_COMPOSITION, c, nUM_TIME_STEPS_PER_SNAP, nUM_SNAP);
							runnable[rId].setBaseDir(baseDir);

							for (int i = RUNNABLE_OFFSET; i < RUNNABLE_OFFSET
									+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; i++) {

								String key = POP_PROP_INIT_PREFIX + Integer.toString(i);
								if (PROP.containsKey(key)) {
									runnable[rId].getRunnable_fields()[i - RUNNABLE_OFFSET] = PropValUtils
											.propStrToObject(PROP.getProperty(key),
													runnable[rId].getRunnable_fields()[i - RUNNABLE_OFFSET].getClass());
								}

							}

							runnable[rId].setSimSetting(1); // No output
							setOptParamInRunnable(runnable[rId], PROP, point, c == null);
							runnable[rId].initialse();
							runnable[rId].fillRiskCatMap(riskGrpArr);
							runnable[rId].allocateSeedInfection(sEED_INFECTION, sTART_TIME);
							rId++;
						}

						if (oPT_TARGET == null || oPT_TARGET.length == 0) {
							System.out.println("OPT_TARGET missing. Printing out runnable fields instead.");
							Runnable_ClusterModel_Transmission showRunnable = runnable[0];
							for (int i = 0; i < showRunnable.getRunnable_fields().length; i++) {
								System.out.printf("Runnable field #%d:\n", i + RUNNABLE_OFFSET);
								Object obj = showRunnable.runnable_fields[i];
								System.out.println("Entry:");
								if (obj == null) {
									System.out.println("null");
								} else if (obj instanceof Object[]) {
									System.out.println(Arrays.deepToString((Object[]) obj));
								} else if (obj instanceof float[]) {
									System.out.println(Arrays.toString((float[]) obj));
								} else if (obj instanceof double[]) {
									System.out.println(Arrays.toString((double[]) obj));
								} else {
									System.out.println(obj.toString());
								}
							}
							System.exit(1);
							return Double.NaN;

						} else {

							if (rId == 1 || nUM_THREADS <= 1) {
								for (int r = 0; r < rId; r++) {
									runnable[r].run();
								}
							} else {
								exec = Executors.newFixedThreadPool(nUM_THREADS);
								for (int r = 0; r < rId; r++) {
									exec.submit(runnable[r]);
								}
								exec.shutdown();
								try {
									if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
										throw new InterruptedException("Time out");
									}
								} catch (InterruptedException e) {
									e.printStackTrace(System.err);
									for (int r = 0; r < rId; r++) {
										runnable[r].run();
									}
								}
							}

							StringBuilder pt_str = new StringBuilder();
							for (double pt : point) {
								if (pt_str.length() != 0) {
									pt_str.append(',');
								}
								pt_str.append(String.format("%.5f", pt));
							}

							double sqSumTotal = 0;

							int start_k = 2;

							Integer[] keys = new Integer[nUM_SNAP];
							keys[0] = nUM_TIME_STEPS_PER_SNAP;
							for (int k = 1; k < keys.length; k++) {
								keys[k] = keys[k - 1] + nUM_TIME_STEPS_PER_SNAP;
							}

							for (int r = 0; r < rId; r++) {
								@SuppressWarnings("unchecked")
								HashMap<Integer, int[][]> infectious_count_map = (HashMap<Integer, int[][]>) runnable[r]
										.getSim_output()
										.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT);

								@SuppressWarnings("unchecked")
								HashMap<Integer, int[]> cumul_treatment_map = (HashMap<Integer, int[]>) runnable[r]
										.getSim_output()
										.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON);

								long cm_seed = bASE_CONTACT_MAP_SEED[r];

								StringBuilder str_disp = new StringBuilder();

								double sqSum = calculateOptFitness(point, oPT_TARGET, pOP_COMPOSITION,
										infectious_count_map, cumul_treatment_map, keys, start_k,
										String.format("CM_Seed = %d, sim_seed = %d", cm_seed, sim_seed), str_disp);

								// Display trends
								try {
									File opt_output_file = new File(baseDir,
											String.format(oPT_RESULT_FILENAME_FORMAT, cm_seed, sim_seed));
									boolean newFile = !opt_output_file.exists();

									FileWriter fWri = new FileWriter(opt_output_file, true);
									PrintWriter pWri = new PrintWriter(fWri);

									if (newFile) {
										pWri.printf("Opt target (%s in total):\n", oPT_TARGET.length);
										for (float[] opt_target_ent : oPT_TARGET) {
											pWri.println(Arrays.toString(opt_target_ent));
										}
										pWri.println();
									}

									pWri.println(pt_str);
									pWri.println(str_disp.toString());
									pWri.close();
									fWri.close();

								} catch (IOException e) {
									e.printStackTrace(System.err);
									System.out.println(str_disp.toString());

								}

								sqSumTotal += sqSum;

							}

							String outMsg;

							if (runnable.length == 1) {
								outMsg = String.format(
										"P = [%s], V = %.2e, map_seed= %d, sim_seed = %d, Time req = %.3fs\n",
										pt_str.toString(), sqSumTotal, runnable[0].getcMap_seed(),
										runnable[0].getSim_seed(), (System.currentTimeMillis() - tic) / 1000f);

							} else {
								outMsg = String.format("P = [%s], V = %f, Time req = %.3fs\n", pt_str.toString(),
										sqSumTotal, (System.currentTimeMillis() - tic) / 1000f);
							}

							try {
								File opt_output_file = new File(baseDir, FILENAME_OPT_RESULT);
								FileWriter fWri = new FileWriter(opt_output_file, true);
								PrintWriter pWri = new PrintWriter(fWri);
								pWri.print(outMsg);
								pWri.close();
								fWri.close();

							} catch (IOException ex) {
								ex.printStackTrace(System.err);
							}

							System.out.println(outMsg);

							return sqSumTotal;
						}

					}

				};
				runSimplex(func, init_param, boundaries, numEval);
			}
		}

	}

	private static void runSimplex(MultivariateFunction func, final double[] param_init,
			final double[][] param_boundaries, int numEval) {
		final double RELATIVE_TOLERANCE = 1e-5;
		final double ABSOLUTE_TOLERANCE = 1e-10;

		MultivariateFunctionMappingAdapter wrapper = new MultivariateFunctionMappingAdapter(func, param_boundaries[0],
				param_boundaries[1]);

		ObjectiveFunction objFunc = new ObjectiveFunction(wrapper);

		InitialGuess initial_guess;
		final NelderMeadSimplex simplex;

		initial_guess = new InitialGuess(wrapper.boundedToUnbounded(param_init));

		simplex = new NelderMeadSimplex(param_init.length);

		SimplexOptimizer optimizer = new SimplexOptimizer(RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);

		try {
			PointValuePair pV;
			pV = optimizer.optimize(objFunc, simplex, GoalType.MINIMIZE, initial_guess, new MaxEval(numEval));
			double[] point = wrapper.unboundedToBounded(pV.getPoint());

			StringBuilder pt_str = new StringBuilder();
			for (double pt : point) {
				if (pt_str.length() != 0) {
					pt_str.append(',');
				}
				pt_str.append(String.format("%.5f", pt));
			}

			System.out.printf("Optimisation Completed.\nP = [%s], V = %f\n", pt_str.toString(), pV.getValue());

		} catch (org.apache.commons.math3.exception.TooManyEvaluationsException ex) {
			System.out.printf("Eval limit of %d reached.\nSimplex (bounded):\n", numEval);

			PointValuePair[] res = simplex.getPoints();
			Arrays.sort(res, new Comparator<PointValuePair>() {
				@Override
				public int compare(PointValuePair o1, PointValuePair o2) {
					return Double.compare(o1.getValue(), o2.getValue());
				}

			});

			for (PointValuePair pV : res) {
				double[] point = wrapper.unboundedToBounded(pV.getPoint());

				StringBuilder pt_str = new StringBuilder();
				for (double pt : point) {
					if (pt_str.length() != 0) {
						pt_str.append(',');
					}
					pt_str.append(String.format("%.5f", pt));
				}

				System.out.printf("P = [%s], V = %f\n", pt_str.toString(), pV.getValue());

			}

		}
	}

	private static int extractContactMap(final ContactMap[] BASE_CONTACT_MAP, final long[] BASE_CONTACT_MAP_SEED,
			File[] preGenClusterFiles, final int NUM_THREADS) {
		Pattern pattern_baseCMap_filename = Pattern.compile(CMAP_REGEX_STR);
		if (NUM_THREADS > 1 || preGenClusterFiles.length > 1) {

			ExecutorService exec = null;
			@SuppressWarnings("unchecked")
			Future<ContactMap>[] cm_futures = new Future[preGenClusterFiles.length];

			exec = Executors.newFixedThreadPool(NUM_THREADS);
			int mapPt = 0;
			for (File cMap_file : preGenClusterFiles) {
				Matcher m = pattern_baseCMap_filename.matcher(cMap_file.getName());
				m.matches();
				BASE_CONTACT_MAP_SEED[mapPt] = Long.parseLong(m.group(1));

				Callable<ContactMap> cm_read_callable = Abstract_Runnable_ClusterModel
						.generateContactMapCallable(cMap_file);
				cm_futures[mapPt] = exec.submit(cm_read_callable);
				mapPt++;
			}

			try {
				exec.shutdown();
				if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
					throw new InterruptedException("Time out");
				}

				for (int c = 0; c < cm_futures.length; c++) {
					BASE_CONTACT_MAP[c] = cm_futures[c].get();
				}
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace(System.err);
			}
		}

		// Non-thread version
		int cMap_count = 0;
		for (int c = 0; c < BASE_CONTACT_MAP.length; c++) {
			if (BASE_CONTACT_MAP[c] == null) {
				Callable<ContactMap> cm_read_callable = Abstract_Runnable_ClusterModel
						.generateContactMapCallable(preGenClusterFiles[c]);
				try {
					BASE_CONTACT_MAP[c] = cm_read_callable.call();
				} catch (Exception e) {
					e.printStackTrace(System.err);
				}
			}
			if (BASE_CONTACT_MAP[c] != null) {
				Matcher m = pattern_baseCMap_filename.matcher(preGenClusterFiles[c].getName());
				m.matches();
				BASE_CONTACT_MAP_SEED[c] = Long.parseLong(m.group(1));
				BASE_CONTACT_MAP[c].setId(BASE_CONTACT_MAP_SEED[c]);
				cMap_count++;
			}
		}
		return cMap_count;
	}

	public static void partnerDistributionFit_Simplex(File baseDir, final double[] init_val_bounded,
			final double[][] boundaries) throws FileNotFoundException, IOException, InvalidPropertiesFormatException {

		final File propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);
		final File outputFile = new File(baseDir, String.format("Opt_output_%d.txt", System.currentTimeMillis()));

		if (propFile.exists()) {

			FileInputStream fIS = new FileInputStream(propFile);
			Properties prop = new Properties();
			prop.loadFromXML(fIS);
			fIS.close();

			System.out.println(String.format("Properties file < %s > loaded.", propFile.getAbsolutePath()));

			Simulation_ClusterModelGeneration sim = new Simulation_ClusterModelGeneration();
			sim.setBaseDir(baseDir);
			sim.loadProperties(prop);

			final long seed;
			final int numSnap;
			final int snapFreq;
			final int[] sampleRange;

			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
				seed = Long
						.parseLong(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
			} else {
				RandomGenerator rngBase = new MersenneTwisterRandomGenerator(System.currentTimeMillis());
				seed = rngBase.nextLong();
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])) {
				numSnap = Integer
						.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
			} else {
				numSnap = 1;
			}
			if (prop.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
				snapFreq = Integer
						.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ]));
			} else {
				snapFreq = AbstractIndividualInterface.ONE_YEAR_INT;
			}
			if (prop.containsKey("POP_PROP_INIT_PREFIX_12")) {
				sampleRange = (int[]) util.PropValUtils.propStrToObject(prop.getProperty("POP_PROP_INIT_PREFIX_12"),
						int[].class);
			} else {
				sampleRange = new int[] { 0, snapFreq * numSnap };
			}

			MultivariateFunction func = new MultivariateFunction() {

				@Override
				public double value(double[] point) {

					long tic = System.currentTimeMillis();

					Population_Bridging population = new Population_Bridging(seed);
					for (int f = 0; f < Population_Bridging.LENGTH_FIELDS_BRIDGING_POP; f++) {
						if (sim.simFields[f] != null) {
							population.getFields()[f] = sim.simFields[f];
						}
					}
					population.setPrintStatus(null);

					HashMap<String, Object> stepwiseOutput = new HashMap<>();

					double[] var_by_gender = new double[Population_Bridging.LENGTH_GENDER];
					Arrays.fill(var_by_gender, 1);
					switch (point.length) {
					case 1:
						// All
						for (int i = 0; i < var_by_gender.length; i++) {
							var_by_gender[i] = point[0];
						}
						break;
					case 2:
						// Hetro female and all msm
						var_by_gender[Population_Bridging.GENDER_FEMALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_HETRO_MALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_MSMO] = point[1];
						var_by_gender[Population_Bridging.GENDER_MSMW] = point[1];
						break;
					case 3:
						// Hetro female, MSMO and MSMW
						var_by_gender[Population_Bridging.GENDER_FEMALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_HETRO_MALE] = point[0];
						var_by_gender[Population_Bridging.GENDER_MSMO] = point[1];
						var_by_gender[Population_Bridging.GENDER_MSMW] = point[2];
						break;
					default:
						var_by_gender = Arrays.copyOf(point, Population_Bridging.LENGTH_GENDER);
					}

					int[] numInPop = (int[]) population.getFields()[Population_Bridging.FIELD_POP_COMPOSITION];
					float[] field_mean_num_partners = (float[]) population
							.getFields()[Population_Bridging.FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS];
					float[] default_mean_num_partners = Arrays.copyOf(field_mean_num_partners,
							field_mean_num_partners.length);

					int numCat = field_mean_num_partners.length / (1 + Population_Bridging.LENGTH_GENDER);

					if (field_mean_num_partners.length == Population_Bridging.LENGTH_GENDER) {
						for (int g = 0; g < field_mean_num_partners.length; g++) {
							field_mean_num_partners[g] = Math
									.round(field_mean_num_partners[g] * (float) var_by_gender[g]);
						}
					} else {
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							for (int c = 0; c < numCat; c++) {
								int pdIndex = numCat + g * numCat + c;
								field_mean_num_partners[pdIndex] = (float) (field_mean_num_partners[pdIndex]
										* var_by_gender[g]);
							}
						}
					}

					// population.setPrintStatus(System.out);
					population.initialise();

					int[] population_num_partner_in_last_12_months_total = new int[field_mean_num_partners.length];
					int sampleCount = 0;
					int startTime = Math.max(0, sampleRange[0] - 1);

					for (int s = 0; s < numSnap; s++) {
						for (int f = 0; f < snapFreq; f++) {
							if (population.getGlobalTime() == startTime) {
								population.setStepwise_output(stepwiseOutput);
							} else if (population.getGlobalTime() >= sampleRange[1]) {
								population.setStepwise_output(null);
							}

							population.advanceTimeStep(1);

							if (population.getGlobalTime() >= sampleRange[0]
									&& population.getGlobalTime() <= sampleRange[1]) {
								int[] population_num_partner_in_last_12_months = (int[]) stepwiseOutput
										.get(Population_Bridging.STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS);
								for (int i = numCat; i < population_num_partner_in_last_12_months.length; i++) {
									population_num_partner_in_last_12_months_total[i] += population_num_partner_in_last_12_months[i];
								}
								sampleCount++;
							}
						}
					}

					double residue = 0;
					double diff;
					double[] gender_weight = new double[Population_Bridging.LENGTH_GENDER];
					Arrays.fill(gender_weight, 1);
					gender_weight[Population_Bridging.GENDER_HETRO_MALE] = 0;

					for (int g = 0; g < gender_weight.length; g++) {
						if (field_mean_num_partners.length == Population_Bridging.LENGTH_GENDER) {
							diff = gender_weight[g]
									* (((double) population_num_partner_in_last_12_months_total[g]) / sampleCount
											- Math.round(default_mean_num_partners[g] * numInPop[g]));
							residue += Math.pow(diff, 2);
						} else {
							for (int c = 0; c < numCat; c++) {
								int pdIndex = numCat + g * numCat + c;
								diff = gender_weight[g]
										* (((double) population_num_partner_in_last_12_months_total[pdIndex])
												/ sampleCount
												- Math.round(default_mean_num_partners[pdIndex] * numInPop[g]));
								residue += Math.pow(diff, 2);
							}
						}
					}

					String outputString = String.format("x = %s, R = %.4f, Time req. = %.3fs", Arrays.toString(point),
							residue, (System.currentTimeMillis() - tic) / 1000f);

					System.out.println(outputString);
					try {
						PrintWriter pWri = new PrintWriter(new FileWriter(outputFile, true));
						pWri.println(outputString);
						pWri.close();

					} catch (IOException ex) {
						ex.printStackTrace(System.err);
					}

					return residue;
				}

			};

			MultivariateFunctionMappingAdapter wrapper = new MultivariateFunctionMappingAdapter(func, boundaries[0],
					boundaries[1]);
			SimplexOptimizer optimizer = new SimplexOptimizer(1e-5, 1e-10);
			NelderMeadSimplex simplex = new NelderMeadSimplex(init_val_bounded.length);
			ObjectiveFunction objFunc = new ObjectiveFunction(wrapper);

			PointValuePair var = optimizer.optimize(new MaxEval(100), objFunc, simplex, GoalType.MINIMIZE,
					new InitialGuess(wrapper.boundedToUnbounded(init_val_bounded)));

			String outputString = "Optimised value = " + Arrays.toString(wrapper.unboundedToBounded(var.getPoint()));
			System.out.println(outputString);
			try {
				PrintWriter pWri = new PrintWriter(new FileWriter(outputFile, true));
				pWri.println(outputString);
				pWri.close();

			} catch (IOException ex) {
				ex.printStackTrace(System.err);
			}

		}
	}

	public static void extractParamToOpt(Properties prop, String optParmaStr) {
		String[] parameter_settings = optParmaStr.split(",");
		String[] parameter_disp = new String[parameter_settings.length];

		for (int param_arr_index = 0; param_arr_index < parameter_settings.length; param_arr_index++) {
			String param_setting = parameter_settings[param_arr_index];
			param_setting = param_setting.replaceAll("\\s", "");
			String[] param_setting_arr = param_setting.split("_");
			int param_name_index = Integer.parseInt(param_setting_arr[0]);

			String propName = String.format("%s%d", POP_PROP_INIT_PREFIX, param_name_index);
			String propType = String.format("%s%d", Simulation_ClusterModelTransmission.POP_PROP_INIT_PREFIX_CLASS,
					param_name_index);

			String dispVal = null;
			try {
				Object paraObj;

				paraObj = PropValUtils.propStrToObject(prop.getProperty(propName),
						Class.forName(prop.getProperty(propType)));

				for (int i = 1; i < param_setting_arr.length - 1; i++) {
					// Search for the first one
					int arraySel = Integer.parseInt(param_setting_arr[i]);
					Object[] testArr = ((Object[]) paraObj);
					paraObj = null;
					for (int j = 0; j < testArr.length && paraObj == null; j++) {
						if ((arraySel & 1 << j) != 0) {
							paraObj = testArr[j];
						}
					}
				}

				// Last entry - group indices
				if (paraObj instanceof int[] || paraObj instanceof float[] || paraObj instanceof double[]) {

					int arraySel = Integer.parseInt(param_setting_arr[param_setting_arr.length - 1]);

					if (paraObj instanceof int[]) {
						int[] val_int_array = (int[]) paraObj;
						for (int i = 0; i < val_int_array.length; i++) {
							if ((arraySel & 1 << i) != 0) {
								String val = Integer.toString(val_int_array[i]);
								if (dispVal == null) {
									dispVal = val;
								} else {
									if (!dispVal.equals(val)) {
										System.err.printf("Warning - Mistmatch param values between %s and %s\n",
												dispVal, val);
									}
								}
							}
						}
					} else if (paraObj instanceof float[]) {
						float[] val_float_array = (float[]) paraObj;
						for (int i = 0; i < val_float_array.length; i++) {
							if ((arraySel & 1 << i) != 0) {
								String val = Float.toString(val_float_array[i]);
								if (dispVal == null) {
									dispVal = val;
								} else {
									if (!dispVal.equals(val)) {
										System.err.printf("Warning - Mistmatch param values between %s and %s\n",
												dispVal, val);
									}
								}
							}
						}
					} else {
						double[] val_double_array = (double[]) paraObj;
						for (int i = 0; i < val_double_array.length; i++) {
							if ((arraySel & 1 << i) != 0) {
								String val = Double.toString(val_double_array[i]);
								if (dispVal == null) {
									dispVal = val;
								} else {
									if (!dispVal.equals(val)) {
										System.err.printf("Warning - Mistmatch param values between %s and %s\n",
												dispVal, val);
									}
								}
							}
						}
					}
				} else {
					dispVal = paraObj.toString();
				}
			} catch (ClassNotFoundException ex) {
				dispVal = String.format("Not defined due to %s", ex.toString());
			}

			parameter_disp[param_arr_index] = dispVal;

		}

		System.out.printf("Param_Setting: %s\n", Arrays.toString(parameter_settings).replaceAll("\\s", ""));
		System.out.printf("Param_Display: %s\n", Arrays.toString(parameter_disp).replaceAll("\\s", ""));

	}

	public static void setOptParamInRunnable(Abstract_Runnable_ClusterModel_Transmission target_runnable, Properties prop,
			double[] point, boolean display_only) {
		String[] parameter_settings = null;
		if (prop.containsKey(OptTrendFittingFunction.POP_PROP_OPT_PARAM_FIT_SETTING)) {
			parameter_settings = prop.getProperty(OptTrendFittingFunction.POP_PROP_OPT_PARAM_FIT_SETTING).split(",");
		}
		setOptParamInRunnable(target_runnable, parameter_settings, point, display_only);
	}

	public static void setOptParamInRunnable(Abstract_Runnable_ClusterModel_Transmission target_runnable,
			String[] parameter_settings, double[] point, boolean display_only) {
		if (target_runnable instanceof Runnable_ClusterModel_Transmission) {
			setOptParamInRunnableSingleTransmission((Runnable_ClusterModel_Transmission) target_runnable,
					parameter_settings, point, display_only);
		} else {
			System.err.printf("setOptParamInRunnable for %s not supported yet. Exiting...\n",
					target_runnable.getClass().getName());
			System.exit(1);
		}
	}

	public static void setOptParamInRunnableSingleTransmission(Runnable_ClusterModel_Transmission target_runnable,
			String[] parameter_settings, double[] point, boolean display_only) {

		HashMap<Integer, Object> modified_param = new HashMap<>();

		if (parameter_settings == null || parameter_settings.length != point.length) {
			// Backward compatibility.
			System.out.println("Warning Parameter setting not used as it mismatches with number or parameters.");
			setOptParamInRunnable(target_runnable, point, display_only);
		} else {
			for (int param_arr_index = 0; param_arr_index < parameter_settings.length; param_arr_index++) {
				String param_setting = parameter_settings[param_arr_index];
				param_setting = param_setting.replaceAll("\\s", "");
				String[] param_setting_arr = param_setting.split("_");
				int param_name_index = Integer.parseInt(param_setting_arr[0]);
				Object val = target_runnable.getRunnable_fields()[param_name_index - RUNNABLE_OFFSET];
				if (val != null) {
					int setting_level = 1;
					switch (param_name_index - RUNNABLE_OFFSET) {
					case Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD:
						double[][] inf_dur = (double[][]) val;
						int site_key = Integer.parseInt(param_setting_arr[1]);
						for (int s = 0; s < inf_dur.length; s++) {
							if ((1 << s & site_key) != 0) {
								double org_mean = inf_dur[s][0];
								inf_dur[s][0] = point[param_arr_index];
								// Adjust SD based on ratio from mean
								inf_dur[s][1] = (inf_dur[s][0] / org_mean) * inf_dur[s][1];
							}
						}
						break;
					case Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM:
						double[] sought_test_param = (double[]) val;
						int index_key = Integer.parseInt(param_setting_arr[1]);
						for (int s = 0; s < sought_test_param.length; s++) {
							if ((1 << s & index_key) != 0) {
								double org_dur = sought_test_param[s];
								sought_test_param[s] = point[param_arr_index];
								if (s % 2 == 0) {
									// Adjust SD based on ratio from mean
									sought_test_param[s + 1] = (point[9] / org_dur) * sought_test_param[s + 1];
								}
							}
						}
						break;
					default:
						recursiveRunnableFieldReplace(val, param_arr_index, point, param_setting_arr, setting_level);

					}
					// Special modification for
					modified_param.put(param_name_index, val);

				} else {
					System.err.printf("Setting of parameter not supported (wrong param number of %d?). Exiting.\n",
							param_name_index);

				}
			}

			if (display_only) {
				System.out.println("Opt. parameter display:");
				Integer[] param = modified_param.keySet().toArray(new Integer[modified_param.size()]);
				Arrays.sort(param);
				for (Integer pI : param) {
					System.out.printf("POP_PROP_INIT_PREFIX_%d:\n", pI);
					Object val = modified_param.get(pI);
					System.out.println(PropValUtils.objectToPropStr(val, val.getClass()));
					System.out.println();
				}
				System.exit(0);

			}

		}
	}

	private static void recursiveRunnableFieldReplace(Object runnableField, int param_index, double[] param_val_all,
			String[] param_setting_all, int setting_level) {
		int arraySel = Integer.parseInt(param_setting_all[setting_level]);

		double offset = 0;

		if (runnableField instanceof int[] || runnableField instanceof float[] || runnableField instanceof double[]) {
			Matcher m = OptTrendFittingFunction.POP_PROP_OPT_PARAM_FIT_SETTING_DIFF_FORMAT
					.matcher(param_setting_all[param_setting_all.length - 1]);
			if (m.find()) {
				int offsetIndex = Integer.parseInt((m.group(1)));
				offset = param_val_all[offsetIndex];
			}

			if (runnableField instanceof int[]) {
				int[] val_int_array = (int[]) runnableField;
				for (int i = 0; i < val_int_array.length; i++) {
					if ((arraySel & 1 << i) != 0) {
						val_int_array[i] = (int) Math.round(offset + param_val_all[param_index]);
					}
				}
			} else if (runnableField instanceof float[]) {
				float[] val_float_array = (float[]) runnableField;
				for (int i = 0; i < val_float_array.length; i++) {
					if ((arraySel & 1 << i) != 0) {
						val_float_array[i] = (float) (offset + param_val_all[param_index]);
					}
				}
			} else if (runnableField instanceof double[]) {
				double[] val_double_array = (double[]) runnableField;
				for (int i = 0; i < val_double_array.length; i++) {
					if ((arraySel & 1 << i) != 0) {
						val_double_array[i] = offset + param_val_all[param_index];
					}
				}
			}

		} else {
			if (runnableField.getClass().isArray()) {
				Object[] obj_array = (Object[]) runnableField;
				for (int i = 0; i < obj_array.length; i++) {
					if ((arraySel & 1 << i) != 0) {
						recursiveRunnableFieldReplace(obj_array[i], param_index, param_val_all, param_setting_all,
								setting_level + 1);

					}
				}
			} else {
				System.err.printf("Class contructor for %s not supported (wrong param number?). Exiting.\n",
						runnableField.getClass().getName());
				System.exit(1);
			}

		}
	}

	private static void setOptParamInRunnable(Abstract_Runnable_ClusterModel_Transmission target_runnable,
			double[] point, boolean display_only) {
		double[][][] transmission_rate = (double[][][]) target_runnable
				.getRunnable_fields()[Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE];

		double[] sym_test_rate = (double[]) target_runnable
				.getRunnable_fields()[Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM];

		double[][] inf_dur = (double[][]) target_runnable
				.getRunnable_fields()[Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD];

		float[][] sym_rate = (float[][]) target_runnable
				.getRunnable_fields()[Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SYM_RATE];

		switch (point.length) {
		case 8:
			// TRANS_P2R, TRANS_R2P
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[0];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[1];
			// TRANS_P2O, TRANS_O2P
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[3];
			// TRANS_R2O, TRANS_O2R
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM] = new double[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[4];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[5];
			// TRANS_O2O
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[6];

			// SYM_TEST_PERIOD
			sym_test_rate[0] = point[7];
			// Adjust SD based on ratio from mean
			sym_test_rate[1] = (point[7] / 3) * 0.86 * Math.sqrt(3 * 0.86 * 0.86);

			break;
		case 10:
		case 14:
		case 15:
		case 16:
			double org_mean;
			// TRANS_P2V, TRANS_V2P
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][Abstract_Runnable_ClusterModel_Transmission.SITE_VAGINA][0] = point[0];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][Abstract_Runnable_ClusterModel_Transmission.SITE_VAGINA][1] = 0;
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_VAGINA][Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[1];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_VAGINA][Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][1] = 0;
			// TRANS_P2R, TRANS_R2P
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[3];
			// TRANS_P2O, TRANS_O2P
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[4];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[5];
			// TRANS_R2O, TRANS_O2R
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM] = new double[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[6];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[7];
			// TRANS_O2O
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Abstract_Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[8];
			// SYM_TEST_PERIOD
			org_mean = sym_test_rate[0];
			sym_test_rate[0] = point[9];
			// Adjust SD based on ratio from mean
			sym_test_rate[1] = (point[9] / org_mean) * sym_test_rate[1];

			if (point.length >= 14) {
				// Duration by site
				for (int s = 0; s < Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
					org_mean = inf_dur[s][0];
					inf_dur[s][0] = point[s + 10];
					// Adjust SD based on ratio from mean
					inf_dur[s][1] = (inf_dur[s][0] / org_mean) * inf_dur[s][1];
				}

				// Sym test adjustment for hetrosexual male
				if (point.length >= 15) {
					// Backward compatibility to single mean-sd option
					if (sym_test_rate.length < 2 * Population_Bridging.LENGTH_GENDER) {
						sym_test_rate = Arrays.copyOf(sym_test_rate,
								sym_test_rate.length * Population_Bridging.LENGTH_GENDER);
						target_runnable
								.getRunnable_fields()[Abstract_Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM] = sym_test_rate;
						for (int g = 1; g < Population_Bridging.LENGTH_GENDER; g++) {
							sym_test_rate[2 * g] = sym_test_rate[0];
							sym_test_rate[2 * g + 1] = sym_test_rate[1];
						}
					}

					// Hetro_male
					sym_test_rate[2] = point[14];
					sym_test_rate[3] = (point[14] / sym_test_rate[0]) * sym_test_rate[1];

				}

				// Sym rate for urethral infection for male
				if (point.length >= 16) {
					for (int g = 1; g < Population_Bridging.LENGTH_GENDER; g++) {
						sym_rate[g][Abstract_Runnable_ClusterModel_Transmission.SITE_PENIS] = (float) point[15];
					}
				}
			}

			break;

		default:
			System.err.printf("Optimisation: Parameter intrepretation %s not defined. Exiting...\n",
					Arrays.toString(point));
			System.exit(-1);

		}

		if (display_only) {
			System.out.println("Opt. parameter display:");

			System.out.println("RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE");
			System.out.println(Arrays.deepToString(transmission_rate));
			System.out.println();

			System.out.println("RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD");
			System.out.println(Arrays.deepToString(inf_dur));
			System.out.println();

			System.out.println("RUNNABLE_FIELD_TRANSMISSION_SYM_RATE");
			System.out.println(Arrays.deepToString(sym_rate));
			System.out.println();

			System.out.println("RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM");
			System.out.println(Arrays.toString(sym_test_rate));
			System.out.println();

			System.out.println("Opt. parameter display completed.");
			System.exit(0);
		}
	}

	public static double calculateOptFitness(double[] parameters, final float[][] opt_target,
			final int[] pop_composition, HashMap<Integer, int[][]> infectious_count_map,
			HashMap<Integer, int[]> cumul_treatment_map, Integer[] map_keys, int start_key, String simIdentifier,
			StringBuilder str_disp) {

		double sqSum = 0;

		for (int k = start_key; k < map_keys.length; k++) {

			if (str_disp != null) {
				str_disp.append(map_keys[k]);
			}
			// Number of infections
			int[][] inf_count = null;
			if (infectious_count_map != null) {
				inf_count = infectious_count_map.get(map_keys[k]);
			} else {
				System.err.printf("Warning: infection count map not defined for Sim #[%s] under parameter of %s\n",
						simIdentifier, Arrays.toString(parameters));
			}

			// Number of treatment / DX
			int[] current_treatment_count = null;
			int[] pre_treatment_count = null;
			if (cumul_treatment_map != null) {
				current_treatment_count = cumul_treatment_map.get(map_keys[k]);
				pre_treatment_count = cumul_treatment_map.get(map_keys[k - 1]);
			} else {
				System.err.printf("Warning: treatment count map not defined for Sim #[%s] under parameter of %s\n",
						simIdentifier, Arrays.toString(parameters));
			}

			if (pre_treatment_count == null && current_treatment_count != null) {
				pre_treatment_count = new int[current_treatment_count.length];
			}

			for (float[] opt_target_ent : opt_target) {
				switch ((int) opt_target_ent[OPT_TARGET_STABLE_DATA_FITTING_TYPE]) {
				case OPT_TARGET_STABLE_FITTING_TYPE_NUM_INFECTED_BY_SITE:
					if (inf_count != null) {
						float[] inf_count_total_by_site = new float[Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE];
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							if ((1 << g & (int) opt_target_ent[OPT_TARGET_STABLE_GROUPS_TO_INCLUDE]) > 0) {
								for (int s = 0; s < Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
									inf_count_total_by_site[s] += inf_count[g][s];
								}
							}
						}
						for (int s = 0; s < Abstract_Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
							float numInfected = inf_count_total_by_site[s];
							sqSum += opt_target_ent[OPT_TARGET_STABLE_OPT_WEIGHTING]
									* Math.pow(numInfected - opt_target_ent[OPT_TARGET_STABLE_TARGET_VALUES + s], 2);
							if (str_disp != null) {
								str_disp.append(',');
								str_disp.append(numInfected);
							}
						}
					} else {
						// Extinction where it shouldn't be
						boolean target_extinct = true;
						for (int s = OPT_TARGET_STABLE_TARGET_VALUES; s < opt_target_ent.length
								&& target_extinct; s++) {
							target_extinct &= opt_target_ent[OPT_TARGET_STABLE_TARGET_VALUES + s] == 0;
						}
						if (!target_extinct) {
							sqSum += Double.POSITIVE_INFINITY;
						}
					}
					break;
				case OPT_TARGET_STABLE_FITTING_TYPE_NOTIFICATIONS_BY_PERSON:
					if (current_treatment_count != null) {
						float treatment_count = 0;
						float treatment_count_denom = 0;
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							if ((1 << g & (int) opt_target_ent[OPT_TARGET_STABLE_GROUPS_TO_INCLUDE]) > 0) {
								treatment_count += current_treatment_count[g] - pre_treatment_count[g];
								treatment_count_denom += pop_composition[g];
							}
						}
						float treatment_rate = treatment_count / treatment_count_denom;
						sqSum += opt_target_ent[OPT_TARGET_STABLE_OPT_WEIGHTING]
								* Math.pow(treatment_rate - opt_target_ent[OPT_TARGET_STABLE_TARGET_VALUES], 2);

						if (str_disp != null) {
							str_disp.append(',');
							str_disp.append(treatment_rate);
						}
					} else if (opt_target_ent[OPT_TARGET_STABLE_TARGET_VALUES] != 0) {
						// Extinction where it shouldn't be
						sqSum += Double.POSITIVE_INFINITY;

					}
					break;
				default:
					System.err.printf("Warning: Opt fitting of type = %f not defined. Fitting of %s ignored/n",
							opt_target_ent[OPT_TARGET_STABLE_DATA_FITTING_TYPE], Arrays.toString(opt_target_ent));

				}

			}

			if (str_disp != null) {
				str_disp.append('\n');
			}
		}
		return sqSum;
	}

}
