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
import sim.Runnable_ClusterModel_ContactMap_Generation;
import sim.Runnable_ClusterModel_Transmission;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelGeneration;
import sim.Simulation_ClusterModelTransmission;
import util.ArrayUtilsRandomGenerator;
import util.PropValUtils;

public class Optimisation_Factory {

	private static final String POP_PROP_INIT_PREFIX = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX;
	private static final String FILENAME_FORMAT_ALL_CMAP = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP;

	public static final String POP_PROP_OPT_TARGET = "POP_PROP_OPT_TARGET";

	// Format:
	// float[][] opt_target = new float[NUM_DATA_TO_FIT]
	// { DATA_FITTING_TYPE, GROUPS_TO_INCLUDE, OPT_WEIGHTING, OPT_TARGET_VALUES_0,
	// ...}

	private static final int OPT_TARGET_FITTING_TYPE_NUM_INFECTED_BY_SITE = 0;
	private static final int OPT_TARGET_FITTING_TYPE_NOTIFICATIONS_BY_PERSON = 1;

	private static final int OPT_TARGET_FITTING_TYPE = 0;
	private static final int OPT_TARGET_GROUPS_TO_INCLUDE = OPT_TARGET_FITTING_TYPE + 1;
	private static final int OPT_TARGET_OPT_WEIGHTING = OPT_TARGET_GROUPS_TO_INCLUDE + 1;
	private static final int OPT_TARGET_TARGET_VALUES_0 = OPT_TARGET_OPT_WEIGHTING + 1;

	public static void stable_prevalence_by_tranmission_fit_GA(String[] args)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException, InterruptedException {
		final String USAGE_INFO = "Usage: PROP_FILE_DIRECTORY INIT_PARAM_VALUE (double[]) BOUNDARIES (double[][]) GA_POP_SIZE (int) "
				+ "<-ta TOURNAMENT_ARITY (int)> <-mr MUTATION_RATE (float)> <-mg MAX_GENERATION (int)> <-exportAll true|false> "
				+ "<-useInitalValues false|true> <-useCMapMapping true|false> <-useCMapEdgeList true|false>"
				+ "<-singleExecutor false|true>";

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
			boolean useSingleExecutor = false;

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
				if("-singleExecutor".equals(args[a])) {
					useSingleExecutor = Boolean.parseBoolean(args[a + 1]);
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
			final File GA_POP_FILE = new File(baseDir, "OptRes_GA_Pop.csv");
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
				
				System.out.printf("# processors available = %d. # threads used = %d.\n", Runtime.getRuntime().availableProcessors(), numThreads);

				// Check for contact cluster generated

				final String filename_cmap_regex = FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");
				final Pattern pattern_baseCMap_filename = Pattern.compile(filename_cmap_regex);

				File[] preGenClusterFiles = contactMapDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isFile() && Pattern.matches(filename_cmap_regex, pathname.getName());

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

				if (GA_POP_FILE.exists()) {
					// Fill previous pop file
					BufferedReader reader = new BufferedReader(new FileReader(GA_POP_FILE));
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
						if(updateGA) {
							PrintWriter pWri = new PrintWriter(GA_POP_FILE);
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
					if (ga_population.get(GA_MAX_POP_SIZE - 1)[GA_ENT_FITNESS].doubleValue()
							- ga_population.get(0)[GA_ENT_FITNESS].doubleValue() < ABSOLUTE_TOLERANCE) {
						System.out.printf(
								"The difference between best and worst fit in GA_POPULATION is less than ABSOLUTE_TOLERANCE of %.f.\n",
								ABSOLUTE_TOLERANCE);
						exiting = true;
					} else if (num_gen > max_num_generation) {
						System.out.printf("Maximum number of generations (%d) reached.\n", (int) max_num_generation);
						exiting = true;
					} else {
						// Populate GA_POPULATION fitness

						ExecutorService exec = Executors.newFixedThreadPool(numThreads);
						int threadCount = 0;

						for (Number[] ga_ent : ga_population) {
							if (Double.isNaN(ga_ent[GA_ENT_FITNESS].doubleValue())) {

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

											Runnable_ClusterModel_Transmission runnable = new Runnable_ClusterModel_Transmission(
													(long) ga_ent[GA_ENT_CMAP_SEED], (long) ga_ent[GA_ENT_SIM_SEED],
													POP_COMPOSITION, cmap, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP);
											runnable.setBaseDir(baseDir);

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
													runnable.getRunnable_fields()[i - runnnable_offset] = PropValUtils
															.propStrToObject(PROP.getProperty(key),
																	runnable.getRunnable_fields()[i - runnnable_offset]
																			.getClass());
												}
											}

											runnable.setSimSetting(1); // No output
											Number[] param_number = Arrays.copyOfRange(ga_ent, GA_ENT_PARM_START,
													ga_ent.length);
											double[] param_double = new double[param_number.length];
											for (int i = 0; i < param_number.length; i++) {
												param_double[i] = param_number[i].doubleValue();
											}

											setOptParamInRunnable(runnable, param_double, false);
											runnable.initialse();
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
													POP_COMPOSITION, infectious_count_map, cumul_treatment_map, keys,
													start_k,
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
														FileWriter export_all_fWri = new FileWriter(GA_ALL_FILE, true);
														PrintWriter export_all_pWri = new PrintWriter(export_all_fWri);
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
								
								if(!useSingleExecutor) {
									if(threadCount == numThreads) {
										exec.shutdown();
										if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
											System.err.println("Thread time-out!");
										}										
										exec = Executors.newFixedThreadPool(numThreads);
										threadCount = 0;
									}																		
									
								}
								
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

						// Fill the rest with next generation
						int next_gen_start = tournament_arity;

						for (int g = next_gen_start; g < GA_MAX_POP_SIZE; g++) {
							Number[] ga_ent = ga_population.get(g);
							ga_ent[GA_ENT_FITNESS] = Double.NaN;
							ga_ent[GA_ENT_CMAP_SEED] = BASE_CONTACT_MAP_SEED[RNG.nextInt(BASE_CONTACT_MAP_SEED.length)];
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

						// Export GA_POPULATION
						PrintWriter pWri = new PrintWriter(GA_POP_FILE);
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
								dateFormat.format(new Date(System.currentTimeMillis())), GA_POP_FILE.getAbsolutePath());

						num_gen++;
					}
				}

			}

		}

	}

	public static void stable_prevalence_by_tranmission_fit_Simplex(String[] args)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException {

		final String USAGE_INFO = "Usage: PROP_FILE_DIRECTORY INIT_PARAM_VALUE (double[]) BOUNDARIES (double[][])  <optional: NUM_EVAL (int)>";

		int numEval = 100;
		if (args.length < 3) {
			System.out.println(USAGE_INFO);
			System.exit(0);
		} else {
			File baseDir = new File(args[0]);
			double[] init_transmissionProb = (double[]) PropValUtils.propStrToObject(args[1], double[].class);
			double[][] boundaries = (double[][]) PropValUtils.propStrToObject(args[2], double[][].class);
			if (args.length > 3) {
				numEval = (int) Integer.parseInt(args[3]);
			}

			stable_prevalence_by_tranmission_fit_Simplex(baseDir, init_transmissionProb, boundaries, numEval);
		}

	}

	public static void stable_prevalence_by_tranmission_fit_Simplex(File baseDir, final double[] init_transmissionProb,
			final double[][] boundaries, int numEval)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException {

		final double RELATIVE_TOLERANCE = 1e-5;
		final double ABSOLUTE_TOLERANCE = 1e-10;
		final File propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);

		final String PROP_SEED_INFECTION = POP_PROP_INIT_PREFIX
				+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
						+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
						+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
						+ Simulation_ClusterModelTransmission.SIM_FIELD_SEED_INFECTION); // "POP_PROP_INIT_PREFIX_14";

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
				contactMapDir = new File(prop.getProperty(Simulation_ClusterModelTransmission.PROP_CONTACT_MAP_LOC));
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

			// Check for contact cluster generated

			final String REGEX_STR = FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");

			File[] preGenClusterFiles = contactMapDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isFile() && Pattern.matches(REGEX_STR, pathname.getName());

				}
			});

			numThreads = Math.min(numThreads, preGenClusterFiles.length);

			final RandomGenerator RNG;
			final int NUM_TIME_STEPS_PER_SNAP;
			final int NUM_SNAP;
			final int[] POP_COMPOSITION;
			final int NUM_THREADS;
			final ContactMap[] BASE_CONTACT_MAP;
			final long[] BASE_CONTACT_MAP_SEED;
			final float[][] OPT_TARGET;
			final int START_TIME;
			final int[][] SEED_INFECTION;
			final Properties PROP;

			RNG = new MersenneTwisterRandomGenerator(seed);
			NUM_TIME_STEPS_PER_SNAP = num_time_steps_per_snap;
			NUM_SNAP = numSnap;
			POP_COMPOSITION = pop_composition;
			NUM_THREADS = numThreads;
			OPT_TARGET = opt_target;
			START_TIME = contact_map_start_time;
			SEED_INFECTION = seed_infection;
			PROP = prop;

			BASE_CONTACT_MAP = new ContactMap[preGenClusterFiles.length];
			BASE_CONTACT_MAP_SEED = new long[BASE_CONTACT_MAP.length];
			Pattern pattern_baseCMap_filename = Pattern.compile(REGEX_STR);

			long tic = System.currentTimeMillis();

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

			System.out.printf("%d ContactMap(s) from %s loaded. Time req. = %.3fs\n", cMap_count,
					contactMapDir.getAbsolutePath(), (System.currentTimeMillis() - tic) / 1000f);

			MultivariateFunction func = new MultivariateFunction() {
				@Override
				public double value(double[] point) {
					long sim_seed = RNG.nextLong();

					ContactMap[] cMap = BASE_CONTACT_MAP;
					long[] cMap_seed = BASE_CONTACT_MAP_SEED;

					if (BASE_CONTACT_MAP.length == 0) {
						cMap = new ContactMap[] { null }; // Special case for opt parameter value disp.
						cMap_seed = new long[] { 0 };
					}

					Runnable_ClusterModel_Transmission[] runnable = new Runnable_ClusterModel_Transmission[cMap.length];
					int rId = 0;

					ExecutorService exec = null;

					long tic = System.currentTimeMillis();

					for (ContactMap c : cMap) {
						runnable[rId] = new Runnable_ClusterModel_Transmission(cMap_seed[rId], sim_seed,
								POP_COMPOSITION, c, NUM_TIME_STEPS_PER_SNAP, NUM_SNAP);
						runnable[rId].setBaseDir(baseDir);

						int runnnable_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
								+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
								+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
								+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD;

						for (int i = runnnable_offset; i < runnnable_offset
								+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; i++) {

							String key = POP_PROP_INIT_PREFIX + Integer.toString(i);
							if (PROP.containsKey(key)) {
								runnable[rId].getRunnable_fields()[i - runnnable_offset] = PropValUtils.propStrToObject(
										PROP.getProperty(key),
										runnable[rId].getRunnable_fields()[i - runnnable_offset].getClass());
							}

						}

						runnable[rId].setSimSetting(1); // No output
						setOptParamInRunnable(runnable[rId], point, c == null);
						runnable[rId].initialse();
						runnable[rId].allocateSeedInfection(SEED_INFECTION, START_TIME);

						rId++;
					}
					if (rId == 1 || NUM_THREADS <= 1) {
						for (int r = 0; r < rId; r++) {
							runnable[r].run();
						}
					} else {
						exec = Executors.newFixedThreadPool(NUM_THREADS);
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

					Integer[] keys = new Integer[NUM_SNAP];
					keys[0] = NUM_TIME_STEPS_PER_SNAP;
					for (int k = 1; k < keys.length; k++) {
						keys[k] = keys[k - 1] + NUM_TIME_STEPS_PER_SNAP;
					}

					for (int r = 0; r < rId; r++) {
						@SuppressWarnings("unchecked")
						HashMap<Integer, int[][]> infectious_count_map = (HashMap<Integer, int[][]>) runnable[r]
								.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT);

						@SuppressWarnings("unchecked")
						HashMap<Integer, int[]> cumul_treatment_map = (HashMap<Integer, int[]>) runnable[r]
								.getSim_output()
								.get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_CUMUL_TREATMENT_BY_PERSON);

						long cm_seed = BASE_CONTACT_MAP_SEED[r];

						StringBuilder str_disp = new StringBuilder();

						double sqSum = calculateOptFitness(point, OPT_TARGET, POP_COMPOSITION, infectious_count_map,
								cumul_treatment_map, keys, start_k,
								String.format("CM_Seed = %d, sim_seed = %d", cm_seed, sim_seed), str_disp);

						// Display trends
						try {
							File opt_output_file = new File(baseDir, String.format("Opt_trend_%d.txt", r));
							boolean newFile = !opt_output_file.exists();

							FileWriter fWri = new FileWriter(opt_output_file, true);
							PrintWriter pWri = new PrintWriter(fWri);

							if (newFile) {
								pWri.printf("Opt target (%s in total):\n", OPT_TARGET.length);
								for (float[] opt_target_ent : OPT_TARGET) {
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
						outMsg = String.format("P = [%s], V = %.2e, map_seed= %d, sim_seed = %d, Time req = %.3fs\n",
								pt_str.toString(), sqSumTotal, runnable[0].getcMap_seed(), runnable[0].getSim_seed(),
								(System.currentTimeMillis() - tic) / 1000f);

					} else {
						outMsg = String.format("P = [%s], V = %.2e, Time req = %.3fs\n", pt_str.toString(), sqSumTotal,
								(System.currentTimeMillis() - tic) / 1000f);
					}

					try {
						File opt_output_file = new File(baseDir, "Opt_res.txt");
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

			};

			MultivariateFunctionMappingAdapter wrapper = new MultivariateFunctionMappingAdapter(func, boundaries[0],
					boundaries[1]);

			ObjectiveFunction objFunc = new ObjectiveFunction(wrapper);

			InitialGuess initial_guess;
			final NelderMeadSimplex simplex;

			initial_guess = new InitialGuess(wrapper.boundedToUnbounded(init_transmissionProb));

			simplex = new NelderMeadSimplex(init_transmissionProb.length);

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

				System.out.printf("Optimisation Completed.\nP = [%s], V = %.2e\n", pt_str.toString(), pV.getValue());

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

					System.out.printf("P = [%s], V = %.2e\n", pt_str.toString(), pV.getValue());

				}

			}

		}

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

	public static void setOptParamInRunnable(Runnable_ClusterModel_Transmission targer_runnable, double[] point,
			boolean display_only) {
		double[][][] transmission_rate = (double[][][]) targer_runnable
				.getRunnable_fields()[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE];

		double[] sym_test_rate = (double[]) targer_runnable
				.getRunnable_fields()[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM];

		double[][] inf_dur = (double[][]) targer_runnable
				.getRunnable_fields()[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD];

		float[][] sym_rate = (float[][]) targer_runnable
				.getRunnable_fields()[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SYM_RATE];

		switch (point.length) {
		case 8:
			// TRANS_P2R, TRANS_R2P
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_PENIS][Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[0];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_RECTUM][Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[1];
			// TRANS_P2O, TRANS_O2P
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_PENIS][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[3];
			// TRANS_R2O, TRANS_O2R
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_RECTUM][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_RECTUM] = new double[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_RECTUM][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[4];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[5];
			// TRANS_O2O
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[6];

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
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_PENIS][Runnable_ClusterModel_Transmission.SITE_VAGINA][0] = point[0];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_PENIS][Runnable_ClusterModel_Transmission.SITE_VAGINA][1] = 0;
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_VAGINA][Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[1];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_VAGINA][Runnable_ClusterModel_Transmission.SITE_PENIS][1] = 0;
			// TRANS_P2R, TRANS_R2P
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_PENIS][Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_RECTUM][Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[3];
			// TRANS_P2O, TRANS_O2P
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_PENIS][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[4];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_PENIS][0] = point[5];
			// TRANS_R2O, TRANS_O2R
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_RECTUM][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_RECTUM] = new double[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_RECTUM][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[6];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_RECTUM][0] = point[7];
			// TRANS_O2O
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX] = new double[2];
			transmission_rate[Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][Runnable_ClusterModel_Transmission.SITE_OROPHARYNX][0] = point[8];
			// SYM_TEST_PERIOD
			org_mean = sym_test_rate[0];
			sym_test_rate[0] = point[9];
			// Adjust SD based on ratio from mean
			sym_test_rate[1] = (point[9] / org_mean) * sym_test_rate[1];

			if (point.length >= 14) {
				// Duration by site
				for (int s = 0; s < Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
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
						targer_runnable.runnable_fields[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM] = sym_test_rate;
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
						sym_rate[g][Runnable_ClusterModel_Transmission.SITE_PENIS] = (float) point[15];
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
				switch ((int) opt_target_ent[OPT_TARGET_FITTING_TYPE]) {
				case OPT_TARGET_FITTING_TYPE_NUM_INFECTED_BY_SITE:
					if (inf_count != null) {
						float[] inf_count_total_by_site = new float[Runnable_ClusterModel_Transmission.LENGTH_SITE];
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							if ((1 << g & (int) opt_target_ent[OPT_TARGET_GROUPS_TO_INCLUDE]) > 0) {
								for (int s = 0; s < Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
									inf_count_total_by_site[s] += inf_count[g][s];
								}
							}
						}
						for (int s = 0; s < Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
							float numInfected = inf_count_total_by_site[s];
							sqSum += opt_target_ent[OPT_TARGET_OPT_WEIGHTING]
									* Math.pow(numInfected - opt_target_ent[OPT_TARGET_TARGET_VALUES_0 + s], 2);
							if (str_disp != null) {
								str_disp.append(',');
								str_disp.append(numInfected);
							}
						}
					} else {
						// Extinction where it shouldn't be
						boolean target_extinct = true;
						for (int s = OPT_TARGET_TARGET_VALUES_0; s < opt_target_ent.length && target_extinct; s++) {
							target_extinct &= opt_target_ent[OPT_TARGET_TARGET_VALUES_0 + s] == 0;
						}
						if (!target_extinct) {
							sqSum += Double.POSITIVE_INFINITY;
						}
					}
					break;
				case OPT_TARGET_FITTING_TYPE_NOTIFICATIONS_BY_PERSON:
					if (current_treatment_count != null) {
						float treatment_count = 0;
						float treatment_count_denom = 0;
						for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
							if ((1 << g & (int) opt_target_ent[OPT_TARGET_GROUPS_TO_INCLUDE]) > 0) {
								treatment_count += current_treatment_count[g] - pre_treatment_count[g];
								treatment_count_denom += pop_composition[g];
							}
						}
						float treatment_rate = treatment_count / treatment_count_denom;
						sqSum += opt_target_ent[OPT_TARGET_OPT_WEIGHTING]
								* Math.pow(treatment_rate - opt_target_ent[OPT_TARGET_TARGET_VALUES_0], 2);

						if (str_disp != null) {
							str_disp.append(',');
							str_disp.append(treatment_rate);
						}
					} else if (opt_target_ent[OPT_TARGET_TARGET_VALUES_0] != 0) {
						// Extinction where it shouldn't be
						sqSum += Double.POSITIVE_INFINITY;

					}
					break;
				default:
					System.err.printf("Warning: Opt fitting of type = %f not defined. Fitting of %s ignored/n",
							opt_target_ent[OPT_TARGET_FITTING_TYPE], Arrays.toString(opt_target_ent));

				}

			}

			if (str_disp != null) {
				str_disp.append('\n');
			}
		}
		return sqSum;
	}

}
