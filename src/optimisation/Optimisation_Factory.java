package optimisation;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.InvalidPropertiesFormatException;
import java.util.Properties;
import java.util.concurrent.Callable;
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
import util.PropValUtils;

public class Optimisation_Factory {

	private static final String POP_PROP_INIT_PREFIX = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX;
	private static final String FILENAME_FORMAT_ALL_CMAP = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP;

	public static void stable_prevalence_by_tranmission_fit_Simplex(String[] args)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException {

		final String USAGE_INFO = "Usage: PROP_FILE_DIRECTORY (double[]) INIT_TRANSMISSION_VALUE (double[][]) BOUNDARIES <optional: NUM_EVAL (int)>";

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

		final String TARGET_PREVAL_STR = POP_PROP_INIT_PREFIX
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
			int snapFreq = 1;
			int[] pop_composition = new int[] { 500000, 500000, 20000, 20000 };
			int numThreads = Runtime.getRuntime().availableProcessors();
			int[][] target_infected = new int[Population_Bridging.LENGTH_GENDER][Runnable_ClusterModel_Transmission.LENGTH_SITE];
			int contact_map_start_time = 365;

			if (prop.containsKey(TARGET_PREVAL_STR)) {
				target_infected = (int[][]) PropValUtils.propStrToObject(prop.getProperty(TARGET_PREVAL_STR),
						int[][].class);

				// Fitting MSMO and MSMW only
				for (int g : new int[] { Population_Bridging.GENDER_FEMALE, Population_Bridging.GENDER_HETRO_MALE }) {
					for (int s = 0; s < target_infected[g].length; s++) {
						target_infected[g][s] = -1;
					}
				}

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
				snapFreq = Integer
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

			File[] preGenClusterFiles = baseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.isFile() && Pattern.matches(REGEX_STR, pathname.getName());

				}
			});

			numThreads = Math.min(numThreads, preGenClusterFiles.length);

			final RandomGenerator RNG;
			final int NUM_TIME_STEPS_PER_SNAP;
			final int SNAP_FREQ;
			final int[] POP_COMPOSITION;
			final int NUM_THREADS;
			final ContactMap[] BASE_CONTACT_MAP;
			final long[] BASE_CONTACT_MAP_SEED;
			final int[][] TARGET_INFECTED;
			final int START_TIME;

			RNG = new MersenneTwisterRandomGenerator(seed);
			NUM_TIME_STEPS_PER_SNAP = snapFreq;
			SNAP_FREQ = numSnap;
			POP_COMPOSITION = pop_composition;
			NUM_THREADS = numThreads;
			TARGET_INFECTED = target_infected;
			START_TIME = contact_map_start_time;

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
					cMap_count++;
				}
			}

			System.out.printf("%d ContactMap(s) from %s loaded. Time req. = %.3fs\n", cMap_count,
					baseDir.getAbsolutePath(), (System.currentTimeMillis() - tic) / 1000f);

			MultivariateFunction func = new MultivariateFunction() {
				@Override
				public double value(double[] point) {
					long sim_seed = RNG.nextLong();
					Runnable_ClusterModel_Transmission[] runnable = new Runnable_ClusterModel_Transmission[BASE_CONTACT_MAP.length];
					int rId = 0;

					ExecutorService exec = null;

					long tic = System.currentTimeMillis();
					for (ContactMap c : BASE_CONTACT_MAP) {
						if (c != null) {
							runnable[rId] = new Runnable_ClusterModel_Transmission(BASE_CONTACT_MAP_SEED[rId], sim_seed,
									POP_COMPOSITION, c, NUM_TIME_STEPS_PER_SNAP, SNAP_FREQ);
							runnable[rId].setBaseDir(baseDir);
							runnable[rId].setSimSetting(1); // No output

							double[][][] transmission_rate = (double[][][]) runnable[rId]
									.getRunnable_fields()[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE];

							double[] sym_test_rate = (double[]) runnable[rId]
									.getRunnable_fields()[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM];

							double[][] inf_dur = (double[][]) runnable[rId]
									.getRunnable_fields()[Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD];

							switch (point.length) {
							case 10:
							case 14:
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
								}

								break;

							default:
								System.err.printf("Optimisation: Parameter intrepretation %s not defined. Exiting...\n",
										Arrays.toString(point));
								System.exit(-1);

							}

							runnable[rId].initialse();
							runnable[rId].allocateSeedInfection(TARGET_INFECTED, START_TIME);

						}
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

					double sqSum = 0;

					int start_k = 2;

					Integer[] keys = new Integer[SNAP_FREQ];
					keys[0] = NUM_TIME_STEPS_PER_SNAP;
					for (int k = 1; k < keys.length; k++) {
						keys[k] = keys[k - 1] + NUM_TIME_STEPS_PER_SNAP;
					}

					for (int r = 0; r < rId; r++) {
						@SuppressWarnings("unchecked")
						HashMap<Integer, int[][]> infectious_count_map = (HashMap<Integer, int[][]>) runnable[r]
								.getSim_output().get(Runnable_ClusterModel_Transmission.SIM_OUTPUT_INFECTIOUS_COUNT);

						StringBuilder str_disp = new StringBuilder();

						for (int k = start_k; k < keys.length; k++) {
							str_disp.append(keys[k]);
							int[][] inf_count = infectious_count_map.get(keys[k]);

							for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
								for (int s = 0; s < Runnable_ClusterModel_Transmission.LENGTH_SITE; s++) {
									int val = 0;
									if (inf_count != null) {
										val = inf_count[g][s];
									}
									if (TARGET_INFECTED[g][s] >= 0) {
										sqSum += Math.pow(val - TARGET_INFECTED[g][s], 2);
									}
									str_disp.append(',');
									str_disp.append(val);
								}
							}
							str_disp.append('\n');
						}

						// Display trends
						try {
							File opt_output_file = new File(baseDir, String.format("Opt_trend_%d.txt", r));
							boolean newFile = !opt_output_file.exists();

							FileWriter fWri = new FileWriter(opt_output_file, true);
							PrintWriter pWri = new PrintWriter(fWri);

							if (newFile) {
								pWri.println("Target = " + Arrays.deepToString(TARGET_INFECTED));
							}

							pWri.println(pt_str);
							pWri.println(str_disp.toString());
							pWri.close();
							fWri.close();

						} catch (IOException e) {
							e.printStackTrace(System.err);
							System.out.println(str_disp.toString());

						}

					}
					
					String outMsg = String.format("P = [%s], V = %.2e, Time req = %.3fs\n", pt_str.toString(), sqSum,
							(System.currentTimeMillis() - tic) / 1000f);
					
					try {
						File opt_output_file = new File(baseDir, "Opt_res.txt");						
						FileWriter fWri = new FileWriter(opt_output_file, true);
						PrintWriter pWri = new PrintWriter(fWri);					
						pWri.println(pt_str);
						pWri.println(outMsg);
						pWri.close();
						fWri.close();
						
					}catch(IOException ex) {
						ex.printStackTrace(System.err);
					}
					
					System.out.println(outMsg);

					return sqSum;

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

}
