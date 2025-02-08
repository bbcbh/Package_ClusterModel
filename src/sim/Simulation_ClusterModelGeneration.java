package sim;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.InvalidPropertiesFormatException;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import population.Population_Bridging_NetworkDensity;
import population.Population_Bridging_Scheduled;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import util.PersonClassifier;
import util.PropValUtils;

public class Simulation_ClusterModelGeneration implements SimulationInterface {

	public static final Object[] DEFAULT_BRIDGING_MAP_GEN_SIM_FIELDS = {};
	public static final int LENGTH_SIM_MAP_GEN_FIELD = 0;

	public static final String POP_PROP_INIT_PREFIX = "POP_PROP_INIT_PREFIX_";
	public static final String POP_PROP_INIT_PREFIX_CLASS = "POP_PROP_INIT_PREFIX_CLASS_";

	public Object[] simFields = new Object[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP + LENGTH_SIM_MAP_GEN_FIELD
			+ Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD];
	public Class<?>[] simFieldClass = new Class[simFields.length];

	protected File baseDir = null;
	protected boolean stopNextTurn = false;
	protected Properties loadedProperties = null; // From .prop file, if any
	protected ArrayList<Long> skipSeeds = null;
	protected ArrayList<Long> useSeeds = null;
	// protected transient Map<Long, ContactMap[]> contactMapSet = null;
	protected transient Map<Long, Abstract_Runnable_ClusterModel_ContactMap_Generation> runnablesMap = null;

	protected boolean printOutput = false;
	protected boolean space_save = false;

	public static final String FILENAME_FORMAT_ALL_CMAP = "All_ContactMap_%d_%d.csv";
	public static final String FILENAME_FORMAT_OUTPUT = "Output_%d.txt";

	public static final String ARG_PRINTOUTPUT = "-printOutput";
	public static final String ARG_USE_EXIST_POP = "-useExistingPop";
	public static final String ARG_SPACE_SAVE = "-space_save";
	public static final String ARG_PRE_SEED = "-preSeed=";

	public void setSpaceSave(boolean space_save) {
		this.space_save = space_save;
	}

	public void setPrintOutput(boolean printOutput) {
		this.printOutput = printOutput;
	}

	public void setSkipSeeds(ArrayList<Long> skipSeeds) {
		this.skipSeeds = skipSeeds;
		Collections.sort(this.skipSeeds);
	}

	public void setUseSeeds(ArrayList<Long> useSeeds) {
		this.useSeeds = useSeeds;
		Collections.sort(this.useSeeds);
	}

	public Simulation_ClusterModelGeneration() {

		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP + LENGTH_SIM_MAP_GEN_FIELD
				+ Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD;
		for (int i = 0; i < simFields.length; i++) {
			// All simulation levels
			if (i >= sim_offset) {
				simFields[i] = DEFAULT_BRIDGING_MAP_GEN_SIM_FIELDS[i - sim_offset];
				simFieldClass[i] = DEFAULT_BRIDGING_MAP_GEN_SIM_FIELDS[i - sim_offset].getClass();
			}
		}
	}

	public Map<Long, Abstract_Runnable_ClusterModel_ContactMap_Generation> getRunnables() {
		return runnablesMap;
	}

	/*
	 * public Map<Long, ContactMap[]> getContactMapSet() { return contactMapSet; }
	 */

	@Override
	public void loadProperties(Properties prop) {

		loadedProperties = prop;

		for (int i = 0; i < simFields.length; i++) {
			String propName = String.format("%s%d", POP_PROP_INIT_PREFIX, i);
			if (prop.containsKey(propName)) {
				String objStr = prop.getProperty(propName);
				if (simFieldClass[i] == null) {
					String str = String.format("%s%d", POP_PROP_INIT_PREFIX_CLASS, i);
					if (prop.containsKey(str)) {
						str = prop.getProperty(str);
						try {
							simFieldClass[i] = Class.forName(str);
						} catch (ClassNotFoundException e) {
							e.printStackTrace();
							simFieldClass[i] = Object.class;
						}
					}
				}
				simFields[i] = PropValUtils.propStrToObject(objStr, simFieldClass[i]);
			}

		}
	}

	@Override
	public Properties generateProperties() {
		Properties prop;
		if (loadedProperties == null) {
			prop = new Properties();
		} else {
			prop = loadedProperties;
		}

		prop.put("PROP_GENERATED_AT", Long.toString(System.currentTimeMillis()));
		for (int i = 0; i < simFields.length; i++) {
			String propName = String.format("%s%d", POP_PROP_INIT_PREFIX, i);
			String className = String.format("%s%d", POP_PROP_INIT_PREFIX_CLASS, i);
			if (simFields[i] != null) {
				prop.put(propName, PropValUtils.objectToPropStr(simFields[i], simFields[i].getClass()));
				prop.put(className, simFields[i].getClass().getName());
			}

		}
		return prop;
	}

	@Override
	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;

	}

	@Override
	public void setStopNextTurn(boolean stopNextTurn) {
		this.stopNextTurn = stopNextTurn;
	}

	@Override
	public void setSnapshotSetting(PersonClassifier[] snapshotCountClassifier, boolean[] snapshotCountAccum) {
		throw new UnsupportedOperationException("Method setSnapshotSetting not support yet.");

	}

	@Override
	public void generateOneResultSet() throws IOException, InterruptedException {
		int numThreads = Runtime.getRuntime().availableProcessors();
		int numSim = 1;
		long seed = System.currentTimeMillis();
		int numSnap = 1;
		int snapFreq = AbstractIndividualInterface.ONE_YEAR_INT;

		if (loadedProperties != null) {
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL])) {
				numThreads = Integer.parseInt(loadedProperties
						.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET])) {
				numSim = Integer.parseInt(loadedProperties
						.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
				seed = Long.parseLong(
						loadedProperties.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])) {
				numSnap = Integer.parseInt(
						loadedProperties.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
				snapFreq = Integer.parseInt(loadedProperties
						.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ]));
			}
		}

		RandomGenerator rngBase = new MersenneTwisterRandomGenerator(seed);

		boolean useParallel = numThreads > 1 && numSim > 1;

		ExecutorService executor = null;
		int numInPool = 0;
		long tic = System.currentTimeMillis();

		Abstract_Runnable_ClusterModel_ContactMap_Generation[] runnables = new Abstract_Runnable_ClusterModel_ContactMap_Generation[numSim];
		runnablesMap = new HashMap<Long, Abstract_Runnable_ClusterModel_ContactMap_Generation>();

		if (useSeeds != null) {
			numSim = useSeeds.size();
		}

		for (int i = 0; i < numSim; i++) {
			long popSeed;

			if (useSeeds != null) {
				popSeed = useSeeds.get(i);
			} else {
				popSeed = rngBase.nextLong();
			}

			if (skipSeeds == null || Collections.binarySearch(skipSeeds, popSeed) < 0) {

				Population_Bridging population = null;

				if (Population_Bridging_Scheduled.class.getName().equals(
						loadedProperties.get(SimulationInterface.PROP_NAME[SimulationInterface.PROP_POP_TYPE]))) {
					population = new Population_Bridging_Scheduled(popSeed);
					((Population_Bridging_Scheduled) population).setSpace_save(space_save);
				} else if (Population_Bridging_NetworkDensity.class.getName().equals(
						loadedProperties.get(SimulationInterface.PROP_NAME[SimulationInterface.PROP_POP_TYPE]))) {
					population = new Population_Bridging_NetworkDensity(popSeed);
					((Population_Bridging_NetworkDensity) population).setSpace_save(space_save);

				} else {
					population = null;
				}

				Abstract_Runnable_ClusterModel_ContactMap_Generation r;

				if (population != null) {
					for (int f = 0; f < Population_Bridging.LENGTH_FIELDS_BRIDGING_POP; f++) {
						if (simFields[f] != null) {
							population.getFields()[f] = simFields[f];
						}
					}
					r = new Runnable_ClusterModel_ContactMap_Generation_BridgingPop(population.getSeed());
					((Runnable_ClusterModel_ContactMap_Generation_BridgingPop) r).setPopulation(population);
				} else {
					String popType = loadedProperties
							.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_POP_TYPE]);

					if ("MutiMap".equals(popType)) {
						r = new Runnable_ClusterModel_ContactMap_Generation_MultiMap(popSeed);

					} else {
						r = null;
						System.err.printf("Error: GenMap Pop type <%s> not defined. Exiting.\n");
						System.exit(-1);
					}

				}

				runnables[i] = r;
				r.setSpace_save(space_save);
				r.setNumSnaps(numSnap);
				r.setSnapFreq(snapFreq);
				r.setBaseDir(baseDir);
				if (population != null) {
					population.setBaseDir(baseDir);
				}
				
				r.setRunnable_fields(simFields);

				/*
				for (int f = 0; f < Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD; f++) {
					if (simFields[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP + LENGTH_SIM_MAP_GEN_FIELD
							+ f] != null) {
						r.getRunnable_fields()[f] = simFields[f + LENGTH_SIM_MAP_GEN_FIELD
								+ Population_Bridging.LENGTH_FIELDS_BRIDGING_POP];
					}
				}
				*/

				runnablesMap.put(popSeed, r);

				if (printOutput) {
					PrintStream[] outputPS;

					if (numThreads < 0) {
						outputPS = new PrintStream[2];
						System.out.println("Debug: Output print to console called by PROP_USE_PARALLEL < 0");
						outputPS[1] = System.out;

					} else {
						outputPS = new PrintStream[1];
					}
					File outputFile = new File(baseDir, String.format(FILENAME_FORMAT_OUTPUT, r.getMapSeed()));
					FileOutputStream fOut = new FileOutputStream(outputFile, true);
					outputPS[0] = new PrintStream(fOut);

					outputPS[0].println(String.format("Seed = %d", r.getMapSeed()));
					runnables[i].setPrintStatus(outputPS);
				}

				if (!useParallel) {
					runnables[i].run();
				} else {
					if (executor == null) {
						executor = Executors.newFixedThreadPool(numThreads);
					}
					executor.submit(runnables[i]);
					numInPool++;
					if (numInPool == numThreads) {
						executor.shutdown();
						if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
							showStrStatus("Thread time-out!");
						}
						numInPool = 0;
						executor = null;
						System.gc();
					}
				}

				runnables[i] = null; // Free up memory?

			} else {
				showStrStatus(String.format("Contact cluster with seed of %d skipped", popSeed));
			}

		}
		if (useParallel && executor != null) {
			executor.shutdown();
			if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
				showStrStatus("Thread time-out!");
			}
			numInPool = 0;
			executor = null;
			System.gc();
		}

		/*
		 * contactMapSet = new HashMap<Long, ContactMap[]>();
		 * 
		 * for (int r = 0; r < runnables.length; r++) { if (runnables[r] != null) {
		 * contactMapSet.put(runnables[r].getPopulation().getSeed(),
		 * runnables[r].getGen_cMap()); } }
		 */

		showStrStatus(String.format("Simulation time required = %.3f s", (System.currentTimeMillis() - tic) / 1000f));
	}

	private void showStrStatus(String string) {
		System.out.println(string);

	}

	public static void launch(String[] args)
			throws InvalidPropertiesFormatException, IOException, InterruptedException {
		File baseDir = null;
		File propFile = null;
		File[] preGenClusterFile = new File[0];
		ArrayList<Long> preGenClusterSeed = new ArrayList<>();

		boolean printOut = false;
		boolean useExistingPop = false;
		boolean space_save = false;
		ArrayList<Long> preSeedList = new ArrayList<>();

		if (args.length > 0) {
			baseDir = new File(args[0]);

			for (int i = 1; i < args.length; i++) {
				if (args[i].equals(ARG_PRINTOUTPUT)) {
					printOut = true;
				}
				if (args[i].equals(ARG_USE_EXIST_POP)) {
					useExistingPop = true;
				}
				if (args[i].equals(ARG_SPACE_SAVE)) {
					space_save = true;
				}
				if (args[i].startsWith(ARG_PRE_SEED)) {
					String[] arg_seeds = args[i].substring(ARG_PRE_SEED.length()).split(",");
					for (String arg_s : arg_seeds) {
						preSeedList.add(Long.parseLong(arg_s));
					}
				}

			}

		} else {
			System.out.println(String.format("Usage: java %s PROP_FILE_DIRECTORY",
					Simulation_ClusterModelGeneration.class.getName()));
			System.exit(0);
		}

		if (baseDir != null) {
			if (baseDir.isDirectory()) {
				// Reading of PROP file
				propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);
				FileInputStream fIS = new FileInputStream(propFile);
				Properties prop = new Properties();
				prop.loadFromXML(fIS);
				fIS.close();
				System.out.println(String.format("Properties file < %s > loaded.", propFile.getAbsolutePath()));

				// Check for cluster generated previously
				final String REGEX_STR = FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");

				preGenClusterFile = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isFile() && Pattern.matches(REGEX_STR, pathname.getName());

					}
				});

				Pattern p = Pattern.compile(FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));

				for (int i = 0; i < preGenClusterFile.length; i++) {
					Matcher m = p.matcher(preGenClusterFile[i].getName());
					if (m.matches()) {
						long seed = Long.parseLong(m.group(1));
						preGenClusterSeed.add(seed);
						System.out.println(String.format("ContactMap for seed #%d located at %s", seed,
								preGenClusterFile[i].getAbsolutePath()));
					}
				}

				Simulation_ClusterModelGeneration sim = new Simulation_ClusterModelGeneration();
				sim.setBaseDir(baseDir);
				sim.setSkipSeeds(preGenClusterSeed);
				sim.loadProperties(prop);
				sim.setPrintOutput(printOut);
				sim.setSpaceSave(space_save);

				if (useExistingPop) {
					ArrayList<Long> useSeeds = new ArrayList<>();
					final String regEx_PopSnap = Abstract_Runnable_ClusterModel_ContactMap_Generation.EXPORT_POP_FILENAME
							.replaceAll("%d", "(0|-{0,1}(?!0)\\\\d+)");
					final Pattern pattern_pop_snap = Pattern.compile(regEx_PopSnap);

					File[] existingPops = baseDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							return pathname.isFile() && Pattern.matches(regEx_PopSnap, pathname.getName());

						}
					});

					for (File existPop : existingPops) {
						Matcher m = pattern_pop_snap.matcher(existPop.getName());
						if (m.matches()) {
							long seed = Long.parseLong(m.group(1));
							useSeeds.add(seed);
							System.out.printf("Attempt to continue population of seed #%d located at %s\n", seed,
									existPop.getName());
						}
					}
					sim.setUseSeeds(useSeeds);
				}

				if (!preSeedList.isEmpty()) {
					sim.setUseSeeds(preSeedList);
				}

				sim.generateOneResultSet();

			}

		}

	}

}
