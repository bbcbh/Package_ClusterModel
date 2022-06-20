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
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PersonClassifier;
import util.PropValUtils;

public class Simulation_ClusterModelGeneration implements SimulationInterface {

	public static final Object[] DEFAULT_BRIDGING_MAP_GEN_SIM_FIELDS = {};
	public static final int LENGTH_BRIDGING_MAP_GEN_SIM_FIELD = 0;

	public static final String POP_PROP_INIT_PREFIX = "POP_PROP_INIT_PREFIX_";
	public static final String POP_PROP_INIT_PREFIX_CLASS = "POP_PROP_INIT_PREFIX_CLASS_";

	public Object[] simFields = new Object[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Runnable_ContactMapGeneration.LENGTH_RUNNABLE_MAP_GEN_FIELD + LENGTH_BRIDGING_MAP_GEN_SIM_FIELD];
	public Class<?>[] simFieldClass = new Class[simFields.length];

	protected File baseDir = null;
	protected boolean stopNextTurn = false;
	protected Properties loadedProperties = null; // From .prop file, if any
	protected ArrayList<Long> skipSeeds = null;
	protected transient Map<Long, ContactMap[]> contactMapSet = null;
	protected transient Map<Long, Runnable_ContactMapGeneration> runnablesMap = null;

	protected boolean printOutput = false;

	public static final String FILENAME_FORMAT_ALL_CMAP = "All_ContactMap_%d_%d.csv";

	public void setPrintOutput(boolean printOutput) {
		this.printOutput = printOutput;
	}

	public void setSkipSeeds(ArrayList<Long> skipSeeds) {
		this.skipSeeds = skipSeeds;
		Collections.sort(this.skipSeeds);
	}

	public Simulation_ClusterModelGeneration() {

		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ContactMapGeneration.LENGTH_RUNNABLE_MAP_GEN_FIELD;
		for (int i = 0; i < simFields.length; i++) {
			// All simulation levels
			if (i >= sim_offset) {
				simFields[i] = DEFAULT_BRIDGING_MAP_GEN_SIM_FIELDS[i - sim_offset];
				simFieldClass[i] = DEFAULT_BRIDGING_MAP_GEN_SIM_FIELDS[i - sim_offset].getClass();
			}
		}
	}

	public Map<Long, Runnable_ContactMapGeneration> getRunnables() {
		return runnablesMap;
	}

	public Map<Long, ContactMap[]> getContactMapSet() {
		return contactMapSet;
	}

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

		Runnable_ContactMapGeneration[] runnables = new Runnable_ContactMapGeneration[numSim];
		runnablesMap = new HashMap<Long, Runnable_ContactMapGeneration>();

		for (int i = 0; i < numSim; i++) {
			long popSeed = rngBase.nextLong();

			if (skipSeeds == null || Collections.binarySearch(skipSeeds, popSeed) < 0) {

				Population_Bridging population = new Population_Bridging(popSeed);
				for (int f = 0; f < Population_Bridging.LENGTH_FIELDS_BRIDGING_POP; f++) {
					if (simFields[f] != null) {
						population.getFields()[f] = simFields[f];
					}
				}

				Runnable_ContactMapGeneration r = new Runnable_ContactMapGeneration();

				runnables[i] = r;
				r.setPopulation(population);
				r.setNumSnaps(numSnap);
				r.setSnapFreq(snapFreq);
				r.setBaseDir(baseDir);

				for (int f = 0; f < Runnable_ContactMapGeneration.LENGTH_RUNNABLE_MAP_GEN_FIELD; f++) {
					if (simFields[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP + f] != null) {
						r.getRunnable_fields()[f] = simFields[f + Population_Bridging.LENGTH_FIELDS_BRIDGING_POP];
					}
				}

				runnablesMap.put(r.getPopulation().getSeed(), r);

				if (printOutput) {
					PrintStream outputPS;
					if (numThreads < 0) { 
						System.out.println("Debug: Output print to console called by PROP_USE_PARALLEL < 0");
						outputPS = System.out;

					} else {
						File outputFile = new File(baseDir,
								String.format("Output_%d.txt", r.getPopulation().getSeed()));
						FileOutputStream fOut = new FileOutputStream(outputFile, true);
						outputPS = new PrintStream(fOut);
					}
					outputPS.println(String.format("Seed = %d", r.getPopulation().getSeed()));
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
					}
				}
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
		}

		contactMapSet = new HashMap<Long, ContactMap[]>();

		for (int r = 0; r < runnables.length; r++) {
			contactMapSet.put(runnables[r].getPopulation().getSeed(), runnables[r].getGen_cMap());
		}

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

		if (args.length > 0) {
			baseDir = new File(args[0]);

			for (int i = 1; i < args.length; i++) {
				if (args[i].equals("-printOutput")) {
					printOut = true;
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
				final String REGEX_STR = FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "-{0,1}(?!0)\\\\d+");

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
				sim.generateOneResultSet();

			}

		}

	}

}
