package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PersonClassifier;
import util.PropValUtils;

public class Simulation_ClusterModelTransmission implements SimulationInterface {

	protected Properties loadedProperties = null; // From .prop file, if any
	protected boolean stopNextTurn = false;
	protected File baseDir = null;
	protected ContactMap baseContactMap;

	public static final String POP_PROP_INIT_PREFIX = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX;
	public static final String POP_PROP_INIT_PREFIX_CLASS = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX_CLASS;

	public static final String FILENAME_FORMAT_ALL_CMAP = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP;

	public static final Object[] DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS = {
			// BRIDGING_MAP_TRANS_SIM_FIELD_SEED_INFECTION
			// Usage: [gender_type][site] = probability (out of total)
			new float[][] {},

	};

	public static final int BRIDGING_MAP_TRANS_SIM_FIELD_SEED_INFECTION = 0;
	public static final int LENGTH_BRIDGING_MAP_TRANS_SIM_FIELD = BRIDGING_MAP_TRANS_SIM_FIELD_SEED_INFECTION  +1;

	public Object[] simFields = new Object[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Runnable_ContactMapGeneration.LENGTH_RUNNABLE_MAP_GEN_FIELD
			+ Simulation_ClusterModelGeneration.LENGTH_BRIDGING_MAP_GEN_SIM_FIELD
			+ LENGTH_BRIDGING_MAP_TRANS_SIM_FIELD];
	public Class<?>[] simFieldClass = new Class[simFields.length];

	public Simulation_ClusterModelTransmission() {
		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ContactMapGeneration.LENGTH_RUNNABLE_MAP_GEN_FIELD
				+ Simulation_ClusterModelGeneration.LENGTH_BRIDGING_MAP_GEN_SIM_FIELD;
		for (int i = 0; i < simFields.length; i++) {
			// All simulation levels
			if (i >= sim_offset) {
				simFields[i] = DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS[i - sim_offset];
				simFieldClass[i] = DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS[i - sim_offset].getClass();
			}
		}
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
	public void setStopNextTurn(boolean stopNextTurn) {
		this.stopNextTurn = stopNextTurn;

	}

	@Override
	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
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
		return prop;
	}

	@Override
	public void setSnapshotSetting(PersonClassifier[] snapshotCountClassifier, boolean[] snapshotCountAccum) {
		throw new UnsupportedOperationException("Method setSnapshotSetting not support yet.");

	}

	public void setBaseContactMap(ContactMap baseContactMap) {
		this.baseContactMap = baseContactMap;
	}

	@Override
	public void generateOneResultSet() throws IOException, InterruptedException {

		int numThreads = Runtime.getRuntime().availableProcessors();
		int numSim = 1;
		long seed = System.currentTimeMillis();
		int numSnap = 1;
		int snapFreq = 1;
		int[] pop_composition = new int[] { 500000, 500000, 20000, 20000 };
		int startTime = 2781;

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

			String popCompositionKey = POP_PROP_INIT_PREFIX
					+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
			if (loadedProperties.containsKey(popCompositionKey)) {
				pop_composition = (int[]) PropValUtils.propStrToObject(loadedProperties.getProperty(popCompositionKey),
						int[].class);
			}
			
			String contactMapRange =  POP_PROP_INIT_PREFIX
					+ Integer.toString(Runnable_ContactMapGeneration.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);
			
			if (loadedProperties.containsKey(contactMapRange)) {
				startTime = ((int[]) PropValUtils.propStrToObject(loadedProperties.getProperty(contactMapRange),
						int[].class))[0];
			}
		}

		// Map stat
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] personStat = new ArrayList[Population_Bridging.LENGTH_GENDER];
		int[] cumul_pop = new int[Population_Bridging.LENGTH_GENDER];

		int offset = 0;

		for (int g = 0; g < cumul_pop.length; g++) {
			cumul_pop[g] = offset + pop_composition[g];
			offset += pop_composition[g];
			personStat[g] = new ArrayList<Integer>();
		}

		for (Integer v : baseContactMap.vertexSet()) {
			int g = Runnable_ContactMapTransmission.getGenderType(v, cumul_pop);
			personStat[g].add(v);
		}
		
		int[] personCount = new int[personStat.length];		
		for(int i = 0; i < personCount.length; i++) {
			personCount[i] = personStat[i].size();
		}
		System.out.println(String.format("Number of indivduals in contact map = %s", Arrays.toString(personCount)));

		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ContactMapGeneration.LENGTH_RUNNABLE_MAP_GEN_FIELD
				+ Simulation_ClusterModelGeneration.LENGTH_BRIDGING_MAP_GEN_SIM_FIELD;
		
		float[][] seedInfectParam = (float[][]) simFields[sim_offset + BRIDGING_MAP_TRANS_SIM_FIELD_SEED_INFECTION];
		int[][] seedInfectNum = new int[Population_Bridging.LENGTH_GENDER][Runnable_ContactMapTransmission.SITE_LENGTH];

		for (int g = 0; g < seedInfectParam.length; g++) {
			for (int s = 0; s < seedInfectParam[g].length; s++) {
				float seedProp = seedInfectParam[g][s];
				if (seedProp < 1) {
					seedInfectNum[g][s] = Math.round(personStat[g].size() * seedProp);
				}else {
					seedInfectNum[g][s] = Math.round(seedInfectParam[g][s]);
				}
			}
		}

		RandomGenerator rngBase = new MersenneTwisterRandomGenerator(seed);
		boolean useParallel = numThreads > 1 && numSim > 1;

		ExecutorService exec = null;
		int inExec = 0;

		Runnable_ContactMapTransmission[] runnable = new Runnable_ContactMapTransmission[numSim];

		for (int s = 0; s < numSim; s++) {
			long simSeed = rngBase.nextLong();

			runnable[s] = new Runnable_ContactMapTransmission(simSeed, pop_composition, baseContactMap,
					numSnap * snapFreq);

			runnable[s].initialse();

			// Add infected
			for (int gender = 0; gender < Population_Bridging.LENGTH_GENDER; gender++) {
				for (int site = 0; site < Runnable_ContactMapTransmission.SITE_LENGTH; site++) {
					if (seedInfectNum[gender][site] > 0) {						
						Integer[] seedInf = util.ArrayUtilsRandomGenerator.randomSelect(
								personStat[gender].toArray(new Integer[personStat[gender].size()]),
								seedInfectNum[gender][site], rngBase);
						
						for(Integer inf: seedInf) {
							runnable[s].addInfected(inf, site, startTime + 180);													
						}
						System.out.println(String.format("Seeding %s of gender #%d at site #%d", 
								Arrays.toString(seedInf), gender, site));						
					}
				}
			}
			if (useParallel) {
				if (exec == null) {
					exec = Executors.newFixedThreadPool(numThreads);
				}
				exec.submit(runnable[s]);

				inExec++;
				if (inExec == numThreads) {
					exec.shutdown();
					if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
						System.err.println("Thread time-out!");
					}
					inExec = 0;
					exec = null;
				}
			} else {
				runnable[s].run();
			}
		}
		if (exec != null && inExec != 0) {
			exec.shutdown();
			if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
				System.err.println("Thread time-out!");
			}
			inExec = 0;
			exec = null;
		}
	}

	public static void launch(String[] args) throws IOException, InterruptedException {
		File baseDir = null;
		File propFile = null;

		File[] preGenClusterMap = new File[0];

		if (args.length > 0) {
			baseDir = new File(args[0]);
		} else {
			System.out.println(String.format("Usage: java %s PROP_FILE_DIRECTORY",
					Simulation_ClusterModelTransmission.class.getName()));
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

				// Check for contact cluster generated

				final String REGEX_STR = FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "\\\\d+");

				preGenClusterMap = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isFile() && Pattern.matches(REGEX_STR, pathname.getName());

					}
				});

				for (int i = 0; i < preGenClusterMap.length; i++) {
					System.out.println(String.format("Running simulation set based on ContactMap located at %s",
							preGenClusterMap[i].getAbsolutePath()));
					StringWriter cMap_str = new StringWriter();
					PrintWriter pWri = new PrintWriter(cMap_str);

					BufferedReader reader = new BufferedReader(new FileReader(preGenClusterMap[i]));
					String line;

					while ((line = reader.readLine()) != null) {
						pWri.println(line);
					}

					pWri.close();
					reader.close();

					Simulation_ClusterModelTransmission sim = new Simulation_ClusterModelTransmission();
					sim.setBaseDir(baseDir);
					sim.loadProperties(prop);
					sim.setBaseContactMap(ContactMap.ContactMapFromFullString(cMap_str.toString()));
					sim.generateOneResultSet();

				}

			}
		}

	}
}
