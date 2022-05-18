package sim;

import java.io.File;
import java.io.IOException;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

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

	
	public void setBaseContactMap(ContactMap baseContactMap) {
		this.baseContactMap = baseContactMap;
	}

	@Override
	public void loadProperties(Properties prop) {
		loadedProperties = prop;
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

	@Override
	public void generateOneResultSet() throws IOException, InterruptedException {

		int numThreads = Runtime.getRuntime().availableProcessors();
		int numSim = 1;
		long seed = System.currentTimeMillis();
		int numSnap = 1;
		int snapFreq = 1;
		int[] pop_composition = new int[] {500000, 500000, 20000, 20000};

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
			
			String popCompositionKey = POP_PROP_INIT_PREFIX + Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);			
			if(loadedProperties.containsKey(popCompositionKey)){
				pop_composition = (int[]) PropValUtils
						.propStrToObject(loadedProperties.getProperty(popCompositionKey), int[].class);
			}
		}

		RandomGenerator rngBase = new MersenneTwisterRandomGenerator(seed);
		boolean useParallel = numThreads > 1 && numSim > 1;

		ExecutorService exec = null;
		int inExec = 0;

		Runnable_ContactMapTransmission[] runnable = new Runnable_ContactMapTransmission[numSim];

		for (int s = 0; s < numSim; s++) {
			long simSeed = rngBase.nextLong();

			runnable[s] = new Runnable_ContactMapTransmission(simSeed, pop_composition, 
					baseContactMap, numSnap * snapFreq);

			runnable[s].initialse();
			
			//TODO: Add infected

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
}
