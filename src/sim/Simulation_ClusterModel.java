package sim;

import java.io.File;
import java.io.IOException;
import java.util.Properties;
import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;

import optimisation.AbstractParameterOptimiser;
import optimisation.AbstractResidualFunc;
import optimisation.LineSearchOptimisier;
import person.AbstractIndividualInterface;
import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import util.PersonClassifier;
import util.PropValUtils;

public class Simulation_ClusterModel implements SimulationInterface {

	public static final Object[] DEFAULT_BRIDGING_SIM_FIELDS = {};
	public static final int LENGTH_SIM_CLUSTER_MODEL_FIELD =  Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD;

	public static final String POP_PROP_INIT_PREFIX = "POP_PROP_INIT_PREFIX_";
	public static final String POP_PROP_INIT_PREFIX_CLASS = "POP_PROP_INIT_PREFIX_CLASS_";

	public Object[] simFields = new Object[LENGTH_SIM_CLUSTER_MODEL_FIELD];
	public Class<?>[] simFieldClass = new Class[LENGTH_SIM_CLUSTER_MODEL_FIELD];

	protected File baseDir = null;
	protected boolean stopNextTurn = false;
	protected Properties loadedProperties = null; // From .prop file, if any

	public Simulation_ClusterModel() {
		for (int i = 0; i < simFields.length; i++) {
			// All simulation levels
			if (i >= Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD) {
				simFields[i] = DEFAULT_BRIDGING_SIM_FIELDS[i
						- Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD];
				simFieldClass[i] = DEFAULT_BRIDGING_SIM_FIELDS[i
						- Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD].getClass();
			}
		}
	}


	@Override
	public void loadProperties(Properties prop) {

		loadedProperties = prop;

		for (int i = 0; i < LENGTH_SIM_CLUSTER_MODEL_FIELD; i++) {
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
		for (int i = 0; i < LENGTH_SIM_CLUSTER_MODEL_FIELD; i++) {
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

		Runnable_ClusterModel[] runnable = new Runnable_ClusterModel[numSim];

		if (useParallel) {
			executor = Executors.newFixedThreadPool(numThreads);
		}

		for (int i = 0; i < numSim; i++) {
			Runnable_ClusterModel r = new Runnable_ClusterModel();
			runnable[i] = r;

			Population_Bridging population = new Population_Bridging(rngBase.nextLong());
			r.setPopulation(population);
			r.setNumSnaps(numSnap);
			r.setSnapFreq(snapFreq);

			if (!useParallel) {
				runnable[i].run();
			} else {
				executor.submit(runnable[i]);
				numInPool++;
				if (numInPool == numThreads) {
					executor.shutdown();
					if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
						showStrStatus("Thread time-out!");
					}
					numInPool = 0;
				}
			}
		}
		if (useParallel && numInPool != 0) {
			executor.shutdown();
			if (!executor.awaitTermination(2, TimeUnit.DAYS)) {
				showStrStatus("Thread time-out!");
			}
		}

		showStrStatus(String.format("Time required = %.3f s", (System.currentTimeMillis() - tic) / 1000f));

	}

	private void showStrStatus(String string) {
		System.out.println(string);

	}

}
