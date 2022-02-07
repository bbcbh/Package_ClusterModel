package sim;

import java.io.File;
import java.io.IOException;
import java.util.Properties;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;

import util.PersonClassifier;
import util.PropValUtils;

public class Simulation_ClusterModel implements SimulationInterface {

	public static final Object[] DEFAULT_BRIDGING_SIM_FIELDS = {
			// FIELD_MAX_CASUAL_PARTNER_NUM_PROB
			// float[GENDER]{CUMUL_PROB_1, CUMUL_PROB_2, ... , MAX_PARTNER_1 , MAX_PARTNER_2
			// ...}}
			// Default:
			// National Centre in HIV Epidemiology and Clinical Research. Phase A of the
			// National Gay Men’s Syphilis Action Plan:
			// Modelling evidence and research on acceptability of interventions for
			// controlling syphilis in Australia.
			// Sydney: National Centre in HIV Epidemiology and Clinical Research, 2009.
			// Fogarty A, Mao L, Zablotska Manos I, et al.
			// The Health in Men and Positive Health cohorts: A comparison of trends in the
			// health and sexual behaviour of
			// HIV-negative and HIV-positive gay men, 2002-2005. Sydney: National Centre in
			// HIV Social Research, 2006.d
			new float[][] { { 1, 0 }, { 1, 0 }, { 0.51f, 10 }, { 0.51f, 10 } }, };

	public static final int FIELD_MAX_CASUAL_PARTNER_NUM_PROB = Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD;
	public static final int LENGTH_SIM_CLUSTER_MODEL_FIELD = FIELD_MAX_CASUAL_PARTNER_NUM_PROB + 1;

	public static final String POP_PROP_INIT_PREFIX = "POP_PROP_INIT_PREFIX_";
	public static final String POP_PROP_INIT_PREFIX_CLASS = "POP_PROP_INIT_PREFIX_CLASS_";

	public Object[] simFields = new Object[LENGTH_SIM_CLUSTER_MODEL_FIELD];
	public Class<?>[] simFieldClass = new Class[LENGTH_SIM_CLUSTER_MODEL_FIELD];
	
	
	protected File baseDir = null;
	protected boolean stopNextTurn = false; 

	public Simulation_ClusterModel() {
		for (int i = 0; i < simFields.length; i++) {
			// All simulation levels
			if (i > Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD) {
				simFields[i] = DEFAULT_BRIDGING_SIM_FIELDS[Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD
						+ i];
				simFieldClass[i] = DEFAULT_BRIDGING_SIM_FIELDS[Runnable_ClusterModel.LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD
						+ i].getClass();
			}
		}
	}

	public static PoissonDistribution PoissonDistributionFit(float targetCDF, int targetVal) {
		double init_lamda = targetVal;
		final double RELATIVE_TOLERANCE = 0.005;
		final double ABSOLUTE_TOLERANCE = 0.001;

		final BrentOptimizer OPTIMIZER = new BrentOptimizer(RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);
		final SearchInterval interval = new SearchInterval(1, 100, init_lamda);
		final UnivariateObjectiveFunction func = new UnivariateObjectiveFunction(new UnivariateFunction() {
			@Override
			public double value(double x) {
				PoissonDistribution func_dist = new PoissonDistribution(x);
				return func_dist.cumulativeProbability(targetVal) - targetCDF;
			}
		});

		double bestFit = OPTIMIZER.optimize(func, GoalType.MINIMIZE, interval).getPoint();
		PoissonDistribution dist = new PoissonDistribution(bestFit);
		return dist;
	}

	@Override
	public void loadProperties(Properties prop) {
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
		Properties prop = new Properties();
		for (int i = 0; i < LENGTH_SIM_CLUSTER_MODEL_FIELD; i++) {
			String propName = String.format("%s%d", POP_PROP_INIT_PREFIX, i);
			String className = String.format("%s%d", POP_PROP_INIT_PREFIX_CLASS, i);
			if (simFields[i] != null) {
				prop.put(propName, PropValUtils.objectToPropStr(simFields[i], simFields[i].getClass()));
				prop.put(className, simFields[i].getClass().toString());
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
		// TODO Auto-generated method stub

	}

	@Override
	public void generateOneResultSet() throws IOException, InterruptedException {
		// TODO Auto-generated method stub

	}

}
