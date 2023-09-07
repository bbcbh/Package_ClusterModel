package optimisation;

import org.apache.commons.math3.analysis.MultivariateFunction;

import sim.Runnable_ClusterModel_Transmission;

public abstract class OptFittingFunction implements MultivariateFunction {

	private Runnable_ClusterModel_Transmission[] runnables;
	private double[] bestResidue_by_runnable;		

	public Runnable_ClusterModel_Transmission[] getRunnables() {
		return runnables;
	}

	public double[] getBestResidue_by_runnable() {
		return bestResidue_by_runnable;
	}

}