package optimisation;

import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;

import sim.Runnable_ClusterModel_Transmission;

public final class OptFittingFunctionWrapper extends MultivariateFunctionMappingAdapter {

	private OptFittingFunction bounded;

	public OptFittingFunctionWrapper(OptFittingFunction bounded, double[] lower, double[] upper) {
		super(bounded, lower, upper);
		this.bounded = bounded;
	}

	public Runnable_ClusterModel_Transmission[] getRunnables() {
		return bounded.getRunnables();
	}

	public double[] getBestResidue_by_runnable() {
		return bounded.getBestResidue_by_runnable();
	}

	public OptFittingFunction getBoundedFunc() {
		return bounded;
	}

}