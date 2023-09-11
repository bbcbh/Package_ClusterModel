package optimisation;

import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;

public final class OptFittingFunctionWrapper extends MultivariateFunctionMappingAdapter {

	private OptFittingFunction bounded;

	public OptFittingFunctionWrapper(OptFittingFunction bounded, double[] lower, double[] upper) {
		super(bounded, lower, upper);
		this.bounded = bounded;
	}	

	public OptFittingFunction getBoundedFunc() {
		return bounded;
	}

}