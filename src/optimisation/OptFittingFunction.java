package optimisation;

import java.util.Properties;

import org.apache.commons.math3.analysis.MultivariateFunction;

public abstract class OptFittingFunction implements MultivariateFunction {   
	
	public abstract long[] getSim_seeds();	
	public abstract long[] getCMap_seeds();	
	public abstract Properties getProperties();	
	public abstract double[] getBestResidue_by_runnable();

}