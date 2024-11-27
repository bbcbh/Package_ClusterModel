package optimisation;

import java.util.Properties;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.MultivariateFunction;

public abstract class OptFittingFunction implements MultivariateFunction {   
	
	// POP_PROP_OPT_PARAM_FIT_SETTING_DIFF_FORMAT
	// Eg. 18_Diff1 means Param[18] = opt. parameter value + value of opt parameter with index 1
	public static final Pattern POP_PROP_OPT_PARAM_FIT_SETTING_DIFF_FORMAT = Pattern.compile("Diff(\\d+)");	
	
	// POP_PROP_OPT_PARAM_TRANSFORM
	// Format: String[][] {popPropInitPrefix_IncIndices,
	// popPropInitPrefix_IncIndices_src_1, ratio_1 ...}
	// e.g. [A, *B, 0.2, C, 0.1, Const ] => A = A * 0.2 B + 0.1 C + Const
	public static final String POP_PROP_OPT_PARAM_TRANSFORM = "POP_PROP_OPT_PARAM_TRANSFORM";
	
	public abstract long[] getSim_seeds();	
	public abstract long[] getCMap_seeds();	
	public abstract Properties getProperties();	
	public abstract double[] getBestResidue_by_runnable();

}