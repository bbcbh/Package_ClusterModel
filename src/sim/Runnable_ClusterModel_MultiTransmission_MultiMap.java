package sim;

import java.util.Properties;
import java.util.regex.Pattern;

/**
 * Not used atm
*/

public class Runnable_ClusterModel_MultiTransmission_MultiMap extends Runnable_ClusterModel_MultiTransmission {
	
	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("MultiTransmission_MultiMap_(\\d+)_(\\d+)_(\\d+)");
	
	protected String[][] cmaps_str;
	protected int[] cmaps_pt;

	public Runnable_ClusterModel_MultiTransmission_MultiMap(
			long cMap_seed, long sim_seed, 
			String[][] cmaps_str,
			Properties prop, 
			int num_inf, int num_site, int num_act) {
		super(cMap_seed, sim_seed, null, prop, num_inf, num_site, num_act);
		this.cmaps_str = cmaps_str;
		this.cmaps_pt = new int[this.cmaps_str.length];
	}

	
}
