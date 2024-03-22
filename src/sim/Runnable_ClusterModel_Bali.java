package sim;

import java.util.regex.Pattern;

import relationship.ContactMap;

public class Runnable_ClusterModel_Bali extends Runnable_ClusterModel_MultiTransmission {
	
	// 0: TP, 1: NG, 2: CT, 3: HIV
	private static final int num_inf = 4; 
	// 0: ANY, 1: URETHAL, 2: RECTAL
	private static final int num_site = 3;
	 //0 = ANY, 1 = ANAL
	private static final int num_act = 2;
	
	
	public static final Pattern PROP_TYPE_PATTERN = Pattern
			.compile("Bali_Model");
	

	public Runnable_ClusterModel_Bali(long cMap_seed, long sim_seed, int[] pop_composition, ContactMap base_cMap,
			int numTimeStepsPerSnap, int numSnap) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap, num_inf, num_site, num_act);		
	}

}
