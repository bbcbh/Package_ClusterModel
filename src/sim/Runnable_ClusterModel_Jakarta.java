package sim;

import java.util.Properties;
import java.util.regex.Pattern;

public class Runnable_ClusterModel_Jakarta extends Runnable_ClusterModel_MultiTransmission {

	// 0: TP, 1: NG, 2: CT, 3: HIV
	private static final int num_inf = 4;
	// 0: ANY, 1: URETHAL, 2: RECTAL
	private static final int num_site = 3;
	// 0 = ANY, 1 = ANAL
	private static final int num_act = 2;
	
	
	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("Jakarta_Model");

	public Runnable_ClusterModel_Jakarta(long cMap_seed, long sim_seed, Properties prop) {
		super(cMap_seed, sim_seed, null, prop, num_inf, num_site, num_act);		
	}
	
	
	
	@Override
	protected void postSimulation() {
		super.postSimulation();
		//TODO: To be added
		
		System.out.printf("%s.postSimulation() to be implemented.\n", this.getClass().getName());
		
		if (print_progress != null && runnableId != null) {
			try {
				print_progress.printf("Post simulation file generation for Thread <%s> completed. Timestamp = %tc.\n",
						runnableId, System.currentTimeMillis());
			} catch (Exception ex) {
				System.err.printf("Post simulation file generation for Thread <%s> completed.\n", runnableId);
			}
		}
	}

}
