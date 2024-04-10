package sim;

import java.util.Properties;
import java.util.regex.Pattern;

import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_Prophylaxis extends Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis {
	
	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("ClusterModel_Prophylaxis");

	private static final int num_inf = 3; // TP, NG and CT
	private static final int num_site = 4;
	private static final int num_act = 3;
	
	protected float prophylaxis_adherence;
	protected float prophylaxis_uptake_HIV_PrEP;
	protected float prophylaxis_uptake_last_TP;
	protected float prophylaxis_uptake_last_STI;	
	protected RandomGenerator rng_PEP;
	
	public static final String PROP_PEP_ADHERENCE = "PROP_PEP_ADHERENCE";
	public static final String PROP_PEP_UPTAKE = "PROP_PEP_UPTAKE";	

	public Runnable_ClusterModel_Prophylaxis(long cMap_seed, long sim_seed, ContactMap base_cMap, Properties prop) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);
		this.prophylaxis_starts_at = Integer.parseInt(prop.getProperty(PROP_PEP_START_AT, "-1"));
		this.prophylaxis_adherence = Float.parseFloat(prop.getProperty(PROP_PEP_ADHERENCE,"0"));
		float[] update_rate = (float[]) PropValUtils.propStrToObject(prop.getProperty(PROP_PEP_UPTAKE,"[0.0,0.0,0.0]"), float[].class);
		this.prophylaxis_uptake_HIV_PrEP = update_rate[0];
		this.prophylaxis_uptake_last_TP = update_rate[1];
		this.prophylaxis_uptake_last_STI = update_rate[2];
		
		rng_PEP = new MersenneTwisterRandomGenerator(sim_seed);
	}


	@Override
	protected void postTimeStep(int currentTime) {		
		super.postTimeStep(currentTime);
		if(currentTime == this.prophylaxis_starts_at) {
			// TODO: Allocate PEP base on rate
			
			
		}
		
	}
	
	
	
	
	
	
	
}
