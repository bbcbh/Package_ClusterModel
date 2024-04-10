package sim;

import java.util.HashMap;
import java.util.Properties;

import relationship.ContactMap;

public class Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis
		extends Runnable_ClusterModel_MultiTransmission {

	protected int prophylaxis_duration = 3;
	protected float prophylaxis_efficiency = 1; // Assume 100 efficiency atm
	protected int prophylaxis_starts_at = -1;

	protected transient HashMap<Integer, int[]> prophylaxis_record;

	protected static final int PROPHYLAXIS_REC_LAST_OFFER_AT = 0;
	protected static final int PROPHYLAXIS_REC_LAST_UPTAKE_AT = PROPHYLAXIS_REC_LAST_OFFER_AT + 1;
	protected static final int PROPHYLAXIS_REC_DOSAGE = PROPHYLAXIS_REC_LAST_UPTAKE_AT + 1;
	protected static final int PROPHYLAXIS_REC_PROTECT_UNTIL = PROPHYLAXIS_REC_DOSAGE + 1;
	protected static final int LENGTH_PROPHYLAXIS_REC = PROPHYLAXIS_REC_PROTECT_UNTIL + 1;
	
	public static final String PROP_PEP_START_AT = "PROP_PEP_START_AT";
	public static final String PROP_PEP_DURATION_PER_DOSE = "PROP_PEP_DURATION_PER_DOSE";
	public static final String PROP_PEP_EFFICIENCY = "PROP_PEP_EFFICIENCY";
	
	public Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis(long cMap_seed, long sim_seed,
			ContactMap base_cMap, Properties prop, int num_inf,	int num_site, int num_act) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);
		
		if(prop.containsKey(PROP_PEP_DURATION_PER_DOSE)) {
			prophylaxis_duration = Integer.parseInt(prop.getProperty(PROP_PEP_DURATION_PER_DOSE));			
		}
		if(prop.containsKey(PROP_PEP_EFFICIENCY)) {
			prophylaxis_efficiency = Float.parseFloat(prop.getProperty(PROP_PEP_EFFICIENCY));
		}			
		
	}

	public Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis(long cMap_seed, long sim_seed,
			int[] pop_composition, ContactMap base_cMap, int numTimeStepsPerSnap, int numSnap, int num_inf,
			int num_site, int num_act) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap, num_inf, num_site,
				num_act);
	}

	@Override
	public void initialse() {
		super.initialse();
		prophylaxis_record = new HashMap<>();
	}

	@Override
	protected double getTransmissionProb(int currentTime, int inf_id, int pid_inf_src, int pid_inf_tar,
			int partnershipDuration, int actType, int src_site, int tar_site) {
		double transProb = super.getTransmissionProb(currentTime, inf_id, pid_inf_src, pid_inf_tar, partnershipDuration,
				actType, src_site, tar_site);
		for (int pid : new int[] { pid_inf_tar }) { // PREP only effect susceptibility not transmission
			int[] prop_rec = prophylaxis_record.get(pid);
			if (prop_rec != null) {				
				// Limited prophylaxis dosage options 
				if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] < currentTime && prop_rec[PROPHYLAXIS_REC_DOSAGE] > 0
						&& prop_rec[PROPHYLAXIS_REC_DOSAGE] != Integer.MAX_VALUE) {
					prop_rec[PROPHYLAXIS_REC_DOSAGE]--;
					prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] = currentTime + prophylaxis_duration;
				}
				if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] == Integer.MAX_VALUE 
						|| prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] >= currentTime) {
					transProb *= 1-prophylaxis_efficiency; 
				}
			}
		}
		return transProb;
	}

}
