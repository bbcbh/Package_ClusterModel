package sim;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Properties;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PropValUtils;

public class Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis
		extends Runnable_ClusterModel_MultiTransmission {

	protected int prophylaxis_duration_per_dose = 3;
	protected float[] prophylaxis_efficacy; // Assume 100 efficiency atm
	protected int prophylaxis_starts_at = -1;

	protected transient HashMap<Integer, int[]> prophylaxis_record;
	protected transient HashMap<String, Float> prophylaxis_resist_profile; // key="INF,ID,SITE"
	protected RandomGenerator rng_PEP;

	protected static final int PROPHYLAXIS_REC_LAST_OFFER_AT = 0;
	protected static final int PROPHYLAXIS_REC_LAST_USE_AT = PROPHYLAXIS_REC_LAST_OFFER_AT + 1;
	protected static final int PROPHYLAXIS_REC_DOSAGE = PROPHYLAXIS_REC_LAST_USE_AT + 1;
	protected static final int PROPHYLAXIS_REC_PROTECT_UNTIL = PROPHYLAXIS_REC_DOSAGE + 1;
	protected static final int LENGTH_PROPHYLAXIS_REC = PROPHYLAXIS_REC_PROTECT_UNTIL + 1;

	public static final String PROP_PEP_START_AT = "PROP_PEP_START_AT";
	public static final String PROP_PEP_DURATION_PER_DOSE = "PROP_PEP_DURATION_PER_DOSE";
	public static final String PROP_PEP_EFFICACY = "PROP_PEP_EFFICACY";

	public static final int SIM_SETTING_KEY_GEN_PEP_USGAGE_RECORD = Simulation_ClusterModelTransmission.SIM_SETTING_KEY_TREATMENT_ON_INFECTIOUS_ONLY
			+ 1;

	public Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis(long cMap_seed, long sim_seed,
			ContactMap base_cMap, Properties prop, int num_inf, int num_site, int num_act) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);

		if (prop.containsKey(PROP_PEP_DURATION_PER_DOSE)) {
			prophylaxis_duration_per_dose = Integer.parseInt(prop.getProperty(PROP_PEP_DURATION_PER_DOSE));
		}
		if (prop.containsKey(PROP_PEP_EFFICACY)) {
			prophylaxis_efficacy = (float[]) PropValUtils.propStrToObject(prop.getProperty(PROP_PEP_EFFICACY),
					float[].class);
		}
		rng_PEP = new MersenneTwisterRandomGenerator(sim_seed);

	}

	public Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis(long cMap_seed, long sim_seed,
			int[] pop_composition, ContactMap base_cMap, int numTimeStepsPerSnap, int numSnap, int num_inf,
			int num_site, int num_act) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap, num_inf, num_site,
				num_act);
	}

	public float[][] generateCurrentPopulationInfectiousResistProfile() {
		ArrayList<float[]> res = new ArrayList<float[]>();		
		
		for (Entry<String, Float> valSet : prophylaxis_resist_profile.entrySet()) {
			String[] profileEnt = valSet.getKey().split(","); // src = 
															  //  String.format("%d,%d,%d", infection_id, pid,site_id);
			int inf = Integer.parseInt(profileEnt[0]);
			int pid = Integer.parseInt(profileEnt[1]);
			int site = Integer.parseInt(profileEnt[2]);
			int[][] inf_stage_all = map_currrent_infection_stage.get(pid);

			if (inf_stage_all == null) {
				System.err.printf("Warning: infection stage for #%d not found "
						+ "despite having prophylaxis_resist_profile of %s\n", pid, valSet.toString());
			} else {
				float[] ent = new float[] { inf, site, pid, inf_stage_all[inf][site], valSet.getValue()};
				int pt = Collections.binarySearch(res, ent, new Comparator<float[]>() {
					@Override
					public int compare(float[] o1, float[] o2) {
						int res = 0;
						int pt = 0;						
						while(res == 0 && pt < o1.length) {
							res = Float.compare(o1[pt], o2[pt]);
							pt++;
						}											
						return res;
					}
				});				
				res.add(~pt, ent);
			}
		}
		return res.toArray(new float[0][]);
	}

	@Override
	public void initialse() {
		super.initialse();
		prophylaxis_record = new HashMap<>();
		prophylaxis_resist_profile = new HashMap<>();
	}

	@Override
	protected int updateInfectionStage(Integer pid, int infection_id, int site_id, int current_infection_stage,
			int current_time, int[][] current_stage_arr, int[][] infection_state_switch, int state_duration_preset) {
		int res = super.updateInfectionStage(pid, infection_id, site_id, current_infection_stage, current_time,
				current_stage_arr, infection_state_switch, state_duration_preset);
		if (!prophylaxis_resist_profile.isEmpty()) {
			// Update doxy-pep resistence profile if needed
			if (infection_state_switch[infection_id][site_id] == AbstractIndividualInterface.INFECT_S) {
				String key_resist_profile = String.format("%d,%d,%d", infection_id, pid, site_id);
				prophylaxis_resist_profile.remove(key_resist_profile);
			}
		}
		return res;
	}

	@Override
	protected void simulate_transmission_success_act(int currentTime, int inf_id, Integer pid_inf_src, int pid_inf_tar,
			int src_site, int tar_site) {
		super.simulate_transmission_success_act(currentTime, inf_id, pid_inf_src, pid_inf_tar, src_site, tar_site);
		if (!prophylaxis_resist_profile.isEmpty()) {
			String key_resist_profile_src = String.format("%d,%d,%d", inf_id, pid_inf_src, src_site);
			Float val = prophylaxis_resist_profile.get(key_resist_profile_src);
			if (val != null) {
				String key_resist_profile_tar = String.format("%d,%d,%d", inf_id, pid_inf_tar, tar_site);
				prophylaxis_resist_profile.put(key_resist_profile_tar, val.floatValue());
			}
		}
	}

	@Override
	protected double getTransmissionProb(int currentTime, int inf_id, int pid_inf_src, int pid_inf_tar,
			int partnershipDuration, int actType, int src_site, int tar_site) {
		double transProb = super.getTransmissionProb(currentTime, inf_id, pid_inf_src, pid_inf_tar, partnershipDuration,
				actType, src_site, tar_site);
		for (int pid : new int[] { pid_inf_tar }) { // PREP only effect susceptibility not transmission
			int[] prop_rec = prophylaxis_record.get(pid);
			if (prop_rec != null) {
				float prophylaxis_efficacy_adjusted = prophylaxis_efficacy.length == NUM_INF
						? prophylaxis_efficacy[inf_id]
						: prophylaxis_efficacy[inf_id * NUM_SITE + tar_site];

				// Implementation of negative prophylaxis_efficacy_adjusted
				if (prophylaxis_efficacy_adjusted < 0) {
					String key_resist_profile = String.format("%d,%d,%d", inf_id, pid_inf_src, src_site);
					Float val = prophylaxis_resist_profile.get(key_resist_profile);
					if (val == null) {
						val = rng_PEP.nextFloat() < -prophylaxis_efficacy_adjusted ? 1f : 0f;
						prophylaxis_resist_profile.put(key_resist_profile, val);
					}
					prophylaxis_efficacy_adjusted *= val.floatValue();
				}

				if (prop_rec[PROPHYLAXIS_REC_LAST_USE_AT] == currentTime) {
					if (prop_rec[PROPHYLAXIS_REC_DOSAGE] != Integer.MAX_VALUE) {
						// Limited prophylaxis dosage
						if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] < currentTime
								&& prop_rec[PROPHYLAXIS_REC_DOSAGE] > 0) {
							prop_rec[PROPHYLAXIS_REC_DOSAGE]--;
							prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] = currentTime + prophylaxis_duration_per_dose;
							transProb *= 1 - prophylaxis_efficacy_adjusted;
						}

					} else {
						// Unlimited dosage
						if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] == Integer.MAX_VALUE
								|| prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] >= currentTime) {
							if (prophylaxis_efficacy != null) {
								transProb *= 1 - prophylaxis_efficacy_adjusted;
							}
						}

					}
				}

			}
		}
		return transProb;
	}

}
