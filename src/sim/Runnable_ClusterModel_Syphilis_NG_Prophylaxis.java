package sim;

import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;

import relationship.ContactMap;

public class Runnable_ClusterModel_Syphilis_NG_Prophylaxis extends Runnable_ClusterModel_MultiTransmission {

	float prophylaxis_uptake_per_treatment;
	int prophylaxis_dosage_max;
	int prophylaxis_starts_at;
	int prophylaxis_duration = 3;

	private static final int num_inf = 2;
	private static final int num_site = 4;
	private static final int num_act = 3;

	public static final Pattern PROP_TYPE_PATTERN = Pattern
			.compile("Syphilis_NG_Prophylaxis_(-?\\d+\\.\\d+)_(\\d+)_(-?\\d+)");

	protected transient HashMap<Integer, int[]> prophylaxis_record; // K=PID
	protected static final int PROPHYLAXIS_REC_LAST_OFFER_AT = 0;
	protected static final int PROPHYLAXIS_REC_LAST_UPTAKE_AT = PROPHYLAXIS_REC_LAST_OFFER_AT + 1;
	protected static final int PROPHYLAXIS_REC_DOSAGE = PROPHYLAXIS_REC_LAST_UPTAKE_AT + 1;
	protected static final int PROPHYLAXIS_REC_PROTECT_UNTIL = PROPHYLAXIS_REC_DOSAGE + 1;
	protected static final int LENGTH_PROPHYLAXIS_REC = PROPHYLAXIS_REC_PROTECT_UNTIL + 1;

	public Runnable_ClusterModel_Syphilis_NG_Prophylaxis(long cMap_seed, long sim_seed, int[] pop_composition,
			ContactMap base_cMap, int numTimeStepsPerSnap, int numSnap, float prophylaxis_uptake_per_treatment,
			int prophylaxis_dosage, int prophylaxis_starts_at) {
		super(cMap_seed, sim_seed, pop_composition, base_cMap, numTimeStepsPerSnap, numSnap, num_inf, num_site,
				num_act);

		this.prophylaxis_uptake_per_treatment = prophylaxis_uptake_per_treatment;
		this.prophylaxis_dosage_max = prophylaxis_dosage;
		this.prophylaxis_starts_at = prophylaxis_starts_at;

	}

	@Override
	public void initialse() {
		super.initialse();
		prophylaxis_record = new HashMap<>();
	}

	@Override
	protected void applyTreatment(int currentTime, int infId, int pid, int[][] inf_stage) {

		if (prophylaxis_starts_at > 0 && currentTime >= prophylaxis_starts_at) {
			int[] prop_rec = prophylaxis_record.get(pid);
			if (prop_rec == null) {
				prop_rec = new int[LENGTH_PROPHYLAXIS_REC];
				Arrays.fill(prop_rec, -1);
				prophylaxis_record.put(pid, prop_rec);
			}
			if (RNG.nextFloat() < prophylaxis_uptake_per_treatment) {
				if (prop_rec[PROPHYLAXIS_REC_LAST_OFFER_AT] < currentTime) {
					prop_rec[PROPHYLAXIS_REC_DOSAGE] = prophylaxis_dosage_max;
					prop_rec[PROPHYLAXIS_REC_LAST_UPTAKE_AT] = currentTime;

				}
			}
			prop_rec[PROPHYLAXIS_REC_LAST_OFFER_AT] = currentTime;
		}

		super.applyTreatment(currentTime, infId, pid, inf_stage);

	}

	@Override
	protected void simulate_non_infectious_act(int currentTime, ContactMap cMap, HashMap<String, int[]> acted_today) {
		for (Integer pid_inf : prophylaxis_record.keySet()) {
			int[] prop_rec = prophylaxis_record.get(pid_inf);

			if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] < currentTime && prop_rec[PROPHYLAXIS_REC_DOSAGE] > 0) {
				if (cMap.containsVertex(pid_inf)) {
					int g_s = getGenderType(pid_inf);
					Integer[][] edges = cMap.edgesOf(pid_inf).toArray(new Integer[0][]);
					for (int i = 0; i < edges.length && prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] < currentTime; i++) {
						Integer[] e = edges[i];
						int pid_inf_tar = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1].equals(pid_inf)
								? e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]
								: e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];

						int[] partners = new int[] { pid_inf, pid_inf_tar };
						int g_t = getGenderType(pid_inf_tar);
						Arrays.sort(partners);

						int[] hasActed = acted_today.get(Arrays.toString(partners));
						if (hasActed == null) {
							boolean acted = false;
							for (int a = 0; a < NUM_ACT && !acted; a++) {
								double[] fieldEntry = table_act_frequency[a][g_s][g_t];
								if (RNG.nextDouble() < fieldEntry[FIELD_ACT_FREQ_ACT_PER_DAY]) {
									acted = true;
								}
							}
							if (acted && prop_rec[PROPHYLAXIS_REC_DOSAGE] > 0) {
								prop_rec[PROPHYLAXIS_REC_DOSAGE]--;
								prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] = currentTime + prophylaxis_duration;
							}

						}

					}
				}
			}

		}
	}

	@Override
	protected double getTransmissionProb(int currentTime, int inf_id, int pid_inf_src, int pid_inf_tar,
			int partnershipDuration, int actType, int src_site, int tar_site) {
		double transProb = super.getTransmissionProb(currentTime, inf_id, pid_inf_src, pid_inf_tar, partnershipDuration,
				actType, src_site, tar_site);
		for (int pid : new int[] { pid_inf_src, pid_inf_tar }) {
			int[] prop_rec = prophylaxis_record.get(pid);
			if (prop_rec != null) {
				if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] < currentTime && prop_rec[PROPHYLAXIS_REC_DOSAGE] > 0) {
					prop_rec[PROPHYLAXIS_REC_DOSAGE]--;
					prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] = currentTime + prophylaxis_duration;
				}
				if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] >= currentTime && pid == pid_inf_tar) { // Doesn't prevent
																									// transmission
					transProb *= 0;
				}
			}
		}
		return transProb;
	}

}
