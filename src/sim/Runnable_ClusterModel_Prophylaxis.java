package sim;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.RealDistribution;

import person.AbstractIndividualInterface;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_Prophylaxis extends Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis {

	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("ClusterModel_Prophylaxis");

	private static final int num_inf = 3; // TP, NG and CT
	private static final int num_site = 4;
	private static final int num_act = 5;
	private RealDistribution persistenceDist;

	private static final int UPTAKE_HIV_PrEP = 0;
	private static final int UPTAKE_DX_TP = UPTAKE_HIV_PrEP + 1;
	private static final int UPTAKE_DX_STI = UPTAKE_DX_TP + 1;
	private static final int UPTAKE_PARTNERS = UPTAKE_DX_STI + 1;
	private static final int UPTAKE_PARTNERS_LIMIT = UPTAKE_PARTNERS + 1;
	private static final int UPTAKE_RATE_INDIVDUAL_ADJ_UPTAKE = UPTAKE_PARTNERS_LIMIT + 1;
	private static final int UPTAKE_RATE_INDIVDUAL_ADJ_NON_UPTAKE = UPTAKE_RATE_INDIVDUAL_ADJ_UPTAKE + 1;

	protected double[] prophylaxis_persistence_adherence;
	protected float[] prophylaxis_uptake;
	protected HashMap<Integer, int[]> dx_last_12_months;
	protected HashMap<Integer, int[][]> sexual_contact_last_12_months;

	protected HashMap<Integer, Float> pep_uptake_individual_rate_adj;
	protected HashMap<Integer, ArrayList<Integer>> pep_usage_record;

	public static final String PROP_PEP_PERSISTENCE_ADHERENCE = "PROP_PEP_PERSISTENCE_ADHERENCE";
	public static final String PROP_PEP_UPTAKE = "PROP_PEP_UPTAKE";

	protected static final int PEP_AVAIL_ANY_STI = 1; // More than once (including current)

	// FILENAME_PREVALENCE_PERSON, FILENAME_CUMUL_INCIDENCE_PERSON
	private static final int[] COL_SEL_INF_GENDER = new int[] { 2, 6, 10 };
	// FILENAME_CUMUL_INCIDENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE = new int[] { 8, 25, 26, 27, 41, 42, 43 };
	// "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_SITE = new int[] { 0, 5, 6, 7, 9, 10, 11 };
	// "Infected_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE_AT = null; // new int[] { 64, 65, 192, 194, 196, 198, 200, 202,
																	// 204 160, 162, 164, 166, 224, 225 };

	private static final String SIM_OUTPUT_KEY_PEP_RESIST_PROFILE = "SIM_OUTPUT_KEY_PEP_RESIST_PROFILE";
	private static final String SIM_OUTPUT_KEY_PEP_COVERAGE = "SIM_OUTPUT_PEP_COVERAGE";
	// ent: Ever PEP, Currently on PEP, # script offered
	protected int count_PEP_OFFERED = 0;
	protected int count_PEP_USED = 0;
	protected int count_PEP_USED_EFFECTIVE = 0;

	public Runnable_ClusterModel_Prophylaxis(long cMap_seed, long sim_seed, ContactMap base_cMap, Properties prop) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);
		this.prophylaxis_starts_at = Integer.parseInt(prop.getProperty(PROP_PEP_START_AT, "-1"));
		this.prophylaxis_persistence_adherence = (double[]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PEP_PERSISTENCE_ADHERENCE, "[0,0,1]"), double[].class);
		prophylaxis_uptake = (float[]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PEP_UPTAKE, "[0.0,0.0,0.0,0.0,0.0,0.0,0.0]"), float[].class);

		dx_last_12_months = new HashMap<>();
		sexual_contact_last_12_months = new HashMap<>();
		pep_uptake_individual_rate_adj = new HashMap<>();
		pep_usage_record = new HashMap<>();

		double[] persistence = Arrays.copyOf(prophylaxis_persistence_adherence, 2);

		if (persistence[1] == 0) {
			persistenceDist = generateNonDistribution(rng_PEP, persistence);
		} else if (persistence[1] < 0) {
			double[] param = Arrays.copyOf(persistence, 2);
			param[1] = Math.abs(param[1]);
			persistenceDist = generateUniformDistribution(rng_PEP, param);
		} else {
			persistenceDist = generateGammaDistribution(rng_PEP, persistence);
		}
	}

	@Override
	protected void simulate_non_infectious_act(int currentTime, ContactMap cMap, HashMap<String, int[]> acted_today) {
		super.simulate_non_infectious_act(currentTime, cMap, acted_today);

		if (this.prophylaxis_starts_at > 0
				&& currentTime > prophylaxis_starts_at - AbstractIndividualInterface.ONE_YEAR_INT) {

			// Update act status for all
			// (only called if PEP is offered based on act history)

			for (Integer pid : cMap.vertexSet()) {
				for (Integer[] edge : cMap.edgesOf(pid)) {
					Integer[] partners = Arrays.copyOf(edge, 2);
					Arrays.sort(partners);
					String key = Arrays.toString(partners);
					int[] hasActed = acted_today.get(key);
					boolean acted = false;
					if (hasActed == null) { // Only count if acts between partnership is not included previously
						hasActed = new int[NUM_ACT];
						for (int a = 0; a < NUM_ACT; a++) {
							double[] fieldEntry = table_act_frequency[a][getGenderType(partners[0])][getGenderType(
									partners[1])];
							if (RNG.nextDouble() < fieldEntry[FIELD_ACT_FREQ_ACT_PER_DAY]) {
								hasActed[a] = 1;
							}
						}
						acted_today.put(key, hasActed);

						for (int a = 0; a < NUM_ACT && !acted; a++) {
							acted |= hasActed[a] != 0;
						}
						if (acted) {
							int recIndec = currentTime % AbstractIndividualInterface.ONE_YEAR_INT;
							int partnerId = partners[0].equals(pid) ? partners[1] : partners[0];
							int[][] contact_hist = sexual_contact_last_12_months.get(pid);
							if (contact_hist == null) {
								contact_hist = new int[AbstractIndividualInterface.ONE_YEAR_INT][];
								sexual_contact_last_12_months.put(pid, contact_hist);
							}
							if (contact_hist[recIndec] == null) {
								contact_hist[recIndec] = new int[] { partnerId };
							} else {
								contact_hist[recIndec] = Arrays.copyOf(contact_hist[recIndec],
										contact_hist[recIndec].length + 1);
								contact_hist[recIndec][contact_hist[recIndec].length - 1] = partnerId;
							}

							for (int pId : partners) {
								int[] prop_rec = prophylaxis_record.get(pId);
								if (prop_rec != null) {
									if (Math.abs(prop_rec[PROPHYLAXIS_REC_LAST_USE_AT]) != currentTime) { // if time <
																											// 0,
																											// then PrEP
																											// is
																											// not used
										double adherence = prophylaxis_persistence_adherence[prophylaxis_persistence_adherence.length
												- 1];
										prop_rec[PROPHYLAXIS_REC_LAST_USE_AT] = (adherence >= 1
												|| rng_PEP.nextDouble() < adherence) ? currentTime : -currentTime;

										if (prop_rec[PROPHYLAXIS_REC_LAST_USE_AT] == currentTime) {
											if ((simSetting
													& 1 << Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis.SIM_SETTING_KEY_GEN_PEP_USGAGE_RECORD) != 0) {
												ArrayList<Integer> ent = pep_usage_record.get(pId);
												if (ent == null) {
													ent = new ArrayList<>();
													pep_usage_record.put(pId, ent);
												}
												ent.add(-currentTime);
											}
											count_PEP_USED++;
										}
									}

								}
							}

						}
					}

				}

			} // End of for (Integer pid : cMap.vertexSet()) {

		}
	}

	@Override
	protected void testPerson(int currentTime, int pid, int infIncl, int siteIncl, int[][] cumul_treatment_by_person) {
		super.testPerson(currentTime, pid, infIncl, siteIncl, cumul_treatment_by_person);

		// Doxy-PEP allocation by test
		if (this.prophylaxis_starts_at > 0 && currentTime >= this.prophylaxis_starts_at) {

			boolean allocatePEP = false;
			boolean offeredPEP = false;

			float individual_uptake_rate_adj = 1;
			float pAlloc;

			if (pep_uptake_individual_rate_adj.containsKey(pid)) {
				individual_uptake_rate_adj = pep_uptake_individual_rate_adj.get(pid);
			}

			// Risk group based PEP
			if (!allocatePEP && prophylaxis_uptake[UPTAKE_HIV_PrEP] > 0) {
				if (!allocatePEP && risk_cat_map.get(pid).intValue() == 0 || risk_cat_map.get(pid).intValue() == 1) {
					pAlloc = individual_uptake_rate_adj * prophylaxis_uptake[UPTAKE_HIV_PrEP];
					offeredPEP |= true;
					allocatePEP |= pAlloc >= 1 || rng_PEP.nextFloat() < pAlloc;
				}
			}

			// Treatment rate based PEP
			if (!allocatePEP && (prophylaxis_uptake[UPTAKE_DX_TP] > 0 || prophylaxis_uptake[UPTAKE_DX_STI] > 0)) {
				int[] dx_hist = dx_last_12_months.get(pid);
				if (dx_hist != null) {
					int current_dx = dx_hist[currentTime % dx_hist.length];
					if (current_dx != 0) { // Has a positive dx today
						// Has TP DX
						if (!allocatePEP && prophylaxis_uptake[UPTAKE_DX_TP] > 0 && ((current_dx & 1) != 0)) {
							pAlloc = individual_uptake_rate_adj * prophylaxis_uptake[UPTAKE_DX_TP];
							offeredPEP |= true;
							allocatePEP |= pAlloc >= 1 || rng_PEP.nextFloat() < pAlloc;
						}
						// Check for STI DX
						if (!allocatePEP && prophylaxis_uptake[UPTAKE_DX_STI] > 0) {
							int num_dx_12_month_any = 0;
							for (int i = 0; i < dx_hist.length && !(num_dx_12_month_any > PEP_AVAIL_ANY_STI); i++) {
								if (dx_hist[i] != 0) {
									num_dx_12_month_any++;
								}
							}
							if (num_dx_12_month_any > PEP_AVAIL_ANY_STI) {
								pAlloc = individual_uptake_rate_adj * prophylaxis_uptake[UPTAKE_DX_STI];
								offeredPEP |= true;
								allocatePEP |= pAlloc >= 1 || rng_PEP.nextFloat() < pAlloc;
							}
						}
					}
				}
			} // End of treatment based PEP

			// Partner based PEP
			if (!allocatePEP && prophylaxis_uptake[UPTAKE_PARTNERS] > 0) {
				int[][] partner_hist = sexual_contact_last_12_months.get(pid);
				if (partner_hist != null) {
					ArrayList<Integer> partnerList = new ArrayList<>();
					for (int i = 0; i < partner_hist.length; i++) {
						if (partner_hist[i] != null) {
							for (Integer part : partner_hist[i]) {
								int pt = Collections.binarySearch(partnerList, part);
								if (pt < 0) {
									partnerList.add(~pt, part);
								}
							}
						}
					}
					if (partnerList.size() > prophylaxis_uptake[UPTAKE_PARTNERS_LIMIT]) {
						pAlloc = individual_uptake_rate_adj * prophylaxis_uptake[UPTAKE_PARTNERS];
						offeredPEP |= true;
						allocatePEP |= pAlloc >= 1 || rng_PEP.nextFloat() < pAlloc;
					}
				}
			}

			if (allocatePEP) {
				allocateProphylaxis(currentTime, pid);
				if (prophylaxis_uptake[UPTAKE_RATE_INDIVDUAL_ADJ_UPTAKE] > 0) {
					pep_uptake_individual_rate_adj.put(pid,
							individual_uptake_rate_adj * prophylaxis_uptake[UPTAKE_RATE_INDIVDUAL_ADJ_UPTAKE]);
				}
			} else {
				if (offeredPEP && prophylaxis_uptake[UPTAKE_RATE_INDIVDUAL_ADJ_NON_UPTAKE] > 0) {
					pep_uptake_individual_rate_adj.put(pid,
							individual_uptake_rate_adj * prophylaxis_uptake[UPTAKE_RATE_INDIVDUAL_ADJ_NON_UPTAKE]);

				}
			}

		}

	}

	@Override
	protected void applyTreatment(int currentTime, int infId, int pid, int[][] inf_stage) {
		super.applyTreatment(currentTime, infId, pid, inf_stage);
		int[] dx_hist = dx_last_12_months.get(pid);
		if (dx_hist == null) {
			dx_hist = new int[AbstractIndividualInterface.ONE_YEAR_INT];
			dx_last_12_months.put(pid, dx_hist);
		}
		dx_hist[currentTime % dx_hist.length] |= 1 << infId;
	}

	@Override
	protected void postTimeStep(int currentTime) {
		super.postTimeStep(currentTime);

		// Treatment history
		for (Integer treated_previously : dx_last_12_months.keySet()) {
			int[] dx_hist = dx_last_12_months.get(treated_previously);
			if (dx_hist != null) {
				// Reset record for the next day
				dx_hist[(currentTime + 1) % dx_hist.length] = 0;
			}
		}

		for (Integer hasPartnership : sexual_contact_last_12_months.keySet()) {
			int[][] partner_hist = sexual_contact_last_12_months.get(hasPartnership);
			if (partner_hist != null) {
				partner_hist[(currentTime + 1) % partner_hist.length] = null;
			}
		}

		if (currentTime % nUM_TIME_STEPS_PER_SNAP == 0 && prophylaxis_record.size() > 0) {
			@SuppressWarnings("unchecked")
			HashMap<Integer, int[]> pep_coverage = (HashMap<Integer, int[]>) sim_output
					.get(SIM_OUTPUT_KEY_PEP_COVERAGE);
			if (pep_coverage == null) {
				pep_coverage = new HashMap<>();
				sim_output.put(SIM_OUTPUT_KEY_PEP_COVERAGE, pep_coverage);
			}
			int[] stat = new int[5]; // Ever PEP, Currently on PEP, # script offered, PEP used (any), PEP used
										// (effective)
			stat[0] = prophylaxis_record.size();
			for (int pid : prophylaxis_record.keySet()) {
				int[] prop_rec = prophylaxis_record.get(pid);
				if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] > currentTime) {
					stat[1]++;
				}

			}
			stat[2] = count_PEP_OFFERED;
			stat[3] = count_PEP_USED;
			stat[4] = count_PEP_USED_EFFECTIVE;
			pep_coverage.put(currentTime, stat);

			if (!prophylaxis_efficacy_adjust.isEmpty()) {
				@SuppressWarnings("unchecked")
				HashMap<Integer, int[]> pep_resist_by_inf_site = (HashMap<Integer, int[]>) sim_output
						.get(SIM_OUTPUT_KEY_PEP_RESIST_PROFILE);
				if (pep_resist_by_inf_site == null) {
					pep_resist_by_inf_site = new HashMap<>();
					sim_output.put(SIM_OUTPUT_KEY_PEP_RESIST_PROFILE, pep_resist_by_inf_site);
				}

				int[] count_pep_resist = new int[3 * NUM_INF * NUM_SITE];

				for (int i = 0; i < NUM_INF; i++) {
					for (int s = 0; s < NUM_SITE; s++) {
						ArrayList<Integer> res = map_currently_infectious.get(String.format("%d,%d", i, s));
						if (res != null) {
							count_pep_resist[i * NUM_SITE + s] = res.size();
							for (Integer inf_pid : res) {															
								// 0 = Resist to DoxyPEP, 1 = Sensitive to DoxyPEP
								Integer pep_efficiency_adj = prophylaxis_efficacy_adjust
										.get(String.format("%d,%d,%d", i, inf_pid, s));
								if (pep_efficiency_adj != null) {
									int offset = NUM_INF * NUM_SITE;
									if (pep_efficiency_adj < 1) {
										offset *= 2;
									}
									count_pep_resist[offset + i * NUM_SITE + s]++;
								}
							}
						}

					}
				}

				pep_resist_by_inf_site.put(currentTime, count_pep_resist);
			}

		}

	}

	protected void allocateProphylaxis(int currentTime, Integer pid) {
		int[] prop_rec = prophylaxis_record.get(pid);
		if (prop_rec == null) {
			prop_rec = new int[LENGTH_PROPHYLAXIS_REC];
			Arrays.fill(prop_rec, -1);
			prophylaxis_record.put(pid, prop_rec);
		}
		prop_rec[PROPHYLAXIS_REC_LAST_OFFER_AT] = currentTime;
		prop_rec[PROPHYLAXIS_REC_LAST_USE_AT] = 0;
		prop_rec[PROPHYLAXIS_REC_DOSAGE] = Integer.MAX_VALUE; // Infinite in this model
		prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] = currentTime + (int) Math.round(persistenceDist.sample());

		count_PEP_OFFERED++;
	}

	@Override
	protected double getTransmissionProb(int currentTime, int inf_id, int pid_inf_src, int pid_inf_tar,
			int partnershipDuration, int actType, int src_site, int tar_site) {
		int[] prop_rec = prophylaxis_record.get(pid_inf_tar);
		if (prop_rec != null) {
			if (Math.abs(prop_rec[PROPHYLAXIS_REC_LAST_USE_AT]) != currentTime) { // if time < 0, then PrEP is not used
				double adherence = prophylaxis_persistence_adherence[prophylaxis_persistence_adherence.length - 1];
				prop_rec[PROPHYLAXIS_REC_LAST_USE_AT] = (adherence >= 1 || rng_PEP.nextDouble() < adherence)
						? currentTime
						: -currentTime;

				if (prop_rec[PROPHYLAXIS_REC_LAST_USE_AT] == currentTime) {
					if ((simSetting
							& 1 << Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis.SIM_SETTING_KEY_GEN_PEP_USGAGE_RECORD) != 0) {
						ArrayList<Integer> ent = pep_usage_record.get(pid_inf_tar);
						if (ent == null) {
							ent = new ArrayList<>();
							pep_usage_record.put(pid_inf_tar, ent);
						}
						ent.add(currentTime);
					}

					count_PEP_USED_EFFECTIVE++;
					count_PEP_USED++;
				}
			}
		}

		return super.getTransmissionProb(currentTime, inf_id, pid_inf_src, pid_inf_tar, partnershipDuration, actType,
				src_site, tar_site);

	}

	@SuppressWarnings("unchecked")
	@Override
	protected void postSimulation() {
		super.postSimulation();
		String key, fileName;
		HashMap<Integer, int[]> countMap;
		String filePrefix = this.getRunnableId() == null ? "" : this.getRunnableId();

		// PEP Usage
		countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_PEP_COVERAGE);
		if (countMap != null) {
			fileName = String.format(filePrefix + "PEP_Stat_%d_%d.csv", cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "PEP_User_Type_%d", new int[] { 5 });
		}

		// PEP resist
		countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_PEP_RESIST_PROFILE);
		if (countMap != null) {
			fileName = String.format(filePrefix + "PEP_Resist_Profile_%d_%d.csv", cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "PEP_ResistId_%d_Inf_%d_Site_%d", new int[] { 3, NUM_INF, NUM_SITE });
			// NG only: new int[] {5,6,7,17,18,19,29,30,31}
		}

		if ((simSetting
				& 1 << Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis.SIM_SETTING_KEY_GEN_PEP_USGAGE_RECORD) != 0) {
			if (!pep_usage_record.isEmpty()) {
				fileName = String.format(filePrefix + "PEP_Usage_Record_%d_%d.csv", cMAP_SEED, sIM_SEED);
				try {
					PrintWriter pWri = new PrintWriter(new java.io.File(baseDir, fileName));
					pWri.println("PID,TIME_PEP_TAKEN");
					for (Entry<Integer, ArrayList<Integer>> valSet : pep_usage_record.entrySet()) {
						pWri.print(valSet.getKey());
						for (Integer timeStep : valSet.getValue()) {
							pWri.print(',');
							pWri.print(timeStep);
						}
						pWri.println();
					}
					pWri.close();
				} catch (IOException e) {
					e.printStackTrace(System.err);
				}
			}
		}

		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE) != 0) {
			key = String.format(SIM_OUTPUT_KEY_CUMUL_TREATMENT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_TREATMENT_FILE);

			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_TREATMENT_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);

		}
		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE) != 0) {

			key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);

			key = String.format(SIM_OUTPUT_KEY_CUMUL_INCIDENCE_SITE,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_INCIDENCE_FILE);

			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(filePrefix + Simulation_ClusterModelTransmission.FILENAME_CUMUL_INCIDENCE_SITE,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d_Site_%d", new int[] { NUM_INF, NUM_GENDER, NUM_SITE },
					COL_SEL_INF_GENDER_SITE);

		}

		if ((simSetting & 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE) != 0) {

			key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_GENDER_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(
					filePrefix + "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_PERSON,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d", new int[] { NUM_INF, NUM_GENDER },
					COL_SEL_INF_GENDER);

			key = String.format(SIM_OUTPUT_KEY_INFECTIOUS_SITE_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(
					filePrefix + "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE,
					cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Site_%d", new int[] { NUM_INF, NUM_SITE }, COL_SEL_INF_SITE);

			key = String.format(SIM_OUTPUT_KEY_INFECTED_AT_GENDER_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			countMap = (HashMap<Integer, int[]>) sim_output.get(key);
			fileName = String.format(
					filePrefix + "Infected_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE, cMAP_SEED,
					sIM_SEED);
			printCountMap(countMap, fileName, "Inf_%d_Gender_%d_Infected_SiteInc_%d",
					new int[] { NUM_INF, NUM_GENDER, 1 << (NUM_SITE + 1) }, COL_SEL_INF_GENDER_SITE_AT);

			key = String.format(SIM_OUTPUT_KEY_INFECTED_SITE_STAGE_COUNT,
					Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_PREVAL_FILE);
			HashMap<Integer, HashMap<String, Integer>> infected_site_stage_count = (HashMap<Integer, HashMap<String, Integer>>) sim_output
					.get(key);

			fileName = String.format(
					filePrefix + "Infected_All_Stages_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE,
					cMAP_SEED, sIM_SEED);
			try {
				PrintWriter pWri = new PrintWriter(new FileWriter(new java.io.File(baseDir, fileName)));
				Integer[] timeArr = infected_site_stage_count.keySet().toArray(new Integer[0]);
				Arrays.sort(timeArr);
				ArrayList<String[]> col_names = new ArrayList<>();

				for (Integer time : timeArr) {
					HashMap<String, Integer> ent = infected_site_stage_count.get(time);
					for (String col_name : ent.keySet()) {
						String[] cs = col_name.split(",");
						int pt = Collections.binarySearch(col_names, cs, new Comparator<String[]>() {
							@Override
							public int compare(String[] s1, String[] s2) {
								int res = 0;
								for (int i = 0; i < s1.length; i++) {
									res = Integer.compare(Integer.parseInt(s1[i]), Integer.parseInt(s2[i]));
									if (res != 0) {
										return res;
									}
								}
								return res;
							}
						});
						if (pt < 0) {
							col_names.add(~pt, cs);
						}
					}

				}

				// Header
				pWri.print("Time");
				for (String[] col_name_split : col_names) {
					pWri.print(',');
					pWri.printf("Inf_%s_Site_%s_Stage_%s", col_name_split[0], col_name_split[1], col_name_split[2]);
				}
				pWri.println();

				for (Integer time : timeArr) {
					HashMap<String, Integer> ent = infected_site_stage_count.get(time);
					pWri.print(time);
					for (int i = 0; i < col_names.size(); i++) {
						pWri.print(',');
						String ckey = String.format("%s,%s,%s", col_names.get(i)[0], col_names.get(i)[1],
								col_names.get(i)[2]);
						Integer val = ent.get(ckey);
						if (val == null) {
							val = 0;
						}
						pWri.print(val);
					}
					pWri.println();
				}
				pWri.close();

			} catch (IOException e) {
				e.printStackTrace(System.err);
			}

		}

		if (print_progress != null && runnableId != null) {
			try {
				print_progress.printf("Post simulation file generation for Thread <%s> completed. Timestamp = %tc.\n",
						runnableId, System.currentTimeMillis());
			} catch (Exception ex) {
				System.err.printf("Post simulation file generation for Thread <%s> completed.\n", runnableId);
			}
		}

	}

	@Override
	public ArrayList<Integer> loadOptParameter(String[] parameter_settings, double[] point, int[][] seedInfectNum,
			boolean display_only) {

		ArrayList<String> common_parameter_name = new ArrayList<>();
		ArrayList<Double> common_parameter_value = new ArrayList<>();

		for (int i = 0; i < parameter_settings.length; i++) {
			if (PROP_PEP_START_AT.equals(parameter_settings[i])) {
				this.prophylaxis_starts_at = (int) point[i];
			} else if (parameter_settings[i].startsWith(PROP_PEP_UPTAKE)) {
				Matcher m = Pattern.compile(PROP_PEP_UPTAKE + "_(\\d+)").matcher(parameter_settings[i]);
				if (m.matches()) {
					prophylaxis_uptake[Integer.parseInt(m.group(1))] = (float) point[i];
				}
			} else if (parameter_settings[i].startsWith(PROP_PEP_PERSISTENCE_ADHERENCE)) {
				Matcher m = Pattern.compile(PROP_PEP_PERSISTENCE_ADHERENCE + "_(\\d+)").matcher(parameter_settings[i]);
				if (m.matches()) {
					prophylaxis_persistence_adherence[Integer.parseInt(m.group(1))] = point[i];
					double[] persistence = Arrays.copyOf(prophylaxis_persistence_adherence, 2);
					if (persistence[1] == 0) {
						persistenceDist = generateNonDistribution(rng_PEP, persistence);
					} else if (persistence[1] < 0) {
						double[] param = Arrays.copyOf(persistence, 2);
						param[1] = Math.abs(param[1]);
						persistenceDist = generateUniformDistribution(rng_PEP, param);
					} else {
						persistenceDist = generateGammaDistribution(rng_PEP, persistence);
					}
				}

			} else {
				common_parameter_name.add(parameter_settings[i]);
				common_parameter_value.add(point[i]);
			}
		}

		Double[] common_parameter_val_obj = common_parameter_value.toArray(new Double[common_parameter_value.size()]);

		double[] common_parameter_val = new double[common_parameter_value.size()];
		for (int i = 0; i < common_parameter_val.length; i++) {
			common_parameter_val[i] = common_parameter_val_obj[i].doubleValue();
		}

		return super.loadOptParameter(common_parameter_name.toArray(new String[common_parameter_name.size()]),
				common_parameter_val, seedInfectNum, display_only);
	}

}
