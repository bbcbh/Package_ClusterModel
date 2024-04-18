package sim;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.RealDistribution;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PropValUtils;

public class Runnable_ClusterModel_Prophylaxis extends Abstract_Runnable_ClusterModel_MultiTransmission_Prophylaxis {

	public static final Pattern PROP_TYPE_PATTERN = Pattern.compile("ClusterModel_Prophylaxis");

	private static final int num_inf = 3; // TP, NG and CT
	private static final int num_site = 4;
	private static final int num_act = 3;
	private RealDistribution adherenceDist;

	protected double[] prophylaxis_adherence;
	protected float prophylaxis_uptake_HIV_PrEP;
	protected float prophylaxis_uptake_last_TP;
	protected float prophylaxis_uptake_last_STI;
	protected RandomGenerator rng_PEP;
	protected HashMap<Integer, int[]> dx_last_12_months;

	public static final String PROP_PEP_ADHERENCE = "PROP_PEP_ADHERENCE";
	public static final String PROP_PEP_UPTAKE = "PROP_PEP_UPTAKE";
	
	// FILENAME_PREVALENCE_PERSON, FILENAME_CUMUL_INCIDENCE_PERSON
	private static final int[] COL_SEL_INF_GENDER = new int[] { 2, 6, 10 };	
	// FILENAME_CUMUL_INCIDENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE = new int[] { 8, 25, 26, 27, 41 , 42, 43 };
	// "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_SITE = new int[] { 0, 5, 6, 7, 9, 10, 11};
	// "Infected_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE_AT =  null; //new int[] { 64, 65, 192, 194, 196, 198, 200, 202, 204  160, 162, 164, 166,	224, 225 };

	public Runnable_ClusterModel_Prophylaxis(long cMap_seed, long sim_seed, ContactMap base_cMap, Properties prop) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);
		this.prophylaxis_starts_at = Integer.parseInt(prop.getProperty(PROP_PEP_START_AT, "-1"));
		this.prophylaxis_adherence = (double[]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PEP_ADHERENCE, "[0,0]"), double[].class);
		float[] update_rate = (float[]) PropValUtils.propStrToObject(prop.getProperty(PROP_PEP_UPTAKE, "[0.0,0.0,0.0]"),
				float[].class);
		this.prophylaxis_uptake_HIV_PrEP = update_rate[0];
		this.prophylaxis_uptake_last_TP = update_rate[1];
		this.prophylaxis_uptake_last_STI = update_rate[2];

		rng_PEP = new MersenneTwisterRandomGenerator(sim_seed);
		dx_last_12_months = new HashMap<>();

		if (prophylaxis_adherence.length == 1) {
			adherenceDist = generateNonDistribution(prophylaxis_adherence);
		} else if (prophylaxis_adherence[1] < 0) {
			double[] param = Arrays.copyOf(prophylaxis_adherence, 2);
			param[1] = Math.abs(param[1]);
			adherenceDist = generateUniformDistribution(param);
		} else {
			adherenceDist = generateGammaDistribution(prophylaxis_adherence);
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

		// Risk group based PEP
		if (currentTime == this.prophylaxis_starts_at && prophylaxis_uptake_HIV_PrEP > 0) {
			for (Integer pid : bASE_CONTACT_MAP.vertexSet()) {
				if (risk_cat_map.get(pid).intValue() == 0 || risk_cat_map.get(pid).intValue() == 1) {
					if (prophylaxis_uptake_HIV_PrEP >= 1 || rng_PEP.nextFloat() < prophylaxis_uptake_HIV_PrEP) {
						allocateProphylaxis(currentTime, pid);
					}
				}
			}
		}

		// Treatment rate based PEP

		for (Integer treated_previously : dx_last_12_months.keySet()) {
			int[] dx_hist = dx_last_12_months.get(treated_previously);
			if (dx_hist != null) {

				if (currentTime == this.prophylaxis_starts_at
						&& (prophylaxis_uptake_last_TP > 0 || prophylaxis_uptake_last_STI > 0)) {
					int[] num_dx_12_months = new int[num_inf];
					int num_dx_12_month_any = 0;
					for (int dx_daily_record : dx_hist) {
						for (int i = 0; i < num_dx_12_months.length; i++) {
							if ((dx_daily_record & 1 << i) != 0) {
								num_dx_12_months[i]++;
								num_dx_12_month_any++;
							}
						}
					}

					boolean allocateProphylaxis = prophylaxis_record.containsKey(treated_previously);

					if (!allocateProphylaxis && num_dx_12_months[0] > 0 && prophylaxis_uptake_last_TP > 0) {
						allocateProphylaxis |= prophylaxis_uptake_last_TP >= 1
								|| rng_PEP.nextFloat() < prophylaxis_uptake_last_TP;
					}

					if (!allocateProphylaxis && num_dx_12_month_any > 2 && prophylaxis_uptake_last_STI > 0) {
						allocateProphylaxis |= prophylaxis_uptake_last_STI >= 1
								|| rng_PEP.nextFloat() < prophylaxis_uptake_last_STI;
					}

					if (!allocateProphylaxis && prophylaxis_uptake_last_STI < 0) { // Special case with no other STI
																					// transmission
						final float[][] sti_incident_all = new float[][] { // HIV-, HIV+
								new float[] { 23.9f, 31.5f }, new float[] { 29.2f, 40.9f } };
						float[] sti_incident = risk_cat_map.get(treated_previously).intValue() == 0
								? sti_incident_all[1]
								: sti_incident_all[0];
						for (int i = 0; i < sti_incident.length && num_dx_12_month_any <= 2; i++) {
							if (rng_PEP.nextFloat() < sti_incident[i] / 100f) {
								num_dx_12_month_any++;
							}
						}
						allocateProphylaxis = num_dx_12_month_any > 2;
					}
					if (allocateProphylaxis) {
						allocateProphylaxis(currentTime, treated_previously);
					}
				}
				// Reset record
				dx_hist[currentTime % dx_hist.length] = 0;
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
		prop_rec[PROPHYLAXIS_REC_LAST_UPTAKE_AT] = currentTime; // Not used
		prop_rec[PROPHYLAXIS_REC_DOSAGE] = Integer.MAX_VALUE; // Infinite in this model
		prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] = currentTime + (int) Math.round(adherenceDist.sample());
	}
	
	@SuppressWarnings("unchecked")
	protected void postSimulation() {
		String key, fileName;
		HashMap<Integer, int[]> countMap;
		String filePrefix = this.getRunnableId() == null ? "" : this.getRunnableId();
		
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
		
		
		/*
		String[] zipFormat = new String[] { "Infected_All_Stages_Prevalence_Site_", 
				"Infected_Prevalence_Site_", "Infectious_Prevalence_Person_", "Infectious_Prevalence_Site_"};		
		String[] stat_file_name = new String[] { null, null, null, null};
		
		try {
			Simulation_ClusterModelTransmission.output_analysis_csv(baseDir, zipFormat, stat_file_name);
		} catch (IOException e) {			
			e.printStackTrace(System.err);
		}
		*/
		
		
	}

}
