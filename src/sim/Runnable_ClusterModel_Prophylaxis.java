package sim;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.regex.Matcher;
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
	private static final int num_act = 5;
	private RealDistribution adherenceDist;

	protected double[] prophylaxis_adherence;
	protected float prophylaxis_uptake_HIV_PrEP;
	protected float prophylaxis_uptake_last_TP;
	protected float prophylaxis_uptake_last_STI;
	protected float prophylaxis_uptake_num_partners;
	protected float prophylaxis_uptake_num_partners_limit;

	protected RandomGenerator rng_PEP;
	protected HashMap<Integer, int[]> dx_last_12_months;
	protected HashMap<Integer, int[][]> partners_last_12_months;

	public static final String PROP_PEP_ADHERENCE = "PROP_PEP_ADHERENCE";
	public static final String PROP_PEP_UPTAKE = "PROP_PEP_UPTAKE";

	protected static final int PEP_AVAIL_TP = 1; // At least once
	protected static final int PEP_AVAIL_ANY_STI = 2; // More than twice

	// FILENAME_PREVALENCE_PERSON, FILENAME_CUMUL_INCIDENCE_PERSON
	private static final int[] COL_SEL_INF_GENDER = new int[] { 2, 6, 10 };
	// FILENAME_CUMUL_INCIDENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE = new int[] { 8, 25, 26, 27, 41, 42, 43 };
	// "Infectious_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_SITE = new int[] { 0, 5, 6, 7, 9, 10, 11 };
	// "Infected_" + Simulation_ClusterModelTransmission.FILENAME_PREVALENCE_SITE
	private static final int[] COL_SEL_INF_GENDER_SITE_AT = null; // new int[] { 64, 65, 192, 194, 196, 198, 200, 202,
																	// 204 160, 162, 164, 166, 224, 225 };

	private static final String SIM_OUTPUT_KEY_PEP_COVERAGE = "SIM_OUTPUT_PEP_COVERAGE";

	public Runnable_ClusterModel_Prophylaxis(long cMap_seed, long sim_seed, ContactMap base_cMap, Properties prop) {
		super(cMap_seed, sim_seed, base_cMap, prop, num_inf, num_site, num_act);
		this.prophylaxis_starts_at = Integer.parseInt(prop.getProperty(PROP_PEP_START_AT, "-1"));
		this.prophylaxis_adherence = (double[]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PEP_ADHERENCE, "[0,0]"), double[].class);
		float[] update_rate = (float[]) PropValUtils
				.propStrToObject(prop.getProperty(PROP_PEP_UPTAKE, "[0.0,0.0,0.0,0.0,0.0]"), float[].class);
		this.prophylaxis_uptake_HIV_PrEP = update_rate[0];
		this.prophylaxis_uptake_last_TP = update_rate[1];
		this.prophylaxis_uptake_last_STI = update_rate[2];
		this.prophylaxis_uptake_num_partners = update_rate[3];
		this.prophylaxis_uptake_num_partners_limit = update_rate[4];

		rng_PEP = new MersenneTwisterRandomGenerator(sim_seed);
		dx_last_12_months = new HashMap<>();
		partners_last_12_months = new HashMap<>();

		if (prophylaxis_adherence.length == 1) {
			adherenceDist = generateNonDistribution(rng_PEP, prophylaxis_adherence);
		} else if (prophylaxis_adherence[1] < 0) {
			double[] param = Arrays.copyOf(prophylaxis_adherence, 2);
			param[1] = Math.abs(param[1]);
			adherenceDist = generateUniformDistribution(rng_PEP, param);
		} else {
			adherenceDist = generateGammaDistribution(rng_PEP, prophylaxis_adherence);
		}
	}

	@Override
	protected void addPartnership(ContactMap cMap, Integer[] edge) {
		super.addPartnership(cMap, edge);

		if (prophylaxis_uptake_num_partners > 0 && prophylaxis_uptake_num_partners_limit > 0) {
			// Only keep partner history if needed
			for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
					Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
				int pid = edge[index];
				int partnerId = edge[index == Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1
						? Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2
						: Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1];

				int recIndec = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
						% AbstractIndividualInterface.ONE_YEAR_INT;

				int[][] partner_hist = partners_last_12_months.get(pid);
				if (partner_hist == null) {
					partner_hist = new int[AbstractIndividualInterface.ONE_YEAR_INT][];
					partners_last_12_months.put(pid, partner_hist);
				}

				if (partner_hist[recIndec] == null) {
					partner_hist[recIndec] = new int[] { partnerId };
				} else {
					partner_hist[recIndec] = Arrays.copyOf(partner_hist[recIndec], partner_hist[recIndec].length + 1);
					partner_hist[recIndec][partner_hist[recIndec].length - 1] = partnerId;
				}
			}
		}
	}

	@Override
	protected void testPerson(int currentTime, int pid, int infIncl, int siteIncl, int[][] cumul_treatment_by_person) {
		super.testPerson(currentTime, pid, infIncl, siteIncl, cumul_treatment_by_person);

		// Doxy-PEP allocation by test
		if (this.prophylaxis_starts_at > 0 && currentTime >= this.prophylaxis_starts_at) {

			boolean allocatePEP = false;

			// Risk group based PEP
			if (!allocatePEP && prophylaxis_uptake_HIV_PrEP > 0) {
				if (!allocatePEP && risk_cat_map.get(pid).intValue() == 0 || risk_cat_map.get(pid).intValue() == 1) {
					allocatePEP |= rng_PEP.nextFloat() < prophylaxis_uptake_HIV_PrEP;
				}
			}

			// TODO: Partner based PEP
			if (!allocatePEP && prophylaxis_uptake_num_partners > 0) {
				int[][] partner_hist = partners_last_12_months.get(pid);
				if (partner_hist != null) {
					ArrayList<Integer> partnerList = new ArrayList<>();
					for (int i = 0; i < partner_hist.length; i++) {
						if(partner_hist[i] != null) {
							for(Integer part : partner_hist[i]) {
								int pt = Collections.binarySearch(partnerList, part);
								if(pt < 0) {
									partnerList.add(~pt, part);
								}																
							}
						}
					}					
					if(partnerList.size() > prophylaxis_uptake_num_partners_limit) {
						allocatePEP |= rng_PEP.nextFloat() < prophylaxis_uptake_num_partners;
					}					
				}
			}

			// Treatment rate based PEP
			if (!allocatePEP && (prophylaxis_uptake_last_TP > 0 || prophylaxis_uptake_last_STI > 0)) {
				int[] dx_hist = dx_last_12_months.get(pid);
				if (dx_hist != null) {
					if (dx_hist[currentTime % dx_hist.length] != 0) { // Has a positive dx today
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
						if (!allocatePEP && num_dx_12_months[0] >= PEP_AVAIL_TP && prophylaxis_uptake_last_TP > 0) {
							allocatePEP |= prophylaxis_uptake_last_TP >= 1
									|| rng_PEP.nextFloat() < prophylaxis_uptake_last_TP;
						}
						if (!allocatePEP && num_dx_12_month_any > PEP_AVAIL_ANY_STI
								&& prophylaxis_uptake_last_STI > 0) {
							allocatePEP |= prophylaxis_uptake_last_STI >= 1
									|| rng_PEP.nextFloat() < prophylaxis_uptake_last_STI;
						}
						// Special case with no other STI transmission
						if (!allocatePEP && prophylaxis_uptake_last_STI < 0) {
							final float[][] sti_incident_all = new float[][] { // HIV-, HIV+
									new float[] { 23.9f, 31.5f }, new float[] { 29.2f, 40.9f } };
							float[] sti_incident = risk_cat_map.get(pid).intValue() == 0 ? sti_incident_all[1]
									: sti_incident_all[0];
							for (int i = 0; i < sti_incident.length && num_dx_12_month_any <= 2; i++) {
								if (rng_PEP.nextFloat() < sti_incident[i] / 100f) {
									num_dx_12_month_any++;
								}
							}
							allocatePEP = num_dx_12_month_any > 2;
						}
					}
				}

			} // End of treatment based PEP

			if (allocatePEP) {
				allocateProphylaxis(currentTime, pid);
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

		for (Integer hasPartnership : partners_last_12_months.keySet()) {
			int[][] partner_hist = partners_last_12_months.get(hasPartnership);
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
			int[] stat = new int[2]; // Ever PEP, Currently on PEP
			stat[0] = prophylaxis_record.size();
			for (int pid : prophylaxis_record.keySet()) {
				int[] prop_rec = prophylaxis_record.get(pid);
				if (prop_rec[PROPHYLAXIS_REC_PROTECT_UNTIL] > currentTime) {
					stat[1]++;
				}

			}
			pep_coverage.put(currentTime, stat);

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

		// PEP Usage
		countMap = (HashMap<Integer, int[]>) sim_output.get(SIM_OUTPUT_KEY_PEP_COVERAGE);
		if (countMap != null) {
			fileName = String.format(filePrefix + "PEP_Stat_%d_%d.csv", cMAP_SEED, sIM_SEED);
			printCountMap(countMap, fileName, "PEP_User_Type_%d", new int[] { 2 });
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
	public ArrayList<Integer> loadOptParamter(String[] parameter_settings, double[] point, int[][] seedInfectNum,
			boolean display_only) {

		ArrayList<String> common_parameter_name = new ArrayList<>();
		ArrayList<Double> common_parameter_value = new ArrayList<>();

		for (int i = 0; i < parameter_settings.length; i++) {
			if (PROP_PEP_START_AT.equals(parameter_settings[i])) {
				this.prophylaxis_starts_at = (int) point[i];
			} else if (parameter_settings[i].startsWith(PROP_PEP_UPTAKE)) {
				Matcher m = Pattern.compile(PROP_PEP_UPTAKE + "_(\\d+)").matcher(parameter_settings[i]);
				if (m.matches()) {
					switch (Integer.parseInt(m.group(1))) {
					case 0:
						prophylaxis_uptake_HIV_PrEP = (float) point[i];
						break;
					case 1:
						prophylaxis_uptake_last_TP = (float) point[i];
						break;
					case 2:
						prophylaxis_uptake_last_STI = (float) point[i];
						break;
					case 3:
						prophylaxis_uptake_num_partners = (float) point[i];
						break;
					case 4:
						prophylaxis_uptake_num_partners_limit = (float) point[i];
						break;
					default:
						System.err.printf("PROP_PEP_UPTAKE of value %s not defined.\n", parameter_settings[i]);
					}
				}
			} else if (parameter_settings[i].startsWith(PROP_PEP_ADHERENCE)) {
				Matcher m = Pattern.compile(PROP_PEP_ADHERENCE + "_(\\d+)").matcher(parameter_settings[i]);
				if (m.matches()) {
					int adhereId = Integer.parseInt(m.group(1));
					if (adhereId >= prophylaxis_adherence.length) {
						prophylaxis_adherence = Arrays.copyOf(prophylaxis_adherence, adhereId + 1);
					}
					prophylaxis_adherence[adhereId] = point[i];

					if (prophylaxis_adherence.length == 1) {
						adherenceDist = generateNonDistribution(rng_PEP, prophylaxis_adherence);
					} else if (prophylaxis_adherence[1] < 0) {
						double[] param = Arrays.copyOf(prophylaxis_adherence, 2);
						param[1] = Math.abs(param[1]);
						adherenceDist = generateUniformDistribution(rng_PEP, param);
					} else {
						adherenceDist = generateGammaDistribution(rng_PEP, prophylaxis_adherence);
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

		return super.loadOptParamter(common_parameter_name.toArray(new String[common_parameter_name.size()]),
				common_parameter_val, seedInfectNum, display_only);
	}

}
