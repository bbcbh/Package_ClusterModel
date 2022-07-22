package sim;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import population.person.Person_Bridging_Pop;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;

public class Runnable_ClusterModel_Transmission extends Abstract_Runnable_ClusterModel {

	public static final int ACT_INDEX_GENITAL = 0;
	public static final int ACT_INDEX_ANAL = ACT_INDEX_GENITAL + 1;
	public static final int ACT_INDEX_FELLATIO = ACT_INDEX_ANAL + 1;

	public static final int SITE_VAGINA = 0;
	public static final int SITE_PENIS = SITE_VAGINA + 1;
	public static final int SITE_RECTUM = SITE_PENIS + 1;
	public static final int SITE_OROPHARYNX = SITE_RECTUM + 1;
	public static final int LENGTH_SITE = SITE_OROPHARYNX + 1;

	public static final String FILENAME_FORMAT_INDEX_CASE_LIST = "Seed_%d_IndexCases.txt";

	// Transmission
	private double[] DEFAULT_TRANS_V2P = new double[] { 0.4, 0.10 };
	private double[] DEFAULT_TRANS_P2V = new double[] { 0.2, 0.05 };
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private double[] DEFAULT_TRANS_P2R = new double[] { 0.63, 0 };
	private double[] DEFAULT_TRANS_R2P = new double[] { 0.02, 0 };
	private double[] DEFAULT_TRANS_P2O = new double[] { 0.44, 0 };
	private double[] DEFAULT_TRANS_O2P = new double[] { 0.01, 0 };

	// Duration
	private double[] DEFAULT_INFECTIOUS_PERIOD_VAGINA = new double[] { 15 * 7, 5 * 7 };
	private double[] DEFAULT_INFECTIOUS_PERIOD_PENIS = new double[] { 15 * 7, 5 * 7 };
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private double[] DEFAULT_INFECTIOUS_PERIOD_RECTUM = new double[] { 307.2, 5.54 };
	private double[] DEFAULT_INFECTIOUS_PERIOD_OROPHARYNX = new double[] { 80.4, 5.67 };

	// Incubation
	private double[] DEFAULT_INCUBATION_RANGE = new double[] { 3, 6 }; // 3 - 5 days

	// Sought test period
	// From Qibin's paper of Gamma(3, 0.86)
	private double[] DEFAULT_SYM_TEST = new double[] { 3 * 0.86, Math.sqrt(3 * 0.86 * 0.86) };

	// Act frequency
	// ASHR2: Those who were in heterosexual relationships had had
	// sex on average 1.84 times a week in the past 4 weeks,
	private float DEFAULT_ACT_GENITAL_FREQ = 1.84f / 7;
	// ASHR2 : only 1% of men and 0.4% of women had anal intercourse
	// the last time they had sex with a partner of the other sex
	private float DEFAULT_ACT_ANAL_FREQ_HETRO = 0.1f / 100;
	// ASHR2: oral sex was reported in only approximately one in four encounters
	private float DEFAULT_ACT_FELLATIO_FREQ_HETRO = 0.25f;
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private float DEFAULT_ACT_ANAL_FREQ_MSM = (1.6f + 2.4f) / 2 / 7;
	private float DEFAULT_ACT_FELLATIO_FREQ_MSM = (1.6f + 2.4f) / 2 / 7;

	private float[] DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM = new float[] { 20 };
	private float[] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_LOW_RISK = new float[] { 0.375f, 1, 720, 360, 90 };
	private float[] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_HIGH_RISK = new float[] { 0.05f, 0.49f, 0.71f, 0.99f, 720,
			360, 180, 120, 90 };
	private float[][] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM = new float[][] {
			DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_LOW_RISK, DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_HIGH_RISK };

	private float[] DEFAULT_SYM_RATE_BY_SITE_MALE = new float[] { Float.NaN, 0.9f, 0.12f, 0 };
	private float[] DEFAULT_SYM_RATE_BY_SITE_FEMALE = new float[] { 0.4f, Float.NaN, 0.12f, 0 };

	public static final int RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ = 0;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE = RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ + 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD = RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_INCUBATION_PERIOD = RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SYM_RATE = RUNNABLE_FIELD_TRANSMISSION_INCUBATION_PERIOD + 1;

	public static final int RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS = RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES = RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM = RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES
			+ 1;
	public static final int LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD = RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
			+ 1;

	public Object[] runnable_fields = {
			// RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ
			// double[ACT_TYPE][GENDER_FROM][GENDER_TO]
			new float[][][] {
					// ACT_INDEX_GENITAL
					new float[][] { new float[] { 0, DEFAULT_ACT_GENITAL_FREQ, 0, DEFAULT_ACT_GENITAL_FREQ },
							new float[] { DEFAULT_ACT_GENITAL_FREQ, 0, 0, 0 }, new float[] { 0, 0, 0, 0 },
							new float[] { DEFAULT_ACT_GENITAL_FREQ, 0, 0, 0 }, },
					// ACT_INDEX_ANAL
					new float[][] { new float[] { 0, DEFAULT_ACT_ANAL_FREQ_HETRO, 0, DEFAULT_ACT_ANAL_FREQ_HETRO },
							new float[] { DEFAULT_ACT_ANAL_FREQ_HETRO, DEFAULT_ACT_ANAL_FREQ_HETRO, 0, 0 },
							new float[] { 0, 0, DEFAULT_ACT_ANAL_FREQ_MSM, DEFAULT_ACT_ANAL_FREQ_MSM },
							new float[] { DEFAULT_ACT_ANAL_FREQ_HETRO, 0, DEFAULT_ACT_ANAL_FREQ_MSM,
									DEFAULT_ACT_ANAL_FREQ_MSM }, },
					// ACT_INDEX_FELLATIO
					new float[][] {
							new float[] { 0, DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, DEFAULT_ACT_FELLATIO_FREQ_HETRO },
							new float[] { DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, 0, 0 },
							new float[] { 0, 0, DEFAULT_ACT_FELLATIO_FREQ_MSM, DEFAULT_ACT_FELLATIO_FREQ_MSM },
							new float[] { DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, DEFAULT_ACT_FELLATIO_FREQ_MSM,
									DEFAULT_ACT_FELLATIO_FREQ_MSM }, },

			},
			// RUNNABLE_FIELD_TRANSMISSION_MAP_TRANSMISSION_RATE
			// double[SITE_FROM][SITE_TO][] =
			new double[][][] {
					// FROM SITE_VAGINA
					new double[][] { null, DEFAULT_TRANS_V2P, null, null },
					// FROM SITE_PENIS
					new double[][] { DEFAULT_TRANS_P2V, null, DEFAULT_TRANS_P2R, DEFAULT_TRANS_P2O },
					// FROM SITE_RECTUM
					new double[][] { null, DEFAULT_TRANS_R2P, null, null },
					// FROM SITE_OROPHARYNX
					new double[][] { null, DEFAULT_TRANS_O2P, null, null }, },

			// RUNNABLE_FIELD_INFECTIOUS_PERIOD
			// double[SITE]
			new double[][] { DEFAULT_INFECTIOUS_PERIOD_VAGINA, DEFAULT_INFECTIOUS_PERIOD_PENIS,
					DEFAULT_INFECTIOUS_PERIOD_RECTUM, DEFAULT_INFECTIOUS_PERIOD_OROPHARYNX, },
			// RUNNABLE_FIELD_INCUBATION_PERIOD
			// double[SITE]
			new double[][] { DEFAULT_INCUBATION_RANGE, DEFAULT_INCUBATION_RANGE, DEFAULT_INCUBATION_RANGE,
					DEFAULT_INCUBATION_RANGE },
			// RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			new float[][] { DEFAULT_SYM_RATE_BY_SITE_FEMALE, DEFAULT_SYM_RATE_BY_SITE_MALE,
					DEFAULT_SYM_RATE_BY_SITE_MALE, DEFAULT_SYM_RATE_BY_SITE_MALE },
			// RUNNABLE_FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			new float[][] { null, null, DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM,
					DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM },
			// RUNNABLE_FIELD_TESTING_RATE_BY_RISK_CATEGORIES
			new float[][][] { null, null, DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM,
					DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM },
			// RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
			DEFAULT_SYM_TEST,

	};

	// FIELD_POP_COMPOSITION
	// int[] {NUM_FEMALE, NUM_MALE, NUM_MSMO, NUM_MSMW}
	final int[] cumulative_pop_composition;
	final ContactMap BASE_CONTACT_MAP;
	final int NUM_TIME_STEPS_PER_SNAP;
	final int SNAP_FREQ;
	final long seed;

	protected RandomGenerator RNG;

	protected transient RealDistribution[][] tranmissionMatrix = new RealDistribution[LENGTH_SITE][LENGTH_SITE];
	protected transient RealDistribution[] infectious_period = new RealDistribution[LENGTH_SITE];
	protected transient RealDistribution[] incubation_period = new RealDistribution[LENGTH_SITE];
	protected transient RealDistribution sym_test_period;

	protected transient ArrayList<Integer>[] currently_infectious;

	// For schedule, key are day and entries are list of person id, ordered
	protected transient HashMap<Integer, ArrayList<Integer>>[] schedule_incubation;
	protected transient HashMap<Integer, ArrayList<Integer>>[] schedule_recovery;
	protected transient HashMap<Integer, ArrayList<Integer>> schedule_testing;
	protected transient HashMap<Integer, Integer[]> mapping_infection_schedule;

	protected transient HashMap<Integer, double[][]> trans_prob;

	protected transient HashMap<Integer, Integer> risk_cat_map;
	protected transient int firstSeedTime = Integer.MAX_VALUE;
	protected transient HashMap<String, Object> sim_output = null;

	public static final String SIM_OUTPUT_INFECTIOUS_COUNT = "SIM_OUTPUT_PREVALENCE";

	public Runnable_ClusterModel_Transmission(long seed, int[] POP_COMPOSITION, ContactMap BASE_CONTACT_MAP,
			int NUM_TIME_STEPS_PER_SNAP, int SNAP_FREQ) {
		super();

		this.cumulative_pop_composition = new int[POP_COMPOSITION.length];
		int offset = 0;

		for (int g = 0; g < this.cumulative_pop_composition.length; g++) {
			this.cumulative_pop_composition[g] = offset + POP_COMPOSITION[g];
			offset += POP_COMPOSITION[g];
		}

		this.BASE_CONTACT_MAP = BASE_CONTACT_MAP;
		this.NUM_TIME_STEPS_PER_SNAP = NUM_TIME_STEPS_PER_SNAP;
		this.SNAP_FREQ = SNAP_FREQ;
		this.seed = seed;

	}

	@Override
	public Object[] getRunnable_fields() {
		return runnable_fields;
	}

	@SuppressWarnings("unchecked")
	public void initialse() {
		RNG = new MersenneTwisterRandomGenerator(seed);

		// Transmission
		double[][][] tranParm = (double[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE];
		for (int sf = 0; sf < tranmissionMatrix.length; sf++) {
			for (int st = 0; st < tranmissionMatrix[sf].length; st++) {
				double[] param = tranParm[sf][st];
				if (param != null) {
					if (param[1] == 0) {
						tranmissionMatrix[sf][st] = generateNonDistribution(param);
					} else {
						tranmissionMatrix[sf][st] = generateBetaDistribution(param);
					}
				}
			}
		}

		// Duration & incubation
		double[][] durParam = (double[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD];
		double[][] incParam = (double[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_INCUBATION_PERIOD];
		for (int s = 0; s < infectious_period.length; s++) {
			infectious_period[s] = generateGammaDistribution(durParam[s]);
			incubation_period[s] = new UniformRealDistribution(RNG, incParam[s][0], incParam[s][1]);
		}
		sym_test_period = generateGammaDistribution(
				(double[]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM]);

		// Lists
		currently_infectious = new ArrayList[LENGTH_SITE];
		schedule_incubation = new HashMap[LENGTH_SITE];
		schedule_recovery = new HashMap[LENGTH_SITE];

		for (int s = 0; s < LENGTH_SITE; s++) {
			currently_infectious[s] = new ArrayList<>();
			schedule_incubation[s] = new HashMap<>();
			schedule_recovery[s] = new HashMap<>();

		}

		trans_prob = new HashMap<>();
		sim_output = new HashMap<>();
		schedule_testing = new HashMap<>();
		mapping_infection_schedule = new HashMap<>();
		risk_cat_map = new HashMap<>();
	}

	public HashMap<String, Object> getSim_output() {
		return sim_output;
	}

	public void scheduleNextTest(Integer personId, int lastTestTime) {
		int genderType = getGenderType(personId);

		float[][] testRate = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES])[genderType];
		if (testRate != null) {
			int riskCat = Math.max(0, getRiskCategories(personId, genderType));
			float[] testRateByCat = testRate[riskCat];
			int divder = (testRateByCat.length - 1) / 2;

			float p = RNG.nextFloat();
			int pI = Arrays.binarySearch(testRateByCat, 0, divder, p);
			if (pI < 0) {
				pI = ~pI;
			}

			if (pI < divder) {

				double testGapTime = (testRateByCat[divder + pI] + testRateByCat[divder + pI + 1]) / 2;
				testGapTime *= 1 + RNG.nextGaussian() / 10;

				Integer testTime = (int) Math.round(lastTestTime + (testGapTime));

				ArrayList<Integer> testEnt = schedule_testing.get(testTime);
				if (testEnt == null) {
					testEnt = new ArrayList<>();
					schedule_testing.put(testTime, testEnt);
				}
				int key = Collections.binarySearch(testEnt, personId);
				if (key < 0) {
					testEnt.add(~key, personId);
				}
			}

		}

	}

	private int getRiskCategories(Integer personId, int genderType) {

		if (risk_cat_map.containsKey(personId)) {
			return risk_cat_map.get(personId);
		} else {
			int riskCat = -1;
			float[] riskCatList = ((float[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS])[genderType];
			if (riskCatList != null) {
				int numCasual = 0;
				Set<Integer[]> edges = BASE_CONTACT_MAP.edgesOf(personId);
				for (Integer[] e : edges) {
					if (e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME] >= firstSeedTime
							&& e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME] < firstSeedTime
									+ NUM_TIME_STEPS_PER_SNAP
							&& e[Population_Bridging.CONTACT_MAP_EDGE_DURATION] <= 1) {
						numCasual++;
					}
				}

				float numCasual1Year = ((float) AbstractIndividualInterface.ONE_YEAR_INT) * numCasual
						/ NUM_TIME_STEPS_PER_SNAP;
				riskCat = Arrays.binarySearch(riskCatList, numCasual1Year);
				if (riskCat < 0) {
					riskCat = ~riskCat;
				}

			}

			risk_cat_map.put(personId, riskCat);

			return riskCat;
		}

	}

	public void allocateSeedInfection(int[][] num_infectioned_by_gender_site, int time) {

		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] pid_collection = new ArrayList[Population_Bridging.LENGTH_GENDER];
		for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
			pid_collection[g] = new ArrayList<>();
		}

		for (Integer v : BASE_CONTACT_MAP.vertexSet()) {
			pid_collection[getGenderType(v)].add(v);
		}

		for (int g = 0; g < Population_Bridging.LENGTH_GENDER; g++) {
			Integer[] candidate = pid_collection[g].toArray(new Integer[pid_collection[g].size()]);
			for (int s = 0; s < LENGTH_SITE; s++) {
				int numInfect = num_infectioned_by_gender_site[g][s];
				for (int p = 0; p < candidate.length && numInfect > 0; p++) {
					if (RNG.nextInt(candidate.length - p) < numInfect) {
						addInfectious(candidate[p], s, time, time + (int) Math.round(infectious_period[s].sample()));
						numInfect--;
					}
				}
			}
		}
	}

	public int addInfectious(Integer infectedId, int site, int infectious_time, int recoveredAt) {
		int key = Collections.binarySearch(currently_infectious[site], infectedId);
		if (key < 0) {
			currently_infectious[site].add(~key, infectedId);
			firstSeedTime = Math.min(firstSeedTime, infectious_time);

			// Recovery
			ArrayList<Integer> sch = schedule_recovery[site].get(recoveredAt);
			if (sch == null) {
				sch = new ArrayList<>();
				schedule_recovery[site].put(recoveredAt, sch);
			}

			int keyR = Collections.binarySearch(sch, infectedId);
			if (keyR < 0) {
				sch.add(~keyR, infectedId);
				updateScheduleMap(infectedId, LENGTH_SITE + site, recoveredAt);
			}

			int gender = getGenderType(infectedId);

			// Transmission probability
			if (!trans_prob.containsKey(infectedId)) {
				double[][] trans = new double[LENGTH_SITE][LENGTH_SITE];
				for (int sf = 0; sf < LENGTH_SITE; sf++) {
					boolean sample = (gender == Person_Bridging_Pop.GENDER_TYPE_FEMALE) ? sf != SITE_PENIS
							: sf != SITE_VAGINA;
					if (sample) {
						for (int st = 0; st < LENGTH_SITE; st++) {
							if (tranmissionMatrix[sf][st] != null) {
								trans[sf][st] = tranmissionMatrix[sf][st].sample();
							}
						}
					}
				}
				trans_prob.put(infectedId, trans);
			}

			// Schedule sym test if needed.
			float sym_rate = ((float[][]) getRunnable_fields()[RUNNABLE_FIELD_TRANSMISSION_SYM_RATE])[gender][site];
			if (sym_rate > 0) {
				if (RNG.nextFloat() < sym_rate) {
					int sym_test_day = (int) Math.round(infectious_time + sym_test_period.sample());
					sch = schedule_testing.get(sym_test_day);

					if (sch == null) {
						sch = new ArrayList<>();
						schedule_testing.put(sym_test_day, sch);
					}

					int keyS = Collections.binarySearch(sch, infectedId);
					if (keyS < 0) {
						keyS = Collections.binarySearch(sch, -infectedId);
						if (keyS < 0) {
							sch.add(~keyS, -infectedId);
						}
					}

				}

			}

		}
		return key;
	}

	public int removeInfected(Integer infectedId, int site) {
		int key = Collections.binarySearch(currently_infectious[site], infectedId);
		if (key >= 0) {
			currently_infectious[site].remove(key);
		}
		return key;
	}

	public void removeInfected(Integer infectedId) {
		for (int s = 0; s < LENGTH_SITE; s++) {
			removeInfected(infectedId, s);
		}
	}

	public int getGenderType(Integer personId) {
		return getGenderType(personId, cumulative_pop_composition);
	}

	public static int getGenderType(Integer personId, int[] cumul_pop_comp) {
		int index = Arrays.binarySearch(cumul_pop_comp, personId);

		if (index < 0) {
			return ~index;
		} else {
			return index;
		}
	}

	@Override
	public void run() {
		int startTime = firstSeedTime;

		if (startTime < Integer.MAX_VALUE) {

			Object[] simulation_store = preSimulation();
			
			// Make a copy of contact map for interaction
			ContactMap cMap;
			try {
				cMap = ContactMap.ContactMapFromFullString(BASE_CONTACT_MAP.toFullString());
			} catch (IOException e1) {
				cMap = BASE_CONTACT_MAP;
				e1.printStackTrace(System.err);
			}

			HashSet<Integer[]> removeEdges = new HashSet<>();
			// Schedule testing
			for (Integer personId : cMap.vertexSet()) {
				scheduleNextTest(personId, startTime);
			}

			int snap_index = 0;
			for (int currentTime = startTime; currentTime < startTime
					+ NUM_TIME_STEPS_PER_SNAP * SNAP_FREQ; currentTime++) {

				for (int site_src = 0; site_src < LENGTH_SITE; site_src++) {
					// Update infectious
					ArrayList<Integer> becomeInfectiousToday = schedule_incubation[site_src].remove(currentTime);

					if (becomeInfectiousToday != null) {
						for (Integer toInfectiousId : becomeInfectiousToday) {
							int recoveredAt = (int) Math.round(infectious_period[site_src].sample()) + currentTime;
							addInfectious(toInfectiousId, site_src, currentTime, recoveredAt);
							Integer[] schMap = mapping_infection_schedule.get(toInfectiousId);
							if (schMap != null) {
								schMap[site_src] = null;
							}

						}
					}

					// Update recovery
					ArrayList<Integer> recoveredToday = schedule_recovery[site_src].remove(currentTime);
					if (recoveredToday != null) {
						for (Integer toRecoveredId : recoveredToday) {
							removeInfected(toRecoveredId, site_src);
							mapping_infection_schedule.remove(toRecoveredId);
						}
					}

					// Transmission
					for (Integer infectious : currently_infectious[site_src]) {
						if (cMap.containsVertex(infectious)) {
							Set<Integer[]> edges = cMap.edgesOf(infectious);
							double[][] trans = trans_prob.get(infectious);

							for (Integer[] e : edges) {
								int startIndex = Population_Bridging.CONTACT_MAP_EDGE_START_TIME;

								while (startIndex < e.length) {
									int durationIndex = startIndex + 1;

									if (currentTime >= e[startIndex]
											&& currentTime < (e[startIndex] + e[durationIndex])) {

										int partner = e[Population_Bridging.CONTACT_MAP_EDGE_P1].equals(infectious)
												? e[Population_Bridging.CONTACT_MAP_EDGE_P2]
												: e[Population_Bridging.CONTACT_MAP_EDGE_P1];

										for (int site_target = 0; site_target < LENGTH_SITE; site_target++) {
											if (trans[site_src][site_target] != 0) {
												// Transmission is possible
												if (Collections.binarySearch(currently_infectious[site_target],
														partner) < 0) {

													int g_s = getGenderType(infectious);
													int g_t = getGenderType(partner);
													int actType;
													// Determine act type
													switch (site_src) {
													case SITE_VAGINA:
														actType = ACT_INDEX_GENITAL;
														break;
													case SITE_OROPHARYNX:
														actType = ACT_INDEX_FELLATIO;
														break;
													case SITE_RECTUM:
														actType = ACT_INDEX_ANAL;
														break;
													default:
														switch (site_target) {
														case SITE_VAGINA:
															actType = ACT_INDEX_GENITAL;
															break;
														case SITE_OROPHARYNX:
															actType = ACT_INDEX_FELLATIO;
															break;
														case SITE_RECTUM:
															actType = ACT_INDEX_ANAL;
															break;
														default:
															actType = -1;
														}
													}
													boolean transmitted = actType != -1;

													if (transmitted) {
														float actProb = ((float[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ])[actType][g_s][g_t];
														double transProb = trans[site_src][site_target];
														transmitted &= actProb > 0;
														if (transmitted) {
															transmitted &= RNG.nextDouble() < (actProb * transProb);
														}
													}

													if (transmitted) {
														transmission_success(currentTime, infectious, partner,
																site_target, actType, simulation_store);
													}
												}
											}
										}
									} else if (currentTime >= (e[startIndex] + e[durationIndex])
											&& startIndex + 2 > e.length) {
										removeEdges.add(e);
									}
									startIndex += 2;
								}
							}
						}
					}
				}

				// Testing
				ArrayList<Integer> testToday = schedule_testing.remove(currentTime);
				if (testToday != null) {
					for (Integer tId : testToday) {
						Integer test_pid = Math.abs(tId);
						// Remove from infectious
						for (int site = 0; site < LENGTH_SITE; site++) {
							int key = Collections.binarySearch(currently_infectious[site], test_pid);
							if (key >= 0) {
								currently_infectious[site].remove(key);
							}
						}
						// Remove from incubation and recovery
						Integer[] infection_schMap = mapping_infection_schedule.get(test_pid);
						if (infection_schMap != null) {
							for (int i = 0; i < infection_schMap.length; i++) {
								if (infection_schMap[i] != null) {
									Integer dateEnt = infection_schMap[i];
									int site = i % LENGTH_SITE;
									ArrayList<Integer> schArr;
									if (i < LENGTH_SITE) {
										schArr = schedule_incubation[site].get(dateEnt);
									} else {
										schArr = schedule_recovery[site].get(dateEnt);
									}
									int key = Collections.binarySearch(schArr, test_pid);
									if (key >= 0) {
										schArr.remove(key);
									}

								}
							}
							Arrays.fill(infection_schMap, null);
						}

						if (tId < 0) {
							// Schedule next test
							scheduleNextTest(test_pid, currentTime);
						}
					}
				}

				for (Integer[] e : removeEdges) {
					cMap.removeEdge(e);
				}

				// Storing snapshot infected in sim_output

				if (snap_index == 0) {
					@SuppressWarnings("unchecked")
					// K = time, V= int [gender][site]
					HashMap<Integer, int[][]> infectious_count_map = (HashMap<Integer, int[][]>) sim_output
							.get(SIM_OUTPUT_INFECTIOUS_COUNT);
					if (infectious_count_map == null) {
						infectious_count_map = new HashMap<>();
						sim_output.put(SIM_OUTPUT_INFECTIOUS_COUNT, infectious_count_map);
					}
					int[][] infectious_count = new int[cumulative_pop_composition.length][LENGTH_SITE];
					for (int site = 0; site < LENGTH_SITE; site++) {
						for (Integer infected_id : currently_infectious[site]) {
							int gender_type = getGenderType(infected_id);
							infectious_count[gender_type][site]++;
						}
					}
					infectious_count_map.put(currentTime, infectious_count);
				}

				snap_index = (snap_index + 1) % NUM_TIME_STEPS_PER_SNAP;

			}
			// End of simulations
			if (runnableId != null) {
				System.out.println(String.format("Thread <%s> completed.", runnableId));
			}
			postSimulation(simulation_store);

		}

	}

	protected Object[] preSimulation() {
		// Do nothing by default		
		return new Object[0];

	}

	/**
	 * Procedure that are called after simulations. *
	 * 
	 */
	protected void postSimulation(Object[] simulation_store) {
		// Do nothing by default
	}


	protected void transmission_success(int currentTime, Integer infectious, int partner, int site_target,
			int actType, Object[] simulation_store) {
		Integer incubation_end_at = currentTime + (int) incubation_period[site_target].sample();
		ArrayList<Integer> ent = schedule_incubation[site_target].get(incubation_end_at);
		if (ent == null) {
			ent = new ArrayList<>();
			schedule_incubation[site_target].put(incubation_end_at, ent);
		}
		int key = Collections.binarySearch(ent, partner);
		if (key < 0) {
			ent.add(~key, partner);
			updateScheduleMap(partner, site_target, incubation_end_at);
		}

	}

	

	protected void updateScheduleMap(int personId, int schMap_index, Integer schMap_ent) {

		Integer[] schMap = mapping_infection_schedule.get(personId);

		// key = person_id, ent = {incubation by site_0, ... recovery_by_site_0 ...}

		if (schMap == null) {
			schMap = new Integer[2 * LENGTH_SITE];
			Arrays.fill(schMap, null);
			mapping_infection_schedule.put(personId, schMap);
		}
		schMap[schMap_index] = schMap_ent;
	}

	protected AbstractRealDistribution generateNonDistribution(double[] input) {
		return new AbstractRealDistribution(RNG) {
			private static final long serialVersionUID = -4946118496555960005L;

			@Override
			public boolean isSupportUpperBoundInclusive() {
				return true;
			}

			@Override
			public boolean isSupportLowerBoundInclusive() {
				return true;
			}

			@Override
			public boolean isSupportConnected() {
				return true;
			}

			@Override
			public double getSupportUpperBound() {
				return input[0];
			}

			@Override
			public double getSupportLowerBound() {
				return input[0];
			}

			@Override
			public double getNumericalVariance() {
				return 0;
			}

			@Override
			public double getNumericalMean() {
				return input[0];
			}

			@Override
			public double density(double x) {
				return x == input[0] ? 1 : 0;
			}

			@Override
			public double cumulativeProbability(double x) {
				return x < input[0] ? 0 : 1;
			}

			@Override
			public double sample() {
				return input[0];
			}

		};

	}

	protected GammaDistribution generateGammaDistribution(double[] input) {
		// For Gamma distribution
		// shape = alpha = mean*mean / variance = mean / scale
		// scale = 1/ (beta or lambda or rate) = variance / mean;
		double[] res = new double[2];
		double var = input[1] * input[1];
		// rate or 1/ beta or lambda
		res[1] = var / input[0];
		// shape or alpha
		res[0] = input[0] / res[1];
		return new GammaDistribution(RNG, res[0], res[1]);
	}

	protected BetaDistribution generateBetaDistribution(double[] input) {
		// For Beta distribution,
		// alpha = mean*(mean*(1-mean)/variance - 1)
		// beta = (1-mean)*(mean*(1-mean)/variance - 1)
		double[] res = new double[2];
		double var = input[1] * input[1];
		double rP = input[0] * (1 - input[0]) / var - 1;
		// alpha
		res[0] = rP * input[0];
		// beta
		res[1] = rP * (1 - input[0]);
		return new BetaDistribution(RNG, res[0], res[1]);
	}

}
