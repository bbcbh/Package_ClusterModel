package sim;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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

public class Runnable_ContactMapTransmission extends Abstract_Runnable_ContactMap {

	public static final int ACT_INDEX_GENITAL = 0;
	public static final int ACT_INDEX_ANAL = ACT_INDEX_GENITAL + 1;
	public static final int ACT_INDEX_FELLATIO = ACT_INDEX_ANAL + 1;

	public static final int SITE_VAGINA = 0;
	public static final int SITE_PENIS = SITE_VAGINA + 1;
	public static final int SITE_RECTUM = SITE_PENIS + 1;
	public static final int SITE_OROPHARYNX = SITE_RECTUM + 1;
	public static final int LENGTH_SITE = SITE_OROPHARYNX + 1;

	public static final String FILENAME_FORMAT_TRANSMISSION_CMAP = "Seed_%s_TransmissionMap_%d.csv";
	public static final String DIRNAME_FORMAT_TRANSMISSION_CMAP = "TransMap_%d";
	public static final String FILENAME_FORMAT_INDEX_CASE_LIST = "Seed_%d_IndexCases.txt";

	public static final int TRANSMAP_EDGE_INFECTIOUS = 0;
	public static final int TRANSMAP_EDGE_SUSCEPTIBLE = TRANSMAP_EDGE_INFECTIOUS + 1;
	public static final int TRANSMAP_EDGE_START_TIME = TRANSMAP_EDGE_SUSCEPTIBLE + 1;
	public static final int TRANSMAP_EDGE_ACT_INVOLVED = TRANSMAP_EDGE_START_TIME + 1;
	public static final int LENGTH_TRANSMAP_EDGE = TRANSMAP_EDGE_ACT_INVOLVED + 1;

	// Transmission
	private static double[] DEFAULT_TRANS_V2P = new double[] { 0.4, 0.10 };
	private static double[] DEFAULT_TRANS_P2V = new double[] { 0.2, 0.05 };
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private static double[] DEFAULT_TRANS_P2R = new double[] { 0.63, 0 };
	private static double[] DEFAULT_TRANS_R2P = new double[] { 0.02, 0 };
	private static double[] DEFAULT_TRANS_P2O = new double[] { 0.44, 0 };
	private static double[] DEFAULT_TRANS_O2P = new double[] { 0.01, 0 };

	// Duration
	private static double[] DEFAULT_INFECTIOUS_PERIOD_VAGINA = new double[] { 15 * 7, 5 * 7 };
	private static double[] DEFAULT_INFECTIOUS_PERIOD_PENIS = new double[] { 15 * 7, 5 * 7 };
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private static double[] DEFAULT_INFECTIOUS_PERIOD_RECTUM = new double[] { 307.2, 5.54 };
	private static double[] DEFAULT_INFECTIOUS_PERIOD_OROPHARYNX = new double[] { 80.4, 5.67 };

	// Incubation
	private static double[] DEFAULT_INCUBATION_RANGE = new double[] { 3, 6 }; // 3 - 5 days

	// Act frequency
	// ASHR2: Those who were in heterosexual relationships had had
	// sex on average 1.84 times a week in the past 4 weeks,
	private static double DEFAULT_ACT_GENITAL_FREQ = 1.84 / 7;
	// ASHR2 : only 1% of men and 0.4% of women had anal intercourse
	// the last time they had sex with a partner of the other sex
	private static double DEFAULT_ACT_ANAL_FREQ_HETRO = 0.1 / 100;
	// ASHR2: oral sex was reported in only approximately one in four encounters
	private static double DEFAULT_ACT_FELLATIO_FREQ_HETRO = 0.25;
	// From Qibin's paper
	// 10.1371/journal.pcbi.1009385
	private static double DEFAULT_ACT_ANAL_FREQ_MSM = (1.6 + 2.4) / 2 / 7;
	private static double DEFAULT_ACT_FELLATIO_FREQ_MSM = (1.6 + 2.4) / 2 / 7;
	private static double[] DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM = new double[] { 20 };
	private static double[] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_LOW_RISK = new double[] { 0.375, 1, 720, 360, 90 };
	private static double[] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_HIGH_RISK = new double[] { 0.05, 0.49, 0.71, 1, 720,
			360, 180, 120, 90 };
	private static double[][] DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM = new double[][] {
			DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_LOW_RISK, DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM_HIGH_RISK };

	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ = 0;
	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_TRANSMISSION_RATE = RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_INFECTIOUS_PERIOD = RUNNABLE_FIELD_TRANSMISSION_MAP_TRANSMISSION_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_INCUBATION_PERIOD = RUNNABLE_FIELD_TRANSMISSION_MAP_INFECTIOUS_PERIOD
			+ 1;
	public static final int RUNNABLE_FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS = RUNNABLE_FIELD_TRANSMISSION_MAP_INFECTIOUS_PERIOD
			+ 1;
	public static final int RUNNABLE_FIELD_TESTING_BY_RISK_CATEGORIES = RUNNABLE_FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			+ 1;

	public Object[] runnable_fields = {
			// RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ
			// double[ACT_TYPE][GENDER_FROM][GENDER_TO]
			new double[][][] {
					// ACT_INDEX_GENITAL
					new double[][] { new double[] { 0, DEFAULT_ACT_GENITAL_FREQ, 0, DEFAULT_ACT_GENITAL_FREQ },
							new double[] { DEFAULT_ACT_GENITAL_FREQ, 0, 0, 0 }, new double[] { 0, 0, 0, 0 },
							new double[] { DEFAULT_ACT_GENITAL_FREQ, 0, 0, 0 }, },
					// ACT_INDEX_ANAL
					new double[][] { new double[] { 0, DEFAULT_ACT_ANAL_FREQ_HETRO, 0, DEFAULT_ACT_ANAL_FREQ_HETRO },
							new double[] { DEFAULT_ACT_ANAL_FREQ_HETRO, DEFAULT_ACT_ANAL_FREQ_HETRO, 0, 0 },
							new double[] { 0, 0, DEFAULT_ACT_ANAL_FREQ_MSM, DEFAULT_ACT_ANAL_FREQ_MSM },
							new double[] { DEFAULT_ACT_ANAL_FREQ_HETRO, 0, DEFAULT_ACT_ANAL_FREQ_MSM,
									DEFAULT_ACT_ANAL_FREQ_MSM }, },
					// ACT_INDEX_FELLATIO
					new double[][] {
							new double[] { 0, DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, DEFAULT_ACT_FELLATIO_FREQ_HETRO },
							new double[] { DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, 0, 0 },
							new double[] { 0, 0, DEFAULT_ACT_FELLATIO_FREQ_MSM, DEFAULT_ACT_FELLATIO_FREQ_MSM },
							new double[] { DEFAULT_ACT_FELLATIO_FREQ_HETRO, 0, DEFAULT_ACT_FELLATIO_FREQ_MSM,
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
			// RUNNABLE_FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			new double[][] { null, null, DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM,
					DEFAULT_RISK_CATEGORIES_CASUAL_PARNTERS_MSM },
			// RUNNABLE_FIELD_TESTING_BY_RISK_CATEGORIES
			new double[][][] { null, null, DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM,
					DEFAULT_TESTING_RATE_BY_CATEGORIES_MSM },

	};

	// FIELD_POP_COMPOSITION
	// int[] {NUM_FEMALE, NUM_MALE, NUM_MSMO, NUM_MSMW}
	final int[] cumulative_pop_composition;
	final ContactMap BASE_CONTACT_MAP;
	final int NUM_TIME_STEPS;
	final long seed;

	private ContactMap transmissionMap;
	private RandomGenerator RNG;

	protected transient RealDistribution[][] tranmissionMatrix = new RealDistribution[LENGTH_SITE][LENGTH_SITE];
	protected transient RealDistribution[] infectious_period = new RealDistribution[LENGTH_SITE];
	protected transient RealDistribution[] incubation_period = new RealDistribution[LENGTH_SITE];

	protected transient ArrayList<Integer>[] currently_infectious;
	protected transient HashMap<Integer, ArrayList<Integer>>[] incubation_schedule;
	protected transient HashMap<Integer, ArrayList<Integer>>[] recovery_schedule;

	protected transient HashMap<Integer, ArrayList<Integer>> testing_schedule;
	protected transient HashMap<Integer, double[][]> trans_prob;

	protected transient HashMap<Integer, Integer> risk_cat_map;
	protected transient int firstSeedTime = Integer.MAX_VALUE;

	private transient HashMap<String, Object> sim_output = null;

	public static final String SIM_OUTPUT_TRANMISSION_MAP = "SIM_OUTPUT_TRANMISSION_MAP";
	public static final String SIM_OUTPUT_CLUSTERS = "SIM_OUTPUT_CLUSTERS";

	public Runnable_ContactMapTransmission(long seed, int[] POP_COMPOSITION, ContactMap BASE_CONTACT_MAP,
			int NUM_TIME_STEPS) {
		super();

		this.cumulative_pop_composition = new int[POP_COMPOSITION.length];
		int offset = 0;

		for (int g = 0; g < this.cumulative_pop_composition.length; g++) {
			this.cumulative_pop_composition[g] = offset + POP_COMPOSITION[g];
			offset += POP_COMPOSITION[g];
		}

		this.BASE_CONTACT_MAP = BASE_CONTACT_MAP;
		this.NUM_TIME_STEPS = NUM_TIME_STEPS;
		this.seed = seed;

	}

	@Override
	public Object[] getRunnable_fields() {
		return runnable_fields;
	}

	@SuppressWarnings("unchecked")
	public void initialse() {
		transmissionMap = new ContactMap();

		RNG = new MersenneTwisterRandomGenerator(seed);

		// Transmission
		double[][][] tranParm = (double[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_MAP_TRANSMISSION_RATE];
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
		double[][] durParam = (double[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_MAP_INFECTIOUS_PERIOD];
		double[][] incParam = (double[][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_MAP_INCUBATION_PERIOD];
		for (int s = 0; s < infectious_period.length; s++) {
			infectious_period[s] = generateGammaDistribution(durParam[s]);
			incubation_period[s] = new UniformRealDistribution(RNG, incParam[s][0], incParam[s][1]);
		}

		// Lists
		currently_infectious = new ArrayList[LENGTH_SITE];
		incubation_schedule = new HashMap[LENGTH_SITE];
		recovery_schedule = new HashMap[LENGTH_SITE];

		for (int s = 0; s < LENGTH_SITE; s++) {
			currently_infectious[s] = new ArrayList<>();
			incubation_schedule[s] = new HashMap<>();
			recovery_schedule[s] = new HashMap<>();

		}

		trans_prob = new HashMap<>();
		sim_output = new HashMap<>();
		testing_schedule = new HashMap<>();
		risk_cat_map = new HashMap<>();
	}

	public HashMap<String, Object> getSim_output() {
		return sim_output;
	}

	public void scheduleNextTest(Integer personId, int lastTestTime) {
		// TODO: Check method
		int genderType = getGenderType(personId);

		double[][] testRate = ((double[][][]) runnable_fields[RUNNABLE_FIELD_TESTING_BY_RISK_CATEGORIES])[genderType];
		if (testRate != null) {
			int riskCat = Math.max(0, getRiskCategories(personId, genderType));
			double[] testRateByCat = testRate[riskCat];
			int divder = (testRateByCat.length - 1) / 2;

			double p = RNG.nextDouble();
			int pI = Arrays.binarySearch(testRateByCat, 0, divder, p);
			if (pI < 0) {
				pI = ~pI;
			}

			double testGapTime = (testRateByCat[divder + pI] + testRateByCat[divder + pI + 1]) / 2;
			testGapTime *= 1 + RNG.nextGaussian() / 10;

			ArrayList<Integer> testEnt = testing_schedule.get((int) Math.round(lastTestTime + (testGapTime)));
			if(testEnt == null) {
				testEnt = new ArrayList<>();
				testing_schedule.put(personId, testEnt);
			}
			testEnt.add(personId);		

		}

	}

	private int getRiskCategories(Integer personId, int genderType) {

		if (risk_cat_map.containsKey(personId)) {
			return risk_cat_map.get(personId);
		} else {

			int riskCat = -1;

			double[] riskCatList = ((double[][]) runnable_fields[RUNNABLE_FIELD_RISK_CATEGORIES_BY_CASUAL_PARTNERS])[genderType];
			if (riskCatList != null) {
				int numCasual = 0;
				Set<Integer[]> edges = BASE_CONTACT_MAP.edgesOf(personId);
				for (Integer[] e : edges) {
					if (e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME] >= firstSeedTime
							&& e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME] < firstSeedTime + NUM_TIME_STEPS
							&& e[Population_Bridging.CONTACT_MAP_EDGE_DURATION] == 1) {
						numCasual++;
					}
				}
				double numCasual1Year = ((double) AbstractIndividualInterface.ONE_YEAR_INT) * numCasual
						/ NUM_TIME_STEPS;
				riskCat = Arrays.binarySearch(riskCatList, numCasual1Year);
			}

			risk_cat_map.put(personId, riskCat);

			return riskCat;
		}

	}

	public int addInfected(Integer infectedId, int site, int firstContactTime, int recoveredAt) {
		int key = Collections.binarySearch(currently_infectious[site], infectedId);
		if (key < 0) {
			currently_infectious[site].add(~key, infectedId);
			firstSeedTime = Math.min(firstSeedTime, firstContactTime);

			// Recovery
			ArrayList<Integer> sch = recovery_schedule[site].get(recoveredAt);
			if (sch == null) {
				sch = new ArrayList<>();
				recovery_schedule[site].put(recoveredAt, sch);
			}
			sch.add(infectedId);

			// Transmission probability
			if (!trans_prob.containsKey(infectedId)) {
				double[][] trans = new double[LENGTH_SITE][LENGTH_SITE];
				int gender = getGenderType(infectedId);
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

	public ContactMap getTransmissionMap() {
		return transmissionMap;
	}

	@Override
	public void run() {

		int startTime = firstSeedTime;

		int[][] seedInfected = new int[currently_infectious.length][];

		// Store initially infected and set initial start time
		for (int site = 0; site < currently_infectious.length; site++) {
			ArrayList<Integer> currently_infectious_by_site = currently_infectious[site];
			seedInfected[site] = new int[currently_infectious_by_site.size()];
			int c = 0;
			for (Integer infectious : currently_infectious_by_site) {
				seedInfected[site][c] = infectious;
				c++;
			}
		}

		if (startTime < Integer.MAX_VALUE) {

			StringBuilder seedInfectedStr = new StringBuilder();
			for (int site = 0; site < seedInfected.length; site++) {
				for (int i = 0; i < seedInfected[site].length; i++) {
					seedInfectedStr.append(site);
					seedInfectedStr.append(',');
					seedInfectedStr.append(seedInfected[site][i]);
					seedInfectedStr.append('\n');
				}
			}

			ContactMap cMap;
			try {
				cMap = ContactMap.ContactMapFromFullString(BASE_CONTACT_MAP.toFullString());
			} catch (IOException e1) {
				cMap = BASE_CONTACT_MAP;
				e1.printStackTrace(System.err);
			}

			HashSet<Integer[]> removeEdges = new HashSet<>();

			for (int currentTime = startTime; currentTime < startTime + NUM_TIME_STEPS; currentTime++) {

				for (int site_src = 0; site_src < LENGTH_SITE; site_src++) {
					// Update infectious
					ArrayList<Integer> becomeInfectiousToday = incubation_schedule[site_src].remove(currentTime);

					if (becomeInfectiousToday != null) {
						for (Integer i : becomeInfectiousToday) {
							int recoveredAt = (int) Math.round(infectious_period[site_src].sample()) + currentTime;
							addInfected(i, site_src, currentTime, recoveredAt);
						}
					}

					// Update recovery
					ArrayList<Integer> recoveredToday = recovery_schedule[site_src].remove(currentTime);
					if (recoveredToday != null) {
						for (Integer i : recoveredToday) {
							removeInfected(i, site_src);
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
													boolean tranmitted = actType != -1;

													if (tranmitted) {

														double actProb = ((double[][][]) runnable_fields[RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ])[actType][g_s][g_t];
														double transProb = trans[site_src][site_target];
														tranmitted &= actProb > 0;
														if (tranmitted) {
															tranmitted &= RNG.nextDouble() < (actProb * transProb);
														}
													}

													if (tranmitted) {
														Integer incubation_end_at = currentTime
																+ (int) incubation_period[site_target].sample();
														ArrayList<Integer> ent = incubation_schedule[site_target]
																.get(incubation_end_at);
														if (ent == null) {
															ent = new ArrayList<>();
															incubation_schedule[site_target].put(incubation_end_at,
																	ent);
														}
														ent.add(partner);

														if (!transmissionMap.containsVertex(infectious)) {
															transmissionMap.addVertex(infectious);
														}
														if (!transmissionMap.containsVertex(partner)) {
															transmissionMap.addVertex(partner);
														}

														Integer[] existEdge = transmissionMap.getEdge(infectious,
																partner);

														if (existEdge == null) {
															existEdge = new Integer[LENGTH_TRANSMAP_EDGE];
															existEdge[TRANSMAP_EDGE_INFECTIOUS] = infectious;
															existEdge[TRANSMAP_EDGE_SUSCEPTIBLE] = partner;
															existEdge[TRANSMAP_EDGE_START_TIME] = currentTime;
															existEdge[TRANSMAP_EDGE_ACT_INVOLVED] = 1 << actType;
															transmissionMap.addEdge(infectious, partner, existEdge);
														} else {
															existEdge[TRANSMAP_EDGE_ACT_INVOLVED] |= 1 << actType;
														}

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
				
				// TODO: Testing

				for (Integer[] e : removeEdges) {
					cMap.removeEdge(e);
				}

			}
			// End of simulations

			sim_output.put(SIM_OUTPUT_TRANMISSION_MAP, transmissionMap);
			Set<ContactMap> clustersSet = transmissionMap.getContactCluster();

			sim_output.put(SIM_OUTPUT_CLUSTERS, clustersSet);

			if (runnableId != null) {
				System.out.println(String.format("Thread <%s> completed.", runnableId));
			}

			// Display clusters as CSV

			ContactMap[] clusters = clustersSet.toArray(new ContactMap[clustersSet.size()]);

			Arrays.sort(clusters, new Comparator<ContactMap>() {
				@Override
				public int compare(ContactMap o1, ContactMap o2) {
					return Integer.compare(o1.vertexSet().size(), o2.vertexSet().size());
				}
			});

			File clusterExport = new File(baseDir, String.format(DIRNAME_FORMAT_TRANSMISSION_CMAP, this.seed));
			clusterExport.mkdirs();

			File printFile;
			PrintWriter expWri;

			try {

				printFile = new File(clusterExport, String.format(FILENAME_FORMAT_INDEX_CASE_LIST, this.seed));

				expWri = new PrintWriter(printFile);
				expWri.println(seedInfectedStr.toString());
				expWri.close();

				for (int cI = 0; cI < clusters.length; cI++) {
					ContactMap c = clusters[cI];

					printFile = new File(clusterExport,
							String.format(FILENAME_FORMAT_TRANSMISSION_CMAP, Long.toString(this.seed), cI));

					expWri = new PrintWriter(printFile);
					expWri.println(c.toFullString());
					expWri.close();

				}
			} catch (IOException ex) {
				ex.printStackTrace(System.err);
				System.out.println("Index case:");
				System.out.println(seedInfectedStr.toString());

				for (int cI = 0; cI < clusters.length; cI++) {
					ContactMap c = clusters[cI];
					System.out.println(String.format("Transmission map <%d, %d>", this.seed, cI));
					System.out.println(c.toFullString());
					System.out.println();
				}

			}

		}

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
