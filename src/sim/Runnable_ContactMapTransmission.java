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

import population.Population_Bridging;
import population.person.Person_Bridging_Pop;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;

public class Runnable_ContactMapTransmission implements Runnable {

	public static final int ACT_INDEX_GENITAL = 0;
	public static final int ACT_INDEX_ANAL = ACT_INDEX_GENITAL + 1;
	public static final int ACT_INDEX_FELLATIO = ACT_INDEX_ANAL + 1;

	public static final int SITE_VAGINA = 0;
	public static final int SITE_PENIS = SITE_VAGINA + 1;
	public static final int SITE_RECTUM = SITE_PENIS + 1;
	public static final int SITE_OROPHARYNX = SITE_RECTUM + 1;
	public static final int SITE_LENGTH = SITE_OROPHARYNX + 1;

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

	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ = 0;
	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_TRANSMISSION_RATE = RUNNABLE_FIELD_TRANSMISSION_MAP_ACT_FREQ
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_INFECTIOUS_PERIOD = RUNNABLE_FIELD_TRANSMISSION_MAP_TRANSMISSION_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_MAP_INCUBATION_PERIOD = RUNNABLE_FIELD_TRANSMISSION_MAP_INFECTIOUS_PERIOD
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

	};

	// FIELD_POP_COMPOSITION
	// int[] {NUM_FEMALE, NUM_MALE, NUM_MSMO, NUM_MSMW}
	final int[] cumulative_pop_composition;
	final ContactMap BASE_CONTACT_MAP;
	final int NUM_TIME_STEPS;
	final long seed;

	private ContactMap transmissionMap;
	private RandomGenerator RNG;

	protected transient RealDistribution[][] tranmissionMatrix = new RealDistribution[SITE_LENGTH][SITE_LENGTH];
	protected transient RealDistribution[] infectious_period = new RealDistribution[SITE_LENGTH];
	protected transient RealDistribution[] incubation_period = new RealDistribution[SITE_LENGTH];

	protected transient ArrayList<Integer>[] currently_infectious;
	protected transient HashMap<Integer, ArrayList<Integer>>[] incubation_schedule;
	protected transient HashMap<Integer, ArrayList<Integer>>[] recovery_schedule;
	protected transient HashMap<Integer, double[][]> trans_prob;

	private transient HashMap<String, Object> sim_output = null;

	public static final String SIM_OUTPUT_TRANMISSION_MAP = "SIM_OUTPUT_TRANMISSION_MAP";
	public static final String SIM_OUTPUT_CLUSTERS = "SIM_OUTPUT_CLUSTERS";

	public String runnableId = null;

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

	public String getRunnableId() {
		return runnableId;
	}

	public void setRunnableId(String runnableId) {
		this.runnableId = runnableId;
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
		currently_infectious = new ArrayList[SITE_LENGTH];
		incubation_schedule = new HashMap[SITE_LENGTH];
		recovery_schedule = new HashMap[SITE_LENGTH];

		for (int s = 0; s < SITE_LENGTH; s++) {
			currently_infectious[s] = new ArrayList<>();
			incubation_schedule[s] = new HashMap<>();
			recovery_schedule[s] = new HashMap<>();
		}

		trans_prob = new HashMap<>();
		sim_output = new HashMap<>();
	}

	public HashMap<String, Object> getSim_output() {
		return sim_output;
	}

	public int addInfected(Integer infectedId, int site, int recoveredAt) {
		int key = Collections.binarySearch(currently_infectious[site], infectedId);
		if (key < 0) {
			currently_infectious[site].add(~key, infectedId);

			// Recovery
			ArrayList<Integer> sch = recovery_schedule[site].get(recoveredAt);
			if (sch == null) {
				sch = new ArrayList<>();
				recovery_schedule[site].put(recoveredAt, sch);
			}
			sch.add(infectedId);

			// Transmission probability
			if (!trans_prob.containsKey(infectedId)) {
				double[][] trans = new double[SITE_LENGTH][SITE_LENGTH];
				int gender = getGenderType(infectedId);
				for (int sf = 0; sf < SITE_LENGTH; sf++) {
					boolean sample = (gender == Person_Bridging_Pop.GENDER_TYPE_FEMALE) ? sf != SITE_PENIS
							: sf != SITE_VAGINA;
					if (sample) {
						for (int st = 0; st < SITE_LENGTH; st++) {
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
		for (int s = 0; s < SITE_LENGTH; s++) {
			removeInfected(infectedId, s);
		}
	}

	public int getGenderType(Integer personId) {

		int index = Arrays.binarySearch(cumulative_pop_composition, personId);

		if (index < 0) {
			return ~index;
		} else {
			return index + 1;
		}

	}

	public ContactMap getTransmissionMap() {
		return transmissionMap;
	}

	@Override
	public void run() {

		int startTime = Integer.MAX_VALUE;

		// Set initial start time
		for (ArrayList<Integer> currently_infectious_by_site : currently_infectious) {
			for (Integer infectious : currently_infectious_by_site) {
				if (BASE_CONTACT_MAP.containsVertex(infectious)) {
					Set<Integer[]> edges = BASE_CONTACT_MAP.edgesOf(infectious);
					for (Integer[] e : edges) {
						startTime = Math.min(startTime, e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME]);
					}
				}
			}
		}

		if (startTime < Integer.MAX_VALUE) {

			ContactMap cMap;
			try {
				cMap = ContactMap.ContactMapFromFullString(BASE_CONTACT_MAP.toFullString());
			} catch (IOException e1) {
				cMap = BASE_CONTACT_MAP;
				e1.printStackTrace(System.err);
			}

			HashSet<Integer[]> removeEdges = new HashSet<>();

			for (int currentTime = startTime; currentTime < startTime + NUM_TIME_STEPS; currentTime++) {

				for (int site_src = 0; site_src < SITE_LENGTH; site_src++) {
					// Update infectious
					ArrayList<Integer> becomeInfectiousToday = incubation_schedule[site_src].remove(currentTime);

					if (becomeInfectiousToday != null) {
						for (Integer i : becomeInfectiousToday) {
							int recoveredAt = (int) Math.round(infectious_period[site_src].sample()) + currentTime;
							addInfected(i, site_src, recoveredAt);
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
									int durationIndex = startIndex+1;

									if (currentTime >= e[startIndex]
											&& currentTime < (e[startIndex]
													+ e[durationIndex])) {

										int partner = e[Population_Bridging.CONTACT_MAP_EDGE_P1].equals(infectious)
												? e[Population_Bridging.CONTACT_MAP_EDGE_P2]
												: e[Population_Bridging.CONTACT_MAP_EDGE_P1];

										for (int site_target = 0; site_target < SITE_LENGTH; site_target++) {
											if (trans[site_src][site_target] != 0) {
												// Transmission is possible
												if (Collections.binarySearch(currently_infectious[site_target],
														partner) < 0) {

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
														int g_s = getGenderType(infectious);
														int g_t = getGenderType(partner);
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
																infectious);

														if (existEdge == null) {
															existEdge = new Integer[] { infectious, partner,
																	currentTime, 1 << actType };
															transmissionMap.addEdge(infectious, partner, existEdge);
														} else {
															existEdge[existEdge.length - 1] |= 1 << actType;
														}

													}

												}

											}

										}

									} else if (currentTime >= (e[startIndex] + e[durationIndex]) 
											&& startIndex+2 > e.length) {
										removeEdges.add(e);
									}
									startIndex += 2;																	
								}
							}
						}
					}
				}

				for (Integer[] e : removeEdges) {
					cMap.removeEdge(e);
				}

			}
			// End of simulations

			sim_output.put(SIM_OUTPUT_TRANMISSION_MAP, transmissionMap);
			Set<ContactMap> clusters = transmissionMap.getContactCluster();

			sim_output.put(SIM_OUTPUT_CLUSTERS, clusters);

			if (runnableId != null) {
				System.out.println(String.format("Thread <%s> completed.", runnableId));
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
