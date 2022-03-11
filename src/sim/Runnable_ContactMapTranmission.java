package sim;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;

public class Runnable_ContactMapTranmission implements Runnable {

	ContactMap transmissionMap;
	ArrayList<Integer> currently_infectious;
	RandomGenerator RNG;

	// FIELD_POP_COMPOSITION
	// int[] {NUM_FEMALE, NUM_MALE, NUM_MSMO, NUM_MSMW}
	final int[] POP_COMPOSITION;
	final ContactMap BASE_CONTACT_MAP;
	final int NUM_TIME_STEPS;
	final long seed;

	double[] INCUBATION_RANGE = new double[] { 3, 3 }; // 3 - 5 days
	double[] INFECTIOUS_DIST_PARAM = new double[] { 15 * 7, 5 * 7 };
	double[] TRANS_MF = new double[] { 0.2, 0.05 };
	double[] TRANS_FM = new double[] { 0.4, 0.10 };
	double[] TRANS_MM = new double[] { 0.4, 0.10 }; // Assumption

	protected transient AbstractRealDistribution infectious_period;
	protected transient AbstractRealDistribution tran_MF;
	protected transient AbstractRealDistribution tran_FM;
	protected transient AbstractRealDistribution tran_MM;

	protected transient HashMap<Integer, ArrayList<Integer>> incubation_schedule;
	protected transient HashMap<Integer, ArrayList<Integer>> recovery_schedule;
	protected transient HashMap<Integer, double[]> trans_prob;

	public Runnable_ContactMapTranmission(long seed, int[] POP_COMPOSITION, ContactMap BASE_CONTACT_MAP,
			int NUM_TIME_STEPS) {
		super();
		this.POP_COMPOSITION = POP_COMPOSITION;
		this.BASE_CONTACT_MAP = BASE_CONTACT_MAP;
		this.NUM_TIME_STEPS = NUM_TIME_STEPS;
		this.seed = seed;

	}

	public void initialse() {
		transmissionMap = new ContactMap();
		currently_infectious = new ArrayList<>();
		RNG = new MersenneTwisterRandomGenerator(seed);

		infectious_period = generateGammaDistribution(INFECTIOUS_DIST_PARAM);
		tran_FM = generateBetaDistribution(TRANS_FM);
		tran_MF =  generateBetaDistribution(TRANS_MF);
		tran_MM = generateBetaDistribution(TRANS_MM);

		incubation_schedule = new HashMap<>();
		recovery_schedule = new HashMap<>();
		trans_prob = new HashMap<>();

	}
	
	protected GammaDistribution generateGammaDistribution(double[] input) {
		return  new GammaDistribution(RNG, input[0], 1.0f / input[1]);
	}
	
	protected BetaDistribution generateBetaDistribution(double[] input) {
        // For Beta distribution, 
        // alpha = mean*(mean*(1-mean)/variance - 1)
        // beta = (1-mean)*(mean*(1-mean)/variance - 1)
        double[] res = new double[2];
        double var = input[1] * input[1];
        double rP = input[0] * (1 - input[0]) / var - 1;
        //alpha
        res[0] = rP * input[0];
        //beta
        res[1] = rP * (1 - input[0]);        
        return new BetaDistribution(RNG, res[0], res[1]);       
    }

	public int addInfected(Integer infectedId, int recoveredAt) {
		int key = Collections.binarySearch(currently_infectious, infectedId);
		if (key < 0) {
			currently_infectious.add(~key, infectedId);

			// Recovery

			ArrayList<Integer> sch = recovery_schedule.get(recoveredAt);
			if (sch == null) {
				sch = new ArrayList<>();
				recovery_schedule.put(recoveredAt, sch);
			}
			sch.add(infectedId);

			// Trans prob
			if (!trans_prob.containsKey(infectedId)) {
				double[] trans = new double[] { tran_FM.sample(), tran_MF.sample(), tran_MM.sample() };
				trans_prob.put(infectedId, trans);
			}

		}
		return key;
	}

	public int removeInfected(Integer infectedId) {
		int key = Collections.binarySearch(currently_infectious, infectedId);
		if (key >= 0) {
			currently_infectious.remove(key);
		}
		return key;
	}

	public int getGenderType(Integer personId) {
		for (int i = 0; i < Population_Bridging.LENGTH_GENDER; i++) {
			if (personId < POP_COMPOSITION[i]) {
				return i;
			}
		}
		return -1;
	}

	public ContactMap getTransmissionMap() {
		return transmissionMap;
	}

	@Override
	public void run() {

		int startTime = Integer.MAX_VALUE;

		// Set initial start time
		for (Integer infectious : currently_infectious) {
			if (BASE_CONTACT_MAP.containsVertex(infectious)) {
				Set<Integer[]> edges = BASE_CONTACT_MAP.edgesOf(infectious);
				for (Integer[] e : edges) {
					startTime = Math.min(startTime, e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME]);
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

				// Update infectious
				ArrayList<Integer> becomeInfectiousToday = incubation_schedule.remove(currentTime);

				if (becomeInfectiousToday != null) {
					for (Integer i : becomeInfectiousToday) {
						int recoveredAt = (int) Math.round(infectious_period.sample()) + currentTime;
						addInfected(i, recoveredAt);
					}
				}
				ArrayList<Integer> recoveredToday = recovery_schedule.remove(currentTime);
				if (recoveredToday != null) {
					for (Integer i : recoveredToday) {
						removeInfected(i);
					}
				}

				for (Integer infectious : currently_infectious) {
					if (cMap.containsVertex(infectious)) {
						Set<Integer[]> edges = cMap.edgesOf(infectious);

						for (Integer[] e : edges) {														
							
							if (currentTime <=  e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME] 
									&&  currentTime< (e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME] 
											+  e[Population_Bridging.CONTACT_MAP_EDGE_DURATION])) {

								int partner = e[Population_Bridging.CONTACT_MAP_EDGE_P1].equals(infectious)
										? e[Population_Bridging.CONTACT_MAP_EDGE_P2]
										: e[Population_Bridging.CONTACT_MAP_EDGE_P1];

								if (Collections.binarySearch(currently_infectious, partner) < 0) {

									boolean tranmitted = false;
									double[] trans = trans_prob.get(infectious);

									if (getGenderType(infectious) == Population_Bridging.GENDER_HETRO_FEMALE) {
										tranmitted = RNG.nextDouble() < trans[0];
									} else if (getGenderType(partner) == Population_Bridging.GENDER_HETRO_FEMALE) {
										tranmitted = RNG.nextDouble() < trans[1];
									} else {
										tranmitted = RNG.nextDouble() < trans[2];
									}

									if (tranmitted) {
										Integer[] tranmission_edge = new Integer[] { infectious, partner, currentTime };
										if (!transmissionMap.containsVertex(infectious)) {
											transmissionMap.addVertex(infectious);
										}
										if (!transmissionMap.containsVertex(partner)) {
											transmissionMap.addVertex(partner);
										}
										transmissionMap.addEdge(infectious, partner, tranmission_edge);
									}

								}

							} else {
								removeEdges.add(e);
							}
						}
					}
				}

				for (Integer[] e : removeEdges) {
					cMap.removeEdge(e);
				}

			}
		}

	}

}
