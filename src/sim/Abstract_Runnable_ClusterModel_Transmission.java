package sim;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;

import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;

public abstract class Abstract_Runnable_ClusterModel_Transmission extends Abstract_Runnable_ClusterModel {

	public static final int ACT_INDEX_GENITAL = 0;
	public static final int ACT_INDEX_ANAL = ACT_INDEX_GENITAL + 1;
	public static final int ACT_INDEX_FELLATIO = ACT_INDEX_ANAL + 1;
	public static final int ACT_INDEX_RIMMING = ACT_INDEX_FELLATIO + 1;
	public static final int ACT_INDEX_KISSING = ACT_INDEX_RIMMING + 1;
	public static final int LENGTH_ACT = ACT_INDEX_KISSING + 1;

	public static final int SITE_VAGINA = 0;
	public static final int SITE_PENIS = SITE_VAGINA + 1;
	public static final int SITE_RECTUM = SITE_PENIS + 1;
	public static final int SITE_OROPHARYNX = SITE_RECTUM + 1;
	public static final int LENGTH_SITE = SITE_OROPHARYNX + 1;

	public static final int RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ = 0;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE = RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ + 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD = RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD = RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SYM_RATE = RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD + 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS = RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES = RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM = RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_DX_TEST_ACCURACY = RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
			+ 1;

	protected final int[] cUMULATIVE_POP_COMPOSITION;
	protected final ContactMap bASE_CONTACT_MAP;
	protected final int nUM_TIME_STEPS_PER_SNAP;
	protected final int nUM_SNAP;
	protected final long sIM_SEED;
	protected final long cMAP_SEED;

	protected PrintStream print_progress = null;
	protected int simSetting = 1;
	protected RandomGenerator RNG;
	protected ArrayList<Integer[]> edges_list;
	protected HashMap<Integer, HashMap<Integer, String>> propSwitch_map;
	protected transient HashMap<Integer, Integer> risk_cat_map;
	protected transient int firstSeedTime = Integer.MAX_VALUE;
	protected transient HashMap<String, Object> sim_output = null;

	public Abstract_Runnable_ClusterModel_Transmission(long cMap_seed, long sim_seed, int[] pop_composition,
			ContactMap base_cMap, int numTimeStepsPerSnap, int numSnap) {
		super();
		this.cUMULATIVE_POP_COMPOSITION = new int[pop_composition.length];
		int offset = 0;

		for (int g = 0; g < this.cUMULATIVE_POP_COMPOSITION.length; g++) {
			this.cUMULATIVE_POP_COMPOSITION[g] = offset + pop_composition[g];
			offset += pop_composition[g];
		}

		this.bASE_CONTACT_MAP = base_cMap;
		this.nUM_TIME_STEPS_PER_SNAP = numTimeStepsPerSnap;
		this.nUM_SNAP = numSnap;
		this.cMAP_SEED = cMap_seed;
		this.sIM_SEED = sim_seed;

		RNG = new MersenneTwisterRandomGenerator(sIM_SEED);
	}
	
	
	public abstract void initialse();	
	public abstract void allocateSeedInfection(int[][] num_infected_count, int time);
	public abstract int addInfectious(Integer infectedId, int infId, int site, int infectious_time, int recoveredAt);

	public long getSim_seed() {
		return sIM_SEED;
	}

	public long getcMap_seed() {
		return cMAP_SEED;
	}

	public void setPrint_progress(PrintStream print_progess) {
		this.print_progress = print_progess;
	}

	public int getSimSetting() {
		return simSetting;
	}

	public void setSimSetting(int simSetting) {
		this.simSetting = simSetting;
	}

	protected final AbstractRealDistribution generateNonDistribution(double[] input) {
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

	protected AbstractRealDistribution generateGammaDistribution(double[] input) {

		if (input[1] != 0) {
			// For Gamma distribution
			// GammaDistribution(RandomGenerator rng, double shape, double scale)
			// shape = mean / scale i.e. mean / (var / mean)
			// scale = var / mean
			double[] res = new double[2];
			double var = input[1] * input[1];
			// scale
			res[1] = var / input[0];
			// shape
			res[0] = input[0] / res[1];
			return new GammaDistribution(RNG, res[0], res[1]);
		} else {
			return generateNonDistribution(input);

		}
	}

	protected AbstractRealDistribution generateBetaDistribution(double[] input) {
		if (input[1] != 0) {

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
		} else {
			return generateNonDistribution(input);
		}
	}

	public void setEdges_list(ArrayList<Integer[]> edges_list) {
		this.edges_list = edges_list;
	}

	public void setPropSwitch_map(HashMap<Integer, HashMap<Integer, String>> propSwitch_map) {
		this.propSwitch_map = propSwitch_map;
	}


	public void fillRiskCatMap(ArrayList<Number[]> prealloactedRiskGrpArr) {
		if (prealloactedRiskGrpArr != null) {
			for (Number[] preAllocRisk : prealloactedRiskGrpArr) {
				if (risk_cat_map == null) {
					risk_cat_map = new HashMap<>();
				}
	
				risk_cat_map.put(
						(Integer) preAllocRisk[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_PID],
						(Integer) preAllocRisk[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP]);
			}
		}
	}

}
