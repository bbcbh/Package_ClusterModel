package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Properties;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PropValUtils;

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
	public static final int RUNNABLE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES = RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
			+ 1;
	protected static final int RUNNABLE_OFFSET = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
			+ +Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
			+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD;

	public static int getGenderType(Integer personId, int[] cumul_pop_comp) {
		int index = Arrays.binarySearch(cumul_pop_comp, personId);

		if (index < 0) {
			return ~index;
		} else {
			return index;
		}
	}

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
	public static final int PROPSWITCH_MAP_KEPT_EDGE_BY_DURATION_INDEX = -1;

	protected transient HashMap<Integer, Integer> risk_cat_map;
	protected transient int firstSeedTime = Integer.MAX_VALUE;
	protected transient HashMap<String, Object> sim_output = null;

	protected float[][] cMap_Kept_Edge_By_Duration_Setting = null;
	private static final int cMAP_KEPT_EDGE_BY_DURATION_DUR_LIST = 0;
	private static final int cMAP_KEPT_EDGE_BY_DURATION_PROB_LIST = cMAP_KEPT_EDGE_BY_DURATION_DUR_LIST + 1;

	protected static final String popCompositionKey = Simulation_ClusterModelTransmission.POP_PROP_INIT_PREFIX
			+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);

	public Abstract_Runnable_ClusterModel_Transmission(long cMap_seed, long sim_seed, ContactMap base_cMap,
			Properties prop) {

		this(cMap_seed, sim_seed,
				(int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey), int[].class), base_cMap,
				Integer.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])),
				Integer.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])));
	}

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

	public abstract int addInfectious(Integer infectedId, int infId, int site, int stage_id, int infectious_time,
			int recoveredAt);

	public abstract void scheduleNextTest(Integer personId, int lastTestTime);

	public abstract void refreshField(int fieldId, boolean clearAll);

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

	protected static final AbstractRealDistribution generateNonDistribution(RandomGenerator RNG, double[] input) {
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

	protected static AbstractRealDistribution generateGammaDistribution(RandomGenerator RNG, double[] input) {
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
			return generateNonDistribution(RNG, input);

		}
	}

	protected static AbstractRealDistribution generateBetaDistribution(RandomGenerator RNG, double[] input) {
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
			return generateNonDistribution(RNG, input);
		}
	}

	protected static AbstractRealDistribution generateUniformDistribution(RandomGenerator RNG, double[] input) {
		if (input[1] != 0) {
			return new UniformRealDistribution(RNG, input[0], input[1]);
		} else {
			return generateNonDistribution(RNG, input);
		}
	}

	public void setEdges_list(ArrayList<Integer[]> edges_list) {
		this.edges_list = edges_list;
	}

	public void setPropSwitch_map(HashMap<Integer, HashMap<Integer, String>> propSwitch_map) {
		this.propSwitch_map = propSwitch_map;
	}

	protected void loadExtraPropSwitchSetting(Integer index, String entry) {
		switch (index.intValue()) {
		case PROPSWITCH_MAP_KEPT_EDGE_BY_DURATION_INDEX:
			cMap_Kept_Edge_By_Duration_Setting = entry.trim().length() > 0
					? (float[][]) PropValUtils.propStrToObject(entry, float[][].class)
					: null;
			break;
		default:
			System.err.printf("Warning: Index of %d in simSpecificSwitch.prop not supported.\n", index);
		}
	}

	public int getGenderType(Integer personId) {
		return getGenderType(personId, cUMULATIVE_POP_COMPOSITION);
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

	protected Integer[][] getEdgesArrayFromBaseConctactMap() {
		if (edges_list == null) {

			try {
				edges_list = generateMapEdgeArray(bASE_CONTACT_MAP).call();
			} catch (Exception e) {
				e.printStackTrace(System.err);
				System.err.println("Error in generating edge list from BASE_CONTACT_MAP. Exiting...");
				edges_list = new ArrayList<>();
				System.exit(-1);

			}
		}
		Integer[][] edges_array = edges_list.toArray(new Integer[edges_list.size()][]);
		return edges_array;
	}

	protected int initaliseCMap(ContactMap cMap, Integer[][] edges_array, int edges_array_pt, int startTime,
			HashMap<Integer, ArrayList<Integer[]>> removeEdges) {
		ArrayList<Integer[]> toRemove;

		// Skip invalid edges
		while (edges_array_pt < edges_array.length
				&& edges_array[edges_array_pt][Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] < startTime) {
			Integer[] edge = edges_array[edges_array_pt];
			int edge_start_time = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
			int expireAt = edge_start_time + edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];

			if (expireAt > startTime) {
				toRemove = removeEdges.get(expireAt);
				if (toRemove == null) {
					toRemove = new ArrayList<>();
					removeEdges.put(expireAt, toRemove);
				}
				toRemove.add(edge);

				for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
						Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
					if (!cMap.containsVertex(edge[index])) {
						cMap.addVertex(edge[index]);
					}
				}
				addPartnership(cMap, edge);
			}

			edges_array_pt++;
		}
		return edges_array_pt;
	}

	protected void setPreAllocatedRiskFromFile() {
		File pre_allocate_risk_file = new File(baseDir,
				String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, cMAP_SEED));

		if (pre_allocate_risk_file.isFile()) {
			try {
				BufferedReader reader = new BufferedReader(new FileReader(pre_allocate_risk_file));
				String line;

				if (risk_cat_map == null) {
					risk_cat_map = new HashMap<>();

				}

				while ((line = reader.readLine()) != null) {
					String[] lineSp = line.split(",");
					risk_cat_map.put(
							Integer.parseInt(
									lineSp[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_PID]),
							Integer.parseInt(
									lineSp[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP]));
				}

				reader.close();
			} catch (Exception e) {
				e.printStackTrace(System.err);
			}
		}
	}

	protected int updateCMap(ContactMap cMap, int currentTime, Integer[][] edges_array, int edges_array_pt,
			HashMap<Integer, ArrayList<Integer[]>> edgesToRemove) {
		ArrayList<Integer[]> toRemove;
		// Remove expired edges
		toRemove = edgesToRemove.get(currentTime);
		if (toRemove != null) {
			for (Integer[] edge : toRemove) {
				cMap.removeEdge(edge);
			}
		}

		// Add new edges and update removal schedule
		while (edges_array_pt < edges_array.length
				&& edges_array[edges_array_pt][Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] <= currentTime) {
			boolean addEdge = cMap_Kept_Edge_By_Duration_Setting == null;
			Integer[] edge = edges_array[edges_array_pt];

			if (cMap_Kept_Edge_By_Duration_Setting != null) {				
				int index = Arrays.binarySearch(cMap_Kept_Edge_By_Duration_Setting[cMAP_KEPT_EDGE_BY_DURATION_DUR_LIST],
						Math.max(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION], 1));
				if (index < 0) {
					index = ~index;
				}
				addEdge = index >= cMap_Kept_Edge_By_Duration_Setting[cMAP_KEPT_EDGE_BY_DURATION_PROB_LIST].length;
				if (!addEdge) {
					float keptEdgeProb = cMap_Kept_Edge_By_Duration_Setting[cMAP_KEPT_EDGE_BY_DURATION_PROB_LIST][index];
					if(keptEdgeProb > 0) {
						addEdge = keptEdgeProb >= 1 ? true : RNG.nextFloat() < keptEdgeProb;
					}
				}
			}
			if (addEdge) {
				Integer expireAt = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
						+ Math.max(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION], 1);

				toRemove = edgesToRemove.get(expireAt);
				if (toRemove == null) {
					toRemove = new ArrayList<>();
					edgesToRemove.put(expireAt, toRemove);
				}
				toRemove.add(edge);

				for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
						Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
					if (!cMap.containsVertex(edge[index])) {
						cMap.addVertex(edge[index]);
					}
				}

				addPartnership(cMap, edge);
			}

			edges_array_pt++;
		}
		return edges_array_pt;
	}

	protected void addPartnership(ContactMap cMap, Integer[] edge) {
		cMap.addEdge(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
				edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], edge);
	}

	public HashMap<String, Object> getSim_output() {
		return sim_output;
	}

	public abstract ArrayList<Integer> loadOptParameter(String[] parameter_settings, double[] point,
			int[][] seedInfectNum, boolean display_only);

}
