package sim;

import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;

public class Runnable_ClusterModel_ContactMap_Generation_MultiMap
		extends Abstract_Runnable_ClusterModel_ContactMap_Generation {

	// Runnable fields
	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP = LENGTH_RUNNABLE_MAP_GEN_FIELD;
	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_GRP_CANDIDATE = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP
			+ 1;

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_GRP_CANDIDATE
			+ 1;

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST
			+ 1;

	public static final int LENGTH_RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP
			+ 1;

	public static final Object[] DEFAULT_RUNNABLE_MAP_GEN_MULTIMAP_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP
			new int[] { 1000 },
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_GRP_CANDIDATE
			// 1 << Grp Number
			new int[] { 1 },
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST
			// double[GRP_NUMBER]{dist}
			new double[][] { new double[] { 18 * AbstractIndividualInterface.ONE_YEAR_INT,
					34 * AbstractIndividualInterface.ONE_YEAR_INT, 0, 1 }, },
			
			//RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP
			// double[GRP_NUMBER][PROB_, 
			new double[][] {
						new double[] {}
			},
	
	};

	// Population Index
	private static final int POP_INDEX_GRP = 0;
	private static final int POP_INDEX_ENTER_POP_AGE = POP_INDEX_GRP + 1;
	private static final int POP_INDEX_ENTER_POP_AT = POP_INDEX_ENTER_POP_AGE + 1;
	private static final int POP_INDEX_EXIT_POP_AT = POP_INDEX_ENTER_POP_AT + 1;
	private static final int LENGTH_POP_ENTRIES = POP_INDEX_EXIT_POP_AT + 1;

	public static final String MAPFILE_FORMAT = "ContactMap_Type_%d_%d.csv"; // Type, Seed,
	public static final String POPSTAT_FORMAT = "POP_STAT_%d"; // Seed

	private RandomGenerator RNG;

	public Runnable_ClusterModel_ContactMap_Generation_MultiMap(long mapSeed) {
		super(mapSeed);

		Object[] newFields = Arrays.copyOf(runnable_fields, LENGTH_RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP);

		for (int i = LENGTH_RUNNABLE_MAP_GEN_FIELD; i < LENGTH_RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP; i++) {
			newFields[i] = DEFAULT_RUNNABLE_MAP_GEN_MULTIMAP_FIELDS[i - LENGTH_RUNNABLE_MAP_GEN_FIELD];
		}

		runnable_fields = newFields;
		RNG = new MersenneTwisterRandomGenerator(mapSeed);

	}

	@Override
	public void run() {
		HashMap<Integer, Object[]> population = new HashMap<>();
		int nextId = 1;
		int[] numInGrp = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP];

		UnivariateInterpolator polator = new LinearInterpolator();

		// Initialise pop
		for (int g = 0; g < numInGrp.length; g++) {
			double[] ageDist = (double[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST];

			UnivariateFunction ageDistFunc = polator.interpolate(Arrays.copyOf(ageDist, ageDist.length / 2),
					Arrays.copyOfRange(ageDist, ageDist.length / 2, ageDist.length));

			for (int i = 0; i < numInGrp[g]; i++) {
				Object[] newPerson = new Object[LENGTH_POP_ENTRIES];
				newPerson[POP_INDEX_GRP] = g;
				newPerson[POP_INDEX_ENTER_POP_AT] = 1;
				newPerson[POP_INDEX_ENTER_POP_AGE] = (int) Math.round(ageDistFunc.value(RNG.nextDouble()));
				newPerson[POP_INDEX_EXIT_POP_AT] = ((int) ageDist[ageDist.length / 2]
						- (int) newPerson[POP_INDEX_ENTER_POP_AGE]) + (int) newPerson[POP_INDEX_ENTER_POP_AT];
				// Undefined

				population.put(nextId, newPerson);

				nextId++;
			}

		}

	}

}
