package sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

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
	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP+1
			+ 1;

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST
			+ 1;

	public static final int LENGTH_RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP
			+ 1;

	public static final Object[] DEFAULT_RUNNABLE_MAP_GEN_MULTIMAP_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP
			new int[] { 1000 },
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST
			// double[GRP_NUMBER]{dist}
			new double[][] { new double[] { 18 * AbstractIndividualInterface.ONE_YEAR_INT,
					34 * AbstractIndividualInterface.ONE_YEAR_INT, 0, 1 }, },

			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP		
			
			new double[][] { new double[] { 1, 1, 1, 1, 0.7, 0.7, -2.8 }, },

	};

	// Population Index
	private static final int POP_INDEX_GRP = 0;
	private static final int POP_INDEX_ENTER_POP_AGE = POP_INDEX_GRP + 1;
	private static final int POP_INDEX_ENTER_POP_AT = POP_INDEX_ENTER_POP_AGE + 1;
	private static final int POP_INDEX_EXIT_POP_AT = POP_INDEX_ENTER_POP_AT + 1;
	private static final int LENGTH_POP_ENTRIES = POP_INDEX_EXIT_POP_AT + 1;
	
	// MAPSETTING
	// If DUR_FREQ < 0, then it is number of one partnership to form within snapshot
	// else it is the duration of partnership
	private static final int MAPSETTING_MAP_TYPE = 0;
	private static final int MAPSETTING_SNAP_FREQ = MAPSETTING_MAP_TYPE+1;
	private static final int MAPSETTING_GRP_INDEX_P1 = MAPSETTING_SNAP_FREQ+1;
	private static final int MAPSETTING_GRP_INDEX_P2 = MAPSETTING_GRP_INDEX_P1+1;
	private static final int MAPSETTING_PROB_HAS_PARTNERSHIP_P1 =  MAPSETTING_GRP_INDEX_P2+1;
	private static final int MAPSETTING_PROB_HAS_PARTNERSHIP_P2 =  MAPSETTING_PROB_HAS_PARTNERSHIP_P1 + 1;
	private static final int MAPSETTING_DUR_FREQ = MAPSETTING_GRP_INDEX_P2+1;
	
	

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
		HashMap<Integer, ArrayList<Integer>> active_in_pop = new HashMap<>();

		int nextId = 1;
		int popTime = 0;
		int[] numInGrp = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP];
		double[][] ageDist = (double[][]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST];

		UnivariateInterpolator polator = new LinearInterpolator();

		// Reuse variables
		ArrayList<Integer> active_by_grp;

		// Initialise pop
		for (int g = 0; g < numInGrp.length; g++) {			
			UnivariateFunction ageDistFunc = polator.interpolate(Arrays.copyOf(ageDist[g], ageDist[g].length / 2),
					Arrays.copyOfRange(ageDist[g], ageDist[g].length / 2, ageDist[g].length));

			active_by_grp = new ArrayList<>();
			for (int i = 0; i < numInGrp[g]; i++) {
				Object[] newPerson = new Object[LENGTH_POP_ENTRIES];
				newPerson[POP_INDEX_GRP] = g;
				newPerson[POP_INDEX_ENTER_POP_AT] = 1;
				newPerson[POP_INDEX_ENTER_POP_AGE] = (int) Math.round(ageDistFunc.value(RNG.nextDouble()));
				newPerson[POP_INDEX_EXIT_POP_AT] = ((int) ageDist[g][ageDist[g].length / 2]
						- (int) newPerson[POP_INDEX_ENTER_POP_AGE]) + (int) newPerson[POP_INDEX_ENTER_POP_AT];
				population.put(nextId, newPerson);

				active_by_grp.add(nextId);
				nextId++;
			}
			active_in_pop.put(g, active_by_grp);

		}

		for (int snapC = 0; snapC < numSnaps; snapC++) {
			popTime += snapFreq;
			
			for (int g = 0; g < numInGrp.length; g++) {
				active_by_grp = active_in_pop.get(g);
				int numRemoved = 0;

				Iterator<Integer> iter = active_by_grp.iterator();
				while (iter.hasNext()) {
					int pId = iter.next();
					Object[] perStat = population.get(pId);
					// Remove age out person
					if ((int) perStat[POP_INDEX_EXIT_POP_AT] <= popTime) {
						iter.remove();
						numRemoved++;
					}
					// Add new person
					while(numRemoved > 0) {						
						Object[] newPerson = new Object[LENGTH_POP_ENTRIES];
						newPerson[POP_INDEX_GRP] = g;
						newPerson[POP_INDEX_ENTER_POP_AT] = popTime - RNG.nextInt(snapFreq) ;
						newPerson[POP_INDEX_ENTER_POP_AGE] = (int) ageDist[g][0];
						newPerson[POP_INDEX_EXIT_POP_AT] = ((int) ageDist[g][ageDist.length / 2]
								- (int) newPerson[POP_INDEX_ENTER_POP_AGE]) + (int) newPerson[POP_INDEX_ENTER_POP_AT];
						population.put(nextId, newPerson);
						active_by_grp.add(nextId);
						nextId++;													
						numRemoved--;
					}
				}																	
				
			}
			
			for(double[] map_setting :  (double[][]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP]) {				
				// TODO: Generate map
				int map_type = (int) map_setting[MAPSETTING_MAP_TYPE];
			}
			
			
			
		

		}

	}

}
