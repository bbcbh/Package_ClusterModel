package sim;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import relationship.ContactMap;

public class Runnable_ClusterModel implements Runnable {

	public static final Object[] DEFAULT_BRIDGING_RUNNABLE_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_VALID_RANGE
			new int[] { 30, 30 + AbstractIndividualInterface.ONE_YEAR_INT },

	};

	public static final int RUNNABLE_FIELD_CONTACT_MAP_VALID_RANGE = 0;
	public static final int LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD = RUNNABLE_FIELD_CONTACT_MAP_VALID_RANGE + 1;

	private Population_Bridging population;
	private int numSnaps;
	private int snapFreq;
	private Object[] runnable_fields = new Object[LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD];
	private ContactMap[] gen_cMap = null;
	private boolean printStatus = false;

	public Runnable_ClusterModel() {
		super();
		for (int i = 0; i < DEFAULT_BRIDGING_RUNNABLE_FIELDS.length; i++) {
			runnable_fields[i] = DEFAULT_BRIDGING_RUNNABLE_FIELDS[i];
		}
	}

	public Object[] getRunnable_fields() {
		return runnable_fields;
	}

	public void setNumSnaps(int numSnaps) {
		this.numSnaps = numSnaps;
	}

	public ContactMap[] getGen_cMap() {
		return gen_cMap;
	}

	public void setSnapFreq(int snapFreq) {
		this.snapFreq = snapFreq;
	}

	public void setPrintStatus(boolean printStatus) {
		this.printStatus = printStatus;
	}

	public Population_Bridging getPopulation() {
		return population;
	}

	public void setPopulation(Population_Bridging population) {
		this.population = population;
	}

	@Override
	public void run() {
		
		if(printStatus) {
			population.setPrintStatus(printStatus);
		}
		population.initialise();
		if (printStatus) {
			System.out.println(population.printCurrentPartnershipStatus());
		}

		int[] contactMapValidRange = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_VALID_RANGE];

		for (int s = 0; s < numSnaps; s++) {
			for (int f = 0; f < snapFreq; f++) {
				if (population.getGlobalTime() + 1 == contactMapValidRange[0]) {
					gen_cMap = (ContactMap[]) population.getFields()[Population_Bridging.FIELD_CONTACT_MAP];
					for (int i = 0; i < gen_cMap.length; i++) {
						gen_cMap[i] = new ContactMap();
					}
				} else if (population.getGlobalTime() + 1 == contactMapValidRange[1]) {
					ContactMap[] cMapR = (ContactMap[]) population.getFields()[Population_Bridging.FIELD_CONTACT_MAP];
					for (int i = 0; i < cMapR.length; i++) {
						cMapR[i] = null;
					}
				}
				population.advanceTimeStep(1);

			}
			if (printStatus) {
				System.out.println(population.printCurrentPartnershipStatus());
			}

		}

	}

}
