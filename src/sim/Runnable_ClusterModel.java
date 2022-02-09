package sim;

import population.Population_Bridging;

public class Runnable_ClusterModel implements Runnable {

	public static final int LENGTH_RUNNABLE_CLUSTER_MODEL_FIELD = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP;

	Population_Bridging population;
	int numSnaps;
	int snapFreq;

	public void setNumSnaps(int numSnaps) {
		this.numSnaps = numSnaps;
	}

	public void setSnapFreq(int snapFreq) {
		this.snapFreq = snapFreq;
	}

	public Population_Bridging getPopulation() {
		return population;
	}

	public void setPopulation(Population_Bridging population) {
		this.population = population;
	}

	@Override
	public void run() {
		// TODO Runnable method
		population.initialise();		
		
		for (int s = 0; s < numSnaps; s++) {
			for (int f = 0; f < snapFreq; f++) {
				
				population.advanceTimeStep(1);
			}

		}

	}

}
