package population.availability;

import java.util.Arrays;

import availability.AbstractAvailability;
import person.AbstractIndividualInterface;
import random.RandomGenerator;
import util.ArrayUtilsRandomGenerator;

public class Availability_Bridging_Population extends AbstractAvailability {

	/**
	 * 
	 */
	private static final long serialVersionUID = 8014779566813031704L;

	protected AbstractIndividualInterface[][] available;
	protected AbstractIndividualInterface[][] pairing = null;
	protected boolean isBipartitieMapping = true;

	public static final String BIPARTITE_MAPPING = "BIPARTITE_MAPPING";

	public Availability_Bridging_Population(RandomGenerator RNG) {
		super(RNG);

	}

	@Override
	public int generatePairing() {
		int numPairs = 0;
		if (available != null) {
			numPairs = isBipartitieMapping ? Integer.MAX_VALUE : 0;
			for (int i = 0; i < available.length; i++) {
				if (isBipartitieMapping) {
					numPairs = Math.min(numPairs, available[i].length);
				} else {
					numPairs += available[i].length;
				}
			}

		}

		AbstractIndividualInterface[][] candidate = new AbstractIndividualInterface[2][];

		if (isBipartitieMapping) {
			for (int i = 0; i < candidate.length; i++) {
				if (available[i].length == numPairs) {
					candidate[i] = Arrays.copyOf(candidate[i], candidate[i].length);
				} else {
					candidate[i] = ArrayUtilsRandomGenerator.randomSelect(available[i], numPairs, getRNG());
				}

			}
		} else {
			numPairs = numPairs / 2;
			candidate[0] = ArrayUtilsRandomGenerator.randomSelect(available[0], 2 * numPairs, getRNG());
			candidate[1] = Arrays.copyOfRange(candidate[0], numPairs, candidate[0].length);
			candidate[0] = Arrays.copyOfRange(candidate[0], 0, numPairs);
		}

		for (int i = 0; i < candidate.length; i++) {
			ArrayUtilsRandomGenerator.shuffleArray(candidate[i], getRNG());
		}

		pairing = new AbstractIndividualInterface[numPairs][2];
		for (int i = 0; i < pairing.length; i++) {
			for (int p = 0; p < pairing.length; p++) {
				pairing[i][p] = candidate[p][i];
			}

		}

		return numPairs;
	}

	@Override
	public AbstractIndividualInterface[][] getPairing() {
		return pairing;
	}

	@Override
	public boolean removeMemberAvailability(AbstractIndividualInterface p) {
		if (available == null) {
			return false;
		}

		// Brutal method at this stage
		for (int a = 0; a < available.length; a++) {
			for (int m = 0; m < available[a].length; m++) {
				if (p.getId() == available[a][m].getId()) {
					System.arraycopy(available[a], m, available[a], m + 1, available[a].length - m);
					available[a] = Arrays.copyOf(available[a], available[a].length - 1);
					return true;
				}
			}
		}

		return false;
	}

	@Override
	public boolean memberAvailable(AbstractIndividualInterface p) {
		if (available == null) {
			return false;
		}
		// Brutal method at this stage
		for (int a = 0; a < available.length; a++) {
			for (int m = 0; m < available[a].length; m++) {
				if (p.getId() == available[a][m].getId()) {
					return true;
				}
			}
		}
		return false;
	}

	@Override
	public void setAvailablePopulation(AbstractIndividualInterface[][] available) {
		this.available = available;

	}

	@Override
	public boolean setParameter(String id, Object value) {
		if (BIPARTITE_MAPPING.equals(id)) {
			isBipartitieMapping = (Boolean) value;
			return true;
		} else {
			return false;
		}

	}

	@Override
	public Object getParameter(String id) {
		if (BIPARTITE_MAPPING.equals(id)) {
			return isBipartitieMapping;
		} else {
			return null;
		}
	}

}
