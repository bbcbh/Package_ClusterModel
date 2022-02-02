package population;

import java.util.Arrays;

import person.AbstractIndividualInterface;
import population.person.Person_Bridging_Pop;
import relationship.ContactMap;
import relationship.RelationshipMap;
import relationship.SingleRelationship;

public class Bridging_Population extends AbstractFieldsArrayPopulation {

	/**
	 * 
	 */
	private static final long serialVersionUID = 4789823681005424649L;

	public final Object[] DEFAULT_BRIDGING_POP_FIELDS = {
			// FIELD_POP_COMPOSITION
			// int[] {NUM_FEMALE, NUM_MALE, NUM_MSMO, NUM_MSMW}
			new int[] { 5000, 5000, 200, 50 },
			// FIELD_CONTACT_MAP		
			new ContactMap[3],

	};

	public static final int FIELD_POP_COMPOSITION = AbstractFieldsArrayPopulation.LENGTH_FIELDS;
	public static final int FIELD_CONTACT_MAP = FIELD_POP_COMPOSITION + 1;
	public static final int LENGTH_FIELDS_BRIDGING_POP = FIELD_CONTACT_MAP + 1;

	public static final int CONTACT_MAP_ALL = 0;
	public static final int CONTACT_MAP_HETRO = CONTACT_MAP_ALL + 1;
	public static final int CONTACT_MAP_MSM = CONTACT_MAP_HETRO + 1;

	public Bridging_Population(long seed) {
		setSeed(seed);
		setRNG(new random.MersenneTwisterRandomGenerator(seed));

		Object[] newFields = Arrays.copyOf(fields, LENGTH_FIELDS_BRIDGING_POP);
		for (int i = AbstractFieldsArrayPopulation.LENGTH_FIELDS; i < newFields.length; i++) {
			newFields[i] = DEFAULT_BRIDGING_POP_FIELDS[i - AbstractFieldsArrayPopulation.LENGTH_FIELDS];
		}

		super.setFields(newFields);
	}

	@Override
	protected SingleRelationship formRelationship(AbstractIndividualInterface[] pair, RelationshipMap relMap, int d,
			int mapType) {

		ContactMap cMapSpec = null;
		Integer[] link = new Integer[pair.length];
		SingleRelationship rel;
		ContactMap cMapAll = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL];

		// 0b00 = Female, 0b01 = Hetro_male, Ob10 = MSMO, 0b11 = MSMW
		int contactType = ((Person_Bridging_Pop) pair[0]).getGenderType()
				& ((Person_Bridging_Pop) pair[1]).getGenderType();

		if (contactType == 0) { // Hetro sex involving female
			cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_HETRO];
		} else if (contactType == 0b10) { // Involve MSM
			cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_MSM];
		}

		for (int i = 0; i < pair.length; i++) {
			if (!relMap.containsVertex(pair[i].getId())) {
				relMap.addVertex(pair[i].getId());
			}
			if (!cMapAll.containsVertex(pair[i].getId())) {
				cMapAll.addVertex(pair[i].getId());
			}
			if (!cMapSpec.containsVertex(pair[i].getId())) {
				cMapSpec.addVertex(pair[i].getId());
			}
			link[i] = pair[i].getId();
		}

		Arrays.sort(link);

		rel = new SingleRelationship(link);
		if (relMap.addEdge(link[0], link[1], rel)) {

			Integer[] c = Arrays.copyOf(link, link.length);
			if (!cMapAll.containsEdge(c)) {
				cMapAll.addEdge(c[0], c[1], c);
			}
			if (!cMapSpec.containsEdge(c)) {
				cMapSpec.addEdge(c[0], c[1], c);
			}

			return rel;
		} else {
			return null;
		}
	}

	@Override
	public void initialise() {
		int[] popSizes = (int[]) getFields()[FIELD_POP_COMPOSITION];
		int popSizeTotal = 0;
		for(int i = 0; i < popSizes.length; i++) {
			popSizeTotal += popSizes[i];
		}
		
		AbstractIndividualInterface[] pop = new AbstractIndividualInterface[popSizeTotal];		
		int popPt = 0;
		for(int i = 0; i < popSizes.length; i++) {
			for(int g = 0; g < popSizes[i]; g++) {								
				pop[popPt] = new Person_Bridging_Pop(popPt+1, i, 
						18*AbstractIndividualInterface.ONE_YEAR_INT, getGlobalTime(), // Initial age - might not be used.
						3); // 3 sites
			}
			
		}
		
		
		// TODO Auto-generated method stub

	}

	@Override
	public void advanceTimeStep(int deltaT) {
		// TODO Auto-generated method stub

	}

}
