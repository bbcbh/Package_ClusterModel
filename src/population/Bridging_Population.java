package population;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import availability.AbstractAvailability;
import person.AbstractIndividualInterface;
import population.availability.Bridging_Population_Availability;
import population.person.Person_Bridging_Pop;
import relationship.ContactMap;
import relationship.RelationshipMap;
import relationship.SingleRelationship;
import util.ArrayUtilsRandomGenerator;

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
			// FIELD_PARTNER_TYPE_PROB
			// float[GENDER]{CUMUL_REG_ONLY, CUMUL_CAS_ONLY}
			// Default:
			// Lee E, Mao L, Broady T, et al. Gay Community Periodic Survey: Melbourne 2018.
			// Sydney: Centre for Social Research in Health, UNSW Sydney, 2018.
			new float[][] { { 1, 1 }, { 1, 1 }, { 0.33f, 0.33f + 0.26f }, { 0.33f, 0.33f + 0.26f } },
			// FIELD_MEAN_CASUAL_PARTNER_NUM
			// Currently bases on a PossionDistribution
			new float[] { 0, 0, 10, 10 },
			// FIELD_MEAN_MSM_REG_PARTNERSHIP_DURATION
			4.0 * AbstractIndividualInterface.ONE_YEAR,

	};

	public static final int FIELD_POP_COMPOSITION = AbstractFieldsArrayPopulation.LENGTH_FIELDS;
	public static final int FIELD_CONTACT_MAP = FIELD_POP_COMPOSITION + 1;
	public static final int FIELD_PARTNER_TYPE_PROB = FIELD_CONTACT_MAP + 1;
	public static final int FIELD_MEAN_CASUAL_PARTNER_NUM = FIELD_PARTNER_TYPE_PROB + 1;
	public static final int FIELD_MEAN_MSM_REG_PARTNERSHIP_DURATION = FIELD_PARTNER_TYPE_PROB + 1;
	public static final int LENGTH_FIELDS_BRIDGING_POP = FIELD_MEAN_MSM_REG_PARTNERSHIP_DURATION + 1;

	public static final int CONTACT_MAP_ALL = 0;
	public static final int CONTACT_MAP_HETRO = CONTACT_MAP_ALL + 1;
	public static final int CONTACT_MAP_MSM = CONTACT_MAP_HETRO + 1;

	public static final int RELMAP_HETRO = 0;
	public static final int RELMAP_MSM = 1;
	public static final int RELMAP_TOTAL = 2;

	private AbstractIntegerDistribution msm_reg_partner_duration;

	public Bridging_Population(long seed) {
		setSeed(seed);
		setRNG(new random.MersenneTwisterRandomGenerator(seed));

		Object[] newFields = Arrays.copyOf(fields, LENGTH_FIELDS_BRIDGING_POP);
		for (int i = AbstractFieldsArrayPopulation.LENGTH_FIELDS; i < newFields.length; i++) {
			newFields[i] = DEFAULT_BRIDGING_POP_FIELDS[i - AbstractFieldsArrayPopulation.LENGTH_FIELDS];
		}

		super.setFields(newFields);
	}

	public Class<? extends Object> getFieldsClass(int fieldIndex) {
		return this.getFields()[fieldIndex].getClass();
	}

	@Override
	protected SingleRelationship formRelationship(AbstractIndividualInterface[] pair, RelationshipMap relMap,
			int duration, int mapType) {

		ContactMap cMapSpec = null;
		Integer[] link = new Integer[pair.length];
		SingleRelationship rel;
		ContactMap cMapAll = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL];

		if (mapType == 0) { // Hetro sex involving female
			cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_HETRO];
		} else { // Involve MSM
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

			rel.setDurations(duration);

			return rel;
		} else {
			return null;
		}
	}

	protected boolean seekingCasualToday(Person_Bridging_Pop person) {
		int maxCasual = (int) person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_6_MONTHS));
		boolean seekCasual = false;

		if (maxCasual > 0) {
			int numCasual = person.getNumCasualInRecord();
			seekCasual = getRNG().nextInt(6 * 30) < (maxCasual - numCasual);
		}
		return seekCasual;
	}

	protected void formCasualPartnership(ArrayList<Person_Bridging_Pop> casualListMSM) {
		ContactMap[] cMaps = new ContactMap[] { ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL],
				((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_MSM] };

		Person_Bridging_Pop[] casualMSMArr = casualListMSM.toArray(new Person_Bridging_Pop[casualListMSM.size()]);
		ArrayUtilsRandomGenerator.shuffleArray(casualMSMArr, getRNG());

		for (int i = 0; (i + 1) < casualMSMArr.length; i += 2) {
			Person_Bridging_Pop[] pair = new Person_Bridging_Pop[] { casualMSMArr[i], casualMSMArr[i + 1] };

			// For consistency
			if (pair[0].getId() > pair[1].getId()) {
				pair[0] = casualMSMArr[i + 1];
				pair[1] = casualMSMArr[i];
			}

			// Form casual pairing if there are not already in a regular partnership
			if ((this.getRelMap()[RELMAP_MSM]).containsEdge(pair[0].getId(), pair[1].getId())) {
				for (int c = 0; c < cMaps.length; c++) {
					for (int p = 0; p < pair.length; p++) {
						if (!cMaps[c].containsVertex(pair[p].getId())) {
							cMaps[c].addVertex(pair[p].getId());
						}
					}

					cMaps[c].addEdge(pair[0].getId(), pair[1].getId(),
							new Integer[] { pair[0].getId(), pair[1].getId() });
				}

			}

		}

	}

	@Override
	public void initialise() {

		// Initialise relationship map
		RelationshipMap[] relMaps = new RelationshipMap[RELMAP_TOTAL];
		relMaps[RELMAP_HETRO] = new RelationshipMap();
		relMaps[RELMAP_MSM] = new RelationshipMap();
		this.setRelMap(relMaps);

		// Initialise population
		int[] popSizes = (int[]) getFields()[FIELD_POP_COMPOSITION];
		int popSizeTotal = 0;
		for (int i = 0; i < popSizes.length; i++) {
			popSizeTotal += popSizes[i];
		}

		AbstractIndividualInterface[] pop = new AbstractIndividualInterface[popSizeTotal];

		int popPt = 0;
		float[][] partner_type = (float[][]) getFields()[FIELD_PARTNER_TYPE_PROB];
		float[] mean_casual_partners = (float[]) getFields()[FIELD_MEAN_CASUAL_PARTNER_NUM];

		ArrayList<Person_Bridging_Pop> casualListMSM = new ArrayList<>();

		for (int i = 0; i < popSizes.length; i++) {

			// Shared among same gender
			float[] p_type_prob = partner_type[i];

			AbstractIntegerDistribution dist = null;

			for (int g = 0; g < popSizes[i]; g++) {
				pop[popPt] = new Person_Bridging_Pop(popPt + 1, i, 18 * AbstractIndividualInterface.ONE_YEAR_INT,
						getGlobalTime(), // Initial age - might not be used.
						3); // 3 sites

				Person_Bridging_Pop person = (Person_Bridging_Pop) pop[popPt];

				// FIELD_PARTNER_TYPE_PROB
				// float[GENDER]{CUMUL_REG_ONLY, CUMUL_CAS_ONLY}

				// Initially set all partner type
				person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_6_MONTHS), 1);
				person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER), 1);

				boolean canCasual = true;
				float prob = getRNG().nextFloat();

				for (int p = 0; p < p_type_prob.length; p++) {

					if (prob < p_type_prob[p]) {
						if (p == 0) { // CUMUL_REG_ONLY
							person.setParameter(
									Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_6_MONTHS), 0);
							canCasual = false;
						}
						if (p == 1) { // CUMUL_CAS_ONLY
							person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER), 0);

						}
					}
				}
				// Set max. number of casual in 6 months
				if (canCasual) {
					if (dist == null) {
						dist = new PoissonDistribution(mean_casual_partners[i]);
					}
					person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_6_MONTHS),
							dist.sample());

					if (seekingCasualToday(person) && (0b10 & person.getGenderType()) != 0) {
						casualListMSM.add(person);
					}
				}

				if (i != Person_Bridging_Pop.GENDER_TYPE_MSMO) {
					relMaps[RELMAP_HETRO].addAvailablePerson(pop[popPt]);
				}
				if (i == Person_Bridging_Pop.GENDER_TYPE_MSMO || i == Person_Bridging_Pop.GENDER_TYPE_MSMW) {
					relMaps[RELMAP_MSM].addAvailablePerson(pop[popPt]);
				}
			}
		}

		this.setPop(pop);

		// Set availability

		AbstractAvailability[] avail = new AbstractAvailability[RELMAP_TOTAL];

		avail[RELMAP_HETRO] = new Bridging_Population_Availability(getRNG());
		avail[RELMAP_HETRO].setParameter(Bridging_Population_Availability.BIPARTITE_MAPPING, true);
		avail[RELMAP_MSM] = new Bridging_Population_Availability(getRNG());
		avail[RELMAP_MSM].setParameter(Bridging_Population_Availability.BIPARTITE_MAPPING, false);
		
		// TODO: Form regular partnership - Hetro

		// Form regular partnership - MSM

		msm_reg_partner_duration = new PoissonDistribution(
				(double) getFields()[FIELD_MEAN_MSM_REG_PARTNERSHIP_DURATION]);
		avail[RELMAP_MSM].setAvailablePopulation(relMaps[RELMAP_MSM].getPersonsAvailable(null, this.getLocalData()));

		int numMSM_Pair = avail[RELMAP_MSM].generatePairing();
		AbstractIndividualInterface[][] pairs = avail[RELMAP_MSM].getPairing();
		for (int p = 0; p < numMSM_Pair; p++) {
			formRelationship(pairs[p], relMaps[RELMAP_MSM], msm_reg_partner_duration.sample(), 1);

		}

		// Form casual partnership - MSM
		
		formCasualPartnership(casualListMSM);

	}

	@Override
	public void advanceTimeStep(int deltaT) {
		incrementTime(deltaT);

		// TODO Auto-generated method stub

		updatePairs();

	}

	protected void updatePairs() {

		for (int map = 0; map < getRelMap().length; map++) {

			RelationshipMap relMap = getRelMap()[map];

			// Update existing
			SingleRelationship[] relArr = relMap.getRelationshipArray();

			if (relMap.edgeSet().size() != relArr.length) {
				relArr = relMap.edgeSet().toArray(new SingleRelationship[relMap.edgeSet().size()]);
			}

			for (SingleRelationship relArr1 : relArr) {
				if (relArr1.incrementTime(1) <= 0) {
					relMap.removeEdge(relArr1);
				}
			}

			// No need to be sorted
			getAvailability()[map].setAvailablePopulation(getRelMap()[map].getPersonsAvailable(null, getLocalData()));

			// Generate new pairing
			int pairNum = getAvailability()[map].generatePairing();
			AbstractIndividualInterface[][] pairs = getAvailability()[map].getPairing();
			for (int pairId = 0; pairId < pairNum; pairId++) {
				int duration = -1;
				if (map == RELMAP_MSM) {
					duration = msm_reg_partner_duration.sample();
				}
				formRelationship(pairs[pairId], getRelMap()[map], duration, map);

			}
		}
	}

}
