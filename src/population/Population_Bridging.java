package population;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import availability.AbstractAvailability;
import person.AbstractIndividualInterface;
import population.availability.Availability_Bridging_Population;
import population.person.Person_Bridging_Pop;
import relationship.ContactMap;
import relationship.RelationshipMap;
import relationship.SingleRelationship;
import util.ArrayUtilsRandomGenerator;

public class Population_Bridging extends AbstractFieldsArrayPopulation {

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
			// Hetro: ASHR
			// Rissel C, Badcock PB, Smith AMA, et al. Heterosexual experience and recent
			// heterosexual
			// encounters among Australian adults: the Second Australian Study of Health and
			// Relationships.
			// Sexual health 2014; 11:416-26.
			// MSM:
			// Lee E, Mao L, Broady T, et al. Gay Community Periodic Survey: Melbourne 2018.
			// Sydney: Centre for Social Research in Health, UNSW Sydney, 2018.
			new float[][] { { 1 - 0.027f, 1 - 0.027f }, { 1 - 0.053f, 1 - 0.053f }, { 0.33f, 0.33f + 0.26f },
					{ 0.33f, 0.33f + 0.26f } },
			// FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS
			// float[GENDER][MEAN_PARTNER_IN_12_MONTHS]
			// Default: Hetro ASHR2 (16 - 30) , MSM - ASHR2
			new float[] { (1.0f + 1.1f) / 2, (1.4f + 1.4f) / 2, 6.8f, 6.8f },
			// FIELD_MEAN_REG_PARTNERSHIP_DUR
			// float[]{HETRO, MSM}
			// Default:
			// Hetro (weight ave from ASHR)
			// Badcock PB, Smith AMA, Richters J, Rissel C, de Visser RO, Simpson JM, et al.
			// Characteristics of heterosexual regular relationships among a representative
			// sample of adults:
			// the Second Australian Study of Health and Relationships. Sexual health.
			// 2014;11(5):427-38.(see ACCEPt appendix)
			// MSM:
			// Prestage GP, Hudson J, Bradley J, et al. TOMS - Three or More Study.
			// Sydney: National Centre in HIV Epidemiology and Clinical Research, University
			// of New South Wales, 2008. (see MSM tech appendix)
			new float[] { 19.2f * AbstractIndividualInterface.ONE_YEAR_INT,
					4f * AbstractIndividualInterface.ONE_YEAR_INT },
			// FIELD_MEAN_CASUAL_PARTNER_IN_12_MONTHS
			// Currently bases on a PossionDistribution from MSM
			// Assumption for Hetro
			new float[] { 12f, 12f, 20.55f, 20.55f },

	};

	public static final int FIELD_POP_COMPOSITION = AbstractFieldsArrayPopulation.LENGTH_FIELDS;
	public static final int FIELD_CONTACT_MAP = FIELD_POP_COMPOSITION + 1;
	public static final int FIELD_PARTNER_TYPE_PROB = FIELD_CONTACT_MAP + 1;
	public static final int FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS = FIELD_PARTNER_TYPE_PROB + 1;
	public static final int FIELD_MEAN_REG_PARTNERSHIP_DUR = FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS + 1;
	public static final int FIELD_MEAN_NUM_CASUAL_PARTNER_IN_12_MONTHS = FIELD_MEAN_REG_PARTNERSHIP_DUR + 1;
	public static final int LENGTH_FIELDS_BRIDGING_POP = FIELD_MEAN_NUM_CASUAL_PARTNER_IN_12_MONTHS + 1;

	public static final int GENDER_HETRO_FEMALE = 0;
	public static final int GENDER_HETRO_MALE = GENDER_HETRO_FEMALE + 1;
	public static final int GENDER_MSMO = GENDER_HETRO_MALE + 1;
	public static final int GENDER_MSMW = GENDER_MSMO + 1;
	public static final int LENGTH_GENDER = GENDER_MSMW + 1;

	public static final int CONTACT_MAP_ALL = 0;
	public static final int CONTACT_MAP_HETRO = CONTACT_MAP_ALL + 1;
	public static final int CONTACT_MAP_MSM = CONTACT_MAP_HETRO + 1;

	public static final int RELMAP_HETRO = 0;
	public static final int RELMAP_MSM = 1;
	public static final int RELMAP_TOTAL = 2;

	@SuppressWarnings("unchecked")
	private transient ArrayList<Person_Bridging_Pop>[] canSeekRelPartners = new ArrayList[LENGTH_GENDER];
	@SuppressWarnings("unchecked")
	private transient ArrayList<Person_Bridging_Pop>[] hasRelPartners = new ArrayList[LENGTH_GENDER];
	@SuppressWarnings("unchecked")
	private transient ArrayList<Person_Bridging_Pop>[] canSeekCasPartners = new ArrayList[LENGTH_GENDER];

	private transient Person_Bridging_Pop[][] casualPartnerFormed = null;

	private AbstractIntegerDistribution[] regPartDuration = new AbstractIntegerDistribution[RELMAP_TOTAL];

	private boolean printStatus = false;

	public void setPrintStatus(boolean printStatus) {
		this.printStatus = printStatus;
	}

	public Population_Bridging(long seed) {
		setSeed(seed);
		setRNG(new random.MersenneTwisterRandomGenerator(seed));

		Object[] newFields = Arrays.copyOf(fields, LENGTH_FIELDS_BRIDGING_POP);
		for (int i = AbstractFieldsArrayPopulation.LENGTH_FIELDS; i < newFields.length; i++) {
			newFields[i] = DEFAULT_BRIDGING_POP_FIELDS[i - AbstractFieldsArrayPopulation.LENGTH_FIELDS];
		}

		super.setFields(newFields);

		for (int g = 0; g < LENGTH_GENDER; g++) {
			canSeekRelPartners[g] = new ArrayList<Person_Bridging_Pop>();
			hasRelPartners[g] = new ArrayList<Person_Bridging_Pop>();
			canSeekCasPartners[g] = new ArrayList<Person_Bridging_Pop>();
		}

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

			link[i] = pair[i].getId();
		}

		checkContactMaps(link, new ContactMap[] { cMapAll, cMapSpec });

		rel = new SingleRelationship(link);

		if (relMap.addEdge(link[0], link[1], rel)) {
			Integer[] c = Arrays.copyOf(link, link.length);

			rel.setDurations(duration);
			for (int p = 0; p < pair.length; p++) {
				((Person_Bridging_Pop) pair[p]).addRegularPartner(pair[(p + 1) % 2]);
			}

			return rel;
		} else {
			return null;
		}
	}

	protected boolean seekingRegularToday(Person_Bridging_Pop person) {
		int maxReg = (Integer) person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS));

		int numReg = person.getNumRegularInRecord();
		boolean seekRegular = maxReg > numReg;

		if (seekRegular) {
			int noPWin = AbstractIndividualInterface.ONE_YEAR_INT;
			if ((int) person.getField(Person_Bridging_Pop.FIELD_LAST_REGULAR_PARTNER_AT_AGE) > 0) {
				noPWin = Math.max(1, noPWin - ((int) person.getAge()
						- (int) person.getField(Person_Bridging_Pop.FIELD_LAST_REGULAR_PARTNER_AT_AGE)));
			}
			seekRegular = getRNG().nextInt(noPWin) < (maxReg - numReg);
		}

		return seekRegular;
	}

	protected boolean seekingCasualToday(Person_Bridging_Pop person) {
		int maxCasual = (int) person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));
		int numCasual = person.getNumCasualInRecord();
		boolean seekCasual = maxCasual > numCasual;

		if (seekCasual) {
			int noPWin = AbstractIndividualInterface.ONE_YEAR_INT;
			if ((int) person.getField(Person_Bridging_Pop.FIELD_LAST_REGULAR_PARTNER_AT_AGE) > 0) {
				noPWin = Math.max(1, noPWin - ((int) person.getAge()
						- (int) person.getField(Person_Bridging_Pop.FIELD_LAST_REGULAR_PARTNER_AT_AGE)));
			}
			seekCasual = getRNG().nextInt(noPWin) < (maxCasual - numCasual);
		}
		return seekCasual;
	}

	protected void formCasualPartnership(int[] pop_diff_num_partner_12_months) {

		casualPartnerFormed = new Person_Bridging_Pop[0][];
		int partPt = 0;
		ContactMap[] cMaps;

		Person_Bridging_Pop[][] casualCandidate = new Person_Bridging_Pop[LENGTH_GENDER][0];

		for (int g = 0; g < casualCandidate.length; g++) {
			if (pop_diff_num_partner_12_months[g] > 0) {
				casualCandidate[g] = canSeekCasPartners[g]
						.toArray(new Person_Bridging_Pop[canSeekCasPartners[g].size()]);
				if (casualCandidate[g].length > pop_diff_num_partner_12_months[g]) {
					casualCandidate[g] = ArrayUtilsRandomGenerator.randomSelect(casualCandidate[g],
							pop_diff_num_partner_12_months[g], getRNG());
				}

			}
		}

		// MSM
		if (casualCandidate[GENDER_MSMO].length + casualCandidate[GENDER_MSMW].length > 1) {

			Person_Bridging_Pop[] casualMSMArr = Arrays.copyOf(casualCandidate[GENDER_MSMO],
					casualCandidate[GENDER_MSMO].length + casualCandidate[GENDER_MSMW].length);

			System.arraycopy(casualCandidate[GENDER_MSMW], 0, casualMSMArr, casualCandidate[GENDER_MSMO].length,
					casualCandidate[GENDER_MSMW].length);

			int numMSMCasualToFormed = casualMSMArr.length / 2;

			if (numMSMCasualToFormed > 0) {
				cMaps = new ContactMap[] { ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL],
						((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_MSM] };

				ArrayUtilsRandomGenerator.shuffleArray(casualMSMArr, getRNG());

				casualPartnerFormed = Arrays.copyOf(casualPartnerFormed,
						casualPartnerFormed.length + numMSMCasualToFormed);

				for (int i = 0; (i + 1) < casualMSMArr.length && partPt < casualPartnerFormed.length; i += 2) {
					Person_Bridging_Pop[] pair = new Person_Bridging_Pop[] { casualMSMArr[i], casualMSMArr[i + 1] };

					// For consistency
					if (pair[0].getId() > pair[1].getId()) {
						pair[0] = casualMSMArr[i + 1];
						pair[1] = casualMSMArr[i];
					}

					// Form casual pairing if there are not already in a regular partnership
					if (!(this.getRelMap()[RELMAP_MSM]).containsEdge(pair[0].getId(), pair[1].getId())) {
						for (int p = 0; p < pair.length; p++) {
							pair[p].addCasualPartner(pair[(p + 1) % 2]);
						}
						casualPartnerFormed[partPt] = new Person_Bridging_Pop[] { pair[0], pair[1] };
						checkContactMaps(new Integer[] { pair[0].getId(), pair[1].getId() }, cMaps);
						partPt++;
					}
				}
			}
		}

		// Heterosexual
		if (casualCandidate[GENDER_HETRO_FEMALE].length > 0
				&& (casualCandidate[GENDER_HETRO_MALE].length + casualCandidate[GENDER_MSMW].length) > 0) {

			int numHetroCasualToFormed = Math.min(casualCandidate[GENDER_HETRO_FEMALE].length,
					casualCandidate[GENDER_HETRO_MALE].length + casualCandidate[GENDER_MSMW].length);

			if (numHetroCasualToFormed > 0) {

				cMaps = new ContactMap[] { ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL],
						((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_HETRO] };

				Person_Bridging_Pop[] casualFemale = casualCandidate[GENDER_HETRO_FEMALE];
				Person_Bridging_Pop[] casualMale = Arrays.copyOf(casualCandidate[GENDER_HETRO_MALE],
						casualCandidate[GENDER_HETRO_MALE].length + casualCandidate[GENDER_MSMW].length);

				System.arraycopy(casualCandidate[GENDER_MSMW], 0, casualMale, casualCandidate[GENDER_HETRO_MALE].length,
						casualCandidate[GENDER_MSMW].length);

				ArrayUtilsRandomGenerator.shuffleArray(casualFemale, getRNG());
				ArrayUtilsRandomGenerator.shuffleArray(casualMale, getRNG());

				casualPartnerFormed = Arrays.copyOf(casualPartnerFormed,
						casualPartnerFormed.length + numHetroCasualToFormed);

				int arrPt = 0;
				while (partPt < casualPartnerFormed.length
						&& arrPt < Math.min(casualMale.length, casualFemale.length)) {
					Person_Bridging_Pop[] pair = new Person_Bridging_Pop[] { casualMale[arrPt], casualFemale[arrPt] };

					arrPt++;
					if (!(this.getRelMap()[RELMAP_HETRO]).containsEdge(pair[0].getId(), pair[1].getId())) {
						for (int p = 0; p < pair.length; p++) {
							pair[p].addCasualPartner(pair[(p + 1) % 2]);
						}
						casualPartnerFormed[partPt] = new Person_Bridging_Pop[] { pair[0], pair[1] };
						checkContactMaps(new Integer[] { pair[0].getId(), pair[1].getId() }, cMaps);
						partPt++;
					}
				}

			}
		}

		if (partPt < casualPartnerFormed.length) {
			casualPartnerFormed = Arrays.copyOf(casualPartnerFormed, partPt);
		}
	}

	@Override
	public void initialise() {

		// Initialise relationship map
		RelationshipMap[] relMaps = new RelationshipMap[RELMAP_TOTAL];
		relMaps[RELMAP_HETRO] = new RelationshipMap();
		relMaps[RELMAP_MSM] = new RelationshipMap();
		this.setRelMap(relMaps);

		float[] meanRelDur = (float[]) getFields()[FIELD_MEAN_REG_PARTNERSHIP_DUR];
		regPartDuration[RELMAP_HETRO] = new PoissonDistribution(getRNG(), meanRelDur[RELMAP_HETRO],
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		regPartDuration[RELMAP_MSM] = new PoissonDistribution(getRNG(), meanRelDur[RELMAP_MSM],
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);

		// Initialise population
		int[] popSizes = (int[]) getFields()[FIELD_POP_COMPOSITION];
		int popSizeTotal = 0;
		for (int i = 0; i < popSizes.length; i++) {
			popSizeTotal += popSizes[i];
		}

		AbstractIndividualInterface[] pop = new AbstractIndividualInterface[popSizeTotal];

		int popPt = 0;
		float[][] partner_type = (float[][]) getFields()[FIELD_PARTNER_TYPE_PROB];
		float[] mean_casual_partners = (float[]) getFields()[FIELD_MEAN_NUM_CASUAL_PARTNER_IN_12_MONTHS];
		int[] population_num_partner_in_last_12_months = new int[LENGTH_GENDER];

		for (int genderGrpIndex = 0; genderGrpIndex < popSizes.length; genderGrpIndex++) {

			// Shared among same gender
			float[] p_type_prob = partner_type[genderGrpIndex];

			AbstractIntegerDistribution dist = null;

			for (int g = 0; g < popSizes[genderGrpIndex]; g++) {
				pop[popPt] = new Person_Bridging_Pop(popPt + 1, genderGrpIndex,
						18 * AbstractIndividualInterface.ONE_YEAR_INT, getGlobalTime(), // Initial age - might not be
																						// used.
						3); // 3 sites

				Person_Bridging_Pop person = (Person_Bridging_Pop) pop[popPt];

				// FIELD_PARTNER_TYPE_PROB
				// float[GENDER]{CUMUL_REG_ONLY, CUMUL_CAS_ONLY}

				// Initially set all partner type
				person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS), 1);
				person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS), 1);

				boolean canCasual = true;
				float prob = getRNG().nextFloat();

				int prob_pt = 0;
				while (prob_pt < p_type_prob.length && prob > p_type_prob[prob_pt]) {
					prob_pt++;
				}

				if (prob_pt == 0) { // CUMUL_REG_ONLY
					person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS), 0);
					canCasual = false;
				} else if (prob_pt == 1) { // CUMUL_CAS_ONLY
					person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS), 0);
				}

				// Set max. number of casual in 12 months
				if (canCasual) {
					if (dist == null) {
						dist = new PoissonDistribution(getRNG(), mean_casual_partners[genderGrpIndex],
								PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
					}
					person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS),
							Math.max(1, dist.sample()));

					if (seekingCasualToday(person)) {
						canSeekCasPartners[person.getGenderType()].add(person);
					}
				}

				if (seekingRegularToday(person)) {
					canSeekRelPartners[person.getGenderType()].add(person);
				}
				popPt++;
			}
		}

		this.setPop(pop);
		this.getFields()[FIELDS_NEXT_ID] = popPt;

		// Set availability

		AbstractAvailability[] avail = new AbstractAvailability[RELMAP_TOTAL];

		avail[RELMAP_HETRO] = new Availability_Bridging_Population(getRNG());
		avail[RELMAP_HETRO].setParameter(Availability_Bridging_Population.BIPARTITE_MAPPING, true);
		avail[RELMAP_MSM] = new Availability_Bridging_Population(getRNG());
		avail[RELMAP_MSM].setParameter(Availability_Bridging_Population.BIPARTITE_MAPPING, false);
		this.setAvailability(avail);

		formRegularCasualPartnerships(population_num_partner_in_last_12_months);

	}

	public void formRegularCasualPartnerships(int[] population_num_partner_in_last_12_months) {
		int[] pop_diff_num_partner_12_months = cal_pop_diff_num_partner(population_num_partner_in_last_12_months);
		int[] reg_part_formed = formRegularPartnership(pop_diff_num_partner_12_months);
		for (int g = 0; g < pop_diff_num_partner_12_months.length; g++) {
			pop_diff_num_partner_12_months[g] -= reg_part_formed[g];
		}

		formCasualPartnership(pop_diff_num_partner_12_months);
	}

	@Override
	public void advanceTimeStep(int deltaT) {
		incrementTime(deltaT);

		int[] population_num_partner_in_last_12_months = new int[LENGTH_GENDER];

		// Clear partnership status

		for (int g = 0; g < LENGTH_GENDER; g++) {
			canSeekRelPartners[g].clear();
			hasRelPartners[g].clear();
			canSeekCasPartners[g].clear();
		}

		for (AbstractIndividualInterface p : this.getPop()) {
			Person_Bridging_Pop person = (Person_Bridging_Pop) p;
			person.incrementTime(deltaT, getInfList());
			int genderType = person.getGenderType();
			int numReg12Months = person.getNumRegularInRecord();
			int numCas12Months = person.getNumCasualInRecord();

			population_num_partner_in_last_12_months[genderType] += numReg12Months + numCas12Months;

			if (seekingRegularToday(person)) {
				canSeekRelPartners[genderType].add(person);
			}

			RelationshipMap[] relMap = getRelMap();
			boolean hasRel;

			switch (genderType) {
			case Person_Bridging_Pop.GENDER_TYPE_MSMO:
				hasRel = relMap[RELMAP_MSM].degreeOf(person.getId()) > 0;
				break;
			case Person_Bridging_Pop.GENDER_TYPE_MSMW:
				hasRel = relMap[RELMAP_MSM].degreeOf(person.getId()) > 0
						|| relMap[RELMAP_HETRO].degreeOf(person.getId()) > 0;
				break;
			default:
				hasRel = relMap[RELMAP_HETRO].degreeOf(person.getId()) > 0;
			}
			if (hasRel) {
				hasRelPartners[genderType].add(person);
			}

			if (seekingCasualToday(person)) {
				canSeekCasPartners[person.getGenderType()].add(person);
			}

		}

		if (printStatus) {
			System.out.println(String.format("# partners in last 12 months = %s - Day %d",
					Arrays.toString(population_num_partner_in_last_12_months), getGlobalTime()));
		}
		updateRelRelationshipMap();
		formRegularCasualPartnerships(population_num_partner_in_last_12_months);

	}

	private int[] cal_pop_diff_num_partner(int[] population_num_partner_in_last_12_months) {
		float[] mean_target = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int[] numInGrp = (int[]) (getFields()[FIELD_POP_COMPOSITION]);
		int[] pop_diff = new int[LENGTH_GENDER];

		for (int g = 0; g < mean_target.length; g++) {
			pop_diff[g] = Math.round(mean_target[g] * numInGrp[g] - population_num_partner_in_last_12_months[g]);
		}
		return pop_diff;
	}

	public String printCurrentPartnershipStatus() {

		StringWriter strWri = new StringWriter();
		PrintWriter pWri = new PrintWriter(strWri);

		pWri.println(String.format("# relationship - Day %d", getGlobalTime()));

		for (int r = 0; r < getRelMap().length; r++) {
			pWri.println(String.format("RelMap #%d: %d", r, getRelMap()[r].edgeSet().size()));
		}

		if (casualPartnerFormed != null) {
			pWri.println(String.format("# casual partnership - Day %d", getGlobalTime()));
			for (int c = 0; c < casualPartnerFormed.length; c++) {
				pWri.println(String.format("#%d: <%d (%d), %d (%d)>", c, casualPartnerFormed[c][0].getId(),
						casualPartnerFormed[c][0].getGenderType(), casualPartnerFormed[c][1].getId(),
						casualPartnerFormed[c][1].getGenderType()));
			}
		}

		return strWri.toString();
	}

	protected int[] formRegularPartnership(int[] pop_diff_num_partner_12_months) {

		AbstractIndividualInterface[][] candidates = new AbstractIndividualInterface[LENGTH_GENDER][];

		int[] partformed = new int[LENGTH_GENDER];

		for (int g = 0; g < pop_diff_num_partner_12_months.length; g++) {

			if (pop_diff_num_partner_12_months[g] > 0) {
				candidates[g] = canSeekRelPartners[g]
						.toArray(new AbstractIndividualInterface[canSeekRelPartners[g].size()]);
				if (pop_diff_num_partner_12_months[g] < canSeekRelPartners[g].size()) {
					candidates[g] = ArrayUtilsRandomGenerator.randomSelect(candidates[g],
							pop_diff_num_partner_12_months[g], getRNG());
				}

			} else if (pop_diff_num_partner_12_months[g] <= 0) {

				candidates[g] = hasRelPartners[g].toArray(new AbstractIndividualInterface[hasRelPartners[g].size()]);

				if (-pop_diff_num_partner_12_months[g] < candidates[g].length) {
					candidates[g] = ArrayUtilsRandomGenerator.randomSelect(candidates[g],
							-pop_diff_num_partner_12_months[g], getRNG());
				}

				for (AbstractIndividualInterface bR : candidates[g]) {
					Person_Bridging_Pop breakRelPerson = (Person_Bridging_Pop) bR;
					removeSingleRelationship(breakRelPerson);
				}
				// Removed and no longer able to form new partnership this turn.
				candidates[g] = new AbstractIndividualInterface[0];

			}

		}

		for (int map = 0; map < getRelMap().length; map++) {

			AbstractIndividualInterface[][] availiablePerson = new AbstractIndividualInterface[2][];

			if (map == RELMAP_HETRO) {
				availiablePerson[0] = candidates[GENDER_HETRO_MALE];
				availiablePerson[1] = candidates[GENDER_HETRO_FEMALE];
			} else {
				availiablePerson[0] = candidates[GENDER_MSMO];
				availiablePerson[1] = new AbstractIndividualInterface[0];
			}
			if (candidates[GENDER_MSMW].length > 0) {
				int orgLen = availiablePerson[0].length;
				availiablePerson[0] = Arrays.copyOf(availiablePerson[0],
						availiablePerson[0].length + candidates[GENDER_MSMW].length);
				System.arraycopy(candidates[GENDER_MSMW], 0, availiablePerson[0], orgLen,
						candidates[GENDER_MSMW].length);
			}

			// No need to be sorted
			getAvailability()[map].setAvailablePopulation(availiablePerson);
			// Generate new pairing
			int pairNum = getAvailability()[map].generatePairing();
			AbstractIndividualInterface[][] pairs = getAvailability()[map].getPairing();
			for (int pairId = 0; pairId < pairNum; pairId++) {
				// Assume no concurrency
				for (int p = 0; p < pairs[pairId].length; p++) {
					removeSingleRelationship((Person_Bridging_Pop) pairs[pairId][p]);
				}

				int duration = regPartDuration[map].sample();

				if (formRelationship(pairs[pairId], getRelMap()[map], duration, map) != null) {
					for (int p = 0; p < pairs[pairId].length; p++) {
						partformed[((Person_Bridging_Pop) pairs[pairId][p]).getGenderType()]++;
					}

				}

			}
		}

		return partformed;
	}

	protected void updateRelRelationshipMap() {
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
				} else {
					Integer[] link = relArr1.getLinks();

					ContactMap cMapAll = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL];
					ContactMap cMapSpec;

					if (map == 0) { // Hetro sex involving female
						cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_HETRO];
					} else { // Involve MSM
						cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_MSM];
					}

					checkContactMaps(link, new ContactMap[] { cMapAll, cMapSpec });

				}
			}
		}
	}

	private void checkContactMaps(Integer[] link, ContactMap[] cMaps) {
		for (ContactMap c : cMaps) {
			if (c != null) {
				for (Integer p : link) {
					if (!c.containsVertex(p)) {
						c.addVertex(p);
					}
				}
				if (!c.containsEdge(link[0], link[1])) {
					c.addEdge(link[0], link[1], link);
				}
			}

		}
	}

	private void removeSingleRelationship(Person_Bridging_Pop breakRelPerson) {
		RelationshipMap tarMap;

		if (breakRelPerson.getGenderType() == Person_Bridging_Pop.GENDER_TYPE_MSMW) {
			int degHetro = getRelMap()[RELMAP_HETRO].containsVertex(breakRelPerson.getId())
					? getRelMap()[RELMAP_HETRO].degreeOf(breakRelPerson.getId())
					: 0;
			int degMSM = getRelMap()[RELMAP_MSM].containsVertex(breakRelPerson.getId())
					? getRelMap()[RELMAP_MSM].degreeOf(breakRelPerson.getId())
					: 0;

			if (degHetro > 0 && degMSM > 0) {
				// 50:50 chance
				if (getRNG().nextFloat() < 0.5f) {
					tarMap = getRelMap()[RELMAP_MSM];
				} else {
					tarMap = getRelMap()[RELMAP_HETRO];
				}
			} else if (degMSM > 0) {
				tarMap = getRelMap()[RELMAP_MSM];
			} else {
				tarMap = getRelMap()[RELMAP_HETRO];
			}

		} else if (breakRelPerson.getGenderType() == Person_Bridging_Pop.GENDER_TYPE_MSMO) {
			tarMap = getRelMap()[RELMAP_MSM];
		} else {
			tarMap = getRelMap()[RELMAP_HETRO];
		}

		if (tarMap.containsVertex(breakRelPerson.getId())) {
			Set<SingleRelationship> edgeSet = tarMap.edgesOf(breakRelPerson.getId());
			if (!edgeSet.isEmpty()) {
				int rId = 0;
				if (edgeSet.size() > 1) {
					rId = getRNG().nextInt(edgeSet.size());
				}
				SingleRelationship relRemove = (SingleRelationship) edgeSet.toArray()[rId];
				tarMap.removeEdge(relRemove);
			}
		}

	}

}
