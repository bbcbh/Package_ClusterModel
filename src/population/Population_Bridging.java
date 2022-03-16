package population;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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
import relationship.SingleRelationshipTimeStamp;
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
			// float[]{MEAN_PARTNER_IN_12_MONTHS_BY_GENDER_TYPE}
			// Default: Hetro ASHR2 (16 - 30) , MSM - ASHR2
			// new float[] { (1.0f + 1.1f) / 2, (1.4f + 1.4f) / 2, 6.8f, 6.8f },
			// Alt:
			// float[]{CATEGORIES_LIMIT ..... ,
			// DISR_1, DISR_2...}
			new float[] { //
					0, 1, 2, 9, 49, 50, // Categories limit
					0.175f, 0.760f, 0.034f, 0.029f, 0.002f, 0.000f, // Female
					0.131f, 0.743f, 0.047f, 0.066f, 0.013f, 0.001f, // Male
					0.181f, 0.320f, 0.063f, 0.180f, 0.240f, 0.016f, // MSMO
					0.527f, 0.210f, 0.041f, 0.171f, 0.051f, 0.000f, }, // MSMW
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
			new float[] { 20.55f, 20.55f, 20.55f, 20.55f },

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
	public static final int LENGTH_RELMAP = 2;

	public static final int CONTACT_MAP_EDGE_P1 = 0;
	public static final int CONTACT_MAP_EDGE_P2 = CONTACT_MAP_EDGE_P1 + 1;
	public static final int CONTACT_MAP_EDGE_START_TIME = CONTACT_MAP_EDGE_P2 + 1;
	public static final int CONTACT_MAP_EDGE_DURATION = CONTACT_MAP_EDGE_START_TIME + 1;
	public static final int LENGTH_CONTACT_MAP_EDGE = CONTACT_MAP_EDGE_DURATION + 1;

	@SuppressWarnings("unchecked")
	private transient ArrayList<Person_Bridging_Pop>[] canSeekRelPartners = new ArrayList[LENGTH_GENDER];
	private transient int[][] canSeekRelPartners_catorgories_offset = new int[LENGTH_GENDER][];
	@SuppressWarnings("unchecked")
	private transient ArrayList<Person_Bridging_Pop>[] hasRelPartners = new ArrayList[LENGTH_GENDER];
	private transient int[][] hasRelPartners_catorgories_offset = new int[LENGTH_GENDER][];
	@SuppressWarnings("unchecked")
	private transient ArrayList<Person_Bridging_Pop>[] canSeekCasPartners = new ArrayList[LENGTH_GENDER];
	private transient int[][] canSeekCasPartners_catorgories_offset = new int[LENGTH_GENDER][];

	private transient Person_Bridging_Pop[][] casualPartnerFormed = null;

	private AbstractIntegerDistribution[] regPartDuration = new AbstractIntegerDistribution[LENGTH_RELMAP];

	private PrintStream printStatus = null;
	
	// TODO: Debug print - might remove in the final version
	private PrintStream[] partnerDistPrint = new PrintStream[LENGTH_GENDER];
	private java.io.File partnerDistPrintFolder = null;

	public void setPrintStatus(PrintStream printStatus) {
		this.printStatus = printStatus;
	}	
	public void setPartnerDistPrintFolder(java.io.File partnerDistPrintFolder) {
		this.partnerDistPrintFolder = partnerDistPrintFolder;
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
		Integer[] link = new Integer[LENGTH_CONTACT_MAP_EDGE];
		SingleRelationshipTimeStamp rel;
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

		link[CONTACT_MAP_EDGE_START_TIME] = getGlobalTime();
		link[CONTACT_MAP_EDGE_DURATION] = duration;

		checkContactMaps(link, new ContactMap[] { cMapAll, cMapSpec });

		rel = new SingleRelationshipTimeStamp(new Integer[] { link[0], link[1] });
		rel.setRelStartTime(getGlobalTime());

		if (relMap.addEdge(link[0], link[1], rel)) {
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

		int numReg = getNumRegularPartnersCurrently(person);

		boolean seekRegular = maxReg > 0; // numReg;

		if (((float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS])).length == LENGTH_GENDER) {
			// Mean (instead of distribution) version
			if (seekRegular) {
				int noPWin = AbstractIndividualInterface.ONE_YEAR_INT;
				if ((int) person.getField(Person_Bridging_Pop.FIELD_LAST_REGULAR_PARTNER_AT_AGE) > 0) {
					noPWin = Math.max(1, noPWin - ((int) person.getAge()
							- (int) person.getField(Person_Bridging_Pop.FIELD_LAST_REGULAR_PARTNER_AT_AGE)));
				}
				seekRegular = getRNG().nextInt(noPWin) < (maxReg - numReg);
			}
		}

		return seekRegular;
	}

	protected boolean seekingCasualToday(Person_Bridging_Pop person) {
		int maxCasual = (int) person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));
		int numCasual = person.getNumCasualInRecord();
		boolean seekCasual = maxCasual > numCasual;

		if (((float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS])).length == LENGTH_GENDER) {
			// Mean (instead of distribution) version
			if (seekCasual) {
				int noPWin = AbstractIndividualInterface.ONE_YEAR_INT;
				if ((int) person.getField(Person_Bridging_Pop.FIELD_LAST_CASUAL_PARNTER_AT_AGE) > 0) {
					noPWin = Math.max(1, noPWin - ((int) person.getAge()
							- (int) person.getField(Person_Bridging_Pop.FIELD_LAST_CASUAL_PARNTER_AT_AGE)));
				}
				seekCasual = getRNG().nextInt(noPWin) < (maxCasual - numCasual);
			}
		}

		return seekCasual;
	}

	protected void formCasualPartnership(int[] pop_diff_num_partner_12_months) {

		casualPartnerFormed = new Person_Bridging_Pop[0][];
		int partPt = 0;
		ContactMap[] cMaps;

		Person_Bridging_Pop[][] casualCandidate = candidatesForPartnerships(canSeekCasPartners,
				canSeekCasPartners_catorgories_offset, pop_diff_num_partner_12_months, false);

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
						checkContactMaps(new Integer[] { pair[0].getId(), pair[1].getId(), getGlobalTime(), 1 }, cMaps);
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
						checkContactMaps(new Integer[] { pair[0].getId(), pair[1].getId(), getGlobalTime(), 1 }, cMaps);
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
		RelationshipMap[] relMaps = new RelationshipMap[LENGTH_RELMAP];
		relMaps[RELMAP_HETRO] = new RelationshipMap();
		relMaps[RELMAP_MSM] = new RelationshipMap();
		this.setRelMap(relMaps);

		float[] meanRelDur = (float[]) getFields()[FIELD_MEAN_REG_PARTNERSHIP_DUR];
		regPartDuration[RELMAP_HETRO] = new PoissonDistribution(getRNG(), meanRelDur[RELMAP_HETRO],
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		regPartDuration[RELMAP_MSM] = new PoissonDistribution(getRNG(), meanRelDur[RELMAP_MSM],
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		
		// Initalise dist print folder (might remove later)
		
		if (partnerDistPrintFolder != null) {
			for (int i = 0; i < partnerDistPrint.length; i++) {
				try {
					partnerDistPrint[i] = new PrintStream(
							new java.io.File(partnerDistPrintFolder, "Gender_" + i + ".csv"));
				} catch (FileNotFoundException e) {
					e.printStackTrace(System.err);
				}
			}
		}

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
		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int[] population_num_partner_in_last_12_months = new int[field_mean_number_partner.length];

		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		if (field_mean_number_partner.length != LENGTH_GENDER) {
			int[] numInGrp = (int[]) (getFields()[FIELD_POP_COMPOSITION]);
			StringWriter wri = new StringWriter();
			PrintWriter pri = new PrintWriter(wri);
			pri.println("Ideal partnership distribution");
			for (int g = 0; g < LENGTH_GENDER; g++) {
				pri.print(String.format(" %d: [", g));
				for (int c = 0; c < numCat; c++) {
					if (c != 0) {
						pri.print(',');
					}
					pri.print(Math.round(field_mean_number_partner[numCat + g * numCat + c] * numInGrp[g]));
				}
				pri.println("]");
			}

			pri.close();
			printStatus.println(wri.toString());
			System.out.println(wri.toString());
		}

		if (numCat > 0) {
			for (int g = 0; g < LENGTH_GENDER; g++) {
				canSeekRelPartners_catorgories_offset[g] = new int[numCat];
				hasRelPartners_catorgories_offset[g] = new int[numCat];
				canSeekCasPartners_catorgories_offset[g] = new int[numCat];
			}

		}

		for (int genderGrpIndex = 0; genderGrpIndex < popSizes.length; genderGrpIndex++) {

			// Shared among same gender
			float[] p_type_prob = partner_type[genderGrpIndex];

			AbstractIntegerDistribution dist = null;

			// By distribution
			// All of them have no partner initially
			if (population_num_partner_in_last_12_months.length != LENGTH_GENDER) {
				population_num_partner_in_last_12_months[numCat + genderGrpIndex * numCat] = popSizes[genderGrpIndex];
			}

			for (int g = 0; g < popSizes[genderGrpIndex]; g++) {
				pop[popPt] = new Person_Bridging_Pop(popPt + 1, genderGrpIndex,
						18 * AbstractIndividualInterface.ONE_YEAR_INT, // Initial age - might not be used
						getGlobalTime(), 3); // 3 sites

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
						canSeekCasPartners_catorgories_offset[person.getGenderType()][0]++;
					}
				}

				if (seekingRegularToday(person)) {
					canSeekRelPartners[person.getGenderType()].add(person);
					canSeekRelPartners_catorgories_offset[person.getGenderType()][0]++;
				}
				popPt++;

			}
		}

		this.setPop(pop);
		this.getFields()[FIELDS_NEXT_ID] = popPt;

		// Set availability

		AbstractAvailability[] avail = new AbstractAvailability[LENGTH_RELMAP];

		avail[RELMAP_HETRO] = new Availability_Bridging_Population(getRNG());
		avail[RELMAP_HETRO].setParameter(Availability_Bridging_Population.BIPARTITE_MAPPING, true);
		avail[RELMAP_MSM] = new Availability_Bridging_Population(getRNG());
		avail[RELMAP_MSM].setParameter(Availability_Bridging_Population.BIPARTITE_MAPPING, false);
		this.setAvailability(avail);

		formPartnerships(population_num_partner_in_last_12_months);

	}

	public void formPartnerships(int[] population_num_partner_in_last_12_months) {
		int[] pop_diff_num_partner_12_months = cal_pop_diff_num_partner(population_num_partner_in_last_12_months);
		int[] reg_part_formed = formRegularPartnership(pop_diff_num_partner_12_months);

		if (pop_diff_num_partner_12_months.length == LENGTH_GENDER) {
			for (int g = 0; g < pop_diff_num_partner_12_months.length; g++) {
				pop_diff_num_partner_12_months[g] -= reg_part_formed[g];
			}
		} else {
			float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
			int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);
			for (int g = 0; g < LENGTH_GENDER; g++) {
				for (int c = 0; c < numCat; c++) {
					int index = numCat + g * numCat + c;
					if (reg_part_formed[index] > 0 && (c >= 1)
							&& field_mean_number_partner[c] - 1 <= field_mean_number_partner[c - 1]) {
						pop_diff_num_partner_12_months[index] -= reg_part_formed[index];
						pop_diff_num_partner_12_months[index - 1] += reg_part_formed[index];
					}
				}
			}

		}
		formCasualPartnership(pop_diff_num_partner_12_months);
	}

	@Override
	public void advanceTimeStep(int deltaT) {
		incrementTime(deltaT);

		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int[] population_num_partner_in_last_12_months = new int[field_mean_number_partner.length];
		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		// Clear partnership status

		for (int g = 0; g < LENGTH_GENDER; g++) {
			canSeekRelPartners[g].clear();
			hasRelPartners[g].clear();
			canSeekCasPartners[g].clear();
			if (numCat > 0) {
				Arrays.fill(canSeekRelPartners_catorgories_offset[g], 0);
				Arrays.fill(hasRelPartners_catorgories_offset[g], 0);
				Arrays.fill(canSeekCasPartners_catorgories_offset[g], 0);
			}
		}

		for (AbstractIndividualInterface p : this.getPop()) {
			Person_Bridging_Pop person = (Person_Bridging_Pop) p;
			person.incrementTime(deltaT, getInfList());
			int genderType = person.getGenderType();
			int numReg12Months = getNumRegularPartnersCurrently(person);
			int numCas12Months = person.getNumCasualInRecord();
			int numPart12Months = numReg12Months + numCas12Months;
			int cPt = -1;

			if (population_num_partner_in_last_12_months.length == LENGTH_GENDER) {
				population_num_partner_in_last_12_months[genderType] += numPart12Months;
			} else {
				cPt = findPartnershipCatogories(field_mean_number_partner, numCat, numPart12Months);
				population_num_partner_in_last_12_months[numCat + genderType * numCat + cPt]++;
			}

			if (seekingRegularToday(person)) {
				canSeekRelPartners[genderType].add(person);
				if (cPt != -1) {
					canSeekRelPartners_catorgories_offset[genderType][cPt]++;

				}
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
				if (cPt != -1) {
					hasRelPartners_catorgories_offset[genderType][cPt]++;
				}
			}

			if (seekingCasualToday(person)) {
				canSeekCasPartners[person.getGenderType()].add(person);
				if (cPt != -1) {
					canSeekCasPartners_catorgories_offset[genderType][cPt]++;
				}
			}

		}

		if (numCat > 0) {
			Comparator<Person_Bridging_Pop> cmp = new Comparator<Person_Bridging_Pop>() {
				@Override
				public int compare(Person_Bridging_Pop o1, Person_Bridging_Pop o2) {
					int cP1 = findPartnershipCatogories(field_mean_number_partner, numCat,
							getNumRegularPartnersCurrently(o1) + o1.getNumCasualInRecord());
					int cP2 = findPartnershipCatogories(field_mean_number_partner, numCat,
							getNumRegularPartnersCurrently(o2) + o2.getNumCasualInRecord());

					if (cP1 == cP2) {
						return Integer.compare(o1.getId(), o2.getId());
					} else {
						return Integer.compare(cP1, cP2);
					}
				}

			};
			for (int g = 0; g < LENGTH_GENDER; g++) {
				Collections.sort(canSeekRelPartners[g], cmp);
				Collections.sort(hasRelPartners[g], cmp);
				Collections.sort(canSeekCasPartners[g], cmp);
			}
		}

		if (printStatus != null) {
			if (population_num_partner_in_last_12_months.length == LENGTH_GENDER) {
				printStatus.println(String.format("# partners in last 12 months = %s - Day %d",
						Arrays.toString(population_num_partner_in_last_12_months), getGlobalTime()));

				System.out.println(String.format("# partners in last 12 months = %s - Day %d",
						Arrays.toString(population_num_partner_in_last_12_months), getGlobalTime()));
			} else {
				StringWriter wri = new StringWriter();
				PrintWriter pri = new PrintWriter(wri);
				pri.println(String.format("# partners in last 12 months -  Day %d", getGlobalTime()));
				for (int g = 0; g < LENGTH_GENDER; g++) {
					pri.println(String.format(" %d: %s", g,
							Arrays.toString(Arrays.copyOfRange(population_num_partner_in_last_12_months,
									numCat + g * numCat, 2 * numCat + g * numCat))));
					
					
					partnerDistPrint[g].print(getGlobalTime());
					
					if(partnerDistPrintFolder != null) {
					for(int n = 0; n < numCat; n++) {
						partnerDistPrint[g].print(',');
						partnerDistPrint[g].print(population_num_partner_in_last_12_months[numCat + g * numCat + n]);
					}
					partnerDistPrint[g].println();
					}
				}

				pri.close();
				printStatus.println(wri.toString());
				System.out.println(wri.toString());
			}

			printStatus.println("Has Regular Partners:");
			for (int i = 0; i < hasRelPartners.length; i++) {
				printStatus.println(String.format(" %d: %d", i, hasRelPartners[i].size()));
			}
			printStatus.println("Seeking Casual Partners:");
			for (int i = 0; i < canSeekCasPartners.length; i++) {
				printStatus.println(String.format(" %d: %d", i, canSeekCasPartners[i].size()));
			}

		}
		updateRelRelationshipMap();
		formPartnerships(population_num_partner_in_last_12_months);

	}

	private int getNumRegularPartnersCurrently(Person_Bridging_Pop person) {
		int numReg12Currently = person.getNumRegularInRecord();
		// Has current partner with partnership longer for 12 months
		for (int r = 0; r < LENGTH_RELMAP; r++) {
			if (getRelMap()[r].containsVertex(person.getId()) && getRelMap()[r].degreeOf(person.getId()) > 0) {
				for (SingleRelationship rel : getRelMap()[r].edgesOf(person.getId())) {
					if (((SingleRelationshipTimeStamp) rel).getRelStartTime() < getGlobalTime()
							- AbstractIndividualInterface.ONE_YEAR_INT) {
						numReg12Currently++;
					}
				}
			}
		}

		return numReg12Currently;
	}

	private int findPartnershipCatogories(float[] field_mean_number_partner, int numCat, int numPart12Months) {
		int cPt = Arrays.binarySearch(field_mean_number_partner, 0, numCat, numPart12Months);
		if (cPt < 0) {
			cPt = ~cPt;
			if (cPt >= numCat) { // Last one is +
				cPt = numCat - 1;
			}
		}
		return cPt;
	}

	private int[] cal_pop_diff_num_partner(int[] population_num_partner_in_last_12_months) {
		float[] field_mean_target = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int[] numInGrp = (int[]) (getFields()[FIELD_POP_COMPOSITION]);
		int[] pop_diff;

		if (field_mean_target.length == LENGTH_GENDER) {
			pop_diff = new int[LENGTH_GENDER];

			for (int g = 0; g < field_mean_target.length; g++) {
				pop_diff[g] = Math
						.round(field_mean_target[g] * numInGrp[g] - population_num_partner_in_last_12_months[g]);
			}
		} else {
			// float[]{CATORGORIES_LIMIT ..... ,
			// DISR_1, DISR_2...}

			pop_diff = new int[field_mean_target.length];
			int numCat = field_mean_target.length / (1 + LENGTH_GENDER);

			for (int g = 0; g < LENGTH_GENDER; g++) {
				for (int c = 0; c < numCat; c++) {
					int pdIndex = numCat + g * numCat + c;

					pop_diff[pdIndex] = Math.round(field_mean_target[pdIndex] * numInGrp[g]
							- population_num_partner_in_last_12_months[pdIndex]);

				}

			}
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

		int[] partformed = new int[pop_diff_num_partner_12_months.length];

		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		AbstractIndividualInterface[][] candidates = candidatesForPartnerships(canSeekRelPartners,
				canSeekRelPartners_catorgories_offset, pop_diff_num_partner_12_months, true);

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
						Person_Bridging_Pop person = (Person_Bridging_Pop) pairs[pairId][p];

						if (partformed.length == LENGTH_GENDER) {
							partformed[person.getGenderType()]++;
						} else {
							int numPartIn12Months = person.getNumCasualInRecord()
									+ getNumRegularPartnersCurrently(person);
							int cPt = 0;
							for (int cI = 0; cI < numCat; cI++) {
								if (numPartIn12Months > field_mean_number_partner[1 + cI]) {
									cPt++;
								}
							}
							partformed[numCat + person.getGenderType() * numCat + cPt]++;
						}
					}

				}

			}
		}

		return partformed;
	}

	protected Person_Bridging_Pop[][] candidatesForPartnerships(ArrayList<Person_Bridging_Pop>[] candidateCollection,
			int[][] categories_offset, int[] pop_diff_num_partner_12_months, boolean allowedBreakExistingPartnership) {

		Person_Bridging_Pop[][] candidates = new Person_Bridging_Pop[LENGTH_GENDER][];
		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		if (pop_diff_num_partner_12_months.length == LENGTH_GENDER) {

			for (int g = 0; g < pop_diff_num_partner_12_months.length; g++) {

				if (pop_diff_num_partner_12_months[g] > 0) {
					candidates[g] = candidateCollection[g]
							.toArray(new Person_Bridging_Pop[candidateCollection[g].size()]);
					if (pop_diff_num_partner_12_months[g] < candidateCollection[g].size()) {
						candidates[g] = ArrayUtilsRandomGenerator.randomSelect(candidates[g],
								pop_diff_num_partner_12_months[g], getRNG());
					}

				} else if (pop_diff_num_partner_12_months[g] <= 0 && allowedBreakExistingPartnership) {

					candidates[g] = hasRelPartners[g].toArray(new Person_Bridging_Pop[hasRelPartners[g].size()]);

					if (-pop_diff_num_partner_12_months[g] < candidates[g].length) {
						candidates[g] = ArrayUtilsRandomGenerator.randomSelect(candidates[g],
								-pop_diff_num_partner_12_months[g], getRNG());
					}

					for (AbstractIndividualInterface bR : candidates[g]) {
						Person_Bridging_Pop breakRelPerson = (Person_Bridging_Pop) bR;
						removeSingleRelationship(breakRelPerson);
					}
					// Removed and no longer able to form new partnership this turn.
					candidates[g] = new Person_Bridging_Pop[0];
				}

			}
		} else {

			for (int g = 0; g < LENGTH_GENDER; g++) {

				int[] diffByCatgories = Arrays.copyOfRange(pop_diff_num_partner_12_months, numCat + g * numCat,
						2 * numCat + g * numCat);

				ArrayList<Person_Bridging_Pop> seekingList = new ArrayList<>();
				int[] catOffset = Arrays.copyOf(categories_offset[g], categories_offset[g].length);


				for (int c = numCat - 1; c > 0; c--) {

					int pickedForCatgoriesToday = diffByCatgories[c];
					
					// TODO: At this stage trial and errors
					int perDaySeekfactor = AbstractIndividualInterface.ONE_YEAR_INT;					
					pickedForCatgoriesToday = ((int) (pickedForCatgoriesToday/perDaySeekfactor));
					int remainder = diffByCatgories[c] % perDaySeekfactor;

					if(remainder != 0) {
						if(getRNG().nextInt(perDaySeekfactor) < remainder) {
							pickedForCatgoriesToday++;
						}
					}

					while (pickedForCatgoriesToday > 0) {
						// Add candidate from next smallest categories
						int selCat = c - 1;
						while (selCat >= 0 && catOffset[selCat] == 0) {
							selCat--;
						}

						if (selCat >= 0) {

							int startOff = 0;

							for (int i = 0; i < selCat; i++) {
								startOff += catOffset[i];
							}

							int selCandidate = startOff + getRNG().nextInt(catOffset[selCat]);

							Person_Bridging_Pop candidate = candidateCollection[g].remove(selCandidate);
							seekingList.add(candidate);

							if (allowedBreakExistingPartnership) {
								removeSingleRelationship(candidate);
							}
							pickedForCatgoriesToday--;
							catOffset[selCat]--;
						} else {
							break;
						}
					}
				}

				candidates[g] = seekingList.toArray(new Person_Bridging_Pop[seekingList.size()]);
			}

		}
		return candidates;
	}

	protected void updateRelRelationshipMap() {
		for (int map = 0; map < getRelMap().length; map++) {
			RelationshipMap relMap = getRelMap()[map];

			// Update existing
			SingleRelationship[] relArr = relMap.getRelationshipArray();

			if (relMap.edgeSet().size() != relArr.length) {
				relArr = relMap.edgeSet().toArray(new SingleRelationship[relMap.edgeSet().size()]);
			}

			ContactMap cMapAll = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL];
			ContactMap cMapSpec;

			if (map == 0) { // Hetro sex involving female
				cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_HETRO];
			} else { // Involve MSM
				cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_MSM];
			}

			for (SingleRelationship rel : relArr) {
				Integer[] link = rel.getLinks();
				link = Arrays.copyOf(link, LENGTH_CONTACT_MAP_EDGE);
				link[CONTACT_MAP_EDGE_START_TIME] = getGlobalTime();
				if (rel.incrementTime(1) <= 0) {
					relMap.removeEdge(rel);
					link[CONTACT_MAP_EDGE_DURATION] = -1;
				} else {
					link[CONTACT_MAP_EDGE_DURATION] = (int) rel.getDurations();
				}
				checkContactMaps(link, new ContactMap[] { cMapAll, cMapSpec });
			}
		}
	}

	private void checkContactMaps(Integer[] link, ContactMap[] cMaps) {
		for (ContactMap c : cMaps) {
			if (c != null) {
				if (!c.containsVertex(link[CONTACT_MAP_EDGE_P1])) {
					c.addVertex(link[CONTACT_MAP_EDGE_P1]);
				}
				if (!c.containsVertex(link[CONTACT_MAP_EDGE_P2])) {
					c.addVertex(link[CONTACT_MAP_EDGE_P2]);
				}

				if (!c.containsEdge(link[CONTACT_MAP_EDGE_P1], link[CONTACT_MAP_EDGE_P2])
						&& link[CONTACT_MAP_EDGE_DURATION] > 0) {
					c.addEdge(link[CONTACT_MAP_EDGE_P1], link[CONTACT_MAP_EDGE_P2], link);
				}
				if (c.containsEdge(link[CONTACT_MAP_EDGE_P1], link[CONTACT_MAP_EDGE_P2])
						&& link[CONTACT_MAP_EDGE_DURATION] > 0) {
					Integer[] e = c.getEdge(link[CONTACT_MAP_EDGE_P1], link[CONTACT_MAP_EDGE_P2]);
					e[CONTACT_MAP_EDGE_DURATION] += link[CONTACT_MAP_EDGE_DURATION];
				}

				if (c.containsEdge(link[CONTACT_MAP_EDGE_P1], link[CONTACT_MAP_EDGE_P2])
						&& link[CONTACT_MAP_EDGE_DURATION] < 0) {
					Integer[] e = c.getEdge(link[CONTACT_MAP_EDGE_P1], link[CONTACT_MAP_EDGE_P2]);
					e[CONTACT_MAP_EDGE_DURATION] = link[CONTACT_MAP_EDGE_START_TIME] - e[CONTACT_MAP_EDGE_START_TIME];

				}
			}
		}

	}

	private void removeSingleRelationship(Person_Bridging_Pop breakRelPerson) {

		int tarMapIndex = -1;

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
					tarMapIndex = RELMAP_MSM;
				} else {
					tarMapIndex = RELMAP_HETRO;
				}
			} else if (degMSM > 0) {
				tarMapIndex = RELMAP_MSM;
			} else {
				tarMapIndex = RELMAP_HETRO;
			}

		} else if (breakRelPerson.getGenderType() == Person_Bridging_Pop.GENDER_TYPE_MSMO) {
			tarMapIndex = RELMAP_MSM;
		} else {
			tarMapIndex = RELMAP_HETRO;
		}

		RelationshipMap tarMap = getRelMap()[tarMapIndex];
		ContactMap cMapAll = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL];
		ContactMap cMapSpec;

		if (tarMapIndex == RELMAP_HETRO) { // Hetro sex involving female
			cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_HETRO];
		} else { // Involve MSM
			cMapSpec = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_MSM];
		}

		if (tarMap.containsVertex(breakRelPerson.getId())) {
			Set<SingleRelationship> edgeSet = tarMap.edgesOf(breakRelPerson.getId());
			if (!edgeSet.isEmpty()) {

				SingleRelationship[] repArr = edgeSet.toArray(new SingleRelationship[edgeSet.size()]);

				if (repArr.length > 1) {
					Arrays.sort(repArr, new Comparator<SingleRelationship>() {
						@Override
						public int compare(SingleRelationship o1, SingleRelationship o2) {
							return Double.compare(o1.incrementTime(0), o2.incrementTime(0));
						}
					});
				}
				SingleRelationship relRemove = repArr[0];

				// Update contact map
				Integer[] link = relRemove.getLinks();
				link = Arrays.copyOf(link, LENGTH_CONTACT_MAP_EDGE);
				link[CONTACT_MAP_EDGE_START_TIME] = getGlobalTime();
				link[CONTACT_MAP_EDGE_DURATION] = -1;
				checkContactMaps(link, new ContactMap[] { cMapAll, cMapSpec });

				tarMap.removeEdge(relRemove);

			}
		}

	}

}
