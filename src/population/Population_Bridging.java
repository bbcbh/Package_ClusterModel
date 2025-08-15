package population;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import availability.AbstractAvailability;
import infection.AbstractInfection;
import person.AbstractIndividualInterface;
import population.availability.Availability_Bridging_Population;
import population.person.Person_Bridging_Pop;
import relationship.ContactMap;
import relationship.RelationshipMap;
import relationship.SingleRelationship;
import relationship.SingleRelationshipTimeStamp;
import sim.Abstract_Runnable_ClusterModel;
import util.ArrayUtilsRandomGenerator;

public class Population_Bridging extends AbstractFieldsArrayPopulation {

	/**
	 * 
	 */
	private static final long serialVersionUID = 202504100918L;

	public final Object[] DEFAULT_BRIDGING_POP_FIELDS = {
			// FIELD_POP_COMPOSITION
			// int[] {NUM_FEMALE, NUM_MALE, NUM_MSMO, NUM_MSMW}
			new int[] { 5000, 5000, 200, 50 },
			// FIELD_CONTACT_MAP
			new ContactMap[1],
			// FIELD_PARTNER_TYPE_PROB
			// float[GENDER]{CUMUL_REG_ONLY, CUMUL_CAS_ONLY}
			// Alt:
			// float[GENDER]{CUMUL_REG_ONLY, CUMUL_CAS_ONLY, ASSOCATIVINESS}
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

	public static final int PARTNER_TYPE_INDEX_REG_ONLY = 0;
	public static final int PARTNER_TYPE_INDEX_CAS_ONLY = PARTNER_TYPE_INDEX_REG_ONLY + 1;
	public static final int PARTNER_TYPE_ASSORTATIVITY = PARTNER_TYPE_INDEX_CAS_ONLY + 1;
	public static final int PARTNER_TYPE_NON_MAPPED_ENCOUNTER_PROB = PARTNER_TYPE_ASSORTATIVITY + 1;
	public static final int PARTNER_TYPE_NON_MAPPED_ENCOUNTER_TARGET_GENDER = PARTNER_TYPE_NON_MAPPED_ENCOUNTER_PROB
			+ 1;

	public static final int GENDER_FEMALE = 0;
	public static final int GENDER_HETRO_MALE = GENDER_FEMALE + 1;
	public static final int GENDER_MSMO = GENDER_HETRO_MALE + 1;
	public static final int GENDER_MSMW = GENDER_MSMO + 1;
	public static final int LENGTH_GENDER = GENDER_MSMW + 1;

	public static final int CONTACT_MAP_ALL = 0;

	public static final int RELMAP_HETRO = 0;
	public static final int RELMAP_MSM = 1;
	public static final int LENGTH_RELMAP = 2;

	public static final String STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS = "STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS";

	protected float prob_no_bridge = -1f;

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
	protected transient HashMap<String, Object> stepwise_output = null;
	protected transient File baseDir = null;

	private static final Comparator<Person_Bridging_Pop> COMPARATOR_CASUAL_MIXING = new Comparator<Person_Bridging_Pop>() {
		@Override
		public int compare(Person_Bridging_Pop o1, Person_Bridging_Pop o2) {
			int m1 = (int) o1.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));
			int m2 = (int) o2.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));
			return Integer.compare(m1, m2);
		}
	};

	private static final int CASUAL_MIXING_SD_RANGE = 2;

	protected AbstractIntegerDistribution[] regPartDuration = new AbstractIntegerDistribution[LENGTH_RELMAP];
	protected PrintStream[] printStatus = null;

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

	public static Population_Bridging decodeFromStream(java.io.ObjectInputStream inStr)
			throws IOException, ClassNotFoundException {
		int globalTime = inStr.readInt();
		AbstractInfection[] infList = (AbstractInfection[]) inStr.readObject();
		Object[] decoded_fields = (Object[]) inStr.readObject();

		Population_Bridging pop = new Population_Bridging((long) decoded_fields[0]);

		if (decoded_fields.length != pop.getFields().length) {
			int oldLen = decoded_fields.length;
			decoded_fields = Arrays.copyOf(decoded_fields, pop.getFields().length);
			for (int i = oldLen; i < decoded_fields.length; i++) {
				decoded_fields[i] = pop.getFields()[i];
			}
		}
		pop.setGlobalTime(globalTime);
		pop.setInfList(infList);
		pop.setFields(decoded_fields);

		// Initiation of transient fields
		pop.initialiseTransientFields();

		return pop;

	}

	public void setProb_no_bridge(float prob_no_bridge) {
		this.prob_no_bridge = prob_no_bridge;
	}

	protected void initialiseTransientFields() {
		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		if (numCat > 0) {
			for (int g = 0; g < LENGTH_GENDER; g++) {
				canSeekRelPartners_catorgories_offset[g] = new int[numCat];
				hasRelPartners_catorgories_offset[g] = new int[numCat];
				canSeekCasPartners_catorgories_offset[g] = new int[numCat];
			}

		}

		float[] meanRelDur = (float[]) getFields()[FIELD_MEAN_REG_PARTNERSHIP_DUR];
		regPartDuration[RELMAP_HETRO] = new PoissonDistribution(getRNG(), meanRelDur[RELMAP_HETRO],
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		regPartDuration[RELMAP_MSM] = new PoissonDistribution(getRNG(), meanRelDur[RELMAP_MSM],
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
	}

	public Class<? extends Object> getFieldsClass(int fieldIndex) {
		return this.getFields()[fieldIndex].getClass();
	}

	@Override
	protected SingleRelationship formRelationship(AbstractIndividualInterface[] pair, RelationshipMap relMap,
			int duration, int mapType) {

		Integer[] link = new Integer[Abstract_Runnable_ClusterModel.LENGTH_CONTACT_MAP_EDGE];
		SingleRelationshipTimeStamp rel;
		ContactMap cMapAll = ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL];

		for (int i = 0; i < pair.length; i++) {
			if (!relMap.containsVertex(pair[i].getId())) {
				relMap.addVertex(pair[i].getId());
			}

			link[i] = pair[i].getId();
		}

		link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] = getGlobalTime();
		link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] = duration;

		checkContactMaps(link, new ContactMap[] { cMapAll });

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
				cMaps = new ContactMap[] { ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL] };

				casualPartnerFormed = Arrays.copyOf(casualPartnerFormed,
						casualPartnerFormed.length + numMSMCasualToFormed);

				// Assortative mixing (MSM)
				Arrays.sort(casualMSMArr, COMPARATOR_CASUAL_MIXING);

				ArrayList<Person_Bridging_Pop> casualMSM_Person_List = new ArrayList<>(List.of(casualMSMArr));
				ArrayList<Integer> casualMSM_MaxP_List = new ArrayList<>();

				for (Person_Bridging_Pop msm : casualMSM_Person_List) {
					casualMSM_MaxP_List.add((int) msm
							.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS)));
				}

				partPt = casualPartnerPairing(partPt, casualMSM_Person_List, casualMSM_MaxP_List, casualMSM_Person_List,
						casualMSM_MaxP_List, this.getRelMap()[RELMAP_MSM], false, numMSMCasualToFormed, cMaps);
			}
		}

		Person_Bridging_Pop[] bridging_candidate = casualCandidate[GENDER_MSMW];

		if (prob_no_bridge > 0) {
			int numMSM_bridge = Math.round(bridging_candidate.length * Math.max(0, 1 - prob_no_bridge));
			if (numMSM_bridge > 0) {
				bridging_candidate = ArrayUtilsRandomGenerator.randomSelect(bridging_candidate, numMSM_bridge,
						getRNG());
			} else {
				bridging_candidate = new Person_Bridging_Pop[0];
			}

		}

		// Heterosexual
		if (casualCandidate[GENDER_FEMALE].length > 0
				&& (casualCandidate[GENDER_HETRO_MALE].length + bridging_candidate.length) > 0)

		{

			int numHetroCasualToFormed = Math.min(casualCandidate[GENDER_FEMALE].length,
					casualCandidate[GENDER_HETRO_MALE].length + bridging_candidate.length);

			if (numHetroCasualToFormed > 0) {

				cMaps = new ContactMap[] { ((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL] };

				Person_Bridging_Pop[] casualFemale = casualCandidate[GENDER_FEMALE];
				Person_Bridging_Pop[] casualMale = Arrays.copyOf(casualCandidate[GENDER_HETRO_MALE],
						casualCandidate[GENDER_HETRO_MALE].length + bridging_candidate.length);

				if (bridging_candidate.length > 0) {
					System.arraycopy(bridging_candidate, 0, casualMale, casualCandidate[GENDER_HETRO_MALE].length,
							bridging_candidate.length);
				}

				// Assortative mixing (Hetrosexual)
				Arrays.sort(casualFemale, COMPARATOR_CASUAL_MIXING);
				Arrays.sort(casualMale, COMPARATOR_CASUAL_MIXING);

				ArrayList<Person_Bridging_Pop> f_list = new ArrayList<>(List.of(casualFemale));
				ArrayList<Person_Bridging_Pop> m_list = new ArrayList<>(List.of(casualMale));

				ArrayList<Integer> f_max = new ArrayList<>();
				for (Person_Bridging_Pop f : casualFemale) {
					f_max.add((int) f
							.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS)));
				}
				ArrayList<Integer> m_max = new ArrayList<>();
				for (Person_Bridging_Pop m : casualMale) {
					m_max.add((int) m
							.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS)));
				}

				casualPartnerFormed = Arrays.copyOf(casualPartnerFormed,
						casualPartnerFormed.length + numHetroCasualToFormed);

				boolean src_is_female = f_list.size() < m_list.size();

				ArrayList<Person_Bridging_Pop> src_list = src_is_female ? f_list : m_list;
				ArrayList<Integer> src_max_list = src_list == f_list ? f_max : m_max;
				ArrayList<Person_Bridging_Pop> tar_base_list = src_is_female ? m_list : f_list;
				ArrayList<Integer> tar_base_max_list = src_is_female ? m_max : f_max;

				partPt = casualPartnerPairing(partPt, src_list, src_max_list, tar_base_list, tar_base_max_list,
						this.getRelMap()[RELMAP_HETRO], true, numHetroCasualToFormed, cMaps);

			}
		}

		if (partPt < casualPartnerFormed.length) {
			casualPartnerFormed = Arrays.copyOf(casualPartnerFormed, partPt);
		}

		// System.out.printf("# casual = %d at t = %d\n", casualPartnerFormed.length,
		// getGlobalTime());
	}

	private int casualPartnerPairing(int partershipPt, ArrayList<Person_Bridging_Pop> src_list,
			ArrayList<Integer> src_max_list, ArrayList<Person_Bridging_Pop> tar_base_list,
			ArrayList<Integer> tar_base_max_list, RelationshipMap relMapCheck, boolean rel_order_by_gender,
			int numCasualToFormed, ContactMap[] existingContactMap) {

		while (numCasualToFormed > 0 && src_list.size() > 0) {

			int src_index = getRNG().nextInt(src_list.size());
			Person_Bridging_Pop src_person = src_list.remove(src_index);
			Integer src_max = src_max_list.remove(src_index);

			Person_Bridging_Pop tar_person = null;

			boolean randMix = true;

			float[] part_type = ((float[][]) getFields()[FIELD_PARTNER_TYPE_PROB])[src_person.getGenderType()];
			if (PARTNER_TYPE_ASSORTATIVITY < part_type.length) {
				float mixProb = part_type[PARTNER_TYPE_ASSORTATIVITY];
				if (mixProb >= 1) {
					System.err.printf("Warning: PARTNER_TYPE_ASSORTATIVITY = %f > 1\n", mixProb);
				}
				randMix = mixProb < getRNG().nextFloat();
			}

			ArrayList<Person_Bridging_Pop> tar_list = new ArrayList<>(tar_base_list);
			ArrayList<Integer> tar_max_partner = new ArrayList<>(tar_base_max_list);

			while (tar_person == null && !tar_list.isEmpty()) {

				tar_person = findPartnerFromList(randMix, src_person, src_max, tar_list, tar_max_partner, relMapCheck);

				if (tar_person != null) {
					Person_Bridging_Pop[] pair;

					// For consistency
					if (rel_order_by_gender) {
						if (src_person.getGenderType() == Population_Bridging.GENDER_FEMALE) {
							pair = new Person_Bridging_Pop[] { tar_person, src_person };
						} else {
							pair = new Person_Bridging_Pop[] { src_person, tar_person };
						}
					} else {
						if (src_person.getId() > tar_person.getId()) {
							pair = new Person_Bridging_Pop[] { tar_person, src_person };
						} else {
							pair = new Person_Bridging_Pop[] { src_person, tar_person };
						}
					}

					// Add casual partner
					for (int p = 0; p < pair.length; p++) {
						pair[p].addCasualPartner(pair[(p + 1) % 2]);
					}

					casualPartnerFormed[partershipPt] = new Person_Bridging_Pop[] { pair[0], pair[1] };
					checkContactMaps(new Integer[] { pair[0].getId(), pair[1].getId(), getGlobalTime(), 1 },
							existingContactMap);

					// Clean up relevant list
					int removeIndex = -1;
					for (int r = 0; r < tar_base_list.size(); r++) {
						Person_Bridging_Pop removePerson = tar_base_list.get(r);
						if (removePerson.getId() == tar_person.getId()) {
							removeIndex = r;
							break;
						}
					}
					if (removeIndex >= 0) {
						tar_base_list.remove(removeIndex);
						tar_base_max_list.remove(removeIndex);
					}

					partershipPt++;
					numCasualToFormed--;
				}
			}

		}
		return partershipPt;
	}

	private Person_Bridging_Pop findPartnerFromList(boolean randomMixing, Person_Bridging_Pop src_person,
			Integer src_behavior, ArrayList<Person_Bridging_Pop> list_tar, ArrayList<Integer> list_tar_behavior,
			RelationshipMap relMap_check) {
		int tar_index;
		Person_Bridging_Pop tar_person;
		if (randomMixing || list_tar.size() == 1) {
			tar_index = getRNG().nextInt(list_tar.size());
		} else {
			int[] indices = new int[] {
					Collections.binarySearch(list_tar_behavior, src_behavior - CASUAL_MIXING_SD_RANGE),
					Collections.binarySearch(list_tar_behavior, src_behavior),
					Collections.binarySearch(list_tar_behavior, src_behavior + CASUAL_MIXING_SD_RANGE), };

			for (int i = 0; i < indices.length; i++) {
				if (indices[i] < 0) {
					indices[i] = ~indices[i];
				}
			}

			if (indices[2] > indices[0]) {
				float sd = (indices[2] - indices[0]) / 2f;
				tar_index = -1;
				// Resample until have a valid match
				while (tar_index < 0 || tar_index >= list_tar.size()) {
					tar_index = (int) Math.round(getRNG().nextGaussian() * sd + indices[1]);
				}
			} else {
				// Go back to random if SD = 0 (e.g. all the same);
				tar_index = getRNG().nextInt(list_tar.size());
			}
		}

		tar_person = list_tar.remove(tar_index);
		list_tar_behavior.remove(tar_index);

		if (relMap_check != null) {
			// For consistency
			boolean swapPerson;

			if ((src_person.getGenderType() & tar_person.getGenderType()) == Population_Bridging.GENDER_FEMALE) {
				swapPerson = src_person.getGenderType() == Population_Bridging.GENDER_FEMALE;
			} else {
				swapPerson = src_person.getId() > tar_person.getId();
			}

			if (relMap_check.containsEdge(swapPerson ? tar_person.getId() : src_person.getId(),
					!swapPerson ? tar_person.getId() : src_person.getId())) {
				tar_person = null;
			}
		}
		return tar_person;
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

		if (field_mean_number_partner.length != LENGTH_GENDER && printStatus != null) {
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
			if (printStatus != null) {
				for (PrintStream out : printStatus) {
					out.println(wri.toString());
				}
			}

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

			AbstractIntegerDistribution casual_partner_dist_by_gender_grp = null;

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

				float prob = getRNG().nextFloat();

				int prob_pt = 0;
				while (prob_pt < PARTNER_TYPE_ASSORTATIVITY && prob > p_type_prob[prob_pt]) {
					prob_pt++;
				}

				if (prob_pt == PARTNER_TYPE_INDEX_REG_ONLY) { // CUMUL_REG_ONLY
					person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS), 0);

				} else if (prob_pt == PARTNER_TYPE_INDEX_CAS_ONLY) { // CUMUL_CAS_ONLY
					person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS), 0);
				}

				// Set max. number of casual in 12 months
				if (casual_partner_dist_by_gender_grp == null) {
					casual_partner_dist_by_gender_grp = new PoissonDistribution(getRNG(),
							mean_casual_partners[genderGrpIndex], PoissonDistribution.DEFAULT_EPSILON,
							PoissonDistribution.DEFAULT_MAX_ITERATIONS);
				}

				initialiseNumPartners(person, genderGrpIndex, g, casual_partner_dist_by_gender_grp);
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

		if (stepwise_output != null) {
			stepwise_output.put(STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS, population_num_partner_in_last_12_months);
		}

		formPartnerships(population_num_partner_in_last_12_months);

	}

	protected void initialiseNumPartners(Person_Bridging_Pop person, int gender_grp_index, int gender_grp_count,
			AbstractIntegerDistribution casual_partner_dist_by_gender_grp) {

		if (!person.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS)).equals(0)) {
			person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS),
					Math.max(1, casual_partner_dist_by_gender_grp.sample()));
			if (seekingCasualToday(person)) {
				canSeekCasPartners[person.getGenderType()].add(person);
				canSeekCasPartners_catorgories_offset[person.getGenderType()][0]++;
			}
		}

		if (seekingRegularToday(person)) {
			canSeekRelPartners[person.getGenderType()].add(person);
			canSeekRelPartners_catorgories_offset[person.getGenderType()][0]++;
		}

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
				for (PrintStream out : printStatus) {
					out.printf("# partners in last 12 months = %s - Day %d\n",
							Arrays.toString(population_num_partner_in_last_12_months), getGlobalTime());
				}

			} else {
				StringWriter wri = new StringWriter();
				PrintWriter pri = new PrintWriter(wri);
				pri.println(String.format("# partners in last 12 months -  Day %d", getGlobalTime()));
				for (int g = 0; g < LENGTH_GENDER; g++) {
					pri.println(String.format(" %d: %s", g,
							Arrays.toString(Arrays.copyOfRange(population_num_partner_in_last_12_months,
									numCat + g * numCat, 2 * numCat + g * numCat))));

				}
				pri.close();

				for (PrintStream out : printStatus) {
					out.println(wri.toString());
				}

			}

		}
		updateRelRelationshipMap();

		if (stepwise_output != null) {
			stepwise_output.put(STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS, population_num_partner_in_last_12_months);
		}
		formPartnerships(population_num_partner_in_last_12_months);

	}

	protected int getNumRegularPartnersCurrently(Person_Bridging_Pop person) {
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

	protected int findPartnershipCatogories(float[] field_mean_number_partner, int numCat, int numPart12Months) {
		int cPt = Arrays.binarySearch(field_mean_number_partner, 0, numCat, numPart12Months);
		if (cPt < 0) {
			cPt = ~cPt;
			if (cPt >= numCat) { // Last one is +
				cPt = numCat - 1;
			}
		}
		return cPt;
	}

	protected int[] cal_pop_diff_num_partner(int[] population_num_partner_in_last_12_months) {
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
				availiablePerson[1] = candidates[GENDER_FEMALE];
			} else {
				availiablePerson[0] = candidates[GENDER_MSMO];
				availiablePerson[1] = new AbstractIndividualInterface[0];
			}

			if (candidates[GENDER_MSMW].length > 0) {
				AbstractIndividualInterface[] candidate_bridge = candidates[GENDER_MSMW];
				if (prob_no_bridge > 0) {
					int numMSM_bridge = Math.round(candidates[GENDER_MSMW].length * Math.max(0, 1 - prob_no_bridge));
					if (numMSM_bridge > 0) {
						candidate_bridge = ArrayUtilsRandomGenerator.randomSelect(candidate_bridge, numMSM_bridge,
								getRNG());
					} else {
						candidate_bridge = new AbstractIndividualInterface[0];
					}

				}

				if (candidate_bridge.length > 0) {
					int orgLen = availiablePerson[0].length;
					availiablePerson[0] = Arrays.copyOf(availiablePerson[0],
							availiablePerson[0].length + candidate_bridge.length);
					System.arraycopy(candidate_bridge, 0, availiablePerson[0], orgLen, candidate_bridge.length);
				}
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

					int perDaySeekfactor = AbstractIndividualInterface.ONE_YEAR_INT;
					pickedForCatgoriesToday = ((int) (pickedForCatgoriesToday / perDaySeekfactor));
					int remainder = diffByCatgories[c] % perDaySeekfactor;

					if (remainder != 0) {
						if (getRNG().nextInt(perDaySeekfactor) < remainder) {
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

			for (SingleRelationship rel : relArr) {
				Integer[] link = rel.getLinks();
				link = Arrays.copyOf(link, Abstract_Runnable_ClusterModel.LENGTH_CONTACT_MAP_EDGE);
				link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] = ((SingleRelationshipTimeStamp) rel)
						.getRelStartTime();
				if (rel.incrementTime(1) <= 0) {
					relMap.removeEdge(rel);
					link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] = -1;
				} else {
					link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] = 0;
				}
				checkContactMaps(link, new ContactMap[] { cMapAll });
			}
		}
	}

	protected void checkContactMaps(Integer[] link, ContactMap[] cMaps) {
		for (ContactMap c : cMaps) {
			if (c != null) {

				if (!c.containsVertex(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1])) {
					c.addVertex(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]);
				}
				if (!c.containsVertex(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2])) {
					c.addVertex(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
				}

				if (!c.containsEdge(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
						link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2])) {
					if (link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] >= 0) {
						c.addEdge(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
								link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], link);
					} else {
						// System.err.println(String.format("checkContactMaps: Edge of %d duration not
						// added.",
						// link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]));
					}
				} else {
					// Edge already existed
					Integer[] e = c.getEdge(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
							link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
					int lastRelStart = e[e.length - 2].intValue();
					if (link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] < 0) {
						e[e.length - 1] = getGlobalTime() - lastRelStart;
					} else {
						if (link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] > 0) {
							// Extend edge
							c.removeEdge(e);
							e = Arrays.copyOf(e, e.length + 2);
							e[e.length - 1] = link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];
							e[e.length - 2] = link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
							c.addEdge(link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
									link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], e);
						} else {
							// Update edge
							e[e.length - 1] = Math.max(e[e.length - 1], getGlobalTime() - lastRelStart);

						}
					}

				}

			}
		}

	}

	protected void removeSingleRelationship(Person_Bridging_Pop breakRelPerson) {

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
				link = Arrays.copyOf(link, Abstract_Runnable_ClusterModel.LENGTH_CONTACT_MAP_EDGE);
				link[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] = -1;
				checkContactMaps(link, new ContactMap[] { cMapAll });

				tarMap.removeEdge(relRemove);

			}
		}

	}

	public File getBaseDir() {
		return baseDir;
	}

	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
	}

	public void setPrintStatus(PrintStream[] printStatus) {
		this.printStatus = printStatus;
	}

	public HashMap<String, Object> getStepwise_output() {
		return stepwise_output;
	}

	public void setStepwise_output(HashMap<String, Object> stepwise_output) {
		this.stepwise_output = stepwise_output;
	}

	public static final Integer[][] orderedEdgesFromContactMap(ContactMap cm) {
		Integer[][] edge_list = cm.edgeSet().toArray(new Integer[cm.edgeSet().size()][]);
		Arrays.sort(edge_list, new Comparator<Integer[]>() {
			@Override
			public int compare(Integer[] o1, Integer[] o2) {
				int cmp = Integer.compare(o1[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME],
						o2[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
				if (cmp == 0) {
					cmp = Integer.compare(o1[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION],
							o2[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]);
					if (cmp == 0) {
						cmp = Integer.compare(o1[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
								o2[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1]);
						if (cmp == 0) {
							cmp = Integer.compare(o1[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2],
									o2[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2]);
						}
					}
				}
				return cmp;
			}
		});

		return edge_list;
	}

}
