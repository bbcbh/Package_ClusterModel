package population;

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;

import availability.AbstractAvailability;
import person.AbstractIndividualInterface;
import population.availability.Bridging_Population_Availability;
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
			// FIELD_PARTNER_TYPE_PROB
			// float[GENDER]{CUMUL_REG_ONLY, CUMUL_CAS_ONLY}
			// Default:
			// Lee E, Mao L, Broady T, et al. Gay Community Periodic Survey: Melbourne 2018.
			// Sydney: Centre for Social Research in Health, UNSW Sydney, 2018.
			new float[][] { { 1, 1 }, { 1, 1 }, { 0.33f, 0.33f + 0.26f }, { 0.33f, 0.33f + 0.26f } },
			// FIELD_MAX_CASUAL_PARTNER_NUM_PROB
			// float[GENDER]{CUMUL_PROB_1, CUMUL_PROB_2, ... , MAX_PARTNER_1 , MAX_PARTNER_2
			// ...}}
			// Default:
			// National Centre in HIV Epidemiology and Clinical Research. Phase A of the
			// National Gay Men’s Syphilis Action Plan:
			// Modelling evidence and research on acceptability of interventions for
			// controlling syphilis in Australia.
			// Sydney: National Centre in HIV Epidemiology and Clinical Research, 2009.
			// Fogarty A, Mao L, Zablotska Manos I, et al.
			// The Health in Men and Positive Health cohorts: A comparison of trends in the
			// health and sexual behaviour of
			// HIV-negative and HIV-positive gay men, 2002-2005. Sydney: National Centre in
			// HIV Social Research, 2006.
			new float[][] { { 1, 0 }, { 1, 0 }, { 0.51f, 10 }, { 0.51f, 10 } },

	};

	public static final int FIELD_POP_COMPOSITION = AbstractFieldsArrayPopulation.LENGTH_FIELDS;
	public static final int FIELD_CONTACT_MAP = FIELD_POP_COMPOSITION + 1;
	public static final int FIELD_PARTNER_TYPE_PROB = FIELD_CONTACT_MAP + 1;
	public static final int FIELD_MAX_CASUAL_PARTNER_NUM_PROB = FIELD_PARTNER_TYPE_PROB + 1;
	public static final int LENGTH_FIELDS_BRIDGING_POP = FIELD_PARTNER_TYPE_PROB + 1;

	public static final int CONTACT_MAP_ALL = 0;
	public static final int CONTACT_MAP_HETRO = CONTACT_MAP_ALL + 1;
	public static final int CONTACT_MAP_MSM = CONTACT_MAP_HETRO + 1;

	public static final int RELMAP_HETRO = 0;
	public static final int RELMAP_MSM = 1;
	public static final int RELMAP_TOTAL = 2;

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
		float[][] max_casual_partners = (float[][]) getFields()[FIELD_MAX_CASUAL_PARTNER_NUM_PROB];

		for (int i = 0; i < popSizes.length; i++) {

			// Shared among same gender
			float[] p_type_prob = partner_type[i];

			AbstractIntegerDistribution dist = null;

			for (int g = 0; g < popSizes[i]; g++) {
				pop[popPt] = new Person_Bridging_Pop(popPt + 1, i, 18 * AbstractIndividualInterface.ONE_YEAR_INT,
						getGlobalTime(), // Initial age - might not be used.
						3); // 3 sites

				// TODO Set individual behavior

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
						// float[GENDER]{CUMUL_PROB, MAX_PARTNER}}
						float[] p_casual_prob = max_casual_partners[i];

						final float targetCDF = p_casual_prob[0];
						final int targetVal = (int) p_casual_prob[1];
						final float init_lamda = targetVal;

						final double RELATIVE_TOLERANCE = 0.005;
						final double ABSOLUTE_TOLERANCE = 0.001;
						final BrentOptimizer OPTIMIZER = new BrentOptimizer(RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE);
						final SearchInterval interval = new SearchInterval(1, 100, init_lamda);
						final UnivariateObjectiveFunction func = new UnivariateObjectiveFunction(
								new UnivariateFunction() {
									@Override
									public double value(double x) {
										PoissonDistribution func_dist = new PoissonDistribution(x);
										return func_dist.cumulativeProbability(targetVal) - targetCDF;
									}
								});

						double bestFit = OPTIMIZER.optimize(func, GoalType.MINIMIZE, interval).getPoint();
						dist = new PoissonDistribution(bestFit);

					}

					person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_6_MONTHS),
							dist.sample());

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

	}

	@Override
	public void advanceTimeStep(int deltaT) {
		// TODO Auto-generated method stub

	}

}
