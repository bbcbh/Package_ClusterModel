package population;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

import infection.AbstractInfection;
import person.AbstractIndividualInterface;
import population.person.Person_Bridging_Pop;
import relationship.ContactMap;

public class Population_Bridging_Scheduled extends Population_Bridging {
	/**
	 * 
	 */
	private static final long serialVersionUID = 7076023161972931764L;

	protected HashMap<Integer, ArrayList<Integer[]>> schedule_partnership;

	protected transient int lastPartnershipScheduling = -1;

	public static final int SCHEDULE_PARTNERSHIP_P1 = 0;
	public static final int SCHEDULE_PARTNERSHIP_P2 = SCHEDULE_PARTNERSHIP_P1 + 1;
	public static final int SCHEDULE_PARTNERSHIP_DUR = SCHEDULE_PARTNERSHIP_P2 + 1;
	public static final int LENGTH_SCHEDULE_PARTNERSHIP = SCHEDULE_PARTNERSHIP_DUR + 1;

	public Population_Bridging_Scheduled(long seed) {
		super(seed);
		schedule_partnership = new HashMap<>();
		Object[] newFields = Arrays.copyOf(super.getFields(), super.getFields().length + 1);
		newFields[newFields.length - 1] = schedule_partnership;
		super.setFields(newFields);
	}

	@Override
	protected void initialiseTransientFields() {
		super.initialiseTransientFields();
		lastPartnershipScheduling = 0;
	}

	@SuppressWarnings("unchecked")
	public static Population_Bridging_Scheduled decodeFromStream(java.io.ObjectInputStream inStr)
			throws IOException, ClassNotFoundException {

		int globalTime = inStr.readInt();
		AbstractInfection[] infList = (AbstractInfection[]) inStr.readObject();
		Object[] decoded_fields = (Object[]) inStr.readObject();

		Population_Bridging_Scheduled pop = new Population_Bridging_Scheduled((long) decoded_fields[0]);

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
		pop.schedule_partnership = (HashMap<Integer, ArrayList<Integer[]>>) (decoded_fields[decoded_fields.length - 1]);
		return pop;
	}

	@Override

	public void advanceTimeStep(int deltaT) {
		incrementTime(deltaT);

		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int[] population_num_partner_in_last_12_months = new int[field_mean_number_partner.length];
		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		boolean reportPartnerStat = getGlobalTime() == lastPartnershipScheduling
				+ AbstractIndividualInterface.ONE_YEAR_INT;

		reportPartnerStat = true; // TODO: Debug for advanceTimeStep (Schedule)

		for (AbstractIndividualInterface p : this.getPop()) {
			Person_Bridging_Pop person = (Person_Bridging_Pop) p;
			person.incrementTime(deltaT, getInfList());

			if (reportPartnerStat) {

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
			}
		}

		if (reportPartnerStat) {
			if (printStatus != null) {
				if (population_num_partner_in_last_12_months.length == LENGTH_GENDER) {
					printStatus.printf("# partners in last 12 months = %s - Day %d",
							Arrays.toString(population_num_partner_in_last_12_months), getGlobalTime());
				} else {
					StringWriter wri = new StringWriter();
					PrintWriter pri = new PrintWriter(wri);
					pri.println(String.format("# partners in last 12 months -  Day %d", getGlobalTime()));
					for (int g = 0; g < LENGTH_GENDER; g++) {
						pri.printf(" %d: %s", g,
								Arrays.toString(Arrays.copyOfRange(population_num_partner_in_last_12_months,
										numCat + g * numCat, 2 * numCat + g * numCat)));
					}
					pri.close();
					printStatus.println(wri.toString());

				}

			}
		}

		updateRelRelationshipMap();

		if (stepwise_output != null) {
			stepwise_output.put(STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS, population_num_partner_in_last_12_months);
		}
		formPartnerships(population_num_partner_in_last_12_months);

	}

	@Override
	public void formPartnerships(int[] population_num_partner_in_last_12_months) {

		boolean reqPartnerScheduling = (lastPartnershipScheduling < 0)
				|| (getGlobalTime() == lastPartnershipScheduling + AbstractIndividualInterface.ONE_YEAR_INT);

		if (reqPartnerScheduling) {
			lastPartnershipScheduling = getGlobalTime();

			int numCat = population_num_partner_in_last_12_months.length / (1 + LENGTH_GENDER);
			float[] cat_value = Arrays.copyOf((float[]) getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS], numCat);

			int[][] candidates_any = new int[getPop().length][LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];

			int[] gender_end = new int[LENGTH_GENDER];

			int pI = 0;

			for (AbstractIndividualInterface absPerson : getPop()) {
				Person_Bridging_Pop person = (Person_Bridging_Pop) absPerson;
				int[] rc = getNumPartnerSought(person);
				candidates_any[pI][COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID] = person.getId();
				candidates_any[pI][COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_GENDER] = person.getGenderType();
				candidates_any[pI][COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_NUM_TO_SOUGHT] = rc[0] + rc[1];
				gender_end[person.getGenderType()]++;
				pI++;
			}

			Arrays.sort(candidates_any, COMPARATOR_BY_PARTNERSHIP_SOUGHT);
		
			HashMap<Integer, int[]> numSoughtThisRound = new HashMap<>();

			

			int[] pop_diff_num_partner_12_months = cal_pop_diff_num_partner(population_num_partner_in_last_12_months);
			for (int c = numCat - 1; c >= 0; c++) {
				for (int g = 0; g < LENGTH_GENDER; g++) {
					int g_start = g > 0 ? gender_end[g - 1] : 0;
					int g_end = gender_end[g];
					int pdIndex = numCat + g * numCat + c;
					while (pop_diff_num_partner_12_months[pdIndex] > 0) {
						int numPartToSought = (int) cat_value[c];
						int maxPart = (c + 1 < cat_value.length) ? (int) cat_value[c + 1] : 60;
						numPartToSought += getRNG().nextInt(maxPart - numPartToSought);

						int src_start = ~Arrays.binarySearch(candidates_any, g_start, g_end,
								new int[] { -1, g, numPartToSought }, COMPARATOR_BY_PARTNERSHIP_SOUGHT);

						int src_end = ~Arrays.binarySearch(candidates_any, g_start, g_end,
								new int[] { Integer.MAX_VALUE, g, numPartToSought }, COMPARATOR_BY_PARTNERSHIP_SOUGHT);
						
						int candidate_index = src_start + getRNG().nextInt(src_end-src_start);
						int[] cmp_ent = candidates_any[candidate_index];
						
						Person_Bridging_Pop src_person = (Person_Bridging_Pop)
								getLocalData().get(cmp_ent[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID]);															
						
						int[] numSoughted_src = numSoughtThisRound.get(cmp_ent[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID]);
						
						if(numSoughted_src == null) {
							numSoughted_src = new int[1];
							numSoughtThisRound.put(cmp_ent[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID], numSoughted_src);
						}

						while (numPartToSought > 0) {
							
							// TODO: Current progress - Form partnership on a 12 months schedule
							
							
							
							numSoughted_src[0]++;

							numPartToSought--;
						}
						
						
						int adjC = Arrays.binarySearch(cat_value, numSoughted_src[0]);
						if(adjC < 0) {
							adjC = ~adjC;
						}																		
						int adjpdIndex =  numCat + g * numCat + adjC;
						

						pop_diff_num_partner_12_months[adjpdIndex]--;
						
					}

				}

			}
		}

		ArrayList<Integer[]> partnerships = schedule_partnership.remove(getGlobalTime());

		// TODO: Check actual formation of partnership
		if (partnerships != null) {
			for (Integer[] edge : partnerships) {
				Person_Bridging_Pop[] pair = new Person_Bridging_Pop[] {
						(Person_Bridging_Pop) getLocalData().get(edge[SCHEDULE_PARTNERSHIP_P1]),
						(Person_Bridging_Pop) getLocalData().get(edge[SCHEDULE_PARTNERSHIP_P2]), };

				int duration = edge[SCHEDULE_PARTNERSHIP_DUR];

				if (duration > 1) { // Regular
					int mapType;
					if (pair[1].getGenderType() == GENDER_HETRO_FEMALE
							|| pair[0].getGenderType() == GENDER_HETRO_FEMALE) {
						mapType = RELMAP_HETRO;
					} else {
						mapType = RELMAP_MSM;
					}

					// Assume no concurrency
					for (int p = 0; p < pair.length; p++) {
						removeSingleRelationship(pair[p]);
					}

					formRelationship(pair, getRelMap()[mapType], duration, mapType);

				} else { // Casual
					ContactMap[] cMaps = new ContactMap[] {
							((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL] };

					// Add casual partner
					for (int p = 0; p < pair.length; p++) {
						pair[p].addCasualPartner(pair[(p + 1) % 2]);
					}

					checkContactMaps(new Integer[] { pair[0].getId(), pair[1].getId(), getGlobalTime(), duration },
							cMaps);
				}

			}
		}
	}

	private int[] getNumPartnerSought(Person_Bridging_Pop person) {
		int[] rc = new int[2];
		rc[0] = (int) person.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS));
		rc[1] = (int) person.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));

		rc[0] = Math.max(0, rc[0] - getNumRegularPartnersCurrently(person));
		rc[1] = Math.max(0, rc[1] - person.getNumCasualInRecord());
		return rc;
	}

	private static int COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID = 0;
	private static int COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_GENDER = COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID;
	private static int COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_NUM_TO_SOUGHT = COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_GENDER
			+ 1;
	private static int LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT = COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_NUM_TO_SOUGHT
			+ 1;

	private final Comparator<int[]> COMPARATOR_BY_PARTNERSHIP_SOUGHT = new Comparator<int[]>() {

		@Override
		public int compare(int[] o1, int[] o2) {
			int cmp;

			cmp = Float.compare(o1[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_GENDER],
					o2[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_GENDER]);
			if (cmp == 0) {
				cmp = Float.compare(o1[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_NUM_TO_SOUGHT],
						o2[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_NUM_TO_SOUGHT]);
				if (cmp == 0) {
					cmp = Float.compare(o1[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID],
							o2[COMPARATOR_BY_PARTNERSHIP_SOUGHT_INDEX_ID]);
				}
			}

			return cmp;
		}

	};

}
