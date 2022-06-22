package population;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

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

	protected transient int lastPartnershipScheduling = 0;

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

		boolean reqPartnerScheduling = getGlobalTime() == lastPartnershipScheduling
				+ AbstractIndividualInterface.ONE_YEAR_INT;

		if (reqPartnerScheduling) {
			lastPartnershipScheduling = getGlobalTime();

			ArrayList<Person_Bridging_Pop> candidates_any = new ArrayList<>(getPop().length);
			ArrayList<Person_Bridging_Pop> candidates_reg = new ArrayList<>(getPop().length);
			ArrayList<Person_Bridging_Pop> candidates_cas = new ArrayList<>(getPop().length);

			int[] gender_count_any = new int[LENGTH_GENDER];
			int[] gender_count_reg = new int[LENGTH_GENDER];
			int[] gender_count_cas = new int[LENGTH_GENDER];

			HashMap<Integer, int[]> rel_status = new HashMap<>();

			for (AbstractIndividualInterface absPerson : getPop()) {
				Person_Bridging_Pop person = (Person_Bridging_Pop) absPerson;
				int[] rc = getNumPartnerSought(person);

				if (rc[0] + rc[1] > 0) {
					candidates_any.add(person);
					gender_count_any[person.getGenderType()]++;
					rel_status.put(person.getId(), rc);
				}
				if (rc[0] > 0) {
					candidates_reg.add(person);
					gender_count_reg[person.getGenderType()]++;
				}
				if (rc[1] > 0) {
					candidates_cas.add(person);
					gender_count_cas[person.getGenderType()]++;
				}
			}

			candidates_any.sort(COMPARATOR_BY_PARTNERSHIP_SOUGHT);
			candidates_reg.sort(COMPARATOR_BY_PARTNERSHIP_SOUGHT);
			candidates_cas.sort(COMPARATOR_BY_PARTNERSHIP_SOUGHT);

			// TODO: Current progress - Form partnership on a 12 months schedule

			int[] pop_diff_num_partner_12_months = cal_pop_diff_num_partner(population_num_partner_in_last_12_months);
			int numCat = pop_diff_num_partner_12_months.length / (1 + LENGTH_GENDER);

			for (int g = 0; g < LENGTH_GENDER; g++) {
				for (int c = numCat - 1; c >= 0; c++) {
					int pdIndex = numCat + g * numCat + c;
					int numToSeek = pop_diff_num_partner_12_months[pdIndex];

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

	private final Comparator<AbstractIndividualInterface> COMPARATOR_BY_PARTNERSHIP_SOUGHT = new Comparator<AbstractIndividualInterface>() {
		@Override
		public int compare(AbstractIndividualInterface o1, AbstractIndividualInterface o2) {
			int cmp;
			if (o1 instanceof Person_Bridging_Pop && o2 instanceof Person_Bridging_Pop) {
				cmp = Integer.compare(((Person_Bridging_Pop) o1).getGenderType(),
						((Person_Bridging_Pop) o2).getGenderType());
				if (cmp == 0) {
					int[] rc1 = getNumPartnerSought((Person_Bridging_Pop) o1);
					int[] rc2 = getNumPartnerSought((Person_Bridging_Pop) o2);
					cmp = Integer.compare(rc1[0] + rc1[1], rc2[0] + rc2[1]);
					if (cmp == 0) {
						cmp = Integer.compare(o1.getId(), o2.getId());
					}
				}
			} else {
				int g1 = o1.isMale() ? 1 : 0;
				int g2 = o2.isMale() ? 1 : 0;
				cmp = Integer.compare(g1, g2);
				if (cmp == 0) {
					cmp = Integer.compare(o1.getId(), o2.getId());
				}
			}
			return cmp;
		}

	};

}
