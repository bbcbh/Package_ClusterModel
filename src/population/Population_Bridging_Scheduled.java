package population;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;

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
	public static final int SCHEDULE_PARTNERSHIP_TYPE = SCHEDULE_PARTNERSHIP_P2 + 1;
	public static final int LENGTH_SCHEDULE_PARTNERSHIP = SCHEDULE_PARTNERSHIP_TYPE + 1;

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
	protected void initialiseNumPartners(Person_Bridging_Pop person, int gender_grp_index, int gender_grp_count,
			AbstractIntegerDistribution casual_partner_dist_by_gender_grp) {

		
		int[] popSizes = (int[]) getFields()[FIELD_POP_COMPOSITION];
		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		int offset = numCat + numCat * gender_grp_index;

		int current_cat_probability_index = offset;
		float prob = ((float) gender_grp_count) / popSizes[gender_grp_index];
		float cumul_prob = field_mean_number_partner[current_cat_probability_index];

		while (cumul_prob < prob && current_cat_probability_index < offset + numCat) {
			current_cat_probability_index++;
			cumul_prob += field_mean_number_partner[current_cat_probability_index];
		}

		current_cat_probability_index -= offset;

		// Group 0 with 1 (for now)		
		int num_12_months = Math.max(1, (int) field_mean_number_partner[current_cat_probability_index]); 
		int range = 0;
		if (current_cat_probability_index > 0 ) {
			range = num_12_months - (int) field_mean_number_partner[current_cat_probability_index - 1];
		}
		if (range > 2) {
			num_12_months -= getRNG().nextInt(range);
		}

		// Seek regular
		boolean seekReg = !person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS)).equals(0);

		// Seek casual
		boolean seekCas = !person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS)).equals(0);

		if (seekReg && seekCas) {
			person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS), num_12_months);
			person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS), num_12_months);
		} else if (seekReg) {
			person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS),
					num_12_months);
		} else {
			person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS),
					num_12_months);
		}

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

			int[] schedule_current_length = new int[population_num_partner_in_last_12_months.length];
			int[][] schedule_current_cat_id_arr = new int[population_num_partner_in_last_12_months.length][getPop().length];
			HashMap<Integer, Integer> schedule_current_id_cat_map = new HashMap<>(); // Id, catergories

			// 0: ANY, 1: REG , 2: CAS
			final int SOUGHT_ANY = 0;
			final int SOUGHT_REG = SOUGHT_ANY + 1;
			final int SOUGHT_CAS = SOUGHT_REG + 1;
			final int LENGTH_SOUGHT = SOUGHT_CAS + 1;

			// Set to candidates arrays
			int[][][] candidates = new int[LENGTH_SOUGHT][getPop().length][];
			ComparatorByPartnershipSought[] comparators = new ComparatorByPartnershipSought[] {
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY),
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG),
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS) };

			int[] gender_end = new int[LENGTH_GENDER];

			int pI = 0;

			for (AbstractIndividualInterface absPerson : getPop()) {
				Person_Bridging_Pop person = (Person_Bridging_Pop) absPerson;
				int[] rc = getNumPartnerSought(person);
				int[] ent = new int[ComparatorByPartnershipSought.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];
				ent[ComparatorByPartnershipSought.INDEX_ID] = person.getId();
				ent[ComparatorByPartnershipSought.INDEX_GENDER] = person.getGenderType();
				ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] = Math.max(rc[0],rc[1]);
				ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG] = rc[0];
				ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS] = rc[1];
				ent[ComparatorByPartnershipSought.INDEX_MAX_REG] = rc[2];
				ent[ComparatorByPartnershipSought.INDEX_MAX_CAS] = rc[3];

				for (int a = 0; a < candidates.length; a++) {
					candidates[a][pI] = ent;
				}
				gender_end[person.getGenderType()]++;

				pI++;
			}

			for (int a = 0; a < candidates.length; a++) {
				Arrays.sort(candidates[a], comparators[a]);
			}

			int[] pop_diff_num_partner_12_months = cal_pop_diff_num_partner(population_num_partner_in_last_12_months);
			for (int c = numCat - 1; c >= 0; c++) {
				for (int g = 0; g < LENGTH_GENDER; g++) {
					int g_start = g > 0 ? gender_end[g - 1] : 0;
					int g_end = gender_end[g];
					int pdIndex = numCat + g * numCat + c;

					while (pop_diff_num_partner_12_months[pdIndex] > schedule_current_length[pdIndex]) {

						int numPartToSought_min = (int) cat_value[c];
						int numPartToSought_max = (c + 1 < cat_value.length) ? (int) cat_value[c + 1] :							
						(int) ((float[])getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS])[numCat-1];

						int[] numPartToSoughtAdj = new int[] {
								numPartToSought_min + getRNG().nextInt(numPartToSought_max - numPartToSought_min),
								numPartToSought_min + getRNG().nextInt(numPartToSought_max - numPartToSought_min), };
						Arrays.sort(numPartToSoughtAdj);

						int[] src_index = new int[2];

						int[] binary_key = new int[ComparatorByPartnershipSought.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];
						binary_key[ComparatorByPartnershipSought.INDEX_GENDER] = g;

						while (src_index[1] == src_index[0] && numPartToSoughtAdj != null) {

							binary_key[ComparatorByPartnershipSought.INDEX_ID] = -1;
							binary_key[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] = numPartToSoughtAdj[0];
							src_index[0] = ~Arrays.binarySearch(candidates[0], g_start, g_end, binary_key,
									comparators[SOUGHT_ANY]);

							binary_key[ComparatorByPartnershipSought.INDEX_ID] = Integer.MAX_VALUE;
							binary_key[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] = numPartToSoughtAdj[1];
							src_index[1] = ~Arrays.binarySearch(candidates[0], g_start, g_end, binary_key,
									comparators[SOUGHT_ANY]);

							numPartToSoughtAdj[0]--;
							numPartToSoughtAdj[1]++;

							if (numPartToSoughtAdj[0] < numPartToSought_min
									&& numPartToSoughtAdj[1] > numPartToSought_max) {
								numPartToSoughtAdj = null;
							} else {
								numPartToSoughtAdj[0] = Math.max(numPartToSought_min, numPartToSoughtAdj[0]);
								numPartToSoughtAdj[1] = Math.min(numPartToSought_max, numPartToSoughtAdj[1]);
							}
						}

						int src_candidate_index = src_index[0];
						if (src_index[1] > src_index[0]) {
							src_candidate_index += getRNG().nextInt(src_index[1] - src_index[0]);
						}

						int[] src_candidate_cmp_ent = candidates[SOUGHT_ANY][src_candidate_index];
						boolean onlySoughtReg = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_MAX_CAS] == 0;
						boolean onlySoughtCas = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_MAX_REG] == 0;
						int maxToSought = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY];

						int[] tar_possible_gender;

						switch (src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_GENDER]) {
						case GENDER_FEMALE:
							tar_possible_gender = new int[] { GENDER_HETRO_MALE, GENDER_MSMW };
							break;
						case GENDER_HETRO_MALE:
							tar_possible_gender = new int[] { GENDER_FEMALE };
							break;
						case GENDER_MSMO:
							tar_possible_gender = new int[] { GENDER_MSMO, GENDER_MSMW };
							break;
						default:
							tar_possible_gender = new int[] { GENDER_FEMALE, GENDER_MSMO, GENDER_MSMW };
						}

						final int tar_sough_partner_type_index;
						final int tar_binary_key_num_sought_partner_type_index;

						if (onlySoughtReg) {
							tar_sough_partner_type_index = SOUGHT_REG;
							tar_binary_key_num_sought_partner_type_index = ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG;

						} else if (onlySoughtCas) {
							tar_sough_partner_type_index = SOUGHT_CAS;
							tar_binary_key_num_sought_partner_type_index = ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS;
						} else {
							tar_sough_partner_type_index = SOUGHT_ANY;
							tar_binary_key_num_sought_partner_type_index = ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY;
						}

						int[][] candidateRangeByGender = new int[tar_possible_gender.length][2];
						int[] totalLength = new int[tar_possible_gender.length];

						Arrays.fill(binary_key, 0);

						for (int cG = 0; cG < candidateRangeByGender.length; cG++) {
							int cg_start = cG > 0 ? gender_end[cG - 1] : 0;
							int cg_end = gender_end[cG];

							// At least seek one partner
							binary_key[ComparatorByPartnershipSought.INDEX_GENDER] = cG;
							binary_key[tar_binary_key_num_sought_partner_type_index] = 1;
							binary_key[ComparatorByPartnershipSought.INDEX_ID] = -1;

							candidateRangeByGender[cG][0] = ~Arrays.binarySearch(
									candidates[tar_sough_partner_type_index], cg_start, cg_end, binary_key,
									comparators[tar_sough_partner_type_index]);

							binary_key[tar_binary_key_num_sought_partner_type_index] = Integer.MAX_VALUE;
							binary_key[ComparatorByPartnershipSought.INDEX_ID] = Integer.MAX_VALUE;

							candidateRangeByGender[cG][1] = ~Arrays.binarySearch(
									candidates[tar_sough_partner_type_index], cg_start, cg_end, binary_key,
									comparators[tar_sough_partner_type_index]);

							totalLength[cG] = candidateRangeByGender[cG][1] - candidateRangeByGender[cG][0];

							if (cG > 0) {
								totalLength[cG] += totalLength[cG - 1];
							}
						}

						int[][] target_candidates = new int[totalLength[totalLength.length]][];

						int offset = 0;
						for (int cG = 0; cG < candidateRangeByGender.length; cG++) {
							System.arraycopy(candidates[tar_sough_partner_type_index], candidateRangeByGender[cG][0],
									target_candidates, offset, totalLength[cG]);
							offset += totalLength[cG];
						}

						Arrays.sort(target_candidates, new Comparator<int[]>() {
							@Override
							public int compare(int[] o1, int[] o2) {
								return Integer.compare(o1[tar_binary_key_num_sought_partner_type_index],
										o2[tar_binary_key_num_sought_partner_type_index]);
							}
						});

						int[] cumul_weight = new int[target_candidates.length];
						cumul_weight[0] = target_candidates[0][tar_binary_key_num_sought_partner_type_index];
						for (int i = 1; i < target_candidates.length; i++) {
							cumul_weight[i] = cumul_weight[i - 1]
									+ target_candidates[i][tar_binary_key_num_sought_partner_type_index];
						}

						// Choose partnership based on number sought
						int numTar = Math.min(maxToSought, target_candidates.length);
						int[][] partnered_with = new int[numTar][];

						for (int i = 0; i < cumul_weight[cumul_weight.length - 1] && numTar > 0; i++) {
							if (getRNG().nextInt(cumul_weight[cumul_weight.length - 1] - i) < numTar) {
								int tI = Arrays.binarySearch(cumul_weight, i);
								if (tI < 0) {
									tI = ~tI;
								}
								partnered_with[partnered_with.length - numTar] = target_candidates[tI];
								numTar--;
							}
						}

						for (int i = 0; i < partnered_with.length; i++) {

							int[] tar_candidate_cmp_ent = partnered_with[i];

							Integer partner_form_time = getGlobalTime()
									+ getRNG().nextInt(AbstractIndividualInterface.ONE_YEAR_INT);

							// Determine if it can be a regular or casual partnership
							int partnership_type = tar_sough_partner_type_index;
							if (tar_sough_partner_type_index == SOUGHT_ANY) {
								if (tar_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG] > 0) {
									partnership_type = SOUGHT_REG;
								} else {
									partnership_type = SOUGHT_CAS;
								}
							}

							Integer[] schedule_partnership_ent = new Integer[LENGTH_SCHEDULE_PARTNERSHIP];

							schedule_partnership_ent[SCHEDULE_PARTNERSHIP_P1] = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_ID];
							schedule_partnership_ent[SCHEDULE_PARTNERSHIP_P2] = tar_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_ID];
							schedule_partnership_ent[SCHEDULE_PARTNERSHIP_TYPE] = partnership_type;

							ArrayList<Integer[]> partnerships = schedule_partnership.get(partner_form_time);
							if (partnerships == null) {
								partnerships = new ArrayList<>();
								schedule_partnership.put(partner_form_time, partnerships);
							}
							partnerships.add(schedule_partnership_ent);

							tar_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY]--;

							switch (tar_sough_partner_type_index) {
							case SOUGHT_REG:
								tar_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG]--;
							case SOUGHT_CAS:
								tar_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS]--;
							}

							// TODO: Current progress update schedule_current_length

						}

						// Update schedule current

						int adjC = Arrays.binarySearch(cat_value, partnered_with.length);
						if (adjC < 0) {
							adjC = ~adjC;
						}
						int adjpdIndex = numCat + g * numCat + adjC;

						int src_id = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_ID];

						schedule_current_cat_id_arr[adjpdIndex][schedule_current_length[adjpdIndex]] = src_id;
						schedule_current_length[adjpdIndex]++;

						Arrays.sort(schedule_current_cat_id_arr[adjpdIndex], 0, schedule_current_length[adjpdIndex]);
						schedule_current_id_cat_map.put(src_id, adjpdIndex);

						// Update candidate arrays
						for (int s = 0; s < LENGTH_SOUGHT; s++) {
							Arrays.sort(candidates[s], comparators[s]);
						}

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

				// For consistency
				if (pair[1].getGenderType() == Person_Bridging_Pop.GENDER_TYPE_FEMALE
						|| pair[1].getId() < pair[0].getId()) {
					Person_Bridging_Pop temp = pair[0];
					pair[0] = pair[1];
					pair[1] = temp;
				}

				int partnershipType = edge[SCHEDULE_PARTNERSHIP_TYPE];

				// TODO: create partnership based on type

				if (partnershipType > 1) { // Regular
					int mapType;
					if (pair[1].getGenderType() == GENDER_FEMALE || pair[0].getGenderType() == GENDER_FEMALE) {
						mapType = RELMAP_HETRO;
					} else {
						mapType = RELMAP_MSM;
					}

					// Assume no concurrency
					for (int p = 0; p < pair.length; p++) {
						removeSingleRelationship(pair[p]);
					}

					formRelationship(pair, getRelMap()[mapType], partnershipType, mapType);

				} else { // Casual
					ContactMap[] cMaps = new ContactMap[] {
							((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL] };

					// Add casual partner
					for (int p = 0; p < pair.length; p++) {
						pair[p].addCasualPartner(pair[(p + 1) % 2]);
					}

					checkContactMaps(
							new Integer[] { pair[0].getId(), pair[1].getId(), getGlobalTime(), partnershipType },
							cMaps);
				}

			}
		}
	}

	private int[] getNumPartnerSought(Person_Bridging_Pop person) {
		int[] rc = new int[4]; // 0 - 1 : To be sought, 2 - 3: Max in 12 month
		rc[2] = (int) person.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS));
		rc[3] = (int) person.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));

		rc[0] = Math.max(0, rc[2] - getNumRegularPartnersCurrently(person));
		rc[1] = Math.max(0, rc[3] - person.getNumCasualInRecord());

		return rc;
	}

	private static final class ComparatorByPartnershipSought implements Comparator<int[]> {

		private static int INDEX_ID = 0;
		private static int INDEX_GENDER = INDEX_ID + 1;
		private static int INDEX_NUM_TO_SOUGHT_ANY = INDEX_GENDER + 1;
		private static int INDEX_NUM_TO_SOUGHT_REG = INDEX_NUM_TO_SOUGHT_ANY + 1;
		private static int INDEX_NUM_TO_SOUGHT_CAS = INDEX_NUM_TO_SOUGHT_REG + 1;
		private static int INDEX_MAX_REG = INDEX_NUM_TO_SOUGHT_CAS + 1;
		private static int INDEX_MAX_CAS = INDEX_MAX_REG + 1;
		private static int LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT = INDEX_MAX_CAS + 1;

		final int cmpMethod;

		public ComparatorByPartnershipSought(int cmpMethod) {
			super();
			this.cmpMethod = cmpMethod;
		}

		@Override
		public int compare(int[] o1, int[] o2) {
			int cmp;
			cmp = Float.compare(o1[INDEX_GENDER], o2[INDEX_GENDER]);
			if (cmp == 0) {
				cmp = Float.compare(o1[cmpMethod], o2[cmpMethod]);
				if (cmp == 0) {
					cmp = Float.compare(o1[INDEX_ID], o2[INDEX_ID]);
				}
			}
			return cmp;
		}
	}
}
