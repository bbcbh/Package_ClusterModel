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
import relationship.RelationshipMap;

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
		if (current_cat_probability_index > 0) {
			range = num_12_months - ((int) field_mean_number_partner[current_cat_probability_index - 1] + 1);
		}
		if (range > 1) {
			num_12_months -= getRNG().nextInt(range);
		}

		// Seek regular
		boolean seekReg = !person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS)).equals(0);

		// Seek casual
		boolean seekCas = !person
				.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS)).equals(0);

		if (seekReg && seekCas) {
			person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS),
					num_12_months);
			person.setParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS),
					num_12_months);
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

		boolean reportPartnerStat = lastPartnershipScheduling < AbstractIndividualInterface.ONE_YEAR_INT
				? (getGlobalTime() == AbstractIndividualInterface.ONE_YEAR_INT)
				: getGlobalTime() == lastPartnershipScheduling + AbstractIndividualInterface.ONE_YEAR_INT;

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
					printStatus.printf("# partners in last 12 months = %s - Day %d\n",
							Arrays.toString(population_num_partner_in_last_12_months), getGlobalTime());
				} else {
					StringWriter wri = new StringWriter();
					PrintWriter pri = new PrintWriter(wri);
					pri.println(String.format("# partners in last 12 months -  Day %d\n", getGlobalTime()));
					for (int g = 0; g < LENGTH_GENDER; g++) {
						pri.printf(" %d: %s\n", g,
								Arrays.toString(Arrays.copyOfRange(population_num_partner_in_last_12_months,
										numCat + g * numCat, 2 * numCat + g * numCat)));
					}
					printStatus.println(wri.toString());
					pri.close();

				}

			}
		}

		updateRelRelationshipMap();

		if (stepwise_output != null) {
			stepwise_output.put(STEPWISE_OUTPUT_NUM_PARTNERS_IN_12_MONTHS, population_num_partner_in_last_12_months);
		}
		formPartnerships(population_num_partner_in_last_12_months);

	}

	private static final int CANDIDATE_ARRAY_SCHEDULE_LIMIT = 0;
	private static final int CANDIDATE_ARRAY_SOUGHT_ANY = CANDIDATE_ARRAY_SCHEDULE_LIMIT + 1;
	private static final int CANDIDATE_ARRAY_SOUGHT_REG = CANDIDATE_ARRAY_SOUGHT_ANY + 1;
	private static final int CANDIDATE_ARRAY_SOUGHT_CAS = CANDIDATE_ARRAY_SOUGHT_REG + 1;
	private static final int CANDIDATE_ARRAY_CURRENTLY_HAS = CANDIDATE_ARRAY_SOUGHT_CAS + 1;

	private final boolean scheduling_debug = !true;

	@Override
	public void formPartnerships(int[] population_num_partner_in_last_12_months) {

		boolean reqPartnerScheduling = (lastPartnershipScheduling < 0)
				|| (getGlobalTime() == lastPartnershipScheduling + AbstractIndividualInterface.ONE_YEAR_INT);

		int schedule_range = AbstractIndividualInterface.ONE_YEAR_INT;

		if (reqPartnerScheduling) {

			lastPartnershipScheduling = getGlobalTime();

			int numCat = population_num_partner_in_last_12_months.length / (1 + LENGTH_GENDER);
			float[] categories_values = Arrays.copyOf((float[]) getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS],
					numCat);

			// Set to candidates arrays
			int[][] candidates = new int[getPop().length][];
			int[] gender_end = new int[LENGTH_GENDER];
			ComparatorByPartnershipSought[] comparators = new ComparatorByPartnershipSought[] {
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_SCHEDULE_LIMIT),
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY),
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG),
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS),
					new ComparatorByPartnershipSought(ComparatorByPartnershipSought.INDEX_CURRENT_ANY), };

			fillCandidateList(candidates, gender_end, schedule_range == 0);

			int[] pop_demand_num_partner_12_months = cal_pop_diff_num_partner(
					new int[population_num_partner_in_last_12_months.length]);
			int[] addressed_demand_so_far = new int[population_num_partner_in_last_12_months.length];
			int[] binary_key = new int[ComparatorByPartnershipSought.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];

			for (int src_cat_index = numCat - 1; src_cat_index > 0; src_cat_index--) {
				for (int src_gender = 0; src_gender < LENGTH_GENDER; src_gender++) {
					int src_gender_index_start = src_gender > 0 ? gender_end[src_gender - 1] : 0;
					int src_gender_index_end = gender_end[src_gender];
					int src_pdIndex = numCat + src_gender * numCat + src_cat_index;

					determine_src_candidate: while (pop_demand_num_partner_12_months[src_pdIndex] > addressed_demand_so_far[src_pdIndex]) {
						// i.e. has pop demand > 0

						int numPartToSought_max = (int) categories_values[src_cat_index];
						int numPartToSought_min = ((src_cat_index == 0) ? 0
								: (int) categories_values[src_cat_index - 1]) + 1;

						Arrays.fill(binary_key, 0);
						binary_key[ComparatorByPartnershipSought.INDEX_GENDER] = src_gender;

						int[] src_candidate_cmp_ent;

						if (schedule_range == 0) {
							Arrays.sort(candidates, src_gender_index_start, src_gender_index_end,
									comparators[CANDIDATE_ARRAY_CURRENTLY_HAS]);

							int[][] src_candidates = new int[0][];

							for (int src_sub_cat_index = src_cat_index
									- 1; src_sub_cat_index >= 0; src_sub_cat_index--) {
								int sub_cat_pdIndex = numCat + src_gender * numCat + src_sub_cat_index;
								if (pop_demand_num_partner_12_months[sub_cat_pdIndex] < addressed_demand_so_far[sub_cat_pdIndex]) {
									// i.e. has excess (i.e. demand < 0)

									int subCat_numPartToSought_min = (src_sub_cat_index == 0) ? 0
											: (int) categories_values[src_sub_cat_index - 1] + 1;

									int subCat_numPartToSought_max = (int) categories_values[src_sub_cat_index];

									int[] subCat_sel_range = new int[2];
									binary_key[ComparatorByPartnershipSought.INDEX_ID] = -1;
									binary_key[ComparatorByPartnershipSought.INDEX_CURRENT_ANY] = subCat_numPartToSought_min;

									subCat_sel_range[0] = ~Arrays.binarySearch(candidates, src_gender_index_start,
											src_gender_index_end, binary_key,
											comparators[CANDIDATE_ARRAY_CURRENTLY_HAS]);

									binary_key[ComparatorByPartnershipSought.INDEX_ID] = Integer.MAX_VALUE;
									binary_key[ComparatorByPartnershipSought.INDEX_CURRENT_ANY] = subCat_numPartToSought_max;

									subCat_sel_range[1] = ~Arrays.binarySearch(candidates, src_gender_index_start,
											src_gender_index_end, binary_key,
											comparators[CANDIDATE_ARRAY_CURRENTLY_HAS]);

									if (subCat_sel_range[1] > subCat_sel_range[0]) {

										int[][] withinCat = Arrays.copyOfRange(candidates, subCat_sel_range[0],
												subCat_sel_range[1]);

										Arrays.sort(withinCat, comparators[CANDIDATE_ARRAY_SOUGHT_ANY]);

										binary_key[ComparatorByPartnershipSought.INDEX_ID] = Integer.MAX_VALUE;
										binary_key[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] = 0;

										int seekingPartners = ~Arrays.binarySearch(withinCat, 0, withinCat.length,
												binary_key, comparators[CANDIDATE_ARRAY_SOUGHT_ANY]);

										binary_key[ComparatorByPartnershipSought.INDEX_ID] = -1;
										binary_key[ComparatorByPartnershipSought.INDEX_SCHEDULE_LIMIT] = numPartToSought_min;

										subCat_sel_range[0] = ~Arrays.binarySearch(withinCat, seekingPartners,
												withinCat.length, binary_key,
												comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);

										binary_key[ComparatorByPartnershipSought.INDEX_ID] = Integer.MAX_VALUE;
										binary_key[ComparatorByPartnershipSought.INDEX_SCHEDULE_LIMIT] = numPartToSought_max;

										subCat_sel_range[1] = ~Arrays.binarySearch(withinCat, 0, withinCat.length,
												binary_key, comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);

										if (subCat_sel_range[1] > subCat_sel_range[0]) {
											int src_last = src_candidates.length;
											src_candidates = Arrays.copyOf(src_candidates,
													subCat_sel_range[1] - subCat_sel_range[0]);
											System.arraycopy(withinCat, subCat_sel_range[0], src_candidates, src_last,
													subCat_sel_range[1] - subCat_sel_range[0]);
										}
									}
								}
							}
							if (src_candidates.length == 0) {
								// No more candidate that can sought required number of partners
								break determine_src_candidate;
							} else {
								src_candidate_cmp_ent = src_candidates[getRNG().nextInt(src_candidates.length)];
							}
						} else {
							// Non-daily version
							int[] numPartToSoughtAdj = new int[] { numPartToSought_min, numPartToSought_max, };
							int[] src_index_range = new int[2];

							Arrays.fill(binary_key, 0);
							binary_key[ComparatorByPartnershipSought.INDEX_GENDER] = src_gender;

							Arrays.sort(candidates, src_gender_index_start, src_gender_index_end,
									comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);

							while (src_index_range[1] == src_index_range[0] && numPartToSoughtAdj != null) {

								binary_key[ComparatorByPartnershipSought.INDEX_ID] = -1;
								binary_key[ComparatorByPartnershipSought.INDEX_SCHEDULE_LIMIT] = numPartToSoughtAdj[0];
								src_index_range[0] = ~Arrays.binarySearch(candidates, src_gender_index_start,
										src_gender_index_end, binary_key, comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);

								binary_key[ComparatorByPartnershipSought.INDEX_ID] = Integer.MAX_VALUE;
								binary_key[ComparatorByPartnershipSought.INDEX_SCHEDULE_LIMIT] = numPartToSoughtAdj[1];
								src_index_range[1] = ~Arrays.binarySearch(candidates, src_gender_index_start,
										src_gender_index_end, binary_key, comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);

								if (src_index_range[1] == src_index_range[0]) {
									// No more candidate that can sought required number of partners
									break determine_src_candidate;

								}
							}

							int src_candidate_index = src_index_range[0];
							int src_candidate_range = src_index_range[1] - src_index_range[0];
							if (src_candidate_range > 1) {
								src_candidate_index += getRNG().nextInt(src_candidate_range);
							}

							src_candidate_cmp_ent = candidates[src_candidate_index];
						}

						boolean onlySoughtReg = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_MAX_CAS] == 0;
						boolean onlySoughtCas = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_MAX_REG] == 0;
						int numPartnerToSought = src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY];

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

						final int tar_sought_candidate_partner_type_index;
						final int tar_sought_comparator_partner_type_index;

						if (onlySoughtReg) {
							tar_sought_candidate_partner_type_index = CANDIDATE_ARRAY_SOUGHT_REG;
							tar_sought_comparator_partner_type_index = ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG;
						} else if (onlySoughtCas) {
							tar_sought_candidate_partner_type_index = CANDIDATE_ARRAY_SOUGHT_CAS;
							tar_sought_comparator_partner_type_index = ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS;
						} else {
							tar_sought_candidate_partner_type_index = CANDIDATE_ARRAY_SOUGHT_ANY;
							tar_sought_comparator_partner_type_index = ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY;
						}

						int[][] candidateRangeByGender = new int[tar_possible_gender.length][2];
						int[] totalLength = new int[tar_possible_gender.length];

						Arrays.fill(binary_key, 0);

						for (int cG = 0; cG < tar_possible_gender.length; cG++) {
							int tar_gender = tar_possible_gender[cG];
							int cg_start = tar_gender > 0 ? gender_end[tar_gender - 1] : 0;
							int cg_end = gender_end[tar_gender];

							Arrays.sort(candidates, cg_start, cg_end,
									comparators[tar_sought_candidate_partner_type_index]);

							// At least seek one partner
							binary_key[ComparatorByPartnershipSought.INDEX_GENDER] = tar_gender;
							binary_key[tar_sought_comparator_partner_type_index] = 1;
							binary_key[ComparatorByPartnershipSought.INDEX_ID] = -1;

							candidateRangeByGender[cG][0] = ~Arrays.binarySearch(candidates, cg_start, cg_end,
									binary_key, comparators[tar_sought_candidate_partner_type_index]);

							binary_key[tar_sought_comparator_partner_type_index] = Integer.MAX_VALUE;
							binary_key[ComparatorByPartnershipSought.INDEX_ID] = Integer.MAX_VALUE;

							candidateRangeByGender[cG][1] = ~Arrays.binarySearch(candidates, cg_start, cg_end,
									binary_key, comparators[tar_sought_candidate_partner_type_index]);

							totalLength[cG] = candidateRangeByGender[cG][1] - candidateRangeByGender[cG][0];

							if (cG > 0) {
								totalLength[cG] += totalLength[cG - 1];
							}
						}

						int[][] target_candidates = new int[totalLength[totalLength.length - 1]][];

						if (target_candidates.length == 0) {
							break determine_src_candidate;
						}

						int offset = 0;
						for (int cG = 0; cG < candidateRangeByGender.length; cG++) {
							int copyLength = candidateRangeByGender[cG][1] - candidateRangeByGender[cG][0];
							System.arraycopy(candidates, candidateRangeByGender[cG][0], target_candidates, offset,
									copyLength);
							offset += copyLength;
						}

						Arrays.sort(target_candidates, new Comparator<int[]>() {
							@Override
							public int compare(int[] o1, int[] o2) {
								return Integer.compare(o1[tar_sought_comparator_partner_type_index],
										o2[tar_sought_comparator_partner_type_index]);
							}
						});

						int[] cumul_weight = new int[target_candidates.length];
						cumul_weight[0] = target_candidates[0][tar_sought_comparator_partner_type_index];
						for (int i = 1; i < target_candidates.length; i++) {
							cumul_weight[i] = cumul_weight[i - 1]
									+ target_candidates[i][tar_sought_comparator_partner_type_index];
						}

						// Choose partnership based on number sought
						numPartnerToSought = Math.min(numPartnerToSought, target_candidates.length);
						int[][] partnered_with = new int[numPartnerToSought][];
						int cumul_weight_lastIndex = cumul_weight.length - 1;
						int tar_candidate_index = -1;
						int num_partner_found = 0;

						// Choose partnership based on number sought
						for (int i = 0; i < partnered_with.length; i++) {
							if (tar_candidate_index != -1) {
								if (tar_candidate_index >= cumul_weight.length) {
									// Already removed?

								} else {
									int adjWeight = cumul_weight[tar_candidate_index];
									if (tar_candidate_index > 0) {
										adjWeight -= cumul_weight[tar_candidate_index - 1];
									}
									for (int adj = tar_candidate_index; adj < cumul_weight_lastIndex; adj++) {
										target_candidates[adj] = target_candidates[adj + 1];
										cumul_weight[adj] = cumul_weight[adj + 1] - adjWeight;
									}
								}

							}
							int prob_weight = getRNG().nextInt(cumul_weight[cumul_weight_lastIndex]);
							tar_candidate_index = Arrays.binarySearch(cumul_weight, prob_weight);
							if (tar_candidate_index < 0) {
								tar_candidate_index = ~tar_candidate_index;
							}
							partnered_with[i] = target_candidates[tar_candidate_index];
							num_partner_found++;

							// Re-pick if the same person was picked
							if (partnered_with[i][ComparatorByPartnershipSought.INDEX_ID] == src_candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_ID]) {
								num_partner_found--;
								for (int rI = tar_candidate_index; rI < target_candidates.length - 1; rI++) {
									target_candidates[rI] = target_candidates[rI + 1];
									cumul_weight[rI] = target_candidates[rI][tar_sought_comparator_partner_type_index];
									if (rI > 0) {
										cumul_weight[rI] += cumul_weight[rI - 1];
									}
								}
								target_candidates = Arrays.copyOf(target_candidates, target_candidates.length - 1);
								cumul_weight = Arrays.copyOf(cumul_weight, cumul_weight.length - 1);
								if (target_candidates.length == 0) {
									break;
								}
								i--;
							}
							cumul_weight_lastIndex--;
						}

						if (num_partner_found < partnered_with.length) {
							partnered_with = Arrays.copyOf(partnered_with, num_partner_found);
						}

						for (int i = 0; i < partnered_with.length; i++) {
							int[] tar_candidate_cmp_ent = partnered_with[i];
							Integer partner_form_time = getGlobalTime();

							if (schedule_range > 1) {
								partner_form_time += getRNG().nextInt(AbstractIndividualInterface.ONE_YEAR_INT);
							}

							// Determine if it can be a regular or casual partnership
							int partnership_type = tar_sought_candidate_partner_type_index;

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

							updateCandidateComparatorEntry(tar_candidate_cmp_ent);
							updateCandidateComparatorEntry(src_candidate_cmp_ent);

							updateScheduledPopDiff(addressed_demand_so_far, src_candidate_cmp_ent, categories_values);
							updateScheduledPopDiff(addressed_demand_so_far, tar_candidate_cmp_ent, categories_values);
						}

					}

				}

			}

			// Debug statement
			if (scheduling_debug) {
				System.out.printf("Schedule Partnership at Day %d:\n", getGlobalTime());
				for (int g = 0; g < LENGTH_GENDER; g++) {
					int pdI = numCat + g * numCat;
					System.out.printf(" %d : %s\n", g,
							Arrays.toString(Arrays.copyOfRange(addressed_demand_so_far, pdI, pdI + numCat)));
				}
				System.out.println();
			}
		}

		ArrayList<Integer[]> partnerships = schedule_partnership.remove(getGlobalTime());

		// The formation of partnership
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

				int relmapType;
				if (pair[1].getGenderType() == GENDER_FEMALE || pair[0].getGenderType() == GENDER_FEMALE) {
					relmapType = RELMAP_HETRO;
				} else {
					relmapType = RELMAP_MSM;
				}

				int partnershipType = edge[SCHEDULE_PARTNERSHIP_TYPE];

				if (partnershipType == CANDIDATE_ARRAY_SOUGHT_ANY) {
					RelationshipMap relMap = getRelMap()[relmapType];

					boolean hasExistingRegPartner = false;
					for (int p = 0; p < 2; p++) {
						hasExistingRegPartner |= relMap.containsVertex(pair[p].getId())
								&& relMap.degreeOf(pair[p].getId()) > 0;
					}

					if (hasExistingRegPartner) {
						partnershipType = CANDIDATE_ARRAY_SOUGHT_CAS;
					} else {
						partnershipType = CANDIDATE_ARRAY_SOUGHT_REG;
					}

				}

				if (partnershipType == CANDIDATE_ARRAY_SOUGHT_REG) { // Regular

					// Assume no concurrency
					for (int p = 0; p < pair.length; p++) {
						removeSingleRelationship(pair[p]);
					}

					int duration = regPartDuration[relmapType].sample();

					formRelationship(pair, getRelMap()[relmapType], duration, relmapType);

				} else { // Casual
					ContactMap[] cMaps = new ContactMap[] {
							((ContactMap[]) getFields()[FIELD_CONTACT_MAP])[CONTACT_MAP_ALL] };

					// Add casual partner
					for (int p = 0; p < pair.length; p++) {
						pair[p].addCasualPartner(pair[(p + 1) % 2]);
					}

					// Casual partnership of duration of 1 day
					checkContactMaps(new Integer[] { pair[0].getId(), pair[1].getId(), getGlobalTime(), 1 }, cMaps);
				}

			}
		}
	}

	private void fillCandidateList(int[][] candidates, int[] gender_end, boolean inclCurrent) {
		int pI = 0;

		float[] field_mean_target = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int[] numInGrp = (int[]) (getFields()[FIELD_POP_COMPOSITION]);
		int numCat = field_mean_target.length / (1 + LENGTH_GENDER);

		int[] numNotSoughting = new int[LENGTH_GENDER];
		int[][][] singleSoughtCandidate = new int[LENGTH_GENDER][][];
		int[] singleSoughtEnd = new int[LENGTH_GENDER];

		for (int g = 0; g < numNotSoughting.length; g++) {
			numNotSoughting[g] = (int) (field_mean_target[(numCat + g * numCat)] * numInGrp[g]); // Floor
			singleSoughtCandidate[g] = new int[numInGrp[g]][];
		}

		for (AbstractIndividualInterface absPerson : getPop()) {
			Person_Bridging_Pop person = (Person_Bridging_Pop) absPerson;
			int[] ent = new int[ComparatorByPartnershipSought.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];
			ent[ComparatorByPartnershipSought.INDEX_ID] = person.getId();
			ent[ComparatorByPartnershipSought.INDEX_GENDER] = person.getGenderType();

			ent[ComparatorByPartnershipSought.INDEX_MAX_REG] = (int) person
					.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS));
			ent[ComparatorByPartnershipSought.INDEX_MAX_CAS] = (int) person
					.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));

			if (inclCurrent) {
				ent[ComparatorByPartnershipSought.INDEX_CURRENT_REG] = getNumRegularPartnersCurrently(person);
				ent[ComparatorByPartnershipSought.INDEX_CURRENT_CAS] = person.getNumCasualInRecord();

			} else {
				ent[ComparatorByPartnershipSought.INDEX_CURRENT_REG] = 0;
				ent[ComparatorByPartnershipSought.INDEX_CURRENT_CAS] = 0;
			}

			ent[ComparatorByPartnershipSought.INDEX_CURRENT_ANY] = ent[ComparatorByPartnershipSought.INDEX_CURRENT_REG]
					+ ent[ComparatorByPartnershipSought.INDEX_CURRENT_CAS];

			ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG] = Math.max(0,
					ent[ComparatorByPartnershipSought.INDEX_MAX_REG]
							- ent[ComparatorByPartnershipSought.INDEX_CURRENT_ANY]);
			ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS] = Math.max(0,
					ent[ComparatorByPartnershipSought.INDEX_MAX_CAS]
							- ent[ComparatorByPartnershipSought.INDEX_CURRENT_ANY]);
			ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] = Math.max(
					ent[ComparatorByPartnershipSought.INDEX_MAX_REG], ent[ComparatorByPartnershipSought.INDEX_MAX_CAS])
					- ent[ComparatorByPartnershipSought.INDEX_CURRENT_ANY];

			ent[ComparatorByPartnershipSought.INDEX_SCHEDULE_LIMIT] = ent[ComparatorByPartnershipSought.INDEX_CURRENT_ANY]
					+ ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY];

			candidates[pI] = ent;
			gender_end[person.getGenderType()]++;

			if (ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] == 1) {
				singleSoughtCandidate[person.getGenderType()][singleSoughtEnd[person.getGenderType()]] = ent;
				singleSoughtEnd[person.getGenderType()]++;
			} else if (ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] == 0) {
				numNotSoughting[person.getGenderType()]--;
			}
			pI++;
		}
		for (int g = 1; g < gender_end.length; g++) {
			gender_end[g] += gender_end[g - 1];
		}

		// Randomly pick those who will non be sought this year round

		for (int g = 0; g < LENGTH_GENDER; g++) {
			for (int i = 0; i < singleSoughtEnd[g] && numNotSoughting[g] > 0; i++) {
				if (getRNG().nextInt(singleSoughtEnd[g] - i) < numNotSoughting[g]) {
					int[] ent = singleSoughtCandidate[g][i];
					ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG] = 0;
					ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS] = 0;
					ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] = 0;
					numNotSoughting[g]--;
				}

			}
		}

	}

	private void updateScheduledPopDiff(int[] scheduled_pop_diff_so_far, int[] candidate_cmp_ent,
			float[] categories_values) {
		int num_partner_currently = Math.max(candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_MAX_REG],
				candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_MAX_CAS])
				- candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY];

		int catAdj_prev = Arrays.binarySearch(categories_values, num_partner_currently - 1);

		if (catAdj_prev < 0) {
			catAdj_prev = ~catAdj_prev;
		}

		int catAdj = Arrays.binarySearch(categories_values, num_partner_currently);
		if (catAdj < 0) {
			catAdj = ~catAdj;
		}

		if (catAdj != catAdj_prev) {
			int adj_partner_dist_index = categories_values.length
					+ candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_GENDER] * categories_values.length + catAdj;
			scheduled_pop_diff_so_far[adj_partner_dist_index]++;

			if (catAdj_prev > 0) {
				// Subtract previous if it already schedule previously
				adj_partner_dist_index = categories_values.length
						+ candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_GENDER] * categories_values.length
						+ catAdj_prev;
				scheduled_pop_diff_so_far[adj_partner_dist_index]--;

			}
		}
	}

	private void updateCandidateComparatorEntry(int[] candidate_cmp_ent) {
		candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY]--;
		if (candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG] > 0) {
			candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_REG]--;
		}
		if (candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS] > 0) {
			candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_CAS]--;
		}
		if (candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_NUM_TO_SOUGHT_ANY] == 0) {
			// No need to sought anymore under current schedule
			candidate_cmp_ent[ComparatorByPartnershipSought.INDEX_SCHEDULE_LIMIT] = 0;
		}

	}

	private static final class ComparatorByPartnershipSought implements Comparator<int[]> {

		private static int INDEX_ID = 0;
		private static int INDEX_GENDER = INDEX_ID + 1;
		private static int INDEX_MAX_REG = INDEX_GENDER + 1;
		private static int INDEX_MAX_CAS = INDEX_MAX_REG + 1;
		private static int INDEX_CURRENT_REG = INDEX_MAX_CAS + 1;
		private static int INDEX_CURRENT_CAS = INDEX_CURRENT_REG + 1;
		private static int INDEX_CURRENT_ANY = INDEX_CURRENT_CAS + 1;
		private static int INDEX_NUM_TO_SOUGHT_REG = INDEX_CURRENT_ANY + 1;
		private static int INDEX_NUM_TO_SOUGHT_CAS = INDEX_NUM_TO_SOUGHT_REG + 1;
		private static int INDEX_NUM_TO_SOUGHT_ANY = INDEX_NUM_TO_SOUGHT_CAS + 1;
		private static int INDEX_SCHEDULE_LIMIT = INDEX_NUM_TO_SOUGHT_ANY + 1;
		private static int LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT = INDEX_SCHEDULE_LIMIT + 1;

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
