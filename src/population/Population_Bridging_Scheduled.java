package population;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.AbstractIntegerDistribution;

import infection.AbstractInfection;
import person.AbstractIndividualInterface;
import population.person.Person_Bridging_Pop;
import relationship.ContactMap;
import relationship.RelationshipMap;
import sim.Abstract_Runnable_ClusterModel_ContactMap_Generation;
import sim.Simulation_ClusterModelGeneration;

public class Population_Bridging_Scheduled extends Population_Bridging {
	/**
	 * 
	 */
	private static final long serialVersionUID = 7076023161972931765L;

	protected HashMap<Integer, ArrayList<Integer[]>> schedule_partnership;

	protected transient int lastPartnershipScheduling = -1;
	protected transient boolean fitHighActFirst = false;

	private final boolean schedule_debug = !true;
	private boolean space_save = false;

	protected long export_period_form_partnership_progress = 5 * 60 * 1000l;
	public static final String FORMAT_FORM_PARTNERSHIP_PROGRESS_PREFIX = "FormPartnership_Progess_%d";
	public static final String FORMAT_FORM_PARTNERSHIP_PROGRESS = FORMAT_FORM_PARTNERSHIP_PROGRESS_PREFIX + "_%d.obj";

	public static final int SCHEDULE_PARTNERSHIP_P1 = 0;
	public static final int SCHEDULE_PARTNERSHIP_P2 = SCHEDULE_PARTNERSHIP_P1 + 1;
	public static final int SCHEDULE_PARTNERSHIP_TYPE = SCHEDULE_PARTNERSHIP_P2 + 1;
	public static final int LENGTH_SCHEDULE_PARTNERSHIP = SCHEDULE_PARTNERSHIP_TYPE + 1;

	protected static final int CANDIDATE_ARRAY_SCHEDULE_LIMIT = 0;
	protected static final int CANDIDATE_ARRAY_SOUGHT_ANY = CANDIDATE_ARRAY_SCHEDULE_LIMIT + 1;
	protected static final int CANDIDATE_ARRAY_SOUGHT_REG = CANDIDATE_ARRAY_SOUGHT_ANY + 1;
	protected static final int CANDIDATE_ARRAY_SOUGHT_CAS = CANDIDATE_ARRAY_SOUGHT_REG + 1;

	protected final Comparator<Integer[]> schedule_partnership_comparator = new Comparator<Integer[]>() {
		@Override
		public int compare(Integer[] o1, Integer[] o2) {
			int res = 0;
			int pt = 0;
			while (res == 0 && pt < o1.length) {
				res = Integer.compare(o1[pt], o2[pt]);
				pt++;
			}
			return res;
		}
	};

	public Population_Bridging_Scheduled(long seed) {
		super(seed);
		schedule_partnership = new HashMap<>();
		Object[] newFields = Arrays.copyOf(super.getFields(), super.getFields().length + 1);
		newFields[newFields.length - 1] = schedule_partnership;
		super.setFields(newFields);
	}

	public void setSpace_save(boolean space_save) {
		this.space_save = space_save;
	}

	public boolean isSpace_save() {
		return this.space_save;
	}

	@Override
	public void initialise() {
		// Check to see if fitting of low activity group first
		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);

		for (int i = 0; i < field_mean_number_partner.length; i++) {
			fitHighActFirst |= field_mean_number_partner[i] < 0;
			field_mean_number_partner[i] = Math.abs(field_mean_number_partner[i]);
		}

		if (fitHighActFirst && printStatus != null) {
			for (PrintStream out : printStatus) {
				out.println("Parntership will be allocated to high activity groups first.");
			}
		}
		super.initialise();
	}

	@Override
	protected void initialiseTransientFields() {

		// Check to see if fitting of low activity group first
		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);

		for (int i = 0; i < field_mean_number_partner.length; i++) {
			fitHighActFirst &= field_mean_number_partner[i] < 0;
			field_mean_number_partner[i] = Math.abs(field_mean_number_partner[i]);
		}

		if (fitHighActFirst && printStatus != null) {
			for (PrintStream out : printStatus) {
				out.println("Parntership will be allocated to high activity groups first.");
			}
		}

		super.initialiseTransientFields();
	}

	@SuppressWarnings("unchecked")
	public static Population_Bridging decodeFromStream(java.io.ObjectInputStream inStr)
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
		pop.lastPartnershipScheduling = (globalTime / AbstractIndividualInterface.ONE_YEAR_INT)
				* AbstractIndividualInterface.ONE_YEAR_INT;

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

		float[] field_mean_number_partner = (float[]) (getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS]);
		int[] population_num_partner_in_last_12_months = new int[field_mean_number_partner.length];
		int numCat = field_mean_number_partner.length / (LENGTH_GENDER + 1);

		// Check if it was in the middle of form partnership progress
		final File progressFile = new File(getBaseDir(),
				String.format(FORMAT_FORM_PARTNERSHIP_PROGRESS, this.getSeed(), getGlobalTime()));

		if (progressFile.exists()) {
			lastPartnershipScheduling = -1;
			formPartnerships(population_num_partner_in_last_12_months);
		}

		incrementTime(deltaT);

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
			// Export contact map
			ContactMap cMap = ((ContactMap[]) getFields()[Population_Bridging.FIELD_CONTACT_MAP])[0];

			if (cMap != null) {
				File allContactFile = new File(baseDir, String.format(
						Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP, getSeed(), getGlobalTime()));

				BufferedWriter fileWriAll;

				try {
					fileWriAll = new BufferedWriter(new FileWriter(allContactFile));
					fileWriAll.append(cMap.toFullString());
					fileWriAll.close();
				} catch (IOException e) {
					e.printStackTrace(System.err);
					System.out.println(cMap.toFullString());
				}

				if (space_save) {

					String str_pattern_contactMap = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP
							.replaceFirst("%d", Long.toString(getSeed()));
					str_pattern_contactMap = str_pattern_contactMap.replaceFirst("%d", "(-{0,1}\\\\d+)");
					Pattern pattern_contactMap = Pattern.compile(str_pattern_contactMap);
					File[] oldContactMap = baseDir.listFiles(new FileFilter() {
						@Override
						public boolean accept(File pathname) {
							boolean res;
							Matcher m = pattern_contactMap.matcher(pathname.getName());
							res = m.matches();
							if (res) {
								res = !pathname.getName().equals(allContactFile.getName());
							}
							return res;
						}
					});

					for (File oldMap : oldContactMap) {
						try {
							FileUtils.delete(oldMap);
						} catch (IOException ex) {
							ex.printStackTrace(System.err);
							oldMap.deleteOnExit();
						}
					}

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
					pri.println(String.format("# partners in last 12 months -  Day %d\n", getGlobalTime()));
					for (int g = 0; g < LENGTH_GENDER; g++) {
						pri.printf(" %d: %s\n", g,
								Arrays.toString(Arrays.copyOfRange(population_num_partner_in_last_12_months,
										numCat + g * numCat, 2 * numCat + g * numCat)));
					}
					for (PrintStream out : printStatus) {
						out.println(wri.toString());
					}
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

	@SuppressWarnings("unchecked")
	@Override
	public void formPartnerships(int[] population_num_partner_in_last_12_months) {

		boolean reqPartnerScheduling = (lastPartnershipScheduling == -1) ? true
				: (getGlobalTime() == lastPartnershipScheduling + AbstractIndividualInterface.ONE_YEAR_INT);

		int schedule_range = AbstractIndividualInterface.ONE_YEAR_INT;

		if (reqPartnerScheduling) {

			long last_export_at = System.currentTimeMillis();

			lastPartnershipScheduling = getGlobalTime();

			int numCat = population_num_partner_in_last_12_months.length / (1 + LENGTH_GENDER);
			float[] categories_values = Arrays.copyOf((float[]) getFields()[FIELD_MEAN_NUM_PARTNER_IN_12_MONTHS],
					numCat);

			// Candidates arrays

			int[][] candidates_all = new int[getPop().length][];

			Comparator_Candidate_Entry[] candidates_array_comparators = new Comparator_Candidate_Entry[] {
					new Comparator_Candidate_Entry(Comparator_Candidate_Entry.INDEX_SCHEDULE_LIMIT),
					generateTargetCandidatesComparator(Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY),
					generateTargetCandidatesComparator(Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_REG),
					generateTargetCandidatesComparator(Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_CAS) };

			int[] candidates_array_num_sought_indices = new int[] { CANDIDATE_ARRAY_SOUGHT_ANY,
					CANDIDATE_ARRAY_SOUGHT_REG, CANDIDATE_ARRAY_SOUGHT_CAS };

			int[] binary_key = new int[Comparator_Candidate_Entry.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];

			// Imported from file

			int[][][] candidates_array_by_partnership_type = new int[candidates_array_comparators.length][][];
			int[] gender_end = new int[LENGTH_GENDER];

			// Population demand and demand addressed
			int[] pop_demand_num_partner_12_months = cal_pop_diff_num_partner(
					new int[population_num_partner_in_last_12_months.length]);
			int[] addressed_demand_so_far = new int[population_num_partner_in_last_12_months.length];

			// Categories order for source
			int[] src_cat_order = new int[numCat - 1];

			int[] completed_src_pdIndex = new int[src_cat_order.length * LENGTH_GENDER];
			Arrays.fill(completed_src_pdIndex, -1);
			int next_completed_src_pdIndex_pt = 0;

			int progressing_src_pdIndex = -1;
			ArrayList<int[]> progressing_src_candidate_list = null;
			int completed_src_candidate_index = 0;

			// End of Imported from file

			final File progressFile = new File(getBaseDir(),
					String.format(FORMAT_FORM_PARTNERSHIP_PROGRESS, this.getSeed(), getGlobalTime()));
			boolean importSuccess = false;

			if (progressFile.isFile()) {
				ObjectInputStream objIn;
				try {
					objIn = new ObjectInputStream(new FileInputStream(progressFile));
					schedule_partnership = (HashMap<Integer, ArrayList<Integer[]>>) objIn.readObject();
					candidates_array_by_partnership_type = (int[][][]) objIn.readObject();
					gender_end = (int[]) objIn.readObject();
					addressed_demand_so_far = (int[]) objIn.readObject();
					completed_src_pdIndex = (int[]) objIn.readObject();
					next_completed_src_pdIndex_pt = objIn.readInt();
					progressing_src_pdIndex = objIn.readInt();
					progressing_src_candidate_list = (ArrayList<int[]>) objIn.readObject();
					completed_src_candidate_index = objIn.readInt();
					objIn.close();

					candidates_all = candidates_array_by_partnership_type[0];
					Arrays.sort(completed_src_pdIndex, 0, next_completed_src_pdIndex_pt);

					importSuccess = true;

				} catch (ClassNotFoundException | IOException e) {
					File tempFile = new File(baseDir, String.format("%s_temp", progressFile.getName()));

					if (tempFile.exists()) {

						try {
							objIn = new ObjectInputStream(new FileInputStream(tempFile));
							schedule_partnership = (HashMap<Integer, ArrayList<Integer[]>>) objIn.readObject();
							candidates_array_by_partnership_type = (int[][][]) objIn.readObject();
							gender_end = (int[]) objIn.readObject();
							addressed_demand_so_far = (int[]) objIn.readObject();
							completed_src_pdIndex = (int[]) objIn.readObject();
							next_completed_src_pdIndex_pt = objIn.readInt();
							progressing_src_pdIndex = objIn.readInt();
							progressing_src_candidate_list = (ArrayList<int[]>) objIn.readObject();
							completed_src_candidate_index = objIn.readInt();
							objIn.close();

							candidates_all = candidates_array_by_partnership_type[0];
							Arrays.sort(completed_src_pdIndex, 0, next_completed_src_pdIndex_pt);

							importSuccess = true;
						} catch (ClassNotFoundException | IOException e1) {
							e1.printStackTrace(System.err);
							importSuccess = false;

						}
					} else {
						e.printStackTrace(System.err);
						importSuccess = false;
					}
				}
			}

			if (!importSuccess) {
				candidates_array_by_partnership_type[0] = candidates_all;
				fillCandidateList(candidates_all, gender_end, schedule_range == 0);
				Arrays.sort(candidates_all, candidates_array_comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);
			}

			Arrays.fill(src_cat_order, -1);
			int src_cat_order_pt = 0;

			if (fitHighActFirst) {
				// Order by activity (highest to lowest)
				for (int src_cat_index = numCat - 1; src_cat_index > 0; src_cat_index--) {
					src_cat_order[src_cat_order_pt] = src_cat_index;
					src_cat_order_pt++;
				}
			} else {
				for (int src_cat_index = 1; src_cat_index < numCat; src_cat_index++) {
					src_cat_order[src_cat_order_pt] = src_cat_index;
					src_cat_order_pt++;
				}
			}

			for (int src_cat_index : src_cat_order) {
				for (int src_gender = 0; src_gender < LENGTH_GENDER; src_gender++) {
					int src_gender_index_start = src_gender > 0 ? gender_end[src_gender - 1] : 0;
					int src_gender_index_end = gender_end[src_gender];
					int src_pdIndex = numCat + src_gender * numCat + src_cat_index;

					if (Arrays.binarySearch(completed_src_pdIndex, 0, next_completed_src_pdIndex_pt, src_pdIndex) < 0) {

						if (pop_demand_num_partner_12_months[src_pdIndex] > addressed_demand_so_far[src_pdIndex]) {

							if (progressing_src_pdIndex != src_pdIndex || progressing_src_candidate_list == null) {

								int[][] src_candidate_cmp_ent_array;
								src_candidate_cmp_ent_array = new int[pop_demand_num_partner_12_months[src_pdIndex]
										- addressed_demand_so_far[src_pdIndex]][];

								int numPartToSought_max = (int) categories_values[src_cat_index];
								int numPartToSought_min = ((src_cat_index == 0) ? 0
										: (int) categories_values[src_cat_index - 1]) + 1;

								Arrays.fill(binary_key, -1);
								binary_key[Comparator_Candidate_Entry.INDEX_GENDER] = src_gender;

								int[] numPartToSoughtAdj = new int[] { numPartToSought_min, numPartToSought_max, };
								int[] src_index_range = new int[2];

								Arrays.fill(binary_key, 0);
								binary_key[Comparator_Candidate_Entry.INDEX_GENDER] = src_gender;

								binary_key[Comparator_Candidate_Entry.INDEX_ID] = -1;
								binary_key[Comparator_Candidate_Entry.INDEX_SCHEDULE_LIMIT] = numPartToSoughtAdj[0];
								src_index_range[0] = ~Arrays.binarySearch(candidates_all, src_gender_index_start,
										src_gender_index_end, binary_key,
										candidates_array_comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);

								binary_key[Comparator_Candidate_Entry.INDEX_ID] = Integer.MAX_VALUE;
								binary_key[Comparator_Candidate_Entry.INDEX_SCHEDULE_LIMIT] = numPartToSoughtAdj[1];
								src_index_range[1] = ~Arrays.binarySearch(candidates_all, src_gender_index_start,
										src_gender_index_end, binary_key,
										candidates_array_comparators[CANDIDATE_ARRAY_SCHEDULE_LIMIT]);

								int next_index_src_candidate_cmp_ent_array = 0;

								for (int i = src_index_range[0]; i < src_index_range[1]
										&& next_index_src_candidate_cmp_ent_array < src_candidate_cmp_ent_array.length; i++) {

									if (getRNG().nextInt(src_index_range[1] - i) < src_candidate_cmp_ent_array.length
											- next_index_src_candidate_cmp_ent_array) {
										src_candidate_cmp_ent_array[next_index_src_candidate_cmp_ent_array] = candidates_all[i];
										next_index_src_candidate_cmp_ent_array++;
									}
								}

								if (next_index_src_candidate_cmp_ent_array < src_candidate_cmp_ent_array.length) {
									src_candidate_cmp_ent_array = Arrays.copyOf(src_candidate_cmp_ent_array,
											next_index_src_candidate_cmp_ent_array);
								}

								progressing_src_candidate_list = new ArrayList<>(List.of(src_candidate_cmp_ent_array));

								if (progressing_src_candidate_list.size() > 0) {
									int[] tar_possible_gender;
									switch (src_gender) {
									case GENDER_FEMALE:
										tar_possible_gender = prob_no_bridge < 0
												? new int[] { GENDER_HETRO_MALE, GENDER_MSMW }
												: new int[] { GENDER_HETRO_MALE };
										break;
									case GENDER_HETRO_MALE:
										tar_possible_gender = new int[] { GENDER_FEMALE };
										break;
									case GENDER_MSMO:
										tar_possible_gender = new int[] { GENDER_MSMO, GENDER_MSMW };
										break;
									default:
										tar_possible_gender = prob_no_bridge < 0
												? new int[] { GENDER_FEMALE, GENDER_MSMO, GENDER_MSMW }
												: new int[] { GENDER_MSMO, GENDER_MSMW };
									}

									Collections.shuffle(progressing_src_candidate_list,
											new Random(getRNG().nextLong()));

									for (int candidates_array_num_sought_index : candidates_array_num_sought_indices) {
										int[][] target_candidates = candidates_array_by_partnership_type[candidates_array_num_sought_index];
										target_candidates = generateTargetCandidateArray(candidates_all,
												candidates_array_comparators, gender_end, tar_possible_gender,
												candidates_array_num_sought_index);
										Arrays.sort(target_candidates,
												candidates_array_comparators[candidates_array_num_sought_index]);
										candidates_array_by_partnership_type[candidates_array_num_sought_index] = target_candidates;

									}

								}
							}

							int src_candidate_counter = 0;

							for (int[] src_candidate_cmp_ent : progressing_src_candidate_list) {

								if (progressing_src_pdIndex != src_pdIndex
										|| src_candidate_counter > completed_src_candidate_index) {

									boolean onlySoughtReg = src_candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_MAX_CAS] == 0;
									boolean onlySoughtCas = src_candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_MAX_REG] == 0;
									int numPartnerToSought = src_candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY];

									final int tar_sought_candidate_partner_type_index;
									if (onlySoughtReg) {
										tar_sought_candidate_partner_type_index = CANDIDATE_ARRAY_SOUGHT_REG;

									} else if (onlySoughtCas) {
										tar_sought_candidate_partner_type_index = CANDIDATE_ARRAY_SOUGHT_CAS;

									} else {
										tar_sought_candidate_partner_type_index = CANDIDATE_ARRAY_SOUGHT_ANY;
									}

									int[][] target_candidates = candidates_array_by_partnership_type[tar_sought_candidate_partner_type_index];

									Comparator_Candidate_Entry target_candidate_comparator = candidates_array_comparators[tar_sought_candidate_partner_type_index];

									int tar_sought_comparator_partner_type_index = target_candidate_comparator
											.getCmpMethod();
									if (tar_sought_comparator_partner_type_index < 0) {
										tar_sought_comparator_partner_type_index = ~tar_sought_comparator_partner_type_index;
									}

									ArrayList<int[]> src_specific_candidate_list = new ArrayList<>(
											List.of(target_candidates));

									src_specific_candidate_list = new ArrayList<>(
											cleanUpCandidateList(src_specific_candidate_list,
													target_candidate_comparator, src_candidate_cmp_ent));

									numPartnerToSought = Math.min(numPartnerToSought, target_candidates.length);
									int[][] partnered_with = new int[numPartnerToSought][];

									int tar_candidate_index = -1;
									int num_partner_found = 0;

									// Assortativity
									boolean assort_mix = false;
									float[] part_type = ((float[][]) getFields()[FIELD_PARTNER_TYPE_PROB])[src_candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_GENDER]];

									int assort_index = -1; // ASSORTATIVITY_INC_INDEX.ASSORTATIVITY_PROB

									if (PARTNER_TYPE_ASSORTATIVITY < part_type.length
											&& src_specific_candidate_list.size() > 0) {
										assort_index = (int) part_type[PARTNER_TYPE_ASSORTATIVITY];
										assort_mix = getRNG()
												.nextFloat() < (part_type[PARTNER_TYPE_ASSORTATIVITY] - assort_index);

									}

									for (int partnerIndex = 0; partnerIndex < partnered_with.length
											&& src_specific_candidate_list.size() > 0; partnerIndex++) {
										float[] cumul_weight = new float[src_specific_candidate_list.size()];
										int cumul_weight_pt = 0;
										int maxAssortWeightDiff = -1;

										if (assort_mix) {
											for (int[] candidate : src_specific_candidate_list) {
												int total_weight = 0;
												for (int i = 0; i < src_candidate_cmp_ent.length; i++) {
													if ((assort_index & 1 << i) != 0) {
														total_weight += Math
																.abs(src_candidate_cmp_ent[i] - candidate[i]);
													}
												}
												maxAssortWeightDiff = Math.max(maxAssortWeightDiff, total_weight + 1);
											}

										}

										for (int[] target_candidate_cmp_ent : src_specific_candidate_list) {

											// Choose partnership based on number sought (by default)
											// but with additional adjustment to assortativity
											int partner_weight = target_candidate_cmp_ent[tar_sought_comparator_partner_type_index];

											if (partner_weight != 0 && assort_mix
													&& src_specific_candidate_list.size() > 1) {
												int weight_diff = 0;
												for (int i = 0; i < src_candidate_cmp_ent.length; i++) {
													if ((assort_index & 1 << i) != 0) {
														weight_diff += Math.abs(
																src_candidate_cmp_ent[i] - target_candidate_cmp_ent[i]);
													}
												}
												partner_weight = (maxAssortWeightDiff - weight_diff);
											}

											if (cumul_weight_pt == 0) {
												cumul_weight[cumul_weight_pt] = partner_weight;
											} else {
												cumul_weight[cumul_weight_pt] = cumul_weight[cumul_weight_pt - 1]
														+ partner_weight;
											}
											cumul_weight_pt++;
										}

										if (cumul_weight[cumul_weight.length - 1] == 0) {
											tar_candidate_index = getRNG().nextInt(cumul_weight.length);
										} else {
											float prob_weight = getRNG().nextFloat()
													* (cumul_weight[cumul_weight.length - 1]);
											tar_candidate_index = Arrays.binarySearch(cumul_weight, prob_weight);
											if (tar_candidate_index < 0) {
												tar_candidate_index = ~tar_candidate_index;
											} else {
												if (tar_candidate_index + 1 < cumul_weight.length) {
													tar_candidate_index++;
												}
											}
										}

										partnered_with[partnerIndex] = src_specific_candidate_list
												.remove(tar_candidate_index);

										// Update target array
										int[] tar_candidate_cmp_ent = partnered_with[partnerIndex];
										int[][] org_target_range = generateCandidateArraySearchRange(
												tar_candidate_cmp_ent, candidates_array_by_partnership_type,
												candidates_array_comparators, candidates_array_num_sought_indices,
												gender_end);
										updateCandidateComparatorEntry(tar_candidate_cmp_ent);
										updateCandidateArrayAll(tar_candidate_cmp_ent,
												candidates_array_by_partnership_type, candidates_array_comparators,
												candidates_array_num_sought_indices, org_target_range);

										// Update source array
										int[][] org_src_range = generateCandidateArraySearchRange(src_candidate_cmp_ent,
												candidates_array_by_partnership_type, candidates_array_comparators,
												candidates_array_num_sought_indices, gender_end);
										updateCandidateComparatorEntry(src_candidate_cmp_ent);
										updateCandidateArrayAll(src_candidate_cmp_ent,
												candidates_array_by_partnership_type, candidates_array_comparators,
												candidates_array_num_sought_indices, org_src_range);

										updateScheduledPopDiff(addressed_demand_so_far, src_candidate_cmp_ent,
												categories_values);
										updateScheduledPopDiff(addressed_demand_so_far, tar_candidate_cmp_ent,
												categories_values);

										num_partner_found++;

										if (partnerIndex + 1 < partnered_with.length) {

											src_specific_candidate_list = new ArrayList<>(List.of(
													candidates_array_by_partnership_type[tar_sought_candidate_partner_type_index]));

											src_specific_candidate_list = new ArrayList<>(
													cleanUpCandidateList(src_specific_candidate_list,
															target_candidate_comparator, src_candidate_cmp_ent));

										}

									}

									for (int partnerIndex = 0; partnerIndex < num_partner_found; partnerIndex++) {
										int[] tar_candidate_cmp_ent = partnered_with[partnerIndex];
										Integer partner_form_time = getGlobalTime();

										if (schedule_range > 1) {
											partner_form_time += getRNG()
													.nextInt(AbstractIndividualInterface.ONE_YEAR_INT);
										}

										// Determine if it can be a regular or casual partnership
										int partnership_type = tar_sought_candidate_partner_type_index;

										Integer[] schedule_partnership_ent = new Integer[LENGTH_SCHEDULE_PARTNERSHIP];

										schedule_partnership_ent[SCHEDULE_PARTNERSHIP_P1] = src_candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_ID];
										schedule_partnership_ent[SCHEDULE_PARTNERSHIP_P2] = tar_candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_ID];
										schedule_partnership_ent[SCHEDULE_PARTNERSHIP_TYPE] = partnership_type;

										ArrayList<Integer[]> partnerships = schedule_partnership.get(partner_form_time);
										if (partnerships == null) {
											partnerships = new ArrayList<>();
											schedule_partnership.put(partner_form_time, partnerships);
										}

										int pt = Collections.binarySearch(partnerships, schedule_partnership_ent,
												schedule_partnership_comparator);
										if (pt < 0) {
											partnerships.add(~pt, schedule_partnership_ent);
										}

									}

								}
								src_candidate_counter++;

								if (export_period_form_partnership_progress > 0 && System.currentTimeMillis()
										- last_export_at > export_period_form_partnership_progress) {

									if (!space_save) {
										exportPartnerFormationProgress(progressFile,
												candidates_array_by_partnership_type, gender_end,
												addressed_demand_so_far, completed_src_pdIndex,
												next_completed_src_pdIndex_pt, src_pdIndex,
												progressing_src_candidate_list, src_candidate_counter);
									}

									last_export_at = System.currentTimeMillis();

									// Debug statement
									if (schedule_debug) {
										System.out.printf("Time = %tF %<tT\n", new Date());
										System.out.printf("Completed src_pd_index = %s\n", Arrays.toString(
												Arrays.copyOf(completed_src_pdIndex, next_completed_src_pdIndex_pt)));
										System.out.printf("Current src_pd_index = %d\n", src_pdIndex);
										System.out.printf("Completed src candidate = %d out of %d\n",
												src_candidate_counter, progressing_src_candidate_list.size());

										System.out.printf("Schedule Partnership at Day %d:\n", getGlobalTime());
										for (int g = 0; g < LENGTH_GENDER; g++) {
											int pdI = numCat + g * numCat;
											System.out.printf(" %d : %s\n", g, Arrays.toString(
													Arrays.copyOfRange(addressed_demand_so_far, pdI, pdI + numCat)));
										}
										System.out.println();
									}
								}
							}

						}

						completed_src_pdIndex[next_completed_src_pdIndex_pt] = src_pdIndex;
						next_completed_src_pdIndex_pt++;
						progressing_src_candidate_list = null;

					}
				}
			}

			// Final export to complete the process

			if (!space_save) {
				exportPartnerFormationProgress(progressFile, candidates_array_by_partnership_type, gender_end,
						addressed_demand_so_far, completed_src_pdIndex, next_completed_src_pdIndex_pt, -1, null, 0);
			}
			last_export_at = System.currentTimeMillis();

			final String oldProgressFile_prefix = String.format(FORMAT_FORM_PARTNERSHIP_PROGRESS_PREFIX, getSeed());
			final String oldExportPop_prefix = String
					.format(Abstract_Runnable_ClusterModel_ContactMap_Generation.EXPORT_POP_FILENAME_PRRFIX, getSeed());

			// Delete old export file
			File[] oldExportFile = baseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {

					boolean isOldFile = pathname.getName().startsWith(oldProgressFile_prefix);

					isOldFile |= !pathname.getName()
							.equals(String.format(
									Abstract_Runnable_ClusterModel_ContactMap_Generation.EXPORT_POP_FILENAME, getSeed(),
									getGlobalTime()))
							&& pathname.getName().startsWith(oldExportPop_prefix);

					return isOldFile;
				}
			});

			for (File toBeDeleted : oldExportFile) {
				try {
					FileUtils.delete(toBeDeleted);
				} catch (IOException ex) {
					ex.printStackTrace(System.err);
					toBeDeleted.deleteOnExit();
				}
			}

			// Debug statement
			if (schedule_debug) {
				System.out.printf("Time = %tF %<tT\n", new Date());
				System.out.printf("Completed src_pd_index = %s\n",
						Arrays.toString(Arrays.copyOf(completed_src_pdIndex, next_completed_src_pdIndex_pt)));
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

	private List<int[]> cleanUpCandidateList(List<int[]> candidate_list,
			Comparator_Candidate_Entry candidate_comparator, int[] remove_candidate_cmp_ent) {

		int[] binary_key = new int[Comparator_Candidate_Entry.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];
		Arrays.fill(binary_key, 0);

		int sought_comparator_partner_type_index = candidate_comparator.getCmpMethod();
		if (sought_comparator_partner_type_index < 0) {
			sought_comparator_partner_type_index = ~sought_comparator_partner_type_index;
		}

		/// Remove all that do not sought partners
		binary_key[Comparator_Candidate_Entry.INDEX_ID] = Integer.MAX_VALUE;
		binary_key[sought_comparator_partner_type_index] = 0;
		int excl_index = ~Collections.binarySearch(candidate_list, binary_key, candidate_comparator);
		candidate_list = candidate_list.subList(excl_index, candidate_list.size());

		if (remove_candidate_cmp_ent != null) {
			binary_key[Comparator_Candidate_Entry.INDEX_ID] = remove_candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_ID];
			binary_key[sought_comparator_partner_type_index] = remove_candidate_cmp_ent[sought_comparator_partner_type_index];

			// Remove self (if found)
			excl_index = Collections.binarySearch(candidate_list, binary_key, candidate_comparator);

			if (excl_index >= 0) {
				candidate_list.remove(excl_index);
			}
		}
		return candidate_list;
	}

	private int[][] generateCandidateArraySearchRange(int[] candidate_cmp_ent,
			int[][][] candidates_array_by_partnership_type, Comparator_Candidate_Entry[] candidates_array_comparators,
			int[] candidates_array_num_sought_indices, int[] gender_end) {
		int[][] candidate_range = new int[candidates_array_by_partnership_type.length][2];
		int candidate_gender = candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_GENDER];
		int candidate_gender_index_start = candidate_gender > 0 ? gender_end[candidate_gender - 1] : 0;
		int candidate_gender_index_end = gender_end[candidate_gender];

		candidate_range[0][0] = candidate_gender_index_start;
		candidate_range[0][1] = searchCandidateIndices(candidate_cmp_ent, candidates_array_by_partnership_type[0],
				candidates_array_comparators[0],
				new int[] { candidate_gender_index_start, candidate_gender_index_end });

		for (int candidates_array_num_sought_index : candidates_array_num_sought_indices) {
			candidate_range[candidates_array_num_sought_index][0] = 0;
			candidate_range[candidates_array_num_sought_index][1] = searchCandidateIndices(candidate_cmp_ent,
					candidates_array_by_partnership_type[candidates_array_num_sought_index],
					candidates_array_comparators[candidates_array_num_sought_index], null);
		}
		return candidate_range;
	}

	private void updateCandidateArrayAll(int[] cmp_entry, int[][][] candidates_array_by_partnership_type,
			Comparator_Candidate_Entry[] candidates_array_comparators, int[] candidates_array_num_sought_indices,
			int[][] ranges) {

		if (cmp_entry[Comparator_Candidate_Entry.INDEX_SCHEDULE_LIMIT] == 0) {
			updateCandidateArray(cmp_entry, candidates_array_by_partnership_type[0], candidates_array_comparators[0],
					ranges[0]);
		}

		for (int candidates_array_num_sought_index : candidates_array_num_sought_indices) {
			if (ranges[candidates_array_num_sought_index][1] > 0) {
				updateCandidateArray(cmp_entry, candidates_array_by_partnership_type[candidates_array_num_sought_index],
						candidates_array_comparators[candidates_array_num_sought_index],
						ranges[candidates_array_num_sought_index]);
			}
		}
	}

	private static void updateCandidateArray(int[] candidate_cmp_ent, int[][] candidates_array,
			Comparator_Candidate_Entry candidates_cmp, int[] search_range) {
		if (search_range != null && search_range[1] > 0) {
			int new_tar_index = searchCandidateIndices(candidate_cmp_ent, candidates_array, candidates_cmp,
					search_range);

			int replace_util_index = search_range[1];
			System.arraycopy(candidates_array, ~new_tar_index, candidates_array, ~new_tar_index + 1,
					replace_util_index - ~new_tar_index);
			candidates_array[~new_tar_index] = candidate_cmp_ent;

		}
	}

	private void exportPartnerFormationProgress(final File progressFile,
			final int[][][] f_candidates_array_by_partnership_type, final int[] f_gender_end,
			final int[] f_addressed_demand_so_far, final int[] f_completed_src_pdIndex,
			final int f_next_completed_src_pdIndex_pt, final int f_progressing_src_pdIndex,
			final ArrayList<int[]> f_src_candidate_list, final int f_completed_src_candidate_index) {

		Runnable exportThread = new Runnable() {

			@Override
			public void run() {
				try {
					File tempFile = null;
					File exportPopFile = new File(baseDir,
							String.format(Abstract_Runnable_ClusterModel_ContactMap_Generation.EXPORT_POP_FILENAME,
									getSeed(), getGlobalTime()));

					if (exportPopFile.isFile()) {
						tempFile = new File(baseDir, String.format("%s_temp", exportPopFile.getName()));
						Files.move(exportPopFile.toPath(), tempFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

					}

					ObjectOutputStream outStream = new ObjectOutputStream(
							new BufferedOutputStream(new FileOutputStream(exportPopFile)));
					encodePopToStream(outStream);
					outStream.close();

					if (tempFile != null) {
						FileUtils.delete(tempFile);
					}

					File newProgresFile = new File(baseDir, progressFile.getName());

					if (progressFile.isFile()) {
						tempFile = new File(baseDir, String.format("%s_temp", progressFile.getName()));
						Files.move(progressFile.toPath(), tempFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
					}

					// Only write progress if there are still src_pdIndex to go through
					if (f_next_completed_src_pdIndex_pt < f_completed_src_pdIndex.length) {
						ObjectOutputStream objOut = new ObjectOutputStream(new FileOutputStream(newProgresFile));
						objOut.writeObject(schedule_partnership);
						objOut.writeObject(f_candidates_array_by_partnership_type);
						objOut.writeObject(f_gender_end);
						objOut.writeObject(f_addressed_demand_so_far);
						objOut.writeObject(f_completed_src_pdIndex);
						objOut.writeInt(f_next_completed_src_pdIndex_pt);
						objOut.writeInt(f_progressing_src_pdIndex);
						objOut.writeObject(f_src_candidate_list);
						objOut.writeInt(f_completed_src_candidate_index);
						objOut.close();
					}

					if (tempFile != null) {
						FileUtils.delete(tempFile);
					}
				} catch (IOException e) {
					e.printStackTrace(System.err);
				}

			}
		};

		ExecutorService exec_export = Executors.newSingleThreadExecutor();
		exec_export.submit(exportThread);
		exec_export.shutdown();

	}

	private static int searchCandidateIndices(int[] candidate_ent, int[][] candidates_array,
			Comparator_Candidate_Entry candidates_cmp, int[] search_range) {
		int index;
		int search_start = 0;
		int search_end = candidates_array.length;

		if (search_range != null) {
			search_start = search_range[0];
			search_end = search_range[1];
		}

		index = Arrays.binarySearch(candidates_array, search_start, search_end, candidate_ent, candidates_cmp);
		return index;
	}

	private int[][] generateTargetCandidateArray(int[][] candidates_src, Comparator_Candidate_Entry[] comparators,
			int[] gender_end, int[] tar_possible_gender, final int tar_sought_candidate_partner_type_index) {

		int[][] candidateRangeByGender = new int[tar_possible_gender.length][2];
		int[] totalCandidateLength = new int[tar_possible_gender.length];

		int[] binary_key = new int[Comparator_Candidate_Entry.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];

		int[][] candidates = Arrays.copyOf(candidates_src, candidates_src.length);

		int sought_comparator_partner_type_index = comparators[tar_sought_candidate_partner_type_index].getCmpMethod();

		if (sought_comparator_partner_type_index < 0) {
			sought_comparator_partner_type_index = ~sought_comparator_partner_type_index;
		}

		for (int cG = 0; cG < tar_possible_gender.length; cG++) {
			int tar_gender = tar_possible_gender[cG];
			int cg_start = tar_gender > 0 ? gender_end[tar_gender - 1] : 0;
			int cg_end = gender_end[tar_gender];

			Arrays.sort(candidates, cg_start, cg_end, comparators[tar_sought_candidate_partner_type_index]);

			// At least seek one partner
			binary_key[Comparator_Candidate_Entry.INDEX_GENDER] = tar_gender;
			binary_key[sought_comparator_partner_type_index] = 1;
			binary_key[Comparator_Candidate_Entry.INDEX_ID] = -1;

			candidateRangeByGender[cG][0] = ~Arrays.binarySearch(candidates, cg_start, cg_end, binary_key,
					comparators[tar_sought_candidate_partner_type_index]);

			binary_key[sought_comparator_partner_type_index] = Integer.MAX_VALUE;
			binary_key[Comparator_Candidate_Entry.INDEX_ID] = Integer.MAX_VALUE;

			candidateRangeByGender[cG][1] = ~Arrays.binarySearch(candidates, cg_start, cg_end, binary_key,
					comparators[tar_sought_candidate_partner_type_index]);

			totalCandidateLength[cG] = candidateRangeByGender[cG][1] - candidateRangeByGender[cG][0];

			if (cG > 0) {
				totalCandidateLength[cG] += totalCandidateLength[cG - 1];
			}
		}

		int[][] target_candidates = new int[totalCandidateLength[totalCandidateLength.length - 1]][];

		int offset = 0;
		for (int cG = 0; cG < candidateRangeByGender.length; cG++) {
			int copyLength = candidateRangeByGender[cG][1] - candidateRangeByGender[cG][0];
			System.arraycopy(candidates, candidateRangeByGender[cG][0], target_candidates, offset, copyLength);
			offset += copyLength;
		}

		return target_candidates;
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
			int[] ent = new int[Comparator_Candidate_Entry.LENGTH_COMPARATOR_BY_PARTNERSHIP_SOUGHT];
			ent[Comparator_Candidate_Entry.INDEX_ID] = person.getId();
			ent[Comparator_Candidate_Entry.INDEX_GENDER] = person.getGenderType();

			ent[Comparator_Candidate_Entry.INDEX_MAX_REG] = (int) person
					.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_REGULAR_PARTNER_12_MONTHS));
			ent[Comparator_Candidate_Entry.INDEX_MAX_CAS] = (int) person
					.getParameter(Integer.toString(Person_Bridging_Pop.FIELD_MAX_CASUAL_PARTNERS_12_MONTHS));

			if (inclCurrent) {
				ent[Comparator_Candidate_Entry.INDEX_CURRENT_REG] = getNumRegularPartnersCurrently(person);
				ent[Comparator_Candidate_Entry.INDEX_CURRENT_CAS] = person.getNumCasualInRecord();

			} else {
				ent[Comparator_Candidate_Entry.INDEX_CURRENT_REG] = 0;
				ent[Comparator_Candidate_Entry.INDEX_CURRENT_CAS] = 0;
			}

			ent[Comparator_Candidate_Entry.INDEX_CURRENT_ANY] = ent[Comparator_Candidate_Entry.INDEX_CURRENT_REG]
					+ ent[Comparator_Candidate_Entry.INDEX_CURRENT_CAS];

			ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_REG] = Math.max(0,
					ent[Comparator_Candidate_Entry.INDEX_MAX_REG] - ent[Comparator_Candidate_Entry.INDEX_CURRENT_ANY]);
			ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_CAS] = Math.max(0,
					ent[Comparator_Candidate_Entry.INDEX_MAX_CAS] - ent[Comparator_Candidate_Entry.INDEX_CURRENT_ANY]);
			ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY] = Math
					.max(ent[Comparator_Candidate_Entry.INDEX_MAX_REG], ent[Comparator_Candidate_Entry.INDEX_MAX_CAS])
					- ent[Comparator_Candidate_Entry.INDEX_CURRENT_ANY];

			ent[Comparator_Candidate_Entry.INDEX_SCHEDULE_LIMIT] = ent[Comparator_Candidate_Entry.INDEX_CURRENT_ANY]
					+ ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY];

			candidates[pI] = ent;
			gender_end[person.getGenderType()]++;

			if (ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY] == 1) {
				singleSoughtCandidate[person.getGenderType()][singleSoughtEnd[person.getGenderType()]] = ent;
				singleSoughtEnd[person.getGenderType()]++;
			} else if (ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY] == 0) {
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
					ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_REG] = 0;
					ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_CAS] = 0;
					ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY] = 0;
					ent[Comparator_Candidate_Entry.INDEX_SCHEDULE_LIMIT] = 0;
					numNotSoughting[g]--;
				}

			}
		}

	}

	private void updateScheduledPopDiff(int[] scheduled_pop_diff_so_far, int[] candidate_cmp_ent,
			float[] categories_values) {
		int num_partner_currently = Math.max(candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_MAX_REG],
				candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_MAX_CAS])
				- candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY];

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
					+ candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_GENDER] * categories_values.length + catAdj;
			scheduled_pop_diff_so_far[adj_partner_dist_index]++;

			if (catAdj_prev > 0) {
				// Subtract previous if it already schedule previously
				adj_partner_dist_index = categories_values.length
						+ candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_GENDER] * categories_values.length
						+ catAdj_prev;
				scheduled_pop_diff_so_far[adj_partner_dist_index]--;

			}
		}
	}

	private void updateCandidateComparatorEntry(int[] candidate_cmp_ent) {
		candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY]--;

		if (candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_REG] > 0) {
			candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_REG]--;
		}
		if (candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_CAS] > 0) {
			candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_CAS]--;
		}
		if (candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_NUM_TO_SOUGHT_ANY] == 0) {
			// No need to sought anymore under current schedule
			candidate_cmp_ent[Comparator_Candidate_Entry.INDEX_SCHEDULE_LIMIT] = 0;
		}

	}

	public void setExport_period_form_partnership_progress(long export_period_form_partnership_progress) {
		this.export_period_form_partnership_progress = export_period_form_partnership_progress;
	}

	private static Comparator_Candidate_Entry generateTargetCandidatesComparator(
			final int comparator_partner_type_index) {
		Comparator_Candidate_Entry target_candidates_comparator = new Comparator_Candidate_Entry(
				~comparator_partner_type_index);

		return target_candidates_comparator;
	}

	private static final class Comparator_Candidate_Entry implements Comparator<int[]> {

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

		final int cmpMethod; // Gender free compare of ~cmpMethod if compMethod < 0

		public Comparator_Candidate_Entry(int cmpMethod) {
			super();
			this.cmpMethod = cmpMethod;
		}

		public int getCmpMethod() {
			return cmpMethod;
		}

		@Override
		public int compare(int[] o1, int[] o2) {
			int cmp;
			if (this.cmpMethod >= 0) {
				cmp = Float.compare(o1[INDEX_GENDER], o2[INDEX_GENDER]);
				if (cmp == 0) {
					cmp = Float.compare(o1[cmpMethod], o2[cmpMethod]);
					if (cmp == 0) {
						cmp = Float.compare(o1[INDEX_ID], o2[INDEX_ID]);
					}
				}
			} else {
				cmp = Integer.compare(o1[~this.cmpMethod], o2[~this.cmpMethod]);
				if (cmp == 0) {
					cmp = Integer.compare(o1[INDEX_ID], o2[INDEX_ID]);
				}
			}
			return cmp;
		}
	}

}
