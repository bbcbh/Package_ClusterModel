package sim;

import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.distribution.PoissonDistribution;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;

public class Runnable_ClusterModel_ContactMap_Generation_MultiMap
		extends Abstract_Runnable_ClusterModel_ContactMap_Generation {

	// Runnable fields
	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP = LENGTH_RUNNABLE_MAP_GEN_FIELD;
	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP
			+ 1;

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST
			+ 1;

	public static final int LENGTH_RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP
			+ 1;

	public static final Object[] DEFAULT_RUNNABLE_MAP_GEN_MULTIMAP_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP
			new int[] { 1000 },
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST
			// double[GRP_NUMBER]{dist}
			new double[][] { new double[] { 18 * AbstractIndividualInterface.ONE_YEAR_INT,
					34 * AbstractIndividualInterface.ONE_YEAR_INT, 0, 1 }, },
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP
			new double[][] { new double[] { 1, 1, 1, 1, 0.7, 0.7, -2.8 }, },

	};

	// Population Index
	private static final int POP_INDEX_GRP = 0;
	private static final int POP_INDEX_ENTER_POP_AGE = POP_INDEX_GRP + 1;
	private static final int POP_INDEX_ENTER_POP_AT = POP_INDEX_ENTER_POP_AGE + 1;
	private static final int POP_INDEX_EXIT_POP_AT = POP_INDEX_ENTER_POP_AT + 1;
	private static final int POP_INDEX_HAS_REG_PARTNER_UNTIL = POP_INDEX_EXIT_POP_AT + 1;
	private static final int LENGTH_POP_ENTRIES = POP_INDEX_HAS_REG_PARTNER_UNTIL + 1;

	// MAPSETTING
	// If DUR_FREQ > 0, then it is number of one partnership to form within snapshot
	// else it is the duration of partnership
	private static final int MAPSETTING_MAP_TYPE = 0;
	private static final int MAPSETTING_SNAP_FREQ = MAPSETTING_MAP_TYPE + 1;
	private static final int MAPSETTING_GRP_INDEX_P1 = MAPSETTING_SNAP_FREQ + 1;
	private static final int MAPSETTING_GRP_INDEX_P2 = MAPSETTING_GRP_INDEX_P1 + 1;
	private static final int MAPSETTING_PROB_HAS_PARTNERSHIP_P1 = MAPSETTING_GRP_INDEX_P2 + 1;
	private static final int MAPSETTING_PROB_HAS_PARTNERSHIP_P2 = MAPSETTING_PROB_HAS_PARTNERSHIP_P1 + 1;
	private static final int MAPSETTING_PARTNERSHIP_FREQ = MAPSETTING_PROB_HAS_PARTNERSHIP_P2 + 1;

	public static final String MAPFILE_FORMAT = "ContactMap_Type_%d_%d.csv"; // Type, Seed,
	public static final String POPSTAT_FORMAT = "POP_STAT_%d"; // Seed

	private RandomGenerator RNG;

	public Runnable_ClusterModel_ContactMap_Generation_MultiMap(long mapSeed) {
		super(mapSeed);

		Object[] newFields = Arrays.copyOf(runnable_fields, LENGTH_RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP);

		for (int i = LENGTH_RUNNABLE_MAP_GEN_FIELD; i < LENGTH_RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP; i++) {
			newFields[i] = DEFAULT_RUNNABLE_MAP_GEN_MULTIMAP_FIELDS[i - LENGTH_RUNNABLE_MAP_GEN_FIELD];
		}

		runnable_fields = newFields;
		RNG = new MersenneTwisterRandomGenerator(mapSeed);

	}

	@Override
	public void run() {
		HashMap<Integer, Object[]> population = new HashMap<>();
		HashMap<Integer, ArrayList<Integer>> active_in_pop = new HashMap<>();

		int nextId = 1;
		int popTime = 0;
		int[] numInGrp = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP];
		double[][] ageDist = (double[][]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST];

		UnivariateInterpolator polator = new LinearInterpolator();

		// Reuse variables
		ArrayList<Integer> active_by_grp;

		// Initialise pop
		for (int g = 0; g < numInGrp.length; g++) {
			UnivariateFunction ageDistFunc = polator.interpolate(Arrays.copyOf(ageDist[g], ageDist[g].length / 2),
					Arrays.copyOfRange(ageDist[g], ageDist[g].length / 2, ageDist[g].length));

			active_by_grp = new ArrayList<>();
			for (int i = 0; i < numInGrp[g]; i++) {
				Object[] newPerson = new Object[LENGTH_POP_ENTRIES];
				newPerson[POP_INDEX_GRP] = g;
				newPerson[POP_INDEX_ENTER_POP_AT] = 1;
				newPerson[POP_INDEX_ENTER_POP_AGE] = (int) Math.round(ageDistFunc.value(RNG.nextDouble()));
				newPerson[POP_INDEX_EXIT_POP_AT] = ((int) ageDist[g][ageDist[g].length / 2]
						- (int) newPerson[POP_INDEX_ENTER_POP_AGE]) + (int) newPerson[POP_INDEX_ENTER_POP_AT];
				newPerson[POP_INDEX_HAS_REG_PARTNER_UNTIL] = -1;
				population.put(nextId, newPerson);

				active_by_grp.add(nextId);
				nextId++;
			}
			active_in_pop.put(g, active_by_grp);

		}

		for (int snapC = 0; snapC < numSnaps; snapC++) {
			popTime += snap_dur;
			int maxMapSnap = -1;

			for (int g = 0; g < numInGrp.length; g++) {
				active_by_grp = active_in_pop.get(g);
				int numRemoved = 0;

				Iterator<Integer> iter = active_by_grp.iterator();
				while (iter.hasNext()) {
					int pId = iter.next();
					Object[] perStat = population.get(pId);
					// Remove age out person
					if ((int) perStat[POP_INDEX_EXIT_POP_AT] <= popTime) {
						iter.remove();
						numRemoved++;
					}
					// Add new person
					while (numRemoved > 0) {
						Object[] newPerson = new Object[LENGTH_POP_ENTRIES];
						newPerson[POP_INDEX_GRP] = g;
						newPerson[POP_INDEX_ENTER_POP_AT] = popTime - RNG.nextInt(snap_dur);
						newPerson[POP_INDEX_ENTER_POP_AGE] = (int) ageDist[g][0];
						newPerson[POP_INDEX_EXIT_POP_AT] = ((int) ageDist[g][ageDist.length / 2]
								- (int) newPerson[POP_INDEX_ENTER_POP_AGE]) + (int) newPerson[POP_INDEX_ENTER_POP_AT];
						newPerson[POP_INDEX_HAS_REG_PARTNER_UNTIL] = -1;
						population.put(nextId, newPerson);

						active_by_grp.add(nextId);
						nextId++;
						numRemoved--;
					}
				}

			}

			for (double[] map_setting : (double[][]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP]) {
				// Generate map
				int gen_map_snap_freq = (int) map_setting[MAPSETTING_SNAP_FREQ];
				maxMapSnap = Math.max(gen_map_snap_freq, maxMapSnap);

				if (snapC > gen_map_snap_freq && snapC % gen_map_snap_freq == 0) {

					ArrayList<int[]> partnership_added = new ArrayList<>();

					int map_type = (int) map_setting[MAPSETTING_MAP_TYPE];
					int grp_index_p1 = (int) map_setting[MAPSETTING_GRP_INDEX_P1];
					int grp_index_p2 = (int) map_setting[MAPSETTING_GRP_INDEX_P2];

					PoissonDistribution dur_freq_dist = new PoissonDistribution(RNG,
							Math.abs(map_setting[MAPSETTING_PARTNERSHIP_FREQ]), PoissonDistribution.DEFAULT_EPSILON,
							PoissonDistribution.DEFAULT_MAX_ITERATIONS);

					// Candidate list = int[]{PID, NUMBER_OF_PANTERSHIP_FROM}
					ArrayList<int[]> candidate_p1 = new ArrayList<>();
					ArrayList<int[]> candidate_p2 = candidate_p1; // Shared candidate list

					if (grp_index_p1 != grp_index_p2) {
						candidate_p2 = new ArrayList<>();
					}

					// Fill candidate list
					for (int g = 0; g < numInGrp.length; g++) {
						if ((grp_index_p1 & 1 << g) != 0) {
							active_by_grp = new ArrayList<>(active_in_pop.get(g));

							int numSelectP1 = (int) Math
									.round(active_by_grp.size() * map_setting[MAPSETTING_PROB_HAS_PARTNERSHIP_P1]);
							int numSelectP2 = (int) Math
									.round(active_by_grp.size() * map_setting[MAPSETTING_PROB_HAS_PARTNERSHIP_P2]);

							if (map_setting[MAPSETTING_PARTNERSHIP_FREQ] < 0) {
								// Filter out those already within regular partnership
								for (Integer pid : active_by_grp) {
									int reg_part_until = (int) population.get(pid)[POP_INDEX_HAS_REG_PARTNER_UNTIL];
									if (reg_part_until > popTime - gen_map_snap_freq * snap_dur) {
										numSelectP1--;
										numSelectP2--;
									}
								}
							}

							for (int i = 0; i < active_by_grp.size() && numSelectP1 > 0; i++) {
								if (RNG.nextInt(active_by_grp.size() - i) < numSelectP1) {
									int dur_freq_val = dur_freq_dist.sample();
									if (dur_freq_val > 0) {
										candidate_p1.add(new int[] { active_by_grp.get(i), dur_freq_val });
										numSelectP1--;
									}
								}
							}

							if (grp_index_p1 != grp_index_p2) {
								for (int i = 0; i < active_by_grp.size() && numSelectP1 > 0; i++) {
									if (RNG.nextInt(active_by_grp.size() - i) < numSelectP2) {
										int dur_freq_val = dur_freq_dist.sample();
										if (dur_freq_val > 0) {
											candidate_p2.add(new int[] { active_by_grp.get(i), dur_freq_val });
											numSelectP2--;
										}
									}
								}
							}
						}
					}

					// Form partnerships
					while (candidate_p1.size() > 0 && candidate_p2.size() > 0
							&& candidate_p1.size() + candidate_p2.size() > 1) {

						int[][] selected_candidates = new int[2][];

						selected_candidates[0] = candidate_p1.remove(RNG.nextInt(candidate_p1.size()));
						selected_candidates[1] = candidate_p2.remove(RNG.nextInt(candidate_p2.size()));

						if (map_setting[MAPSETTING_PARTNERSHIP_FREQ] > 0) {
							partnership_added.add(new int[] { selected_candidates[0][0], selected_candidates[1][0],
									popTime - RNG.nextInt(gen_map_snap_freq), 1 });
							for (int p = 0; p < selected_candidates.length; p++) {
								selected_candidates[p][1]--;

							}
						} else {
							int dur = Math.min(selected_candidates[0][1], selected_candidates[0][1]);
							partnership_added.add(new int[] { selected_candidates[0][0], selected_candidates[1][0],
									popTime - RNG.nextInt(gen_map_snap_freq), dur });
							for (int p = 0; p < selected_candidates.length; p++) {
								selected_candidates[p][1] = 0;
								population.get(selected_candidates[p][0])[POP_INDEX_HAS_REG_PARTNER_UNTIL] = dur
										+ popTime;
							}
						}
						// Add candidate back if needed
						if (selected_candidates[0][1] > 0) {
							candidate_p1.add(selected_candidates[0]);
						}
						if (selected_candidates[1][1] > 0) {
							candidate_p2.add(selected_candidates[1]);
						}
					}

					// Export partnership
					int[][] partnership_added_arr = partnership_added.toArray(new int[0][]);
					partnership_added.sort(new Comparator<int[]>() {
						@Override
						public int compare(int[] o1, int[] o2) {
							int res = Integer.compare(o1[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME],
									o2[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);

							if (res == 0) {
								for (int i = 0; i < Math.min(o1.length, o2.length) && res == 0; i++) {
									res = Integer.compare(o1[i], o2[i]);
								}
							}
							return res;
						}
					});

					try {
						String genMapFilename = String.format(MAPFILE_FORMAT, map_type, getMapSeed());
						File genMapFile = getTargetFile(genMapFilename);
						PrintWriter pWri = new PrintWriter(new FileWriter(genMapFile, true));
						for (int[] partnership_added_ent : partnership_added_arr) {
							pWri.println(String.format("%d,%d,%d,%d\n",
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2],
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME],
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]));

						}
						pWri.close();

					} catch (IOException ex) {
						ex.printStackTrace(System.err);
					}

				}

			}

			if (snapC % maxMapSnap == 0) {
				exportPopulationToFile(population);
			}

		}
		// Finishing
		exportPopulationToFile(population);
	}

	public void exportPopulationToFile(HashMap<Integer, Object[]> population) {
		try {
			String exportFileName = String.format(POPSTAT_FORMAT, mapSeed);
			File exportPopFile = getTargetFile(exportFileName);
			PrintWriter pWri_exportPop = new PrintWriter(exportPopFile);

			pWri_exportPop.println("ID,GRP,ENTER_POP_AGE,ENTER_POP_AT,EXIT_POP_AT,HAS_REG_PARTNER_UNTIL");
			for (Integer id : population.keySet()) {
				pWri_exportPop.printf("%d", id);
				Object[] ent = population.get(id);
				for (Object obj : ent) {
					pWri_exportPop.print(',');
					pWri_exportPop.print(obj.toString());
				}
				pWri_exportPop.println();
			}
			pWri_exportPop.close();
		} catch (IOException ex) {
			ex.printStackTrace(System.err);
		}
	}

	public File getTargetFile(String inputFilename) throws IOException {
		File genMapFile = new File(baseDir, inputFilename);

		if (genMapFile.exists()) {
			File[] delFile = baseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.getName().startsWith(inputFilename + "_");
				}
			});
			for (File d : delFile) {
				Files.delete(d.toPath());
			}
			genMapFile.renameTo(new File(baseDir, String.format("%s_%d", inputFilename, System.currentTimeMillis())));
			genMapFile = new File(baseDir, inputFilename);
		}
		return genMapFile;
	}

}
