package sim;

import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.distribution.PoissonDistribution;

import person.AbstractIndividualInterface;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import util.Util_7Z_CSV_Entry_Extract_Callable;

public class Runnable_ClusterModel_ContactMap_Generation_MultiMap
		extends Abstract_Runnable_ClusterModel_ContactMap_Generation {

	// Runnable fields
	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP = Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD;
	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGE_DIST = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP
			+ 1;

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGE_DIST
			+ 1;

	public static final int LENGTH_RUNNABLE_MAP_GEN_FIELD = RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP
			+ 1;

	public static final Object[] DEFAULT_RUNNABLE_MAP_GEN_MULTIMAP_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP
			new int[] { 1000 },
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGEING_DIST
			// double[GRP_NUMBER]{dist}
			new double[][] { new double[] { 18 * AbstractIndividualInterface.ONE_YEAR_INT,
					34 * AbstractIndividualInterface.ONE_YEAR_INT, 0, 1 }, },
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP
			new double[][] { new double[] { 1, 30, 1, 1, 0.7, 0.7, -2.8 }, },

	};

	// MAPSETTING
	// If DUR_FREQ > 0, then it is number of one partnership to form within snapshot
	// else it is the duration of partnership
	private static final int MAPSETTING_MAP_TYPE = 0;
	private static final int MAPSETTING_FREQ = MAPSETTING_MAP_TYPE + 1;
	private static final int MAPSETTING_GRP_INDEX_P1 = MAPSETTING_FREQ + 1;
	private static final int MAPSETTING_GRP_INDEX_P2 = MAPSETTING_GRP_INDEX_P1 + 1;
	private static final int MAPSETTING_PROB_HAS_PARTNERSHIP_P1 = MAPSETTING_GRP_INDEX_P2 + 1;
	private static final int MAPSETTING_PROB_HAS_PARTNERSHIP_P2 = MAPSETTING_PROB_HAS_PARTNERSHIP_P1 + 1;
	private static final int MAPSETTING_PARTNERSHIP_FREQ = MAPSETTING_PROB_HAS_PARTNERSHIP_P2 + 1;

	public static final String MAPFILE_FORMAT = "ContactMap_Type_%d_%d.csv"; // Type, Seed,
	public static final String POPSTAT_FORMAT = "POP_STAT_%d.csv"; // Seed

	private static String OUTPUTMSG_FORMAT = "%d:RelMap #%d generated with %d new partnerships added. Time req. = %.3f seconds\n";

	protected RandomGenerator RNG;

	public Runnable_ClusterModel_ContactMap_Generation_MultiMap(long mapSeed) {
		super(mapSeed);

		Object[] newFields = Arrays.copyOf(runnable_fields, LENGTH_RUNNABLE_MAP_GEN_FIELD);

		for (int i = Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD; i < LENGTH_RUNNABLE_MAP_GEN_FIELD; i++) {
			newFields[i] = DEFAULT_RUNNABLE_MAP_GEN_MULTIMAP_FIELDS[i
					- Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD];
		}

		runnable_fields = newFields;
		RNG = new MersenneTwisterRandomGenerator(mapSeed);

	}

	@Override
	public void run() {
		HashMap<Integer, Object[]> population = new HashMap<>();
		HashMap<Integer, ArrayList<Integer>> active_in_pop = new HashMap<>();
		// Fields
		int[] contactMapValidRange = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE];
		final float noBridge = (float) runnable_fields[RUNNABLE_FIELD_NO_BRIDGE];
		int[] numInGrp = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_NUMBER_OF_GRP];
		double[][] ageDist = (double[][]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_AGE_DIST];

		int nextId = 1;
		int popTime = contactMapValidRange[0];
		long lastSnapTime = -1;
		UnivariateInterpolator polator = new LinearInterpolator();

		// Reuse variables
		ArrayList<Integer> active_by_grp;

		// Check for previous simulation
		File outputFile = new File(baseDir,
				String.format(Simulation_ClusterModelGeneration.FILENAME_FORMAT_OUTPUT, mapSeed));

		if (outputFile.exists()) {
			try {
				String patternStr = "(\\d+):.*";
				Pattern lastLinePattern = Pattern.compile(patternStr);
				String[] lines = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(outputFile);
				int lastlinePt = lines.length - 1;

				Matcher m = lastLinePattern.matcher(lines[lastlinePt]);

				while (!m.matches() && lastlinePt > 0) {
					lastlinePt--;
					m = lastLinePattern.matcher(lines[lastlinePt]);
				}

				if (m.matches()) {
					popTime = Integer.parseInt(m.group(1));

					// Update population
					File popFile = new File(baseDir, String.format(POPSTAT_FORMAT, mapSeed));
					if (popFile.exists()) {
						String[] popLines = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(popFile);
						popFile.renameTo(new File(baseDir,
								String.format(popFile.getName() + "_%d", System.currentTimeMillis())));
						for (int i = 1; i < popLines.length; i++) {
							String[] ent = popLines[i].split(",");
							Object[] newPerson = new Object[Abstract_Runnable_ClusterModel.LENGTH_POP_ENTRIES];
							int pid = Integer.parseInt(ent[0].toString());
							for (int j = 0; j < newPerson.length; j++) {
								newPerson[j] = Integer.parseInt(ent[j + 1]);
							}

							if (((Integer) newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AT]) <= popTime) {
								population.put(pid, newPerson);

								active_by_grp = active_in_pop
										.get(newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_GRP]);
								if (active_by_grp == null) {
									active_by_grp = new ArrayList<Integer>();
									active_in_pop.put((Integer) newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_GRP],
											active_by_grp);
								}
								active_by_grp.add(pid);
								nextId = Math.max(nextId, pid + 1);
							}

						}
					}

					// Remove edges from map
					int numMap = ((double[][]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP]).length;
					for (int mapId = 0; mapId < numMap; mapId++) {
						File mapFile = new File(baseDir, String.format(MAPFILE_FORMAT, mapId, mapSeed));
						if (mapFile.exists()) {
							String[] edges = Util_7Z_CSV_Entry_Extract_Callable.extracted_lines_from_text(mapFile);
							mapFile.renameTo(new File(baseDir,
									String.format(mapFile.getName() + "_%d", System.currentTimeMillis())));

							// New map file
							mapFile = new File(baseDir, String.format(MAPFILE_FORMAT, mapId, mapSeed));
							PrintWriter pWri = new PrintWriter(mapFile);

							for (String edge : edges) {
								String[] ent = edge.split(",");
								int edgeFormAt = Integer
										.parseInt(ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
								if (edgeFormAt <= popTime) {
									int edgeDur = Integer
											.parseInt(ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]);
									if (edgeDur > 1) {
										for (int pIdPt : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
												Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
											int pId = Integer.parseInt(ent[pIdPt]);
											Object[] person = population.get(pId);
											person[Abstract_Runnable_ClusterModel.POP_INDEX_HAS_REG_PARTNER_UNTIL] = Math
													.max((Integer) person[Abstract_Runnable_ClusterModel.POP_INDEX_HAS_REG_PARTNER_UNTIL],
															edgeFormAt + edgeDur);
										}
									}

									pWri.println(edge);
								}
							}
							pWri.close();
						}
					}

				}

			} catch (IOException e) {
				e.printStackTrace(System.err);
			}

		}

		if (population.isEmpty()) {
			// Initialise pop
			for (int g = 0; g < numInGrp.length; g++) {
				UnivariateFunction ageDistFunc = polator.interpolate(
						Arrays.copyOfRange(ageDist[g], ageDist[g].length / 2, ageDist[g].length),
						Arrays.copyOf(ageDist[g], ageDist[g].length / 2));

				active_by_grp = new ArrayList<>();
				for (int i = 0; i < numInGrp[g]; i++) {
					Object[] newPerson = new Object[Abstract_Runnable_ClusterModel.LENGTH_POP_ENTRIES];
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_GRP] = g;
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AT] = 1;
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AGE] = (int) Math
							.round(ageDistFunc.value(RNG.nextDouble()));
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_EXIT_POP_AT] = ((int) ageDist[g][ageDist[g].length
							/ 2 - 1] - (int) newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AGE])
							+ (int) newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AT];
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_HAS_REG_PARTNER_UNTIL] = -1;
					population.put(nextId, newPerson);

					active_by_grp.add(nextId);
					nextId++;
				}
				active_in_pop.put(g, active_by_grp);

			}
		}
		// First export
		exportPopulationToFile(population);

		int numSnapsSim = Math.max(numSnaps, Math.round(contactMapValidRange[1] / snap_dur));

		for (int snapC = 0; snapC < numSnapsSim; snapC++) {
			popTime += snap_dur;

			// Store who are the person added
			HashMap<Integer, Object[]> justAddedPerson = new HashMap<>();

			for (int g = 0; g < numInGrp.length; g++) {
				active_by_grp = active_in_pop.get(g);
				int numRemoved = 0;
				Iterator<Integer> iter = active_by_grp.iterator();
				// Remove age out person
				while (iter.hasNext()) {
					int pId = iter.next();
					Object[] perStat = population.get(pId);
					if ((int) perStat[Abstract_Runnable_ClusterModel.POP_INDEX_EXIT_POP_AT] <= popTime) {
						iter.remove();
						numRemoved++;
					}
				}

				while (numRemoved > 0) {
					Object[] newPerson = new Object[Abstract_Runnable_ClusterModel.LENGTH_POP_ENTRIES];
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_GRP] = g;
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AT] = popTime - RNG.nextInt(snap_dur);
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AGE] = (int) ageDist[g][0];
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_EXIT_POP_AT] = ((int) ageDist[g][ageDist[g].length
							/ 2 - 1] - (int) newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AGE])
							+ (int) newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AT];
					newPerson[Abstract_Runnable_ClusterModel.POP_INDEX_HAS_REG_PARTNER_UNTIL] = -1;
					population.put(nextId, newPerson);
					justAddedPerson.put(nextId, newPerson);
					active_by_grp.add(nextId);
					nextId++;
					numRemoved--;
				}

			}

			// Add new person
			if (!justAddedPerson.isEmpty()) {
				exportPopulationToFile(justAddedPerson, true);
			}

			for (double[] map_setting : (double[][]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_MULTIMAP_PARTNERSHIP_BY_SNAP]) {
				// Generate map
				int gen_map_freq = (int) map_setting[MAPSETTING_FREQ];

				if ((popTime - contactMapValidRange[0]) >= gen_map_freq
						&& (popTime - contactMapValidRange[0]) % gen_map_freq == 0) {

					ArrayList<int[]> partnership_added = new ArrayList<>();

					int map_type = (int) map_setting[MAPSETTING_MAP_TYPE];
					int grp_index_p1 = (int) map_setting[MAPSETTING_GRP_INDEX_P1];
					int grp_index_p2 = (int) map_setting[MAPSETTING_GRP_INDEX_P2];
					long tic = System.currentTimeMillis();

					// Assume minimum casual of 1 in last month, or 1 day duration for reg
					// partnership
					PoissonDistribution dur_freq_dist = new PoissonDistribution(RNG,
							Math.abs(map_setting[MAPSETTING_PARTNERSHIP_FREQ]) - 1, PoissonDistribution.DEFAULT_EPSILON,
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
							active_by_grp = active_in_pop.get(g);

							int numSelectP1 = (int) Math
									.round(active_by_grp.size() * map_setting[MAPSETTING_PROB_HAS_PARTNERSHIP_P1]);
							int numSelectP2 = (int) Math
									.round(active_by_grp.size() * map_setting[MAPSETTING_PROB_HAS_PARTNERSHIP_P2]);

							if (map_setting[MAPSETTING_PARTNERSHIP_FREQ] < 0) {
								// Filter out those already within regular partnership
								for (Integer pid : active_by_grp) {
									int reg_part_until = (int) population
											.get(pid)[Abstract_Runnable_ClusterModel.POP_INDEX_HAS_REG_PARTNER_UNTIL];
									if (reg_part_until > popTime - gen_map_freq) {
										numSelectP1--;
										numSelectP2--;
									}
								}
							}

							for (int i = 0; i < active_by_grp.size() && numSelectP1 > 0; i++) {
								if (RNG.nextInt(active_by_grp.size() - i) < numSelectP1) {
									int dur_freq_val = dur_freq_dist.sample() + 1;
									candidate_p1.add(new int[] { active_by_grp.get(i), dur_freq_val });
									numSelectP1--;

								}
							}

							if (grp_index_p1 != grp_index_p2) {
								for (int i = 0; i < active_by_grp.size() && numSelectP1 > 0; i++) {
									if (RNG.nextInt(active_by_grp.size() - i) < numSelectP2) {
										int dur_freq_val = dur_freq_dist.sample() + 1;
										candidate_p2.add(new int[] { active_by_grp.get(i), dur_freq_val });
										numSelectP2--;

									}
								}
							}
						}
					}

					// Form partnerships
					while ((grp_index_p1 == grp_index_p2) ? candidate_p1.size() > 1
							: candidate_p1.size() > 0 && candidate_p2.size() > 0) {

						int[][] selected_candidates = new int[2][];

						selected_candidates[0] = candidate_p1.remove(RNG.nextInt(candidate_p1.size()));
						selected_candidates[1] = candidate_p2.remove(RNG.nextInt(candidate_p2.size()));

						if (map_setting[MAPSETTING_PARTNERSHIP_FREQ] > 0) {
							partnership_added.add(new int[] { selected_candidates[0][0], selected_candidates[1][0],
									popTime - RNG.nextInt(gen_map_freq), 1 });
							for (int p = 0; p < selected_candidates.length; p++) {
								selected_candidates[p][1]--;
							}
						} else {
							int dur = Math.min(selected_candidates[0][1], selected_candidates[0][1]) + 1;
							partnership_added.add(new int[] { selected_candidates[0][0], selected_candidates[1][0],
									popTime - RNG.nextInt(gen_map_freq), dur });
							for (int p = 0; p < selected_candidates.length; p++) {
								selected_candidates[p][1] = 0;
								population.get(
										selected_candidates[p][0])[Abstract_Runnable_ClusterModel.POP_INDEX_HAS_REG_PARTNER_UNTIL] = dur
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
						for (int[] partnership_added_ent : partnership_added) {
							pWri.printf(String.format("%d,%d,%d,%d\n",
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2],
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME],
									partnership_added_ent[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]));

						}
						pWri.close();

						if (printStatus != null) {
							for (PrintStream out : printStatus) {

								out.printf(OUTPUTMSG_FORMAT, popTime, map_type, partnership_added.size(),
										(System.currentTimeMillis() - tic) / 1000f);
							}
						}

					} catch (IOException ex) {
						ex.printStackTrace(System.err);
					}

				}

			}
			
		}
		// Finishing
		exportPopulationToFile(population);
	}

	private void exportPopulationToFile(HashMap<Integer, Object[]> population) {
		exportPopulationToFile(population, false);
	}

	private void exportPopulationToFile(HashMap<Integer, Object[]> population, boolean append) {
		try {
			String exportFileName = String.format(POPSTAT_FORMAT, mapSeed);
			File exportPopFile = getTargetFile(exportFileName);
			PrintWriter pWri_exportPop = new PrintWriter(new FileWriter(exportPopFile, append));
			if (!exportPopFile.exists() || !append) {
				pWri_exportPop.println("ID,GRP,ENTER_POP_AGE,ENTER_POP_AT,EXIT_POP_AT,HAS_REG_PARTNER_UNTIL");
			}
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

	protected File getTargetFile(String inputFilename) {
		File genMapFile = new File(baseDir, inputFilename);

		if (genMapFile.exists() && !isSpace_save()) {
			File[] delFile = baseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.getName().startsWith(inputFilename + "_");
				}
			});
			try {
				for (File d : delFile) {
					Files.delete(d.toPath());
				}

				Files.copy(genMapFile.toPath(),
						new File(baseDir, String.format("%s_%d", inputFilename, System.currentTimeMillis())).toPath());
			} catch (IOException ex) {
				ex.printStackTrace(System.err);
			}			
		}
		return genMapFile;
	}

	@Override
	public void setRunnable_fields(Object[] simFields) {
		for (int f = 0; f < getRunnable_fields().length; f++) {
			if (simFields[Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD + f] != null) {
				getRunnable_fields()[f] = simFields[f + Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD];
			}
		}

	}

}
