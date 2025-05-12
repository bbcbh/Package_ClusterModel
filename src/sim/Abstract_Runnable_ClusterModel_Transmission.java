package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.IntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.ArrayUtilsRandomGenerator;
import util.PropValUtils;
import util.Util_7Z_CSV_Entry_Extract_Callable;

public abstract class Abstract_Runnable_ClusterModel_Transmission extends Abstract_Runnable_ClusterModel {

	public static final int ACT_INDEX_GENITAL = 0;
	public static final int ACT_INDEX_ANAL = ACT_INDEX_GENITAL + 1;
	public static final int ACT_INDEX_FELLATIO = ACT_INDEX_ANAL + 1;
	public static final int ACT_INDEX_RIMMING = ACT_INDEX_FELLATIO + 1;
	public static final int ACT_INDEX_KISSING = ACT_INDEX_RIMMING + 1;
	public static final int LENGTH_ACT = ACT_INDEX_KISSING + 1;

	public static final int SITE_VAGINA = 0;
	public static final int SITE_PENIS = SITE_VAGINA + 1;
	public static final int SITE_RECTUM = SITE_PENIS + 1;
	public static final int SITE_OROPHARYNX = SITE_RECTUM + 1;
	public static final int LENGTH_SITE = SITE_OROPHARYNX + 1;

	public static final int RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ = 0;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE = RUNNABLE_FIELD_TRANSMISSION_ACT_FREQ + 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD = RUNNABLE_FIELD_TRANSMISSION_TRANSMISSION_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD = RUNNABLE_FIELD_TRANSMISSION_INFECTIOUS_PERIOD
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SYM_RATE = RUNNABLE_FIELD_TRANSMISSION_STAGE_PERIOD + 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS = RUNNABLE_FIELD_TRANSMISSION_SYM_RATE
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES = RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM = RUNNABLE_FIELD_TRANSMISSION_TESTING_RATE_BY_RISK_CATEGORIES
			+ 1;
	public static final int RUNNABLE_FIELD_TRANSMISSION_DX_TEST_PROPERTIES = RUNNABLE_FIELD_TRANSMISSION_SOUGHT_TEST_PERIOD_BY_SYM
			+ 1;
	protected static final int RUNNABLE_OFFSET = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
			+ +Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
			+ Simulation_ClusterModelTransmission.LENGTH_SIM_MAP_TRANSMISSION_FIELD;

	protected final int[] cUMULATIVE_POP_COMPOSITION;
	protected final ContactMap bASE_CONTACT_MAP;
	protected final int nUM_TIME_STEPS_PER_SNAP;
	protected final int nUM_SNAP;
	protected final long sIM_SEED;
	protected final long cMAP_SEED;

	protected PrintStream print_progress = null;
	protected int simSetting = 1;
	protected RandomGenerator RNG;
	protected ArrayList<Integer[]> edges_list;

	// Multi contact map / large map support
	protected File[] contactMap_files = null;
	protected String[] contactMap_nextString = null;
	protected BufferedReader[] contactMap_reader = null;

	// PopStat
	protected HashMap<Integer, String[]> pop_stat = null;

	protected HashMap<Integer, HashMap<Integer, String>> propSwitch_map;

	protected transient HashMap<Integer, Integer> risk_cat_map;
	protected transient int firstSeedTime = Integer.MAX_VALUE;
	protected transient HashMap<String, Object> sim_output = null;

	// Format:
	// float[][] { edge_duration_to_be_kept, kept_probability }, or empty string to
	// set back to default
	public static final int PROPSWITCH_MAP_KEPT_EDGE_BY_DURATION_INDEX = -1;
	protected float[][] cMap_Kept_Edge_By_Duration_Setting = null;

	private static final int cMAP_KEPT_EDGE_BY_DURATION_DUR_LIST = 0;
	private static final int cMAP_KEPT_EDGE_BY_DURATION_PROB_LIST = cMAP_KEPT_EDGE_BY_DURATION_DUR_LIST + 1;

	protected static final String popCompositionKey = Simulation_ClusterModelTransmission.POP_PROP_INIT_PREFIX
			+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);

	protected Properties baseProp; // From simSpecificSim.prop

	// Non_mapped encounter
	// To be initialise during updateCMap
	protected float[] non_mapped_encounter_prob = null;
	protected int[] non_mapped_encounter_target_gender = null;
	protected ArrayList<Integer>[] non_map_candidate_id_seeker = null;
	protected ArrayList<Integer>[] non_map_candidate_id_target = null;

	protected ArrayList<String[]> non_map_edges_pre_gen = null;
	protected int non_map_edges_pre_gen_last_edge_pt = 0;

	protected Integer[][] non_map_edges_store = null;
	protected RandomGenerator RNG_IMPORT = null;
	protected RandomGenerator RNG_NM = null;
	protected IntegerDistribution[] dist_nm_candidate = null;

	private int NON_MAP_EDGES_STORE_PT = 0;
	private int NON_MAP_EDGES_STORE_MAX = 100000;

	protected static String[] exportFileFormat = new String[] { "Export_RNG_%d_%d_%d.obj",
			"Export_SimOutput_%d_%d_%d.obj", };

	protected static final int EXPORT_RNG = 0;
	protected static final int EXPORT_SIMOUTPUT = EXPORT_RNG + 1;
	protected static final int LENGTH_EXPORT = EXPORT_SIMOUTPUT + 1;	
	protected boolean skipStateGen = false;

	protected int importedAtTime = -1;
	
	protected final int SIM_OFFSET = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
			+ Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD;

	public Abstract_Runnable_ClusterModel_Transmission(long cMap_seed, long sim_seed, ContactMap base_cMap,
			Properties prop) {

		this(cMap_seed, sim_seed,
				(int[]) PropValUtils.propStrToObject(prop.getProperty(popCompositionKey), int[].class), base_cMap,
				Integer.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])),
				Integer.parseInt(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])));

		this.baseProp = prop;
	}

	public Abstract_Runnable_ClusterModel_Transmission(long cMap_seed, long sim_seed, int[] pop_composition,
			ContactMap base_cMap, int numTimeStepsPerSnap, int numSnap) {
		super();
		this.cUMULATIVE_POP_COMPOSITION = new int[pop_composition.length];
		int offset = 0;

		for (int g = 0; g < this.cUMULATIVE_POP_COMPOSITION.length; g++) {
			this.cUMULATIVE_POP_COMPOSITION[g] = offset + pop_composition[g];
			offset += pop_composition[g];
		}

		this.bASE_CONTACT_MAP = base_cMap;
		this.nUM_TIME_STEPS_PER_SNAP = numTimeStepsPerSnap;
		this.nUM_SNAP = numSnap;
		this.cMAP_SEED = cMap_seed;
		this.sIM_SEED = sim_seed;

		RNG = new MersenneTwisterRandomGenerator(sIM_SEED);
		RNG_NM = new MersenneTwisterRandomGenerator(sIM_SEED);
		RNG_IMPORT = new MersenneTwisterRandomGenerator(sIM_SEED);
	}

	public int getImportedAtTime() {
		return importedAtTime;
	}

	@SuppressWarnings("unchecked")
	public void importRunnableTransmission(int time_pt) throws IOException, ClassNotFoundException {
		File objFile;
		ObjectInputStream objStream;
		// RNG
		objFile = new File(baseDir, String.format(exportFileFormat[EXPORT_RNG], cMAP_SEED, sIM_SEED, time_pt));
		objStream = new ObjectInputStream(new FileInputStream(objFile));
		RNG = (RandomGenerator) objStream.readObject();
		objStream.close();

		// SIM_OUPUT
		objFile = new File(baseDir, String.format(exportFileFormat[EXPORT_SIMOUTPUT], cMAP_SEED, sIM_SEED, time_pt));
		objStream = new ObjectInputStream(new FileInputStream(objFile));
		sim_output = (HashMap<String, Object>) objStream.readObject();
		objStream.close();

		importedAtTime = time_pt;
	}

	public void exportRunnableTransmission(int time_pt) throws IOException {
		File objFile, archiveFile;
		ObjectOutputStream objStream;
		// RNG
		objFile = new File(baseDir, String.format(exportFileFormat[EXPORT_RNG], cMAP_SEED, sIM_SEED, time_pt));
		objStream = new ObjectOutputStream(new FileOutputStream(objFile));
		objStream.writeObject(RNG);
		objStream.close();

		archiveFile = new File(baseDir, objFile.getName() + ".7z");
		util.Util_7Z_CSV_Entry_Extract_Callable.zipFile(new File[] { objFile }, archiveFile, false);
		Files.delete(objFile.toPath());

		// SIM_OUPUT
		objFile = new File(baseDir, String.format(exportFileFormat[EXPORT_SIMOUTPUT], cMAP_SEED, sIM_SEED, time_pt));
		objStream = new ObjectOutputStream(new FileOutputStream(objFile));
		objStream.writeObject(sim_output);
		objStream.close();

		archiveFile = new File(baseDir, objFile.getName() + ".7z");
		util.Util_7Z_CSV_Entry_Extract_Callable.zipFile(new File[] { objFile }, archiveFile, false);
		Files.delete(objFile.toPath());
	}		

	public void setSkipStateGen(boolean skipStateGen) {
		this.skipStateGen = skipStateGen;
	}

	public void importExportedTransmissionStates(String[] fileformats) {
		ArrayList<Integer> timePt_Collection = null;
		ArrayList<File> extracted_file = new ArrayList<>();

		for (String fileformat : fileformats) {
			String pat_format = fileformat.replaceFirst("%d", Long.toString(cMAP_SEED));
			pat_format = pat_format.replaceFirst("%d", Long.toString(sIM_SEED));
			pat_format = pat_format.replaceFirst("%d", "(\\\\d+)");
			Pattern pattern_sel_zip = Pattern.compile(pat_format + ".7z");

			File[] candidate = baseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pattern_sel_zip.matcher(pathname.getName()).matches();
				}
			});
			// Replace with the unzip version if exist
			for (int f = 0; f < candidate.length; f++) {
				String rawFileName = candidate[f].getName().substring(0, candidate[f].getName().length() - 3);
				File rawFile = new File(baseDir, rawFileName);
				if (rawFile.exists()) {
					candidate[f] = rawFile;
				}
			}

			Pattern pattern_sel = Pattern.compile(pat_format + ".*");
			ArrayList<Integer> timePt_Collection_fitlered = new ArrayList<>();
			for (File f : candidate) {
				Matcher m = pattern_sel.matcher(f.getName());
				m.matches();
				int time_pt = Integer.parseInt(m.group(1));
				int pos;

				if (timePt_Collection == null || Collections.binarySearch(timePt_Collection, time_pt) >= 0) {
					boolean validFile = f.length() > 0;
					if (f.getName().endsWith(".7z")) {
						try {
							extracted_file.addAll(Util_7Z_CSV_Entry_Extract_Callable.unzipFile(f, baseDir));
						} catch (IOException e) {
							validFile = false;
						}
					}
					if (validFile) {
						pos = Collections.binarySearch(timePt_Collection_fitlered, time_pt);
						if (pos < 0) {
							timePt_Collection_fitlered.add(~pos, time_pt);
						}
					}
				}
			}
			timePt_Collection = timePt_Collection_fitlered;
		}

		// Select best time point

		if (timePt_Collection.size() > 0) {
			int bestTimePt = timePt_Collection.get(timePt_Collection.size() - 1);

			// Remove non-used files
			for (String fileformat : fileformats) {
				String pat_format = fileformat.replaceFirst("%d", Long.toString(cMAP_SEED));
				pat_format = pat_format.replaceFirst("%d", Long.toString(sIM_SEED));
				pat_format = pat_format.replaceFirst("%d", "(\\\\d+)");
				Pattern pattern_sel = Pattern.compile(pat_format + ".*");

				File[] candidate = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						Matcher m = pattern_sel.matcher(pathname.getName());
						boolean res = m.matches();
						if (res) {
							res = Integer.parseInt(m.group(1)) != bestTimePt;
						}
						return res;
					}
				});
				for (File f : candidate) {
					try {
						Files.delete(f.toPath());
					} catch (IOException e) {
						e.printStackTrace(System.err);
					}
				}
			}

			try {
				long tic = System.currentTimeMillis();
				try {
					importRunnableTransmission(bestTimePt);

				} catch (IOException | ClassNotFoundException e) {
					e.printStackTrace(System.err);
				}

				if (print_progress != null && runnableId != null) {
					print_progress.printf("Thread <%s>: Population stat imported at T = %d . Time req. = %.3fs\n",
							runnableId, bestTimePt, (System.currentTimeMillis() - tic) / 1000.0);

				}

			} catch (Exception e) {
				e.printStackTrace(System.err);
			}

		}

		// Remove extracted file
		for (File f : extracted_file) {
			if (f.exists()) {
				try {
					Files.delete(f.toPath());
				} catch (IOException e) {
					e.printStackTrace(System.err);
				}
			}
		}

	}

	public static int getGenderType(Integer personId, int[] cumul_pop_comp) {
		int index = Arrays.binarySearch(cumul_pop_comp, personId);

		if (index < 0) {
			return ~index;
		} else {
			return index;
		}
	}

	public abstract void initialse();

	public Properties getSim_prop() {
		return baseProp;
	}

	public HashMap<Integer, String[]> getPop_stat() {
		return pop_stat;
	}

	public void setPop_stat(HashMap<Integer, String[]> pop_stat) {
		this.pop_stat = pop_stat;
	}

	public Integer[] getCurrentPopulationPId(int time) {
		if (pop_stat == null) {
			if (bASE_CONTACT_MAP != null) {
				// All valid
				return bASE_CONTACT_MAP.vertexSet().toArray(new Integer[bASE_CONTACT_MAP.vertexSet().size()]);
			} else {
				System.err.println("Warning - Population undefined.");
				return null;
			}
		} else {
			ArrayList<Integer> res = new ArrayList<>();
			for (Integer pid : pop_stat.keySet()) {
				String[] popEnt = pop_stat.get(pid);
				if (Integer.parseInt(popEnt[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AT]) <= time
						&& time < Integer.parseInt(popEnt[Abstract_Runnable_ClusterModel.POP_INDEX_EXIT_POP_AT])) {
					res.add(pid);
				}
			}
			return res.toArray(new Integer[0]);
		}
	}
	
	public boolean isValidPerson(int pid, int[] timeRange) {
		if (pop_stat == null) {
			if (bASE_CONTACT_MAP != null) {
				return bASE_CONTACT_MAP.containsVertex(pid);
			}else {
				System.err.println("Warning - Population undefined.");
				return false;
			}			
		}else {			
			String[] popEnt = pop_stat.get(pid);
			boolean res = timeRange[0] < Integer.parseInt(popEnt[Abstract_Runnable_ClusterModel.POP_INDEX_EXIT_POP_AT]);
			if(timeRange.length > 1) {
				res &= Integer.parseInt(popEnt[Abstract_Runnable_ClusterModel.POP_INDEX_ENTER_POP_AT]) < timeRange[1]; 
			}					
			return res;
		}
	}

	public File[] getContactMapFiles() {
		return contactMap_files;
	}

	public String[] getContactMap_nextString() {
		return contactMap_nextString;
	}

	public void setContactMapFiles(File[] contactMapFiles) {
		this.contactMap_files = contactMapFiles;
		this.contactMap_nextString = new String[contactMapFiles.length];
		this.contactMap_reader = new BufferedReader[contactMapFiles.length];
	}

	public void setBaseProp(Properties sim_prop) {
		this.baseProp = sim_prop;
	}

	public abstract void allocateSeedInfection(int[][] num_infected_count, int time);

	public abstract int addInfectious(Integer infectedId, int infId, int site, int stage_id, int infectious_time,
			int recoveredAt);

	public abstract void scheduleNextTest(Integer personId, int lastTestTime);

	public abstract void refreshField(int fieldId, boolean clearAll);

	public long getSim_seed() {
		return sIM_SEED;
	}

	public long getcMap_seed() {
		return cMAP_SEED;
	}

	public void setPrint_progress(PrintStream print_progess) {
		this.print_progress = print_progess;
	}

	public int getSimSetting() {
		return simSetting;
	}

	public void setSimSetting(int simSetting) {
		this.simSetting = simSetting;
	}

	public void setEdges_list(ArrayList<Integer[]> edges_list) {
		this.edges_list = edges_list;
	}

	public void setPropSwitch_map(HashMap<Integer, HashMap<Integer, String>> propSwitch_map) {
		this.propSwitch_map = propSwitch_map;
	}

	protected void loadNonRunnableFieldSetting(Integer index, String entry, int loadTime) {
		
		switch (index.intValue()) {
		case Population_Bridging.FIELD_PARTNER_TYPE_PROB:
			// Reset
			non_mapped_encounter_prob = null;
			non_mapped_encounter_target_gender = null;
			break;
		case SIM_OFFSET + Simulation_ClusterModelTransmission.SIM_FIELD_SEED_INFECTION:
			int[][] num_infected = (int[][]) PropValUtils.propStrToObject(entry, int[][].class);
			allocateSeedInfection(num_infected, loadTime);
			break;
		case PROPSWITCH_MAP_KEPT_EDGE_BY_DURATION_INDEX:
			cMap_Kept_Edge_By_Duration_Setting = entry.trim().length() > 0
					? (float[][]) PropValUtils.propStrToObject(entry, float[][].class)
					: null;
			break;
		default:
			System.err.printf("Warning: Index of %d in simSpecificSwitch.prop not supported.\n", index);
		}
	}

	public int getGenderType(Integer personId) {
		String[] pT = null;
		if (pop_stat != null) {
			pT = pop_stat.get(personId);
		}
		if (pT == null) {
			return getGenderType(personId, cUMULATIVE_POP_COMPOSITION);
		} else {
			return Integer.parseInt(pT[Abstract_Runnable_ClusterModel.POP_INDEX_GRP]);
		}
	}

	public int exitPopAt(Integer personId) {
		int res = Integer.MAX_VALUE;
		if (pop_stat != null) {
			String[] pT = pop_stat.get(personId);
			if (pT != null) {
				res = Integer.parseInt(pT[Abstract_Runnable_ClusterModel.POP_INDEX_EXIT_POP_AT]);
			}
		}
		return res;
	}

	public void fillRiskCatMap(ArrayList<Number[]> prealloactedRiskGrpArr) {
		if (prealloactedRiskGrpArr != null) {
			for (Number[] preAllocRisk : prealloactedRiskGrpArr) {
				if (risk_cat_map == null) {
					risk_cat_map = new HashMap<>();
				}

				risk_cat_map.put(
						(Integer) preAllocRisk[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_PID],
						(Integer) preAllocRisk[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP]);
			}
		}
	}

	protected Integer[][] getEdgesArrayFromBaseConctactMap() {
		if (edges_list == null) {
			if (bASE_CONTACT_MAP != null) {
				try {
					edges_list = generateMapEdgeArray(bASE_CONTACT_MAP).call();
				} catch (Exception e) {
					e.printStackTrace(System.err);
					System.err.println("Error in generating edge list from BASE_CONTACT_MAP. Exiting...");
					edges_list = new ArrayList<>();
					System.exit(-1);
				}
			} else {
				edges_list = new ArrayList<>();
			}
		}
		Integer[][] edges_array = edges_list.toArray(new Integer[edges_list.size()][]);
		return edges_array;
	}

	protected int initaliseCMap(ContactMap cMap, Integer[][] edges_array, int edges_array_pt, int startTime,
			HashMap<Integer, ArrayList<Integer[]>> removeEdges) {
		ArrayList<Integer[]> toRemove;

		// Skip invalid edges
		while (edges_array_pt < edges_array.length
				&& edges_array[edges_array_pt][Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] < startTime) {
			Integer[] edge = edges_array[edges_array_pt];
			int edge_start_time = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
			int expireAt = edge_start_time + edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];

			if (expireAt > startTime) {
				toRemove = removeEdges.get(expireAt);
				if (toRemove == null) {
					toRemove = new ArrayList<>();
					removeEdges.put(expireAt, toRemove);
				}
				toRemove.add(edge);

				for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
						Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
					if (!cMap.containsVertex(edge[index])) {
						cMap.addVertex(edge[index]);
					}
				}
				addPartnership(cMap, edge);
			}

			edges_array_pt++;
		}

		// Add multimap edges if needed
		if (bASE_CONTACT_MAP == null && contactMap_files != null) {
			for (int i = 0; i < contactMap_files.length; i++) {
				updateCMapFromFiles(cMap, startTime, i, removeEdges, null);
			}
		}

		return edges_array_pt;
	}

	protected void setPreAllocatedRiskFromFile() {
		File pre_allocate_risk_file = new File(baseDir,
				String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, cMAP_SEED));

		if (pre_allocate_risk_file.isFile()) {
			try {
				BufferedReader reader = new BufferedReader(new FileReader(pre_allocate_risk_file));
				String line;

				if (risk_cat_map == null) {
					risk_cat_map = new HashMap<>();

				}

				while ((line = reader.readLine()) != null) {
					String[] lineSp = line.split(",");
					risk_cat_map.put(
							Integer.parseInt(
									lineSp[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_PID]),
							Integer.parseInt(
									lineSp[Simulation_ClusterModelTransmission.PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP]));
				}

				reader.close();
			} catch (Exception e) {
				e.printStackTrace(System.err);
			}
		}
	}

	private int updateCMapFromFiles(ContactMap cMap, int currentTime, int mapNumber,
			HashMap<Integer, ArrayList<Integer[]>> removeEdges, ArrayList<Integer> included_pids) {
		String[] edgeSp;
		Integer[] add_e;
		ArrayList<Integer[]> addedEdges = new ArrayList<>();
		int numEdgeAdded = 0;				

		if (contactMap_nextString[mapNumber] != null) {
			edgeSp = contactMap_nextString[mapNumber].split(",");
			add_e = new Integer[edgeSp.length];
			for (int i = 0; i < add_e.length; i++) {
				add_e[i] = Integer.parseInt(edgeSp[i]);
			}
			if (add_e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] <= currentTime
					&& (add_e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
							+ add_e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]) > currentTime) {
				if (addSelectiveEdge(included_pids, add_e)) {
					addPartnership(cMap, add_e);
					addedEdges.add(add_e);
					numEdgeAdded++;					
				}
				contactMap_nextString[mapNumber] = null;
			}						
		}

		// Reading next set of edges
		if (contactMap_nextString[mapNumber] == null) {
			try {
				if (contactMap_reader[mapNumber] == null) {
					contactMap_reader[mapNumber] = new BufferedReader(new FileReader(contactMap_files[mapNumber]));
				}
				while ((contactMap_nextString[mapNumber] = contactMap_reader[mapNumber].readLine()) != null) {
					edgeSp = contactMap_nextString[mapNumber].split(",");
					add_e = new Integer[edgeSp.length];
					for (int i = 0; i < add_e.length; i++) {
						add_e[i] = Integer.parseInt(edgeSp[i]);
					}
					if (add_e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] <= currentTime) {
						if ((add_e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
								+ add_e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION]) > currentTime) {
							addPartnership(cMap, add_e);
							addedEdges.add(add_e);
							numEdgeAdded++;
						}
					} else {
						break;
					}
				}

			} catch (IOException e) {
				e.printStackTrace(System.err);
			}
		}
		// Update remove edges
		ArrayList<Integer[]> toRemove;
		for (Integer[] edge : addedEdges) {
			int edge_start_time = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
			int expireAt = edge_start_time + edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];
			if (expireAt > currentTime) {
				toRemove = removeEdges.get(expireAt);
				if (toRemove == null) {
					toRemove = new ArrayList<>();
					removeEdges.put(expireAt, toRemove);
				}
				toRemove.add(edge);
			}

		}

		return numEdgeAdded;

	}

	protected int[] updateCMap(ContactMap cMap, int currentTime, Integer[][] edges_array, int edges_array_pt,
			HashMap<Integer, ArrayList<Integer[]>> edgesToRemove) {
		return updateCMap(cMap, currentTime, edges_array, edges_array_pt, edgesToRemove, null);
	}

	// return int[]{edges_array_pt, numEdgeAdded}
	protected int[] updateCMap(ContactMap cMap, int currentTime, Integer[][] edges_array, int edges_array_pt,
			HashMap<Integer, ArrayList<Integer[]>> edgesToRemove, ArrayList<Integer> included_pids) {

		int numEdgeAdded = 0;

		// Add multimap and read edges if needed
		if (bASE_CONTACT_MAP == null && contactMap_files != null) {
			for (int i = 0; i < contactMap_files.length; i++) {
				numEdgeAdded += updateCMapFromFiles(cMap, currentTime, i, edgesToRemove, included_pids);
			}
		}

		ArrayList<Integer[]> toRemove;
		// Remove expired edges
		toRemove = edgesToRemove.get(currentTime);
		if (toRemove != null) {
			for (Integer[] edge : toRemove) {
				cMap.removeEdge(edge);
			}
		}

		// Add new edges and update removal schedule
		while (edges_array_pt < edges_array.length
				&& edges_array[edges_array_pt][Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] <= currentTime) {
			boolean addEdge = cMap_Kept_Edge_By_Duration_Setting == null;
			Integer[] edge = edges_array[edges_array_pt];

			if (cMap_Kept_Edge_By_Duration_Setting != null) {
				int index = Arrays.binarySearch(cMap_Kept_Edge_By_Duration_Setting[cMAP_KEPT_EDGE_BY_DURATION_DUR_LIST],
						Math.max(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION], 1));
				if (index < 0) {
					index = ~index;
				}
				addEdge = index >= cMap_Kept_Edge_By_Duration_Setting[cMAP_KEPT_EDGE_BY_DURATION_PROB_LIST].length;
				if (!addEdge) {
					float keptEdgeProb = cMap_Kept_Edge_By_Duration_Setting[cMAP_KEPT_EDGE_BY_DURATION_PROB_LIST][index];
					if (keptEdgeProb > 0) {
						addEdge = keptEdgeProb >= 1 ? true : RNG.nextFloat() < keptEdgeProb;
					}
				}
			}

			// Check for one-off not involving infectious
			if (addEdge) {
				addEdge = addSelectiveEdge(included_pids, edge);
			}
			if (addEdge) {
				Integer expireAt = edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]
						+ Math.max(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION], 1);

				toRemove = edgesToRemove.get(expireAt);
				if (toRemove == null) {
					toRemove = new ArrayList<>();
					edgesToRemove.put(expireAt, toRemove);
				}
				toRemove.add(edge);

				for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
						Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
					if (!cMap.containsVertex(edge[index])) {
						cMap.addVertex(edge[index]);
					}
				}

				addPartnership(cMap, edge);
				numEdgeAdded++;
			}

			edges_array_pt++;
		}

		// Initialise non-mapped encounter
		if (non_mapped_encounter_prob == null) {
			non_mapped_encounter_prob = new float[cUMULATIVE_POP_COMPOSITION.length];
			non_mapped_encounter_target_gender = new int[cUMULATIVE_POP_COMPOSITION.length];
			boolean hasProb = false;

			// Only valid if it is defined in the prop file
			String prop_ent_str = getSim_prop()
					.getProperty(String.format("POP_PROP_INIT_PREFIX_%d", Population_Bridging.FIELD_PARTNER_TYPE_PROB));

			if (prop_ent_str != null) {
				float[][] prop_ent = (float[][]) util.PropValUtils.propStrToObject(prop_ent_str, float[][].class);
				for (int g = 0; g < cUMULATIVE_POP_COMPOSITION.length; g++) {
					float[] ent = prop_ent[g];
					if (ent.length > Population_Bridging.PARTNER_TYPE_NON_MAPPED_ENCOUNTER_PROB) {
						non_mapped_encounter_prob[g] = ent[Population_Bridging.PARTNER_TYPE_NON_MAPPED_ENCOUNTER_PROB];
						non_mapped_encounter_target_gender[g] = (int) ent[Population_Bridging.PARTNER_TYPE_NON_MAPPED_ENCOUNTER_TARGET_GENDER];
						hasProb |= non_mapped_encounter_prob[g] > 0;
					}
				}
			}
			if (!hasProb) {
				non_mapped_encounter_prob = new float[0]; // 0-length = not used
			} else {
				// Direct file
				String fName = String.format(Simulation_ClusterModelTransmission.FILENAME_NON_MAPPED_ENCOUNTER,
						cMAP_SEED, sIM_SEED);
				File preGenFile = new File(baseDir, fName);

				if (preGenFile.exists()) {
					try {
						BufferedReader reader = new BufferedReader(new FileReader(preGenFile));
						non_map_edges_pre_gen = new ArrayList<>();
						String line;
						while ((line = reader.readLine()) != null) {
							non_map_edges_pre_gen.add(line.split(","));
						}
						reader.close();
					} catch (IOException e) {
						e.printStackTrace(System.err);
						non_map_edges_pre_gen = null;
					}
				}
				// Zip file
				String fNameZip = fName.replace(String.format("_%d.csv", sIM_SEED), ".csv.7z");
				preGenFile = new File(baseDir, fNameZip);
				if (preGenFile.exists()) {
					try {
						HashMap<String, ArrayList<String[]>> entMap = util.Util_7Z_CSV_Entry_Extract_Callable
								.extractedLinesFrom7Zip(preGenFile);
						non_map_edges_pre_gen = entMap.get(fName);
					} catch (IOException e) {
						e.printStackTrace(System.err);
						non_map_edges_pre_gen = null;
					}
				}
				if (non_map_edges_pre_gen == null) {
					init_non_mapped_encounter_list();
				}
			}
		}

		if (non_mapped_encounter_prob.length > 0) {
			ArrayList<Integer[]> non_map_edge_to_add_array;
			if (non_map_edges_pre_gen != null) {
				non_map_edge_to_add_array = new ArrayList<>();
				while (non_map_edges_pre_gen_last_edge_pt < non_map_edges_pre_gen.size()) {
					String[] edge_s = non_map_edges_pre_gen.get(non_map_edges_pre_gen_last_edge_pt);
					int non_map_edges_pre_gen_last_edge_time = Integer
							.parseInt(edge_s[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);

					if (non_map_edges_pre_gen_last_edge_time > currentTime) {
						break;
					}
					// Check if the last edge match with current time
					if (non_map_edges_pre_gen_last_edge_time == currentTime) {
						Integer[] edge = new Integer[edge_s.length];
						for (int s = 0; s < edge_s.length; s++) {
							edge[s] = Integer.valueOf(edge_s[s]);
						}
						non_map_edge_to_add_array.add(edge);
					}
					non_map_edges_pre_gen_last_edge_pt++;
				}
			} else {
				non_map_edge_to_add_array = form_non_mapped_edges(cMap, currentTime);
			}

			// Adding of non-map edge
			for (Integer[] nm_edge : non_map_edge_to_add_array) {
				toRemove = edgesToRemove.get(currentTime + 1);
				if (toRemove == null) {
					toRemove = new ArrayList<>();
					edgesToRemove.put(currentTime + 1, toRemove);
				}
				toRemove.add(nm_edge);

				for (int index : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
						Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
					if (!cMap.containsVertex(nm_edge[index])) {
						cMap.addVertex(nm_edge[index]);
					}
				}
				addPartnership(cMap, nm_edge);
				addNonMapEdgeStore(nm_edge);
				numEdgeAdded++;
			}

		}

		return new int[] { edges_array_pt, numEdgeAdded };
	}

	private boolean addSelectiveEdge(ArrayList<Integer> included_pids, Integer[] edge) {
		boolean addEdge = included_pids == null || edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] > 1;
		if (!addEdge) {
			final int[] partnerIndices = new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
					Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 };
			for (int index : partnerIndices) {
				addEdge |= (Collections.binarySearch(included_pids, index) >= 0);
				if (addEdge) {
					return addEdge;
				}
			}
		}
		return addEdge;
	}

	protected void addPartnership(ContactMap cMap, Integer[] edge) {
		for (int i : new int[] { Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1,
				Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2 }) {
			if (!cMap.containsVertex(edge[i])) {
				cMap.addVertex(edge[i]);
			}
		}
		cMap.addEdge(edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P1],
				edge[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_P2], edge);
	}

	public HashMap<String, Object> getSim_output() {
		return sim_output;
	}

	public abstract ArrayList<Integer> loadOptParameter(String[] parameter_settings, double[] point,
			int[][] seedInfectNum, boolean display_only);

	protected ArrayList<Integer[]> form_non_mapped_edges(ContactMap cMap, int currentTime) {
		ArrayList<Integer[]> non_map_edge_to_add_array = new ArrayList<Integer[]>();
		int[] num_non_map_seeker = new int[non_mapped_encounter_prob.length];
		for (int g = 0; g < num_non_map_seeker.length; g++) {
			if (dist_nm_candidate[g] != null) {
				num_non_map_seeker[g] = dist_nm_candidate[g].sample();
			}
		}
		for (int seek_g = 0; seek_g < num_non_map_seeker.length; seek_g++) {
			if (num_non_map_seeker[seek_g] > 0) {
				int num_to_seek = Math.min(num_non_map_seeker[seek_g], Math
						.min(non_map_candidate_id_seeker[seek_g].size(), non_map_candidate_id_target[seek_g].size()));
				if (num_to_seek > 0) {
					int[] nm_seeker_index = ArrayUtilsRandomGenerator.randomSelectIndex(num_to_seek, 0,
							non_map_candidate_id_seeker[seek_g].size(), RNG_NM);
					int[] nm_target_index = ArrayUtilsRandomGenerator.randomSelectIndex(num_to_seek, 0,
							non_map_candidate_id_seeker[seek_g].size(), RNG_NM);
					if (num_to_seek > 1) {
						ArrayUtilsRandomGenerator.shuffleArray(nm_seeker_index, RNG_NM);
						ArrayUtilsRandomGenerator.shuffleArray(nm_target_index, RNG_NM);
					}
					for (int p = 0; p < num_to_seek; p++) {
						int seeker_id = non_map_candidate_id_seeker[seek_g].get(nm_seeker_index[p]);
						int target_id = non_map_candidate_id_target[seek_g].get(nm_target_index[p]);
						int p1 = Math.min(seeker_id, target_id);
						int p2 = p1 == seeker_id ? target_id : seeker_id;
						if (seeker_id != target_id && !cMap.containsEdge(p1, p2)) {
							non_map_edge_to_add_array.add(new Integer[] { p1, p2, currentTime, 1 });
						}
					}
				}
			}
		}
		return non_map_edge_to_add_array;
	}

	@SuppressWarnings("unchecked")
	protected void init_non_mapped_encounter_list() {

		// Seeking non-mapped casual partnership
		non_map_candidate_id_target = new ArrayList[cUMULATIVE_POP_COMPOSITION.length];
		non_map_candidate_id_seeker = new ArrayList[cUMULATIVE_POP_COMPOSITION.length];
		dist_nm_candidate = new IntegerDistribution[cUMULATIVE_POP_COMPOSITION.length];

		for (int seeker_g = 0; seeker_g < cUMULATIVE_POP_COMPOSITION.length; seeker_g++) {
			if (non_mapped_encounter_prob[seeker_g] > 0) {
				non_map_candidate_id_seeker[seeker_g] = new ArrayList<Integer>();
				non_map_candidate_id_target[seeker_g] = new ArrayList<Integer>();
				int id = 1;
				while (id <= cUMULATIVE_POP_COMPOSITION[cUMULATIVE_POP_COMPOSITION.length - 1]) {
					int gender_type = getGenderType(id);
					if (gender_type == seeker_g) {
						non_map_candidate_id_seeker[seeker_g].add(id);
					}
					if (((1 << gender_type) & non_mapped_encounter_target_gender[seeker_g]) > 0) {
						non_map_candidate_id_target[seeker_g].add(id);
					}
					id++;
				}
				dist_nm_candidate[seeker_g] = new PoissonDistribution(RNG_NM,
						non_mapped_encounter_prob[seeker_g] * non_map_candidate_id_seeker[seeker_g].size(),
						PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
			}
		}

		if ((getSimSetting()
				& 1 << Simulation_ClusterModelTransmission.SIM_SETTING_KEY_GEN_NON_MAPPED_ENCOUNTER) != 0) {
			non_map_edges_store = new Integer[NON_MAP_EDGES_STORE_MAX][];
		}

	}

	protected void exportNonMapEdgeStore() {
		if (non_map_edges_store != null) {
			String fName = String.format(Simulation_ClusterModelTransmission.FILENAME_NON_MAPPED_ENCOUNTER, cMAP_SEED,
					sIM_SEED);
			PrintWriter pWri;
			try {
				pWri = new PrintWriter(new FileWriter(new File(baseDir, fName), true));
				for (int i = 0; i < NON_MAP_EDGES_STORE_PT; i++) {
					Integer[] edges = non_map_edges_store[i];
					boolean first_ent = true;
					for (Integer edge : edges) {
						if (!first_ent) {
							pWri.print(',');
						} else {
							first_ent = false;
						}
						pWri.print(edge);
					}
					pWri.println();
				}
				pWri.close();
				NON_MAP_EDGES_STORE_PT = 0;
			} catch (IOException e) {
				e.printStackTrace(System.err);
			}
		}
	}

	private void addNonMapEdgeStore(Integer[] nm_edge) {
		if (non_map_edges_store != null) {
			non_map_edges_store[NON_MAP_EDGES_STORE_PT] = nm_edge;
			NON_MAP_EDGES_STORE_PT++;
			if (NON_MAP_EDGES_STORE_PT >= non_map_edges_store.length) {
				exportNonMapEdgeStore();

			}
		}
	}

	protected void postSimulation() {
		if (contactMap_reader != null) {
			try {
				for (BufferedReader r : contactMap_reader) {
					if (r != null) {
						r.close();
					}
				}
			} catch (IOException e) {
				e.printStackTrace(System.err);
			}
		}
		exportNonMapEdgeStore();
	}

}
