package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.archivers.sevenz.SevenZArchiveEntry;
import org.apache.commons.compress.archivers.sevenz.SevenZFile;
import org.apache.commons.compress.archivers.sevenz.SevenZOutputFile;
import org.apache.commons.io.FileUtils;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import relationship.TransmissionMap;
import util.PersonClassifier;
import util.PropValUtils;

public class Simulation_ClusterModelTransmission implements SimulationInterface {

	protected Properties loadedProperties = null; // From .prop file, if any
	protected boolean stopNextTurn = false;
	protected File baseDir = null;
	protected ContactMap baseContactMap;
	protected long baseContactMapSeed = -1;
	protected ArrayList<Number[]> prealloactedRiskGrpArr = null;
	protected int simSetting = 1 << SIM_SETTING_KEY_TRACK_TRANSMISSION_CLUSTER;
	protected boolean exportSkipBackup = false;
	protected boolean printProgress = false;

	public static final String POP_PROP_INIT_PREFIX = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX;
	public static final String POP_PROP_INIT_PREFIX_CLASS = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX_CLASS;
	public static final String FILENAME_FORMAT_ALL_CMAP = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP;

	public static final String LAUNCH_ARGS_SKIP_BACKUP = "-export_skip_backup";
	public static final String LAUNCH_ARGS_PRINT_PROGRESS = "-printProgress";

	public static final String FILENAME_INDEX_CASE_LIST = "Seed_IndexCases_%d_%d.txt";
	public static final String FILENAME_PREVALENCE_SITE = "Prevalence_Site_%d_%d.csv";
	public static final String FILENAME_CUMUL_INCIDENCE_SITE = "Incidence_Site_%d_%d.csv";
	public static final String FILENAME_PREVALENCE_PERSON = "Prevalence_Person_%d_%d.csv";
	public static final String FILENAME_CUMUL_INCIDENCE_PERSON = "Incidence_Person_%d_%d.csv";
	public static final String FILENAME_CUMUL_TREATMENT_PERSON = "Treatment_Person_%d_%d.csv";
	public static final String FILENAME_CUMUL_POSITIVE_DX_PERSON = "Positive_DX_Person_%d_%d.csv";
	public static final String FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON = "Positive_DX_Sought_Person_%d_%d.csv";
	public static final String FILENAME_INFECTION_HISTORY = "InfectHist_%d_%d.csv";
	public static final String FILENAME_CUMUL_ANTIBIOTIC_USAGE = "Antibiotic_usage_%d_%d.csv";
	public static final String FILENAME_ALL_TRANSMISSION_CMAP = "All_transmap_%d_%d.csv";
	public static final String FILENAME_VACCINE_COVERAGE = "Vaccine_coverage_%d_%d.csv";
	public static final String FILENAME_VACCINE_COVERAGE_PERSON = "Vaccine_coverage_Person_%d_%d.csv";

	public static final String FILENAME_INDEX_CASE_LIST_ZIP = FILENAME_INDEX_CASE_LIST.replaceFirst("_%d", "") + ".7z";
	public static final String FILENAME_PREVALENCE_SITE_ZIP = FILENAME_PREVALENCE_SITE.replaceFirst("_%d", "") + ".7z";
	public static final String FILENAME_PREVALENCE_PERSON_ZIP = FILENAME_PREVALENCE_PERSON.replaceFirst("_%d", "")
			+ ".7z";
	public static final String FILENAME_CUMUL_INCIDENCE_SITE_ZIP = FILENAME_CUMUL_INCIDENCE_SITE.replaceFirst("_%d", "")
			+ ".7z";
	public static final String FILENAME_CUMUL_INCIDENCE_PERSON_ZIP = FILENAME_CUMUL_INCIDENCE_PERSON.replaceFirst("_%d",
			"") + ".7z";
	public static final String FILENAME_CUMUL_POSITIVE_DX_PERSON_ZIP = FILENAME_CUMUL_POSITIVE_DX_PERSON
			.replaceFirst("_%d", "") + ".7z";
	public static final String FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON_ZIP = FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON
			.replaceFirst("_%d", "") + ".7z";

	public static final String FILENAME_CUMUL_TREATMENT_PERSON_ZIP = FILENAME_CUMUL_TREATMENT_PERSON.replaceFirst("_%d",
			"") + ".7z";
	public static final String FILENAME_INFECTION_HISTORY_ZIP = FILENAME_INFECTION_HISTORY.replaceFirst("_%d", "")
			+ ".7z";
	public static final String FILENAME_CUMUL_ANTIBIOTIC_USAGE_ZIP = FILENAME_CUMUL_ANTIBIOTIC_USAGE.replaceFirst("_%d",
			"") + ".7z";
	public static final String FILENAME_ALL_TRANSMISSION_CMAP_ZIP = FILENAME_ALL_TRANSMISSION_CMAP.replaceFirst("_%d",
			"") + ".7z";
	public static final String FILENAME_VACCINE_COVERAGE_ZIP = FILENAME_VACCINE_COVERAGE.replaceFirst("_%d", "")
			+ ".7z";
	public static final String FILENAME_VACCINE_COVERAGE_PERSON_ZIP = FILENAME_VACCINE_COVERAGE_PERSON
			.replaceFirst("_%d", "") + ".7z";

	public static final String FILENAME_PRE_ALLOCATE_RISK_GRP = "RiskGrp_Map_%d.csv";
	public static final int PRE_ALLOCATE_RISK_GRP_INDEX_PID = 0;
	public static final int PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP = PRE_ALLOCATE_RISK_GRP_INDEX_PID + 1;
	public static final int PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATE_TOTAL = PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP + 1;
	public static final int PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATR_FIRSTSNAP = PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATE_TOTAL
			+ 1;

	// Switching parameter
	public static final String FILENAME_PROP_SWITCH = "simSpecificSwitch.prop";
	public static final String POP_PROP_SWITCH_PREFIX = "SWITCH_%d_";
	public static final String POP_PROP_SWITCH_AT = "PROP_SIM_SWITCH_AT";
	protected HashMap<Integer, HashMap<Integer, String>> propSwitch_map = new HashMap<>();

	// Sim setting to indicates what type of simulation need to be run.
	// Sim setting is active if simSetting & 1 << SIM_SETTING_KEY != 0
	public static final String PROP_SIM_SETTING = "PROP_SIM_SETTING";

	public static final int SIM_SETTING_KEY_GLOBAL_TIME_SEED = 0;
	public static final int SIM_SETTING_KEY_GEN_PREVAL_FILE = SIM_SETTING_KEY_GLOBAL_TIME_SEED + 1;
	public static final int SIM_SETTING_KEY_GEN_INCIDENCE_FILE = SIM_SETTING_KEY_GEN_PREVAL_FILE + 1;
	public static final int SIM_SETTING_KEY_GEN_TREATMENT_FILE = SIM_SETTING_KEY_GEN_INCIDENCE_FILE + 1;
	public static final int SIM_SETTING_KEY_TRACK_TRANSMISSION_CLUSTER = SIM_SETTING_KEY_GEN_TREATMENT_FILE + 1;
	public static final int SIM_SETTING_KEY_TRACK_INFECTION_HISTORY = SIM_SETTING_KEY_TRACK_TRANSMISSION_CLUSTER + 1;
	public static final int SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE = SIM_SETTING_KEY_TRACK_INFECTION_HISTORY + 1;
	public static final int SIM_SETTING_KEY_TRACK_VACCINE_COVERAGE = SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE + 1;
	public static final int SIM_SETTING_KEY_TREATMENT_ON_INFECTIOUS_ONLY = SIM_SETTING_KEY_TRACK_VACCINE_COVERAGE + 1;

	public static final String PROP_CONTACT_MAP_LOC = "PROP_CONTACT_MAP_LOC";

	public static final Object[] DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS = {
			// BRIDGING_MAP_TRANS_SIM_FIELD_SEED_INFECTION
			// Usage: [gender_type][site] = probability (out of total)
			new float[][] {},

	};

	public static final int SIM_FIELD_SEED_INFECTION = 0;
	public static final int LENGTH_SIM_MAP_TRANSMISSION_FIELD = SIM_FIELD_SEED_INFECTION + 1;

	public Object[] simFields = new Object[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
			+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
			+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
			+ LENGTH_SIM_MAP_TRANSMISSION_FIELD
			+ Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD];
	public Class<?>[] simFieldClass = new Class[simFields.length];

	public Simulation_ClusterModelTransmission() {
		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
				+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD;
		for (int i = 0; i < simFields.length; i++) {
			// All simulation levels
			if (i >= sim_offset && (i - sim_offset) < DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS.length) {
				simFields[i] = DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS[i - sim_offset];
				simFieldClass[i] = DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS[i - sim_offset].getClass();
			}
		}
	}

	public int getSimSetting() {
		return simSetting;
	}

	public void setBaseContactMapSeed(long baseContactMapSeed) {
		this.baseContactMapSeed = baseContactMapSeed;
		try {
			loadPreallocateRiskGrp(baseContactMapSeed);
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
		}
	}

	public void setExportSkipBackup(boolean exportSkipBackup) {
		this.exportSkipBackup = exportSkipBackup;
	}

	public void setPrintProgress(boolean printProgress) {
		this.printProgress = printProgress;
	}

	@Override
	public void loadProperties(Properties prop) {
		loadedProperties = prop;

		for (int i = 0; i < simFields.length; i++) {
			String propName = String.format("%s%d", POP_PROP_INIT_PREFIX, i);
			if (prop.containsKey(propName)) {
				String objStr = prop.getProperty(propName);
				if (simFieldClass[i] == null) {
					String str = String.format("%s%d", POP_PROP_INIT_PREFIX_CLASS, i);
					if (prop.containsKey(str)) {
						str = prop.getProperty(str);
						try {
							simFieldClass[i] = Class.forName(str);
						} catch (ClassNotFoundException e) {
							e.printStackTrace();
							simFieldClass[i] = Object.class;
						}
					}
				}
				simFields[i] = PropValUtils.propStrToObject(objStr, simFieldClass[i]);
			}
		}

		if (prop.containsKey(PROP_SIM_SETTING)) {
			simSetting = Integer.parseInt(prop.getProperty(PROP_SIM_SETTING));
		}

		File propSwitchFile = new File(baseDir, FILENAME_PROP_SWITCH);

		try {
			if (propSwitchFile.exists()) {
				Properties prop_switch = new Properties();
				FileInputStream fIS_switch = new FileInputStream(propSwitchFile);
				prop_switch.loadFromXML(fIS_switch);
				fIS_switch.close();
				int[] switchTime = (int[]) PropValUtils.propStrToObject(prop_switch.getProperty(POP_PROP_SWITCH_AT),
						int[].class);

				for (int sI = 0; sI < switchTime.length; sI++) {
					int sTime = switchTime[sI];
					String switch_prefix = String.format(POP_PROP_SWITCH_PREFIX, sI);

					for (Object key : prop_switch.keySet()) {
						String keyStr = key.toString();
						if (keyStr.startsWith(switch_prefix)) {
							Integer runnableKey = Integer
									.parseInt(keyStr.substring(switch_prefix.length() + POP_PROP_INIT_PREFIX.length()));
							HashMap<Integer, String> ent = propSwitch_map.get(sTime);
							if (ent == null) {
								ent = new HashMap<>();
								propSwitch_map.put(sTime, ent);
							}
							ent.put(runnableKey, prop_switch.getProperty(keyStr));

						}
					}
				}

				System.out.printf("Properties switch file < %s > loaded.\n", propSwitchFile.getAbsolutePath());
			}
		} catch (IOException ex) {
			ex.printStackTrace(System.err);

		}

	}

	@Override
	public void setStopNextTurn(boolean stopNextTurn) {
		this.stopNextTurn = stopNextTurn;

	}

	@Override
	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
	}

	@Override
	public Properties generateProperties() {
		Properties prop;
		if (loadedProperties == null) {
			prop = new Properties();
		} else {
			prop = loadedProperties;
		}

		prop.put("PROP_GENERATED_AT", Long.toString(System.currentTimeMillis()));
		return prop;
	}

	@Override
	public void setSnapshotSetting(PersonClassifier[] snapshotCountClassifier, boolean[] snapshotCountAccum) {
		throw new UnsupportedOperationException("Method setSnapshotSetting not support yet.");

	}

	public void setBaseContactMap(ContactMap baseContactMap) {
		this.baseContactMap = baseContactMap;
	}

	@Override
	public void generateOneResultSet() throws IOException, InterruptedException {

		int numThreads = Runtime.getRuntime().availableProcessors();
		int numSim = 1;
		long seed = System.currentTimeMillis();
		int num_snap = 1;
		int num_time_steps_per_snap = 1;
		int[] pop_composition = new int[] { 500000, 500000, 20000, 20000 };

		if (loadedProperties != null) {
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL])) {
				numThreads = Integer.parseInt(loadedProperties
						.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_USE_PARALLEL]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET])) {
				numSim = Integer.parseInt(loadedProperties
						.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SIM_PER_SET]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
				seed = Long.parseLong(
						loadedProperties.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP])) {
				num_snap = Integer.parseInt(
						loadedProperties.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
				num_time_steps_per_snap = Integer.parseInt(loadedProperties
						.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ]));
			}

			String popCompositionKey = POP_PROP_INIT_PREFIX
					+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);
			if (loadedProperties.containsKey(popCompositionKey)) {
				pop_composition = (int[]) PropValUtils.propStrToObject(loadedProperties.getProperty(popCompositionKey),
						int[].class);
			}
		}

		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
				+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD;

		final int ent_offset = sim_offset + LENGTH_SIM_MAP_TRANSMISSION_FIELD;

		// Check error in zipped output
		File[] zipFiles = baseDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.getName().endsWith(".csv.7z") || pathname.getName().endsWith(".txt.7z");
			}
		});

		for (File zf : zipFiles) {
			try {
				SevenZFile archive7Z = new SevenZFile(zf);
				archive7Z.close();

			} catch (IOException ex) {
				System.out.printf("Error in reading '%s' --> ", zf.getCanonicalPath());

				Pattern replace_file_pattern = Pattern.compile(zf.getName() + "_(\\d+).7z");

				File[] possible_replace_file = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						Matcher m = replace_file_pattern.matcher(pathname.getName());
						return m.matches();
					}
				});

				boolean replacement_found = false;

				Arrays.sort(possible_replace_file, new Comparator<File>() {
					@Override
					public int compare(File o1, File o2) {
						Matcher m1 = replace_file_pattern.matcher(o1.getName());
						m1.matches();
						Matcher m2 = replace_file_pattern.matcher(o2.getName());
						m2.matches();
						long v1 = Long.parseLong(m1.group(1));
						long v2 = Long.parseLong(m2.group(1));
						return -Long.compare(v1, v2);
					}
				});

				for (int i = 0; i < possible_replace_file.length && !replacement_found; i++) {
					try {
						SevenZFile archive7Z = new SevenZFile(possible_replace_file[i]);
						archive7Z.close();
						FileUtils.copyFile(possible_replace_file[i], zf, StandardCopyOption.REPLACE_EXISTING);
						System.out.printf("Replaced by '%s'\n", possible_replace_file[i].getCanonicalPath());
						replacement_found = true;
					} catch (IOException ex1) {
						// Skip to next available file
					}
				}

				if (!replacement_found) {
					System.out.println("Replacement not found - file removed.");
					FileUtils.delete(zf);
				}
			}
		}
		// Zipping CSV if found
		File[] preCSV = baseDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.getName().endsWith(".csv");
			}
		});
		if (preCSV.length > 0) {
			zipOutputFiles();
		}

		// Map stat
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] personStat = new ArrayList[Population_Bridging.LENGTH_GENDER];
		int[] cumul_pop = new int[Population_Bridging.LENGTH_GENDER];

		int offset = 0;

		for (int g = 0; g < cumul_pop.length; g++) {
			cumul_pop[g] = offset + pop_composition[g];
			offset += pop_composition[g];
			personStat[g] = new ArrayList<Integer>();
		}

		for (Integer v : baseContactMap.vertexSet()) {
			int g = Runnable_ClusterModel_Transmission.getGenderType(v, cumul_pop);
			personStat[g].add(v);
		}

		int[] personCount = new int[personStat.length];
		for (int i = 0; i < personCount.length; i++) {
			personCount[i] = personStat[i].size();
		}
		System.out.println(String.format("Number of indivduals in contact map = %s", Arrays.toString(personCount)));

		float[][] riskCatListAll = ((float[][]) simFields[ent_offset
				+ Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS]);

		int[] contactMapTimeRange = (int[]) simFields[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE];

		int[] cumulative_pop_composition = new int[pop_composition.length];
		int pop_offset = 0;

		for (int g = 0; g < cumulative_pop_composition.length; g++) {
			cumulative_pop_composition[g] = pop_offset + pop_composition[g];
			pop_offset += pop_composition[g];
		}

		if (prealloactedRiskGrpArr == null) {
			prealloactedRiskGrpArr = new ArrayList<>();
			fillRiskGrpArrByCasualPartnership(prealloactedRiskGrpArr, baseContactMap, cumulative_pop_composition,
					riskCatListAll, contactMapTimeRange);

			reallocateRiskGrp(baseContactMapSeed);
		}

		float[][] seedInfectParam = (float[][]) simFields[sim_offset + SIM_FIELD_SEED_INFECTION];
		int[][] seedInfectNum = new int[Population_Bridging.LENGTH_GENDER][Runnable_ClusterModel_Transmission.LENGTH_SITE];

		for (int g = 0; g < seedInfectParam.length; g++) {
			for (int s = 0; s < seedInfectParam[g].length; s++) {
				float seedProp = seedInfectParam[g][s];
				if (seedProp < 1) {
					seedInfectNum[g][s] = Math.round(personStat[g].size() * seedProp);
				} else {
					seedInfectNum[g][s] = Math.round(seedInfectParam[g][s]);
				}
			}
		}

		RandomGenerator rngBase = new MersenneTwisterRandomGenerator(seed);
		boolean useParallel = numThreads > 1 && numSim > 1;

		ExecutorService exec = null;
		int inExec = 0;

		Runnable_ClusterModel_Transmission[] runnable = new Runnable_ClusterModel_Transmission[numSim];
		File seedIndexExport7z = new File(baseDir,
				String.format(Simulation_ClusterModelTransmission.FILENAME_INDEX_CASE_LIST_ZIP, baseContactMapSeed));
		ArrayList<Long> completedSeed = new ArrayList<>();

		if (seedIndexExport7z.exists()) {
			SevenZFile inputZip = new SevenZFile(seedIndexExport7z);
			Pattern p_seedIndexList = Pattern.compile(FILENAME_INDEX_CASE_LIST
					.replaceFirst("%d", Long.toString(baseContactMapSeed)).replaceFirst("%d", "(-{0,1}(?!0)\\\\d+)"));
			SevenZArchiveEntry inputEnt;
			while ((inputEnt = inputZip.getNextEntry()) != null) {
				String entName = inputEnt.getName();
				Matcher m = p_seedIndexList.matcher(entName);
				if (m.matches()) {
					completedSeed.add(Long.parseLong(m.group(1)));
				}
			}
			inputZip.close();
		}

		Collections.sort(completedSeed);

		ArrayList<Integer[]> edge_list = null;
		try {
			edge_list = Abstract_Runnable_ClusterModel.generateMapEdgeArray(baseContactMap).call();
		} catch (Exception e1) {
			e1.printStackTrace(System.err);
		}

		for (int s = 0; s < numSim; s++) {
			long simSeed = rngBase.nextLong();
			boolean runSim = true;
			if (Collections.binarySearch(completedSeed, simSeed) >= 0) {
				System.out.println(String.format("Simulation under seed # %d already generated.", simSeed));
				runSim = false;
			}

			if (runSim) {
				runnable[s] = new Runnable_ClusterModel_Transmission_Map(baseContactMapSeed, simSeed, pop_composition,
						baseContactMap, num_time_steps_per_snap, num_snap);
				runnable[s].setBaseDir(baseDir);
				runnable[s].setEdges_list(edge_list);
				runnable[s].setSimSetting(simSetting);
				runnable[s].setPropSwitch_map(propSwitch_map);

				if (printProgress) {
					runnable[s].setPrint_progress(System.out);
					runnable[s].setRunnableId(String.format("%d,%d", baseContactMapSeed, simSeed));
				}

				for (int f = 0; f < Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; f++) {
					if (simFields[ent_offset + f] != null) {
						runnable[s].getRunnable_fields()[f] = simFields[ent_offset + f];
					}
				}

				if ((simSetting & 1 << SIM_SETTING_KEY_TRACK_TRANSMISSION_CLUSTER) != 0) {
					((Runnable_ClusterModel_Transmission_Map) runnable[s]).setTransmissionMap(new TransmissionMap());
				}

				if ((simSetting & 1 << SIM_SETTING_KEY_TREATMENT_ON_INFECTIOUS_ONLY) != 0) {
					System.out.println("Note: Assuming treatment applied on infectious only");
				}

				runnable[s].initialse();
				runnable[s].fillRiskCatMap(prealloactedRiskGrpArr);

				if ((simSetting & 1 << SIM_SETTING_KEY_GLOBAL_TIME_SEED) != 0) {
					runnable[s].allocateSeedInfection(seedInfectNum, contactMapTimeRange[0]);
				} else {
					// Add infected
					for (int gender = 0; gender < Population_Bridging.LENGTH_GENDER; gender++) {
						for (int site = 0; site < Runnable_ClusterModel_Transmission.LENGTH_SITE; site++) {
							if (seedInfectNum[gender][site] > 0) {
								Integer[] seedInf = util.ArrayUtilsRandomGenerator.randomSelect(
										personStat[gender].toArray(new Integer[personStat[gender].size()]),
										seedInfectNum[gender][site], rngBase);

								int[] seedTime = new int[seedInf.length];

								for (int i = 0; i < seedInf.length; i++) {
									Integer infected = seedInf[i];

									int firstContactTime = contactMapTimeRange[1]; // Start only from valid range

									Set<Integer[]> edgesOfInfected = baseContactMap.edgesOf(infected);

									for (Integer[] e : edgesOfInfected) {
										int edge_start_time = e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME];
										int edge_end_time = edge_start_time
												+ e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION];
										if (edge_end_time > contactMapTimeRange[0]
												&& edge_start_time < contactMapTimeRange[1]) {
											firstContactTime = Math.max(Math.min(firstContactTime, edge_start_time),
													contactMapTimeRange[0]);
										}
									}

									if (firstContactTime == contactMapTimeRange[1]) {
										firstContactTime = contactMapTimeRange[0];
									}

									seedTime[i] = firstContactTime;
									runnable[s].addInfectious(infected, site, firstContactTime, firstContactTime + 180);

								}
								System.out.println(String.format("Seeding %s of gender #%d at site #%d at t = %s",
										Arrays.toString(seedInf), gender, site, Arrays.toString(seedTime)));
							}
						}
					}
				}
				if (useParallel) {
					if (exec == null) {
						exec = Executors.newFixedThreadPool(numThreads);
					}
					exec.submit(runnable[s]);

					inExec++;
					if (inExec == numThreads) {
						exec.shutdown();
						if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
							System.err.println("Thread time-out!");
						}
						inExec = 0;
						exec = null;
						if (numSim > 1) {
							zipOutputFiles();
						}
					}
				} else {
					runnable[s].run();
					if (numSim > 1) {
						zipOutputFiles();
					}
				}
			}

		}
		if (exec != null && inExec != 0) {
			exec.shutdown();
			if (!exec.awaitTermination(2, TimeUnit.DAYS)) {
				System.err.println("Thread time-out!");
			}
			inExec = 0;
			exec = null;
			if (numSim > 1) {
				zipOutputFiles();
			}
		}
	}

	public static void fillRiskGrpArrByCasualPartnership(ArrayList<Number[]> riskGrpArr, ContactMap cMap,
			int[] cumulative_pop_composition, float[][] riskCatListAll, int[] timeStartFrom) {
		// Generated preallocated list
		for (float[] riskCatList : riskCatListAll) {
			if (riskCatList != null && riskCatList[0] < 0) {
				int genderIncl = -(int) riskCatList[0];
				for (int g = 0; g < cumulative_pop_composition.length; g++) {
					if ((genderIncl & 1 << g) > 0) {
						int g_start = 1;
						if (g > 1) {
							g_start = cumulative_pop_composition[g - 1] + 1;
						}
						for (int pid = g_start; pid <= cumulative_pop_composition[g]; pid++) {
							if (cMap.containsVertex(pid)) {
								int numCasual = 0;
								int numCasual1Yr = 0;
								int firstCasualPartnerTime = Integer.MAX_VALUE;
								int lastCasualPartnerTime = 0;
								int firstPartnerTime = Integer.MAX_VALUE;
								int lastPartnerTime = 0;
								Set<Integer[]> edges = cMap.edgesOf(pid);
								for (Integer[] e : edges) {
									firstPartnerTime = Math.min(firstPartnerTime,
											e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
									lastPartnerTime = Math.max(lastPartnerTime,
											e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);

									if (e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] >= timeStartFrom[0]
											&& e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_DURATION] <= 1) {
										firstCasualPartnerTime = Math.min(firstCasualPartnerTime,
												e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
										lastCasualPartnerTime = Math.max(lastCasualPartnerTime,
												e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME]);
										numCasual++;
										if (e[Abstract_Runnable_ClusterModel.CONTACT_MAP_EDGE_START_TIME] < timeStartFrom[0]
												+ AbstractIndividualInterface.ONE_YEAR_INT) {
											numCasual1Yr++;
										}
									}
								}

								float denom = lastCasualPartnerTime - firstCasualPartnerTime;
								if (denom <= 0) {									
									denom = lastPartnerTime - firstPartnerTime;
									if (denom <= 0) {
										denom = timeStartFrom[1] - timeStartFrom[0];
									}
								}

								float numCasualPerYear = (((float) AbstractIndividualInterface.ONE_YEAR_INT)
										* numCasual) / denom;
								float numCasualPerYear_1stYear = (((float) AbstractIndividualInterface.ONE_YEAR_INT)
										* numCasual1Yr) / denom;
								riskGrpArr.add(new Number[] { pid, -1, numCasualPerYear, numCasualPerYear_1stYear });

							}
						}
					}
				}
			}
		}
	}

	private void loadPreallocateRiskGrp(long baseContactMapSeed) throws NumberFormatException, IOException {		
		if(prealloactedRiskGrpArr == null) {
			prealloactedRiskGrpArr = new ArrayList<>();
		}
		boolean reallocate = loadPreallocateRiskGrp(prealloactedRiskGrpArr, baseDir, baseContactMapSeed);
		if (reallocate) {
			reallocateRiskGrp(baseContactMapSeed);
		}
	}

	public static boolean loadPreallocateRiskGrp(ArrayList<Number[]> prealloactedRiskGrpArr, File baseDir,
			long baseContactMapSeed) throws NumberFormatException, IOException {		
		File pre_allocate_risk_file = new File(baseDir,
				String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, baseContactMapSeed));
		boolean reallocate = false;

		if (pre_allocate_risk_file.exists()) {
			BufferedReader reader = new BufferedReader(new FileReader(pre_allocate_risk_file));
			String rLine;
			while ((rLine = reader.readLine()) != null) {
				String[] ent = rLine.split(",");
				Number[] map_ent = new Number[ent.length];
				map_ent[PRE_ALLOCATE_RISK_GRP_INDEX_PID] = Integer.parseInt(ent[PRE_ALLOCATE_RISK_GRP_INDEX_PID]);
				map_ent[PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP] = Integer
						.parseInt(ent[PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP]);
				map_ent[PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATE_TOTAL] = Float
						.parseFloat(ent[PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATE_TOTAL]);
				map_ent[PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATR_FIRSTSNAP] = Float
						.parseFloat(ent[PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATR_FIRSTSNAP]);

				reallocate |= map_ent[PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP].intValue() < 0;

				if (Float.isFinite(map_ent[PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATE_TOTAL].floatValue())) {
					prealloactedRiskGrpArr.add(map_ent);
				}

			}
			reader.close();
		}
		return reallocate;
	}

	private void reallocateRiskGrp(long baseContactMapSeed) {
		final int ent_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
				+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD + LENGTH_SIM_MAP_TRANSMISSION_FIELD;

		final String popCompositionKey = POP_PROP_INIT_PREFIX
				+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION);

		int[] pop_composition = (int[]) PropValUtils.propStrToObject(loadedProperties.getProperty(popCompositionKey),
				int[].class);
		int[] cumulative_pop_composition = new int[pop_composition.length];
		int pop_offset = 0;
		for (int g = 0; g < cumulative_pop_composition.length; g++) {
			cumulative_pop_composition[g] = pop_offset + pop_composition[g];
			pop_offset += pop_composition[g];
		}

		float[][] riskCatListAll = ((float[][]) simFields[ent_offset
				+ Runnable_ClusterModel_Transmission.RUNNABLE_FIELD_TRANSMISSION_RISK_CATEGORIES_BY_CASUAL_PARTNERS]);

		long seed = System.currentTimeMillis();
		if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED])) {
			seed = Long.parseLong(
					loadedProperties.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));
		}

		reallocateRiskGrp(prealloactedRiskGrpArr, baseContactMapSeed, cumulative_pop_composition, riskCatListAll,
				baseDir, seed);
	}

	public static void reallocateRiskGrp(ArrayList<Number[]> prealloactedRiskGrpArr, long baseContactMapSeed,
			int[] cumulative_pop_composition, float[][] riskCatListAll, File baseDir, long rng_seed) {
		File pre_allocate_risk_file = new File(baseDir,
				String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, baseContactMapSeed));

		RandomGenerator rngBase = new MersenneTwisterRandomGenerator(rng_seed);

		// Reset risk group for everyone
		for (Number[] ent : prealloactedRiskGrpArr) {
			ent[PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP] = -1;
		}

		for (float[] riskCatList : riskCatListAll) {
			if (riskCatList != null && riskCatList[0] < 0) {
				// riskCatList: Number[]
				// {-genderIncl_Index, numRiskGrp, numCasualPartnerCat
				// numCasualPartCatUpperRange_0, numCasualPartCatUpperRange_1 ...
				// riskGrp0_Cat0, riskGrp0_Cat1, ...
				// riskGrp1_Cat0, riskGrp1_Cat1, ...}

				int genderIncl = -(int) riskCatList[0];
				int numRiskGrp = (int) riskCatList[1];
				int numCPartCat = (int) riskCatList[2];
				float[] cPartCatUpper = Arrays.copyOfRange(riskCatList, 3, 3 + numCPartCat - 1);
				int unitDistTotal = 0;

				for (int i = 0; i < numRiskGrp * numCPartCat; i++) {
					unitDistTotal += (int) riskCatList[i + 3 + (numCPartCat - 1)];
				}

				int[] riskGrpIndex = new int[unitDistTotal];
				int[] catGrpIndex = new int[unitDistTotal];
				int rPt = 0;
				for (int rG = 0; rG < numRiskGrp; rG++) {
					for (int catG = 0; catG < numCPartCat; catG++) {
						int rep = (int) riskCatList[3 + (numCPartCat - 1) + rG * numCPartCat + catG];
						while (rep > 0) {
							riskGrpIndex[rPt] = rG;
							catGrpIndex[rPt] = catG;
							rPt++;
							rep--;
						}
					}
				}

				ArrayList<Number[]> risk_grp_prealloc_subGrp = new ArrayList<>();
				for (Number[] ent : prealloactedRiskGrpArr) {
					int g = Arrays.binarySearch(cumulative_pop_composition,
							(Integer) ent[PRE_ALLOCATE_RISK_GRP_INDEX_PID]);
					if (g < 0) {
						g = ~g;
					}
					if ((genderIncl & 1 << g) > 0) {
						risk_grp_prealloc_subGrp.add(ent);
					}
				}

				HashMap<Integer, ArrayList<Number[]>> bins = new HashMap<>();

				for (Number[] ent : risk_grp_prealloc_subGrp) {
					float numCasual12Months = (Float) ent[PRE_ALLOCATE_RISK_GRP_INDEX_CASUAL_RATE_TOTAL];
					if (numCasual12Months > 0) { // Only those with casual partners
						int rGrp = Arrays.binarySearch(cPartCatUpper, numCasual12Months);
						if (rGrp < 0) {
							rGrp = ~rGrp;
						}
						ArrayList<Number[]> bin_ent = bins.get(rGrp);
						if (bin_ent == null) {
							bin_ent = new ArrayList<>();
							bins.put(rGrp, bin_ent);
						}
						bin_ent.add(ent);
					}
				}

				boolean endSelection = false;
				Number[][] pickedEntry = new Number[unitDistTotal][];
				int numberAllocated = 0;

				while (!endSelection) {
					Arrays.fill(pickedEntry, null);
					rPt = 0;
					while (rPt < pickedEntry.length && !endSelection) {
						int binSize = bins.get((catGrpIndex[rPt])).size();
						if (binSize > 0) {
							pickedEntry[rPt] = bins.get(catGrpIndex[rPt]).remove(rngBase.nextInt(binSize));
						} else {
							endSelection = true;
						}
						rPt++;
					}

					if (!endSelection) {
						for (int r = 0; r < pickedEntry.length; r++) {
							pickedEntry[r][PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP] = riskGrpIndex[r];
						}
						numberAllocated += pickedEntry.length;
					}

				}

				System.out.printf("CMap_Seed %d:  Risk group allocated for %.1f%% of the population\n",
						baseContactMapSeed, (100.0f * numberAllocated) / prealloactedRiskGrpArr.size());

				// Set the remainder group
				for (Number[] ent : risk_grp_prealloc_subGrp) {
					if (ent[PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP].intValue() < 0) {
						ent[PRE_ALLOCATE_RISK_GRP_INDEX_RISKGRP] = numRiskGrp;
					}
				}

			}
		}
		try {
			PrintWriter pWri_riskGrp = new PrintWriter(pre_allocate_risk_file);
			for (Number[] ent : prealloactedRiskGrpArr) {
				pWri_riskGrp.printf("%d,%d,%f,%f\n", ent[0], ent[1], ent[2], ent[3]);
			}
			pWri_riskGrp.close();
		} catch (IOException ex) {
			ex.printStackTrace(System.err);
		}
	}

	private void zipOutputFiles() throws IOException, FileNotFoundException {

		zipSelectedOutputs(FILENAME_INDEX_CASE_LIST.replaceFirst("%d", Long.toString(baseContactMapSeed)),
				String.format(FILENAME_INDEX_CASE_LIST_ZIP, baseContactMapSeed));

		if ((simSetting & 1 << SIM_SETTING_KEY_GEN_PREVAL_FILE) != 0) {
			zipSelectedOutputs(FILENAME_PREVALENCE_SITE.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_PREVALENCE_SITE_ZIP, baseContactMapSeed));

			zipSelectedOutputs(FILENAME_PREVALENCE_PERSON.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_PREVALENCE_PERSON_ZIP, baseContactMapSeed));

		}
		if ((simSetting & 1 << SIM_SETTING_KEY_GEN_INCIDENCE_FILE) != 0) {
			zipSelectedOutputs(FILENAME_CUMUL_INCIDENCE_SITE.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_CUMUL_INCIDENCE_SITE_ZIP, baseContactMapSeed));
			zipSelectedOutputs(FILENAME_CUMUL_INCIDENCE_PERSON.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_CUMUL_INCIDENCE_PERSON_ZIP, baseContactMapSeed));

		}
		if ((simSetting & 1 << SIM_SETTING_KEY_GEN_TREATMENT_FILE) != 0) {
			zipSelectedOutputs(FILENAME_CUMUL_POSITIVE_DX_PERSON.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_CUMUL_POSITIVE_DX_PERSON_ZIP, baseContactMapSeed));
			zipSelectedOutputs(
					FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON_ZIP, baseContactMapSeed));
			zipSelectedOutputs(FILENAME_CUMUL_TREATMENT_PERSON.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_CUMUL_TREATMENT_PERSON_ZIP, baseContactMapSeed));

		}
		if ((simSetting & 1 << SIM_SETTING_KEY_TRACK_INFECTION_HISTORY) != 0) {
			zipSelectedOutputs(FILENAME_INFECTION_HISTORY.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_INFECTION_HISTORY_ZIP, baseContactMapSeed));

		}
		if ((simSetting & 1 << SIM_SETTING_KEY_TRACK_ANTIBIOTIC_USAGE) != 0) {
			zipSelectedOutputs(FILENAME_CUMUL_ANTIBIOTIC_USAGE.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_CUMUL_ANTIBIOTIC_USAGE_ZIP, baseContactMapSeed));

		}

		if ((simSetting & 1 << SIM_SETTING_KEY_TRACK_TRANSMISSION_CLUSTER) != 0) {
			zipSelectedOutputs(FILENAME_ALL_TRANSMISSION_CMAP.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_ALL_TRANSMISSION_CMAP_ZIP, baseContactMapSeed));
		}

		if ((simSetting & 1 << SIM_SETTING_KEY_TRACK_VACCINE_COVERAGE) != 0) {
			zipSelectedOutputs(FILENAME_VACCINE_COVERAGE.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_VACCINE_COVERAGE_ZIP, baseContactMapSeed));

			zipSelectedOutputs(FILENAME_VACCINE_COVERAGE_PERSON.replaceFirst("%d", Long.toString(baseContactMapSeed)),
					String.format(FILENAME_VACCINE_COVERAGE_PERSON_ZIP, baseContactMapSeed));
		}
	}

	protected void zipSelectedOutputs(String file_name, String zip_file_name)
			throws IOException, FileNotFoundException {
		final Pattern pattern_include_file = Pattern.compile(file_name.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)"));

		File[] files_list = baseDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				Matcher m = pattern_include_file.matcher(pathname.getName());
				return m.matches();
			}
		});

		if (files_list.length > 0) {
			File zipFile = new File(baseDir, zip_file_name);

			File preZip = null;

			if (zipFile.exists()) {
				// Delete old zip
				File[] oldZips = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.getName().startsWith(zip_file_name + "_") && pathname.getName().endsWith(".7z");
					}
				});
				for (File f : oldZips) {
					Files.delete(f.toPath());
				}
				preZip = new File(baseDir, zip_file_name + "_" + Long.toString(System.currentTimeMillis()) + ".7z");
				Files.copy(zipFile.toPath(), preZip.toPath(), StandardCopyOption.COPY_ATTRIBUTES);
			}

			SevenZOutputFile outputZip = new SevenZOutputFile(zipFile);

			// Copy previous entries from zip

			if (preZip != null) {
				SevenZFile inputZip = new SevenZFile(preZip);
				SevenZArchiveEntry inputEnt;
				final int BUFFER = 2048;
				byte[] buf = new byte[BUFFER];
				while ((inputEnt = inputZip.getNextEntry()) != null) {
					outputZip.putArchiveEntry(inputEnt);
					int count;
					while ((count = inputZip.read(buf, 0, BUFFER)) != -1) {
						outputZip.write(Arrays.copyOf(buf, count));
					}
					outputZip.closeArchiveEntry();
				}
				inputZip.close();

				if (exportSkipBackup) {
					preZip.delete();
				}
			}

			// Add new entry to zip
			SevenZArchiveEntry entry;
			FileInputStream fIn;

			for (int fI = 0; fI < files_list.length; fI++) {
				entry = outputZip.createArchiveEntry(files_list[fI], files_list[fI].getName());
				outputZip.putArchiveEntry(entry);
				fIn = new FileInputStream(files_list[fI]);
				outputZip.write(fIn);
				outputZip.closeArchiveEntry();
				fIn.close();
			}

			outputZip.close();

			// Clean up
			for (File f : files_list) {
				f.delete();
			}
		}
	}

	public static void launch(String[] args) throws IOException, InterruptedException {
		File baseDir = null;
		File propFile = null;

		File[] preGenClusterMap = new File[0];

		if (args.length > 0) {
			baseDir = new File(args[0]);
		} else {
			System.out.println(String.format("Usage: java %s PROP_FILE_DIRECTORY <%s>",
					Simulation_ClusterModelTransmission.class.getName(), LAUNCH_ARGS_SKIP_BACKUP));
			System.exit(0);
		}

		if (baseDir != null) {
			if (baseDir.isDirectory()) {
				// Reading of PROP file
				propFile = new File(baseDir, SimulationInterface.FILENAME_PROP);
				FileInputStream fIS = new FileInputStream(propFile);
				Properties prop = new Properties();
				prop.loadFromXML(fIS);
				fIS.close();
				System.out.println(String.format("Properties file < %s > loaded.", propFile.getAbsolutePath()));

				// Check for contact cluster generated

				final String REGEX_STR = FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");

				File contactMapDir = baseDir;

				if (prop.getProperty(PROP_CONTACT_MAP_LOC) != null) {
					contactMapDir = new File(prop.getProperty(PROP_CONTACT_MAP_LOC));
					if (!contactMapDir.exists() || !contactMapDir.isDirectory()) {
						contactMapDir = baseDir;
					}
				}

				preGenClusterMap = contactMapDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isFile() && Pattern.matches(REGEX_STR, pathname.getName());

					}
				});

				Arrays.sort(preGenClusterMap);

				long tic = System.currentTimeMillis();

				for (int i = 0; i < preGenClusterMap.length; i++) {

					System.out.println(String.format("Running simulation set based on ContactMap located at %s",
							preGenClusterMap[i].getAbsolutePath()));

					StringWriter cMap_str = new StringWriter();
					PrintWriter pWri = new PrintWriter(cMap_str);

					BufferedReader reader = new BufferedReader(new FileReader(preGenClusterMap[i]));
					String line;

					while ((line = reader.readLine()) != null) {
						pWri.println(line);
					}

					pWri.close();
					reader.close();

					Simulation_ClusterModelTransmission sim = new Simulation_ClusterModelTransmission();
					sim.setBaseDir(baseDir);
					sim.loadProperties(prop);
					Matcher m = Pattern.compile(REGEX_STR).matcher(preGenClusterMap[i].getName());
					if (m.matches()) {
						long cMap_seed = Long.parseLong(m.group(1));
						sim.setBaseContactMapSeed(cMap_seed);
					}

					sim.setBaseContactMap(ContactMap.ContactMapFromFullString(cMap_str.toString()));

					if (args.length > 1) {
						for (int ai = 1; ai < args.length; ai++) {
							if (LAUNCH_ARGS_SKIP_BACKUP.equals(args[ai])) {
								sim.setExportSkipBackup(true);
							}
							if (LAUNCH_ARGS_PRINT_PROGRESS.equals(args[ai])) {
								sim.setPrintProgress(true);
							}

						}

					}

					sim.generateOneResultSet();

				}

				System.out.println(String.format("%d simulation(s) completed. Runtime (total)= %.2fs",
						preGenClusterMap.length, (System.currentTimeMillis() - tic) / 1000f));

			}
		}

	}

}
