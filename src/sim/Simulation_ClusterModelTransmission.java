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
	public static final String FILENAME_VACCINE_COVERAGE_PERSON_ZIP = FILENAME_VACCINE_COVERAGE_PERSON.replaceFirst("_%d", "")
			+ ".7z";

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

		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD
				+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD;

		float[][] seedInfectParam = (float[][]) simFields[sim_offset + SIM_FIELD_SEED_INFECTION];
		int[][] seedInfectNum = new int[Population_Bridging.LENGTH_GENDER][Runnable_ClusterModel_Transmission.LENGTH_SITE];

		int[] contactMapTimeRange = (int[]) simFields[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE];

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
				
				if(printProgress) {
					runnable[s].setPrint_progress(System.out);
					runnable[s].setRunnableId(String.format("%d,%d",baseContactMapSeed, simSeed));
				}
				

				for (int f = 0; f < Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; f++) {
					int ent_offset = sim_offset + LENGTH_SIM_MAP_TRANSMISSION_FIELD;
					if (simFields[ent_offset + f] != null) {
						runnable[s].getRunnable_fields()[f] = simFields[ent_offset + f];
					}
				}

				if ((simSetting & 1 << SIM_SETTING_KEY_TRACK_TRANSMISSION_CLUSTER) != 0) {
					((Runnable_ClusterModel_Transmission_Map) runnable[s]).setTransmissionMap(new TransmissionMap());
				}
				runnable[s].initialse();

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
			zipSelectedOutputs(FILENAME_CUMUL_POSITIVE_DX_SOUGHT_PERSON.replaceFirst("%d", Long.toString(baseContactMapSeed)),
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
				
				if(prop.getProperty(PROP_CONTACT_MAP_LOC) != null) {
					contactMapDir = new File(prop.getProperty(PROP_CONTACT_MAP_LOC));
					if(!contactMapDir.exists() || !contactMapDir.isDirectory()) {
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
						sim.setBaseContactMapSeed(Long.parseLong(m.group(1)));
					}
					sim.setBaseContactMap(ContactMap.ContactMapFromFullString(cMap_str.toString()));

					if (args.length > 1) {
						for (int ai = 1; ai < args.length; ai++) {
							if (LAUNCH_ARGS_SKIP_BACKUP.equals(args[ai])) {
								sim.setExportSkipBackup(true);
							}
							if(LAUNCH_ARGS_PRINT_PROGRESS.equals(args[ai])) {
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
