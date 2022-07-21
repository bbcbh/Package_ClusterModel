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

import population.Population_Bridging;
import random.MersenneTwisterRandomGenerator;
import random.RandomGenerator;
import relationship.ContactMap;
import util.PersonClassifier;
import util.PropValUtils;

public class Simulation_ClusterModelTransmission implements SimulationInterface {

	protected Properties loadedProperties = null; // From .prop file, if any
	protected boolean stopNextTurn = false;
	protected File baseDir = null;
	protected ContactMap baseContactMap;
	protected long baseContactMapSeed = -1;

	public static final String POP_PROP_INIT_PREFIX = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX;
	public static final String POP_PROP_INIT_PREFIX_CLASS = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX_CLASS;
	public static final String FILENAME_FORMAT_ALL_CMAP = Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP;
	public static final String FILENAME_TRANSMAP_ZIP_PREFIX = "All_transmap_%d";
	public static final String REGEX_TRANSMISSION_CMAP_DIR = Runnable_ClusterModel_Transmission_ContactMap.DIRNAME_FORMAT_TRANSMISSION_CMAP
			.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");
	public static final String REGEX_INDEX_CASE_LIST = Runnable_ClusterModel_Transmission.FILENAME_FORMAT_INDEX_CASE_LIST
			.replaceAll("%d", "(-{0,1}(?!0)\\\\d+)");

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
			+ LENGTH_SIM_MAP_TRANSMISSION_FIELD + Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD];
	public Class<?>[] simFieldClass = new Class[simFields.length];

	public Simulation_ClusterModelTransmission() {
		final int sim_offset = Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
				+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
				+ Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD;
		for (int i = 0; i < simFields.length; i++) {
			// All simulation levels
			if (i >= sim_offset) {
				simFields[i] = DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS[i - sim_offset];
				simFieldClass[i] = DEFAULT_BRIDGING_MAP_TRANS_SIM_FIELDS[i - sim_offset].getClass();
			}
		}
	}

	public void setBaseContactMapSeed(long baseContactMapSeed) {
		this.baseContactMapSeed = baseContactMapSeed;
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
		int numSnap = 1;
		int snapFreq = 1;
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
				numSnap = Integer.parseInt(
						loadedProperties.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_NUM_SNAP]));
			}
			if (loadedProperties.containsKey(SimulationInterface.PROP_NAME[SimulationInterface.PROP_SNAP_FREQ])) {
				snapFreq = Integer.parseInt(loadedProperties
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
		File clusterExport7z = new File(baseDir,
				String.format(FILENAME_TRANSMAP_ZIP_PREFIX, baseContactMapSeed) + ".7z");
		ArrayList<Long> completedSeed = new ArrayList<>();

		if (clusterExport7z.exists()) {
			SevenZFile inputZip = new SevenZFile(clusterExport7z);
			Pattern p_seedIndexList = Pattern.compile(REGEX_INDEX_CASE_LIST);
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

		for (int s = 0; s < numSim; s++) {
			long simSeed = rngBase.nextLong();
			boolean runSim = true;
			if (Collections.binarySearch(completedSeed, simSeed) >= 0) {
				System.out.println(String.format("Transmission map under seed #%d already generated.", simSeed));
				runSim = false;
			}

			if (runSim) {
				runnable[s] = new Runnable_ClusterModel_Transmission_ContactMap(simSeed, pop_composition, baseContactMap,
						numSnap, snapFreq);
				runnable[s].setBaseDir(baseDir);
				
				for(int f = 0; f < Runnable_ClusterModel_Transmission.LENGTH_RUNNABLE_MAP_TRANSMISSION_FIELD; f++) {
					int ent_offset = sim_offset+LENGTH_SIM_MAP_TRANSMISSION_FIELD;
					if(simFields[ent_offset + f] != null) {
						runnable[s].getRunnable_fields()[f] = simFields[ent_offset + f];
					}
				}
				
				
				((Runnable_ClusterModel_Transmission_ContactMap) runnable[s]).setTransmissionMap(new ContactMap());
				runnable[s].initialse();
				

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
									firstContactTime = Math.min(firstContactTime,
											e[Population_Bridging.CONTACT_MAP_EDGE_START_TIME]);
								}
								
								if(firstContactTime == contactMapTimeRange[1]) {
									firstContactTime = contactMapTimeRange[0];
								}

								seedTime[i] = firstContactTime;
								runnable[s].addInfected(infected, site, firstContactTime, firstContactTime + 180);

							}
							System.out.println(String.format("Seeding %s of gender #%d at site #%d at t = %s",
									Arrays.toString(seedInf), gender, site, Arrays.toString(seedTime)));
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
						zipTransmissionMaps();
					}
				} else {
					runnable[s].run();
					zipTransmissionMaps();
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
			zipTransmissionMaps();
		}

	}

	private void zipTransmissionMaps() throws IOException, FileNotFoundException {
		// Zip all transmit output to single file

		File clusterExport7z = new File(baseDir,
				String.format(FILENAME_TRANSMAP_ZIP_PREFIX, baseContactMapSeed) + ".7z");

		File preZip = null;

		if (clusterExport7z.exists()) {
			// Delete old zip
			File[] oldZips = baseDir.listFiles(new FileFilter() {
				@Override
				public boolean accept(File pathname) {
					return pathname.getName()
							.startsWith(String.format(FILENAME_TRANSMAP_ZIP_PREFIX, baseContactMapSeed) + "_")
							&& pathname.getName().endsWith(".7z");
				}
			});
			for (File f : oldZips) {
				Files.delete(f.toPath());
			}
			preZip = new File(baseDir, String.format(FILENAME_TRANSMAP_ZIP_PREFIX, baseContactMapSeed) + "_"
					+ Long.toString(System.currentTimeMillis()) + ".7z");
			Files.copy(clusterExport7z.toPath(), preZip.toPath(), StandardCopyOption.COPY_ATTRIBUTES);

		}

		File[] transMapDirs = baseDir.listFiles(new FileFilter() {

			@Override
			public boolean accept(File pathname) {
				return pathname.isDirectory() && Pattern.matches(REGEX_TRANSMISSION_CMAP_DIR, pathname.getName());

			}
		});

		if (transMapDirs.length > 0) {

			SevenZOutputFile outputZip = new SevenZOutputFile(clusterExport7z);

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
			}

			File[] entFiles;
			SevenZArchiveEntry entry;
			FileInputStream fIn;

			for (int dI = 0; dI < transMapDirs.length; dI++) {
				entFiles = transMapDirs[dI].listFiles();
				for (int fI = 0; fI < entFiles.length; fI++) {
					entry = outputZip.createArchiveEntry(entFiles[fI], entFiles[fI].getName());
					outputZip.putArchiveEntry(entry);
					fIn = new FileInputStream(entFiles[fI]);
					outputZip.write(fIn);
					outputZip.closeArchiveEntry();
					fIn.close();
				}

			}
			outputZip.close();
			// Clean up
			for (int dI = 0; dI < transMapDirs.length; dI++) {
				FileUtils.deleteDirectory(transMapDirs[dI]);
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
			System.out.println(String.format("Usage: java %s PROP_FILE_DIRECTORY",
					Simulation_ClusterModelTransmission.class.getName()));
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

				preGenClusterMap = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.isFile() && Pattern.matches(REGEX_STR, pathname.getName());

					}
				});

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
					sim.generateOneResultSet();

				}

				System.out.println(String.format("%d simulation(s) completed. Runtime (total)= %.2fs",
						preGenClusterMap.length, (System.currentTimeMillis() - tic) / 1000f));

			}
		}

	}
}
