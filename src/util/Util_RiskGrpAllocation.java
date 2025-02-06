package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.InvalidPropertiesFormatException;
import java.util.Properties;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import population.Population_Bridging;
import relationship.ContactMap;
import sim.Abstract_Runnable_ClusterModel_ContactMap_Generation;
import sim.SimulationInterface;
import sim.Simulation_ClusterModelGeneration;
import sim.Simulation_ClusterModelTransmission;

public class Util_RiskGrpAllocation {
	public static final Pattern PATTERN_CMAP_FILE = Pattern
			.compile(Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP.replaceAll("%d", "(-?\\\\d+)"));

	public static void generateRiskGrpAllocationByRiskCat(File cMapFile, float[][] riskCatListAll, int[] map_time_range,
			int[] cumul_pop, long baseSeed) throws FileNotFoundException, IOException {
		Matcher m = PATTERN_CMAP_FILE.matcher(cMapFile.getName());
		m.matches();
		long cMapSeed = Long.parseLong(m.group(1));
		StringWriter outputMsg = new StringWriter();
		PrintWriter pWri = new PrintWriter(outputMsg);
		File riskGrpFile = new File(cMapFile.getParentFile(),
				String.format(Simulation_ClusterModelTransmission.FILENAME_PRE_ALLOCATE_RISK_GRP, cMapSeed));
		if (riskGrpFile.exists()) {
			System.out.printf("Warning: %s already exists. Risk group file not generated.\n", riskGrpFile.getName());
		} else {
			pWri.printf("Reading %s...\n", cMapFile.getName());
			StringBuilder cMap_str = new StringBuilder();
			BufferedReader reader = new BufferedReader(new FileReader(cMapFile));
			String s;
			while ((s = reader.readLine()) != null) {
				cMap_str.append(s);
				cMap_str.append('\n');
			}
			reader.close();

			ContactMap cMap = ContactMap.ContactMapFromFullString(cMap_str.toString());

			pWri.printf("CMap generated from %s.\n", cMapFile.getName());

			ArrayList<Number[]> prealloactedRiskGrpArr = new ArrayList<>();
			Simulation_ClusterModelTransmission.fillRiskGrpArrByCasualPartnership(prealloactedRiskGrpArr, cMap,
					cumul_pop, riskCatListAll, map_time_range);

			pWri.printf("Casual partnership calculated from %s.\n", cMapFile.getName());

			Simulation_ClusterModelTransmission.reallocateRiskGrp(prealloactedRiskGrpArr, cMapSeed, cumul_pop,
					riskCatListAll, cMapFile.getParentFile(), baseSeed);

			pWri.printf("RiskGrp allocation from %s completed.\n\n", cMapFile.getName());
			
			
			pWri.close();
			
			pWri = new PrintWriter(new File(cMapFile.getParent(), String.format("RiskGrp_Output_%s.txt", cMapFile.getName())));
			pWri.println(outputMsg.toString());
			pWri.close();
			

		}
	}

	public static void generateRiskGrpAllocationByRiskCat(File cMapFolder, float[][] riskCatListAll, int numThreads)
			throws FileNotFoundException, IOException, InvalidPropertiesFormatException, InterruptedException {

		File propPath = new File(cMapFolder, "simSpecificSim.prop");
		FileInputStream fIS = new FileInputStream(propPath);
		Properties prop = new Properties();
		prop.loadFromXML(fIS);
		fIS.close();

		int[] pop_composition = (int[]) PropValUtils
				.propStrToObject(prop.getProperty(Simulation_ClusterModelTransmission.POP_PROP_INIT_PREFIX
						+ Integer.toString(Population_Bridging.FIELD_POP_COMPOSITION)), int[].class);

		int[] cumul_pop = Arrays.copyOf(pop_composition, pop_composition.length);
		for (int i = 1; i < cumul_pop.length; i++) {
			cumul_pop[i] = pop_composition[i] + cumul_pop[i - 1];
		}

		final String time_rangeKey = Simulation_ClusterModelGeneration.POP_PROP_INIT_PREFIX
				+ Integer.toString(Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
						+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
						+ Abstract_Runnable_ClusterModel_ContactMap_Generation.RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE);

		int[] map_time_range = (int[]) PropValUtils.propStrToObject(prop.getProperty(time_rangeKey), int[].class);

		File[] cMapFiles = cMapFolder.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				Matcher m = PATTERN_CMAP_FILE.matcher(pathname.getName());
				return m.matches();
			}
		});

		long baseSeed = Long
				.parseLong(prop.getProperty(SimulationInterface.PROP_NAME[SimulationInterface.PROP_BASESEED]));

		ExecutorService execReadMap = Executors
				.newFixedThreadPool(Math.min(numThreads, Runtime.getRuntime().availableProcessors()));

		for (File cMapFile : cMapFiles) {
			if (numThreads <= 1) {
				generateRiskGrpAllocationByRiskCat(cMapFile, riskCatListAll, map_time_range, cumul_pop, baseSeed);
			} else {
				execReadMap.submit(new Runnable() {
					@Override
					public void run() {
						try {
							generateRiskGrpAllocationByRiskCat(cMapFile, riskCatListAll, map_time_range, cumul_pop,
									baseSeed);
						} catch (IOException e) {
							e.printStackTrace(System.err);
						}
					}
				});
			}
		}

		if (numThreads > 1) {
			execReadMap.shutdown();
			if (!execReadMap.awaitTermination(2, TimeUnit.DAYS)) {
				System.err.println("Thread time-out!");
			}
		}

	}

}
