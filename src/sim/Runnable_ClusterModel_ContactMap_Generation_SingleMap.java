package sim;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;

import population.Population_Bridging;
import population.Population_Bridging_Scheduled;
import relationship.ContactMap;

public class Runnable_ClusterModel_ContactMap_Generation_SingleMap
		extends Abstract_Runnable_ClusterModel_ContactMap_Generation {

	private Population_Bridging population;

	public Runnable_ClusterModel_ContactMap_Generation_SingleMap(long mapSeed) {
		super(mapSeed);
	}

	public Population_Bridging getPopulation() {
		return population;
	}

	public void setPopulation(Population_Bridging population) {
		this.population = population;
	}

	private void exportPopSnap(long snapTime) {
		try {

			File exportFile = new File(baseDir,
					String.format(EXPORT_POP_FILENAME, population.getSeed(), population.getGlobalTime()));

			if (!exportFile.exists()) {
				ObjectOutputStream outStream = new ObjectOutputStream(
						new BufferedOutputStream(new FileOutputStream(exportFile)));
				population.encodePopToStream(outStream);
				outStream.close();

				// Remove old snapshot file

				String fileCheck = String.format(EXPORT_POP_FILENAME, population.getSeed(), 0);
				final String fileCheckPrefix = fileCheck.substring(0, fileCheck.indexOf('.') - 1);

				File[] oldPopulationSnapFiles = baseDir.listFiles(new FileFilter() {
					@Override
					public boolean accept(File pathname) {
						return pathname.getName().startsWith(fileCheckPrefix) && !pathname.equals(exportFile);
					}
				});

				for (File df : oldPopulationSnapFiles) {

					try {
						Files.delete(df.toPath());
					} catch (Exception e) {
						System.err.printf("Error in deleteing %s. Skipped deleted", df.toPath().toString());
					}
				}
			}

		} catch (IOException ex) {
			ex.printStackTrace(System.err);

		}
	}

	@Override
	public void run() {

		long tic = System.currentTimeMillis();
		int[] contactMapValidRange = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE];
		// Check for previous snapshot

		String regEx_PopSnap = EXPORT_POP_FILENAME.replaceFirst("%d", Long.toString(population.getSeed()))
				.replaceFirst("%d", "(\\\\d+)");

		final Pattern pattern_pop_snap = Pattern.compile(regEx_PopSnap);

		File[] oldPopulationSnapFiles;
		int skipTimeUntil = -1;
		final float prob_no_bridge = (float) runnable_fields[RUNNABLE_FIELD_NO_BRIDGE];

		oldPopulationSnapFiles = baseDir.listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return Pattern.matches(regEx_PopSnap, pathname.getName());
			}
		});

		if (oldPopulationSnapFiles.length == 0) {
			if (printStatus != null) {
				population.setPrintStatus(printStatus);
			}

			if (contactMapValidRange[0] == 0) {
				gen_cMap = (ContactMap[]) population.getFields()[Population_Bridging.FIELD_CONTACT_MAP];
				for (int i = 0; i < gen_cMap.length; i++) {
					gen_cMap[i] = new ContactMap();
				}
			}

			if (population instanceof Population_Bridging_Scheduled) {
				((Population_Bridging_Scheduled) population).setProb_no_bridge(prob_no_bridge);
			}

			population.initialise();
		} else {

			if (oldPopulationSnapFiles.length > 1) {
				Arrays.sort(oldPopulationSnapFiles, new Comparator<File>() {
					@Override
					public int compare(File o1, File o2) {
						Matcher m1 = pattern_pop_snap.matcher(o1.getName());
						m1.matches();
						long t1 = Long.parseLong(m1.group(1));
						Matcher m2 = pattern_pop_snap.matcher(o2.getName());
						m2.matches();
						long t2 = Long.parseLong(m2.group(1));

						return Long.compare(t1, t2);
					}
				});
			}

			boolean popFileloaded = false;
			for (int fI = oldPopulationSnapFiles.length - 1; fI >= 0 && !popFileloaded; fI--) {
				File prevSnapFile = oldPopulationSnapFiles[fI];

				StringBuilder output = new StringBuilder();
				output.append(String.format("Reusing population snapshot file %s....", prevSnapFile.getName()));

				try {
					ObjectInputStream objIn;
					Population_Bridging snapPop;

					try {
						objIn = new ObjectInputStream(new BufferedInputStream(new FileInputStream(prevSnapFile)));

						if (population instanceof Population_Bridging_Scheduled) {
							snapPop = Population_Bridging_Scheduled.decodeFromStream(objIn);
						} else {
							snapPop = Population_Bridging.decodeFromStream(objIn);
						}
					} catch (Exception e) {
						// See if temp file existed and load that instead
						File tempFile = new File(baseDir, String.format("%s_temp", prevSnapFile.getName()));
						if (tempFile.exists()) {
							output.append(String.format(" FAILED\nRetry using population snapshot file %s....",
									tempFile.getName()));
							objIn = new ObjectInputStream(new BufferedInputStream(new FileInputStream(tempFile)));
							if (population instanceof Population_Bridging_Scheduled) {
								snapPop = Population_Bridging_Scheduled.decodeFromStream(objIn);
							} else {
								snapPop = Population_Bridging.decodeFromStream(objIn);
							}
						} else {
							throw e;
						}
					}

					snapPop.setBaseDir(baseDir);
					this.setPopulation(snapPop);
					objIn.close();
					output.append(String.format(" SUCCESS, with global time set at %d.", population.getGlobalTime()));
					skipTimeUntil = population.getGlobalTime();
					gen_cMap = (ContactMap[]) population.getFields()[Population_Bridging.FIELD_CONTACT_MAP];

					if (population instanceof Population_Bridging_Scheduled) {
						((Population_Bridging_Scheduled) population)
								.setProb_no_bridge(prob_no_bridge);
					}

					if (contactMapValidRange[0] < population.getGlobalTime()
							&& population.getGlobalTime() < contactMapValidRange[1]) {
						for (int i = 0; i < gen_cMap.length; i++) {
							if (gen_cMap[i] == null) {
								gen_cMap[i] = new ContactMap();
							}
						}
					}

					popFileloaded = true;

				} catch (Exception e) {

					output.append(" FAILED. Trying the next snapshot file (if any).");
					FileUtils.deleteQuietly(prevSnapFile);
				}

				System.out.println(output.toString());

			}

			if (!popFileloaded) {
				System.out.println("Reusing population snapshot FAILED. Initialising a new population instead.");

				if (printStatus != null) {
					population.setPrintStatus(printStatus);
				}
				if (contactMapValidRange[0] == 0) {
					gen_cMap = (ContactMap[]) population.getFields()[Population_Bridging.FIELD_CONTACT_MAP];
					for (int i = 0; i < gen_cMap.length; i++) {
						gen_cMap[i] = new ContactMap();
					}
				}

				if (population instanceof Population_Bridging_Scheduled) {
					((Population_Bridging_Scheduled) population).setProb_no_bridge(prob_no_bridge);
				}
				population.initialise();

			}

			if (printStatus != null) {
				population.setPrintStatus(printStatus);
			}
		}

		long lastExportTime = System.currentTimeMillis();

		int stepCount = 0;

		for (int s = 0; s < numSnaps; s++) {
			for (int f = 0; f < snap_dur; f++) {
				if (stepCount >= skipTimeUntil) {
					if (contactMapValidRange[0] != 0 && population.getGlobalTime() == contactMapValidRange[0]) {
						gen_cMap = (ContactMap[]) population.getFields()[Population_Bridging.FIELD_CONTACT_MAP];
						for (int i = 0; i < gen_cMap.length; i++) {
							gen_cMap[i] = new ContactMap();
						}
					} else if (population.getGlobalTime() == contactMapValidRange[1]) {
						// Set to null
						population.setParameter("Population_Bridging.FIELD_CONTACT_MAP",
								Population_Bridging.FIELD_CONTACT_MAP, new ContactMap[gen_cMap.length]);

					}
					population.advanceTimeStep(1);

					boolean exportPop = true;

					if (population instanceof population.Population_Bridging_Scheduled) {
						exportPop &= !((Population_Bridging_Scheduled) population).isSpace_save();
					}

					if (prob_no_bridge > 0 && exportPop) {
						long snapTime = System.currentTimeMillis();
						if (snapTime - lastExportTime > prob_no_bridge) {
							exportPopSnap(snapTime);
							lastExportTime = snapTime;
						}

					}
				}
				stepCount++;

			}
		}

		if (printStatus != null) {
			for (PrintStream out : printStatus) {
				out.printf("Run time (simulation only) = %.3f seconds\n", (System.currentTimeMillis() - tic) / 1000f);
			}
		}

		if (baseDir != null) {

			boolean exportPop = true;

			if (population instanceof population.Population_Bridging_Scheduled) {
				exportPop &= !((Population_Bridging_Scheduled) population).isSpace_save();
			}
			if (exportPop) {
				exportPopSnap(System.currentTimeMillis());
			}

			if (printStatus != null) {
				for (PrintStream out : printStatus) {
					out.printf("Runtime (incl. export pop) = %.3f seconds\n",
							(System.currentTimeMillis() - tic) / 1000f);
				}
			}

		}

		if (printStatus != null) {
			for (PrintStream out : printStatus) {
				out.close();
			}
		}

	}

	@Override
	public void setRunnable_fields(Object[] simFields) {
		for (int f = 0; f < Abstract_Runnable_ClusterModel_ContactMap_Generation.LENGTH_RUNNABLE_MAP_GEN_FIELD; f++) {
			if (simFields[Population_Bridging.LENGTH_FIELDS_BRIDGING_POP
					+ Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD + f] != null) {
				getRunnable_fields()[f] = simFields[f + Simulation_ClusterModelGeneration.LENGTH_SIM_MAP_GEN_FIELD
						+ Population_Bridging.LENGTH_FIELDS_BRIDGING_POP];
			}
		}

	}

}
