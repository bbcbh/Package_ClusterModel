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

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import population.Population_Bridging_Scheduled;
import relationship.ContactMap;

public class Runnable_ClusterModel_ContactMap_Generation extends Abstract_Runnable_ClusterModel {

	public static final Object[] DEFAULT_RUNNABLE_MAP_GEN_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE
			new int[] { 30, 30 + AbstractIndividualInterface.ONE_YEAR_INT },
			// RUNNABLE_FILED_EXPORT_FREQ - in milliseconds
			-1l,

	};

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE = 0;
	public static final int RUNNABLE_FILED_EXPORT_FREQ = RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE + 1;
	public static final int LENGTH_RUNNABLE_MAP_GEN_FIELD = RUNNABLE_FILED_EXPORT_FREQ + 1;

	public static final String EXPORT_POP_FILENAME_PRRFIX = "EXPORT_POP_%d";
	public static final String EXPORT_POP_FILENAME = EXPORT_POP_FILENAME_PRRFIX+ "_%d.obj";

	private Population_Bridging population;
	private int numSnaps;
	private int snapFreq;
	private Object[] runnable_fields = new Object[LENGTH_RUNNABLE_MAP_GEN_FIELD];
	private ContactMap[] gen_cMap = null;
	private PrintStream[] printStatus = null;

	public Runnable_ClusterModel_ContactMap_Generation() {
		super();
		for (int i = 0; i < DEFAULT_RUNNABLE_MAP_GEN_FIELDS.length; i++) {
			runnable_fields[i] = DEFAULT_RUNNABLE_MAP_GEN_FIELDS[i];
		}
	}

	public void setNumSnaps(int numSnaps) {
		this.numSnaps = numSnaps;
	}

	public ContactMap[] getGen_cMap() {
		return gen_cMap;
	}

	public void setSnapFreq(int snapFreq) {
		this.snapFreq = snapFreq;
	}

	public void setPrintStatus(PrintStream[] printStatus) {
		this.printStatus = printStatus;
	}

	public Population_Bridging getPopulation() {
		return population;
	}

	public void setPopulation(Population_Bridging population) {
		this.population = population;
	}

	@Override
	public Object[] getRunnable_fields() {
		return runnable_fields;
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
		final long exportFreq = (long) runnable_fields[RUNNABLE_FILED_EXPORT_FREQ];

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
				((Population_Bridging_Scheduled) population).setExport_period_form_partnership_progress(exportFreq);
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
						((Population_Bridging_Scheduled) population).setExport_period_form_partnership_progress(exportFreq);
					}
					
					if(contactMapValidRange[0] < population.getGlobalTime()
							&&  population.getGlobalTime() < contactMapValidRange[1]) {
						for(int i = 0; i < gen_cMap.length; i++) {
							if(gen_cMap[i] == null) {
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
					((Population_Bridging_Scheduled) population).setExport_period_form_partnership_progress(exportFreq);
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
			for (int f = 0; f < snapFreq; f++) {
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

					if (exportFreq > 0) {
						long snapTime = System.currentTimeMillis();
						if (snapTime - lastExportTime > exportFreq) {
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
			
			exportPopSnap(System.currentTimeMillis());

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

}
