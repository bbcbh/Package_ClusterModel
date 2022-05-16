package sim;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;

import person.AbstractIndividualInterface;
import population.Population_Bridging;
import relationship.ContactMap;

public class Runnable_ContactMapGeneration implements Runnable {

	public static final Object[] DEFAULT_RUNNABLE_MAP_GEN_FIELDS = {
			// RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE
			new int[] { 30, 30 + AbstractIndividualInterface.ONE_YEAR_INT },
			// RUNNABLE_FILED_EXPORT_FREQ - in milliseconds
			-1l,

	};

	public static final int RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE = 0;
	public static final int RUNNABLE_FILED_EXPORT_FREQ = RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE + 1;
	public static final int LENGTH_RUNNABLE_MAP_GEN_FIELD = RUNNABLE_FILED_EXPORT_FREQ + 1;

	public static final String EXPORT_POP_FILENAME = "EXPORT_POP_%d_%d.obj";

	private Population_Bridging population;
	private int numSnaps;
	private int snapFreq;
	private Object[] runnable_fields = new Object[LENGTH_RUNNABLE_MAP_GEN_FIELD];
	private ContactMap[] gen_cMap = null;
	private PrintStream printStatus = null;

	private File baseDir = null;

	public Runnable_ContactMapGeneration() {
		super();
		for (int i = 0; i < DEFAULT_RUNNABLE_MAP_GEN_FIELDS.length; i++) {
			runnable_fields[i] = DEFAULT_RUNNABLE_MAP_GEN_FIELDS[i];
		}
	}

	public void setBaseDir(File baseDir) {
		this.baseDir = baseDir;
	}

	public Object[] getRunnable_fields() {
		return runnable_fields;
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

	public void setPrintStatus(PrintStream printStatus) {
		this.printStatus = printStatus;
	}

	public Population_Bridging getPopulation() {
		return population;
	}

	public void setPopulation(Population_Bridging population) {
		this.population = population;
	}

	@Override
	public void run() {

		if (printStatus != null) {
			population.setPrintStatus(printStatus);
		}
		int[] contactMapValidRange = (int[]) runnable_fields[RUNNABLE_FIELD_CONTACT_MAP_GEN_VALID_RANGE];

		if (contactMapValidRange[0] == 0) {
			gen_cMap = (ContactMap[]) population.getFields()[Population_Bridging.FIELD_CONTACT_MAP];
			for (int i = 0; i < gen_cMap.length; i++) {
				gen_cMap[i] = new ContactMap();
			}

		}

		population.initialise();
		/*
		 * if (printStatus != null) { printStatus.println();
		 * printStatus.println(population.printCurrentPartnershipStatus()); }
		 */

		long lastExportTime = System.currentTimeMillis();
		final long exportFreq = (long) runnable_fields[RUNNABLE_FILED_EXPORT_FREQ];

		for (int s = 0; s < numSnaps; s++) {
			for (int f = 0; f < snapFreq; f++) {
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

						try {
							File exportFile = new File(baseDir,
									String.format(EXPORT_POP_FILENAME, population.getSeed(), snapTime));
							ObjectOutputStream outStream = new ObjectOutputStream(new FileOutputStream(exportFile));														
							population.encodePopToStream(outStream);
							outStream.close();
							
							// Remove old snapshot file							
							
							String fileCheck = 
									String.format(EXPORT_POP_FILENAME,population.getSeed(),0);							
							final String fileCheckPrefix = fileCheck.substring(0, fileCheck.length()-1);
							
							File[] oldSnap = baseDir.listFiles(new FileFilter() {								
								@Override
								public boolean accept(File pathname) {																											
									return pathname.getName().startsWith(fileCheckPrefix) && 
										!pathname.equals(exportFile);
								}
							});
							
							for(File df : oldSnap) {
								Files.delete(df.toPath());
							}
 							
						
							

						} catch (IOException ex) {
							ex.printStackTrace(System.err);

						}

						lastExportTime = snapTime;
					}

				}

			}
		}

		if (baseDir != null) {

			ContactMap cMap = gen_cMap[0];
			File allContactFile = new File(baseDir,
					String.format(Simulation_ClusterModelGeneration.FILENAME_FORMAT_ALL_CMAP, population.getSeed(),
							cMap.vertexSet().size()));

			BufferedWriter fileWriAll;
			try {
				fileWriAll = new BufferedWriter(new FileWriter(allContactFile));
				fileWriAll.append(cMap.toFullString());
				fileWriAll.close();
			} catch (IOException e) {
				e.printStackTrace(System.err);
				System.out.println(cMap.toFullString());
			}

		}

		if (printStatus != null) {
			printStatus.close();
		}

	}

}
